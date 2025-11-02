// pmms.hpp
// Parallel Multiway Merge Sort (templated, header-only).
// + LoserTree merge, PMMS stats, overloads for rvalue/lvalue,
// + Parallel 3-way quicksort (templated) with OpenMP tasks.

#ifndef PMMS_HPP_INCLUDED
#define PMMS_HPP_INCLUDED

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <queue>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#ifdef __linux__
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#endif

namespace pmms {

using std::size_t;
using std::vector;
using hr_clock = std::chrono::high_resolution_clock;

struct PMMSOptions {
    bool usePWayMerge = true;
    bool useSampling = true;
    int oversample = 4;
    bool pinThreads = false;
    int maxSelectIterations = 64;
};

struct PMMSStats {
    double total = 0.0;
    double local_sort = 0.0;
    double splitter_selection = 0.0;
    double partitioning = 0.0;
    double merge = 0.0;
};

inline void tryPinThread(bool doPin) {
#ifdef __linux__
    if (!doPin)
        return;
    int tid = omp_get_thread_num();
    int nprocs = static_cast<int>(sysconf(_SC_NPROCESSORS_ONLN));
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(tid % nprocs, &cpuset);
    pthread_t t = pthread_self();
    pthread_setaffinity_np(t, sizeof(cpuset), &cpuset);
#else
    (void) doPin;
#endif
}

// ---------------- sampling splitters ----------------
template <typename T>
vector<T> sampleSplitters(const vector<vector<T>>& sequences, int p,
                          int oversample) {
    vector<T> samples;
    if (sequences.empty() || p <= 1)
        return samples;
    samples.reserve(static_cast<size_t>(sequences.size()) *
                    std::max(1, oversample * p));
    for (const auto& seq : sequences) {
        int len = static_cast<int>(seq.size());
        if (len == 0)
            continue;
        int s = std::min(len, std::max(1, oversample * p));
        for (int t = 1; t <= s; ++t) {
            int idx =
                static_cast<int>((static_cast<long long>(t) * len) / (s + 1));
            if (idx < 0)
                idx = 0;
            if (idx >= len)
                idx = len - 1;
            samples.push_back(seq[static_cast<size_t>(idx)]);
        }
    }
    if (samples.empty())
        return {};
    std::sort(samples.begin(), samples.end());
    vector<T> splitters;
    splitters.reserve(std::max(0, p - 1));
    size_t Ssz = samples.size();
    for (int j = 1; j < p; ++j) {
        size_t pos = (static_cast<size_t>(j) * Ssz) / p;
        if (pos >= Ssz)
            pos = Ssz - 1;
        splitters.push_back(samples[pos]);
    }
    for (size_t i = 1; i < splitters.size(); ++i)
        if (splitters[i] < splitters[i - 1])
            splitters[i] = splitters[i - 1];
    return splitters;
}

// ---------------- exact selection safe ----------------
template <typename T>
T msSelectExactSafe(const vector<vector<T>>& sequences, long long k,
                    int maxIterations = 64) {
    const int p = static_cast<int>(sequences.size());
    if (p == 0)
        throw std::invalid_argument("msSelectExactSafe: empty input");
    vector<long long> L(p), R(p);
    long long total = 0;
    for (int i = 0; i < p; ++i) {
        L[i] = 0;
        R[i] = static_cast<long long>(sequences[i].size()) - 1;
        if (R[i] >= L[i])
            total += (R[i] - L[i] + 1);
    }
    if (k < 1 || k > total)
        throw std::out_of_range("msSelectExactSafe: k out of range");

    int iter = 0;
    while (true) {
        long long rem = 0;
        for (int i = 0; i < p; ++i)
            if (R[i] >= L[i])
                rem += (R[i] - L[i] + 1);
        if (rem == 1) {
            for (int i = 0; i < p; ++i)
                if (R[i] >= L[i])
                    return sequences[i][static_cast<size_t>(L[i])];
        }
        vector<T> reps;
        reps.reserve(p);
        for (int i = 0; i < p; ++i) {
            if (R[i] >= L[i]) {
                long long mid = (L[i] + R[i]) / 2;
                reps.push_back(sequences[i][static_cast<size_t>(mid)]);
            }
        }
        if (reps.empty())
            throw std::logic_error("msSelectExactSafe: no reps");
        size_t mpos = reps.size() / 2;
        std::nth_element(reps.begin(), reps.begin() + static_cast<long>(mpos),
                         reps.end());
        T pivot = reps[mpos];

        vector<long long> counts(p, 0);
        long long totalLE = 0;
        for (int i = 0; i < p; ++i) {
            if (R[i] < L[i]) {
                counts[i] = 0;
                continue;
            }
            auto itB = sequences[i].begin() + static_cast<size_t>(L[i]);
            auto itE = sequences[i].begin() + static_cast<size_t>(R[i]) + 1;
            auto it = std::upper_bound(itB, itE, pivot);
            counts[i] = static_cast<long long>(it - itB);
            totalLE += counts[i];
        }

        if (totalLE >= k) {
            for (int i = 0; i < p; ++i) R[i] = L[i] + counts[i] - 1;
        } else {
            for (int i = 0; i < p; ++i) L[i] = L[i] + counts[i];
            k -= totalLE;
        }

        ++iter;
        if (iter > maxIterations) {
            vector<T> remElems;
            remElems.reserve(static_cast<size_t>(rem));
            for (int i = 0; i < p; ++i) {
                if (R[i] >= L[i]) {
                    auto b = sequences[i].begin() + static_cast<size_t>(L[i]);
                    auto e =
                        sequences[i].begin() + static_cast<size_t>(R[i]) + 1;
                    remElems.insert(remElems.end(), b, e);
                }
            }
            if ((long long) remElems.size() < k)
                throw std::logic_error(
                    "msSelectExactSafe fallback: size mismatch");
            std::nth_element(remElems.begin(),
                             remElems.begin() + static_cast<size_t>(k - 1),
                             remElems.end());
            return remElems[static_cast<size_t>(k - 1)];
        }
    }
}

// ---------------- LoserTree (segment-tree based) ----------------
template <typename T>
class LoserTree {
    int k;
    int sz;
    vector<int> seg;
    std::vector<std::optional<T>> leaves;
    std::less<T> less;

    inline int minIndex(int a, int b) {
        if (a == -1)
            return b;
        if (b == -1)
            return a;
        return less(leaves[a].value(), leaves[b].value()) ? a : b;
    }

   public:
    LoserTree() : k(0), sz(1) {}
    LoserTree(const std::vector<std::optional<T>>& init) { initTree(init); }

    void initTree(const std::vector<std::optional<T>>& init) {
        leaves = init;
        k = static_cast<int>(leaves.size());
        sz = 1;
        while (sz < k) sz <<= 1;
        seg.assign(2 * sz, -1);
        for (int i = 0; i < k; ++i)
            seg[sz + i] = leaves[i].has_value() ? i : -1;
        for (int i = k; i < sz; ++i) seg[sz + i] = -1;
        for (int node = sz - 1; node >= 1; --node)
            seg[node] = minIndex(seg[2 * node], seg[2 * node + 1]);
    }

    int winner() const {
        if (k == 0)
            return -1;
        return seg.size() > 1 ? seg[1] : -1;
    }

    const std::optional<T>& leafValue(int idx) const {
        return leaves[static_cast<size_t>(idx)];
    }

    void replaceLeaf(int idx, const std::optional<T>& val) {
        if (idx < 0 || idx >= k)
            return;
        leaves[static_cast<size_t>(idx)] = val;
        int pos = sz + idx;
        seg[pos] = val.has_value() ? idx : -1;
        pos >>= 1;
        while (pos >= 1) {
            seg[pos] = minIndex(seg[2 * pos], seg[2 * pos + 1]);
            pos >>= 1;
        }
    }
};

// ---------------- Core merge (LoserTree) ----------------
template <typename T>
vector<T> mergeSortedSequencesInPlace(vector<vector<T>>& sequences, int p,
                                      const PMMSOptions& opt,
                                      PMMSStats* stats = nullptr) {
    const int sCount = static_cast<int>(sequences.size());
    long long total = 0;
    for (const auto& s : sequences) total += static_cast<long long>(s.size());
    if (total == 0)
        return vector<T>{};
    if (p <= 0)
        p = omp_get_max_threads();
    if (p > total)
        p = static_cast<int>(std::min<long long>(total, p));

    auto t_start_total = hr_clock::now();

    // splitter selection
    auto t0 = hr_clock::now();
    vector<T> splitters;
    if (opt.useSampling) {
        splitters =
            sampleSplitters<T>(sequences, p, std::max(1, opt.oversample));
    } else {
        splitters.reserve(static_cast<size_t>(p > 0 ? p - 1 : 0));
        for (int j = 1; j < p; ++j) {
            long long k = ((long long) j * total) / p;
            if (k < 1)
                k = 1;
            splitters.push_back(
                msSelectExactSafe<T>(sequences, k, opt.maxSelectIterations));
        }
    }
    auto t1 = hr_clock::now();
    if (stats)
        stats->splitter_selection =
            std::chrono::duration<double>(t1 - t0).count();

    // partitioning
    t0 = hr_clock::now();
    vector<vector<long long>> pos(
        static_cast<size_t>(sCount),
        vector<long long>(static_cast<size_t>(p + 1), 0));
    vector<long long> bucketTotal(static_cast<size_t>(p), 0LL);
    for (int i = 0; i < sCount; ++i) {
        pos[i][0] = 0;
        size_t curpos = 0;
        for (int j = 0; j < p - 1; ++j) {
            bool hasSplitter = (j < static_cast<int>(splitters.size()));
            auto it = hasSplitter
                          ? std::upper_bound(sequences[i].begin() + curpos,
                                             sequences[i].end(),
                                             splitters[static_cast<size_t>(j)])
                          : sequences[i].end();
            size_t idx = static_cast<size_t>(it - sequences[i].begin());
            pos[i][j + 1] = static_cast<long long>(idx);
            curpos = idx;
        }
        pos[i][p] = static_cast<long long>(sequences[i].size());
        for (int j = 0; j < p; ++j)
            bucketTotal[static_cast<size_t>(j)] += (pos[i][j + 1] - pos[i][j]);
    }
    vector<long long> baseBucket(static_cast<size_t>(p), 0LL);
    for (int j = 1; j < p; ++j)
        baseBucket[static_cast<size_t>(j)] =
            baseBucket[static_cast<size_t>(j - 1)] +
            bucketTotal[static_cast<size_t>(j - 1)];
    t1 = hr_clock::now();
    if (stats)
        stats->partitioning = std::chrono::duration<double>(t1 - t0).count();

    // allocate result
    vector<T> out(static_cast<size_t>(total));

    // merge with LoserTree per bucket
    t0 = hr_clock::now();
#pragma omp parallel for schedule(dynamic) num_threads(p)
    for (int j = 0; j < p; ++j) {
        tryPinThread(opt.pinThreads);
        long long outPos = baseBucket[static_cast<size_t>(j)];
        std::vector<long long> curs(static_cast<size_t>(sCount));
        for (int i = 0; i < sCount; ++i)
            curs[static_cast<size_t>(i)] =
                pos[static_cast<size_t>(i)][static_cast<size_t>(j)];
        std::vector<std::optional<T>> init(static_cast<size_t>(sCount));
        for (int i = 0; i < sCount; ++i) {
            long long b = pos[static_cast<size_t>(i)][static_cast<size_t>(j)];
            long long e =
                pos[static_cast<size_t>(i)][static_cast<size_t>(j + 1)];
            if (b < e)
                init[static_cast<size_t>(i)] =
                    sequences[static_cast<size_t>(i)][static_cast<size_t>(b)];
            else
                init[static_cast<size_t>(i)].reset();
        }
        LoserTree<T> lt(init);
        while (true) {
            int win = lt.winner();
            if (win == -1)
                break;
            const T val = lt.leafValue(win).value();
            out[static_cast<size_t>(outPos++)] = val;
            curs[static_cast<size_t>(win)] += 1;
            long long nextIdx = curs[static_cast<size_t>(win)];
            long long endIdx =
                pos[static_cast<size_t>(win)][static_cast<size_t>(j + 1)];
            if (nextIdx < endIdx)
                lt.replaceLeaf(win, sequences[static_cast<size_t>(win)]
                                             [static_cast<size_t>(nextIdx)]);
            else
                lt.replaceLeaf(win, std::optional<T>());
        }
    }
    t1 = hr_clock::now();
    if (stats)
        stats->merge = std::chrono::duration<double>(t1 - t0).count();

    auto t_end_total = hr_clock::now();
    if (stats)
        stats->total =
            std::chrono::duration<double>(t_end_total - t_start_total).count();

    return out;
}

// Public overloads
template <typename T>
vector<T> parallelMultiwayMergeSort(const vector<vector<T>>& sequences_input,
                                    int p,
                                    const PMMSOptions& opt = PMMSOptions(),
                                    PMMSStats* stats = nullptr) {
    vector<vector<T>> seqs = sequences_input;
    return mergeSortedSequencesInPlace<T>(seqs, p, opt, stats);
}

template <typename T>
vector<T> parallelMultiwayMergeSort(vector<vector<T>>&& sequences_mv, int p,
                                    const PMMSOptions& opt = PMMSOptions(),
                                    PMMSStats* stats = nullptr) {
    return mergeSortedSequencesInPlace<T>(sequences_mv, p, opt, stats);
}

// single-vector API
template <typename T>
vector<T> parallelMultiwayMergeSort(const vector<T>& data, int p,
                                    const PMMSOptions& opt = PMMSOptions(),
                                    PMMSStats* stats = nullptr) {
    const int n = static_cast<int>(data.size());
    if (n == 0) {
        if (stats) {
            stats->total = 0.0;
            stats->local_sort = stats->splitter_selection =
                stats->partitioning = stats->merge = 0.0;
        }
        return {};
    }
    if (p <= 0)
        p = omp_get_max_threads();
    if (p > n)
        p = n;

    vector<vector<T>> sequences(static_cast<size_t>(p));
    int baseChunk = n / p, rem = n % p;
    int cur = 0;
    for (int i = 0; i < p; ++i) {
        int len = baseChunk + (i < rem ? 1 : 0);
        if (len > 0)
            sequences[i].assign(data.begin() + cur, data.begin() + cur + len);
        cur += len;
    }

    auto t_start_total = hr_clock::now();

    // local sorts
    auto t0 = hr_clock::now();
#pragma omp parallel for schedule(static) num_threads(p)
    for (int i = 0; i < p; ++i) {
        tryPinThread(opt.pinThreads);
        std::sort(sequences[i].begin(), sequences[i].end());
    }
    auto t1 = hr_clock::now();
    if (stats)
        stats->local_sort = std::chrono::duration<double>(t1 - t0).count();

    PMMSStats coreStats;
    PMMSStats* corePtr = stats ? &coreStats : nullptr;
    vector<T> result =
        mergeSortedSequencesInPlace<T>(sequences, p, opt, corePtr);

    if (stats) {
        stats->splitter_selection = coreStats.splitter_selection;
        stats->partitioning = coreStats.partitioning;
        stats->merge = coreStats.merge;
        auto t_end_total = hr_clock::now();
        stats->total =
            std::chrono::duration<double>(t_end_total - t_start_total).count();
    }

    return result;
}

/* =========================
   Parallel 3-way quicksort
   Generic templated implementation for any T and Compare.
   Interfaces:
     - sequential_quicksort_3way_rec<T,Compare>(T* a, int lo, int hi, Compare
   comp, int cutoff)
     - parallel_quicksort(vector<T>& v, int num_threads, Compare comp, int
   cutoff)
   ========================= */
template <typename T, typename Compare = std::less<T>>
void sequential_quicksort_3way_rec(T* a, int lo, int hi,
                                   Compare comp = Compare(),
                                   int cutoff = 16384) {
    if (hi - lo <= 1)
        return;
    if (hi - lo <= cutoff) {
        std::sort(a + lo, a + hi, comp);
        return;
    }
    T pivot = a[lo];
    int lt = lo, gt = hi - 1, i = lo + 1;
    while (i <= gt) {
        if (comp(a[i], pivot))
            std::swap(a[lt++], a[i++]);
        else if (comp(pivot, a[i]))
            std::swap(a[i], a[gt--]);
        else
            ++i;
    }
    sequential_quicksort_3way_rec(a, lo, lt, comp, cutoff);
    sequential_quicksort_3way_rec(a, gt + 1, hi, comp, cutoff);
}

template <typename T, typename Compare = std::less<T>>
void parallel_quicksort_task(T* a, int lo, int hi, Compare comp = Compare(),
                             int cutoff = 16384) {
    if (hi - lo <= 1)
        return;
    if (hi - lo <= cutoff) {
        std::sort(a + lo, a + hi, comp);
        return;
    }
    T pivot = a[lo];
    int lt = lo, gt = hi - 1, i = lo + 1;
    while (i <= gt) {
        if (comp(a[i], pivot))
            std::swap(a[lt++], a[i++]);
        else if (comp(pivot, a[i]))
            std::swap(a[i], a[gt--]);
        else
            ++i;
    }
#pragma omp task shared(a) if (lt - lo > cutoff)
    parallel_quicksort_task<T, Compare>(a, lo, lt, comp, cutoff);
#pragma omp task shared(a) if (hi - (gt + 1) > cutoff)
    parallel_quicksort_task<T, Compare>(a, gt + 1, hi, comp, cutoff);
#pragma omp taskwait
}

template <typename T, typename Compare = std::less<T>>
void parallel_quicksort(std::vector<T>& v,
                        int num_threads = omp_get_max_threads(),
                        Compare comp = Compare(), int cutoff = 16384) {
    if (v.size() <= 1)
        return;
    omp_set_num_threads(num_threads);
#pragma omp parallel
    {
#pragma omp single nowait
        {
            parallel_quicksort_task<T, Compare>(
                v.data(), 0, static_cast<int>(v.size()), comp, cutoff);
        }
    }
}

}   // namespace pmms

#endif   // PMMS_HPP_INCLUDED
