// tests.cpp
// Tests and microbench for pmms library with detailed timings.
// Also contains micro-benchmark comparing priority_queue vs LoserTree.
// Build: g++ -O3 -std=c++17 -fopenmp tests.cpp -o tests

#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <queue>
#include <random>
#include <string>

#include "pmms.hpp"

// DRAFT
using namespace pmms;

class Triplet {
   public:
    Triplet(int r, int c, double v) : _row(r), _col(c), _value(v) {}
    Triplet() : _row(0), _col(0), _value(0) {}
    const int& row() const noexcept { return _row; }
    const int& col() const noexcept { return _col; }
    const double& value() const noexcept { return _value; }
    int& row() { return _row; }
    int& col() { return _col; }
    double& value() { return _value; }
    bool operator<(const Triplet& other) const {
        return std::tie(_row, _col) < std::tie(other.row(), other.col());
    }

   private:
    int _row;
    int _col;
    double _value;
};

static bool equalVecsInt(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i])
            return false;
    return true;
}
static bool equalVecsStr(const std::vector<std::string>& a,
                         const std::vector<std::string>& b) {
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i])
            return false;
    return true;
}
template <typename T>
static bool isSortedAsc(const std::vector<T>& v) {
    return std::is_sorted(v.begin(), v.end(), std::less<T>());
}

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& seqs) {
    long long total = 0;
    for (const auto& s : seqs) total += static_cast<long long>(s.size());
    std::vector<T> out;
    out.reserve(static_cast<size_t>(total));
    for (const auto& s : seqs)
        for (const auto& v : s) out.push_back(v);
    return out;
}
template <typename T>
std::vector<std::vector<T>> splitAndSort(const std::vector<T>& data,
                                         int parts) {
    int n = static_cast<int>(data.size());
    if (parts <= 0)
        parts = omp_get_max_threads();
    std::vector<std::vector<T>> seqs(static_cast<size_t>(parts));
    int base = n / parts, rem = n % parts;
    int cur = 0;
    for (int i = 0; i < parts; ++i) {
        int len = base + (i < rem ? 1 : 0);
        if (len > 0) {
            seqs[i].assign(data.begin() + cur, data.begin() + cur + len);
            std::sort(seqs[i].begin(), seqs[i].end());
        }
        cur += len;
    }
    return seqs;
}

// ------------------ Performance benchmark with detailed timings
// ------------------
void perfBenchmarkDetailedInt(int n, int p) {
    std::mt19937 rng(42);
    std::vector<int> data(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) data[i] = static_cast<int>(rng());

    std::cout << "\nDetailed perf benchmark: n=" << n << " p=" << p << "\n";

    // reference std::sort (single-thread)
    auto t0 = std::chrono::high_resolution_clock::now();
    auto ref = data;
    std::sort(ref.begin(), ref.end());
    auto t1 = std::chrono::high_resolution_clock::now();
    double t_std = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "std::sort (single-thread) : " << t_std << " s\n\n";

    // 1) PMMS single-vector, p-way merge, sampling
    PMMSOptions opt;
    opt.useSampling = true;
    opt.usePWayMerge = true;
    opt.oversample = 4;

    PMMSStats stats1;
    t0 = std::chrono::high_resolution_clock::now();
    auto out1 = parallelMultiwayMergeSort<int>(data, p, opt, &stats1);
    t1 = std::chrono::high_resolution_clock::now();
    double wall1 = std::chrono::duration<double>(t1 - t0).count();
    if (!equalVecsInt(out1, ref)) {
        std::cerr << "pmms sampling p-way incorrect\n";
        std::exit(1);
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "PMMS (single-vector, sampling, p-way):\n";
    std::cout << "  wall_total       = " << wall1 << " s\n";
    std::cout << "  stats.total      = " << stats1.total << " s\n";
    std::cout << "  local_sort       = " << stats1.local_sort << " s\n";
    std::cout << "  splitter_select  = " << stats1.splitter_selection << " s\n";
    std::cout << "  partitioning     = " << stats1.partitioning << " s\n";
    std::cout << "  merge            = " << stats1.merge << " s\n";
    std::cout << "  speedup vs std   = " << (t_std / wall1) << "x\n\n";

    // 2) PMMS single-vector, copy+sort (no p-way heap)
    opt.usePWayMerge = false;
    PMMSStats stats2;
    t0 = std::chrono::high_resolution_clock::now();
    auto out2 = parallelMultiwayMergeSort<int>(data, p, opt, &stats2);
    t1 = std::chrono::high_resolution_clock::now();
    double wall2 = std::chrono::duration<double>(t1 - t0).count();
    if (!equalVecsInt(out2, ref)) {
        std::cerr << "pmms sampling copy+sort incorrect\n";
        std::exit(1);
    }

    std::cout << "PMMS (single-vector, sampling, copy+sort):\n";
    std::cout << "  wall_total       = " << wall2 << " s\n";
    std::cout << "  stats.total      = " << stats2.total << " s\n";
    std::cout << "  local_sort       = " << stats2.local_sort << " s\n";
    std::cout << "  splitter_select  = " << stats2.splitter_selection << " s\n";
    std::cout << "  partitioning     = " << stats2.partitioning << " s\n";
    std::cout << "  merge(copy+sort) = " << stats2.merge << " s\n";
    std::cout << "  speedup vs std   = " << (t_std / wall2) << "x\n\n";

    // 3) Merge pre-sorted sequences (rvalue) path (zero-copy)
    auto seqs = splitAndSort<int>(data, p);
    PMMSStats stats3;
    opt.useSampling = true;
    opt.usePWayMerge = true;
    t0 = std::chrono::high_resolution_clock::now();
    auto out3 =
        parallelMultiwayMergeSort<int>(std::move(seqs), p, opt, &stats3);
    t1 = std::chrono::high_resolution_clock::now();
    double wall3 = std::chrono::duration<double>(t1 - t0).count();
    if (!equalVecsInt(out3, ref)) {
        std::cerr << "pmms merge-sorted rvalue incorrect\n";
        std::exit(1);
    }

    std::cout << "PMMS (merge-sorted, rvalue, sampling, p-way):\n";
    std::cout << "  wall_total       = " << wall3 << " s\n";
    std::cout << "  stats.total      = " << stats3.total << " s\n";
    std::cout << "  splitter_select  = " << stats3.splitter_selection << " s\n";
    std::cout << "  partitioning     = " << stats3.partitioning << " s\n";
    std::cout << "  merge            = " << stats3.merge << " s\n";
    std::cout << "  speedup vs std   = " << (t_std / wall3) << "x\n\n";

    // 4) Parallel quick sort
    t0 = std::chrono::high_resolution_clock::now();
    auto out4 = data;
    parallel_quicksort<int>(out4);
    t1 = std::chrono::high_resolution_clock::now();
    double wall4 = std::chrono::duration<double>(t1 - t0).count();
    if (!equalVecsInt(out4, ref)) {
        std::cerr << "parallel quick sort incorrect\n";
        std::exit(1);
    }

    std::cout << "parallel quick sort):\n";
    std::cout << "  wall_total       = " << wall4 << " s\n";
    std::cout << "  speedup vs std   = " << (t_std / wall4) << "x\n\n";
}

// ------------------ Micro-benchmark: priority_queue vs LoserTree
// ------------------ We compute splitters & partitions once, then run two
// different merge backends on the same partitioning so results are directly
// comparable.
static void compute_partitions_for_int(
    const std::vector<std::vector<int>>& sequences, int p,
    const PMMSOptions& opt, std::vector<int>& splitters,
    std::vector<std::vector<long long>>& pos,
    std::vector<long long>& bucketTotal, std::vector<long long>& baseBucket) {
    int sCount = static_cast<int>(sequences.size());
    long long total = 0;
    for (const auto& s : sequences) total += (long long) s.size();
    if (p <= 0)
        p = omp_get_max_threads();
    if (p > total)
        p = static_cast<int>(std::min<long long>(total, p));

    if (opt.useSampling) {
        splitters =
            sampleSplitters<int>(sequences, p, std::max(1, opt.oversample));
    } else {
        splitters.clear();
        for (int j = 1; j < p; ++j) {
            long long k = ((long long) j * total) / p;
            if (k < 1)
                k = 1;
            splitters.push_back(
                msSelectExactSafe<int>(sequences, k, opt.maxSelectIterations));
        }
    }

    pos.assign(static_cast<size_t>(sCount),
               std::vector<long long>(static_cast<size_t>(p + 1), 0));
    bucketTotal.assign(static_cast<size_t>(p), 0LL);
    for (int i = 0; i < sCount; ++i) {
        pos[i][0] = 0;
        size_t curpos = 0;
        for (int j = 0; j < p - 1; ++j) {
            bool hasSplitter = (j < (int) splitters.size());
            auto it = hasSplitter
                          ? std::upper_bound(sequences[i].begin() + curpos,
                                             sequences[i].end(),
                                             splitters[static_cast<size_t>(j)])
                          : sequences[i].end();
            size_t idx = static_cast<size_t>(it - sequences[i].begin());
            pos[i][j + 1] = (long long) idx;
            curpos = idx;
        }
        pos[i][p] = (long long) sequences[i].size();
        for (int j = 0; j < p; ++j)
            bucketTotal[static_cast<size_t>(j)] += (pos[i][j + 1] - pos[i][j]);
    }

    baseBucket.assign(static_cast<size_t>(p), 0LL);
    for (int j = 1; j < p; ++j)
        baseBucket[static_cast<size_t>(j)] =
            baseBucket[static_cast<size_t>(j - 1)] +
            bucketTotal[static_cast<size_t>(j - 1)];
}

// Merge implementation using priority_queue (parallel over buckets)
static std::vector<int> merge_with_priority_queue(
    const std::vector<std::vector<int>>& sequences,
    const std::vector<std::vector<long long>>& pos,
    const std::vector<long long>& baseBucket,
    const std::vector<long long>& bucketTotal) {
    int sCount = (int) sequences.size();
    int p = (int) bucketTotal.size();
    long long total = 0;
    for (auto v : bucketTotal) total += v;
    std::vector<int> out(static_cast<size_t>(total));

#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, p))
    for (int j = 0; j < p; ++j) {
        long long outPos = baseBucket[static_cast<size_t>(j)];
        struct Node {
            int val;
            int seq;
            long long idx;
        };
        struct Cmp {
            bool operator()(Node const& a, Node const& b) const {
                return a.val > b.val;
            }
        };   // min-heap
        std::priority_queue<Node, std::vector<Node>, Cmp> pq;

        for (int i = 0; i < sCount; ++i) {
            long long b = pos[static_cast<size_t>(i)][static_cast<size_t>(j)];
            long long e =
                pos[static_cast<size_t>(i)][static_cast<size_t>(j + 1)];
            if (b < e)
                pq.push(Node{
                    sequences[static_cast<size_t>(i)][static_cast<size_t>(b)],
                    i, b});
        }
        while (!pq.empty()) {
            Node cur = pq.top();
            pq.pop();
            out[static_cast<size_t>(outPos++)] = cur.val;
            long long next = cur.idx + 1;
            int seqIdx = cur.seq;
            if (next <
                pos[static_cast<size_t>(seqIdx)][static_cast<size_t>(j + 1)]) {
                pq.push(Node{sequences[static_cast<size_t>(seqIdx)]
                                      [static_cast<size_t>(next)],
                             seqIdx, next});
            }
        }
    }
    return out;
}

// Merge implementation using LoserTree (parallel over buckets)
static std::vector<int> merge_with_loser_tree(
    const std::vector<std::vector<int>>& sequences,
    const std::vector<std::vector<long long>>& pos,
    const std::vector<long long>& baseBucket,
    const std::vector<long long>& bucketTotal) {
    int sCount = (int) sequences.size();
    int p = (int) bucketTotal.size();
    long long total = 0;
    for (auto v : bucketTotal) total += v;
    std::vector<int> out(static_cast<size_t>(total));

#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, p))
    for (int j = 0; j < p; ++j) {
        long long outPos = baseBucket[static_cast<size_t>(j)];
        // curs per sequence
        std::vector<long long> curs(static_cast<size_t>(sCount));
        for (int i = 0; i < sCount; ++i)
            curs[static_cast<size_t>(i)] =
                pos[static_cast<size_t>(i)][static_cast<size_t>(j)];

        std::vector<std::optional<int>> init(static_cast<size_t>(sCount));
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

        pmms::LoserTree<int> lt(init);
        while (true) {
            int win = lt.winner();
            if (win == -1)
                break;
            int val = lt.leafValue(win).value();
            out[static_cast<size_t>(outPos++)] = val;
            curs[static_cast<size_t>(win)] += 1;
            long long nextIdx = curs[static_cast<size_t>(win)];
            long long endIdx =
                pos[static_cast<size_t>(win)][static_cast<size_t>(j + 1)];
            if (nextIdx < endIdx) {
                lt.replaceLeaf(win, sequences[static_cast<size_t>(win)]
                                             [static_cast<size_t>(nextIdx)]);
            } else {
                lt.replaceLeaf(win, std::optional<int>());
            }
        }
    }
    return out;
}

// micro-benchmark function
void microBenchmarkMergePQvsLT(int n, int p) {
    std::cout << "\nMicro-benchmark: priority_queue vs LoserTree (merge of "
                 "sorted sequences)\n";
    std::mt19937 rng(1234);
    std::vector<int> data(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        data[i] = rng() % 1000;   // moderate domain for duplicates

    auto seqs = splitAndSort<int>(data, p);
    // compute partitions once
    PMMSOptions opt;
    opt.useSampling = true;
    opt.oversample = 4;
    std::vector<int> splitters;
    std::vector<std::vector<long long>> pos;
    std::vector<long long> bucketTotal, baseBucket;
    compute_partitions_for_int(seqs, p, opt, splitters, pos, bucketTotal,
                               baseBucket);

    // warmup
    auto warm1 = merge_with_priority_queue(seqs, pos, baseBucket, bucketTotal);
    auto warm2 = merge_with_loser_tree(seqs, pos, baseBucket, bucketTotal);
    (void) warm1;
    (void) warm2;

    // measure priority_queue
    double best_pq = 1e300;
    for (int r = 0; r < 3; ++r) {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto out_pq =
            merge_with_priority_queue(seqs, pos, baseBucket, bucketTotal);
        auto t1 = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration<double>(t1 - t0).count();
        best_pq = std::min(best_pq, dt);
        // verify sorted
        auto ref = flatten(seqs);
        std::sort(ref.begin(), ref.end());
        if (!equalVecsInt(out_pq, ref)) {
            std::cerr << "micro-bench PQ result mismatch\n";
            std::exit(1);
        }
    }

    // measure loser-tree
    double best_lt = 1e300;
    for (int r = 0; r < 3; ++r) {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto out_lt = merge_with_loser_tree(seqs, pos, baseBucket, bucketTotal);
        auto t1 = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration<double>(t1 - t0).count();
        best_lt = std::min(best_lt, dt);
        auto ref = flatten(seqs);
        std::sort(ref.begin(), ref.end());
        if (!equalVecsInt(out_lt, ref)) {
            std::cerr << "micro-bench LT result mismatch\n";
            std::exit(1);
        }
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  best priority_queue merge : " << best_pq << " s\n";
    std::cout << "  best loser-tree merge     : " << best_lt << " s\n";
    std::cout << "  ratio (PQ / LT)           : " << (best_pq / best_lt)
              << "x\n";
}

// Keep correctness tests (concise)
void runQuickCorrectness() {
    std::cout << "Running quick correctness checks...\n";
    // small int
    std::mt19937 rng(12345);
    std::vector<int> a(1000);
    for (int i = 0; i < 1000; ++i) a[i] = rng() % 200 - 100;
    PMMSOptions opt;
    auto r1 = parallelMultiwayMergeSort<int>(a, 8, opt);
    auto ref = a;
    std::sort(ref.begin(), ref.end());
    if (!equalVecsInt(r1, ref)) {
        std::cerr << "quick check failed\n";
        std::exit(1);
    }
    std::cout << "quick checks OK\n";
}

int main() {
    int p = omp_get_max_threads();
    std::cout << "OpenMP threads available: " << p << "\n";
    std::cout << "Using p = " << p << " partitions\n";

    runQuickCorrectness();

    // detailed perf test
    int n = 20'000'000;   // tune to your machine
    perfBenchmarkDetailedInt(n, p);

    // micro-benchmark (merge backend comparison)
    microBenchmarkMergePQvsLT(50000'000, p);   // moderate n for micro-bench

    std::cout << "Done.\n";
    return 0;
}
