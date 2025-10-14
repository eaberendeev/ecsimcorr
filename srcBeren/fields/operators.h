#pragma once

#include "Mesh.h"
// Функция сравнения для сортировки по ключу (row, col)
inline bool compareTriplets(const Triplet& a, const Triplet& b) {
    return std::tie(a.row(), a.col()) < std::tie(b.row(), b.col());
}

template <typename IndexerX, typename IndexerY, typename IndexerZ,
          typename MatrixType>
void Mesh::convert_block_to_crs_format(MatrixType bmatrix, Operator& mat,
                                       const Domain& domain) {
    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;

    // std::vector<Triplet> trips;
    const auto size = domain.size();
    const int max_i = size.x() - 1;
    const int max_j = size.y() - 1;
    const int max_k = size.z() - 1;
    const int num_threads = std::min(omp_get_max_threads(), 128);

    // Оценка количества итераций и элементов для каждого потока.
    // Если общее число итераций равно:
    //   (max_i - BORDER) * (max_j - BORDER) * (max_k - BORDER)
    // и для каждой итерации добавляется 9 элементов (3x3 окрестность),
    // то можно примерно вычислить estimated_per_thread.
    size_t total_iterations = static_cast<size_t>(max_i - BORDER) *
                              (max_j - BORDER) * (max_k - BORDER);
    size_t estimated_per_thread = (total_iterations / num_threads) * 9;

    // Локальные векторы для каждого потока с предварительным резервированием
    // памяти.
    std::vector<std::vector<Triplet>> local_vectors(num_threads);
    for (int t = 0; t < num_threads; ++t) {
        local_vectors[t].reserve(estimated_per_thread);
    }
    double time1 = omp_get_wtime();
#pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();

        std::vector<Triplet> local_trips;
        local_trips.reserve(144 * 9 * (max_i - BORDER) / omp_get_num_threads());

#pragma omp for collapse(3) schedule(dynamic, 8) nowait
        for (int i = BORDER; i < max_i; ++i)
            for (int j = BORDER; j < max_j; ++j)
                for (int k = BORDER; k < max_k; ++k) {
                    if (!bmatrix.non_zeros[sind(i, j, k)])
                        continue;
                    const auto& block = bmatrix[sind(i, j, k)];

                    // X component
                    processComponent<IndexerX, IndexerX, 0>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<IndexerX, IndexerY, 1>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<IndexerX, IndexerZ, 2>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);

                    // Y component
                    processComponent<IndexerY, IndexerX, 3>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<IndexerY, IndexerY, 4>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<IndexerY, IndexerZ, 5>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);

                    // Z component
                    processComponent<IndexerZ, IndexerX, 6>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<IndexerZ, IndexerY, 7>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<IndexerZ, IndexerZ, 8>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                }

        // Сортировка локального вектора по (row, col)
        auto& vec = local_vectors[tid];
        std::sort(vec.begin(), vec.end(), compareTriplets);

        // Устранение дубликатов: проход по отсортированному вектору, складываем
        // значения для одинаковых ключей
        if (!vec.empty()) {
            size_t index = 0;
            for (size_t j = 1; j < vec.size(); ++j) {
                if (vec[index].row() == vec[j].row() &&
                    vec[index].col() == vec[j].col()) {
                    vec[index].value() += vec[j].value();
                } else {
                    ++index;
                    vec[index] = vec[j];
                }
            }
            vec.resize(index + 1);
        }
    }
    double time2 = omp_get_wtime();

    // Собираем только непустые локальные векторы для дальнейшего слияния
    std::vector<std::vector<Triplet>> non_empty;
    for (auto& v : local_vectors) {
        if (!v.empty()) {
            non_empty.push_back(std::move(v));
        }
    }

    // Многофазное слияние: объединяем пары отсортированных векторов параллельно
    // с предварительным резервированием памяти
    while (non_empty.size() > 1) {
        size_t new_size = (non_empty.size() + 1) / 2;
        std::vector<std::vector<Triplet>> new_vectors(new_size);
        size_t pairs = non_empty.size() / 2;

#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < pairs; ++i) {
            const auto& left = non_empty[2 * i];
            const auto& right = non_empty[2 * i + 1];
            // Резервируем память для слияния двух векторов
            size_t merged_capacity = left.size() + right.size();
            std::vector<Triplet> merged;
            merged.reserve(merged_capacity);

            size_t li = 0, ri = 0;
            while (li < left.size() && ri < right.size()) {
                if (compareTriplets(left[li], right[ri])) {
                    merged.push_back(left[li]);
                    ++li;
                } else if (compareTriplets(right[ri], left[li])) {
                    merged.push_back(right[ri]);
                    ++ri;
                } else {
                    // Если ключи равны – складываем значения
                    Triplet t = left[li];
                    t.value() += right[ri].value();
                    merged.push_back(t);
                    ++li;
                    ++ri;
                }
            }
            while (li < left.size()) {
                merged.push_back(left[li++]);
            }
            while (ri < right.size()) {
                merged.push_back(right[ri++]);
            }
            new_vectors[i] = std::move(merged);
        }
        // Если число векторов нечётное – последний переносим без изменений
        if (non_empty.size() % 2 == 1) {
            new_vectors.back() = std::move(non_empty.back());
        }
        non_empty = std::move(new_vectors);
    }
    double time3 = omp_get_wtime();

    // В non_empty[0] теперь находится глобальный вектор, уже отсортированный и
    // с устранёнными дубликатами.
    const std::vector<Triplet>& trips =
        non_empty.empty() ? std::vector<Triplet>() : non_empty.front();

    // std::vector<Trip> trips2;
    // trips2.reserve(trips.size());
    // for (long i = 0; i < trips.size(); i++){
    //     trips2.emplace_back(trips[i].row, trips[i].col, trips[i].value);
    // }

    // double time4 = omp_get_wtime();

    // mat.setFromTriplets(trips.begin(), trips.end());
    optimizedSetFromTriplets(mat, trips);
    double time5 = omp_get_wtime();
    std::cout << time2 - time1 << " " << time3 - time2 << " " << time5 - time3
              << "\n";
    std::cout << "Matrix L (block) was created."
              << " trips size: " << trips.size() << std::endl;
}
