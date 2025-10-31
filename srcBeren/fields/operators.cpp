#include "Mesh.h"
#include "World.h"
#include "Shape.h"
#include "util.h"

// matrix ColMajor
// (0,0) (0,1) (0,2)
// (1,0) (1,1) (1,2)
// (2,0) (2,1) (2,2)

// Ex(i+/2,j,k), Ey(i,j+1/2,k), Ez(i,j,k+1/2) 
// Bx(i,j+1/2,k+1/2), By(i+1/2,j,k+1/2), Bz(i+/2,j+1/2,k)

 void Mesh::stencil_Imat(Operator &mat, const Domain &domain)
 { 
 // !!!!! needs bound condition and if cases!!!!!!
  std::vector<Trip> trips;
  const auto size = domain.size();
  int totalSize = size.x()*size.y()*size.z();
  trips.reserve(totalSize);

  for(int i = 0; i < size.x(); i++){
    for(int j = 0; j < size.y(); j++){
      for(int k = 0; k < size.z(); k++){

        // i,j,k
        trips.push_back(Trip(vind(i,j,k,0),vind(i,j,k,0),1.0));

        // i,j,k
        trips.push_back(Trip(vind(i,j,k,1),vind(i,j,k,1),1.0));
          
        // i,j,k
        trips.push_back(Trip(vind(i,j,k,2),vind(i,j,k,2),1.0));
          
      }
    }
  }
  std::cout << size << " " << 3*totalSize << " " << Imat.rows() << " " << Imat.cols() <<  "\n";
  mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_Lmat(Operator &mat, const Domain &domain) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    const int rowsCount = 3 * size.x() * size.y() * size.z();
    const size_t totalSize = (size_t)rowsCount * LMAT_MAX_ELEMENTS_PER_ROW;
    std::cout << totalSize << "\n";
    trips.reserve(totalSize);

#pragma omp parallel
    {
        std::vector<Trip> local_trips;
        local_trips.reserve(totalSize / omp_get_max_threads());

#pragma omp for schedule(dynamic, 8)
        for (int row = 0; row < rowsCount; row++) {
            for (const auto &[col, value] : LmatX[row]) {
                if (std::abs(value) > LMAT_VALUE_TOLERANCE)
                    local_trips.emplace_back(row, col, value);
            }
        }

#pragma omp critical
        trips.insert(trips.end(), local_trips.begin(), local_trips.end());
    }
    std::cout << "trips size: " << trips.size() << std::endl;

    mat.setFromTriplets(trips.begin(), trips.end());
}


void Mesh::stencil_Lmat2(Operator& mat, const Domain& domain) {
    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;
    
    //std::vector<Triplet> trips;
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

#pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        // First touch - выделяем память в том же потоке, который будет её
        // использовать
        local_vectors[tid].reserve(estimated_per_thread);
    }
    double time1 = omp_get_wtime();
#pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        
        #pragma omp for collapse(3) schedule(dynamic, 8) nowait
        for (int i = BORDER; i < max_i; ++i)
        for (int j = BORDER; j < max_j; ++j)
        for (int k = BORDER; k < max_k; ++k)
        {
            if (!LmatX2.non_zeros[sind(i, j, k)])
                continue;
            const auto& block = LmatX2[sind(i,j,k)];
            
            // X component
            processComponent<XIndexer, XIndexer, 0>(i,j,k, block, local_vectors[tid], xSize, ySize, zSize, TOL);
            processComponent<XIndexer, YIndexer, 1>(i, j, k, block, local_vectors[tid],
                                                    xSize, ySize, zSize, TOL);
            processComponent<XIndexer, ZIndexer, 2>(i, j, k, block, local_vectors[tid],
                                                    xSize, ySize, zSize, TOL);

            // Y component
            processComponent<YIndexer, XIndexer, 3>(i, j, k, block, local_vectors[tid],
                                                    xSize, ySize, zSize, TOL);
            processComponent<YIndexer, YIndexer, 4>(i, j, k, block, local_vectors[tid],
                                                    xSize, ySize, zSize, TOL);
            processComponent<YIndexer, ZIndexer, 5>(i, j, k, block, local_vectors[tid],
                                                    xSize, ySize, zSize, TOL);

            // Z component
            processComponent<ZIndexer, XIndexer, 6>(i, j, k, block, local_vectors[tid],
                                                    xSize, ySize, zSize, TOL);
            processComponent<ZIndexer, YIndexer, 7>(i, j, k, block, local_vectors[tid],
                                                    xSize, ySize, zSize, TOL);
            processComponent<ZIndexer, ZIndexer, 8>(i, j, k, block, local_vectors[tid],
                                                    xSize, ySize, zSize, TOL);
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

    //mat.setFromTriplets(trips.begin(), trips.end());
    optimizedSetFromTriplets(mat, trips);
        double time5 = omp_get_wtime();
    std::cout << time2-time1  << " " << time3-time2  << " " << time5-time3 << "\n";
    std::cout << "Matrix L (block) was created." << " trips size: " << trips.size() << std::endl;
}

void Mesh::stencil_Lmat2_NGP(Operator &mat, const Domain &domain) {
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
        local_trips.reserve(16 * 9 * (max_i - BORDER) / omp_get_num_threads());

#pragma omp for collapse(3) schedule(dynamic, 8) nowait
        for (int i = BORDER; i < max_i; ++i)
            for (int j = BORDER; j < max_j; ++j)
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX_NGP.non_zeros[sind(i, j, k)])
                        continue;
                    const auto &block = LmatX_NGP[sind(i, j, k)];

                    // X component
                    processComponent<XIndexerNGP, XIndexerNGP, 0>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<XIndexerNGP, YIndexerNGP, 1>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<XIndexerNGP, ZIndexerNGP, 2>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);

                    // Y component
                    processComponent<YIndexerNGP, XIndexerNGP, 3>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<YIndexerNGP, YIndexerNGP, 4>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<YIndexerNGP, ZIndexerNGP, 5>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);

                    // Z component
                    processComponent<ZIndexerNGP, XIndexerNGP, 6>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<ZIndexerNGP, YIndexerNGP, 7>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                    processComponent<ZIndexerNGP, ZIndexerNGP, 8>(
                        i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                        TOL);
                }

        // Сортировка локального вектора по (row, col)
        auto &vec = local_vectors[tid];
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
    for (auto &v : local_vectors) {
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
            const auto &left = non_empty[2 * i];
            const auto &right = non_empty[2 * i + 1];
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
    const std::vector<Triplet> &trips =
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


void Mesh::apply_periodic_boundaries(Operator &LmatX) {
    const auto size = int3(xSize, ySize, zSize);
    constexpr int OVERLAP_SIZE = 3;
    const int last_indx = size.x() - OVERLAP_SIZE;
    const int last_indy = size.y() - OVERLAP_SIZE;
    const int last_indz = size.z() - OVERLAP_SIZE;
    Eigen::SparseMatrix<double, MAJOR> boundaryMatrix(LmatX.rows(),
                                                      LmatX.cols());
    std::vector<Trip> boundaryTrips;
    if (!(bounds.isPeriodic(X) || bounds.isPeriodic(Y) || bounds.isPeriodic(Z))) return;

        boundaryTrips.reserve(size.x() * size.y() * size.z());
    if (bounds.isPeriodic(X)) {
#pragma omp parallel
    {
        std::vector<Trip> localTrips;
        localTrips.reserve(size.x() * size.y() * size.z() /
                           omp_get_num_threads());
#pragma omp for schedule(dynamic, 32)
            for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
                auto ix = pos_vind(i, 0);
                auto iy = pos_vind(i, 1);
                auto iz = pos_vind(i, 2);
                auto id = pos_vind(i, 3);
                if (ix < OVERLAP_SIZE) {
                    auto indBound = vind(last_indx + ix, iy, iz, id);

                    for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(
                             LmatX, i);
                         it; ++it) {
                        auto ind2 = it.col();
                        auto ix1 = pos_vind(ind2, 0);
                        auto iy1 = pos_vind(ind2, 1);
                        auto iz1 = pos_vind(ind2, 2);
                        auto id1 = pos_vind(ind2, 3);

                        if (ix1 < OVERLAP_SIZE) {
                            auto indBound2 =
                                vind(last_indx + ix1, iy1, iz1, id1);
                            localTrips.emplace_back(indBound, indBound2,
                                                    it.value());
                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
                if (ix > last_indx - 1) {
                    auto indBound = vind(ix - last_indx, iy, iz, id);
                    for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(
                             LmatX, i);
                         it; ++it) {
                        auto ind2 = it.col();
                        auto ix1 = pos_vind(ind2, 0);
                        auto iy1 = pos_vind(ind2, 1);
                        auto iz1 = pos_vind(ind2, 2);
                        auto id1 = pos_vind(ind2, 3);
                        if (ix1 > last_indx - 1) {
                            auto indBound2 =
                                vind(ix1 - last_indx, iy1, iz1, id1);
                            localTrips.emplace_back(indBound, indBound2,
                                                    it.value());
                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
            }

#pragma omp critical
        boundaryTrips.insert(boundaryTrips.end(), localTrips.begin(),
                             localTrips.end());
    }
    boundaryMatrix.setFromTriplets(boundaryTrips.begin(), boundaryTrips.end());
    LmatX += boundaryMatrix;
    boundaryTrips.clear();
    }
    if (bounds.isPeriodic(Y)) {
#pragma omp parallel
    {
        std::vector<Trip> localTrips;
        localTrips.reserve(size.x() * size.y() * size.z() /
                           omp_get_num_threads());

#pragma omp for schedule(dynamic, 32)
            for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
                auto iy = pos_vind(i, 1);
                auto ix = pos_vind(i, 0);
                auto iz = pos_vind(i, 2);
                auto id = pos_vind(i, 3);
                if (iy < OVERLAP_SIZE) {
                    auto indBound = vind(ix, last_indy + iy, iz, id);
                    for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(
                             LmatX, i);
                         it; ++it) {
                        auto ind2 = it.col();
                        auto ix1 = pos_vind(ind2, 0);
                        auto iy1 = pos_vind(ind2, 1);
                        auto iz1 = pos_vind(ind2, 2);
                        auto id1 = pos_vind(ind2, 3);

                        if (iy1 < OVERLAP_SIZE) {
                            localTrips.emplace_back(
                                indBound, vind(ix1, last_indy + iy1, iz1, id1),
                                it.value());
                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
                if (iy > last_indy - 1) {
                    auto indBound = vind(ix, iy - last_indy, iz, id);
                    for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(
                             LmatX, i);
                         it; ++it) {
                        auto ind2 = it.col();
                        //auto value = it.value();
                        auto ix1 = pos_vind(ind2, 0);
                        auto iy1 = pos_vind(ind2, 1);
                        auto iz1 = pos_vind(ind2, 2);
                        auto id1 = pos_vind(ind2, 3);

                        if (iy1 > last_indy - 1) {
                            auto indBound2 =
                                vind(ix1, iy1 - last_indy, iz1, id1);
                            localTrips.emplace_back(indBound, indBound2,
                                                    it.value());

                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
            }
#pragma omp critical 
        boundaryTrips.insert(boundaryTrips.end(), localTrips.begin(),
                             localTrips.end());
    }

    boundaryMatrix.setFromTriplets(boundaryTrips.begin(), boundaryTrips.end());
    LmatX += boundaryMatrix;
    boundaryTrips.clear();
    }
    if (bounds.isPeriodic(Z)) {
#pragma omp parallel
    {
        std::vector<Trip> localTrips;
        localTrips.reserve(size.x() * size.y() * size.z() /
                           omp_get_num_threads());
#pragma omp for schedule(dynamic, 32)
            for (int i = 0; i < LmatX.outerSize(); i++) {
                auto iz = pos_vind(i, 2);
                auto ix = pos_vind(i, 0);
                auto iy = pos_vind(i, 1);
                auto id = pos_vind(i, 3);
                if (iz < OVERLAP_SIZE) {
                    auto indBound = vind(ix, iy, last_indz + iz, id);
                    for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(
                             LmatX, i);
                         it; ++it) {
                        auto ind2 = it.col();
                        auto ix1 = pos_vind(ind2, 0);
                        auto iy1 = pos_vind(ind2, 1);
                        auto iz1 = pos_vind(ind2, 2);
                        auto id1 = pos_vind(ind2, 3);

                        if (iz1 < OVERLAP_SIZE) {
                            auto indBound2 =
                                vind(ix1, iy1, last_indz + iz1, id1);
                            localTrips.emplace_back(indBound, indBound2,
                                                    it.value());
                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
                if (iz > last_indz - 1) {
                    auto indBound = vind(ix, iy, iz - last_indz, id);
                    for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(
                             LmatX, i);
                         it; ++it) {
                        auto ind2 = it.col();
                       // auto value = it.value();
                        auto ix1 = pos_vind(ind2, 0);
                        auto iy1 = pos_vind(ind2, 1);
                        auto iz1 = pos_vind(ind2, 2);
                        auto id1 = pos_vind(ind2, 3);

                        if (iz1 > last_indz - 1) {
                            auto indBound2 =
                                vind(ix1, iy1, iz1 - last_indz, id1);
                            localTrips.emplace_back(indBound, indBound2,
                                                    it.value());
                        } else {
                            localTrips.emplace_back(indBound, ind2, it.value());
                        }
                    }
                }
            }

#pragma omp critical
        boundaryTrips.insert(boundaryTrips.end(), localTrips.begin(),
                             localTrips.end());
    }

    boundaryMatrix.setFromTriplets(boundaryTrips.begin(), boundaryTrips.end());
    LmatX += boundaryMatrix;
    boundaryTrips.clear();
    }
}

void Mesh::stencil_curlB(Operator &mat, const Domain &domain) {
    // TO DO: create a different boundary cases
    // NOW X and Y always periodic
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x()*size.y()*size.z()*12;
    trips.reserve(totalSize);
   // if (domain.is_periodic_bound(Z)) {
        stencil_curlB_periodic(trips, domain);
   // } else {
    //    stencil_curlB_openZ(trips, domain);
   // }

    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_curlB_periodic(std::vector<Trip> &trips,
                                  const Domain &domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    const auto size = domain.size();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    auto addTriplet = [](std::vector<Trip> &trips, const Domain &domain,
                         int vindElec, int vindMag, double val) {
        bool onArea = domain.in_region_electric(vindElec);
        if (onArea) {
            trips.push_back(Trip(vindElec, vindMag, val));
        }
    };

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                // const int im = (i != 0) ? i - 1 : size.x() - 4;
                // const int jm = (j != 0) ? j - 1 : size.y() - 4;
                // const int km = (k != 0) ? k - 1 : size.z() - 4;
                const int im = (i != 0) || !domain.is_periodic_bound(X)
                                   ? i - 1
                                   : size.x() - 4;
                const int jm = (j != 0) || !domain.is_periodic_bound(Y)
                                   ? j - 1
                                   : size.y() - 4;
                const int km = (k != 0) || !domain.is_periodic_bound(Z)
                                   ? k - 1
                                   : size.z() - 4;

                const int vindx = vind(i, j, k, 0);
                const int vindy = vind(i, j, k, 1);
                const int vindz = vind(i, j, k, 2);

                // (x)[i+1/2,j,k]
                // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
                double val = 1.0 / dy;
                addTriplet(trips, domain, vindx, vind(i, j, k, 2), val);
                addTriplet(trips, domain, vindx, vind(i, jm, k, 2), -val);
                // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
                val = -1.0 / dz;
                addTriplet(trips, domain, vindx, vind(i, j, k, 1), val);
                addTriplet(trips, domain, vindx, vind(i, j, km, 1), -val);

                // (y)[i,j+1/2,k]
                // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
                val = 1.0 / dz;
                addTriplet(trips, domain, vindy, vind(i, j, k, 0), val);
                addTriplet(trips, domain, vindy, vind(i, j, km, 0), -val);
                // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
                val = -1.0 / dx;
                addTriplet(trips, domain, vindy, vind(i, j, k, 2), val);
                addTriplet(trips, domain, vindy, vind(im, j, k, 2), -val);

                // (z)[i,j,k+1/2]
                // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                val = 1.0 / dx;
                addTriplet(trips, domain, vindz, vind(i, j, k, 1), val);
                addTriplet(trips, domain, vindz, vind(im, j, k, 1), -val);
                // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                val = -1.0 / dy;
                addTriplet(trips, domain, vindz, vind(i, j, k, 0), val);
                addTriplet(trips, domain, vindz, vind(i, jm, k, 0), -val);
            }
        }
    }
}

void Mesh::stencil_curlB_openZ(Operator &mat, const Domain &domain) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 12;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int im = (i != 0) ? i - 1 : size.x() - 4;
                const int jm = (j != 0) ? j - 1 : size.y() - 4;
                const int km = (k != 0) ? k - 1 : size.z() - 4;
                const int vindx = vind(i, j, k, 0);
                const int vindy = vind(i, j, k, 1);
                const int vindz = vind(i, j, k, 2);

                if (k < 2 || k > size.z() - 3) {
                    if (k == 1) {
                        // (z)[i,j,k+1/2]
                        // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                        double val = 1.0 / dx;
                        trips.push_back(Trip(vindz, vind(i, j, k, 1), val));
                        trips.push_back(Trip(vindz, vind(im, j, k, 1), -val));
                        // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                        val = -1.0 / dy;
                        trips.push_back(Trip(vindz, vind(i, j, k, 0), val));
                        trips.push_back(Trip(vindz, vind(i, jm, k, 0), -val));
                    }
                    continue;
                }

                // (x)[i+1/2,j,k]
                // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
                double val = 1.0 / dy;
                trips.push_back(Trip(vindx, vind(i, j, k, 2), val));
                trips.push_back(Trip(vindx, vind(i, jm, k, 2), -val));
                // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
                val = -1.0 / dz;
                trips.push_back(Trip(vindx, vind(i, j, k, 1), val));
                trips.push_back(Trip(vindx, vind(i, j, km, 1), -val));

                // (y)[i,j+1/2,k]
                // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
                val = 1.0 / dz;
                trips.push_back(Trip(vindy, vind(i, j, k, 0), val));
                trips.push_back(Trip(vindy, vind(i, j, km, 0), -val));
                // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
                val = -1.0 / dx;
                trips.push_back(Trip(vindy, vind(i, j, k, 2), val));
                trips.push_back(Trip(vindy, vind(im, j, k, 2), -val));

                // (z)[i,j,k+1/2]
                // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                val = 1.0 / dx;
                trips.push_back(Trip(vindz, vind(i, j, k, 1), val));
                trips.push_back(Trip(vindz, vind(im, j, k, 1), -val));
                // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                val = -1.0 / dy;
                trips.push_back(Trip(vindz, vind(i, j, k, 0), val));
                trips.push_back(Trip(vindz, vind(i, jm, k, 0), -val));
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_curlE(Operator &mat, const Domain &domain) {
    // TO DO: create a different boundary cases
    // NOW X and Y always periodic
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 12;
    trips.reserve(totalSize);
    stencil_curlE_periodic(trips, domain);
    std::cout << trips[0].col() << " " << trips[0].row() << " " << trips.size()
              << " " << totalSize << " " << mat.rows() << std::endl;
    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_curlE_periodic(std::vector<Trip> &trips,
                                  const Domain &domain) {
    const auto size = domain.size();
    double dx = domain.cell_size().x();
    double dy = domain.cell_size().y();
    double dz = domain.cell_size().z();

    auto addTriplet =
        [](std::vector<Trip> &trips, const Domain &domain, int vindMag,
           int vindElec, double val) {
            bool onArea = domain.in_region_magnetic(vindMag);
            if (onArea) {
                trips.push_back(Trip(vindMag, vindElec, val));
                //if (vindMag < 0 || vindElec < 0 || vindMag >= 3* domain.size().x() * domain.size().y() * domain.size().z() ||
                //vindElec >= domain.size().x() * domain.size().y() * domain.size().z())
                 //   std::cout << "(" << vindMag << ", " << vindElec
                  //            << "): " << val << std::endl;
            }
        };

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int ip =
                    (i != size.x() - 1) || !domain.is_periodic_bound(X) ? i + 1
                                                                        : 3;
                const int jp =
                    (j != size.y() - 1) || !domain.is_periodic_bound(Y) ? j + 1
                                                                        : 3;
                const int kp =
                    (k != size.z() - 1) || !domain.is_periodic_bound(Z) ? k + 1
                                                                        : 3;
                const int vindx = vind(i, j, k, 0);
                const int vindy = vind(i, j, k, 1);
                const int vindz = vind(i, j, k, 2);

                // (x)[i,j+1/2,k+1/2]
                // ( Ez[i,j+1,k+1/2] - Ez[i,j,k+1/2] ) / dy
                double val = 1.0 / dy;
                addTriplet(trips, domain, vindx, vind(i, jp, k, 2), val);
                addTriplet(trips, domain, vindx, vind(i, j, k, 2), -val);
                // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
                val = -1.0 / dz;
                addTriplet(trips, domain ,vindx, vind(i, j, kp, 1), val);
                addTriplet(trips, domain, vindx, vind(i, j, k, 1), -val);

                // (y)[i+1/2,j,k+1/2]
                // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
                val = 1.0 / dz;
                addTriplet(trips, domain, vindy, vind(i, j, kp, 0), val);
                addTriplet(trips, domain, vindy, vind(i, j, k, 0), -val);
                // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
                val = -1.0 / dx;
                addTriplet(trips, domain, vindy, vind(ip, j, k, 2), val);
                addTriplet(trips, domain, vindy, vind(i, j, k, 2), -val);

                // (z)[i+1/2,j+1/2,k]
                // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
                val = 1.0 / dx;
                addTriplet(trips, domain, vindz, vind(ip, j, k, 1), val);
                addTriplet(trips, domain, vindz, vind(i, j, k, 1), -val);
                // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
                val = -1.0 / dy;
                addTriplet(trips, domain, vindz, vind(i, jp, k, 0), val);
                addTriplet(trips, domain, vindz, vind(i, j, k, 0), -val);
            }
        }
    }
}


void Mesh::stencil_curlE_openZ(Operator &mat, const Domain &domain) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 12;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z() - 2; k++) {
                const int ip = (i != size.x() - 1) ? i + 1 : 3;
                const int jp = (j != size.y() - 1) ? j + 1 : 3;
                const int kp = k + 1;

                const int vindx = vind(i, j, k, 0);
                const int vindy = vind(i, j, k, 1);
                const int vindz = vind(i, j, k, 2);

                // (x)[i,j+1/2,k+1/2]
                // ( Ez[i,j+1,k+1/2] - Ez[i,j,k+1/2] ) / dy
                double val = 1.0 / dy;
                if (k > 0 && k < size.z() - 2) {
                    trips.push_back(Trip(vindx, vind(i, jp, k, 2), val));
                    trips.push_back(Trip(vindx, vind(i, j, k, 2), -val));
                }
                // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
                val = -1.0 / dz;
                if (kp > 1 && kp < size.z() - 2) {
                    trips.push_back(Trip(vindx, vind(i, j, kp, 1), val));
                }
                if (k > 1 && k < size.z() - 2) {
                    trips.push_back(Trip(vindx, vind(i, j, k, 1), -val));
                }
                // (y)[i+1/2,j,k+1/2]
                // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
                val = 1.0 / dz;
                if (kp > 1 && kp < size.z() - 2) {
                    trips.push_back(Trip(vindy, vind(i, j, kp, 0), val));
                }
                if (k > 1 && k < size.z() - 2) {
                    trips.push_back(Trip(vindy, vind(i, j, k, 0), -val));
                }
                // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
                val = -1.0 / dx;
                if (k > 0 && k < size.z() - 2) {
                    trips.push_back(Trip(vindy, vind(ip, j, k, 2), val));
                    trips.push_back(Trip(vindy, vind(i, j, k, 2), -val));
                }
                // (z)[i+1/2,j+1/2,k]
                // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
                if (k > 1 && k < size.z() - 2) {
                    val = 1.0 / dx;
                    trips.push_back(Trip(vindz, vind(ip, j, k, 1), val));
                    trips.push_back(Trip(vindz, vind(i, j, k, 1), -val));
                    // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
                    val = -1.0 / dy;
                    trips.push_back(Trip(vindz, vind(i, jp, k, 0), val));
                    trips.push_back(Trip(vindz, vind(i, j, k, 0), -val));
                }
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_divE(Operator &mat, const Domain &domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x()*size.y()*size.z()*6;
    trips.reserve(totalSize);
    double dx = domain.cell_size().x();
    double dy = domain.cell_size().y();
    double dz = domain.cell_size().z();
    for(int i = 0; i < size.x(); i++){
      for(int j = 0; j < size.y(); j++){
        for(int k = 0; k < size.z(); k++){
            const int im = (i != 0) ? i - 1 : size.x() - 4;
            const int jm = (j != 0) ? j - 1 : size.y() - 4;
            const int km = (k != 0) ? k - 1 : size.z() - 4;

            const int sindx = sind(i, j, k);

            // [i,j,k]
            // ( Ex[i+1/2,j,k] - Ex[i-1,j,k] ) / dx
            double val = 1.0 / dx;
            trips.push_back(Trip(sindx, vind(i, j, k, 0), val));
            trips.push_back(Trip(sindx, vind(im, j, k, 0), -val));
            // ( Ex[i,j+1/2,k] - Ex[i,j-1/2,k] ) / dy
            val = 1.0 / dy;
            trips.push_back(Trip(sindx, vind(i, j, k, 1), val));
            trips.push_back(Trip(sindx, vind(i, jm, k, 1), -val));
            // ( Ez[i,j,k+1/2] - Ez[i,j,k-1/2] ) / dz
            val = 1.0 / dz;
            trips.push_back(Trip(sindx, vind(i, j, k, 2), val));
            trips.push_back(Trip(sindx, vind(i, j, km, 2), -val));

        }
      }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}
