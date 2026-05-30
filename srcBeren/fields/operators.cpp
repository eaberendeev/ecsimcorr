#include "operators.h"

#include "Mesh.h"
#include "Shape.h"
#include "World.h"
#include "config.h"
#include "pmms.hpp"
#include "timer.h"
#include "util.h"
// matrix ColMajor
// (0,0) (0,1) (0,2)
// (1,0) (1,1) (1,2)
// (2,0) (2,1) (2,2)

// Ex(i+/2,j,k), Ey(i,j+1/2,k), Ez(i,j,k+1/2)
// Bx(i,j+1/2,k+1/2), By(i+1/2,j,k+1/2), Bz(i+/2,j+1/2,k)

void Mesh::stencil_Imat(Operator& mat, const Domain& domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z();
    trips.reserve(totalSize);

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                // i,j,k
                trips.push_back(Trip(vind(i, j, k, 0), vind(i, j, k, 0), 1.0));

                // i,j,k
                trips.push_back(Trip(vind(i, j, k, 1), vind(i, j, k, 1), 1.0));

                // i,j,k
                trips.push_back(Trip(vind(i, j, k, 2), vind(i, j, k, 2), 1.0));
            }
        }
    }
    std::cout << size << " " << 3 * totalSize << " " << Imat.rows() << " " << Imat.cols() << "\n";
    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_Lmat(Operator& mat, const Domain& domain) {
    //     std::vector<Trip> trips;
    //     const auto size = domain.size();
    //     const int rowsCount = 3 * size.x() * size.y() * size.z();
    //     const size_t totalSize = (size_t)rowsCount *
    //     LMAT_MAX_ELEMENTS_PER_ROW; std::cout << totalSize << "\n";
    //     trips.reserve(totalSize);

    // #pragma omp parallel
    //     {
    //         std::vector<Trip> local_trips;
    //         local_trips.reserve(totalSize / omp_get_max_threads());

    // #pragma omp for schedule(dynamic, 8)
    //         for (int row = 0; row < rowsCount; row++) {
    //             for (const auto &[col, value] : LmatX[row]) {
    //                 if (std::abs(value) > LMAT_VALUE_TOLERANCE)
    //                     local_trips.emplace_back(row, col, value);
    //             }
    //         }

    // #pragma omp critical
    //         trips.insert(trips.end(), local_trips.begin(),
    //         local_trips.end());
    //     }
    //     std::cout << "trips size: " << trips.size() << std::endl;

    //     mat.setFromTriplets(trips.begin(), trips.end());
}

static bool equalVecsTriplets(const std::vector<Triplet>& a, const std::vector<Triplet>& b) {
    if (a.size() != b.size()) {
        std::cout << "size a " << a.size() << " size b " << b.size() << "\n";
        return false;
    }
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i].col() != b[i].col() || a[i].row() != b[i].row() || fabs(a[i].value() - b[i].value()) > 1.e-15) {
            std::cout << "row " << a[i].row() << " col " << a[i].col() << " value " << a[i].value() << " row "
                      << b[i].row() << " col " << b[i].col() << " value " << b[i].value() << " "
                      << std::abs(a[i].value() - b[i].value()) << "\n";
            return false;
        }
    return true;
}

std::vector<Triplet> multyPhaseMerge(std::vector<std::vector<Triplet>>& local_vectors) {
    // Собираем только непустые локальные векторы для дальнейшего слияния
    std::vector<std::vector<Triplet>> non_empty;
    for (auto& v : local_vectors) {
        if (!v.empty()) {
            non_empty.push_back(std::move(v));
        }
    }
    if (non_empty.empty()) {
        return {};
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
    return std::move(non_empty[0]);
}

void Mesh::stencil_Lmat2(Operator& mat, const Domain& domain) {
    RECORD_TIMER;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;
    static int check_count = 0;
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
    size_t total_iterations = static_cast<size_t>(max_i - BORDER) * (max_j - BORDER) * (max_k - BORDER);
    size_t estimated_per_thread = (total_iterations / num_threads) * 9;

    // Локальные векторы для каждого потока с предварительным резервированием
    // памяти.
    std::vector<std::vector<Triplet>> local_vectors(num_threads);

    timer::timer timer1("parallel sections");

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
        timer::flatTimer timerLocal("main OMP section");

        int tid = omp_get_thread_num();

#pragma omp for collapse(3) schedule(dynamic, 8) nowait
        for (int i = BORDER; i < max_i; ++i)
            for (int j = BORDER; j < max_j; ++j)
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX2.non_zeros[sind(i, j, k)])
                        continue;
                    const Block& block = LmatX2[sind(i, j, k)];

                    // X component
                    processComponent<XIndexer, XIndexer, 0>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<XIndexer, YIndexer, 1>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<XIndexer, ZIndexer, 2>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);

                    // Y component
                    processComponent<YIndexer, XIndexer, 3>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<YIndexer, YIndexer, 4>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<YIndexer, ZIndexer, 5>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);

                    // Z component
                    processComponent<ZIndexer, XIndexer, 6>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<ZIndexer, YIndexer, 7>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<ZIndexer, ZIndexer, 8>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
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
                if (vec[index].row() == vec[j].row() && vec[index].col() == vec[j].col()) {
                    vec[index].value() += vec[j].value();
                } else {
                    ++index;
                    vec[index] = vec[j];
                }
            }
            vec.resize(index + 1);
        }
    }

    timer1.finish();
    timer::timer timer2("section 2");

    double time02 = omp_get_wtime();
    pmms::PMMSOptions opt;
    opt.useSampling = true;
    opt.usePWayMerge = true;
    opt.oversample = 4;
    std::vector<Triplet> test = pmms::parallelMultiwayMergeSort<Triplet>(std::move(local_vectors), num_threads, opt);
    if (!test.empty()) {
        size_t index = 0;
        for (size_t j = 1; j < test.size(); ++j) {
            if (test[index].row() == test[j].row() && test[index].col() == test[j].col()) {
                test[index].value() += test[j].value();
            } else {
                ++index;
                test[index] = test[j];
            }
        }
        test.resize(index + 1);
    }

    timer2.finish();
    timer::timer timer3("section 3");

    check_count++;
    // В non_empty[0] теперь находится глобальный вектор, уже
    // отсортированный и с устранёнными дубликатами.
    double time2 = omp_get_wtime();
    double time3;
    if (check_count < 20) {
        const std::vector<Triplet> trips = multyPhaseMerge(local_vectors);
        time3 = omp_get_wtime();
        if (equalVecsTriplets(test, trips))
            std::cout << "Equal!\n";
    } else {
        time3 = omp_get_wtime();
    }

    // std::vector<Trip> trips2;
    // trips2.reserve(trips.size());
    // for (long i = 0; i < trips.size(); i++){
    //     trips2.emplace_back(trips[i].row, trips[i].col, trips[i].value);
    // }

    // double time4 = omp_get_wtime();

    // mat.setFromTriplets(trips.begin(), trips.end());

    timer3.finish();

    optimizedSetFromTriplets(mat, test);
    double time5 = omp_get_wtime();
    std::cout << time2 - time1 << " " << time3 - time2 << " " << time5 - time3 << " " << time2 - time02 << "\n";
    std::cout << "Matrix L (block) was created." << " trips size: " << test.size() << std::endl;
}

// TODO implement parallelMultiwayMergeSort and merge to convert block to csr
// (it is dublicated)
void Mesh::stencil_Lmat2_NGP(Operator& mat, const Domain& domain) {
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
    size_t total_iterations = static_cast<size_t>(max_i - BORDER) * (max_j - BORDER) * (max_k - BORDER);
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
                    const auto& block = LmatX_NGP[sind(i, j, k)];

                    // X component
                    processComponent<XIndexerNGP, XIndexerNGP, 0>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<XIndexerNGP, YIndexerNGP, 1>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<XIndexerNGP, ZIndexerNGP, 2>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);

                    // Y component
                    processComponent<YIndexerNGP, XIndexerNGP, 3>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<YIndexerNGP, YIndexerNGP, 4>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<YIndexerNGP, ZIndexerNGP, 5>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);

                    // Z component
                    processComponent<ZIndexerNGP, XIndexerNGP, 6>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<ZIndexerNGP, YIndexerNGP, 7>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<ZIndexerNGP, ZIndexerNGP, 8>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                }

        // Сортировка локального вектора по (row, col)
        auto& vec = local_vectors[tid];
        std::sort(vec.begin(), vec.end(), compareTriplets);

        // Устранение дубликатов: проход по отсортированному вектору, складываем
        // значения для одинаковых ключей
        if (!vec.empty()) {
            size_t index = 0;
            for (size_t j = 1; j < vec.size(); ++j) {
                if (vec[index].row() == vec[j].row() && vec[index].col() == vec[j].col()) {
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
    const std::vector<Triplet>& trips = non_empty.empty() ? std::vector<Triplet>() : non_empty.front();

    // std::vector<Trip> trips2;
    // trips2.reserve(trips.size());
    // for (long i = 0; i < trips.size(); i++){
    //     trips2.emplace_back(trips[i].row, trips[i].col, trips[i].value);
    // }

    // double time4 = omp_get_wtime();

    // mat.setFromTriplets(trips.begin(), trips.end());
    optimizedSetFromTriplets(mat, trips);
    double time5 = omp_get_wtime();
    std::cout << time2 - time1 << " " << time3 - time2 << " " << time5 - time3 << "\n";
    std::cout << "Matrix L (block) was created." << " trips size: " << trips.size() << std::endl;
}

void Mesh::stencil_curlB(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 12;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    bool is_periodic[3] = {bc_handler.is_periodic(0), bc_handler.is_periodic(1), bc_handler.is_periodic(2)};

    auto addRowIfInside = [&](int i, int j, int k, int comp) -> bool {
        return domain.is_inside_node_periodic(i, j, k, FieldType::ELECTRIC, comp, is_periodic);
    };

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int im = bc_handler.wrap_index(i - 1, 0, domain);
                const int jm = bc_handler.wrap_index(j - 1, 1, domain);
                const int km = bc_handler.wrap_index(k - 1, 2, domain);

                const int vindx = domain.vind(i, j, k, 0);   // Ex
                const int vindy = domain.vind(i, j, k, 1);   // Ey
                const int vindz = domain.vind(i, j, k, 2);   // Ez

                // (x)[i+1/2,j,k]
                if (addRowIfInside(i, j, k, 0)) {
                    // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
                    double val = 1.0 / dy;
                    trips.emplace_back(vindx, domain.vind(i, j, k, 2), val);
                    trips.emplace_back(vindx, domain.vind(i, jm, k, 2), -val);
                    // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
                    val = -1.0 / dz;
                    trips.emplace_back(vindx, domain.vind(i, j, k, 1), val);
                    trips.emplace_back(vindx, domain.vind(i, j, km, 1), -val);
                }
                // (y)[i,j+1/2,k]
                if (addRowIfInside(i, j, k, 1)) {
                    // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
                    double val = 1.0 / dz;
                    trips.emplace_back(vindy, domain.vind(i, j, k, 0), val);
                    trips.emplace_back(vindy, domain.vind(i, j, km, 0), -val);
                    // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
                    val = -1.0 / dx;
                    trips.emplace_back(vindy, domain.vind(i, j, k, 2), val);
                    trips.emplace_back(vindy, domain.vind(im, j, k, 2), -val);
                }
                // (z)[i,j,k+1/2]
                if (addRowIfInside(i, j, k, 2)) {
                    // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                    double val = 1.0 / dx;
                    trips.emplace_back(vindz, domain.vind(i, j, k, 1), val);
                    trips.emplace_back(vindz, domain.vind(im, j, k, 1), -val);
                    // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                    val = -1.0 / dy;
                    trips.emplace_back(vindz, domain.vind(i, j, k, 0), val);
                    trips.emplace_back(vindz, domain.vind(i, jm, k, 0), -val);
                }
                //  bc_handler.modify_curlB_stencil(i, j, k, trips, domain);
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_curlE(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 12;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    bool is_periodic[3] = {bc_handler.is_periodic(0), bc_handler.is_periodic(1), bc_handler.is_periodic(2)};

    // Проверяет, принадлежит ли магнитный узел области (с учётом периодичности)
    auto addRowIfInside = [&](int i, int j, int k, int comp) -> bool {
        return domain.is_inside_node_periodic(i, j, k, FieldType::MAGNETIC, comp, is_periodic);
    };

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int ip = bc_handler.wrap_index(i + 1, 0, domain);
                const int jp = bc_handler.wrap_index(j + 1, 1, domain);
                const int kp = bc_handler.wrap_index(k + 1, 2, domain);

                const int vindx = domain.vind(i, j, k, 0);   // Bx
                const int vindy = domain.vind(i, j, k, 1);   // By
                const int vindz = domain.vind(i, j, k, 2);   // Bz

                // (x)[i,j+1/2,k+1/2]
                if (addRowIfInside(i, j, k, 0)) {
                    // ( Ez[i,j+1,k+1/2] - Ez[i,j,k+1/2] ) / dy
                    double val = 1.0 / dy;
                    trips.emplace_back(vindx, domain.vind(i, jp, k, 2), val);
                    trips.emplace_back(vindx, domain.vind(i, j, k, 2), -val);
                    // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
                    val = -1.0 / dz;
                    trips.emplace_back(vindx, domain.vind(i, j, kp, 1), val);
                    trips.emplace_back(vindx, domain.vind(i, j, k, 1), -val);
                }

                // (y)[i+1/2,j,k+1/2]
                if (addRowIfInside(i, j, k, 1)) {
                    // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
                    double val = 1.0 / dz;
                    trips.emplace_back(vindy, domain.vind(i, j, kp, 0), val);
                    trips.emplace_back(vindy, domain.vind(i, j, k, 0), -val);
                    // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
                    val = -1.0 / dx;
                    trips.emplace_back(vindy, domain.vind(ip, j, k, 2), val);
                    trips.emplace_back(vindy, domain.vind(i, j, k, 2), -val);
                }

                // (z)[i+1/2,j+1/2,k]
                if (addRowIfInside(i, j, k, 2)) {
                    // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
                    double val = 1.0 / dx;
                    trips.emplace_back(vindz, domain.vind(ip, j, k, 1), val);
                    trips.emplace_back(vindz, domain.vind(i, j, k, 1), -val);
                    // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
                    val = -1.0 / dy;
                    trips.emplace_back(vindz, domain.vind(i, jp, k, 0), val);
                    trips.emplace_back(vindz, domain.vind(i, j, k, 0), -val);
                }

                //  bc_handler.modify_curlE_stencil(i, j, k, trips, domain);
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_divE(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 6;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    bool is_periodic[3] = {bc_handler.is_periodic(0), bc_handler.is_periodic(1), bc_handler.is_periodic(2)};

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                if (!domain.is_inside_node_periodic(i, j, k, FieldType::DENSITY, 0, is_periodic))
                    continue;

                const int im = bc_handler.wrap_index(i - 1, 0, domain);
                const int jm = bc_handler.wrap_index(j - 1, 1, domain);
                const int km = bc_handler.wrap_index(k - 1, 2, domain);

                const int sindx = domain.sind(i, j, k);

                // [i,j,k]
                // ( Ex[i+1/2,j,k] - Ex[i-1/2,j,k] ) / dx
                double val = 1.0 / dx;
                trips.push_back(Trip(sindx, domain.vind(i, j, k, 0), val));
                trips.push_back(Trip(sindx, domain.vind(im, j, k, 0), -val));
                // ( Ey[i,j+1/2,k] - Ey[i,j-1/2,k] ) / dy
                val = 1.0 / dy;
                trips.push_back(Trip(sindx, domain.vind(i, j, k, 1), val));
                trips.push_back(Trip(sindx, domain.vind(i, jm, k, 1), -val));
                // ( Ez[i,j,k+1/2] - Ez[i,j,k-1/2] ) / dz
                val = 1.0 / dz;
                trips.push_back(Trip(sindx, domain.vind(i, j, k, 2), val));
                trips.push_back(Trip(sindx, domain.vind(i, j, km, 2), -val));
                // TODO: add bc_handler modify stencil
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}
void Mesh::stencil_smooth_1d(Operator& mat, const Domain& domain, int dim) {
    constexpr int COMPONENTS = 3;
    constexpr int STENCIL = 3;
    const auto size = domain.size();
    const int nx = size.x();
    const int ny = size.y();
    const int nz = size.z();

    const int totalCells = nx * ny * nz;
    // checker
    mat.resize(domain.total_size() * 3, domain.total_size() * 3);

    if (dim != Axis::X || dim != Axis::Y || dim != Axis::Z) {
        std::cout << "Error: Invalid dimension." << std::endl;
        return;
    }
    std::vector<Trip> trips;
    trips.reserve(static_cast<size_t>(totalCells * COMPONENTS * STENCIL));

    auto addTriplet = [](std::vector<Trip>& trips, const Domain& domain, int row, int col, double val) {
        bool onArea = domain.is_inside_node(row, FieldType::ELECTRIC);
        if (onArea) {
            trips.push_back(Trip(row, col, val));
        }
    };

    // auto imod = [](int a, int m) {
    //     // % m - 1
    //     if (a > m - 1) return 3;
    //     if (a < 0) return m - 4;
    //     return a;
    // };

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                int im = i - 1, ip = i + 1;
                int jm = j - 1, jp = j + 1;
                int km = k - 1, kp = k + 1;

                // if (dim == Axis::X && domain.is_periodic_bound(X)) {
                //     im = imod(i - 1, nx);
                //     ip = imod(i + 1, nx);
                // } else if (dim == Axis::Y && domain.is_periodic_bound(Y)) {
                //     jm = imod(j - 1, ny);
                //     jp = imod(j + 1, ny);
                // } else if (dim == Axis::Z &&
                //            domain.is_periodic_bound(Z)) {
                //     km = imod(k - 1, nz);
                //     kp = imod(k + 1, nz);
                // }

                for (int ax = 0; ax < COMPONENTS; ++ax) {
                    const int row = vind(i, j, k, ax);

                    int col_m, col_c, col_p;
                    if (dim == Axis::X) {
                        col_m = vind(im, j, k, ax);
                        col_c = vind(i, j, k, ax);
                        col_p = vind(ip, j, k, ax);
                    } else if (dim == Axis::Y) {
                        col_m = vind(i, jm, k, ax);
                        col_c = vind(i, j, k, ax);
                        col_p = vind(i, jp, k, ax);
                    } else {
                        col_m = vind(i, j, km, ax);
                        col_c = vind(i, j, k, ax);
                        col_p = vind(i, j, kp, ax);
                    }

                    addTriplet(trips, domain, row, col_m, 0.25);
                    addTriplet(trips, domain, row, col_c, 0.5);
                    addTriplet(trips, domain, row, col_p, 0.25);
                }
            }
        }
    }

    mat.setFromTriplets(trips.begin(), trips.end());
}
