#pragma once

#include <array>

template <int maxNnz = 12 * 12 * 9>
struct RowBlock {
    RowBlock(int rowIn) : row(rowIn), nnz(0) {
    }
    RowBlock() {
    }

    ~RowBlock() {
    }

    template <int otherNnz>
    RowBlock(int count, const RowBlock<otherNnz>* others) {
        mergeFromOthers(count, others);
    }

    RowBlock& operator=(const RowBlock<maxNnz>& other) {
        nnz = other.nnz;
        row = other.row;
        assert(nnz <= maxNnz);

        std::copy_n(other.values.begin(), nnz, values.begin());
        std::copy_n(other.columns.begin(), nnz, columns.begin());

        return *this;
    }

    RowBlock(const RowBlock& other) {
        nnz = other.nnz;
        row = other.row;

        std::copy_n(other.values.begin(), nnz, values.begin());
        std::copy_n(other.columns.begin(), nnz, columns.begin());
    }

    bool operator!=(const RowBlock& other) {
        if (row != other.row || nnz != other.nnz) {
            return true;
        }
        for (int i = 0; i < nnz; ++i) {
            if (values[i] != other.values[i] || columns[i] != other.columns[i]) {
                return true;
            }
        }
        return false;
    }

    void push_back_value(int col, double val) {
        assert(nnz < maxNnz);
        values[nnz] = val;
        columns[nnz] = col;
        nnz += 1;
    }

    void sort() {
        std::array<std::pair<int, double>, maxNnz> pairs;
        for (int i = 0; i < nnz; ++i) {
            pairs[i] = {columns[i], values[i]};
        }

        std::sort(pairs.begin(), pairs.begin() + nnz,
                  [](const std::pair<int, double>& a, const std::pair<int, double>& b) { return a.first < b.first; });

        for (int i = 0; i < nnz; ++i) {
            columns[i] = pairs[i].first;
            values[i] = pairs[i].second;
        }
    }

    template <int otherNnz>
    void mergeFromOthers(int count, const RowBlock<otherNnz>* others) {
        static int counter = 0;
        counter += 1;

        // if (counter == 13275) {
        //     std::cout << "stack smash at start?" << std::endl;
        // }

        // std::cout << "counter is: " << counter << std::endl;

        int smallestCols[(maxNnz + otherNnz - 1) / otherNnz];

        // Eigen::Vector<int, -1> its(count);
        // its.fill(0);
        int its[count];
        // more cache-friendly access
        int othersNnz[count];
        std::fill_n(its, count, 0);

        // Eigen::Vector<int, -1> othersNnz(count);
        for (int i = 0; i < count; ++i) {
            othersNnz[i] = others[i].nnz;
        }

        for (int i = 1; i < count; ++i) {
            assert(others[i].row == others[0].row);
        }

        row = others[0].row;
        nnz = 0;

        // if (counter == 13275) {
        //     std::cout << "stack smash after?" << std::endl;
        // }

        while (true) {
            int smallestCol = std::numeric_limits<int>::max();
            int smallestColsCount = 0;
            double acc = 0;
            for (int i = 0; i < count; ++i) {
                assert(its[i] >= 0);
                if (its[i] < othersNnz[i]) {
                    const int otherCol = others[i].columns[its[i]];
                    if (otherCol == smallestCol) {
                        acc += others[i].values[its[i]];
                        smallestCols[smallestColsCount] = i;
                        smallestColsCount += 1;
                    } else if (otherCol < smallestCol) {
                        smallestCol = otherCol;
                        acc = others[i].values[its[i]];
                        smallestCols[0] = i;
                        smallestColsCount = 1;
                    }
                }
            }

            if (smallestCol == std::numeric_limits<int>::max()) {
                break;
            }

            for (int i = 0; i < smallestColsCount; ++i) {
                its[smallestCols[i]] += 1;
            }

            assert(nnz < maxNnz);
            columns[nnz] = smallestCol;
            values[nnz] = acc;
            nnz += 1;
        }

        // if (counter == 13275) {
        //     std::cout << "STOP HERE" << std::endl;
        //     std::cout << "nnz: " << nnz << " of " << maxNnz << std::endl;
        //     return;
        // }

        // for (int i = 0; i < count; ++i) {
        //     assert(its[i] >= 0);
        // }

        // std::cout << "END!" << std::endl;
    }

    friend std::ostream& operator<<(std::ostream& os, RowBlock block) {
        os << "Row block for row " << block.row << ", nnz: " << block.nnz << " of " << maxNnz << ", columns-values: ";
        for (int i = 0; i < block.nnz; ++i) {
            os << block.columns[i] << "-" << block.values[i] << ", ";
        }
        return os;
    }

    int row;
    int nnz;
    std::array<double, maxNnz> values;
    std::array<int, maxNnz> columns;
};

template <typename RowIdx, typename ColIdx, int DIR, typename Block_t>
static void blockToRowBlocks(int i_cell, int j_cell, int k_cell, const Block_t& block, [[maybe_unused]] int Nx, int Ny,
                             int Nz, double tolerance, std::vector<RowBlock<12>>& rowBlocks) {
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };
    int addedRows = 0;
    for (int x1 = 0; x1 < RowIdx::size_x; ++x1) {
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1) {
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1) {
                const int row = vind(i_cell + x1 + RowIdx::offset_x, j_cell + y1 + RowIdx::offset_y,
                                     k_cell + z1 + RowIdx::offset_z, RowIdx::dir);

                bool isAddedInThisRow = false;

                const int rowIdx = RowIdx::calculate(x1, y1, z1);

                for (int x2 = 0; x2 < ColIdx::size_x; ++x2) {
                    for (int y2 = 0; y2 < ColIdx::size_y; ++y2) {
                        for (int z2 = 0; z2 < ColIdx::size_z; ++z2) {
                            const int colIdx = ColIdx::calculate(x2, y2, z2);
                            const double val = block(rowIdx, colIdx, DIR);

                            if (std::abs(val) > tolerance) {
                                const int col = vind(i_cell + x2 + ColIdx::offset_x, j_cell + y2 + ColIdx::offset_y,
                                                     k_cell + z2 + ColIdx::offset_z, ColIdx::dir);
                                if (!isAddedInThisRow) [[unlikely]] {
                                    rowBlocks.emplace_back(row);
                                    isAddedInThisRow = true;
                                    addedRows += 1;
                                }
                                rowBlocks.back().push_back_value(col, val);
                            }
                        }
                    }
                }

                // if (row == 1891471 && isAddedInThisRow) {
                //     std::cout << "row " << addedRows << ": " << rowBlocks.back() << std::endl;
                // }
            }
        }
    }
}
