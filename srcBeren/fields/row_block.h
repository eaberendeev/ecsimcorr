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

    // Both mergeFromOthers do the same: one with dynamic-size arrays, another with compile-time arrays
    template <int otherNnz>
    void mergeFromOthers(int count, const RowBlock<otherNnz>* others) {
        if (count == 1) {
            mergeFromOthers<1>(others);
        } else if (count == 2) {
            mergeFromOthers<2>(others);
        } else if (count == 3) {
            mergeFromOthers<3>(others);
        } else if (count == 4) {
            mergeFromOthers<4>(others);
        }

        if (count < 5) {
            return;
        }

        int smallestCols[count];

        int its[count];
        // more cache-friendly access
        int othersNnz[count];
        std::fill_n(its, count, 0);

        for (int i = 0; i < count; ++i) {
            othersNnz[i] = others[i].nnz;
        }

        for (int i = 1; i < count; ++i) {
            assert(others[i].row == others[0].row);
        }

        row = others[0].row;
        nnz = 0;

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
    }

    template <int count, int otherNnz>
    void mergeFromOthers(const RowBlock<otherNnz>* others) {
        if constexpr (count == 0) {
            return;
        } else if constexpr (count == 1) {
            row = others[0].row;
            nnz = others[0].nnz;
            for (int i = 0; i < nnz; ++i) {
                values[i] = others[0].values[i];
                columns[i] = others[0].columns[i];
            }

            return;
        }

        int smallestCols[count];
        int its[count]{0};
        // more cache-friendly access
        int othersNnz[count];
        for (int i = 0; i < count; ++i) {
            othersNnz[i] = others[i].nnz;
        }

        for (int i = 1; i < count; ++i) {
            assert(others[i].row == others[0].row);
        }

        row = others[0].row;
        nnz = 0;

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
