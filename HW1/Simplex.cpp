#include <algorithm>
#include <iostream>
#include <numeric>
#include <ranges>
#include <tuple>
#include <valarray>
#include <vector>
#include <iomanip>

namespace ranges = std::ranges;
using namespace std::string_literals;

void printVec(const auto& slice, std::string sep = " "s) {
    for (const auto v : slice) {
        std::cout << v << sep;
    }
    std::cout << '\n';
}

class Matrix {
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::valarray<double> data;

   public:
    explicit Matrix(std::size_t rows, std::size_t cols, double init = 0.0)
        : rows{rows},
          cols{cols},
          data{std::valarray<double>(init, rows * cols)} {}

    explicit Matrix(std::size_t rows, std::size_t cols,
                    const std::vector<std::vector<double>>& d)
        : rows{rows}, cols{cols}, data{std::valarray<double>(rows * cols)} {
        for (std::size_t i = 0; i != rows; ++i) {
            for (std::size_t j = 0; j != cols; ++j) {
                (*this)(i, j) = d[i][j];
            }
        }
    }

    std::size_t getRows() const { return rows; }

    std::size_t getCols() const { return cols; }

    double& operator()(std::size_t i, std::size_t j) {
        return data[i * cols + j];
    }

    std::slice_array<double> sliceRow(std::size_t row) {
        return data[std::slice(row * cols, cols, 1)];
    }

    std::valarray<double> row(std::size_t row) { return sliceRow(row); }

    std::slice_array<double> sliceCol(std::size_t col) {
        return data[std::slice(col, rows, cols)];
    }

    std::valarray<double> col(std::size_t col) { return sliceCol(col); }

    void print() {
        for (auto i = 0U; i != rows; ++i) {
            printVec(row(i), " \t"s);
        }
        std::cout << '\n';
    }
};

auto isNegative = [](const auto value) { return value < 0; };
auto isPositive = [](const auto value) { return value > 0; };

std::vector<std::size_t> getIndicesIfPred(const auto& vec,
                                          const std::size_t start,
                                          const std::size_t end, auto pred) {
    auto ret = std::vector<std::size_t>{};
    for (std::size_t i = start; i != end; ++i) {
        if (pred(vec[i])) {
            ret.push_back(i);
        }
    }
    return ret;
}

std::vector<std::size_t> initIndices(std::size_t n, int start = 0) {
    auto ret = std::vector<std::size_t>(n);
    std::iota(std::begin(ret), std::end(ret), start);
    return ret;
}

Matrix initData(std::size_t example = 0) {
    // we assume that this is a m x n matrix

    if (example == 1) {
        auto init = std::vector<std::vector<double>>{
            {-2, 1, 3, 1, 0, 0, 4},
            {2, 3, -1, 0, 1, 0, 3},
            {0, 1, 2, 0, 0, 1, 3},
            {1, -1, 0, 0, 0, 0, 0},
        };

        const auto m = init.size();
        const auto n = init[0].size();

        auto ret = Matrix{m, n};
        for (std::size_t i = 0; i != m; ++i) {
            for (std::size_t j = 0; j != n; ++j) {
                ret(i, j) = init[i][j];
            }
        }
        return ret;
    }
    if (example == 2) {
        auto init = std::vector<std::vector<double>>{
            {-1, 1, 1, 0, 2},
            {-2, 1, 0, 1, 1},
            {-1, -2, 0, 0, 0},
        };
        const auto m = init.size();
        const auto n = init[0].size();

        auto ret = Matrix{m, n};
        for (std::size_t i = 0; i != m; ++i) {
            for (std::size_t j = 0; j != n; ++j) {
                ret(i, j) = init[i][j];
            }
        }
        return ret;
    }
    if (example == 3) {
        auto init = std::vector<std::vector<double>>{
            {4, 5, -2, 1, 0, 0, 22},
            {1, -2, 1, 0, 1, 0, 30},
            {1, 1, 4, 0, 0, 1, 24},
            {3, -2, -4, 0, 0, 0, 0},
        };

        const auto m = init.size();
        const auto n = init[0].size();

        auto ret = Matrix{m, n};
        for (std::size_t i = 0; i != m; ++i) {
            for (std::size_t j = 0; j != n; ++j) {
                ret(i, j) = init[i][j];
            }
        }
        return ret;
    }
    auto init = std::vector<std::vector<double>>{
        {2, 1, 1, 0, 0, 100},
        {1, 1, 0, 1, 0, 80},
        {1, 0, 0, 0, 1, 40},
        {-1, -9, 0, 0, 0, 0},
    };

    const auto m = init.size();
    const auto n = init[0].size();

    auto ret = Matrix{m, n};
    for (std::size_t i = 0; i != m; ++i) {
        for (std::size_t j = 0; j != n; ++j) {
            ret(i, j) = init[i][j];
        }
    }
    return ret;
}

void simplexAlgorithm(std::size_t example = 0) {
    auto matrix = initData(example);
    auto initialMatrix = matrix;
    const auto m = matrix.getRows();
    const auto n = matrix.getCols();

    const auto columnIndices = initIndices(n - 1);
    auto rowIndices = initIndices(m - 1, n - m);

    std::cout << "Row indices: ";
    printVec(rowIndices);
    matrix.print();

    auto compareIndices = [&columnIndices](const auto i, const auto j) {
        return columnIndices[i] < columnIndices[j];
    };

    // n - 1 we don't need last column
    for (auto candidatePivots = getIndicesIfPred(matrix.row(m - 1), 0, n - 1, isNegative);
         not candidatePivots.empty();
         candidatePivots = getIndicesIfPred(matrix.row(m - 1), 0, n - 1, isNegative)) {
        const auto pivot =
            *ranges::min_element(candidatePivots, compareIndices);
        std::cout << "Pivot: " << pivot << '\n';
        const auto pivotCol = matrix.col(pivot);

        if (ranges::all_of(pivotCol, isNegative)) {
            printVec(pivotCol);
            matrix.print();
            std::cout << "Problem is unbounded\n";
            return;
        }

        // m - 1 because we don't need last row
        const auto positiveIndices = getIndicesIfPred(pivotCol, 0, m - 1, isPositive);

        auto division = [&](const auto index) {
            // the last element on row is divided by the pivot element on the
            // row
            return matrix(index, n - 1) / matrix(index, pivot);
        };

        auto compare = [&](const auto i, const auto j) {
            // the leaving element must have the minimum division
            // if more than one possible min element, the one who is
            // lexicographically first is chosen
            const auto iVal = division(i);
            const auto jVal = division(j);
            if (iVal == jVal) {
                return rowIndices[i] < rowIndices[j];
            }
            return iVal < jVal;
        };

        const auto leaving = *ranges::min_element(positiveIndices, compare);
        std::cout << "Leaving: " << leaving << '\n';
        const auto tkl = matrix(leaving, pivot);
        std::cout << "tkl: " << tkl << '\n';

        for (std::size_t i = 0; i != m; ++i) {
            if (i == leaving) {
                continue;
            }
            for (std::size_t j = 0; j != n; ++j) {
                if (j == pivot) {
                    continue;
                }
                
                matrix(i, j) = (matrix(i, j) * tkl -
                                matrix(i, pivot) * matrix(leaving, j)) /
                               tkl;
            }
        }

        matrix.sliceCol(pivot) = 0;    // set all column to 0
        matrix(leaving, pivot) = tkl;  // set this value back
        matrix.sliceRow(leaving) = matrix.row(leaving) / tkl;

        // updating indices
        rowIndices[leaving] = pivot;
        std::cout << "Row indices: ";
        printVec(rowIndices);
        matrix.print();
    }
    
    std::cout << "Optimum value: " << -matrix(m - 1, n - 1) << '\n';

    auto variables = std::vector<double>(n - 1);
    for (std::size_t i = 0; i != m - 1; ++i) {
        variables[rowIndices[i]] = matrix(i, n - 1);
    }
    std::cout << "Variables:\n";
    printVec(variables, " \t");
    for (std::size_t i = 0; i != m; ++i) {
        auto row = initialMatrix.row(i);
        printVec(row);
        std::cout << "Inner product: " << std::inner_product(std::begin(variables), std::end(variables), std::begin(row), 0.0) << '\n';
        std::cout << '\n';
    }

}

int main() {
    simplexAlgorithm(1);
    simplexAlgorithm(2);
    simplexAlgorithm(3);

    return 0;
}
