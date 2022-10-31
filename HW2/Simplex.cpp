#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <ranges>
#include <sstream>
#include <tuple>
#include <valarray>
#include <vector>

namespace ranges = std::ranges;
using namespace std::string_literals;

void printVec(auto& stream, const auto& slice, std::string sep = " "s) {
    for (const auto v : slice) {
        stream << v << sep;
    }
    stream << '\n';
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

    void addColumn(std::size_t col, double value)
    {
        // complicated code to add a new column. Also, very ineficient.
        // Faster solution would be to create a new vallaray,
        // copy each row until col inside it, add a 0.0, copy the rest of the row.
        // And each row can be copied in parallel since we know the start indices.
        auto newData = std::valarray<double>(value, rows * (cols + 1));
        for (std::size_t i = 0; i != rows; ++i) {
            for (std::size_t j = 0; j != cols + 1; ++j) {
                if (j == col) {
                    continue;
                }
                if (j > col) {
                    newData[i * (cols + 1) + j] = data[i * cols + j - 1];
                } else {
                    newData[i * (cols + 1) + j] = data[i * cols + j];
                }
            }
        }
        ++cols;
        data = newData;
    }

    void print(auto& stream) {
        for (auto i = 0U; i != rows; ++i) {
            printVec(stream, row(i), " \t"s);
        }
        std::cout << '\n';
    }
};

auto isNegative = [](const auto value) { return value < 0; };
auto isNegativeOrZero = [](const auto value) { return value <= 0; };
auto isPositive = [](const auto value) { return value > 0; };

using Indices = std::vector<std::size_t>;

Indices getIndicesIfPred(const auto& vec, const std::size_t start,
                         const std::size_t end, auto pred) {
    auto ret = Indices{};
    for (std::size_t i = start; i != end; ++i) {
        if (pred(vec[i])) {
            ret.push_back(i);
        }
    }
    return ret;
}

Indices initIndices(std::size_t n, int start = 0) {
    auto ret = Indices(n);
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
    if (example == 4) {
        auto init = std::vector<std::vector<double>>{
            {3, 2, 0, 0, 14},
            {2, -4, -1, 0, 2},
            {4, 3, 0, 1, 0, 19},
            {2, 3, 0, 0, 0, 0},
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
    if (example == 5) {
        auto init = std::vector<std::vector<double>>{
            {2, 1, -1, 0, 7},
            {2, 2, 0, 1, 10},
            {-1, 2, 0, 0, 0},
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

void printSolution(auto& stream, Matrix& matrix, Matrix& initialMatrix,
                   const Indices& rowIndices) {
    const auto m = matrix.getRows();
    const auto n = matrix.getCols();
    stream << "Optimum value: " << -matrix(m - 1, n - 1) << '\n';

    auto variables = std::vector<double>(n - 1);
    for (std::size_t i = 0; i != m - 1; ++i) {
        variables[rowIndices[i]] = matrix(i, n - 1);
    }
    stream << "Variables:\n";
    printVec(stream, variables, " \t");
    for (std::size_t i = 0; i != m; ++i) {
        auto row = initialMatrix.row(i);
        printVec(stream, row);
        stream << "Inner product: "
               << std::inner_product(std::begin(variables), std::end(variables),
                                     std::begin(row), 0.0)
               << '\n';
        stream << '\n';
    }
}

std::pair<Matrix, Indices> simplexAlgorithm(const Matrix& _matrix,
                                            const Indices& _rowIndices,
                                            const bool debug) {
    auto matrix = _matrix;
    const auto m = matrix.getRows();
    const auto n = matrix.getCols();
    auto rowIndices = _rowIndices;
    auto debugStream = std::ostringstream{};
    debugStream << "Row indices: ";
    printVec(debugStream, rowIndices);
    matrix.print(debugStream);

    for (auto candidatePivots =
             getIndicesIfPred(matrix.row(m - 1), 0, n - 1, isNegative);
         not candidatePivots.empty();
         candidatePivots =
             getIndicesIfPred(matrix.row(m - 1), 0, n - 1, isNegative)) {

        const auto pivot = *ranges::min_element(candidatePivots);

        const auto pivotCol = matrix.col(pivot);

        debugStream << "Pivot: " << pivot << '\n';

        if (ranges::all_of(pivotCol, isNegativeOrZero)) {
            printVec(debugStream, pivotCol);
            matrix.print(debugStream);
            std::cout << "Problem is unbounded\n";
            throw std::runtime_error{"Problem is unbounded"};
        }

        // m - 1 because we don't need last row
        const auto positiveIndices =
            getIndicesIfPred(pivotCol, 0, m - 1, isPositive);

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
                return rowIndices[iVal] < rowIndices[jVal];
            }
            return iVal < jVal;
        };

        const auto leaving = *ranges::min_element(positiveIndices, compare);
        const auto tkl = matrix(leaving, pivot);

        debugStream << "Leaving: " << leaving << '\n';
        debugStream << "tkl: " << tkl << '\n';

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
        debugStream << "Row indices: ";
        printVec(debugStream, rowIndices);
        matrix.print(debugStream);
    }

    if (debug) {
        std::cout << debugStream.str();
    }

    return {matrix, rowIndices};
}

Matrix prepareForPhaseOne(const Matrix& _matrix)
{
    auto matrix = _matrix;
    auto artificialSlackVariablesNumber = matrix.getRows() - 1;
    matrix.sliceRow(artificialSlackVariablesNumber) = 0; // setting the last row to 0
    for (std::size_t i = 0; i != artificialSlackVariablesNumber; ++i) {
        const auto columnIndex = matrix.getCols() - 1; // adding column left to RHS
        matrix.addColumn(columnIndex, 0); 
        matrix.sliceRow(artificialSlackVariablesNumber) = matrix.row(artificialSlackVariablesNumber) - matrix.row(i);
        matrix(i, columnIndex) = 1; // setting the artificial slack variable to 1
    }
    return matrix;
}

void doSimplexAlgorithm(std::size_t example = 0, const bool debug = false) {
    auto matrix = initData(example);
    if (debug) {
        matrix.print(std::cout);
    }
   
    try {
        auto phaseOneMatrix = prepareForPhaseOne(matrix);
        const auto m = phaseOneMatrix.getRows();
        const auto n = phaseOneMatrix.getCols();
        auto rowIndices = initIndices(m - 1, n - m);

        auto [newMatrix, newRowIndices] =
            simplexAlgorithm(phaseOneMatrix, rowIndices, debug);

        if (debug) {
            printSolution(std::cout, newMatrix, matrix, newRowIndices);
        }
        auto [newMatrix2, newRowIndices2] =
            simplexAlgorithm(newMatrix, newRowIndices, debug);
    } catch (...) {
    }
}

int main() {
    // doSimplexAlgorithm(1, true);
    // doSimplexAlgorithm(2, true);
    // doSimplexAlgorithm(3, true);
    doSimplexAlgorithm(4, true);
    doSimplexAlgorithm(5, true);

    return 0;
}
