//
// Created by squir on 19.03.2023.
//

#ifndef IOP_MATRIX_H
#define IOP_MATRIX_H

#include <iostream>
#include <vector>
#include <ostream>
#include "Fraction.h"
#include <exception>
#include <iomanip>
#include <set>
#include <cmath>


//генерирует сочетания из n по k
void
combinations_(std::set<std::set<int>> &combinations, std::vector<int> &currentComb, int i, int b, int n, int k);


//генерирует сочетания из n по k
std::set<std::set<int>> combinations(int n, int k);

template<typename T>
requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
std::vector<T> &operator*(std::vector<T> &v, T value) {
    for (auto &x: v) {
        x *= value;
    }
    return v;
}

namespace matrices {
    //класс Матрица
    template<class T> requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    class Matrix : public std::vector<std::vector<T>> {
    public:
        using std::vector<std::vector<T>>::vector;

        //конструктор с входными данными количества строк nRows и столбцов nCols
        Matrix(int nRows, int nCols) : std::vector<std::vector<T>>(nRows, std::vector<T>(nCols)) {}

        //конструктор с входными данными количества строк nRows, столбцов nCols
        // и значением для инициализации всех элементов матрицы value
        Matrix(int nRows, int nCols, T value) : std::vector<std::vector<T>>(nRows, std::vector<T>(nCols, value)) {}

        //конструктор квадратной матрицы размерностью size
        explicit Matrix(int size) : Matrix(size, size) {}

        //возвращает количество строк матрицы
        [[nodiscard]]int nRows() const {
            return this->size();
        }

        //возвращает количество столбцов матрицы
        [[nodiscard]]int nCols() const {
            return this->at(0).size();
        }

        Matrix operator*(Matrix<T> m) {
            if (nCols() != m.nRows())
                throw std::domain_error
                        ("multiplication cannot be performed: matrices of unsuitable dimension");

            Matrix<T> mul(nRows(), m.nCols(), 0);
            for (int i = 0; i < nRows(); ++i)
                for (int j = 0; j < m.nCols(); ++j) {
                    for (int k = 0; k < nCols(); ++k)
                        mul[i][j] += (*this)[i][k] * m[k][j];
                }

            return mul;
        }

        Matrix operator*(T number) {
            Matrix<T> mul(nRows(), nCols());
            for (int j = 0; j < nCols(); ++j) {
                for (int k = 0; k < nCols(); ++k) {
                    mul[j][k] = (*this)[j][k] * number;
                }
            }
            return mul;
        }

        Matrix operator*(std::vector<T> vec) {
            if (nCols() != vec.size())
                throw std::domain_error
                        ("multiplication cannot be performed: matrix and vector of unsuitable dimension");

            Matrix<T> mul(nRows(), vec.size(), 0);
            for (int i = 0; i < nRows(); ++i) {
                for (int k = 0; k < nCols(); ++k) {
                    mul[i][k] += (*this)[i][k] * vec[k];
                }
            }

            return mul;
        }

        Matrix operator+(Matrix<T> m) {
            if (nCols() != m.nCols() || nRows() != m.nRows())
                throw std::domain_error
                        ("addition cannot be performed: matrices have different sizes");

            Matrix<T> sum = *this;
            for (int i = 0; i < nRows(); ++i) {
                for (int j = 0; j < nCols(); ++j) {
                    sum[i][j] += m[i][j];
                }
            }

            return sum;
        }

        Matrix operator-(Matrix<T> m) {
            if (nCols() != m.nCols() || nRows() != m.nRows())
                throw std::domain_error
                        ("subtraction cannot be performed: matrices have different sizes");

            Matrix<T> difference = *this;
            for (int i = 0; i < nRows(); ++i) {
                for (int j = 0; j < nCols(); ++j) {
                    difference[i][j] -= m[i][j];
                }
            }
            return difference;
        }

        //возвращает транспонированную матрицу
        Matrix transpose() {
            Matrix<T> MT(nCols(), nRows());
            for (int i = 0; i < nCols(); ++i) {
                for (int j = 0; j < nRows(); ++j) {
                    MT[i][j] = (*this)[j][i];
                }
            }

            return MT;
        }

        //возвращает ранг матрицы
        int rang() {
            Matrix<T> triangle = GaussJordanMethod();
            int rang = 0;
            for (int i = 0; i < nRows(); ++i) {
                for (int j = 0; j < nCols(); ++j) {
                    if (triangle[i][j] != 0) {
                        rang++;
                        break;
                    }
                }
            }
            return rang;
        }

        //возвращает матрицу, преобразованную прямым ходом метода Гаусса
        //в треугольную
        Matrix directСourseOfGaussMethod() {
            Matrix<T> triangleMatrix = *this;
            for (int i = 0; i < nRows() - 1; ++i) {
                int maxIndex = triangleMatrix.findMaxInPartOfCol(i);
                if (maxIndex != i)
                    std::swap(triangleMatrix[maxIndex], triangleMatrix[i]);
                if (triangleMatrix[i][i] != 0) {
                    for (int j = i + 1; j < nRows(); ++j) {
                        T multiplier = -triangleMatrix[j][i] / triangleMatrix[i][i];

                        for (int k = i; k < nCols(); ++k) {
                            triangleMatrix[j][k] += triangleMatrix[i][k] * multiplier;
                        }
                    }
                }
            }
            return triangleMatrix;
        }


        //возвращает матрицу, преобразованную прямым ходом метода Жордана-Гаусса
        Matrix GaussJordanMethod() {
            Matrix<T> triangleMatrix = *this;
            for (int i = 0; i < nRows(); ++i) {
                int maxIndex = triangleMatrix.findMaxInPartOfCol(i);
                if (maxIndex != i)
                    std::swap(triangleMatrix[maxIndex], triangleMatrix[i]);

                T maxEl = triangleMatrix[i][i];
                if (triangleMatrix[i][i] != 0) {
                    for (int j = i; j < nCols(); ++j) {
                        triangleMatrix[i][j] /= maxEl;
                    }

                    for (int j = 0; j < nRows(); ++j) {

                        if (i != j) {
                            T multiplier = triangleMatrix[j][i];

                            for (int k = i; k < nCols(); ++k) {
                                triangleMatrix[j][k] -= triangleMatrix[i][k] * multiplier;
                            }
                        }
                    }
                }
            }


            return triangleMatrix;
        }

        // умножает строку на число
        void mulRow(int indexRow, T el) {
            for (int i = 0; i < (*this)[indexRow].size(); ++i) {
                (*this)[indexRow][i] *= el;
            }
        }

        // возвращает матрицу всех базисных решений СЛАУ
        Matrix getAllBasicSolutions() {
            Matrix<T> matrix = *this;
            matrix.removeCol(nCols() - 1);

            int r1 = rang();
            int r2 = matrix.rang();
            if (r1 != r2)
                throw std::runtime_error("rang A != rang ~A");
            matrix = GaussJordanMethod();
            for (int i = 0; i <= matrix.size() - r1; ++i) {
                matrix.pop_back();
            }
            int nVars = matrix.nCols() - 1;
            std::set<int> allVars;
            for (int i = 1; i <= nVars; ++i) {
                allVars.insert(i);
            }

            Matrix<T> allBasicSolutions;
            auto combsOfBasicVars = combinations(nVars, r1);
            for (auto &x: combsOfBasicVars) {
                std::set<int> curCombOfFreeVars;
                std::set_difference(allVars.begin(), allVars.end(), x.begin(), x.end(),
                                    std::inserter(curCombOfFreeVars, curCombOfFreeVars.begin()));
                Matrix<T> basicMatrix = matrix;
                int countDeleted = 0;
                for (auto &freeVar: curCombOfFreeVars) {
                    basicMatrix.removeCol(freeVar - 1 - countDeleted);
                    countDeleted++;
                }
                basicMatrix = basicMatrix.GaussJordanMethod();
                auto valuesOfBasicVars = basicMatrix.reverseСourseOfGaussMethod();
                std::vector<T> curBasicSolution(nVars, 0);
                int indexBasicValue = 0;
                for (auto &basicVar: x) {
                    int index = basicVar - 1;
                    curBasicSolution[index] = valuesOfBasicVars[0][indexBasicValue];
                    indexBasicValue++;
                }
                allBasicSolutions.push_back(curBasicSolution);
            }

            return allBasicSolutions;
        }


        //возвращает матрицу содержащую решения системы уравнений
        //с заданной расширенной матрицей коэффициентов,
        // предварительно приведенной к треугольному виду
        Matrix reverseСourseOfGaussMethod() {
            int notFree = nRows();
            int free = nCols() - nRows();
            Matrix<T> solution(free, nRows());
            Matrix<T> triangleMatrix = *this;
            for (int i = nRows() - 1; i >= 0; --i) {
                for (int k = nCols() - 1, s = free - 1; k >= notFree && s >= 0; --k, s--) {
                    for (int j = nCols() - free - 1; j > i; --j) {
                        triangleMatrix[i][k] -= triangleMatrix[i][j] * solution[s][j];
                    }
                }
                for (int k = nCols() - 1, s = free - 1; k >= notFree && s >= 0; --k, s--) {
                    solution[s][i] = triangleMatrix[i][k] / triangleMatrix[i][i];
                }
            }
            return solution;
        }

        //возвращает определитель матрицы, полученный
        //методом Гаусса
        [[nodiscard]]T determinate() {
            if (nRows() != nCols())
                throw std::domain_error("determinant could not be found: matrix is not square");

            auto triangleMatrix = this->directСourseOfGaussMethod();
            T det = 1;
            for (int i = 0; i < nRows(); ++i) {
                det *= triangleMatrix[i][i];
            }

            return det;
        }


        //добавляет в данную матрицу столбцы матрицы add
        void addCols(Matrix<T> &add) {
            if (nRows() != add.nRows())
                throw std::runtime_error("cannot add columns: matrices have different amount of rows\n");

            for (int i = 0; i < nRows(); ++i) {
                (*this)[i].insert((*this)[i].end(), add[i].begin(), add[i].end());
            }
        }

        //удаляет из матрицы столбец с индексом pos
        void removeCol(int pos) {
            for (int i = 0; i < nRows(); ++i) {
                (*this)[i].erase((*this)[i].begin() + pos);
            }
        }

        //возвращает матрицу из одного столбца данной матрицы с индексом pos
        Matrix getCols(int pos) {
            Matrix<T> col(nRows(), 1);
            for (int i = 0; i < nRows(); ++i) {
                col[i][0] = (*this)[i][pos];
            }

            return col;
        }

        //возвращает матрицу из столбцов данной с индексами в промежутке [first,last]
        Matrix getCols(int first, int last) {
            Matrix<T> cols(nRows(), last - first + 1);
            for (int j = 0; j < nRows(); ++j)
                for (int i = first, k = 0; i <= last; ++i, ++k) {
                    cols[j][k] = (*this)[j][i];
                }

            return cols;
        }

    private:
        //возвращает индекс наибольшего элемента в столбце под номером col
        int findMaxInPartOfCol(int col) {
            T max = fabs((*this)[col][col]);
            int maxIndex = col;
            for (int i = col + 1; i < nRows(); ++i) {
                if (max < fabs((*this)[i][col])) {
                    max = fabs((*this)[i][col]);
                    maxIndex = i;
                }
            }
            return maxIndex;
        }
    };


    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    std::istream &operator>>(std::istream &in, Matrix<T> &m) {
        for (auto &row: m) {
            for (auto &value: row) {
                in >> value;
            }
        }

        return in;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    std::ostream &operator<<(std::ostream &out, Matrix<T> &m) {
        for (auto &row: m) {
            for (auto &value: row) {
                out << value << ' ';
            }
            out << '\n';
        }

        return out;
    }

    //возвращает матрицу решений системы уравнений
    // с расширенной матрицей коэффициентов coefficients
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> resolveLinearSystemByGaussMethod(Matrix<T> coefficients) {
        auto triangleMatrix = coefficients.directСourseOfGaussMethod();
        auto solution = triangleMatrix.reverseСourseOfGaussMethod();
        return solution;
    }


    //возвращает опорный план выбранный из базисных решений,
    // представленных матрицей solutions
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> getReferencePlan(Matrix<T> solutions) {
        Matrix<T> referencePlan;
        for (auto &sol: solutions) {
            bool haveNegativeVars = false;
            for (auto &var: sol) {
                if (var < 0) {
                    haveNegativeVars = true;
                    break;
                }
            }
            if (!haveNegativeVars)
                referencePlan.push_back(sol);
        }

        return referencePlan;
    }


    //выводит матрицу solution как решение системы
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void outputSolutionOfSystem(Matrix<T> solution) {
        std::cout << "\nsolutions of system:\n";
        int i = 1;
        for (auto &x: solution) {
            std::cout << i << ") ";
            i++;
            int j = 1;
            for (auto &y: x) {
                std::cout << "X" << j << "=" << y << " ";
            }
            std::cout << "\n";
        }
    }

    //выводит матрицу solution как решение системы
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void outputSolutions(Matrix<T> solution) {
        int i = 1;
        int lastVarIndex = solution.nCols() - 1;
        for (auto &x: solution) {
            std::cout << i << ") (";
            i++;
            for (int j = 0; j < solution.nCols() - 1; j++) {
                std::cout << x[j] << "; ";
            }
            std::cout << x[lastVarIndex] << ")\n";
        }
    }

    //создает единичную матрицу размерностью n
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> identityMatrix(int n) {
        Matrix<T> I(n, n, 0);
        for (int i = 0; i < n; ++i) {
            I[i][i] = 1;
        }
        return I;
    }

    //возвращает матрицу обратную для матрицы m
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> inverseMatrix(Matrix<T> &m) {
        T det = m.determinate();
        if (det == 0)
            throw std::domain_error(
                    "It is impossible to get the inverse matrix, the determinant of the matrix is zero");
        Matrix<T> I = identityMatrix<T>(m.nCols());
        Matrix<T> inverse = m;
        inverse.addCols(I);
        inverse = inverse.directСourseOfGaussMethod();
        inverse = inverse.reverseСourseOfGaussMethod();
        inverse = inverse.transpose();

        return inverse;
    }


    //создание симплекс-таблицы на основе матрицы коэффициентов системы ограничений и
    //коэффициентов целевой функции
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> createSimplexTable(Matrix<T> &matrix, std::vector<T> targetFunc) {
        Matrix<T> m = matrix;
        //приведение коэффициентов базовых переменных к единице
        std::vector<int> basic(m.nRows(), 0);
        for (int i = 0; i < targetFunc.size() - 1; ++i) {
            int count = 0;
            for (int j = 0; j < matrix.nRows(); ++j) {
                if(matrix[j][i] != 0)
                    count++;
            }
            for (int j = 0; j < matrix.nRows(); ++j) {
                if (targetFunc[i] == 0 && matrix[j][i] != 0 && count==1) {
                    T curElToBasic = m[j][i];
                    for (int l = 0; l < m.nCols(); ++l) {
                        m[j][l] /= curElToBasic;
                    }
                    basic[j] = i;
                }
            }
        }

        //заполнение таблицы коэффициентами системы ограничений
        Matrix<T> simplexTable(m.nRows() + 1, m.nCols() + 1, 0);
        for (int i = 0; i < m.nRows(); ++i) {
            for (int j = 0, l = 1; j < m.nCols(); ++j, ++l) {
                simplexTable[i][l] = m[i][j];
            }
        }

        //заполнение строки таблицы коэффициентами целевой функции
        int targetRow = simplexTable.nRows() - 1;
        for (int j = 0, l = 1; j < m.nCols()-1; ++j, ++l) {
            simplexTable[targetRow][l] = -targetFunc[j];
        }
        simplexTable[targetRow].back() = targetFunc.back();
        //заполнение столбца таблицы индексами базовых переменных
        for (int i = 0; i < m.nRows(); ++i) {
            simplexTable[i][0] = basic[i];
        }

        return simplexTable;
    }

    //переход к новой симплекс-таблице
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> substitutionInSimplexTable(Matrix<T> &simplexTable, int generalElementRow, int generalElementCol) {
        T divider = simplexTable[generalElementRow][generalElementCol];
        for (int i = 1; i < simplexTable.nCols(); ++i) {
            simplexTable[generalElementRow][i] /= divider;
        }

        simplexTable[generalElementRow][0] = generalElementCol - 1;

        for (int i = 0; i < simplexTable.nRows(); ++i) {
            if (i != generalElementRow) {
                T mul = -simplexTable[i][generalElementCol];
                for (int j = 1; j < simplexTable.nCols(); ++j) {
                    simplexTable[i][j] += simplexTable[generalElementRow][j] * mul;
                }
            }
        }

        return simplexTable;
    }


    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    T Abs(T el) {
        if (el < 0)
            return -el;
        else
            return el;
    }

    //поиск минимального отрицательного коэффициента в целевой функции
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    int targetFunctionMinNegative(Matrix<T> &simplexTable) {
        int minIndex = 1;
        int funcRow = simplexTable.nRows() - 1;
        T min = Abs(simplexTable[funcRow][minIndex]);
        bool negativeFind = simplexTable[funcRow][minIndex] < 0;
        for (int i = 2; i < simplexTable.nCols() - 1; ++i) {
            if (negativeFind && simplexTable[funcRow][i] < 0 && min < Abs(simplexTable[funcRow][i])) {
                min = Abs(simplexTable[funcRow][i]);
                minIndex = i;
            } else if (!negativeFind && simplexTable[funcRow][i] < 0) {
                min = Abs(simplexTable[funcRow][i]);
                minIndex = i;
                negativeFind = true;
            }
        }

        return minIndex;
    }

    void outputSimplexTable(Matrix<Fraction> &t);

    void outputArtificalSimplexTable(matrices::Matrix<Fraction> &t, int nX, bool method);

    //возвращает конечную симплекс таблицу для задачи ЛП
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> simplexMethod(Matrix<T> &matrix, std::vector<T> targetFunc) {
        Matrix<T> simplexTable = createSimplexTable(matrix, targetFunc);

        int targetRow = simplexTable.nRows() - 1;

        //поиск минимального отрицательного среди коэффициентов целевой функции
        int minIndexInTargetFunc = targetFunctionMinNegative(simplexTable);

        outputSimplexTable(simplexTable);
        while (simplexTable[targetRow][minIndexInTargetFunc] < 0) {

            int generalElementIndexRow = -1;

            //поиск первого положительного элемента в столбце генерального элемента
            for (int i = 0; i < targetRow; ++i) {
                if (simplexTable[i][minIndexInTargetFunc] > 0) {
                    generalElementIndexRow = i;
                    break;
                }
            }

            //исключительная ситуация -  в столбце нет положительных элементов
            if (generalElementIndexRow == -1) {
                throw std::runtime_error("the system has no solutions");
            }

            //поиск генерального элемента
            T relation = simplexTable[generalElementIndexRow][simplexTable.nCols() - 1] /
                         simplexTable[generalElementIndexRow][minIndexInTargetFunc];
            for (int i = generalElementIndexRow + 1; i < targetRow; ++i) {
                if (simplexTable[i][minIndexInTargetFunc] > 0 &&
                    simplexTable[i][simplexTable.nCols() - 1] / simplexTable[i][minIndexInTargetFunc] < relation) {
                    generalElementIndexRow = i;
                    relation = simplexTable[generalElementIndexRow][simplexTable.nCols() - 1] /
                               simplexTable[generalElementIndexRow][minIndexInTargetFunc];
                }
            }

            simplexTable = substitutionInSimplexTable(simplexTable, generalElementIndexRow,
                                                      minIndexInTargetFunc);
            outputSimplexTable(simplexTable);
            //поиск минимального отрицательного среди коэффициентов целевой функции
            minIndexInTargetFunc = targetFunctionMinNegative(simplexTable);
        }

        //возврат конечной симплекс-таблицы
        return simplexTable;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> simplexMethodforArtificalBasis(Matrix<T> &matrix, std::vector<T> targetFunc, int nX, bool method) {
        Matrix<T> simplexTable = createSimplexTable(matrix, targetFunc);

        int targetRow = simplexTable.nRows() - 1;

        //поиск минимального отрицательного среди коэффициентов целевой функции
        int minIndexInTargetFunc = targetFunctionMinNegative(simplexTable);

        outputArtificalSimplexTable(simplexTable, nX, method);
        while (simplexTable[targetRow][minIndexInTargetFunc] < 0) {

            int generalElementIndexRow = -1;

            //поиск первого положительного элемента в столбце генерального элемента
            for (int i = 0; i < targetRow; ++i) {
                if (simplexTable[i][minIndexInTargetFunc] > 0) {
                    generalElementIndexRow = i;
                    break;
                }
            }

            //исключительная ситуация -  в столбце нет положительных элементов
            if (generalElementIndexRow == -1) {
                throw std::runtime_error("the system has no solutions");
            }

            //поиск генерального элемента
            T relation = simplexTable[generalElementIndexRow][simplexTable.nCols() - 1] /
                         simplexTable[generalElementIndexRow][minIndexInTargetFunc];
            for (int i = generalElementIndexRow + 1; i < targetRow; ++i) {
                if (simplexTable[i][minIndexInTargetFunc] > 0 &&
                    simplexTable[i][simplexTable.nCols() - 1] / simplexTable[i][minIndexInTargetFunc] < relation) {
                    generalElementIndexRow = i;
                    relation = simplexTable[generalElementIndexRow][simplexTable.nCols() - 1] /
                               simplexTable[generalElementIndexRow][minIndexInTargetFunc];
                }
            }

            simplexTable = substitutionInSimplexTable(simplexTable, generalElementIndexRow,
                                                      minIndexInTargetFunc);
            outputArtificalSimplexTable(simplexTable, nX, method);
            //поиск минимального отрицательного среди коэффициентов целевой функции
            minIndexInTargetFunc = targetFunctionMinNegative(simplexTable);
        }

        //возврат конечной симплекс-таблицы
        return simplexTable;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> simplexMethod(Matrix<T> &simplexTable, int nX, bool method) {

        int targetRow = simplexTable.nRows() - 1;

        //поиск минимального отрицательного среди коэффициентов целевой функции
        int minIndexInTargetFunc = targetFunctionMinNegative(simplexTable);

        outputArtificalSimplexTable(simplexTable, nX, method);
        while (simplexTable[targetRow][minIndexInTargetFunc] < 0) {

            int generalElementIndexRow = -1;

            //поиск первого положительного элемента в столбце генерального элемента
            for (int i = 0; i < targetRow; ++i) {
                if (simplexTable[i][minIndexInTargetFunc] > 0) {
                    generalElementIndexRow = i;
                    break;
                }
            }

            //исключительная ситуация -  в столбце нет положительных элементов
            if (generalElementIndexRow == -1) {
                throw std::runtime_error("the system has no solutions");
            }

            //поиск генерального элемента
            T relation = simplexTable[generalElementIndexRow][simplexTable.nCols() - 1] /
                         simplexTable[generalElementIndexRow][minIndexInTargetFunc];
            for (int i = generalElementIndexRow + 1; i < targetRow; ++i) {
                if (simplexTable[i][minIndexInTargetFunc] > 0 &&
                    simplexTable[i][simplexTable.nCols() - 1] / simplexTable[i][minIndexInTargetFunc] < relation) {
                    generalElementIndexRow = i;
                    relation = simplexTable[generalElementIndexRow][simplexTable.nCols() - 1] /
                               simplexTable[generalElementIndexRow][minIndexInTargetFunc];
                }
            }

            simplexTable = substitutionInSimplexTable(simplexTable, generalElementIndexRow,
                                                      minIndexInTargetFunc);
            outputArtificalSimplexTable(simplexTable, nX, method);
            //поиск минимального отрицательного среди коэффициентов целевой функции
            minIndexInTargetFunc = targetFunctionMinNegative(simplexTable);
        }

        //возврат конечной симплекс-таблицы
        return simplexTable;
    }

    //выводит оптимальное решение задачи ЛП найденное симплекс-методом в чистом виде
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void outputSolutionSimplexMethod(Matrix<T> simplexTable, int nX) {
        std::vector<T> solution(simplexTable.nCols() - 2, 0);
        int n = simplexTable.nRows() - 1;
        int k = simplexTable.nCols() - 1;
        for (int i = 0; i < n; ++i) {
            int var = static_cast<int>(simplexTable[i][0]);
            solution[var] = simplexTable[i][k];
        }

        k = nX - 1;
        std::cout << "Solution - ( ";
        for (int i = 0; i < k; ++i) {
            std::cout << solution[i] << "; ";
        }
        std::cout << solution[nX - 1] << " )\nz_max = " << simplexTable[n][simplexTable.nCols() - 1];
    }

    //выводит оптимальное решение задачи ЛП найденное симплекс-методом в чистом виде
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void outputSolutionSimplexMethodFractional(Matrix<T> simplexTable, int nX) {
        std::vector<T> solution(simplexTable.nCols() - 2, 0);
        int n = simplexTable.nRows() - 1;
        int k = simplexTable.nCols()-1;
        for (int i = 0; i < n; ++i) {
            int var = static_cast<int>(simplexTable[i][0]);
            solution[var] = simplexTable[i][k];
        }

        k = nX-1 ;
        std::cout << "Solution:\n( ";
        for (int i = 0; i < k; ++i) {
            std::cout << solution[i] << "; ";
        }
        std::cout << solution[nX - 1] <<" )  \n(";

        for (int i = 1; i < k; ++i) {
            std::cout << solution[i]<<" / "<<solution[0] << "; ";
        }
        std::cout << solution[nX - 1]<<" / "<<solution[0] << " )\n";

        std::cout << "\nAnswer - ( ";
        for (int i = 1; i < k; ++i) {
            std::cout << solution[i] / solution[0] << "; ";
        }
        std::cout << solution[nX - 1] / solution[0] << " )\nz_max = " << simplexTable[n][simplexTable.nCols() - 1];
    }
    //выводит оптимальное решение задачи ЛП найденное симплекс-методом в чистом виде
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void outputSolutionSimplexMethodDual(Matrix<T> simplexTable, int nX) {
        std::vector<T> solution(simplexTable.nCols() - 2, 0);
        int n = simplexTable.nRows() - 1;
        int k = simplexTable.nCols() - 1;
        for (int i = 0; i < n; ++i) {
            int var = static_cast<int>(simplexTable[i][0]);
            solution[var] = simplexTable[i][k];
        }


        std::cout << "Solution:\n( ";
        for (int i = 0; i < solution.size() - 1; ++i) {
            std::cout << solution[i] << "; ";
        }
        std::cout << solution.back() << " )\nf_max = " << simplexTable[n][simplexTable.nCols() - 1];


        std::cout << "\n( ";
        for (int i = nX + 1; i < simplexTable.nCols() - 1; ++i) {
            std::cout << simplexTable.back()[i] << "; ";
        }
        for (int i = 1; i < nX; ++i) {
            std::cout << simplexTable.back()[i] << "; ";
        }
        std::cout << simplexTable.back()[nX] << " )\nz_min = " << simplexTable[n][simplexTable.nCols() - 1];

    }

    //выводит оптимальное решение задачи ЛП найденное симплекс-методом в чистом виде
    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void outputSolutionMatrixGame(Matrix<T> simplexTable, Matrix<T> matrix, int nX, bool type) {
        if (!type) {
            std::cout << "Solution:\n";
            std::cout << "\nv = " << static_cast<T>(1) / simplexTable[0][0] << "\nq (";
            for (int i = 0; i < matrix.nRows() - 1; ++i) {
                std::cout << simplexTable[1][i] << "; ";
            }
            std::cout << simplexTable[1][matrix.nRows() - 1] << " )";

            std::cout << "\np (";
            for (int i = 0; i < matrix.nCols() - 1; ++i) {
                std::cout << simplexTable[2][i] << "; ";
            }
            std::cout << simplexTable[2][matrix.nCols() - 1] << " )\n";

        } else {
            std::vector<T> solution(simplexTable.nCols() - 2, 0);
            int n = simplexTable.nRows() - 1;
            int k = simplexTable.nCols() - 1;
            for (int i = 0; i < n; ++i) {
                int var = static_cast<int>(simplexTable[i][0]);
                solution[var] = simplexTable[i][k];
            }


            T v = static_cast<T>(1) / simplexTable[n][simplexTable.nCols() - 1];
            std::cout << "Solution:\n";
            std::cout << "\nv = " << v << "\nq (";
            for (int i = 0; i < matrix.nCols() - 1; ++i) {
                std::cout << solution[i] * v << "; ";
            }
            std::cout << solution.back() * v << " )\n ";


            std::cout << "\np ( ";
            for (int i = nX + 1; i < simplexTable.nCols() - 2; ++i) {
                std::cout << simplexTable.back()[i] * v << "; ";
            }
            std::cout << simplexTable.back()[simplexTable.nCols() - 2] * v << " )";
        }
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void addVar(Matrix<T> &matrix, int index) {
        int newSize = matrix.nCols() + 1;
        for (int i = 0; i < matrix.nRows(); ++i) {
            matrix[i].resize(newSize);
            matrix[i].back() = matrix[i][matrix.nCols() - 2];
            matrix[i][matrix.nCols() - 2] = i == index ? 1 : 0;
        }
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> bringСonditionsToEquality(Matrix<T> &matrix, std::vector<std::string> &signs, int &additionalXIndex) {
        additionalXIndex = matrix.nCols() - 2;
        auto copyMatrix = matrix;
        for (int i = 0; i < signs.size(); ++i) {
            if (signs[i] == "<=") {
                addVar(copyMatrix, i);
                additionalXIndex++;
            } else if (signs[i] == ">=") {
                copyMatrix.mulRow(i, -1);
                addVar(copyMatrix, i);
                additionalXIndex++;
            }
        }
        return copyMatrix;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> addArtificalVars(Matrix<T> matrix, std::vector<std::string> &signs, std::vector<T> targetFunc,
                               int &additionalYIndex,std::vector<bool> &basic) {
        //приведение коэффициентов базовых переменных к единице
        std::vector<bool> basicHave(matrix.nRows(), false);
        for (int i = 0; i < targetFunc.size() - 1; ++i) {
            int count = 0;
            for (int j = 0; j < matrix.nRows()-1; ++j) {
                if(matrix[j][i] != 0)
                    count++;
            }
            for (int j = 0; j < matrix.nRows(); ++j) {
                if (targetFunc[i] == 0 && matrix[j][i] != 0 && count==1) {
                    basicHave[j] = true;
                }
            }
        }

        basic=basicHave;

        additionalYIndex = 0;
        auto copyMatrix = matrix;
        for (int i = 0; i < signs.size(); ++i) {
            if (signs[i] == "=" && !basicHave[i]) {
                addVar(copyMatrix, i);
                additionalYIndex++;
            }
        }
        return copyMatrix;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    std::vector<T> getArtificalTargetFunc(Matrix<T> matrix,  std::vector<std::string> &signs,int nX,std::vector<bool> &basic) {

        std::vector<T> targetFuncNew(matrix.nCols(), 0);
        for (int i = 0; i < signs.size(); ++i) {
            if (signs[i] == "=" && !basic[i]) {
                for (int j = 0; j < nX; ++j) {
                    targetFuncNew[j] += matrix[i][j];
                }
                targetFuncNew.back() += matrix[i].back();
            }
        }

        return targetFuncNew;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    std::vector<T>
    getTargetFuncForMethodOfLargeFines(Matrix<T> matrix, std::vector<T> targetFunc, std::vector<std::string> &signs,
                                       int nX) {
        T mul = 0;
        for (int i = 0; i < matrix.nRows(); ++i) {
            mul += matrix[i].back();
        }
        mul *= 10;
        std::vector<T> targetFuncNew(matrix.nCols(), 0);
        for (int i = 0; i < signs.size(); ++i) {
            if (signs[i] == "=") {
                for (int j = 0; j < nX; ++j) {
                    targetFuncNew[j] += mul * matrix[i][j];
                }
                targetFuncNew.back() += mul * matrix[i].back();
            }
        }

        for (int i = 0; i < targetFunc.size() - 1; ++i) {
            targetFuncNew[i] += targetFunc[i];
        }
        targetFuncNew.back() += targetFunc.back();

        return targetFuncNew;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void getNewBasisFromSimplexTable(Matrix<T> &matrix, int &additionalYIndex) {
        int newSize = matrix.nCols() - additionalYIndex;
        for (int i = 0; i < matrix.nRows(); ++i) {
            matrix[i][newSize - 1] = matrix[i].back();
            matrix[i].resize(newSize);
        }
        matrix.pop_back();
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    std::vector<T> getTargetFunc(Matrix<T> matrix, std::vector<T> targetFunc) {
        std::vector<T> newTargetFunc(matrix.nCols() - 1, 0);
        for (int i = 0; i < targetFunc.size() - 1; ++i) {
            newTargetFunc[i] = targetFunc[i];
        }
        newTargetFunc.back() = targetFunc.back();
        targetFunc = newTargetFunc;
        for (int i = 0; i < matrix.nRows(); ++i) {
            int indexMul = (int) (matrix[i].front());
            T mul = targetFunc[indexMul];
            newTargetFunc[indexMul] = 0;
            for (int j = 0; j < targetFunc.size() - 1; ++j) {
                if (j != indexMul)
                    newTargetFunc[j] -= matrix[i][j + 1] * mul;
            }
            newTargetFunc.back() -= matrix[i].back() * mul;
        }
        return newTargetFunc;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> artificialBasisMethod(Matrix<T> matrix, std::vector<T> targetFunc, std::vector<std::string> &signs) {
        int additionalXIndex = matrix.nCols() - 2;
        auto copyMatrix = bringСonditionsToEquality(matrix, signs, additionalXIndex);
        int additionalYIndex = 0;
        int nX = additionalXIndex + 1;
        std::vector<bool> basic;
        copyMatrix = addArtificalVars(copyMatrix, signs, targetFunc, additionalYIndex, basic);

        auto artificalTargetFunction = getArtificalTargetFunc(copyMatrix, signs, nX, basic);

        auto table = simplexMethodforArtificalBasis(copyMatrix, artificalTargetFunction, nX, false);

        getNewBasisFromSimplexTable(table, additionalYIndex);
        auto newTargetFunc = getTargetFunc(table, targetFunc);

        for (int i = 0; i < table.nRows(); ++i) {
            table[i].erase(table[i].begin());
        }

        table = simplexMethod(table, newTargetFunc);

        return table;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> methodOfLargeFines(Matrix<T> matrix, std::vector<T> targetFunc, std::vector<std::string> &signs) {
        int additionalXIndex = matrix.nCols() - 2;
        auto copyMatrix = bringСonditionsToEquality(matrix, signs, additionalXIndex);
        int additionalYIndex = 0;
        int nX = additionalXIndex + 1;
        std::vector<bool> basic;
        copyMatrix = addArtificalVars(copyMatrix, signs, targetFunc, additionalYIndex, basic);
        auto newTargetFunction = getTargetFuncForMethodOfLargeFines(copyMatrix, targetFunc, signs, nX);

        auto table = simplexMethodforArtificalBasis(copyMatrix, newTargetFunction, nX, true);

        return table;
    }

    template<class T> requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    struct Task {
        Matrix<T> matrix;
        std::vector<std::string> signs;
        std::vector<T> targetFunc;
        std::string typeOfTask;
    };

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Task<T> getDualTask(Task<T> task) {
        Matrix<T> matrix = task.matrix;
        matrix.push_back(task.targetFunc);
        matrix = matrix.transpose();
        Task<T> dualTask;
        dualTask.targetFunc = matrix.back();
        matrix.pop_back();
        dualTask.matrix = matrix;
        std::string sign;
        if (task.typeOfTask == "max") {
            sign = ">=";
            dualTask.typeOfTask = "min";
        } else {
            sign = "<=";
            dualTask.typeOfTask = "max";
        }

        std::vector<std::string> signs(matrix.nRows(), sign);
        dualTask.signs = signs;

        return dualTask;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Task<T> getStandardFormOfTask(Task<T> task) {
        Matrix<T> matrix;
        std::string signOfType = task.typeOfTask == "max" ? "<=" : ">=";
        T mul = static_cast<T>(-1);
        for (int i = 0; i < task.signs.size(); ++i) {
            if (task.signs[i] == "=") {
                matrix.push_back(task.matrix[i]);
                matrix.push_back(task.matrix[i] * mul);
            } else if (signOfType != task.signs[i])
                matrix.push_back(task.matrix[i] * mul);
            else
                matrix.push_back(task.matrix[i]);
        }

        std::vector<std::string> signs(matrix.nRows(), signOfType);
        task.signs = signs;
        task.matrix = matrix;

        return task;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void outputTask(Task<T> task, bool isDual) {
        char letter;
        char funcName;
        if (isDual) {
            std::cout << "DUAL TASK\n";
            letter = 'y';
            funcName = 'f';
        } else {
            std::cout << "INITIAL TASK\n";
            letter = 'x';
            funcName = 'z';
        }

        for (int i = 0; i < task.matrix.nRows(); ++i) {
            for (int j = 0; j < task.matrix.nCols() - 1; ++j) {
                if (task.matrix[i][j] == 0)
                    continue;
                else if (task.matrix[i][j] < 0)
                    std::cout << " - ";
                else if (j != 0)
                    std::cout << " + ";
                if (Abs(task.matrix[i][j]) != 1)
                    std::cout << Abs(task.matrix[i][j]);
                std::cout << letter << '_' << j + 1;
            }
            std::cout << ' ' << task.signs[i];
            if (task.matrix[i].back() < 0)
                std::cout << " - " << Abs(task.matrix[i].back());
            else
                std::cout << ' ' << Abs(task.matrix[i].back());
            std::cout << '\n';
        }

        std::cout << '\n' << funcName << " = ";
        for (int i = 0; i < task.targetFunc.size() - 1; ++i) {
            if (task.targetFunc[i] == 0)
                continue;
            else if (task.targetFunc[i] < 0)
                std::cout << " - ";
            else if (i != 0)
                std::cout << " + ";
            if (Abs(task.targetFunc[i]) != 1)
                std::cout << Abs(task.targetFunc[i]);
            std::cout << letter << '_' << i + 1;
        }

        if (task.targetFunc.back() < 0) {
            std::cout << " - " << Abs(task.targetFunc.back());
        } else if (task.targetFunc.back() > 0) {
            std::cout << " + " << Abs(task.targetFunc.back());
        }

        std::cout << " ---> " << task.typeOfTask << "\n\n";
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    int getMinIndexInFreeCoefficientsCol(Matrix<T> simplexTable) {
        T min = simplexTable[0].back();
        int minIndex = 0;
        for (int i = 1; i < simplexTable.nRows() - 1; ++i) {
            T curEl = simplexTable[i].back();
            if (curEl < min) {
                min = curEl;
                minIndex = i;
            }
        }

        return minIndex;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> methodOfRefiningEstimates(Matrix<T> simplexTable) {

        int targetRow = simplexTable.nRows() - 1;

        //поиск минимального отрицательного среди коэффициентов целевой функции
        int minIndexInFreeCoefficientsCol = getMinIndexInFreeCoefficientsCol(simplexTable);


        outputSimplexTable(simplexTable);
        while (simplexTable[minIndexInFreeCoefficientsCol].back() < 0) {

            int generalElementIndexCol = -1;

            //поиск первого отрицательного элемента в строке генерального элемента
            for (int i = 1; i < simplexTable.nCols() - 1; ++i) {
                if (simplexTable[minIndexInFreeCoefficientsCol][i] < 0) {
                    generalElementIndexCol = i;
                    break;
                }
            }

            //исключительная ситуация - в строке нет отрицательных элементов
            if (generalElementIndexCol == -1) {
                throw std::runtime_error("the system has no solutions");
            }

            //поиск генерального элемента в строке
            T relation = Abs(simplexTable.back()[generalElementIndexCol] /
                             simplexTable[minIndexInFreeCoefficientsCol][generalElementIndexCol]);
            for (int i = generalElementIndexCol + 1; i < simplexTable.nCols() - 1; ++i) {
                if (simplexTable[minIndexInFreeCoefficientsCol][i] < 0 &&
                    Abs(simplexTable.back()[i] /
                        simplexTable[minIndexInFreeCoefficientsCol][i]) < relation) {
                    generalElementIndexCol = i;
                    relation = Abs(simplexTable.back()[generalElementIndexCol] /
                                   simplexTable[minIndexInFreeCoefficientsCol][generalElementIndexCol]);
                }
            }

            simplexTable = substitutionInSimplexTable(simplexTable, minIndexInFreeCoefficientsCol,
                                                      generalElementIndexCol);
            outputSimplexTable(simplexTable);
            //поиск минимального отрицательного среди коэффициентов целевой функции
            minIndexInFreeCoefficientsCol = getMinIndexInFreeCoefficientsCol(simplexTable);
        }

        //возврат конечной симплекс-таблицы
        return simplexTable;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> methodOfRefiningEstimates(Matrix<T> simplexTable, int nX, bool method) {

        int targetRow = simplexTable.nRows() - 1;

        //поиск минимального отрицательного среди коэффициентов целевой функции
        int minIndexInFreeCoefficientsCol = getMinIndexInFreeCoefficientsCol(simplexTable);


        outputArtificalSimplexTable(simplexTable, nX, method);
        while (simplexTable[minIndexInFreeCoefficientsCol].back() < 0) {

            int generalElementIndexCol = -1;

            //поиск первого отрицательного элемента в строке генерального элемента
            for (int i = 1; i < simplexTable.nCols() - 1; ++i) {
                if (simplexTable[minIndexInFreeCoefficientsCol][i] < 0) {
                    generalElementIndexCol = i;
                    break;
                }
            }

            //исключительная ситуация - в строке нет отрицательных элементов
            if (generalElementIndexCol == -1) {
                throw std::runtime_error("the system has no solutions");
            }

            //поиск генерального элемента в строке
            T relation = Abs(simplexTable.back()[generalElementIndexCol] /
                             simplexTable[minIndexInFreeCoefficientsCol][generalElementIndexCol]);
            for (int i = generalElementIndexCol + 1; i < simplexTable.nCols() - 1; ++i) {
                if (simplexTable[minIndexInFreeCoefficientsCol][i] < 0 &&
                    Abs(simplexTable.back()[i] /
                        simplexTable[minIndexInFreeCoefficientsCol][i]) < relation) {
                    generalElementIndexCol = i;
                    relation = Abs(simplexTable.back()[generalElementIndexCol] /
                                   simplexTable[minIndexInFreeCoefficientsCol][generalElementIndexCol]);
                }
            }

            simplexTable = substitutionInSimplexTable(simplexTable, minIndexInFreeCoefficientsCol,
                                                      generalElementIndexCol);
            outputArtificalSimplexTable(simplexTable, nX, method);
            //поиск минимального отрицательного среди коэффициентов целевой функции
            minIndexInFreeCoefficientsCol = getMinIndexInFreeCoefficientsCol(simplexTable);
        }

        //возврат конечной симплекс-таблицы
        return simplexTable;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Task<T> getCanonicalForm(Task<T> task) {
        int additionalXIndex = task.matrix.nCols() - 2;
        auto copyMatrix = bringСonditionsToEquality(task.matrix, task.signs, additionalXIndex);
        task.matrix = copyMatrix;
        std::vector<std::string> signs(copyMatrix.nRows(), "=");
        task.signs = signs;
        std::vector<T> func(copyMatrix.nCols(), 0);
        for (int i = 0; i < task.targetFunc.size() - 1; ++i) {
            func[i] = task.targetFunc[i];
        }
        func.back() = task.targetFunc.back();

        task.targetFunc = func;

        return task;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> dualSimplexMethod(Task<T> task) {
        auto initTask = getStandardFormOfTask(task);
        auto dualTask = getDualTask(initTask);
        outputTask(initTask, false);
        outputTask(dualTask, true);

        Matrix<T> simplexTable;


        std::cout << "Simplex method\n";
        if (initTask.typeOfTask == "max") {
            initTask = getCanonicalForm(initTask);
            simplexTable = simplexMethod(initTask.matrix, initTask.targetFunc);
        } else {
            dualTask = getCanonicalForm(dualTask);
            simplexTable = simplexMethod(dualTask.matrix, dualTask.targetFunc);
        }
        std::cout << "Method of refining estimates\n";
        simplexTable = methodOfRefiningEstimates(simplexTable);

        return simplexTable;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> pureStrategies(Matrix<T> matrix) {
        Matrix transpose = matrix.transpose();

        std::vector<T> minInRows(matrix.nRows());
        std::vector<T> maxInCols(matrix.nCols());

        int bound1 = (matrix.nRows() < matrix.nCols()) ? matrix.nRows() : matrix.nCols();
        int bound2 = (matrix.nRows() < matrix.nCols()) ? matrix.nCols() : matrix.nRows();

        for (int i = 0; i < bound1; i++) {
            minInRows[i] = *(std::min_element(matrix[i].begin(),
                                              matrix[i].end()));
            maxInCols[i] = *(std::max_element(transpose[i].begin(),
                                              transpose[i].end()));
        }

        if (bound1 == matrix.nRows()) {
            for (int i = bound1; i < bound2; i++) {
                maxInCols[i] = *(std::max_element(transpose[i].begin(),
                                                  transpose[i].end()));
            }
        } else {
            for (int i = bound1; i < bound2; i++) {
                minInRows[i] = *(std::min_element(matrix[i].begin(),
                                                  matrix[i].end()));
            }
        }

        T maxMin = *(std::max_element(minInRows.begin(),
                                      minInRows.end()));
        T minMax = *(std::min_element(maxInCols.begin(),
                                      maxInCols.end()));
        Matrix<T> solution(3, (matrix.nRows() < matrix.nCols()) ? matrix.nCols() : matrix.nRows(), 0);
        if (maxMin == minMax) {

            solution[0][0] = maxMin;

            for (int i = 0; i < matrix.nRows(); ++i) {
                auto result = std::find(matrix[i].begin(), matrix[i].end(), maxMin);
                if (result != matrix[i].end())
                    solution[1][i] = 1;
                else
                    solution[1][i] = 0;
            }
            for (int i = 0; i < matrix.nCols(); ++i) {
                auto result = std::find(transpose[i].begin(), transpose[i].end(), maxMin);
                if (result != transpose[i].end())
                    solution[2][i] = 1;
                else
                    solution[2][i] = 0;
            }
            return solution;
        }
        Matrix<T> empty;
        return empty;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> mixedStrategies(Matrix<T> matrix) {
        for (auto &row: matrix) {
            row.push_back(1);
        }
        std::vector<T> targetFunc(matrix.nCols(), 1);
        targetFunc.back() = 0;

        std::vector<std::string> signs(matrix.nRows(), "<=");
        Task<T> task = {matrix, signs, targetFunc, "max"};

        return dualSimplexMethod(task);
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> getSaddlePointOfTheGame(Matrix<T> matrix, bool &type) {
        Matrix<T> solution = pureStrategies(matrix);
        if (!solution.empty()) {
            type = false;
            return solution;
        } else {
            type = true;
            solution = mixedStrategies(matrix);
            return solution;
        }
        return solution;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    bool isIntegerSolution(Matrix<T> table) {
        for (int i = 0; i < table.nRows(); ++i) {
            int digit = static_cast<int>(table[i].back());
            T digitT = static_cast<T>(digit);
            if (table[i].back() != digitT) {
                return false;
            }
        }
        return true;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    int findTheMostFractionIndex(Matrix<T> table) {
        T max = Abs(table[0].back() - static_cast<int>(table[0].back()));
        int index = 0;
        for (int i = 0; i < table.nRows() - 1; ++i) {
            T cur = Abs(table[i].back() - static_cast<int>(table[i].back()));
            if (cur > max) {
                max = cur;
                index = i;
            }
        }
        return index;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void addGomorySection(Matrix<T> &table, int index) {
        int newSize = table.nCols() + 1;
        for (int i = 0; i < table.nRows(); ++i) {
            table[i].resize(newSize);
            table[i].back() = table[i][table.nCols() - 2];
            table[i][table.nCols() - 2] = 0;
        }
        std::vector<T> basic;
        for (int i = 0; i < table.nRows() - 1; ++i) {
            basic.push_back(table[i].front());
        }
        std::vector<T> section(newSize);
        section[0] = newSize - 3;
        for (int i = 1; i < section.size(); ++i) {
            T ind = static_cast<T>(i - 1);
            if (std::find(basic.begin(), basic.end(), ind) != basic.end())
                section[i] = 0;
            else {
                section[i] = (table[index][i] - static_cast<int>(table[index][i])) * -1;
            }
        }
        section[newSize - 2] = 1;

        table.insert(table.begin() + table.nRows() - 1, section);
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T> firstGomoryMethod(Task<T> task) {
        int nX = task.matrix.nCols() - 1;
        std::cout << "\n\n______Simplex method______\n";
        auto table = simplexMethod(task.matrix, task.targetFunc);
        int maxXIndex = findTheMostFractionIndex(table);
        while (table[maxXIndex].back() != static_cast<T> (static_cast<int>( table[maxXIndex].back())) ||
               table.back().back() != static_cast<T> (static_cast<int>( table.back().back()))) {
            addGomorySection(table, maxXIndex);
            std::cout << "______Simplex method______\n";
            table = simplexMethod(table, nX, 0);
            std::cout << "______Method of refining estimates______\n";
            table = methodOfRefiningEstimates(table, nX, 0);
            maxXIndex = findTheMostFractionIndex(table);
        }
        return table;
    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    Matrix<T>
    solvingFractionalLinearProgrammingProblem(Matrix<T> matrix, Matrix<T> targetFunc, std::vector<std::string> &sings) {
        int nX = matrix.nCols() - 1;
        matrix = bringСonditionsToEquality(matrix, sings, nX);

        std::vector<bool> funcVars(matrix.nCols() - 1, false);
        for (int i = 0; i < targetFunc.nRows(); ++i) {
            for (int j = 0; j < targetFunc.nCols(); ++j) {
                if (targetFunc[i][j] != 0)
                    funcVars[j] = true;
            }
        }
        matrix.push_back(targetFunc.back());
        for (int i = 0; i < matrix.nRows(); ++i) {
            matrix[i].insert(matrix[i].begin(), matrix[i].back()*-1);
            matrix[i].back() = 0;
        }
        matrix.back().back() = 1;
        std::vector<T> func = targetFunc.front();
        func.insert(func.begin(), func.back());
        func.back() = 0;


        std::vector<std::string> singsNew(matrix.nRows(), "=");
        return artificialBasisMethod(matrix, func, singsNew);

    }

    template<typename T>
    requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
    void inputSystem(matrices::Matrix<T> &matrix, std::vector<std::string> &signs) {
        for (int i = 0; i < matrix.nRows(); ++i) {
            for (int j = 0; j < matrix.nCols() - 1; ++j) {
                std::cin >> matrix[i][j];
            }
            std::cin >> signs[i];
            std::cin >> matrix[i].back();
        }
    }


}

template<typename T>
requires std::is_arithmetic_v<T> || std::is_same_v<T, Fraction>
std::istream &operator>>(std::istream &in, std::vector<T> &v) {
    for (auto &x: v) {
        in >> x;
    }
    return in;
}

#endif //IOP_MATRIX_H
