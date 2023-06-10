//
// Created by squir on 19.03.2023.
//

#include "Matrix.h"

//генерирует сочетания из n по k
void
combinations_(std::set<std::set<int>> &combinations, std::vector<int> &currentComb, int i, int b, int n, int k) {

    for (int j = b; j <= n - k + i; ++j) {
        currentComb[i - 1] = j;
        if (i == k) {
            std::set<int> curCombSet;
            std::for_each(currentComb.begin(),
                          currentComb.end(),
                          [&](const auto &elem) {
                              curCombSet.insert(elem);
                          });
            combinations.insert(curCombSet);
        } else
            combinations_(combinations, currentComb, i + 1, j + 1, n, k);
    }
}

//генерирует сочетания из n по k
std::set<std::set<int>> combinations(int n, int k) {
    std::set<std::set<int>> combinations;
    std::vector<int> curComb(k);
    combinations_(combinations, curComb, 1, 1, n, k);
    return combinations;
}

//вывод симплекс-таблицы
void matrices::outputSimplexTable(matrices::Matrix<Fraction> &t) {
    std::cout << "Basic Var." << std::setw(7);
    for (int i = 0; i < t.nCols() - 2; ++i) {
        std::cout << "X_" << i << std::setw(8);
    }
    std::cout << std::setw(14) << "Free Term\n";

    for (int i = 0; i < t.nRows() - 1; ++i) {
        std::cout << 'X' << t[i].front() << std::setw(14);
        for (int j = 1; j < t.nCols() - 1; ++j) {
            std::cout << t[i][j];
            int space = t[i][j].denominator() == 1 ? 9 : 7;
            std::cout << std::setw(space);
        }
        std::cout << std::setw(10) << t[i][t.nCols() - 1] << '\n';
    }

    std::cout << 'Z' << std::setw(15);
    for (int i = 1; i < t.nCols() - 1; ++i) {
        std::cout << t[t.nRows() - 1][i];
        int space = t[t.nRows() - 1][i].denominator() == 1 ? 9 : 7;
        std::cout << std::setw(space);
    }
    std::cout << std::setw(10) << t[t.nRows() - 1][t.nCols() - 1] << '\n';
    std::cout << "\n\n";
}


//вывод симплекс-таблицы
void matrices::outputArtificalSimplexTable(matrices::Matrix<Fraction> &t, int nX, bool method) {
    std::cout << "Basic Var." << std::setw(7);
    for (int i = 0; i < nX; ++i) {
        std::cout << "X_" << i << std::setw(8);
    }
    for (int i = 0; i < t.nCols() - 2 - nX; ++i) {
        std::cout << "U_" << i << std::setw(8);
    }

    std::cout << std::setw(14) << "Free Term\n";

    for (int i = 0; i < t.nRows() - 1; ++i) {
        if (t[i].front() < nX)
            std::cout << 'X' << t[i].front() << std::setw(14);
        else
            std::cout << 'U' << t[i].front() - nX << std::setw(14);
        for (int j = 1; j < t.nCols() - 1; ++j) {
            std::cout << t[i][j];
            int space = t[i][j].denominator() == 1 ? 9 : 7;
            std::cout << std::setw(space);
        }
        std::cout << std::setw(10) << t[i][t.nCols() - 1] << '\n';
    }

    if(!method) {
        std::cout << 'Z' << std::setw(15);
    }else
        std::cout << "Z_m" << std::setw(15);
    for (int i = 1; i < t.nCols() - 1; ++i) {
        std::cout << t[t.nRows() - 1][i];
        int space = t[t.nRows() - 1][i].denominator() == 1 ? 9 : 7;
        std::cout << std::setw(space);
    }
    std::cout << std::setw(10) << t[t.nRows() - 1][t.nCols() - 1] << '\n';
    std::cout << "\n\n";
}
