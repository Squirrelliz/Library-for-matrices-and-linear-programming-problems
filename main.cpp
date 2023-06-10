#include <windows.h>
#include <vector>
#include "Matrix.h"


int main() {
    SetConsoleOutputCP(CP_UTF8);
    int nRows, nCols;

    std::cout << "Введите размеры расширенной матрицы коэффициентов системы ограничений\n";

    std::cin >> nRows >> nCols;
    matrices::Matrix<Fraction> matrix(nRows, nCols);
    std::vector<std::string> signs(nRows);
    std::cout
            << "Введите расширенную матрицу коэффициентов системы ограничений вместе со знаком для каждой строки (равенство или неравенство)\n";

    matrices::inputSystem(matrix, signs);

    std::cout << "Введите коэффициенты целевой функции\n";

    std::vector<Fraction> targetFunc(nCols);
    std::cin >> targetFunc;

    std::string type;
    std::cout << "Введите тип задачи\n";
    std::cin >> type;
    matrices::Task<Fraction> task{matrix, signs, targetFunc, type};
    try {
        matrices::Matrix<Fraction> simplexTable = matrices::firstGomoryMethod(task);
        matrices::outputSolutionSimplexMethod(simplexTable,nCols-1);
    } catch (std::runtime_error &s) {
        std::cout << s.what();
    }

    return 0;
}

