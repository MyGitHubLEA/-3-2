
#include "SquareMatrix.h"
#include "Output.h"
#include "UI.h"

bool IfCorrectIndex(ArrayRectangularMatrix<int> matrix, int indexR, int indexC = 0)
{
    if (indexR >= matrix.GetNumberRows() || indexR < 0 || indexC >= matrix.GetNumberColumns() || indexC < 0)
    {
        std::cout << "Incorrect index. Try again!" << std::endl;

        return false;
    }
    return true;
}

void ActionMenu()
{
    std::cout << "1) MultRow";
    std::cout << std::endl;
    std::cout << "2) MultColumn";
    std::cout << std::endl;
    std::cout << "3) AddRowByRow";
    std::cout << std::endl;
    std::cout << "4) AddColumnByColumn";
    std::cout << std::endl;
    std::cout << "5) SwapRows";
    std::cout << std::endl;
    std::cout << "6) SwapColumns";
    std::cout << std::endl;
    std::cout << "7) MultScalar";
    std::cout << std::endl;
    std::cout << "8) GetNorm";
    std::cout << std::endl;
    std::cout << "9) Print";
    std::cout << std::endl;
    std::cout << "10) Back";
    std::cout << std::endl;

}


void MatrixMenu()
{
    int n, m, placeFrom, placeIn, type;
    std::cout << "Select the matrix type: 1) SquareMatrix 2) RectangularMatrix" << std::endl;
    std::cin >> type;
    int** a = nullptr;
    ArrayRectangularMatrix<int>* matrix = nullptr;
    if (type == 2)
    {
        std::cout << "Input size of matrix (row and column)" << std::endl;
        std::cin >> n >> m;
        a = new int* [n];
        std::cout << "Input elements of matrix" << std::endl;

        for (int i = 0; i < n; i++)
        {
            a[i] = new int[m];
            for (int j = 0; j < m; j++)
            {
                std::cin >> a[i][j];
            }
        }
        matrix = new ArrayRectangularMatrix<int>(n, m, a);
    }
    if (type == 1)
    {
        std::cout << "Input size of matrix" << std::endl;
        std::cin >> n;
        a = new int* [n];
        m = n;
        std::cout << "Input elements of matrix" << std::endl;
        for (int i = 0; i < n; i++)
        {
            a[i] = new int[n];
            for (int j = 0; j < n; j++)
            {
                std::cin >> a[i][j];
            }
        }
        matrix = new ArraySquareMatrix<int>(n, a);
    }

    for (int i = 0; i < n; i++)
    {
        delete[] a[i];
    }
    delete[] a;
    int action;
    do
    {
        ActionMenu();
        std::cin >> action;
        switch (action)
        {
        case(1):
            std::cout << "Input coefficient and place:";
            std::cout << std::endl;
            std::cin >> n;
            do
            {
                std::cin >> placeIn;
            } while (!IfCorrectIndex(*matrix, placeIn));
            matrix->MultRow(placeIn, n);
            break;
        case(2):
            std::cout << "Input coefficient and place:";
            std::cout << std::endl;
            std::cin >> n;
            do
            {
                std::cin >> placeIn;
            } while (!IfCorrectIndex(*matrix, placeIn));
            matrix->MultColumn(placeIn, n);
            break;
        case(3):
            std::cout << "Input coefficient and place from, in:";
            std::cout << std::endl;
            std::cin >> n;
            do
            {
                std::cin >> placeFrom >> placeIn;
            } while (!IfCorrectIndex(*matrix, placeIn) && !IfCorrectIndex(*matrix, placeFrom));
            matrix->AddRowByRow(placeIn, placeFrom, n);
            break;
        case(4):
            std::cout << "Input coefficient and place from, in:";
            std::cout << std::endl;
            std::cin >> n;
            do
            {
                std::cin >> placeFrom >> placeIn;
            } while (!IfCorrectIndex(*matrix, placeIn) && !IfCorrectIndex(*matrix, placeFrom));
            matrix->AddColumnByColumn(placeIn, placeFrom, n);
            break;
        case(5):
            std::cout << "Input elements of swap(i, j):";
            std::cout << std::endl;
            do
            {
                std::cin >> placeFrom >> placeIn;
            } while (!IfCorrectIndex(*matrix, placeIn) && !IfCorrectIndex(*matrix, placeFrom));
            matrix->SwapRows(placeIn, placeFrom);
            break;
        case(6):
            std::cout << "Input elements of swap(i, j):";
            std::cout << std::endl;
            do
            {
                std::cin >> placeFrom >> placeIn;
            } while (!IfCorrectIndex(*matrix, placeIn) && !IfCorrectIndex(*matrix, placeFrom));
            matrix->SwapColumns(placeIn, placeFrom);
            break;
        case(7):
            std::cout << "Input coefficient:";
            std::cout << std::endl;
            std::cin >> n;
            matrix->MultScalar(n);
            break;
        case(8):
            std::cout << "Array matrix norm:";
            std::cout << std::endl;
            std::cout << matrix->GetNorm();
            std::cout << std::endl;
            break;
        case(9):
            std::cout << "Array matrix";
            std::cout << std::endl;
            std::cout << *matrix;
            std::cout << std::endl;
            break;
        case(10):
            delete matrix;
            break;
        default:
            std::cout << "Repid exit";
            std::cout << std::endl;
            break;
        }
    } while (action != 10);
}
