#include "ArrayRectangularMatrix.h"
#include <cassert>

int** arr;
const double E = 5.0e-10;

void Inisialisation () 
{
    arr= new int* [4];
    arr[0] = new int[5] {0, 1, 2, 3, 4};
    arr[1] = new int[5] {1, 2, 3, 4, 5};
    arr[2] = new int[5] {2, 3, 4, 5, 6};
    arr[3] = new int[5] {3, 4, 5, 6, 7};

}

void Delete() 
{
    for (int i = 0; i < 4; i++)
    {
        delete[] arr[i];
    }
    delete[] arr;

}

void ReMatrixTestGetNumberColumnsRows()
{
   
    ArrayRectangularMatrix<int> matrix_a (4, 5, arr);
    ArrayRectangularMatrix<int> matrix_b (matrix_a);
   // ArrayRectangularMatrix<int> matrix_c;

    assert(matrix_a.GetNumberRows() == 4);
    assert(matrix_b.GetNumberRows() == 4);
    //assert(matrix_c.GetNumberRows() == 0);
    assert(matrix_a.GetNumberColumns() == 5);
    assert(matrix_b.GetNumberColumns() == 5);
    //assert(matrix_c.GetNumberColumns() == 0);

    std::cout << "GetNumberRows and GetNumberColumns tests passed" << std::endl;
}

void ReMatrixTestGet()
{
   
    ArrayRectangularMatrix<int> matrix (4, 5, arr);
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
           assert(matrix.Get(i, j) == arr[i][j]);
        }
    }
    std::cout << "Get tests passed" << std::endl;
}

void ReMatrixTestMultScalar()
{
   
    ArrayRectangularMatrix<int> matrix (4, 5, arr);
    matrix.MultScalar(2);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 5; ++j)
            //std::cout << matrix.Get(i, j) << " " << arr[i][j] << std::endl;
            assert(matrix.Get(i, j) == 2 * arr[i][j]);
    std::cout << "MultScalar tests passed" << std::endl;
}

void ReMatrixTestGetNorm()
{
   
    ArrayRectangularMatrix<int> matrix (4, 5, arr);
    double norm = matrix.GetNorm();
    assert(17 + E > norm && 17 - E < norm);
    std::cout << "GetNorm test passed" << std::endl;
}

void ReMatrixTestMultColumn()
{
    ArrayRectangularMatrix<int> matrix(4, 5, arr);
    matrix.MultColumn(2, 2);
    for (int i = 0; i < 4; ++i)
        assert(matrix.Get(i, 2) == 2 * arr[i][2]);
    std::cout << "MultColumn test passed" << std::endl;
}

void ReMatrixTestMultRow()
{
    ArrayRectangularMatrix<int> matrix (4, 5, arr);
    matrix.MultRow(2, 2);
    for (int j = 0; j < 5; ++j)
        assert(matrix.Get(2, j) == 2 * arr[2][j]);
    std::cout << "MultRow test passed" << std::endl;
}

void ReMatrixTestAddMatrix()
{
    ArrayRectangularMatrix<int> matrix_a (4, 5, arr);
    ArrayRectangularMatrix<int> matrix_b (matrix_a);
    matrix_a.AddMatrix(&matrix_b);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 5; ++j)
            assert(matrix_a.Get(i, j) == 2 * arr[i][j]);
    std::cout << "AddMatrix tests passed" << std::endl;
}

void ReMatrixTestAddColumnByColumn()
{
    ArrayRectangularMatrix<int> matrix (4, 5, arr);
    matrix.AddColumnByColumn(0, 1, 3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
        {
            if (j == 0)
                assert(matrix.Get(i, j) == arr[i][j] + 3 * arr[i][1]);
            else
                assert(matrix.Get(i, j) == arr[i][j]);
        }
    std::cout << "AddColumnByColumn test passed" << std::endl;
}

void ReMatrixTestAddRowByRow()
{
    ArrayRectangularMatrix<int> matrix(4, 5, arr);
    matrix.AddRowByRow(0, 1, 3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
        {
            if (i == 0)
                assert(matrix.Get(i, j) == arr[i][j] + 3 * arr[1][j]);
            else
                assert(matrix.Get(i, j) ==  arr[i][j]);
        }
    std::cout << "AddRowByRow test passed" << std::endl;
}


void ReMatrixTestSwapRows()
{
    ArrayRectangularMatrix<int> matrix (4, 5, arr);
    matrix.SwapRows(1, 2);
    for (int i = 0; i < 4; ++i)
    {
        assert(matrix.Get(1, i) == arr[2][i]);
        assert(matrix.Get(2, i) == arr[1][i]);
    }
    std::cout << "SwapRows test passed" << std::endl;
}

void ReMatrixTestSwapColumns()
{
    ArrayRectangularMatrix<int> matrix (4, 5, arr);

    matrix.SwapColumns(1, 2);
    for (int i = 0; i < 3; ++i)
    {
        assert(matrix.Get(i, 1) == arr[2][i]);
        assert(matrix.Get(i, 2) == arr[1][i]);
    }
    std::cout << "SwapColumns test passed" << std::endl;
}



void TestarRectangularMatrix()
{
    Inisialisation();
    ReMatrixTestGetNumberColumnsRows();
    ReMatrixTestGet();
    ReMatrixTestMultScalar();
    ReMatrixTestGetNorm();
    ReMatrixTestMultColumn();
    ReMatrixTestMultRow();
    ReMatrixTestAddMatrix();
    ReMatrixTestAddColumnByColumn();
    ReMatrixTestAddRowByRow();
    ReMatrixTestSwapRows();
    ReMatrixTestSwapColumns();
    Delete();
}