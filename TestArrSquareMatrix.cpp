#include "SquareMatrix.h"
#include <cassert>
//------------можно раскомментировать, скомпилируется-------------//
//#include "UI.h"
//#include "UI.h"
int** a;
const double E = 5.0e-10;

void SInisialisation()
{
    a = new int* [4];
    a[0] = new int[4] {0, 1, 2, 3};
    a[1] = new int[4] {1, 2, 3, 4};
    a[2] = new int[4] {2, 3, 4, 5};
    a[3] = new int[4] {3, 4, 5, 6};

}


void SDelete()
{
    for (int i = 0; i < 4; i++)
    {
        delete[] a[i];
    }
    delete[] a;

}


void SMatrixTestGetDimension() 
{
    ArraySquareMatrix<int> matrix_a (4, a);
    ArraySquareMatrix<int> matrix_b (matrix_a);
    ArraySquareMatrix<int> matrix_c;

    assert(matrix_a.GetDimension() == 4);
    assert(matrix_b.GetDimension() == 4);
    assert(matrix_c.GetDimension() == 0);
    // утечки памяти везде, сделать переменные матрицы, завести две статические переменные, инициализатор, 
    //чтобы там явно была видна матрица
    std::cout << "GetDimension tests passed" << std::endl;
}

void SMatrixTestGet()
{
    
    ArraySquareMatrix<int> matrix(4, a);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            assert(matrix.Get(i, j) == a [i][j]); 
    std::cout << "Get tests passed" << std::endl;
}

void SMatrixTestMultScalar() 
{
   
    ArraySquareMatrix<int> matrix (4, a);
    matrix.MultScalar(2);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            assert(matrix.Get(i, j) == 2 * a[i][j]);
    std::cout << "MultScalar tests passed" << std::endl;
}

void SMatrixTestGetNorm()
{
   
    ArraySquareMatrix<int> matrix(4, a);
    double norm = matrix.GetNorm();
    assert( 13 + E > norm && 13 - E < norm);
    std::cout << "GetNorm test passed" << std::endl;
}

void SMatrixTestMultColumn()
{
    
    ArraySquareMatrix<int> matrix (4, a);
    matrix.MultColumn(2, 2);
    for (int i = 0; i < 4; ++i)
       assert(matrix.Get(i, 2) == 2 * a[i][2]);
    std::cout << "MultColumn test passed" << std::endl;
}

void SMatrixTestMultRow()
{
    
    ArraySquareMatrix<int> matrix (4, a);
    matrix.MultRow(2, 2);
    for (int j = 0; j < 4; ++j)
        assert(matrix.Get(2, j) == 2 * a[2][j]);
    std::cout << "MultRow test passed" << std::endl;
}

void SMatrixTestAddMatrix()
{
  
    ArraySquareMatrix<int> matrix_a (4, a);
    ArraySquareMatrix<int> matrix_b (matrix_a);
    matrix_a.AddMatrix(&matrix_b);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            assert(matrix_a.Get(i, j) == 2 * a[i][j]);
    std::cout << "AddMatrix tests passed" << std::endl;
}

void SMatrixTestAddColumnByColumn()
{
    
    ArraySquareMatrix<int> matrix (4, a);
    matrix.AddColumnByColumn(0, 1, 3);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
            if (j == 0)
            assert(matrix.Get(i, j) == a[i][j] + 3 * a[i][1]);
            else 
            assert(matrix.Get(i, j) == a[i] [j]);
        }
    std::cout << "AddColumnByColumn test passed" << std::endl;
}

void SMatrixTestAddRowByRow()
{
    ArraySquareMatrix<int>* matrix = new ArraySquareMatrix<int>(4, a);
    matrix->AddRowByRow(0, 1, 3);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
            if (i == 0) 
            assert(matrix->Get(i, j) == a[i][j] + 3 * a[1][j]);
            else 
            assert(matrix->Get(i, j) ==  a[i][j]);
        }
    std::cout << "AddRowByRow test passed" << std::endl;
}

void SMatrixTestSwapRows()
{

    ArraySquareMatrix<int>* matrix = new ArraySquareMatrix<int>(4, a);
    matrix->SwapRows(1, 2);
    for (int i = 0; i < 4; ++i)
    {
        assert(matrix->Get(1, i) == a[2][i]);
        assert(matrix->Get(2, i) == a[1][i]);
    }
    std::cout << "SwapRows test passed" << std::endl;
}

void SMatrixTestSwapColumns()
{
   

    ArraySquareMatrix<int>* matrix = new ArraySquareMatrix<int>(4, a);
    matrix->SwapColumns(1, 2);
    for (int i = 0; i < 4; ++i)
    {
        assert(matrix->Get(i, 1) == a[2][i]);
        assert(matrix->Get(i, 2) == a[1][i]);
    }
    std::cout << "SwapColumns test passed" << std::endl;
}


void TestarSquareMatrix()
{
    SInisialisation();
    SMatrixTestGetDimension();
    SMatrixTestGet();
    SMatrixTestMultScalar();
    SMatrixTestGetNorm();
    SMatrixTestMultColumn();
    SMatrixTestMultRow();
    SMatrixTestAddMatrix();
    SMatrixTestAddColumnByColumn();
    SMatrixTestAddRowByRow();
    SMatrixTestSwapRows();
    SMatrixTestSwapColumns();
    SDelete();
}