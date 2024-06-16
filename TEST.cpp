#pragma once
#include "ArrayRectangularMatrix.h"
#include "RectangularMatrix.h"
#include <time.h>
#include <random>

void Calculate(clock_t start, clock_t  end) 
{
	double time_taken = double(end - start);
	std::cout << std::endl << "time: " << time_taken;
	std::cout << std::endl << std::endl;
}

int** Inisialisation(int N) 
{
	int** M;
	M = new int* [N];
	for (int i = 0; i < N; i++)
	{
		M[i] = new int[N];
		for (int j = 0; j < N; j++)
		{
			M[i][j] = rand() % 10;
		}
	}
	return M;
}




void TEST()
{
	clock_t start, end;
	//------------ArrayRectangularMatrix - матрица матриц---------------//
	for (int n = 1; n < 250; n += 10)
	{
		int** arr = Inisialisation(n);
		ArrayRectangularMatrix<int> matrixA(n, n, arr);
		//std::cout << matrixA;
		start = clock();
		ArrayRectangularMatrix<int> matrixTransposed(matrixA);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				matrixTransposed.Set(matrixA.Get(i, j), j, i);
			}
		}
		end = clock();
		std::cout << "Size of matrix: " << n;
		Calculate(start, end);
	}
	//------------RectangularMatrix - массив матриц---------------//
	for (int n = 1; n < 250; n += 10)
	{
		int** arr = Inisialisation(n);
		RectangularMatrix<int> matrixB(n, n, arr);
		start = clock();
		RectangularMatrix<int> matrixTransposed(matrixB);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				matrixTransposed.Set(matrixB.Get(i, j), j, i);
			}
		}
		end = clock();
		std::cout <<"Size of matrix: "<< n;
		Calculate(start, end);
	}
}