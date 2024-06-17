#pragma once
#include "Sequence.h"
#include "ArrayRectangularMatrix.h"

template<typename T>
class ArraySquareMatrix : public ArrayRectangularMatrix<T>
{
public:
    ArraySquareMatrix(int dimension, T** items);
    ArraySquareMatrix();
    ArraySquareMatrix(const ArraySquareMatrix<T>& squarearMatrix);
    int GetDimension();
};

template<typename T>
ArraySquareMatrix<T>::ArraySquareMatrix(int dimension, T** items) : ArrayRectangularMatrix<T>(dimension, dimension, items) {}

template<typename T>
ArraySquareMatrix<T>::ArraySquareMatrix() : ArrayRectangularMatrix<T>() {}

template<typename T>
ArraySquareMatrix<T>::ArraySquareMatrix(const ArraySquareMatrix<T>& squarearMatrix) : ArrayRectangularMatrix<T>(squarearMatrix) {}

template<typename T>
int ArraySquareMatrix<T>::GetDimension()
{
    return this->GetNumberRows();
}
