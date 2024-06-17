#pragma once
#include <iostream>
#include "RectangularMatrix.h"
#include "ArrayRectangularMatrix.h"


template<typename T>
std::ostream& operator << (std::ostream& out, const ArrayRectangularMatrix<T>& matrix)
{

    for (int i = 0; i < matrix.GetNumberRows(); ++i)
    {
        for (int j = 0; j < matrix.GetNumberColumns(); ++j)
            out << matrix.Get(i, j) << " ";
        out << std::endl;
    }
    return out;
}


template<typename T>
std::ostream& operator << (std::ostream& out, const RectangularMatrix<T>& matrix)
{

    for (int i = 0; i < matrix.GetNumberRows(); ++i)
    {
        for (int j = 0; j < matrix.GetNumberColumns(); ++j)
            out << matrix.Get(i, j) << " ";
        out << std::endl;
    }
    return out;
}
