#pragma once
#include <iostream>
#include "Sequence.h"
#include <stdexcept>
#include <math.h>

template<typename T>
class ArraySquareMatrix
{
private:
    MutableArraySequence<T>* items;
    //последовательность последовательностей. Тест производительности (1000 на 1000) транспонировать
    // если разницы нет, то миллион на миллион. График зависимости времени исполнения от количества
public:
    ArraySquareMatrix(int dimension, T** items);
    ArraySquareMatrix();
    ArraySquareMatrix(const ArraySquareMatrix<T>& squarearMatrix);
    ~ArraySquareMatrix();

    int GetDimension() const; 
    T Get(int indexRow, int indexColumn) const;


    void MultRow(int numberRow, T scalar); 
    void MultColumn(int numberColumn, T scalar); 
    void AddRowByRow(int indexRowWhereAdd, int indexRowWhicheAdd, T coefficient); 
    void AddColumnByColumn(int indexColumnWhereAdd, int indexColumnWhicheAdd, T coefficient); 
    void SwapRows(int indexFirstRow, int indexSecondRow); 
    void SwapColumns(int indexFirstColumn, int indexSecondColumn); 

    void MultScalar(T scalar); 
    void AddMatrix(ArraySquareMatrix<T>* squarearMatrix); 

    T GetNorm(); 
    ArraySquareMatrix<T>& operator += (const ArraySquareMatrix<T>& squarearMatrix);
    ArraySquareMatrix<T>& operator = (const ArraySquareMatrix<T>& squarearMatrix); 
   // friend ArraySquareMatrix<T>& operator << (const ArraySquareMatrix<T>& squarearMatrix);
    //тетсты для равно
};

template<typename T>
ArraySquareMatrix<T>::ArraySquareMatrix(int dimension, T** items)
{
    this->items = new MutableArraySequence<T>[dimension];
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            (this->items)[i].Append(items[i][j]);
        }
    }
}

template<typename T>
ArraySquareMatrix<T>::ArraySquareMatrix()
{
    this->items = nullptr;
}

template<typename T>
ArraySquareMatrix<T>::ArraySquareMatrix(const ArraySquareMatrix<T>& squarearMatrix)
{
    this->items = new MutableArraySequence<T>[squarearMatrix.GetDimension()];

    for (int i = 0; i < squarearMatrix.GetDimension(); i++)
    {
        for (int j = 0; j < squarearMatrix.GetDimension(); j++)
        {
            (this->items)[i].Append(squarearMatrix.Get(i, j));
        }
    }
}


template<typename T>
ArraySquareMatrix<T>& ArraySquareMatrix<T>::operator = (const ArraySquareMatrix<T>& squarearMatrix)
{
    //удалить items
    this->items = new MutableArraySequence<T>[squarearMatrix.GetDimension()];

    for (int i = 0; i < squarearMatrix.GetDimension(); i++)
    {
        for (int j = 0; j < squarearMatrix.GetDimension(); j++)
        {
            (this->items)[i].Append(squarearMatrix.Get(i, j));
        }
    }
    return *this;
}


template<typename T>
ArraySquareMatrix<T>::~ArraySquareMatrix()
{
    delete[] this->items;
}


template<typename T>
int ArraySquareMatrix<T>::GetDimension() const
{
    if (!(this->items))
    {
        return 0;
    }

    return (this->items)[0].GetLength();
}

template<typename T>
T ArraySquareMatrix<T>::Get(int indexRow, int indexColumn) const
{
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexRow >= this->GetDimension() || indexRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexColumn >= this->GetDimension() || indexColumn < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    return (this->items)[indexRow].Get(indexColumn);
}

template<typename T>
void ArraySquareMatrix<T>::MultRow(int numberRow, T scalar) {
    if (!(this->items)) 
    {
        throw std::domain_error("Empty matrix");
    }

    if (numberRow > this->GetDimension() || numberRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(scalar)) {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetDimension(); i++) {
        this->items[numberRow].Set(scalar * this->Get(numberRow, i), i);
    }
}

template<typename T>
void ArraySquareMatrix<T>::MultColumn(int numberColumn, T scalar) {
    if (!(this->items)) {
        throw std::domain_error("Empty matrix");
    }

    if (numberColumn > this->GetDimension() || numberColumn < 0) {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(scalar)) 
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetDimension(); i++) 
    {
        this->items[i].Set(scalar * this->Get(i, numberColumn), numberColumn);
    }
}

template<typename T>
void ArraySquareMatrix<T>::AddRowByRow(int indexRowWhereAdd, int indexRowWhicheAdd, T coefficient) 
{
    if (!(this->items)) 
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexRowWhereAdd >= this->GetDimension() || indexRowWhereAdd < 0) 
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexRowWhicheAdd >= this->GetDimension() || indexRowWhicheAdd < 0) 
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(coefficient)) {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetDimension(); i++)
    {
        this->items[indexRowWhereAdd].Set(coefficient * this->Get(indexRowWhicheAdd, i) + this->Get(indexRowWhereAdd, i), i);
    }
}

template<typename T>
void ArraySquareMatrix<T>::AddColumnByColumn(int indexColumnWhereAdd, int indexColumnWhichAdd, T coefficient) 
{
    if (!(this->items)) 
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexColumnWhereAdd >= this->GetDimension() || indexColumnWhereAdd < 0) 
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexColumnWhichAdd >= this->GetDimension() || indexColumnWhichAdd < 0) 
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(coefficient)) 
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetDimension(); i++) 
    {
        this->items[i].Set(coefficient * this->Get(i, indexColumnWhichAdd) + this->Get(i, indexColumnWhereAdd), indexColumnWhereAdd);
    }
}

template<typename T>
void ArraySquareMatrix<T>::SwapRows(int indexFirstRow, int indexSecondRow) 
{
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexFirstRow >= this->GetDimension() || indexFirstRow < 0) 
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexSecondRow >= this->GetDimension() || indexSecondRow < 0) 
    {
        throw std::out_of_range("Out of the range of the array");
    }

    MutableArraySequence<T>* buf = new MutableArraySequence<T>(this->items[indexFirstRow]); 
    
    for (int i = 0; i < this->GetDimension(); i++) 
    {
        this->items[indexFirstRow].Set(this->Get(indexSecondRow, i), i);
    }
    for (int i = 0; i < this->GetDimension(); i++) 
    {
        this->items[indexSecondRow].Set(buf->Get(i), i);
    }
    delete buf;
}

template<typename T>
void ArraySquareMatrix<T>::SwapColumns(int indexFirstColumn, int indexSecondColumn)
{
    if (!(this->items)) 
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexFirstColumn >= this->GetDimension() || indexFirstColumn < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexSecondColumn >= this->GetDimension() || indexSecondColumn < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    MutableArraySequence<T> buf;
    for (int i = 0; i < this->GetDimension(); i++) 
    {
        buf.Append(this->Get(i, indexFirstColumn));
    }
    for (int i = 0; i < this->GetDimension(); i++) 
    {
        this->items[i].Set(this->Get(i, indexSecondColumn), indexFirstColumn);
    }
    for (int i = 0; i < this->GetDimension(); i++) 
    {
        this->items[i].Set(buf.Get(i), indexSecondColumn);
    }
}

template<typename T>
void ArraySquareMatrix<T>::MultScalar(T scalar)
{
    if (!(this->items)) 
    {
        throw std::domain_error("Empty matrix");
    }

    if (!(scalar)) 
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetDimension(); i++) 
    {
        for (int j = 0; j < this->GetDimension(); j++) 
        {
            this->items[i].Set(this->Get(i, j) * scalar, j);
        }
    }
}

template<typename T>
void ArraySquareMatrix<T>::AddMatrix(ArraySquareMatrix<T>* squareMatrix) 
{
    if (!(this->items)) 
    {
        throw std::domain_error("Empty matrix");
    }

    if (this->GetDimension() != squareMatrix->GetDimension())
    {
        throw std::invalid_argument("the number of rows in the matrix does not match");
    }

    *this += *squareMatrix;
}

template<typename T>
T ArraySquareMatrix<T>::GetNorm()
{
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    T valueNorm = this->Get(0, 0);
    valueNorm *= valueNorm;

    for (int i = 0; i < this->GetDimension(); i++)
    {
        for (int j = 0; j < this->GetDimension(); j++)
        {
            if (!(i == 0 && j == 0))
            {
                valueNorm += this->Get(i, j) * this->Get(i, j);
            }
        }
    }
    return sqrt(valueNorm);
}

template<typename T>
std::ostream& operator << (std::ostream& out, const ArraySquareMatrix<T>& matrix)
{

    for (int i = 0; i < matrix.GetDimension(); ++i)
    {
        for (int j = 0; j < matrix.GetDimension(); ++j)
            out << matrix.Get(i, j) << " ";
        out << std::endl;
    }
    return out;
}


template<typename T>
ArraySquareMatrix<T>& ArraySquareMatrix<T>::operator += (const ArraySquareMatrix<T>& squarearMatrix)
{

    for (int i = 0; i < squarearMatrix.GetDimension(); i++)
    {
        for (int j = 0; j < squarearMatrix.GetDimension(); j++)
        {
            (this->items)[i].Set(squarearMatrix.Get(i, j) + this->Get(i, j), j);
        }
    }
    return *this;
}


template<typename T>
ArraySquareMatrix<T> operator + (const ArraySquareMatrix<T>& matrix_a, const ArraySquareMatrix<T>& matrix_b)
{
    ArraySquareMatrix<T> res = matrix_a;
    res += matrix_b;
    return res;
}

