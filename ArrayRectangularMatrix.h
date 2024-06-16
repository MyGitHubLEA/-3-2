#pragma once
#include "Sequence.h"
#include <stdexcept>
#include <math.h>
template<typename T>
class ArrayRectangularMatrix 
{
private:
    MutableArraySequence< MutableArraySequence<T>> items;
public:
    ArrayRectangularMatrix(int numberRow, int numberColumn, T** items);
    ArrayRectangularMatrix();
    ArrayRectangularMatrix(const ArrayRectangularMatrix<T>& rectangularMatrix);

    int GetNumberRows() const;
    int GetNumberColumns() const;
    T Get(int indexRow, int indexColumn) const;
    void Set(T val, int indexRow, int indexColumn);

    void MultRow(int numberRow, T scalar);
    void MultColumn(int numberColumn, T scalar);
    void AddRowByRow(int indexRowWhereAdd, int indexRowWhicheAdd, T coefficient);
    void AddColumnByColumn(int indexColumnWhereAdd, int indexColumnWhicheAdd, T coefficient);
    void SwapRows(int indexFirstRow, int indexSecondRow);
    void SwapColumns(int indexFirstColumn, int indexSecondColumn);

    void MultScalar(T scalar);
    void AddMatrix(ArrayRectangularMatrix<T>* rectangularMatrix);

    T GetNorm();

    const ArrayRectangularMatrix<T>& operator = (const ArrayRectangularMatrix<T>& rectangularMatrix);
    ArrayRectangularMatrix<T>& operator += (const ArrayRectangularMatrix<T>& rectangularMatrix);
    template <typename type>
    friend std::ostream& operator<<(std::ostream& out, const ArrayRectangularMatrix<type>& rectanguarMatrix);
};


template<typename T>
ArrayRectangularMatrix<T>::ArrayRectangularMatrix(int numberRow, int numberColumn, T** items) {
    MutableArraySequence<T>* M = new MutableArraySequence<T>[numberRow]();
    for (int i = 0; i < numberRow; i++)
    {
        for (int j = 0; j < numberColumn; j++)
        {
            M[i].Append(items[i][j]);
        }
    }
   
    this->items = MutableArraySequence<MutableArraySequence<T>> (M, numberRow);

    
}

template<typename T>
ArrayRectangularMatrix<T>::ArrayRectangularMatrix() {
   this->items = MutableArraySequence<MutableArraySequence<T>>();
}

template<typename T>
ArrayRectangularMatrix<T>::ArrayRectangularMatrix(const ArrayRectangularMatrix<T>& rectangularMatrix) 
{
    this->items = MutableArraySequence<MutableArraySequence<T>>(rectangularMatrix.items);
}


template<typename T>
const ArrayRectangularMatrix<T>& ArrayRectangularMatrix<T>::operator = (const ArrayRectangularMatrix<T>& rectangularMatrix) 
{
    this->items = MutableArraySequence<MutableArraySequence<T>>(rectangularMatrix->items);
    return *this;
}


template<typename T>
int ArrayRectangularMatrix<T>::GetNumberRows() const 
{
    if ((this->items).GetLength()==0)
    {
        return 0;
    }
    return this->items.GetLength();
}

template<typename T>
int ArrayRectangularMatrix<T>::GetNumberColumns() const 
{
   
    return (this->items).Get(0).GetLength();
}

template<typename T>
T ArrayRectangularMatrix<T>::Get(int indexRow, int indexColumn) const 
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexRow >= this->GetNumberRows() || indexRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexColumn >= this->GetNumberColumns() || indexColumn < 0) 
    {
       throw std::out_of_range("Out of the range of the array");
    }
    return ((this->items).Get(indexRow)).Get(indexColumn);
}

template<typename T>
void ArrayRectangularMatrix<T>::MultRow(int numberRow, T scalar)
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }
    
    if (numberRow > this->GetNumberRows() || numberRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(scalar))
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetNumberColumns(); i++) 
    {

        this->Set(scalar * (this->Get(numberRow, i)), numberRow, i);
    }
}

template<typename T>
void ArrayRectangularMatrix<T>::MultColumn(int numberColumn, T scalar) 
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }
    
    if (numberColumn > this->GetNumberColumns() || numberColumn < 0) 
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(scalar))
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetNumberRows(); i++)
    {
        this->Set(scalar * this->Get(i, numberColumn), i, numberColumn);
    }
}

template<typename T>
void ArrayRectangularMatrix<T>::AddRowByRow(int indexRowWhereAdd, int indexRowWhicheAdd, T coefficient) 
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }
    
    if (indexRowWhereAdd >= this->GetNumberRows() || indexRowWhereAdd < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexRowWhicheAdd >= this->GetNumberColumns() || indexRowWhicheAdd < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(coefficient)) 
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetNumberColumns(); i++) 
    {
        this->Set(coefficient * this->Get(indexRowWhicheAdd, i) + this->Get(indexRowWhereAdd, i),indexRowWhereAdd, i);
    }
}

template<typename T>
void ArrayRectangularMatrix<T>::AddColumnByColumn(int indexColumnWhereAdd, int indexColumnWhichAdd, T coefficient) 
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }
    
    if (indexColumnWhereAdd >= this->GetNumberColumns() || indexColumnWhereAdd < 0) 
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexColumnWhichAdd >= this->GetNumberColumns() || indexColumnWhichAdd < 0) {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(coefficient)) 
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetNumberRows(); i++)
    {
        this->Set(coefficient * this->Get(i, indexColumnWhichAdd) + this->Get(i, indexColumnWhereAdd), i, indexColumnWhereAdd);
    }
}

template<typename T>
void ArrayRectangularMatrix<T>::SwapRows(int indexFirstRow, int indexSecondRow)
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexFirstRow >= this->GetNumberRows() || indexFirstRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexSecondRow >= this->GetNumberRows() || indexSecondRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    MutableArraySequence<T> buf;
    for (int i = 0; i < this->GetNumberColumns(); i++) 
    {
        buf.Append(this->Get(indexFirstRow, i));
    }
    
    this->items.Set(this->items.Get(indexSecondRow), indexFirstRow);
    this->items.Set(buf, indexSecondRow);
}

template<typename T>
void ArrayRectangularMatrix<T>::SwapColumns(int indexFirstColumn, int indexSecondColumn) 
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexFirstColumn >= this->GetNumberColumns() || indexFirstColumn < 0) 
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexSecondColumn >= this->GetNumberColumns() || indexSecondColumn < 0) 
    {
        throw std::out_of_range("Out of the range of the array");
    }

    MutableArraySequence<T> buf;
    for (int i = 0; i < this->GetNumberRows(); i++)
    {
        buf.Append(this->Get(i, indexFirstColumn));
    }
    for (int i = 0; i < this->GetNumberRows(); i++)
    {
        this->Set(this->Get(i, indexSecondColumn), i, indexFirstColumn);
    }
    for (int i = 0; i < this->GetNumberRows(); i++)
    {
        this->Set(buf.Get(i), i, indexSecondColumn);
    }
}

template<typename T>
void ArrayRectangularMatrix<T>::MultScalar(T scalar) 
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }

    if (!(scalar))
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetNumberRows(); i++)
    {
        for (int j = 0; j < this->GetNumberColumns(); j++)
        {
            this->Set(this->Get(i, j) * scalar, i, j);
        }
    }
}


template<typename T>
void ArrayRectangularMatrix<T>::Set(T val, int indexRow, int indexColumn) 
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }
    
    if (indexRow >= this->GetNumberRows() || indexRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexColumn >= this->GetNumberColumns() || indexColumn < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }
    MutableArraySequence<T> newArray(this->items.Get(indexRow));
    newArray.Set(val, indexColumn);
    (this->items).Set(newArray, indexRow);
}


template<typename T>
void ArrayRectangularMatrix<T>::AddMatrix(ArrayRectangularMatrix<T>* rectangularMatrix) 
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }
    
    if (this->GetNumberRows() != rectangularMatrix->GetNumberRows())
    {
        throw std::invalid_argument("the number of rows in the matrix does not match");
    }

    if (this->GetNumberColumns() != rectangularMatrix->GetNumberColumns()) 
    {
        throw std::invalid_argument("the number of columns in the matrix does not match");
    }
    *this += *rectangularMatrix;
}

template<typename T>
T ArrayRectangularMatrix<T>::GetNorm() 
{
    if ((this->items).GetLength() == 0)
    {
        throw std::domain_error("Empty matrix");
    }
    
    T valueNorm = this->Get(0, 0);
    valueNorm *= valueNorm;

    for (int i = 0; i < this->GetNumberRows(); i++) {
        for (int j = 0; j < this->GetNumberColumns(); j++) {
            if (!(i == 0 && j == 0)) {
                valueNorm += this->Get(i, j) * this->Get(i, j);
            }
        }
    }

    return sqrt(valueNorm);
}



template<typename T>
ArrayRectangularMatrix<T>& ArrayRectangularMatrix<T>::operator += (const ArrayRectangularMatrix<T>& rectangularMatrix)
{

    for (int i = 0; i < rectangularMatrix.GetNumberRows(); ++i)
    {
        for (int j = 0; j < rectangularMatrix.GetNumberColumns(); ++j)
        {
            this->Set(rectangularMatrix.Get(i, j) + this->Get(i, j), i, j);
        }
    }
    return *this;
}


template<typename T>
ArrayRectangularMatrix<T> operator + (const ArrayRectangularMatrix<T>& matrix_a, const ArrayRectangularMatrix<T>& matrix_b)
{
    ArrayRectangularMatrix<T> res = matrix_a;
    res += matrix_b;
    return res;
}
