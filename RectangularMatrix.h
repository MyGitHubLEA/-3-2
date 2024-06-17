#pragma once
#include "Sequence.h"
#include <stdexcept>
#include <math.h>
template<typename T>
class RectangularMatrix
{
private:
    int numberRows;
    MutableArraySequence<T>* items;
public:
    RectangularMatrix(int numberRow, int numberColumn, T** items);
    RectangularMatrix();
    RectangularMatrix(const RectangularMatrix<T>& rectangularMatrix);
    ~RectangularMatrix();

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
    void AddMatrix(RectangularMatrix<T>* rectangularMatrix);

    T GetNorm();

    const RectangularMatrix<T>& operator = (const RectangularMatrix<T>& rectangularMatrix);
    RectangularMatrix<T>& operator += (const RectangularMatrix<T>& rectangularMatrix);
    template <typename type>
    friend std::ostream& operator<<(std::ostream& out, const RectangularMatrix<type>& rectanguarMatrix);
};


template<typename T>
RectangularMatrix<T>::RectangularMatrix(int numberRow, int numberColumn, T** items) {
    this->numberRows = numberRow;

    this->items = new MutableArraySequence<T>[numberRow];
    for (int i = 0; i < numberRow; i++) {
        for (int j = 0; j < numberColumn; j++) {
            (this->items)[i].Append(items[i][j]);
        }
    }
}

template<typename T>
RectangularMatrix<T>::RectangularMatrix() {
    this->numberRows = 0;
    this->items = nullptr;
}

template<typename T>
RectangularMatrix<T>::RectangularMatrix(const RectangularMatrix<T>& rectangularMatrix)
{
    this->numberRows = rectangularMatrix.GetNumberRows();
    this->items = new MutableArraySequence<T>[rectangularMatrix.GetNumberRows()];

    for (int i = 0; i < rectangularMatrix.GetNumberRows(); i++) {
        for (int j = 0; j < rectangularMatrix.GetNumberColumns(); j++)
        {
            (this->items)[i].Append(rectangularMatrix.Get(i, j));
        }
    }
}


template<typename T>
const RectangularMatrix<T>& RectangularMatrix<T>::operator = (const RectangularMatrix<T>& rectangularMatrix)
{
    this->numberRows = rectangularMatrix.GetNumberRows();
    this->items = new MutableArraySequence<T>[rectangularMatrix.GetNumberRows()];

    for (int i = 0; i < rectangularMatrix.GetNumberRows(); i++)
    {
        for (int j = 0; j < rectangularMatrix.GetNumberColumns(); j++)
        {
            (this->items)[i].Append(rectangularMatrix.Get(i, j));
        }
    }

    return *this;
}




template<typename T>
RectangularMatrix<T>::~RectangularMatrix()
{
    delete[] this->items;
}





template<typename T>
int RectangularMatrix<T>::GetNumberRows() const
{
    if (!(this->items))
    {
        return 0;
    }
    return this->numberRows;
}

template<typename T>
int RectangularMatrix<T>::GetNumberColumns() const
{
    if (!(this->items))
    {
        return 0;
    }

    return (this->items)[0].GetLength();
}

template<typename T>
T RectangularMatrix<T>::Get(int indexRow, int indexColumn) const
{
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexRow >= this->numberRows || indexRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexColumn >= this->GetNumberColumns() || indexColumn < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    return (this->items)[indexRow].Get(indexColumn);
}



template<typename T>
void RectangularMatrix<T>::MultRow(int numberRow, T scalar)
{
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    if (numberRow > this->numberRows || numberRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(scalar))
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetNumberColumns(); i++)
    {
        this->items[numberRow].Set(scalar * this->Get(numberRow, i), i);
    }
}

template<typename T>
void RectangularMatrix<T>::MultColumn(int numberColumn, T scalar)
{
    if (!(this->items))
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

    for (int i = 0; i < this->numberRows; i++)
    {
        this->items[i].Set(scalar * this->Get(i, numberColumn), numberColumn);
    }
}

template<typename T>
void RectangularMatrix<T>::AddRowByRow(int indexRowWhereAdd, int indexRowWhicheAdd, T coefficient)
{
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexRowWhereAdd >= this->numberRows || indexRowWhereAdd < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexRowWhicheAdd >= this->numberRows || indexRowWhicheAdd < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(coefficient))
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->GetNumberColumns(); i++)
    {
        this->items[indexRowWhereAdd].Set(coefficient * this->Get(indexRowWhicheAdd, i) + this->Get(indexRowWhereAdd, i), i);
    }
}

template<typename T>
void RectangularMatrix<T>::AddColumnByColumn(int indexColumnWhereAdd, int indexColumnWhicheAdd, T coefficient) {
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexColumnWhereAdd >= this->GetNumberColumns() || indexColumnWhereAdd < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexColumnWhicheAdd >= this->GetNumberColumns() || indexColumnWhicheAdd < 0) {
        throw std::out_of_range("Out of the range of the array");
    }

    if (!(coefficient))
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->numberRows; i++)
    {
        this->items[i].Set(coefficient * this->Get(i, indexColumnWhicheAdd) + this->Get(i, indexColumnWhereAdd), indexColumnWhereAdd);
    }
} 

template<typename T>
void RectangularMatrix<T>::SwapRows(int indexFirstRow, int indexSecondRow)
{
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexFirstRow >= this->numberRows || indexFirstRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexSecondRow >= this->numberRows || indexSecondRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    MutableArraySequence<T> buf;
    for (int i = 0; i < this->GetNumberColumns(); i++)
    {
        buf.Append(this->Get(indexFirstRow, i));
    }
    for (int i = 0; i < this->GetNumberColumns(); i++)
    {
        this->items[indexFirstRow].Set(this->Get(indexSecondRow, i), i);
    }
    for (int i = 0; i < this->GetNumberColumns(); i++)
    {
        this->items[indexSecondRow].Set(buf.Get(i), i);
    }
}

template<typename T>
void RectangularMatrix<T>::SwapColumns(int indexFirstColumn, int indexSecondColumn)
{
    if (!(this->items))
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
    for (int i = 0; i < this->numberRows; i++)
    {
        buf.Append(this->Get(i, indexFirstColumn));
    }
    for (int i = 0; i < this->numberRows; i++)
    {
        this->items[i].Set(this->Get(i, indexSecondColumn), indexFirstColumn);
    }
    for (int i = 0; i < this->numberRows; i++)
    {
        this->items[i].Set(buf.Get(i), indexSecondColumn);
    }
}

template<typename T>
void RectangularMatrix<T>::MultScalar(T scalar)
{
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    if (!(scalar))
    {
        throw std::invalid_argument("invalid value");
    }

    for (int i = 0; i < this->numberRows; i++)
    {
        for (int j = 0; j < this->GetNumberColumns(); j++)
        {
            this->items[i].Set(this->Get(i, j) * scalar, j);
        }
    }
}

template<typename T>
void RectangularMatrix<T>::AddMatrix(RectangularMatrix<T>* rectangularMatrix)
{
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    if (this->numberRows != rectangularMatrix->numberRows)
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
T RectangularMatrix<T>::GetNorm()
{
    if (!(this->items))
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
void RectangularMatrix<T>::Set(T val, int indexRow, int indexColumn)
{
    if (!(this->items))
    {
        throw std::domain_error("Empty matrix");
    }

    if (indexRow >= this->numberRows || indexRow < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }

    if (indexColumn >= this->GetNumberColumns() || indexColumn < 0)
    {
        throw std::out_of_range("Out of the range of the array");
    }
    (this->items)[indexRow].Set(val, indexColumn);
}



template<typename T>
RectangularMatrix<T>& RectangularMatrix<T>::operator += (const RectangularMatrix<T>& rectangularMatrix)
{

    for (int i = 0; i < rectangularMatrix.GetNumberRows(); ++i)
    {
        for (int j = 0; j < rectangularMatrix.GetNumberColumns(); ++j)
        {
            (this->items)[i].Set(rectangularMatrix.Get(i, j) + this->Get(i, j), j);
        }
    }
    return *this;
}


template<typename T>
RectangularMatrix<T> operator + (const RectangularMatrix<T>& matrix_a, const RectangularMatrix<T>& matrix_b)
{
    RectangularMatrix<T> res = matrix_a;
    res += matrix_b;
    return res;
}
