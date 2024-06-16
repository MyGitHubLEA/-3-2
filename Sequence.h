#pragma once

#include "DynamicArray.h"
#include <iostream>

template <typename T> class Sequence
{
public:
	virtual T GetFirst() const = 0;
	virtual T GetLast() const = 0;
	virtual T Get(int index) const = 0;
	virtual int GetLength() const = 0;
};

template <typename T> class MutableSequence : public Sequence<T>
{
public:
	virtual void Set(const T& item, int index) = 0;
	virtual void Append(const T& item) = 0;
	virtual void Prepend(const T& item) = 0;
	virtual void InsertAt(const T& item, int index) = 0;
	virtual Sequence<T>* GetSubsequence(int start, int end) const = 0;
	virtual Sequence<T>* Concat(const Sequence<T>* seq) const = 0;
};

template <typename T> class ImmutableSequence : public Sequence<T>
{
public:
	virtual ImmutableSequence<T>* Set(const T& item, int index) const = 0;
	virtual ImmutableSequence<T>* Append(const T& item) const = 0;
	virtual ImmutableSequence<T>* Prepend(const T& item) const = 0;
	virtual ImmutableSequence<T>* InsertAt(const T& item, int index) const = 0;
	virtual Sequence<T>* GetSubsequence(int start, int end) const = 0;
	virtual Sequence<T>* Concat(const Sequence<T>* seq) const = 0;
};



template <typename T> class MutableArraySequence : public MutableSequence<T>
{
private:
	DynamicArray<T>* array;
public:

	MutableArraySequence()
	{
		this->array = new DynamicArray<T>(0);
	}

	MutableArraySequence(const MutableArraySequence<T>& seq)
	{
		this->array = new DynamicArray<T>(*seq.array);
	}

	MutableArraySequence(T* items, int size) 
	{
		this->array = new DynamicArray<T>(items, size);
	}

	~MutableArraySequence()
	{
		delete this->array;
	}

	T GetFirst() const override
	{
		if (this->GetLength() == 0) throw std::out_of_range("Sequence is empty\n");
		return this->array->Get(0);
	}

	T GetLast() const override
	{
		if (this->GetLength() == 0) throw std::out_of_range("Sequence is empty\n");
		return this->array->Get(this->GetLength() - 1);
	}

	T Get(int index) const override
	{
		if (index < 0 || index >= this->GetLength())
			throw std::out_of_range("IndexOutOfRange\n");
		return this->array->Get(index);
	}

	int GetLength() const override 
	{
		return this->array->GetLength();
	}

	void Set(const T& item, int index) override
	{
		if (index < 0 || index >= this->GetLength())
			throw std::out_of_range("IndexOutOfRange\n");
		this->array->Set(item, index);
	}

	MutableArraySequence<T>* GetSubsequence(int start, int end) const override
	{
		if (start < 0 || start > this->GetLength() || start > end || end >= this->GetLength())
			throw std::out_of_range("IndexOutOfRange\n");
		int sizeArray = end - start + 1;
		MutableArraySequence<T>* subSequence = new MutableArraySequence<T>();
		subSequence->array->Resize(sizeArray);
		for (int i = 0; i < sizeArray; i++)
			subSequence->Set(this->Get(i + start), i);
		return subSequence;
	}

	void Append(const T& item) override
	{
		this->array->Resize(this->GetLength() + 1);
		this->array->Set(item, this->GetLength() - 1);
	}

	void Prepend(const T& item) override
	{
		T t1;
		int Len = this->GetLength();
		this->array->Resize(this->GetLength() + 1);
		for (int i = 0; i < Len; i++)
		{
			t1 = this->array->Get(Len - (i + 1));
			this->array->Set(t1, Len - i);
		}
		this->array->Set(item, 0);
	}

	void InsertAt(const T& item, int index) override
	{
		if (index < 0 || index >= this->GetLength()) throw std::out_of_range("IndexOutOfRange\n");
		int Len = this->GetLength();
		this->array->Resize(this->GetLength() + 1);
		T t1;
		for (int i = Len - 1; i >= index; --i)
		{
			t1 = this->array->Get(i);
			this->array->Set(t1, i + 1);
		}
		this->array->Set(item, index);
	}

	MutableArraySequence<T>* Concat(const Sequence<T>* seq) const override
	{
		MutableArraySequence<T>* newArray = new MutableArraySequence<T>(*this);
		for (int i = 0; i < seq->GetLength(); ++i)
			newArray->Append(seq->Get(i));
		return newArray;
	}
	const MutableArraySequence<T>& operator = (const MutableArraySequence<T>& seq)
	{
		delete this->array;
		DynamicArray<T>* newArray = new DynamicArray<T>(*(seq.array));
		this->array = newArray;
		return *this;
	}
};


template <typename T> class ImmutableArraySequence : public ImmutableSequence<T> {
private:
	DynamicArray<T>* array;
public:

	ImmutableArraySequence()
	{
		this->array = new DynamicArray<T>(0);
	}

	ImmutableArraySequence(const ImmutableArraySequence<T>& seq)
	{
		this->array = new DynamicArray<T>(*seq.array);
	}

	ImmutableArraySequence(T* items, int size)
	{
		this->array = new DynamicArray<T>(items, size);
	}

	~ImmutableArraySequence()
	{
		delete this->array;
	}

	T GetFirst() const override
	{
		if (this->GetLength() == 0) throw std::out_of_range("Sequence is empty\n");
		return this->array->Get(0);
	}

	T GetLast() const override
	{
		if (this->GetLength() == 0) throw std::out_of_range("Sequence is empty\n");
		return this->array->Get(this->GetLength() - 1);
	}

	T Get(int index) const override
	{
		if (index < 0 || index >= this->GetLength())
			throw std::out_of_range("IndexOutOfRange\n");
		return this->array->Get(index);
	}

	int GetLength() const override
	{
		return this->array->GetLength();
	}

	ImmutableArraySequence<T>* Set(const T& item, int index) const override
	{
		if (index < 0 || index >= this->GetLength())
			throw std::out_of_range("IndexOutOfRange\n");
		ImmutableArraySequence<T>* newArray = new ImmutableArraySequence<T>(*this);
		newArray->array->Set(item, index);
		return newArray;
	}

	ImmutableArraySequence<T>* GetSubsequence(int start, int end) const override
	{
		if (start < 0 || start >= this->GetLength() || start >= end || end >= this->GetLength())
			throw std::out_of_range("IndexOutOfRange\n");
		int sizeArray = end - start + 1;
		ImmutableArraySequence<T>* subSequence = new ImmutableArraySequence<T>();
		subSequence->array->Resize(sizeArray);
		for (int i = 0; i < sizeArray; i++)
			subSequence->array->Set(this->Get(i + start), i);
		return subSequence;
	}

	ImmutableArraySequence<T>* Append(const T& item) const override
	{
		ImmutableArraySequence<T>* newArray = new ImmutableArraySequence<T>(*this);
		newArray->array->Resize(this->GetLength() + 1);
		newArray->array->Set(item, this->GetLength());
		return newArray;
	}

	ImmutableArraySequence<T>* Prepend(const T& item) const override
	{
		T t1;
		int Len = this->GetLength();
		ImmutableArraySequence<T>* newArray = new ImmutableArraySequence<T>(*this);
		newArray->array->Resize(this->GetLength() + 1);
		for (int i = 0; i < Len; i++)
		{
			t1 = newArray->array->Get(Len - (i + 1));
			newArray->array->Set(t1, Len - i);
		}
		newArray->array->Set(item, 0);
		return newArray;
	}

	ImmutableArraySequence<T>* InsertAt(const T& item, int index) const override
	{
		if (index < 0 || index >= this->GetLength()) throw std::out_of_range("IndexOutOfRange\n");
		ImmutableArraySequence<T>* newArray = new ImmutableArraySequence<T>(*this);
		int Len = newArray->GetLength();
		newArray->array->Resize(newArray->GetLength() + 1);
		T t1;
		for (int i = Len - 1; i >= index; --i)
		{
			t1 = newArray->array->Get(i);
			newArray->array->Set(t1, i + 1);
		}
		newArray->array->Set(item, index);
		return newArray;
	}

	ImmutableArraySequence<T>* Concat(const Sequence<T>* seq) const override
	{
		ImmutableArraySequence<T>* newArray = new ImmutableArraySequence<T>(*this);
		newArray->array->Resize(this->GetLength() + seq->GetLength());
		for (int i = 0; i < seq->GetLength(); ++i)
			newArray->array->Set(seq->Get(i), this->GetLength() + i);
		return newArray;
	}

};

template <typename T>
std::ostream& operator <<(std::ostream& out, const Sequence<T>& Seq)
{
	for (int i = 0; i < Seq.GetLength(); ++i)
	{
		out << Seq.Get(i) << " ";
	}
	out << "" << std::endl;
	return out;
}
