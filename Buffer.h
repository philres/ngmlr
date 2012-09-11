#ifndef __BUFFER_H__
#define __BUFFER_H__

/*
Synchronized Buffer

The buffer can lead to deadlocks if capacity is less than the sum
of the highest values for count in simultaneous Read and Write calls.
(i.e., if count(Write)+count(Read) > capacity at any given time)

This could be prevented by allowing write to write only part of its batch
instead of all-or-nothing, this can however lead to starvation of single
threads.
*/

template <typename T> class Buffer
{
private:
	template <typename U> class _Buffer;
	_Buffer<T> * m_pImpl;
	char const * const m_Name;
public:
	Buffer(char const * const name, int Capacity, int * pCountReport);
	~Buffer();
	int Capacity() const;
	int Count() const;
	int Free() const;
	
	// Fillrate of the buffer
	float Load() const;

	void Register() const;
	void Release() const;

	// Blocking write count elements from data to buffer
	bool Write(T * data, int count, bool block = true);
	bool WriteR(T pData, int count, bool block = true);
	// Non-blocking write count elements from data to buffer
	bool TryWrite(T * data, int count);
	// Blocking read count elements from buffer into data
	int  Read (T * data, int count) const;
};

#endif
