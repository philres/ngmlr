#include "Buffer.h"
#include "NGM.h"

#include <memory.h>

#include "Debug.h"

#undef module_name
#define module_name "BUFFER"

template <typename T>
template <typename U>
class Buffer<T>::_Buffer
{
private:
	int const m_Capacity;
	volatile int m_Count;
	volatile int * m_pCountReport;
	int m_rPos;
	int m_wPos;

	int m_rReq;

	int m_Writers;

	T * m_pData;

public:
	NGMMutex * mutex;
	NGMThreadWait * r_wait;
	NGMThreadWait * w_wait;

	_Buffer(int capacity, int * pCountReport) :
		m_Capacity(capacity),
		m_Count(0), 
		m_rPos(0), 
		m_wPos(0), 
		m_rReq(0),
		m_Writers(0),
		m_pData(new T[capacity]),
		mutex(new NGMMutex()),
		r_wait(new NGMThreadWait()),
		w_wait(new NGMThreadWait())
	{
		if (pCountReport)
			m_pCountReport = pCountReport;
		else
			m_pCountReport = new int();
		NGMInitMutex(mutex);
		NGMInitWait(r_wait);
		NGMInitWait(w_wait);
		Log.Verbose("Buffer internal ctor called");
	}
	~_Buffer()
	{
		delete[] m_pData;
		delete mutex;
		delete r_wait;
		delete w_wait;
	}

	inline T * & Data() { return m_pData; }
	inline int const Capacity() { return m_Capacity; }
	inline int const Count() { return m_Count; }
	inline int const RPos() { return m_rPos; }
	inline int const WPos() { return m_wPos; }
	inline int const RReq() { return m_rReq; }
	inline void ReqR(int const req) { m_rReq = req; }
	inline void ReqW(int const req) {}
	inline int const Free() { return m_Capacity - m_Count; }
	inline void IncR() { ++m_rPos; if (m_rPos >= m_Capacity) m_rPos = 0; }
	inline void IncW() { ++m_wPos; if (m_wPos >= m_Capacity) m_wPos = 0; }
	inline int ModCount(int n) { return *m_pCountReport = (m_Count += n); }

	void Register()
	{ 
		NGMLock(mutex); 
		++m_Writers; 
		NGMUnlock(mutex); 
	}
	void Release() 
	{ 
		NGMLock(mutex);
		--m_Writers; 
		if (m_Writers == 0) 
			NGMSignal(r_wait); 
		NGMUnlock(mutex); 
	}
	int GetWriters() { return m_Writers; }
};


template <typename T> Buffer<T>::Buffer(char const * const name, int capacity, int * pCountReport) : m_Name(name)
{
	Log.Verbose("Creating %i element buffer %s", capacity, name);
	m_pImpl = new _Buffer<T>(capacity, pCountReport);
}
template <typename T> Buffer<T>::~Buffer()
{
	delete m_pImpl;
}

template <typename T> int Buffer<T>::Capacity() const
{
	return m_pImpl->Capacity();
}

template <typename T> int Buffer<T>::Count() const
{
	return m_pImpl->Count();
}

template <typename T> int Buffer<T>::Free() const
{
	return Capacity() - Count();
}

template <typename T> float Buffer<T>::Load() const
{
	return (float)Count() / (float)Capacity();
}

template <typename T> int Buffer<T>::Read(T * data, int count) const
{
	NGMLock(m_pImpl->mutex);
	int req = count;

	if (m_pImpl->GetWriters() == 0)
	{
		if (m_pImpl->Count() == 0)
		{
			NGMUnlock(m_pImpl->mutex);
			return 0;
		}
		else
		{
			// Last (incomplete) batch
			count = m_pImpl->Count();
		}
	}
	while (m_pImpl->Count() < count)
	{
		m_pImpl->ReqR(count);
		Log.Verbose("<%s> waiting for read (req = %i)", m_Name, count);
		NGMWait(m_pImpl->mutex, m_pImpl->r_wait);
		Log.Verbose("<%s> woke up for read", m_Name);

		if (m_pImpl->GetWriters() == 0)
			count = m_pImpl->Count();
	}
	Log.Verbose("<%s> starting to read (c = %i)", m_Name, count);

	if (count > req)
		count = req;

	for (int i = 0; i < count; ++i)
	{
		data[i] = m_pImpl->Data()[m_pImpl->RPos()];
		//m_pImpl->Data()[m_pImpl->RPos()] = 0;
		m_pImpl->IncR();
	}

	m_pImpl->ModCount(-count);
	Log.Verbose("<%s> mod by -%i (req %i) (now %i)", m_Name, count, req, m_pImpl->Count());

	NGMSignal(m_pImpl->w_wait);

	NGMUnlock(m_pImpl->mutex);
	return count;
}

template <typename T> bool Buffer<T>::TryWrite(T * pData, int count)
{
	return Write(pData, count, false);
}

template <typename T> bool Buffer<T>::WriteR(T pData, int count, bool block)
{
	if ( count > m_pImpl->Capacity() || count < 0)
	{
		Log.Error("Illegal write to buffer <%s> (write count %i, capacity %i)", m_Name, count, m_pImpl->Capacity());
		Fatal();
	}
	NGMLock(m_pImpl->mutex);

	while (m_pImpl->Free() < count)
	{
		Log.Verbose("<%s> Buffer full - waiting", m_Name);

		if (block)
		{
			NGMWait(m_pImpl->mutex, m_pImpl->w_wait);
		}
		else
		{
			NGMUnlock(m_pImpl->mutex);
			return false;
		}
	}
	//Log.Message("Writing %i elements (%i free)", count, m_pImpl->Free());

	for (int i = 0; i < count; ++i)
	{
		m_pImpl->Data()[m_pImpl->WPos()] = &pData[i];
		m_pImpl->IncW();
	}

	if (m_pImpl->ModCount(count) >= m_pImpl->RReq())
	{
		m_pImpl->ReqR(0);
		NGMSignal(m_pImpl->r_wait);
	}

	NGMUnlock(m_pImpl->mutex);

	return true;
}

template <typename T> bool Buffer<T>::Write(T * pData, int count, bool block)
{
	if ( count > m_pImpl->Capacity() || count < 0)
	{
		Log.Error("Illegal write to buffer <%s> (write count %i, capacity %i)", m_Name, count, m_pImpl->Capacity());
		Fatal();
	}
	NGMLock(m_pImpl->mutex);

	while (m_pImpl->Free() < count)
	{
		Log.Verbose("<%s> Buffer full - waiting", m_Name);

		if (block)
		{
			NGMWait(m_pImpl->mutex, m_pImpl->w_wait);
		}
		else
		{
			NGMUnlock(m_pImpl->mutex);
			return false;
		}
	}
	//Log.Message("Writing %i elements (%i free)", count, m_pImpl->Free());

	for (int i = 0; i < count; ++i)
	{
		m_pImpl->Data()[m_pImpl->WPos()] = pData[i];
		m_pImpl->IncW();
	}

	if (m_pImpl->ModCount(count) >= m_pImpl->RReq())
	{
		m_pImpl->ReqR(0);
		NGMSignal(m_pImpl->r_wait);
	}

	NGMUnlock(m_pImpl->mutex);

	return true;
}

template <typename T> void Buffer<T>::Register() const
{
	m_pImpl->Register();
	Log.Verbose("Registered writer %s (count = %i)", m_Name, m_pImpl->GetWriters());
}
template <typename T> void Buffer<T>::Release() const
{
	m_pImpl->Release();
	Log.Verbose("Released buffer %s (count = %i)", m_Name, m_pImpl->GetWriters());
}


#include "LocationScore.h"
#include "MappedRead.h"


// explicit template instantiation
template class Buffer<MappedRead*>;
template class Buffer<LocationScore*>;

template class Buffer<MappedRead*>::_Buffer<MappedRead*>;
template class Buffer<LocationScore*>::_Buffer<LocationScore*>;
