#ifndef __TYPES_H__
#define __TYPES_H__

typedef unsigned int uint;

//Type for holding genomic locations
//typedef unsigned long long uloc;
//typedef unsigned long long uloc;

//Temporary, used for identifying illegal conversions
typedef long long loc;



#define ULOC_WRAPPER

#ifdef ULOC_WRAPPER

	class uloc
	{
		unsigned long long _val;

	public:
		explicit uloc() { _val = 0; }
		//explicit uloc(unsigned int v) { _val = v; }
		uloc(const uloc& other) { _val = other._val; }
		const uloc& operator=(const uloc& other) { _val = other._val; return *this; }

		const uloc operator+(const uloc& other) const { uloc sum(*this); sum._val = _val + other._val; return sum; }
		const uloc operator+(const uint& other) const { uloc sum(*this); sum._val = _val + other; return sum; }
		const uloc operator+(const int& other) const { uloc sum(*this); sum._val = _val + other; return sum; }
		const uloc operator-(const uloc& other) const { uloc sum(*this); sum._val = _val - other._val; return sum; }
		const uloc operator-(const uint& other) const { uloc sum(*this); sum._val = _val - other; return sum; }
		const uloc operator-(const int& other) const { uloc sum(*this); sum._val = _val - other; return sum; }

		const uloc operator*(const uloc& other) const { uloc sum(*this); sum._val = _val * other._val; return sum; }
		const uloc operator*(const uint& other) const { uloc sum(*this); sum._val = _val * other; return sum; }
		const uloc operator*(const int& other) const { uloc sum(*this); sum._val = _val * other; return sum; }
		const uloc operator/(const uloc& other) const { uloc sum(*this); sum._val = _val / other._val; return sum; }
		const uloc operator/(const uint& other) const { uloc sum(*this); sum._val = _val / other; return sum; }
		const uloc operator/(const int& other) const { uloc sum(*this); sum._val = _val / other; return sum; }

		const uloc& operator++() { _val++; return *this; }
		const uloc& operator--() { _val--; return *this; }

		const uloc operator++(int) { uloc tmp = *this; _val++; return tmp; }
		const uloc operator--(int) { uloc tmp = *this; _val--; return tmp; }

		const uloc& operator+=(const uloc& other) { _val = _val + other._val; return *this; }
		const uloc& operator+=(const uint& other) { _val = _val + other; return *this; }
		const uloc& operator+=(const int& other) { _val = _val + other; return *this; }
		const uloc& operator-=(const uloc& other) { _val = _val - other._val; return *this; }
		const uloc& operator-=(const uint& other) { _val = _val - other; return *this; }
		const uloc& operator-=(const int& other) { _val = _val - other; return *this; }
		const uloc operator>>(const uint& other) const { uloc sum(*this); sum._val = _val >> other; return sum; }
		const uloc operator<<(const uint& other) const { uloc sum(*this); sum._val = _val << other; return sum; }
		const uloc operator|(const uint& other) const { uloc sum(*this); sum._val = _val | other; return sum; }
		const uloc operator&(const uint& other) const { uloc sum(*this); sum._val = _val & other; return sum; }

		bool operator==(const uloc& other) const { return _val == other._val; }
		bool operator!=(const uloc& other) const { return _val != other._val; }
		bool operator<(const uloc& other) const { return _val > other._val; }
		bool operator>(const uloc& other) const { return _val < other._val; }
		bool operator<=(const uloc& other) const { return _val >= other._val; }
		bool operator>=(const uloc& other) const { return _val <= other._val; }

		bool operator<(const uint& other) const { return _val > other; }
		bool operator>(const uint& other) const { return _val < other; }
		bool operator<=(const uint& other) const { return _val >= other; }
		bool operator>=(const uint& other) const { return _val <= other; }

		bool operator<(const int& other) const { return _val > other; }
		bool operator>(const int& other) const { return _val < other; }
		bool operator<=(const int& other) const { return _val >= other; }
		bool operator>=(const int& other) const { return _val <= other; }

		friend const uloc operator+(int left,const uloc& right) { uloc sum(right); sum._val += left; return sum; }
		friend const uloc operator+(uint left,const uloc& right) { uloc sum(right); sum._val += left; return sum; }
		friend const uloc operator-(int left,const uloc& right) { uloc sum; sum._val = left; sum._val -= right._val; return sum; }
		friend const uloc operator-(uint left,const uloc& right) { uloc sum; sum._val = left; sum._val -= right._val; return sum; }
		friend const uloc operator/(int left,const uloc& right) { uloc sum(right); sum._val /= left; return sum; }
		friend const uloc operator/(uint left,const uloc& right) { uloc sum(right); sum._val /= left; return sum; }
		friend const uloc operator*(int left,const uloc& right) { uloc sum; sum._val = left; sum._val *= right._val; return sum; }
		friend const uloc operator*(uint left,const uloc& right) { uloc sum; sum._val = left; sum._val *= right._val; return sum; }

		friend const bool operator>(int left,const uloc& right) { return left > right._val; }
		friend const bool operator<(int left,const uloc& right) { return left < right._val; }
		friend const bool operator>(uint left,const uloc& right) { return left > right._val; }
		friend const bool operator<(uint left,const uloc& right) { return left < right._val; }
		friend const bool operator>=(int left,const uloc& right) { return left >= right._val; }
		friend const bool operator<=(int left,const uloc& right) { return left <= right._val; }
		friend const bool operator>=(uint left,const uloc& right) { return left >= right._val; }
		friend const bool operator<=(uint left,const uloc& right) { return left <= right._val; }

		static uloc from_uint32(uint v)  { uloc ul; ul._val = v; return ul; }
		static uint to_uint32(const uloc& other) { return other._val; }
		static uloc from_int32(int v)  { uloc ul; ul._val = v; return ul; }
		static int to_int32(const uloc& other) { return other._val; }
		static uloc from_loc(loc v)  { uloc ul; ul._val = v; return ul; }
		static loc to_loc(const uloc& other) { return (loc) other._val; }
		static uloc from_uloc(unsigned long long v)  { uloc ul; ul._val = v; return ul; }
		static unsigned long long to_uloc(const uloc& other) { return (unsigned long long) other._val; }

	};

	#define ULOC_FROM_UINT32 uloc::from_uint32
	#define ULOC_TO_UINT32   uloc::to_uint32

	#define ULOC_FROM_INT32  uloc::from_int32
	#define ULOC_TO_INT32    uloc::to_int32

	#define ULOC_FROM_LOC    uloc::from_loc
	#define ULOC_TO_LOC      uloc::to_loc

	#define ULOC_FROM_ULOC   uloc::from_uloc
	#define ULOC_TO_ULOC     uloc::to_uloc

#else

	typedef unsigned long long uloc;

	#define ULOC_FROM_UINT32
	#define ULOC_TO_UINT32

	#define ULOC_FROM_INT32
	#define ULOC_TO_INT32

	#define ULOC_FROM_LOC
	#define ULOC_TO_LOC

	#define ULOC_FROM_ULOC
	#define ULOC_TO_ULOC

#endif


#ifndef _WIN32
typedef unsigned long ulong;
#endif
#ifdef _WIN32
typedef unsigned long long ulong;
#endif

//#define INSTANCE_COUNTING

#endif
