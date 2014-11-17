#ifndef __TYPES_H__
#define __TYPES_H__

typedef unsigned int uint;

//Type for holding genomic locations
typedef long long loc;
typedef unsigned long long uloc;

#define ULOC_FROM_UINT32
#define ULOC_TO_UINT32

#define ULOC_FROM_INT32
#define ULOC_TO_INT32
#define ULOC_FROM_LOC
#define ULOC_TO_LOC

#define MAKE_ULOC
#define GET_ULOC

#ifndef _WIN32
typedef unsigned long ulong;
#endif
#ifdef _WIN32
typedef unsigned long long ulong;
#endif

//#define INSTANCE_COUNTING

#endif
