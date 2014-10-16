#ifndef __TYPES_H__
#define __TYPES_H__

typedef unsigned int uint;

typedef unsigned long long uint64;

#ifndef _WIN32
typedef unsigned long ulong;
#endif
#ifdef _WIN32
typedef unsigned long long ulong;
#endif

//#define INSTANCE_COUNTING

#endif
