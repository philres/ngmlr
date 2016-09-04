#ifndef __TYPES_H__
#define __TYPES_H__

typedef unsigned int uint;

//Type for holding genomic locations
typedef long long loc;
typedef unsigned long long uloc;

#ifndef _WIN32
typedef unsigned long ulong;
#endif
#ifdef _WIN32
typedef unsigned long long ulong;
#endif

//#define INSTANCE_COUNTING

#endif
