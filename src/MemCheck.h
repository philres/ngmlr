
#ifndef __MEMCHECK_H__
#define __MEMCHECK_H__

#include <cstddef>
#include <stdio.h>
#include <stdlib.h>
#include <new>

void *operator new(size_t size, const char *file, int lineno) {
	void *p = malloc(size);
	if (p == NULL)
		throw std::bad_alloc();

	fprintf(stdout, "MEMCHECK\tNEW\t%s:%d\t%p\t%lu\n", file, lineno, p, size);
	return p;
}

/* ---------------------------------------- operator delete */

void operator delete(void *p) {
	if (p == NULL)
		return;
	fprintf(stdout, "MEMCHECK\tDEL\t-\t%p\t-\n", p);
	free(p);
}

/* ---------------------------------------- operator new[] */

void *operator new[](size_t size, const char *file, int lineno) {
	void *p = malloc(size);
	if (p == NULL)
		throw std::bad_alloc();
	fprintf(stdout, "MEMCHECK\tNEW\t%s:%d\t%p\t%lu\n", file, lineno, p, size);
	return p;
}

/* ---------------------------------------- operator delete[] */

void operator delete[](void *p) {
	if (p == NULL)
		return;
	fprintf(stdout, "MEMCHECK\tDEL\t-\t%p\t-\n", p);
	free(p);

}


#define new new(__FILE__, __LINE__)
//#define new NEW_DEBUG

#endif //__MEMCHECK_H__
