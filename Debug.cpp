#define _CRTDBG_MAPALLOC

#include <stdio.h>

/*void * operator new (size_t size, char const * file, int line)
{
	void * tmp = ::new char[size];
	//if (size != 24)
		printf(">>> created %i byte block at 0x%llx (%s, line %i)\n", size, (unsigned long long)tmp, file, line);
	return tmp;
}

void operator delete(void * p, char const * file, int line)
{
	delete p;
}*/
