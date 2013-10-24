#include <seqan/basic.h>
using namespace seqan;

struct MyClass
{
};

int main()
{
	MyClass* my_class_arr;
	allocate(Default(), my_class_arr, 100);
	arrayConstruct(my_class_arr, my_class_arr + 100);
	arrayDestruct(my_class_arr, my_class_arr + 100);
	deallocate(Default(), my_class_arr, 100);
	Allocator<SimpleAlloc< > > alloc1;
	char * char_array;
	allocate(alloc1, char_array, 300);
	clear(alloc1);
	return 0;
}
