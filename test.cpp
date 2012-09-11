#include "Buffer.h"
#include "CS.h"

#include "PrefixTable.h"
#include "Types.h"

bool CompareRefEntry(RefEntry const * lhs, RefEntry const * rhs);

void CompareTables(PrefixTable * table1, PrefixTable * table2)
{
	for (uint i = 0; i < 1000000; ++i)
	{
		if (!CompareRefEntry(table1->GetRefEntry(i), table2->GetRefEntry(i)))
		{
			Log.Error("Mismatch on prefix 0x%x", i);
		}
	}
	Log.Green("Compare done");
}

bool CompareRefEntry(RefEntry const * lhs, RefEntry const * rhs)
{
	return lhs->refTotal == rhs->refTotal;
}

#include <malloc.h>

void TestMem()
{
#ifdef INSTANCE_COUNTING
	Log.Green("Counts:");
	Log.Message("MappedRead count = %i", MappedRead::sInstanceCount);
	Log.Message("LocationScore count = %i", LocationScore::sInstanceCount);
#endif
	malloc_stats();
}
