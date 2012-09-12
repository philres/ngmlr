#include "NGMStats.h"

#ifndef _WIN32
#include <sys/ipc.h>
#include <sys/shm.h>
#endif

#include <errno.h>
#include <memory.h>

#include "Log.h"

#undef module_name
#define module_name "STATS"

NGMStats * pStats = 0;

NGMStats * NGMStats::InitStats(char const * const AppName)
{
#ifndef _WIN32
	key_t key = ftok(AppName, 's');

	if (key == (key_t)-1)
	{
		Log.Warning("Unable to obtain stats smem key (error %i)", errno);
	}
	else
	{
		int shm_id = shmget( key, sizeof(NGMStats), IPC_CREAT | 0644 );

		if (shm_id == -1)
		{
			Log.Warning("Unable to get stats smem (error %i)", errno);
		}
		else
		{
			void * pmem = shmat(shm_id, 0, 0);
			if (pmem != (void*)-1)
			{
				pStats = (NGMStats*)pmem;
				memset(pmem, 0, sizeof(NGMStats));
				pStats->Cookie = 0x51415001;
				Log.Message("Exporting stats at 0x%x (Key 0x%x, ID 0x%x)", pmem, key, shm_id);
				return pStats;
			}
			else
			{
				Log.Warning("Unable to access stats smem (error %i)", errno);
			}
		}
	}
#endif
	return pStats = new NGMStats();
}

NGMStats::NGMStats()
{
	csTime = 0.0f;
	csLength = 0;
	csOverflows = 0;

	validPairs = 0.0f;
	insertSize = 0.0f;
}

NGMStats::~NGMStats()
{

}
