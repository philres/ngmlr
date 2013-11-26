#ifndef __READSTATUS_H__
#define __READSTATUS_H__

namespace NGMNames
{
	enum ReadStatus
	{
		Unknown			= 0x000,

		NoSrcPair		= 0x010,		
		PairedFail		= 0x020,	
		Empty			= 0x040,
		DeletionPending	= 0x100
	};
}

#endif
