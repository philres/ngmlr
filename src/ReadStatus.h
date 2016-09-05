/**
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Contact: philipp.rescheneder@univie.ac.at
 */

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
