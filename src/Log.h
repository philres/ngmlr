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

#ifndef __LOGGING_H__
#define __LOGGING_H__

#include "ILog.h"

namespace __Log
{
	void Init();
}

class _Log : public ILog
{
public:
	static _Log const & Instance();
	static void Init(char const * logFile, int logLvl);
	static void FilterLevel(int const lvl);
	static void setColor(bool const color);
	static void Cleanup();

	void _Message(int const lvl, char const * const title, char const * const msg, ...) const;
	void _Debug  (int const lvl, char const * const title, char const * const msg, ...) const;
	void Fatal();

private:

	_Log();
	~_Log();

	friend void __Log::Init();
};

#undef Log
#define Log _Log::Instance()
#define Logger _Log::Instance()

#endif
