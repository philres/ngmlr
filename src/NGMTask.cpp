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

#include "NGM.h"
#include "Log.h"

void NGMTask::FinishStage()
{
	if (!m_FinishedStage)
	{
		NGM.FinishStage(m_TID);
		m_FinishedStage = true;
	}
}

void NGMTask::Run() {
	m_FinishedStage = false;

	try {
		DoRun();
	} catch (std::bad_alloc & ex) {
		Log.Error("Exception bad_alloc occured in thread %i. This usually means you ran out of physical or virtual memory (try ulimit -v)", m_TID);
	}
	catch (...)
	{
		Log.Error("Exception in thread %i", m_TID);
		throw;
	}
}
