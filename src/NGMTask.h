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

#ifndef __NGMTASK_H__
#define __NGMTASK_H__

class NGMTask
 {
	int m_ReadProvider;

protected:
	friend class _NGM;
	int m_TID;
	bool m_FinishedStage;

public:
	void Run();
	void FinishStage();
	virtual void DoRun() = 0;
	virtual int GetStage() const = 0;
	const virtual char* GetName() const = 0;

	virtual ~NGMTask() {}
};

#endif
