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

#ifndef _WIN32

#include "NGMThreads.h"

NGMThread NGMCreateThread(pNGMThreadFunc func, int tid)
{
	return NGMCreateThread(func, new int(tid));
}
NGMThread NGMCreateThread(pNGMThreadFunc func, void* data, bool detached)
{
	pthread_t hThread;
	pthread_create(&hThread, 0, func, data);
	if (detached)
		pthread_detach(hThread);
	return hThread;
}
void NGMJoin(NGMThread * thread)
{
	pthread_join(*thread, 0);
}
void NGMJoinAll(NGMThread * thread, int count)
{
	for (int i = 0; i < count; ++i)
		NGMJoin(&thread[i]);
}

void NGMInitWait(NGMThreadWait * wait)
{
	pthread_cond_init(wait, 0);
}
void NGMWait(NGMMutex * mutex, NGMThreadWait * wait)
{
	pthread_cond_wait(wait, mutex);
}
void NGMSignal(NGMThreadWait * wait)
{
	pthread_cond_broadcast(wait);
}

void NGMInitMutex(NGMMutex * mutex)
{
	pthread_mutex_init(mutex, 0);
}
void NGMLock(NGMMutex * mutex)
{
	pthread_mutex_lock(mutex);
}
void NGMUnlock(NGMMutex * mutex)
{
	pthread_mutex_unlock(mutex);
}

void NGMOnce(NGMOnceControl * once_control, void (*init_func)())
{
	pthread_once(once_control, init_func);
}

#endif
