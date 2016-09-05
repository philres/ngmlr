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

#ifndef __UPDATE_CHECK_H__
#define __UPDATE_CHECK_H__

class UpdateCheckInterface
{
public:
	//Remind the user to perform an update-check after the current binary has a certain age
	static void reminder();

	//Perform a version check with the CIBIV server
	static void remoteCheck();
};

#endif//__UPDATE_CHECK_H__
