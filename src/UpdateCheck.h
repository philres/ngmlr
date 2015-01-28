/*
 * UpdateCheck.h
 *
 *  Created on: October 2, 2014
 *      Author: moritz
 */

class UpdateCheckInterface
{
public:
	//Remind the user to perform an update-check after the current binary has a certain age
	static void reminder();

	//Perform a version check with the CIBIV server
	static void remoteCheck();
};
