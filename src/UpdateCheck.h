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
