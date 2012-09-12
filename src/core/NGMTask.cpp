#include "NGM.h"
#include "Debug.h"

void NGMTask::FinishStage()
{
	if (!m_FinishedStage)
	{
		NGM.FinishStage(m_TID);
		m_FinishedStage = true;
	}
}

void NGMTask::Run()
{
	m_FinishedStage = false;

	try
	{
		DoRun();
	}
	catch (std::bad_alloc)
	{
		Log.Error("Exception bad_alloc occured in thread %i. This usually means you ran out of physical or virtual memory (try ulimit -v)", m_TID);
	}
	//catch (...)
	//{
	//	Log.Error("Exception in thread %i", m_TID);
	//	throw;
	//}
}
