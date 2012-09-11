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
