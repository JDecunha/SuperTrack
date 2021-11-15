#include "ThreadAllocation.hh"

//
//Thread Allocation Definitions
//

ThreadAllocation::ThreadAllocation(const Int_t& threadID, const Int_t& randomSeed, const Int_t& nOversamples) 
{
	_threadID = threadID;
	_randomSeed = randomSeed;
	_nOversamples = nOversamples;
	_nTasks = 0;
}

void ThreadAllocation::AddTask(ThreadTask task)
{
	_tasks.push_back(task);
	_nTasks += 1;
}

std::vector<ThreadTask> ThreadAllocation::GetTasks()
{
	return _tasks;
}

Int_t ThreadAllocation::GetRandomSeed()
{
	return _randomSeed;
}

Int_t ThreadAllocation::GetNOversamples()
{
	return _nOversamples;
}

Int_t ThreadAllocation::GetThreadID()
{
	return _threadID;
}