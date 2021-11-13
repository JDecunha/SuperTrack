#include "ThreadAllocation.hh"

//
//Thread Allocation Definitions
//

ThreadAllocation::ThreadAllocation(const Int_t& threadID, const Int_t& randomSeed) 
{
	_threadID = threadID;
	_randomSeed = randomSeed;
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