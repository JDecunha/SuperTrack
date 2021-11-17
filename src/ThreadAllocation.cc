#include "ThreadAllocation.hh"

//
//Thread Allocation Definitions
//

ThreadAllocation::ThreadAllocation(const Int_t& threadID, const Long_t& randomSeed, const Int_t& nOversamples) 
{
	_threadID = threadID;
	_randomSeed = randomSeed;
	_nOversamples = nOversamples;
}

void ThreadAllocation::AddTask(ThreadTask& task)
{
	//Grab the number of tasks before pushing back another
	size_t ntasks = _tasks.size();

	//Push back
	_tasks.push_back(task);

	//Append thread specific properties to the task
	_tasks[ntasks]._randomSeed = _randomSeed;
	_tasks[ntasks]._threadID = _threadID;
	_tasks[ntasks]._nOversamples = _nOversamples;
}

std::vector<ThreadTask> ThreadAllocation::GetTasks()
{
	return _tasks;
}

Long_t ThreadAllocation::GetRandomSeed()
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