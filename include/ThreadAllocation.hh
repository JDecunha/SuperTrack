#pragma once

#include "ThreadTask.hh"
#include <vector>
#include "TROOT.h"


class ThreadAllocation
{
	public:
		ThreadAllocation(const Int_t& threadID, const Int_t& randomSeed);

		void AddTask(ThreadTask task);
		std::vector<ThreadTask> GetTasks();

	private:
		std::vector<ThreadTask> _tasks;
		
		Int_t _threadID;
		Int_t _randomSeed;
		Int_t _nTasks;
};

