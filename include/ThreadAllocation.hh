#pragma once

#include "ThreadTask.hh"
#include <vector>
#include "TROOT.h"


class ThreadAllocation
{
	public:
		ThreadAllocation(const Int_t& threadID, const Int_t& randomSeed, const Int_t& nOversamples);

		void AddTask(ThreadTask task);
		std::vector<ThreadTask> GetTasks();

		Int_t GetRandomSeed();
		Int_t GetThreadID();
		Int_t GetNOversamples();

	private:
		std::vector<ThreadTask> _tasks;
		
		Int_t _threadID;
		Int_t _randomSeed;
		Int_t _nTasks;
		Int_t _nOversamples;
};

