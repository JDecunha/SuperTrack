#pragma once

class ThreadAllocation;
#include "TROOT.h"

class ThreadTask
{
	public:
		ThreadTask(const TString& fileName,const Long_t& entryPoint,const Long_t& exitPoint);

		TString GetFilename() const;
		Long_t GetEntryPoint() const;
		Long_t GetExitPoint() const;
		Long_t GetRandomSeed() const;
		Int_t GetThreadID() const;
		Int_t GetNOversamples() const;

	private:
		friend ThreadAllocation;

		TString _fileName;
		Long_t _entryPoint;
		Long_t _exitPoint;

		//These are pointers to the corresponding fields of ThreadAllocation
		//These pointers are valid because ThreadTasks belong to a ThreadAllocation
		//So when ThreadAllocation goes out of scope so does the task
		Int_t _threadID;
		Long_t _randomSeed;
		Int_t _nOversamples;
};

