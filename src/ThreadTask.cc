#include "ThreadTask.hh"

//
//Thread Tasks Definitions
//

ThreadTask::ThreadTask(const TString& fileName, const Long_t& entryPoint, const Long_t& exitPoint)
{
	_fileName = fileName;
	_entryPoint = entryPoint;
	_exitPoint = exitPoint;
}

TString ThreadTask::GetFilename() const { return _fileName; }
Long_t ThreadTask::GetEntryPoint() const { return _entryPoint; }
Long_t ThreadTask::GetExitPoint() const { return _exitPoint; }
Long_t ThreadTask::GetRandomSeed() const { return _randomSeed; }
Int_t ThreadTask::GetNOversamples() const { return _nOversamples; }
Int_t ThreadTask::GetThreadID() const { return _threadID; }

