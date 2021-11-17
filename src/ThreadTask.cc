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

TString ThreadTask::GetFilename(){ return _fileName; }
Long_t ThreadTask::GetEntryPoint(){ return _entryPoint; }
Long_t ThreadTask::GetExitPoint(){ return _exitPoint; }
Long_t ThreadTask::GetRandomSeed() { return _randomSeed; }
Int_t ThreadTask::GetNOversamples(){ return _nOversamples; }
Int_t ThreadTask::GetThreadID(){ return _threadID; }

