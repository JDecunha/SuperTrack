#pragma once

#include "TROOT.h"

class ThreadTask
{
	public:
		ThreadTask(const TString& fileName,const Long_t& entryPoint,const Long_t& exitPoint);

		TString GetFilename();
		Long_t GetEntryPoint();
		Long_t GetExitPoint();

	private:
		TString _fileName;
		Long_t _entryPoint;
		Long_t _exitPoint;
};

