#pragma once

#include "Histogram.cuh"
#include "SimulationMethod.hh"
#include "ThreadAllocation.hh"
#include <vector>

class SuperTrackManager
{	
	public:
		SuperTrackManager();
		void AddThreadAllocations(std::vector<ThreadAllocation>& allocations);
		void AddHistogram(Histogram* histogram);
		void AddSimulationMethod(SimulationMethod* method);
		void Initialize();
		void Run();

		void TestSuperimpose(std::vector<ThreadAllocation>& allocations, Histogram* histogram, SimulationMethod* method);

	private:
		std::vector<ThreadAllocation> _threadAllocations;
		Histogram* _histogram;
		SimulationMethod* _method;
		bool _bInitialized;
};