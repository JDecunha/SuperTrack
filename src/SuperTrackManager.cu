#include "SuperTrackManager.hh"

SuperTrackManager::SuperTrackManager() { }

void SuperTrackManager::AddThreadAllocations(std::vector<ThreadAllocation>& allocations)
{
	_threadAllocations = std::move(allocations);
}

void SuperTrackManager::AddHistogram(Histogram* histogram)
{
	_histogram = histogram;
}

void SuperTrackManager::AddSimulationMethod(SimulationMethod* method)
{
	_method = method;
}

void SuperTrackManager::TestSuperimpose(std::vector<ThreadAllocation>& allocations, Histogram* histogram, SimulationMethod* method)
{

}