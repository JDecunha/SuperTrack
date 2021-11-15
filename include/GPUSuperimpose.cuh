#pragma once

class TString;
class TH1F;
#include "ThreadAllocation.hh"
#include "TROOT.h"

TH1F score_lineal_GPU(TString filepath, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, Int_t nthreads, Int_t nSamples = 1, Long_t random_seed = time(NULL));
TH1F score_lineal_GPU_New(std::vector<ThreadAllocation> threadAllocations, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter);