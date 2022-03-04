//SuperTrack
#include "Histogram.cuh"
#include "VolumeEdepPair.cuh"
//gTest
#include "gtest/gtest.h"

//Kernel fills volume edep-pair with increasing volume numbers, but an edep of 1 in each
__global__ void FillVolumeEdepPairConstEIncreasingVolume(VolumeEdepPair pairToFill, double edepValue)
{
    //Determine index and stride
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < *(pairToFill.numElements); i+=stride)
    {
        pairToFill.volume[i] = i; //VolumeID increments from 0 to numElements
        pairToFill.edep[i] = edepValue; //edep is just 1 for every value
    }
}

//Kernel fills volume edep-pair with same volume numbers and an edep of 1 in each
__global__ void FillVolumeEdepPairConstEConstVolume(VolumeEdepPair pairToFill, double edepValue)
{
    //Determine index and stride
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < *(pairToFill.numElements); i+=stride)
    {
        pairToFill.volume[i] = 1; //VolumeID increments from 0 to numElements
        pairToFill.edep[i] = edepValue; //edep is just 1 for every value
    }

}

TEST(HistogramTests, HandlesVectorOfOneInDifferentVolumes)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(500,0,500,"lin");

    //Step 2.) Create a test VolumeEdepPair (testPair) and allocate memory for pairs
    int numEdeps = 100;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the testPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the testPair with predefined inputs, a constant edep value each in a different volume
    double edepValue = 1;
    FillVolumeEdepPairConstEIncreasingVolume<<<1,1>>>(testPair,edepValue);
    cudaDeviceSynchronize();

    //Step 5.) SortReduceAndAddToHistogram the testPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Step 6.) Free memory
    testHistogram.FreeTrackProcess(); //This transfers the histogram to CPU, before freeing up GPU memory used
    testPair.Free(); //Free the GPU memory held by the testing volume-edep pair

    //Step 7.) Get the CPU histogram and call assertions
    TH1D output = testHistogram.GetCPUHistogram();

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps);
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }

    }

}

TEST(HistogramTests, HandlesVectorOfFiveInDifferentVolumes)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(500,0,500,"lin");

    //Step 2.) Create a test VolumeEdepPair (testPair) and allocate memory for pairs
    int numEdeps = 100;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the testPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the testPair with predefined inputs, a constant edep value each in a different volume
    double edepValue = 5;
    FillVolumeEdepPairConstEIncreasingVolume<<<1,1>>>(testPair,edepValue);
    cudaDeviceSynchronize();

    //Step 5.) SortReduceAndAddToHistogram the testPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Step 6.) Free memory
    testHistogram.FreeTrackProcess(); //This transfers the histogram to CPU, before freeing up GPU memory used
    testPair.Free(); //Free the GPU memory held by the testing volume-edep pair

    //Step 7.) Get the CPU histogram and call assertions
    TH1D output = testHistogram.GetCPUHistogram();

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps);
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }

    }

}

TEST(HistogramTests, HandlesVectorOfOnesInSameVolume)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(500,0,1000,"lin");

    //Step 2.) Create a test VolumeEdepPair (testPair) and allocate memory for pairs
    int numEdeps = 100;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the testPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the testPair with predefined inputs, a constant edep value each in a different volume
    double edepValue = 5;
    FillVolumeEdepPairConstEConstVolume<<<1,1>>>(testPair,edepValue);
    cudaDeviceSynchronize();

    //Step 5.) SortReduceAndAddToHistogram the testPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Step 6.) Free memory
    testHistogram.FreeTrackProcess(); //This transfers the histogram to CPU, before freeing up GPU memory used
    testPair.Free(); //Free the GPU memory held by the testing volume-edep pair

    //Step 7.) Get the CPU histogram and call assertions
    TH1D output = testHistogram.GetCPUHistogram();

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue*numEdeps && high_edge > edepValue*numEdeps) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, 1) << "Counts equal to " << counts << " rather than " << numEdeps*edepValue << " in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }

    }

}

TEST(HistogramTests, HandlesVectorOfThousandInSameVolume)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(30000,1e-1,1e7,"lin");

    //Step 2.) Create a test VolumeEdepPair (testPair) and allocate memory for pairs
    int numEdeps = 1000;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the testPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the testPair with predefined inputs, a constant edep value each in a different volume
    double edepValue = 1000;
    FillVolumeEdepPairConstEConstVolume<<<1,1>>>(testPair,edepValue);
    cudaDeviceSynchronize();

    //Step 5.) SortReduceAndAddToHistogram the testPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Step 6.) Free memory
    testHistogram.FreeTrackProcess(); //This transfers the histogram to CPU, before freeing up GPU memory used
    testPair.Free(); //Free the GPU memory held by the testing volume-edep pair

    //Step 7.) Get the CPU histogram and call assertions
    TH1D output = testHistogram.GetCPUHistogram();

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue*numEdeps && high_edge > edepValue*numEdeps) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, 1) << "Counts equal to " << counts << " rather than " << numEdeps*edepValue << " in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }

    }

}

TEST(HistogramTests, HandlesVectorOfOneInDifferentVolumesLogarithmic)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(500,0,500,"log");

    //Step 2.) Create a test VolumeEdepPair (testPair) and allocate memory for pairs
    int numEdeps = 100;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the testPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the testPair with predefined inputs, a constant edep value each in a different volume
    double edepValue = 1;
    FillVolumeEdepPairConstEIncreasingVolume<<<1,1>>>(testPair,edepValue);
    cudaDeviceSynchronize();

    //Step 5.) SortReduceAndAddToHistogram the testPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Step 6.) Free memory
    testHistogram.FreeTrackProcess(); //This transfers the histogram to CPU, before freeing up GPU memory used
    testPair.Free(); //Free the GPU memory held by the testing volume-edep pair

    //Step 7.) Get the CPU histogram and call assertions
    TH1D output = testHistogram.GetCPUHistogram();

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps);
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }

    }

}

TEST(HistogramTests, HandlesVectorOfFiveInDifferentVolumesLogarithmic)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(500,0,500,"log");

    //Step 2.) Create a test VolumeEdepPair (testPair) and allocate memory for pairs
    int numEdeps = 100;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the testPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the testPair with predefined inputs, a constant edep value each in a different volume
    double edepValue = 5;
    FillVolumeEdepPairConstEIncreasingVolume<<<1,1>>>(testPair,edepValue);
    cudaDeviceSynchronize();

    //Step 5.) SortReduceAndAddToHistogram the testPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Step 6.) Free memory
    testHistogram.FreeTrackProcess(); //This transfers the histogram to CPU, before freeing up GPU memory used
    testPair.Free(); //Free the GPU memory held by the testing volume-edep pair

    //Step 7.) Get the CPU histogram and call assertions
    TH1D output = testHistogram.GetCPUHistogram();

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps);
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }

    }

}

TEST(HistogramTests, HandlesVectorOfOnesInSameVolumeLogarithmic)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(500,0,1000,"log");

    //Step 2.) Create a test VolumeEdepPair (testPair) and allocate memory for pairs
    int numEdeps = 100;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the testPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the testPair with predefined inputs, a constant edep value each in a different volume
    double edepValue = 5;
    FillVolumeEdepPairConstEConstVolume<<<1,1>>>(testPair,edepValue);
    cudaDeviceSynchronize();

    //Step 5.) SortReduceAndAddToHistogram the testPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Step 6.) Free memory
    testHistogram.FreeTrackProcess(); //This transfers the histogram to CPU, before freeing up GPU memory used
    testPair.Free(); //Free the GPU memory held by the testing volume-edep pair

    //Step 7.) Get the CPU histogram and call assertions
    TH1D output = testHistogram.GetCPUHistogram();

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue*numEdeps && high_edge > edepValue*numEdeps) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, 1) << "Counts equal to " << counts << " rather than " << numEdeps*edepValue << " in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }

    }

}

TEST(HistogramTests, HandlesVectorOfThousandInSameVolumeLogarithmic)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(30000,1e-1,1e7,"log");

    //Step 2.) Create a test VolumeEdepPair (testPair) and allocate memory for pairs
    int numEdeps = 1000;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the testPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the testPair with predefined inputs, a constant edep value each in a different volume
    double edepValue = 1000;
    FillVolumeEdepPairConstEConstVolume<<<1,1>>>(testPair,edepValue);
    cudaDeviceSynchronize();

    //Step 5.) SortReduceAndAddToHistogram the testPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Step 6.) Free memory
    testHistogram.FreeTrackProcess(); //This transfers the histogram to CPU, before freeing up GPU memory used
    testPair.Free(); //Free the GPU memory held by the testing volume-edep pair

    //Step 7.) Get the CPU histogram and call assertions
    TH1D output = testHistogram.GetCPUHistogram();

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue*numEdeps && high_edge > edepValue*numEdeps) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, 1) << "Counts equal to " << counts << " rather than " << numEdeps*edepValue << " in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }

    }
}



