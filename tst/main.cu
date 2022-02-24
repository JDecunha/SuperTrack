//SuperTrack
#include "Histogram.cuh"
#include "VolumeEdepPair.cuh"
//gTest
#include "gtest/gtest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

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

TEST(HistogramTests, HandlesVectorOfOnesInDifferentVolumes)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(256,1e-1,1e1,"lin");

    //Step 2.) Create VolumeEdepPair and allocate for 100 points
    int numEdeps = 100;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the VolumeEdepPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the VolumeEdepPair with a the required inputs. Write a kernel for this? 
    FillVolumeEdepPairConstEIncreasingVolume<<<1,1>>>(testPair,1);
    cudaDeviceSynchronize();

    //Step 5.) Call SortReduceAndAddToHistogram on the VolumeEdepPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Optional: Print GPU Histogram
    testHistogram.Print();

    //Step 6.) Transfer the VolumeEdepPair to the CPU and call assertions on it
}

TEST(HistogramTests, HandlesVectorOfFivesInDifferentVolumes)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(256,1e-1,1e1,"lin");

    //Step 2.) Create VolumeEdepPair and allocate for 100 points
    int numEdeps = 100;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the VolumeEdepPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the VolumeEdepPair with a the required inputs. Write a kernel for this? 
    FillVolumeEdepPairConstEIncreasingVolume<<<1,1>>>(testPair,5);
    cudaDeviceSynchronize();

    //Step 5.) Call SortReduceAndAddToHistogram on the VolumeEdepPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Optional: Print GPU Histogram
    testHistogram.Print();

    //Step 6.) Transfer the VolumeEdepPair to the CPU and call assertions on it
}

TEST(HistogramTests, HandlesVectorOfOnesInSameVolume)
{
    //Step 1.) Create histogram
    Histogram testHistogram = Histogram(256,1e-1,1e2,"lin");

    //Step 2.) Create VolumeEdepPair and allocate for 100 points
    int numEdeps = 100;
    VolumeEdepPair testPair; testPair.Allocate(numEdeps);

    //Step 3.) Allocate the Track Process of the Histogram given the VolumeEdepPair
    testHistogram.AllocateTrackProcess(testPair);

    //Step 4.) Fill the VolumeEdepPair with a the required inputs. Write a kernel for this? 
    FillVolumeEdepPairConstEConstVolume<<<1,1>>>(testPair,1);
    cudaDeviceSynchronize();

    //Step 5.) Call SortReduceAndAddToHistogram on the VolumeEdepPair
    testHistogram.SortReduceAndAddToHistogram(testPair);
    cudaDeviceSynchronize();

    //Optional: Print GPU Histogram
    testHistogram.Print();

    //Step 6.) Transfer the VolumeEdepPair to the CPU and call assertions on it
}

