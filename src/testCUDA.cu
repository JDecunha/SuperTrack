#include <iostream>
#include <stdio.h>
#include "../include/testCUDA.cuh"

__global__ void print_GPU()
{
 printf("This is running on a GPU! \n");
}

void print_test()
{
 std::cout << "Program continues to print_test" << std::endl;
 print_GPU<<<40,32>>>();
 std::cout << "Program leaves kernel" << std::endl;
 cudaError_t cudaerr = cudaDeviceSynchronize();
 if (cudaerr != cudaSuccess)
 {
        printf("kernel launch failed with error \"%s\".\n",cudaGetErrorString(cudaerr));
 }

}