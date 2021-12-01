#pragma once
//CUB
#include <cub/cub.cuh>

struct CUBAddOperator
{
	//Because this is being inlined I think we have to have it in the header file
    template <typename T>
    CUB_RUNTIME_FUNCTION __forceinline__ T operator()(const T &a, const T &b) const
    {
		return a+b;
	}
};