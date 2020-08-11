# About
This implementation of the evolution simulation is the same as [the CUDA implementation](../CUDA), but using NVIDIA CUDA Driver API for native C, instead of CUDA's own syntax.

The kernels are contained in a unique `.cu` file and the header files, separately from the main `.c` program. This `evolution_kernels.cu` file is compiled to `kernels.ptx` using `nvcc`. `evolution_CUDA_g110.c` is compiled using `gcc`.
