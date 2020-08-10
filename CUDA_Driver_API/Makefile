
CFLAGS= -Wall -I./include -O3 -DCUDA_CHECK
CFLAGS_DEBUG= -Wall -I./include -O0 -g -DCUDA_CHECK
LIBS= -lcuda
CC=gcc
PTX=nvcc -ptx --compiler-options

evolution: evolution_CUDA_g110.c include/*.h kernels.ptx
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

debug: evolution_CUDA_g110.c include/*.h kernels_debug.ptx
	$(CC) $(CFLAGS_DEBUG) -DDEBUG $< -o $@ $(LIBS)

kernels.ptx: evolution_kernels.cu include/*h
	$(PTX) $(CFLAGS) $< -o $@

ptx: kernels.ptx

kernels_debug.ptx: evolution_kernels.cu include/*.h
	$(PTX) $(CFLAGS_DEBUG) -DDEBUG $< -o $@

tgz:
	tar -czf evolution_cuda_api.tgz Makefile evolution_CUDA_g110.c evolution_kernels.cu include

clean:
	rm -f evolution debug *.ptx *.tgz

del:
	rm -f evolution debug
