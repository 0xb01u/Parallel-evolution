
CFLAGS= -I./original_src -O3 -DCUDA_CHECK --compiler-options -Wall
CFLAGS_DEBUG= -I./original_src -O0 -G -DCUDA_CHECK --compiler-options -Wall
LIBS=
CC=nvcc

evolution: evolution_CUDA_g110.cu *.h
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

debug: evolution_CUDA_g110.cu *.h
	$(CC) $(CFLAGS_DEBUG) -DDEBUG $< -o $@ $(LIBS)

tgz:
	tar -czf evolution_cuda.tgz Makefile evolution.cu *.h

clean:
	rm -f evolution debug

