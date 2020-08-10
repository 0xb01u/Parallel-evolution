
CFLAGS= -I./original_src -O3
CFLAGS_DEBUG= -I./original_src -g -O0
LIBS=-lm
CC=mpicc

evolution: evolution_MPI_g110.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

debug: evolution_MPI_g110.c
	$(CC) $(CFLAGS_DEBUG) -DDEBUG $< -o $@ $(LIBS)

clean:
	rm -f evolution debug

