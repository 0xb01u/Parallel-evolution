
CFLAGS= -I./original_src -O3
CFLAGS_DEBUG= -I./original_src -g -O0
LIBS=-lm -fopenmp

evolution: evolution_OpenMP_g110.c
	gcc $(CFLAGS) $< -o $@ $(LIBS)

debug: evolution_OpenMP_g110.c
	gcc $(CFLAGS_DEBUG) -DDEBUG $< -o $@ $(LIBS)

clean:
	rm -f evolution debug

