/*
 * Simplified simulation of life evolution
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2019/2020
 *
 * v1.6
 *
 * CHANGES:
 * 1) Float values have been substituted by fixed point arithmetics 
 *	using integers. To simplify, the fixed point arithmetics are done 
 *	with PRECISION in base 10. See precision constant in int_float.h
 * 2) It uses a portable approximation to trigonometric functions using
 *	Taylor polynomials. 
 * 3) nrand48 function has been extracted from glibc source code and 
 *	its internal API simplified to allow its use in the GPU.
 * 4) The kernel functions have been separated from the main program into a 
 *  different file.
 * 5) It now uses CUDA Driver API, instead of regular CUDA syntax, and,
 *  thus, this file can be compiled using gcc.
 *
 * (c) 2020, Arturo Gonzalez Escribano
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#include<cputils.h>
#include<int_float.h>
#include<cuda.h>

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 *	USE THIS SPACE FOR YOUR KERNEL OR DEVICE FUNTIONS
 *
 */

#include "include/cuda_check.h"
#include "include/evolution.h"

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

#ifdef DEBUG
/* 
 * Function: Print the current state of the simulation 
 */
void print_status(int iteration, int rows, int columns, int *culture, int num_cells, Cell *cells, int num_cells_alive, Statistics sim_stat) {
	/* 
	 * You don't need to optimize this function, it is only for pretty printing and debugging purposes.
	 * It is not compiled in the production versions of the program.
	 * Thus, it is never used when measuring times in the leaderboard
	 */
	int i,j;

	printf("Iteration: %d\n", iteration);
	printf("+");
	for(j=0; j<columns; j++) printf("---");
	printf("+\n");
	for(i=0; i<rows; i++) {
		printf("|");
		for(j=0; j<columns; j++) {
			char symbol;
			if (accessMat(culture, i, j) >= 20 * PRECISION) symbol = '+';
			else if (accessMat(culture, i, j) >= 10 * PRECISION) symbol = '*';
			else if (accessMat(culture, i, j) >= 5 * PRECISION) symbol = '.';
			else symbol = ' ';

			int t;
			int counter = 0;
			for(t=0; t<num_cells; t++) {
				int row = (int)(cells[t].pos_row / PRECISION);
				int col = (int)(cells[t].pos_col / PRECISION);
				if (cells[t].alive && row == i && col == j) {
					counter ++;
				}
			}
			if (counter > 9) printf("(M)");
			else if (counter > 0) printf("(%1d)", counter);
			else printf(" %c ", symbol);
		}
		printf("|\n");
	}
	printf("+");
	for(j=0; j<columns; j++) printf("---");
	printf("+\n");
	printf("Num_cells_alive: %04d\nHistory(Cells: %04d, Dead: %04d, Max.alive: %04d, Max.new: %04d, Max.dead: %04d, Max.age: %04d, Max.food: %6f)\n\n", 
		num_cells_alive, 
		sim_stat.history_total_cells, 
		sim_stat.history_dead_cells, 
		sim_stat.history_max_alive_cells, 
		sim_stat.history_max_new_cells, 
		sim_stat.history_max_dead_cells, 
		sim_stat.history_max_age,
		(float)sim_stat.history_max_food / PRECISION
	);
}
#endif

/*
 * Function: Print usage line in stderr
 */
void show_usage(char *program_name) {
	fprintf(stderr,"Usage: %s ", program_name);
	fprintf(stderr,"<rows> <columns> <maxIter> <max_food> <food_density> <food_level> <short_rnd1> <short_rnd2> <short_rnd3> <num_cells>\n");
	fprintf(stderr,"\tOptional arguments for special food spot: [ <row> <col> <size_rows> <size_cols> <density> <level> ]\n");
	fprintf(stderr,"\n");
}


/*
 * MAIN PROGRAM
 */
int main(int argc, char *argv[]) {
	int i,j;

	// Simulation data
	int max_iter;			// Maximum number of simulation steps
	int rows, columns;		// Cultivation area sizes
	int *culture;			// Cultivation area values
	int *culture_cells;		// Ancillary structure to count the number of cells in a culture space

	float max_food;			// Maximum level of food on any position
	float food_density;		// Number of food sources introduced per step
	float food_level;		// Maximum number of food level in a new source

	bool food_spot_active = false;	// Special food spot: Active
	int food_spot_row = 0;		// Special food spot: Initial row
	int food_spot_col = 0;		// Special food spot: Initial row
	int food_spot_size_rows = 0;	// Special food spot: Rows size
	int food_spot_size_cols = 0;	// Special food spot: Cols size
	float food_spot_density = 0.0f;	// Special food spot: Food density
	float food_spot_level = 0.0f;	// Special food spot: Food level

	unsigned short init_random_seq[3];	// Status of the init random sequence
	unsigned short food_random_seq[3];	// Status of the food random sequence
	unsigned short food_spot_random_seq[3];	// Status of the special food spot random sequence

	int	num_cells;		// Number of cells currently stored in the list
	Cell	*cells;			// List to store cells information

	// Statistics
	Statistics sim_stat;	
	sim_stat.history_total_cells = 0;
	sim_stat.history_dead_cells = 0;
	sim_stat.history_max_alive_cells = 0;
	sim_stat.history_max_new_cells = 0;
	sim_stat.history_max_dead_cells = 0;
	sim_stat.history_max_age = 0;
	sim_stat.history_max_food = 0.0f;

	/* 1. Read simulation arguments */
	/* 1.1. Check minimum number of arguments */
	if (argc < 11) {
		fprintf(stderr, "-- Error: Not enough arguments when reading configuration from the command line\n\n");
		show_usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	/* 1.2. Read culture sizes, maximum number of iterations */
	rows = atoi(argv[1]);
	columns = atoi(argv[2]);
	max_iter = atoi(argv[3]);

	/* 1.3. Food data */
	max_food = atof(argv[4]);
	food_density = atof(argv[5]);
	food_level = atof(argv[6]);

	/* 1.4. Read random sequences initializer */
	for(i=0; i<3; i++) {
		init_random_seq[i] = (unsigned short)atoi(argv[7+i]);
	}

	/* 1.5. Read number of cells */
	num_cells = atoi(argv[10]);

	/* 1.6. Read special food spot */
	if (argc > 11) {
		if (argc < 17) {
			fprintf(stderr, "-- Error in number of special-food-spot arguments in the command line\n\n");
			show_usage(argv[0]);
			exit(EXIT_FAILURE);
		}
		else {
			food_spot_active = true;
			food_spot_row = atoi(argv[11]);
			food_spot_col = atoi(argv[12]);
			food_spot_size_rows = atoi(argv[13]);
			food_spot_size_cols = atoi(argv[14]);
			food_spot_density = atof(argv[15]);
			food_spot_level = atof(argv[16]);

			// Check non-used trailing arguments
			if (argc > 17) {
				fprintf(stderr, "-- Error: too many arguments in the command line\n\n");
				show_usage(argv[0]);
				exit(EXIT_FAILURE);
			}
		}
	}

#ifdef DEBUG
	/* 1.7. Print arguments */
	printf("Arguments, Rows: %d, Columns: %d, max_iter: %d\n", rows, columns, max_iter);
	printf("Arguments, Max.food: %f, Food density: %f, Food level: %f\n", max_food, food_density, food_level);
	printf("Arguments, Init Random Sequence: %hu,%hu,%hu\n", init_random_seq[0], init_random_seq[1], init_random_seq[2]);
	if (food_spot_active) {
		printf("Arguments, Food_spot, pos(%d,%d), size(%d,%d), Density: %f, Level: %f\n",
			food_spot_row, food_spot_col, food_spot_size_rows, food_spot_size_cols, food_spot_density, food_spot_level);
	}
	printf("Initial cells: %d\n", num_cells);
#endif // DEBUG


	/* 1.8. Initialize random sequences for food dropping */
	for(i=0; i<3; i++) {
		food_random_seq[i] = (unsigned short)nrand48(init_random_seq);
		food_spot_random_seq[i] = (unsigned short)nrand48(init_random_seq);
	}

	/* 1.9. Initialize random sequences of cells */
	cells = (Cell *)malloc(sizeof(Cell) * (size_t)num_cells);
	if (cells == NULL) {
		fprintf(stderr,"-- Error allocating: %d cells\n", num_cells);
		exit(EXIT_FAILURE);
	}
	for(i=0; i<num_cells; i++) {
		// Initialize the cell ramdom sequences
		for(j=0; j<3; j++) 
			cells[i].random_seq[j] = (unsigned short)nrand48(init_random_seq);
	}


#ifdef DEBUG
	/* 1.10. Print random seed of the initial cells */
	/*
	printf("Initial cells random seeds: %d\n", num_cells);
	for(i=0; i<num_cells; i++)
		printf("\tCell %d, Random seq: %hu,%hu,%hu\n", i, cells[i].random_seq[0], cells[i].random_seq[1], cells[i].random_seq[2]);
	*/
#endif // DEBUG


	// CUDA start
	CUdevice device;
	CUcontext context;

	cudaCheckCall(cuInit(0));
	cudaCheckCall(cuDeviceGet(&device, 0));
	cudaCheckCall(cuCtxCreate(&context, 0, device));

	/* 2. Start global timer */
	double ttotal = cp_Wtime();

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */

/*
 * Block and thread sizes for kernel executions.
 *
 */
#define THREADS1D 128
#define THREADS THREADS1D, 1, 1
#define BLOCK (max(rows*columns, num_cells_alive)/THREADS1D + 1), 1, 1
#define BLOCK_F (max3(rows*columns, num_cells_alive, max_new_sources)/THREADS1D + 1), 1, 1	/* "Food kernels" */
#define BLOCK_C (num_cells_alive)/ THREADS1D + 1, 1, 1	/* "Cell kernels" */
#define BLOCK_P (rows*columns)/THREADS1D + 1, 1, 1		/* "Culture kernels" */


	/* 3. Initialize culture surface and initial cells */
	culture = NULL;
	culture_cells = NULL;

	/* Device equivalents */
	CUdeviceptr culture_d, culture_cells_d;
	cudaCheckCall(cuMemAlloc(&culture_d, sizeof(int) * (size_t)rows * (size_t)columns));
	cudaCheckCall(cuMemAlloc(&culture_cells_d, sizeof(short) * (size_t)rows * (size_t)columns));

	/* Set both surfaces to 0 */
	cudaCheckCall(cuMemsetD32(culture_d, 0, (size_t)rows * (size_t)columns));
	cudaCheckCall(cuMemsetD16(culture_cells_d, 0, (size_t)rows * (size_t)columns));

	/* CUDA streams */
	CUstream alt;
	cudaCheckCall(cuStreamCreate(&alt, CU_STREAM_NON_BLOCKING));

	/* Copy random cell seeds to device */
	unsigned short *random_seqs;
	cudaCheckCall(cuMemAllocHost((void **)&random_seqs, sizeof(unsigned short) * 3 * num_cells));
	CUdeviceptr random_seqs_d;
	cudaCheckCall(cuMemAlloc(&random_seqs_d, sizeof(unsigned short) * 3 * num_cells));

	for (i = 0; i < num_cells; i++)
	{
		random_seqs[3*i] = cells[i].random_seq[0];
		random_seqs[3*i + 1] = cells[i].random_seq[1];
		random_seqs[3*i + 2] = cells[i].random_seq[2];
	}

	cudaCheckCall(cuMemcpyHtoDAsync(random_seqs_d, random_seqs, sizeof(unsigned short) * 3 * num_cells, alt));

	int num_cells_alive = num_cells;

	/* Device cell lists */
	/* They are assigned 2GiB of memory each, so they are never required to realloc. */
	CUdeviceptr cells_d1, cells_d2;
	cudaCheckCall(cuMemAlloc(&cells_d1, (size_t) (1l << 30)));
	cudaCheckCall(cuMemAlloc(&cells_d2, (size_t) (1l << 30)));
	/* Device statistics */
	CUdeviceptr stats_d;
	cudaCheckCall(cuMemAlloc(&stats_d, sizeof(Statistics)));

	/* Device intialization */
	/* Module load */
	static char path[256];
	int c = (int)strlen(argv[0]) - 1;
	while (argv[0][--c] != '/');
	strncpy(path, argv[0], c + 1);

	CUmodule evolution_kernels;
	#ifdef DEBUG
	cudaCheckCall(cuModuleLoad(&evolution_kernels, strcat(path, "kernels_debug.ptx")));
	#else
	cudaCheckCall(cuModuleLoad(&evolution_kernels, strcat(path, "kernels.ptx")));
	#endif // DEBUG

	/* Functions load */
	CUfunction initGPU, initCells, step1, cleanCells, placeFood, step2, recount, step3;
	cudaCheckCall(cuModuleGetFunction(&initGPU, evolution_kernels, "initGPU"));
	cudaCheckCall(cuModuleGetFunction(&initCells, evolution_kernels, "initCells"));
	cudaCheckCall(cuModuleGetFunction(&step1, evolution_kernels, "step1"));
	cudaCheckCall(cuModuleGetFunction(&cleanCells, evolution_kernels, "cleanCells"));
	cudaCheckCall(cuModuleGetFunction(&placeFood, evolution_kernels, "placeFood"));
	cudaCheckCall(cuModuleGetFunction(&step2, evolution_kernels, "step2"));
	cudaCheckCall(cuModuleGetFunction(&recount, evolution_kernels, "recount"));
	cudaCheckCall(cuModuleGetFunction(&step3, evolution_kernels, "step3"));

	void *initGPU_args[8] = { &culture_d, &culture_cells_d, &rows, &columns, &cells_d1, &cells_d2, &num_cells, &stats_d };
	cudaCheckCall(cuLaunchKernel(initGPU, 1, 1, 1, 1, 1, 1, 0, 0, initGPU_args, NULL));
	void *initCells_args[1] = { &random_seqs_d };
	cudaCheckCall(cuLaunchKernel(initCells, BLOCK_C, THREADS, 0, 0, initCells_args, NULL));

	/* 4. Simulation */
	int iter;
	int max_food_int = max_food * PRECISION;

	/* Food generation and placement variables and structures */
	int num_new_sources = (int)(rows * columns * food_density);
	int	num_new_sources_spot = food_spot_active ? (int)(food_spot_size_rows * food_spot_size_cols * food_spot_density) : 0;
	int max_new_sources = num_new_sources + num_new_sources_spot;
	food_t *food_to_place;
	cudaCheckCall(cuMemAllocHost((void **)&food_to_place, sizeof(food_t) * (size_t)max_new_sources));
	CUdeviceptr food_to_place_d;
	cudaCheckCall(cuMemAlloc(&food_to_place_d, sizeof(food_t) * (size_t)max_new_sources));

	/* Food for first iteration. */
	for (i=0; i<num_new_sources; i++) {
		food_to_place[i].pos = int_urand48(rows, food_random_seq)*columns;
		food_to_place[i].pos += int_urand48(columns, food_random_seq);
		food_to_place[i].food = int_urand48(food_level * PRECISION, food_random_seq);
	}
	// In the special food spot
	if (food_spot_active) {
		for (; i<max_new_sources; i++) {
			food_to_place[i].pos = (food_spot_row + int_urand48(food_spot_size_rows, food_spot_random_seq))*columns;
			food_to_place[i].pos += food_spot_col + int_urand48(food_spot_size_cols, food_spot_random_seq);
			food_to_place[i].food = int_urand48(food_spot_level * PRECISION, food_spot_random_seq);
		}
	}

	/* First 10 iterations are different. */
	for(iter=0; iter < min(max_iter, 10) && sim_stat.history_max_food <= max_food_int; iter++) {

		cudaCheckCall(cuLaunchKernel(step1, BLOCK_C, THREADS, sizeof(int) * THREADS1D, 0, NULL, NULL));

		cudaCheckCall(cuMemcpyHtoD(food_to_place_d, food_to_place, sizeof(food_t) * (size_t)max_new_sources));

		/* Steps of the simulation */
		void *placeFood_args[2] = { &food_to_place_d, &max_new_sources };
		cudaCheckCall(cuLaunchKernel(placeFood, BLOCK_F, THREADS, 0, 0, placeFood_args, NULL));
		cudaCheckCall(cuLaunchKernel(step2, BLOCK_F, THREADS, 0, 0, NULL, NULL));
		cudaCheckCall(cuLaunchKernel(recount, 1, 1, 1, 1, 1, 1, 0, 0, NULL, NULL));
		cudaCheckCall(cuLaunchKernel(step3, BLOCK_P, THREADS, sizeof(int) * THREADS1D, 0, NULL, NULL));

		/* 4.1. Spreading new food */
		// Across the whole culture
		for (i=0; i<num_new_sources; i++) {
			food_to_place[i].pos = int_urand48(rows, food_random_seq)*columns;
			food_to_place[i].pos += int_urand48(columns, food_random_seq);
			food_to_place[i].food = int_urand48(food_level * PRECISION, food_random_seq);
		}
		// In the special food spot
		if (food_spot_active) {
			for (; i<max_new_sources; i++) {
				food_to_place[i].pos = (food_spot_row + int_urand48(food_spot_size_rows, food_spot_random_seq))*columns;
				food_to_place[i].pos += food_spot_col + int_urand48(food_spot_size_cols, food_spot_random_seq);
				food_to_place[i].food = int_urand48(food_spot_level * PRECISION, food_spot_random_seq);
			}
		}

		Statistics prev_stats = sim_stat;		
		cudaCheckCall(cuMemcpyDtoH(&sim_stat, stats_d, sizeof(Statistics)));	// For some reason, async calls don't get synchronized
		                                                                    	// even using the cudaStreamSychronize(alt) at the end.

		/* Recalculate number of cells alive */
		if (iter > 0)	// Needed because prev_stats is all zeroes initialy.
			num_cells_alive += (sim_stat.history_total_cells - prev_stats.history_total_cells) - (sim_stat.history_dead_cells - prev_stats.history_dead_cells);
#ifdef DEBUG
		/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
		//print_status(iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat);
#endif // DEBUG
	}

	for(; iter<max_iter && sim_stat.history_max_food <= max_food_int && num_cells_alive > 0; iter++) {

		cudaCheckCall(cuLaunchKernel(step1, BLOCK_C, THREADS,  sizeof(int) * THREADS1D, 0, NULL, NULL));
		cudaCheckCall(cuLaunchKernel(cleanCells, BLOCK_C, THREADS, 0, 0, NULL, NULL));

		cudaCheckCall(cuMemcpyHtoD(food_to_place_d, food_to_place, sizeof(food_t) * (size_t)max_new_sources));

		/* Steps of the simulation */
		void *placeFood_args[2] = { &food_to_place_d, &max_new_sources };
		cudaCheckCall(cuLaunchKernel(placeFood, BLOCK_F, THREADS, 0, 0, placeFood_args, NULL));
		cudaCheckCall(cuLaunchKernel(step2, BLOCK_F, THREADS, 0, 0, NULL, NULL));
		cudaCheckCall(cuLaunchKernel(recount, 1, 1, 1, 1, 1, 1, 0, 0, NULL, NULL));
		cudaCheckCall(cuLaunchKernel(step3, BLOCK_P, THREADS, sizeof(int) * THREADS1D, 0, NULL, NULL));

		/* 4.1. Spreading new food */
		// Across the whole culture
		for (i=0; i<num_new_sources; i++) {
			food_to_place[i].pos = int_urand48(rows, food_random_seq)*columns;
			food_to_place[i].pos += int_urand48(columns, food_random_seq);
			food_to_place[i].food = int_urand48(food_level * PRECISION, food_random_seq);
		}
		// In the special food spot
		if (food_spot_active) {
			for (; i<max_new_sources; i++) {
				food_to_place[i].pos = (food_spot_row + int_urand48(food_spot_size_rows, food_spot_random_seq))*columns;
				food_to_place[i].pos += food_spot_col + int_urand48(food_spot_size_cols, food_spot_random_seq);
				food_to_place[i].food = int_urand48(food_spot_level * PRECISION, food_spot_random_seq);
			}
		}

		Statistics prev_stats = sim_stat;		
		cudaCheckCall(cuMemcpyDtoH(&sim_stat, stats_d, sizeof(Statistics)));	// For some reason, async calls don't get synchronized
																				// even using the cudaStreamSychronize(alt) at the end.

		num_cells_alive += (sim_stat.history_total_cells - prev_stats.history_total_cells) - (sim_stat.history_dead_cells - prev_stats.history_dead_cells);
#ifdef DEBUG
		/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
		//print_status(iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat);
#endif // DEBUG
	}

	cudaCheckCall(cuMemFree(culture_d));
	cudaCheckCall(cuMemFree(culture_cells_d));
	cudaCheckCall(cuStreamSynchronize(alt));

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

	// CUDA stop
	cuCtxDetach(context);

	/* 5. Stop global time */
	ttotal = cp_Wtime() - ttotal;

#ifdef DEBUG
	printf("List of cells at the end of the simulation: %d\n\n", num_cells);
	for(i=0; i<num_cells; i++) {
		printf("Cell %d, Alive: %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
				i,
				cells[i].alive,
				(float)cells[i].pos_row / PRECISION, 
				(float)cells[i].pos_col / PRECISION, 
				(float)cells[i].mov_row / PRECISION, 
				(float)cells[i].mov_col / PRECISION, 
				(float)cells[i].choose_mov[0] / PRECISION, 
				(float)cells[i].choose_mov[1] / PRECISION, 
				(float)cells[i].choose_mov[2] / PRECISION, 
				(float)cells[i].storage / PRECISION,
				cells[i].age);
	}
#endif // DEBUG

	/* 6. Output for leaderboard */
	printf("\n");
	/* 6.1. Total computation time */
	printf("Time: %lf\n", ttotal);

	/* 6.2. Results: Number of iterations and other statistics */
	printf("Result: %d, ", iter);
	printf("%d, %d, %d, %d, %d, %d, %d, %f\n", 
		num_cells_alive, 
		sim_stat.history_total_cells, 
		sim_stat.history_dead_cells, 
		sim_stat.history_max_alive_cells, 
		sim_stat.history_max_new_cells, 
		sim_stat.history_max_dead_cells, 
		sim_stat.history_max_age,
		(float)sim_stat.history_max_food / PRECISION
	);

	/* 7. Free resources */	
	free(culture);
	free(culture_cells);
	free(cells);

	/* 8. End */
	return 0;
}
