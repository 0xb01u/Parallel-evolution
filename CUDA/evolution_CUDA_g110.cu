/*
 * Simplified simulation of life evolution
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2019/2020
 *
 * v1.5
 *
 * CHANGES:
 * 1) Float values have been substituted by fixed point arithmetics 
 *	using integers. To simplify, the fixed point arithmetics are done 
 *	with PRECISION in base 10. See precision constant in int_float.h
 * 2) It uses a portable approximation to trigonometric functions using
 *	Taylor polynomials. 
 * 3) nrand48 function has been extracted from glibc source code and 
 *	its internal API simplified to allow its use in the GPU.
 *
 * (c) 2020, Arturo Gonzalez Escribano
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
#include<cputils.h>
#include<cuda.h>
#include<int_float.h>

/* 
 * Constants: Converted to fixed point with the given PRECISION
 */
#define ENERGY_NEEDED_TO_LIVE		PRECISION / 10	// Equivalent to 0.1
#define ENERGY_NEEDED_TO_MOVE		PRECISION	// Equivalent to 1.0
#define ENERGY_SPENT_TO_LIVE		PRECISION / 5	// Equivalent to 0.2
#define ENERGY_SPENT_TO_MOVE		PRECISION	// Equivalent to 1.0
#define ENERGY_NEEDED_TO_SPLIT		PRECISION * 20	// Equivalent to 20.0


/* Structure to store data of a cell */
typedef struct {
	int pos_row, pos_col;		// Position
	int mov_row, mov_col;		// Direction of movement
	int choose_mov[3];		// Genes: Probabilities of 0 turning-left; 1 advance; 2 turning-right
	int storage;			// Food/Energy stored
	int age;			// Number of steps that the cell has been alive
	unsigned short random_seq[3];	// Status value of its particular random sequence
	bool alive;			// Flag indicating if the cell is still alive
} Cell;


/* Structure for simulation statistics */
typedef struct {
	int history_total_cells;	// Accumulated number of cells created
	int history_dead_cells;		// Accumulated number of dead cells
	int history_max_alive_cells;	// Maximum number of cells alive in a step
	int history_max_new_cells;	// Maximum number of cells created in a step
	int history_max_dead_cells;	// Maximum number of cells died in a step
	int history_max_age;		// Maximum age achieved by a cell
	int history_max_food;		// Maximum food level in a position of the culture
} Statistics;


/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 *	USE THIS SPACE FOR YOUR KERNEL OR DEVICE FUNTIONS
 *
 */

#include "taylor_trig.h"
#include "glibc_nrand48.h"

/*
 * Get an uniformly distributed random number between 0 and max
 * It uses glibc_nrand, that returns a number between 0 and 2^31
 */
#define int_urand48( max, seq )	(int)( (long)(max) * glibc_nrand48( seq ) / 2147483648 )

/* 
 * Macro function to simplify accessing with two coordinates to a flattened array
 * 	This macro-function can be modified by the students if needed
 *
 */
#define accessMat( arr, exp1, exp2 )	arr[ (int)(exp1) * columns + (int)(exp2) ]

/*
 * Function: Choose a new direction of movement for a cell
 * 	This function can be changed and/or optimized by the students
 */
__host__ __device__ void cell_new_direction( Cell *cell ) {
	int angle = int_urand48( INT_2PI, cell->random_seq );
	cell->mov_row = taylor_sin( angle );
	cell->mov_col = taylor_cos( angle );
}

/*
 * Function: Mutation of the movement genes on a new cell
 * 	This function can be changed and/or optimized by the students
 */
__host__ __device__ void cell_mutation( Cell *cell ) {
	/* 1. Select which genes change:
	 	0 Left grows taking part of the Advance part
	 	1 Advance grows taking part of the Left part
	 	2 Advance grows taking part of the Right part
	 	3 Right grows taking part of the Advance part
	*/
	int mutation_type = int_urand48( 4, cell->random_seq );
	/* 2. Select the amount of mutation (up to 50%) */
	int mutation_percentage = int_urand48( PRECISION / 2, cell->random_seq );
	/* 3. Apply the mutation */
	int mutation_value;
	switch( mutation_type ) {
		case 0:
			mutation_value = intfloatMult( cell->choose_mov[1] , mutation_percentage );
			cell->choose_mov[1] -= mutation_value;
			cell->choose_mov[0] += mutation_value;
			break;
		case 1:
			mutation_value = intfloatMult( cell->choose_mov[0] , mutation_percentage );
			cell->choose_mov[0] -= mutation_value;
			cell->choose_mov[1] += mutation_value;
			break;
		case 2:
			mutation_value = intfloatMult( cell->choose_mov[2] , mutation_percentage );
			cell->choose_mov[2] -= mutation_value;
			cell->choose_mov[1] += mutation_value;
			break;
		case 3:
			mutation_value = intfloatMult( cell->choose_mov[1] , mutation_percentage );
			cell->choose_mov[1] -= mutation_value;
			cell->choose_mov[2] += mutation_value;
			break;
	}
	/* 4. Correct potential precision problems */
	cell->choose_mov[2] = PRECISION - cell->choose_mov[1] - cell->choose_mov[0];
}

/*
 * CUDA block reduction
 * Inputs: 
 *	Device pointer to an array of int of any size
 *	Size of the array
 *	Device pointer to an int to store the result
 * 
 * Launching parameters:
 *	One-dimesional grid of any size
 *	Any valid block size
 *	Dynamic shared memory size equal to: sizeof(int) * block size
 *
 * (c) 2020, Arturo Gonzalez-Escribano
 * Simplification for an assignment in a Parallel Computing course,
 * Computing Engineering Degree, Universidad de Valladolid
 * Academic year 2019/2020
 */
__device__ void reductionMax(int *array, int size, int *result)
{
	int tid = threadIdx.x;
	int gid = tid + blockIdx.x * blockDim.x;

	extern __shared__ int buffer[ ];
	if ( gid < size ) { 
		buffer[ tid ] = array[ gid ];
	}
	else buffer[ tid ] = 0;
	__syncthreads();

	for( int step=blockDim.x/2; step>=1; step /= 2 ) {
		if ( tid < step )
			if ( buffer[ tid ] < buffer[ tid + step ] )
				buffer[ tid ] = buffer[ tid + step ];
		if ( step > 32 )
			__syncthreads();
	}

	if ( tid == 0 )
		atomicMax( result, buffer[0] );
}

/* ===================== WE START HERE ===================== */

/*
 * Copy-paste of the other reductionMax, but for cell age.
 *
 */
__device__ void reductionMax(Cell* array, int size, int *result)
{
	int tid = threadIdx.x;
	int gid = tid + blockIdx.x * blockDim.x;

	extern __shared__ int buffer[ ];
	if ( gid < size ) { 
		buffer[ tid ] = array[ gid ].age;
	}
	else buffer[ tid ] = 0;
	__syncthreads();

	for( int step=blockDim.x/2; step>=1; step /= 2 ) {
		if ( tid < step )
			if ( buffer[ tid ] < buffer[ tid + step ] )
				buffer[ tid ] = buffer[ tid + step ];
		if ( step > 32 )
			__syncthreads();
	}

	if ( tid == 0 )
		atomicMax( result, buffer[0] );
}

/*
 * Function to atomically add a value to a short-typed variable.
 *
 */
__device__ int atomicAdd(short *pos, short inc)
{
	if (((long)pos & 0x3) != 0)
	{
		pos--;
		return (short)((atomicAdd((int *)pos, (int)inc * 0x10000) >> 16) & 0xffff);
	}
	return (short)(atomicAdd((int *)pos, (int)inc) & 0xffff);
}

/*
 * Ancillary type for random food generation.
 *
 */
typedef struct {
	int food;
	int pos;
} food_t;

/*
 * Maximum of two or three values
 *
 */
#define max(x, y) (x > y ? x : y)
#define max3(x, y, z) max(max(x, y), z)

/*
 * Minimum of two or values
 *
 */
#define min(x, y) (x < y ? x : y)

/*
 * Global identifier for a device thread
 *
 */
#define GLOBAL_ID threadIdx.x + blockIdx.x * blockDim.x

/*
 * Global device variables.
 * These are the same we'd work with on the CPU.
 *  The names are kept the same so the code is more easily legible.
 *  They are global so the kernels take less arguments.
 */
__device__ int rows = 0;
__device__ int columns = 0;
__device__ int num_cells = 0;
__device__ int *culture = NULL;
__device__ short *culture_cells = NULL;
__device__ Cell *cells = NULL;
__device__ Cell *cells_aux = NULL;
__device__ Statistics *sim_stat;
__device__ int num_cells_alive = 0;
__device__ int step_dead_cells = 0;
__device__ int step_new_cells = 0;
__device__ int free_position = 0;

/*
 * Initialize global device variables.
 *
 */
__global__ void initGPU(int *culture_d, short *culture_cells_d, int rows_d, int columns_d, Cell *cells_d1, Cell *cells_d2, int num_cells_d, Statistics *stats)
{
	rows = rows_d;
	columns = columns_d;
	num_cells = num_cells_d;
	culture = culture_d;
	culture_cells = culture_cells_d;
	cells = cells_d1;
	cells_aux = cells_d2;

	num_cells_alive = num_cells;

	sim_stat = stats;

	sim_stat->history_total_cells = num_cells;
	sim_stat->history_dead_cells = 0;
	sim_stat->history_max_alive_cells = num_cells;
	sim_stat->history_max_new_cells = 0;
	sim_stat->history_max_dead_cells = 0;
	sim_stat->history_max_age = 0;
	sim_stat->history_max_food = 0.0f;
}

/*
 * Initialize cell list on the device.
 *
 */
__global__ void initCells(unsigned short *random_seqs_d)
{
	int gid = GLOBAL_ID;

	if (gid >= num_cells) return;

	Cell *my_cell = &cells[gid];

	my_cell->random_seq[0] = random_seqs_d[3*gid];
	my_cell->random_seq[1] = random_seqs_d[3*gid + 1];
	my_cell->random_seq[2] = random_seqs_d[3*gid + 2];

	my_cell->alive = true;
	// Initial age: Between 1 and 20 
	my_cell->age = 1 + int_urand48( 19, my_cell->random_seq );
	// Initial storage: Between 10 and 20 units
	my_cell->storage = 10 * PRECISION + int_urand48( 10 * PRECISION, my_cell->random_seq );
	// Initial position: Anywhere in the culture arena
	my_cell->pos_row = int_urand48( rows * PRECISION, my_cell->random_seq );
	my_cell->pos_col = int_urand48( columns * PRECISION, my_cell->random_seq );
	// Movement direction: Unity vector in a random direction
	cell_new_direction( my_cell );
	// Movement genes: Probabilities of advancing or changing direction: The sum should be 1.00
	my_cell->choose_mov[0] = PRECISION / 3;
	my_cell->choose_mov[2] = PRECISION / 3;
	my_cell->choose_mov[1] = PRECISION - my_cell->choose_mov[0] - my_cell->choose_mov[2];
}

#ifdef DEBUG
__device__ void print_statusGPU( int rows, int columns, int *culture, int num_cells, Cell *cells, int num_cells_alive, Statistics sim_stat );
#endif // DEBUG

/*
 * Cell movementes.
 *  Section 4.3 of the simulation.
 */
__global__ void step1()
{
	int gid = GLOBAL_ID;

	Cell *my_cell = &cells[gid];

	/* 4.3. Cell movements */
	if (gid < num_cells)
	{
		my_cell->age ++;

		/* 4.3.1. Check if the cell has the needed energy to move or keep alive */
		if ( my_cell->storage < ENERGY_NEEDED_TO_LIVE ) {
			// Cell has died
			my_cell->alive = false;
			atomicAdd(&step_dead_cells, 1);
		}

		if (my_cell->alive)
		{
			if ( my_cell->storage < ENERGY_NEEDED_TO_MOVE ) {
				// Almost dying cell, it cannot move, only if enough food is dropped here it will survive
				my_cell->storage -= ENERGY_SPENT_TO_LIVE;
			}
			else {
				// Consume energy to move
				my_cell->storage -= ENERGY_SPENT_TO_MOVE;
					
				/* 4.3.2. Choose movement direction */
				int prob = int_urand48( PRECISION, my_cell->random_seq );
				if ( prob < my_cell->choose_mov[0] ) {
					// Turn left (90 degrees)
					int tmp = my_cell->mov_col;
					my_cell->mov_col = my_cell->mov_row;
					my_cell->mov_row = -tmp;
				}
				else if ( prob >= my_cell->choose_mov[0] + my_cell->choose_mov[1] ) {
					// Turn right (90 degrees)
					int tmp = my_cell->mov_row;
					my_cell->mov_row = my_cell->mov_col;
					my_cell->mov_col = -tmp;
				}
				// else do not change the direction
				
				/* 4.3.3. Update position moving in the choosen direction*/
				my_cell->pos_row += my_cell->mov_row;
				my_cell->pos_col += my_cell->mov_col;
				// Periodic arena: Left/Rigth edges are connected, Top/Bottom edges are connected
				if ( my_cell->pos_row < 0 ) my_cell->pos_row += rows * PRECISION;
				if ( my_cell->pos_row >= rows * PRECISION) my_cell->pos_row -= rows * PRECISION;
				if ( my_cell->pos_col < 0 ) my_cell->pos_col += columns * PRECISION;
				if ( my_cell->pos_col >= columns * PRECISION) my_cell->pos_col -= columns * PRECISION;
			}
			/* 4.3.4. Annotate that there is one more cell in this culture position */
			short *pos = &accessMat( culture_cells, my_cell->pos_row / PRECISION, my_cell->pos_col / PRECISION );
			int inc = 0x1;
			if (((long)pos) % 4 != 0)
			{
				pos -= 1;
				inc = 0x10000;

			}
			atomicAdd((int *)pos, inc);
		}

	} // End cell movements

	// Statistics: Max age of a cell in the simulation history
	reductionMax(cells, num_cells, &sim_stat->history_max_age);
}

/*
 * Function to clean dead cells from the list of cells.
 *
 */
__global__ void cleanCells()
{
	int gid = GLOBAL_ID;

	if (step_dead_cells > 0 && gid < num_cells)
	{
		Cell *my_cell = &cells[gid];
		if ( my_cell->alive ) {
			cells_aux[atomicAdd(&free_position, 1)] = *my_cell;
		}
		if (gid == 0) num_cells_alive -= step_dead_cells;
	}
}

/*
 * Function to place the randomly-generated in the host food on the
 * device culture structure.
 */
__global__ void placeFood(food_t *food, int num_food)
{
	int gid = GLOBAL_ID;

	if (gid == 0 && step_dead_cells > 0)
	{
		Cell *tmp = cells;
		cells = cells_aux;
		cells_aux = tmp;
	}

	if (gid < num_food) atomicAdd(&culture[food[gid].pos], food[gid].food);
}

/*
 * Cell actions.
 *  Section 4.4 of the simulation.
 */
__global__ void step2()
{
	int gid = GLOBAL_ID;

	Cell *my_cell = &cells[gid];

	/* 4.4.1. Food harvesting */
	if (gid < num_cells_alive)
	{
		int food = accessMat( culture, my_cell->pos_row / PRECISION, my_cell->pos_col / PRECISION );
		int count = accessMat( culture_cells, my_cell->pos_row / PRECISION, my_cell->pos_col / PRECISION );
		int my_food = food / count;
		my_cell->storage += my_food;

		/* 4.4.2. Split cell if the conditions are met: Enough maturity and energy */
		if ( my_cell->age > 30 && my_cell->storage > ENERGY_NEEDED_TO_SPLIT ) {

			// Split: Create new cell
			int pos = atomicAdd(&step_new_cells, 1) + num_cells_alive;

			// Split energy stored and update age in both cells
			my_cell->storage /= 2;
			my_cell->age = 1;

			// New cell is a copy of parent cell
			Cell *new_cell = &cells[pos];
			*new_cell = *my_cell;

			// Random seed for the new cell, obtained using the parent random sequence
			new_cell->random_seq[0] = (unsigned short)glibc_nrand48( my_cell->random_seq );
			new_cell->random_seq[1] = (unsigned short)glibc_nrand48( my_cell->random_seq );
			new_cell->random_seq[2] = (unsigned short)glibc_nrand48( my_cell->random_seq );

			// Both cells start in random directions
			cell_new_direction( my_cell );
			cell_new_direction( new_cell );
		
			// Mutations of the movement genes in both cells
			cell_mutation( my_cell );
			cell_mutation( new_cell );
		} // End cell actions
	}

}

/*
 * Function to correctly recount cells and statistics, on one thread.
 *
 */
__global__ void recount()
{
	// Hello there, Yuri!
	asm(
		"add.s32	%0, %1, %2;"
		: "=r"(num_cells_alive)
		: "r"(num_cells_alive),
		  "r"(step_new_cells)
	);
	asm(
		"add.s32	%0, %6, %5;\n\t"
		"max.s32	%1, %7, %5;\n\t"
		"add.s32	%2, %9, %8;\n\t"
		"max.s32	%3, %10, %8;\n\t"
		"max.s32	%4, %12, %11;"
		: "=r"(sim_stat->history_total_cells),
		  "=r"(sim_stat->history_max_new_cells),
		  "=r"(sim_stat->history_dead_cells),
		  "=r"(sim_stat->history_max_dead_cells),
		  "=r"(sim_stat->history_max_alive_cells)
		: "r"(step_new_cells),
		  "r"(sim_stat->history_total_cells),
		  "r"(sim_stat->history_max_new_cells),
		  "r"(step_dead_cells),
		  "r"(sim_stat->history_dead_cells),
		  "r"(sim_stat->history_max_dead_cells),
		  "r"(num_cells_alive),
		  "r"(sim_stat->history_max_alive_cells)
	);
	asm(
		"mov.s32	%0, 0;\n\t"
		"mov.s32	%1, 0;\n\t"
		"mov.s32	%2,	0;"
		: "=r"(step_dead_cells),
		  "=r"(step_new_cells),
		  "=r"(free_position)
	);
	asm(
		"mov.s32	%0, %1;"
		: "=r"(num_cells)
		: "r"(num_cells_alive)
	);
	// Not gonna lie, we couldn't get it to work in just one asm() call.
}

/*
 * Food decrease and statistics update.
 *  Sections 4.5 and 4.8 of the simulation.
 */
__global__ void step3()
{
	int gid = GLOBAL_ID;

	/* 4.8. Decrease non-harvested food */
	if (gid < rows*columns)
	{
		if (culture_cells[gid] != 0) culture[gid] = 0;
		else culture[gid] -= culture[gid] / 20;
		/* 4.2. Prepare ancillary data structures */	
		/* 4.2.1. Clear ancillary structure of the culture to account alive cells in a position after movement */
		culture_cells[gid] = 0;
	}
	reductionMax(culture, rows*columns, &sim_stat->history_max_food);


#ifdef DEBUG
	if (gid == 0)
	{
		/* In case someone wants to print the debug information for each iteration, this is the place. */
		print_statusGPU(rows, columns, culture, num_cells, cells, num_cells_alive, *sim_stat);
	}
#endif // DEBUG
}

#ifdef DEBUG
/* 
 * Function: Print the current state of the simulation, with verbose information (exact storage and food).
 *  Reconfigured to work on a device thread.
 */
__device__ void print_statusGPU( int rows, int columns, int *culture, int num_cells, Cell *cells, int num_cells_alive, Statistics sim_stat ) {

	int i,j;

	printf("+");
	for( j=0; j<columns; j++ ) printf("---");
	printf("+\n");
	for( i=0; i<rows; i++ ) {
		printf("|");
		for (j = 0; j < columns; j++)
        {
            int t;
            int counter = 0;
            int n = 0;
            for (t = 0; t < num_cells; t++)
            {
                int row = (int)(cells[t].pos_row / PRECISION);
                int col = (int)(cells[t].pos_col / PRECISION);
                if (cells[t].alive && row == i && col == j)
                {
                	n++;
                    counter += cells[t].storage;
                }
            }
            if (counter > 0)
            	if (n > 1)
            		printf("(%06d)%d", counter, n);
            	else
                	printf("(%06d)", counter);
            else
                printf(" %06d ", (accessMat(culture, i, j)));
        }
		printf("|\n");
	}
	printf("+");
	for( j=0; j<columns; j++ ) printf("---");
	printf("+\n");
	printf("Num_cells_alive: %04d\nHistory( Cells: %04d, Dead: %04d, Max.alive: %04d, Max.new: %04d, Max.dead: %04d, Max.age: %04d, Max.food: %6f )\n\n", 
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
#endif // DEBUG

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

#ifdef DEBUG
/* 
 * Function: Print the current state of the simulation 
 */
void print_status( int iteration, int rows, int columns, int *culture, int num_cells, Cell *cells, int num_cells_alive, Statistics sim_stat ) {
	/* 
	 * You don't need to optimize this function, it is only for pretty printing and debugging purposes.
	 * It is not compiled in the production versions of the program.
	 * Thus, it is never used when measuring times in the leaderboard
	 */
	int i,j;

	printf("Iteration: %d\n", iteration );
	printf("+");
	for( j=0; j<columns; j++ ) printf("---");
	printf("+\n");
	for( i=0; i<rows; i++ ) {
		printf("|");
		for( j=0; j<columns; j++ ) {
			char symbol;
			if ( accessMat( culture, i, j ) >= 20 * PRECISION ) symbol = '+';
			else if ( accessMat( culture, i, j ) >= 10 * PRECISION ) symbol = '*';
			else if ( accessMat( culture, i, j ) >= 5 * PRECISION ) symbol = '.';
			else symbol = ' ';

			int t;
			int counter = 0;
			for( t=0; t<num_cells; t++ ) {
				int row = (int)(cells[t].pos_row / PRECISION);
				int col = (int)(cells[t].pos_col / PRECISION);
				if ( cells[t].alive && row == i && col == j ) {
					counter ++;
				}
			}
			if ( counter > 9 ) printf("(M)" );
			else if ( counter > 0 ) printf("(%1d)", counter );
			else printf(" %c ", symbol );
		}
		printf("|\n");
	}
	printf("+");
	for( j=0; j<columns; j++ ) printf("---");
	printf("+\n");
	printf("Num_cells_alive: %04d\nHistory( Cells: %04d, Dead: %04d, Max.alive: %04d, Max.new: %04d, Max.dead: %04d, Max.age: %04d, Max.food: %6f )\n\n", 
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
void show_usage( char *program_name ) {
	fprintf(stderr,"Usage: %s ", program_name );
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
		show_usage( argv[0] );
		exit( EXIT_FAILURE );
	}

	/* 1.2. Read culture sizes, maximum number of iterations */
	rows = atoi( argv[1] );
	columns = atoi( argv[2] );
	max_iter = atoi( argv[3] );

	/* 1.3. Food data */
	max_food = atof( argv[4] );
	food_density = atof( argv[5] );
	food_level = atof( argv[6] );

	/* 1.4. Read random sequences initializer */
	for( i=0; i<3; i++ ) {
		init_random_seq[i] = (unsigned short)atoi( argv[7+i] );
	}

	/* 1.5. Read number of cells */
	num_cells = atoi( argv[10] );

	/* 1.6. Read special food spot */
	if (argc > 11 ) {
		if ( argc < 17 ) {
			fprintf(stderr, "-- Error in number of special-food-spot arguments in the command line\n\n");
			show_usage( argv[0] );
			exit( EXIT_FAILURE );
		}
		else {
			food_spot_active = true;
			food_spot_row = atoi( argv[11] );
			food_spot_col = atoi( argv[12] );
			food_spot_size_rows = atoi( argv[13] );
			food_spot_size_cols = atoi( argv[14] );
			food_spot_density = atof( argv[15] );
			food_spot_level = atof( argv[16] );

			// Check non-used trailing arguments
			if ( argc > 17 ) {
				fprintf(stderr, "-- Error: too many arguments in the command line\n\n");
				show_usage( argv[0] );
				exit( EXIT_FAILURE );
			}
		}
	}

#ifdef DEBUG
	/* 1.7. Print arguments */
	printf("Arguments, Rows: %d, Columns: %d, max_iter: %d\n", rows, columns, max_iter);
	printf("Arguments, Max.food: %f, Food density: %f, Food level: %f\n", max_food, food_density, food_level);
	printf("Arguments, Init Random Sequence: %hu,%hu,%hu\n", init_random_seq[0], init_random_seq[1], init_random_seq[2]);
	if ( food_spot_active ) {
		printf("Arguments, Food_spot, pos(%d,%d), size(%d,%d), Density: %f, Level: %f\n",
			food_spot_row, food_spot_col, food_spot_size_rows, food_spot_size_cols, food_spot_density, food_spot_level );
	}
	printf("Initial cells: %d\n", num_cells );
#endif // DEBUG


	/* 1.8. Initialize random sequences for food dropping */
	for( i=0; i<3; i++ ) {
		food_random_seq[i] = (unsigned short)glibc_nrand48( init_random_seq );
		food_spot_random_seq[i] = (unsigned short)glibc_nrand48( init_random_seq );
	}

	/* 1.9. Initialize random sequences of cells */
	cells = (Cell *)malloc( sizeof(Cell) * (size_t)num_cells );
	if ( cells == NULL ) {
		fprintf(stderr,"-- Error allocating: %d cells\n", num_cells );
		exit( EXIT_FAILURE );
	}
	for( i=0; i<num_cells; i++ ) {
		// Initialize the cell ramdom sequences
		for( j=0; j<3; j++ ) 
			cells[i].random_seq[j] = (unsigned short)glibc_nrand48( init_random_seq );
	}


#ifdef DEBUG
	/* 1.10. Print random seed of the initial cells */
	/*
	printf("Initial cells random seeds: %d\n", num_cells );
	for( i=0; i<num_cells; i++ )
		printf("\tCell %d, Random seq: %hu,%hu,%hu\n", i, cells[i].random_seq[0], cells[i].random_seq[1], cells[i].random_seq[2] );
	*/
#endif // DEBUG


	// CUDA start
	cudaSetDevice(0);
	cudaDeviceSynchronize();

	/* 2. Start global timer */
	double ttotal = cp_Wtime();

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */

#include "cuda_check.h"
/*
 * Simple macro function to check errors on kernel executions.
 * To use alongside cuda_check.h
 *
 */
#define cudaCheckKernel(kernel) { \
	kernel; \
	cudaCheckLast(); \
}

/*
 * Block and thread sizes for kernel executions.
 *
 */
#define THREADS 128
#define BLOCK (max(rows*columns, num_cells_alive)/THREADS + 1)
#define BLOCK_F (max3(rows*columns, num_cells_alive, max_new_sources)/THREADS + 1)	/* "Food kernels" */
#define BLOCK_C (num_cells_alive)/ THREADS + 1	/* "Cell kernels" */
#define BLOCK_P (rows*columns)/THREADS + 1		/* "Culture kernels" */


	/* 3. Initialize culture surface and initial cells */
	culture = NULL;
	culture_cells = NULL;

	/* Device equivalents */
	int *culture_d;
	short *culture_cells_d;
	cudaCheckCall(cudaMalloc(&culture_d, sizeof(int) * (size_t)rows * (size_t)columns));
	cudaCheckCall(cudaMalloc(&culture_cells_d, sizeof(short) * (size_t)rows * (size_t)columns));

	/* Set both surfaces to 0 */
	cudaMemset(culture_d, 0, sizeof(int) * (size_t)rows * (size_t)columns);
	cudaMemset(culture_cells_d, 0, sizeof(short) * (size_t)rows * (size_t)columns);

	/* CUDA streams */
	cudaStream_t alt;
	cudaCheckCall(cudaStreamCreate(&alt));

	/* Copy random cell seeds to device */
	unsigned short *random_seqs;
	cudaCheckCall(cudaMallocHost(&random_seqs, sizeof(unsigned short) * 3 * num_cells));
	unsigned short *random_seqs_d;
	cudaCheckCall(cudaMalloc(&random_seqs_d, sizeof(unsigned short) * 3 * num_cells));

	for (i = 0; i < num_cells; i++)
	{
		random_seqs[3*i] = cells[i].random_seq[0];
		random_seqs[3*i + 1] = cells[i].random_seq[1];
		random_seqs[3*i + 2] = cells[i].random_seq[2];
	}

	cudaCheckCall(cudaMemcpyAsync(random_seqs_d, random_seqs, sizeof(unsigned short) * 3 * num_cells, cudaMemcpyHostToDevice, alt));

	int num_cells_alive = num_cells;

	/* Device cell lists */
	/* They are assigned 2GiB of memory each, so they are never required to realloc. */
	Cell *cells_d1;
	cudaCheckCall(cudaMalloc(&cells_d1, (size_t) (1l << 31)));
	Cell *cells_d2;
	cudaCheckCall(cudaMalloc(&cells_d2, (size_t) (1l << 31)));
	/* Device statistics */
	Statistics *stats_d;
	cudaCheckCall(cudaMalloc(&stats_d, sizeof(Statistics)));

	/* Device intialization */
	initGPU<<<1, 1>>>(culture_d, culture_cells_d, rows, columns, cells_d1, cells_d2, num_cells, stats_d);
	initCells<<<BLOCK_C, THREADS>>>(random_seqs_d);

	/* 4. Simulation */
	int iter;
	int max_food_int = max_food * PRECISION;

	/* Food generation and placement variables and structures */
	int num_new_sources = (int)(rows * columns * food_density);
	int	num_new_sources_spot = food_spot_active ? (int)(food_spot_size_rows * food_spot_size_cols * food_spot_density) : 0;
	int max_new_sources = num_new_sources + num_new_sources_spot;
	food_t *food_to_place;
	cudaCheckCall(cudaMallocHost(&food_to_place, sizeof(food_t) * (size_t)max_new_sources));
	food_t *food_to_place_d;
	cudaCheckCall(cudaMalloc(&food_to_place_d, sizeof(food_t) * (size_t)max_new_sources));

	/* Food for first iteration. */
	for (i=0; i<num_new_sources; i++) {
		food_to_place[i].pos = int_urand48( rows, food_random_seq )*columns;
		food_to_place[i].pos += int_urand48( columns, food_random_seq );
		food_to_place[i].food = int_urand48( food_level * PRECISION, food_random_seq );
	}
	// In the special food spot
	if ( food_spot_active ) {
		for (; i<max_new_sources; i++) {
			food_to_place[i].pos = (food_spot_row + int_urand48( food_spot_size_rows, food_spot_random_seq ))*columns;
			food_to_place[i].pos += food_spot_col + int_urand48( food_spot_size_cols, food_spot_random_seq );
			food_to_place[i].food = int_urand48( food_spot_level * PRECISION, food_spot_random_seq );
		}
	}

	/* First 10 iterations are different. */
	for( iter=0; iter < min(max_iter, 10) && sim_stat.history_max_food <= max_food_int; iter++ ) {

		cudaCheckKernel((step1<<<BLOCK_C, THREADS, sizeof(int) * THREADS>>>()));

		cudaCheckCall(cudaMemcpy(food_to_place_d, food_to_place, sizeof(food_t) * (size_t)max_new_sources, cudaMemcpyHostToDevice));

		/* Steps of the simulation */
		cudaCheckKernel((placeFood<<<BLOCK_F, THREADS>>>(food_to_place_d, max_new_sources)));
		cudaCheckKernel((step2<<<BLOCK_F, THREADS>>>()));
		cudaCheckKernel((recount<<<1, 1>>>()));
		cudaCheckKernel((step3<<<BLOCK_P, THREADS, sizeof(int) * THREADS>>>()));

		/* 4.1. Spreading new food */
		// Across the whole culture
		for (i=0; i<num_new_sources; i++) {
			food_to_place[i].pos = int_urand48( rows, food_random_seq )*columns;
			food_to_place[i].pos += int_urand48( columns, food_random_seq );
			food_to_place[i].food = int_urand48( food_level * PRECISION, food_random_seq );
		}
		// In the special food spot
		if ( food_spot_active ) {
			for (; i<max_new_sources; i++) {
				food_to_place[i].pos = (food_spot_row + int_urand48( food_spot_size_rows, food_spot_random_seq ))*columns;
				food_to_place[i].pos += food_spot_col + int_urand48( food_spot_size_cols, food_spot_random_seq );
				food_to_place[i].food = int_urand48( food_spot_level * PRECISION, food_spot_random_seq );
			}
		}

		Statistics prev_stats = sim_stat;		
		cudaCheckCall((cudaMemcpyAsync(&sim_stat, stats_d, sizeof(Statistics), cudaMemcpyDeviceToHost, alt)));

		/* Recalculate number of cells alive */
		if (iter > 0)	// Needed because prev_stats is all zeroes initialy.
			num_cells_alive += (sim_stat.history_total_cells - prev_stats.history_total_cells) - (sim_stat.history_dead_cells - prev_stats.history_dead_cells);
#ifdef DEBUG
		/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
		//print_status( iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat );
#endif // DEBUG
	}

	for( ; iter<max_iter && sim_stat.history_max_food <= max_food_int && num_cells_alive > 0; iter++ ) {

		cudaCheckKernel((step1<<<BLOCK_C, THREADS, sizeof(int) * THREADS>>>()));
		cudaCheckKernel((cleanCells<<<BLOCK_C, THREADS>>>()));

		cudaCheckCall(cudaMemcpy(food_to_place_d, food_to_place, sizeof(food_t) * (size_t)max_new_sources, cudaMemcpyHostToDevice));

		/* Steps of the simulation */
		cudaCheckKernel((placeFood<<<BLOCK_F, THREADS>>>(food_to_place_d, max_new_sources)));
		cudaCheckKernel((step2<<<BLOCK_F, THREADS>>>()));
		cudaCheckKernel((recount<<<1, 1>>>()));
		cudaCheckKernel((step3<<<BLOCK_P, THREADS, sizeof(int) * THREADS>>>()));

		/* 4.1. Spreading new food */
		// Across the whole culture
		for (i=0; i<num_new_sources; i++) {
			food_to_place[i].pos = int_urand48( rows, food_random_seq )*columns;
			food_to_place[i].pos += int_urand48( columns, food_random_seq );
			food_to_place[i].food = int_urand48( food_level * PRECISION, food_random_seq );
		}
		// In the special food spot
		if ( food_spot_active ) {
			for (; i<max_new_sources; i++) {
				food_to_place[i].pos = (food_spot_row + int_urand48( food_spot_size_rows, food_spot_random_seq ))*columns;
				food_to_place[i].pos += food_spot_col + int_urand48( food_spot_size_cols, food_spot_random_seq );
				food_to_place[i].food = int_urand48( food_spot_level * PRECISION, food_spot_random_seq );
			}
		}

		Statistics prev_stats = sim_stat;		
		cudaCheckCall((cudaMemcpyAsync(&sim_stat, stats_d, sizeof(Statistics), cudaMemcpyDeviceToHost, alt)));

		num_cells_alive += (sim_stat.history_total_cells - prev_stats.history_total_cells) - (sim_stat.history_dead_cells - prev_stats.history_dead_cells);
#ifdef DEBUG
		/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
		//print_status( iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat );
#endif // DEBUG
	}

/* And here goes the ASCII art :D
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkxdddddddddddddddddddddddddddddddddddddddddddddxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkllloxOOOdllldOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxoloxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxdoldkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxc::cccc:odc:cccc::okOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdcdk0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK0kocxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOl,lxkKXXKklcd0XXXOxd::kOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlkkx0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKOxOxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOdooooooxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo;ckXXXXKK0kokXXKXX0;.;coOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkkxxxxxxkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlO0OkoloooooooooooooooooooooooooooooooooolldkO0xlxOOOOOOOOOOOOOOOOOOOOOOOOkxolc:::clloxkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxxl',clk0kxxo:lk0KxoOkoodlcokOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxoc:;;;:coxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKKOcckOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOko;o0KKxlxOOOOOOOOOOOOOOOOOOOOOOOOo:oxO0000kxl:dOOOOOOOOOOOOOOOkdoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooxOOOOOOOOOOOOOOOOOxc::lxOOO0xclolooox00OKXXXK0d':OOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkoldkkxxxkkdlokOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOxl:oOkkkkkl:lkOOOOOOOOOOOOOOxlodddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddlokOOOOOOOOOOOOOkl:lOXXKKKKXXKkxOOOOOod0KXK00d:cdOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkl:lxOOOOOko:lkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOx:dK00O00l:kOOOOOOOOOOOOOOOo:xXKOkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk0KKl:kOOOOOOOOOOOOOk,.oddx0XXKOdolxXXXXKOdlodocccdkOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOl:k0kkkOOccOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOx:dKKKKK0l:kOOOOOOOOOOOOOOOo:xOxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxdxxxxxxxxxxxxxxxoodxxxxxxxxxxxdx0l:xOOOOOOOOOOOOOkl;;;:::cc:..lOKKKKKKKx:::lkOkOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo:kKKKKKO:cOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOkdodkOOkdoxkOOOOOOOOOOOOOOOo:xxlOXKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK0OkkkkkkkkkkkOKKkx0KKKKKKKKKKKxlko:xOOOOxdol:::ldxxoccc:;;:llccdxxkkxxlcok0OOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxloOKKKOdldOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOo;;;:xOOOOOOOOOOOOOOOOOOo:xxlxOOOOOOOOOOOOOOOOkkkOOOOOOOOOOOOOxc............,okdokOOOOOOOOOOOdlOo:kOOkl::cx000d;',oO00kl::cdkdc;,.,;:dkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkdlllllokOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOkd:,;;,cdkOOOOOOOOOOOOOOOOo:xxlxOOOOOOOOOOOOOOOOkoxOOOOOOOOOOkd;.............. .'lkOOOOOOOOOOOOdlOo:kOOd':0XXXKKX0ox0KXKXXXo.;xxO0x:oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkd:,;,;dkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOxcldc;,;olckOOOOOOOOOOOOOOOd:xxlOKKKKKKKKKKKKKKKK0xOKKKKKKKKKKd'... ...............,oOKKKKKKKKKKxlOo:kOOOx:'ldkKXKOOolkKX0xkx;;c:cxd;oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxll:,:,;lldOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOx:oOko:;dl:kOOOOOOOOOOOOOOOo:xxlOKKKKKKKKKKKKKKKK0xkKKKKKKKK0d,....:'.....:'....... .oKKKKKKKKKKxlOo:xOo:::cccdOoloooooxOOokXKKX0d::llokOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOl;xkoc,cx:cOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKO:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOkxolll:,coxkOOOOOOOOOOOOOOOo:xxlxOkkkkkkkkkkkkkkkkxkkkkkkkkOk, ..;xOkl,':xOl',ll'.. .;dkOkkkkkkOolko;llldxOXXXX0kddddddxxx0XXXXXKO:;dc:kOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxlloll;;lldOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOkkkkl;dOOOOOOOOOOOOOOOOOOo:xxlOK00KKKKKKKKKKKKKKKKKKKKK0KOl....c0000OOO00OOOO0d'',..:OK0KKK00Kxlko..lKKXXKKKKXXxoOKKKKxcckXKxlc:okkcckOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkdooo;;okOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKO:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo:xOOOOOOOOOOOOOOOOOOd:xxlOKKKKKKKKKKKKKKKKKKKKKKKKKKx. ...lkl:oO0000d:cx0Odxl. ;0KKKKKKKKxlko..,::cokkko;:lx0XKXKX0o:::c;,dOOkc:kOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkcckOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKO:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOocxOOOOOOOOOOOOOOOOOOo:xxlxOOOOOOOOOOOOOOOOOOOOOOOOOOo.  .lOOkxkOdlokOxkO0Od;.. ;kOOOOOOOOdlOo;dxxxxl;;:;,:xXXXXXXOl:lxxkxloxOkc:kOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkllkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKK0OxdddOKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxdkOOOOOOOOOOOOOOOOOOd:xxlxOOOOOOOOOOOOOOOOOOOOOOOOkOo. ..,oOOO0Oocok0Okl:;.... ,kOOOOOOOOolOo:xOOOOOOOOOkoc:::::::cxOkkkkx:ckkc;lxkkkkkOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkookOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKxcloollxk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo:xxlOKKKKKKKKKKKKKKKKKKKKKKKKKKk;..'...'cdodk00kl'.  ..'..l0KKKKKKKKxlOo:xOOOOOOOOOOOxooooooxOOOxc;:c:;:c::ll:;;;lkOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKd;lddo;ok:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo:xxlOKKKKKKKKKKKKKKKKKKKKKKKKKKK0ddOocccclc.:OO:.:cc:.ckoxKKKKKKKKKKxlko:xOOOOOOOOOOOOOOOOOOOOOOd;,;cl:,:c::c:;;,:kOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkdlcccccccccccccccccccccccclcldkOOOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKOdlclldOk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo:xxcdxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxdc;;;,..,,..,;:cclxxxxxxxxxxxxxllOo:xOOOOOOOOOOOOOOOOOOOOOOd;;:;;;,;;;;;:::,:kOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkdlcooodddddddddddddddddddddddoocldkOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKK00000KO:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo:xOdddddddddddddddddddddddddddddddddl,... ...........:oddddddddddddddxOo:kOOOOOOOOOOOOOOOOOOOOOOd;;::,,::::::::c;:kOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx::ddxkOOOOOOOOOOOOOOOOOOOOOOkxdd:;xOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOd:d0Oxdddddddddddddddddddddddddddddc.. .................,odddddddddddxOOlckOOOOOOOOOOOOOOOOOOOOOOko:;:;;::::::cl:cdkOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx::odkOOOOOOOOOOOOOOOOOOOOOOOOkdd::xOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkxoloooooooooooooooooooooooooooool;.. .............. ....,coloooooooollodkOOOOOOOOOOOOOOOOOOOOOOOOkc;oo:;llccldl;lOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx::ddxOOOOOOOOOOOOOOOOOOOOOOOOkdd::xOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx,  ................... .lOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkc;od:;ldolodl;oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx::odkOOOOOOOOOOOOOOOOOOOOOOOOkddc:xOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKO:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOd,...................... .:xkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkc:ol;,cloool:,oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx::oddkOOOOOOOOOOOOOOOOOOOOOOkxdd::xOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo. ........................;kOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkc;cc;,::clc:;,oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx:,clooxxxxxxxxxxxxxxxxxxxxxxoolc,;xOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo. ....................... ,kOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkl;::;,:::::cc:oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx:,:;;:cccccccccccccccccccccc:;;:,;xOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo. ....................... ,kOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkl;l:;llcclc:dkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkd:;:::;;;;;;;;;;;;;;;;;;;;;;:::;:okOOOOdlOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKO:oKKKxlxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOd..........................;kOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo;lc;ldoloc:xOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxc,,;:;;;;;;;;;;;;;;;;;;;;;;:;,'cxxxxxxooOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKO:oKKKkldxxxddxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxdl,,:........,'........;:';odxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxl;cc;;cooo:;oxxxxxxxkOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxdxxdc',cc,,:llllllllllllllllc;,cc,'cdxddxxddOKK0cc0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKk:oKKKkddxxddddddxddxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxdddddddxxddddlll,'..,c;..;oddddddddxxddxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxddddxl,;:;,;clc;;oxddddlcdOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOl',cc;cdxxxxxxxxxxxxxxxxdc;cl,'lkOkOkOkxdxOkclkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx;lOkxdxkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkxl:,.;cccccccc:;;lxOOOkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo,;:;,;:::,;xOOkxddcoOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo:',cc;ldddddddddoooddddodl;cl,':oooooooooooodxddxxddxxxxxxxxxxxxxxxxxxdxddxdxxxxdoloooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo;,;;..:llollolllc,';;,;loooooooooooooooooooooooooooooooooooooooooooc;::,:c:,:loooo:lkloOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKo',cc;ldddoddddddddddododl;:c,'oKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK0;'c..:llllllllllll,.c,,OKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKO::l;col;lKKKK0clOooOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKo,';,,cdxxxxxxxxxxxxxxxxdc,,;,'lKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK0;.;.'lllllllllllll,.;',OKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKO::c;:oo;lKKKK0clOooOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkl,';:;:cllllllllllllllllc:;:;,'ckkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkx;':.,ddddddddddddd:.;''dOkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkx:;;,;cc;lkkkkx:lOloOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOdoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO0koxOOOOOOOOOOOOOd;,,,;;;,,,,,,,,,,,,,,,,;;;,,,;oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO0xokOOOOOOOOOOOOOOOOOOOOOOOOOkoxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOdok0OOOOOOOOOOOOOOOOkxo';dddxo;,;'cxdx:.lxxOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOddOkl,,;;;oOOOO0k:lOooOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
KKKKkx0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKOdkKKKKKKKKKKKKKKOoloooooooolccooooooooooooolok0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKOdOKKKKKKKKKKKKKKKKKKKKKKKKK0xkKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKkx0KKKKKKKKKKKKKKKKKKKk,.;;;:;.:d'.:;:.'xKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK0xx0Ko'',,;xKKKKK0clOooOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
0000doO00000000000000000000000000000000000000000000000000000kox0000000000000000000000000Ood000000000000000000000000000000000000000000000000000000xok0000000000000000000000000kox0000000000000000000000000000000000000000000000000000000doO0000000000000000000Odl;.    :d'  .,clk00000000000000000000000000000000000000000000od00kddddxO00000OclOloOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
kkkkddkkkkkkkxxkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkxdxkkkkkkkkkkkkkkkkkkkkkkkkkkddkkkkkkxxkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkxdxkkkkkkkkkkkkkkkkkkkkkkkkkkdxkkkkkkxxkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkxdkkkkkkkkkkkkkkkkkkkkkkdc;;;;;oxc;;;;:lxkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkddkkkkkkkkkkkkkkx:lOloOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
KKKKKKKKKKKKOxOKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKOxOKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKkx0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK0xkKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK0clOloOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
OOOOOOOOOOOOkdxOOOOOOOOOOOOOOOOOOOOOOOO0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKOxOKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKkx0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK0xkKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK0clOooOOOOOOOOOOOOOOOOOOOOOOO;                                                                                        
kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkxdkxdddxkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkxoxOkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkdokkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkodkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkx:lOloOOOOOOOOOOOOOOOOOOOOOOO;                                                                                        
xO0kxk0KOxxOKKOxk0OxxxO0Oxk0K0k0OxOkddxkOOOOOdokOOOOOOOOOOOOOOOOOOOOOOOOOxokOOOOOOOOOOOOOOOOOOOOOOkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkodOOOOOOOOOOOOOOOOOOOOOOOOOOdokOOOOOOOOOOOOOOOOOOOOOOkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkoxOOOOOOOOOOOOOOOOOOOOOOOOOOkodOOOOOOOOOOOOOOOOOOOOOOkkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxoxOOOOOOOOOOOOOOOOOOOOOOOOOkodOk:lOloOOOOOOOOOOOOOOOOOOOOOOO;                                                                                        
kkkdox0KOdokK0kodkkkkkkkxodOKOxOOxOkdkkO0KKKKxx0KKKKKKKKKKKKKKKKKKKKKKKKKkdOKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK0xxKKKKKKKKKKKKKKKKKKKKKKKKKKkx0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKOdkKKKKKKKKKKKKKKKKKKKKKKKKKK0xx0KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKOdOKKKKKKKKKKKKKKKKKKKKKKKKKOdkK0clOloOOOOOOOOOOOOOOOOOOOOOOk;                                                                                        
kxdddx0Kkdok00kddoxkkKNKKKXNNNXXXXNXKKKXNNNNNXXNNNNNNNNNNNNNNNNNNNNNNNNNNXXNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNXXNNNNNNNNNNNNNNNNNNNNNNNNNNXXNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNXXNNNNNNNNNNNNNNNNNNNNNNNNNNNXXNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNXXNNNNNNNNNNNNNNNNNNNNNNNNNNNXXNN0KXKKNXXNXXXXXXXXXXXXXXXXXXXOdddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd;                   
odxxO00kddodk0KOxxdodKMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMd                   
kOKKK0kddooddk0KKK0kOXMMM0lcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccOMMMd                   
KKKKOxdodxxddddkKKKKKNMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
KK0Oxddxkkkkxdoxk0KKKNMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
0OxddxkkxdodkkxddxO0KNMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
kdddkOxdodooddOOdodx0NMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
dodkkxdddodoodxkkdodkXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
dox0kododdddddok0xoodKMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
doxOkdodddddddok0xoodKMMMo                                                                                     .:l'                                                                                                                          ,l;.                                                    .cl.                                                                                                                                 lMMMd                   
xodxkkdodddoodkOxdodONMMMo                                     'dkdccdo'                                       'x0c       .ckl.                                                                                      ,kx,    .lkc..         .l0d.                                                    ;O0:       .okc     ;kx'                                                                                                             lMMMd                   
Oxddxkkdoddodxkkdodk0NMMMo                .::.  .:,.           ;XMd..dWK;   .::::;.    .;::::.   .;:::::'   .,::cl'     ':lKMKo:'    ,::::.                .::::;.   .:;.  ':,    '::::'    ':::::;.   .:;.  ,:'  .;:kWWk:;. .OM0dl::,.   .::cl;.    .::::::.     .::::::.                        .;:ccc.     ':oXM0l:.  ,OX:         '::::'                                                                                              lMMMd                   
KOddodxkxddxkkdddk0KKNMMMo                ,dX0dxK0c.           ;XMo  oMX: .cKKl:xXO,  ;OXd:lOk;  :NMk::OXo. .,;dNMd     .:lKMKl:.  'kXk:c0Kl.             :0Ko:dX0;  lMN: .kMO' .oXO::kXx' .kMKl:dXO;  oMX: .OMk. .;;kWWx;;. .OM0dl;xXO,  .:l0MO.    cNWk:cOKo. .lK0c;xWWc                        .,;xWWl     ':oXM0c;.  'l:.       .dXO::xOl.                                                                                            lMMMd                   
K00kdoddxkkxddox0KKKKNMMMo               ,cdXMMMM0lc.          ;XMo  oMX: .kMKc;dNMo  ;0Xd;;c;.  :NWl  oMX;    ;KMd       .kMO.    cWWx;:OMK,            .xMKl;oXMd  ,OXd,c0Xl. ;XMO;;kMN: .kMO. 'dk;  oMX; .OMk.    cWNc    .OMk,. :NWl    .xMO.    cNWc  dMK, '0Mx. cWWc                           :NWl       '0Mx.    ..         'xXk;;c:.                                                                                             lMMMd                   
00KKOkddddddddxOKKK0KNMMMo                .dKkclOO:.           ;XMo  oMX: .kMKdlccc'   .:ccdK0;  :NMl  dMK;    ;XMd       .kMO.    cWWOcccc:.            .xMXdcccl'   'OWMMWl   ;XM0lcclc. .kMO.       oMX; .OMk.    cWNc    .OMk,. :NWl    .kMO.    :NWl  dMK, '0Mx. cWWc                           :NWl       '0Mx.                .,cclOKo.                                                                                            lMMMd                   
xkOO0KOxoddodOK00OkxxKMMMo                .''.  .'.            ;XMk,,xNO, .oX0:'cxo.  'oxc'c0Xl. :NMx,,kNk'  .'oNMk,'.    .lX0c'.  ,0Nd',od;.             cXKc':do'    'kMXo.   'kNk,'ldc. .kMO.       :KXo':0Mk.    ;KNd''  .OMk;. :NWl  .':0MKc'.  cWWl  dMK, .dXO;'dWWc  cOd.                   .'dWWx,'.    .dXO;'.             .cxl',xNk'                                                                                            lMMMd                   
odddxO0OxdddO00xdddodKMMMo                                     .looolo:.    ,oolol.    'loloo,   :NM0ooo;.  .:loooool'      ,ool,   .:oloo;                ,oolol'      ,oc.     .;oolo:.   ;o:.        .lolxXMk.     .colc. .:o;.  .lo'  ,loooool,  'oo'  ,oc.   ;oolOWWc  ;00,                  .cloooooc.      ,ool'              .:ollo;.                                                                                             lMMMd                   
kxddox0KOdokK0koddxkkXMMMo                                                                       :NMl                                                                                                  'od;.;ONd.                                               .:dc..oNK;  ,c,.                                                                                                                                                          lMMMd                   
dkOxox0KOdokK0xodOkddKMMMo                                                                       'dx,                                                                                                   ,oolod:.                                                 .:dolol'                                                                                                                                                                 lMMMd                   
kO0OkO00OkkO00Okk0OkkXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
kkkkxxkkkkkkkkkkkkkkOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
KKK0xd0KKKKKKKKKKKK0KNMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
kkkkdokkkkkkkkkkkkkkkXMMMo                                                    ...        cko.     .....      ......                                                                                                                                                                                                                                                                                                                       lMMMd                   
000000000000Oxk00000KNMMMo                                                   'OXo        ;dc.    .cokXK;    .:lxKXl                                                                                                                                                                                                                                                                                                                       lMMMd                   
KKKKKKKKKKKKOxOKKKKKKNMMMo                                      .colloc.   ;lxNMKoo'  'loddl.       lWN:       ;KMd                'dd'  ;dl.  .:dolol.   :dc. .ld;                                                                                                                                                                                                                                                                       lMMMd                   
OOOOOOOOOOOOkoxOOOOOOXMMMo                                     ,ONx..cdc.  ..cKMk,..  ..,kM0'       lWN:       ;XMd                cWWc .xMK, 'kNk'.dNK; .xM0' ,KMd                                                                                                                                                                                                                                                                       lMMMd                   
kkkkkkkkkkkkkxkkkkkkOXMMMo                                     .c0Ooll;.     ,KMx.      .xM0'       lWN:       ;XMd                cWWc .xM0, ,KMx. lWNc .xM0' ,KMd                                                                                                                                                                                                                                                                       lMMMd                   
KKKKKKKKKKKKKKKKKKKKKNMMMo                                      'c:,,xNk'    ,KMd       .xM0'       lWN:       ;KMo                cWWc .dM0, ,KMd  cWWc .xM0' ,KMd  .'.                                                                                                                                                                                                                                                                  lMMMd                   
KKKKKKKKKKKKKKKKKKKKKNMMMo                                     .l0kccOKo.    .lK0lc.  'cl0MXdc,  .:cOMWkc;. .;cxNM0cc.             'xXkcl0M0, .lK0lckXx'  ;OKdcxNMd  cNX:                                                                                                                                                                                                                                                                 lMMMd                   
lllllllllllllllllllll0MMMo                                       .;;;;.        .;;;.  .;;;;;;;.  .,;;;;;;,.  ';;;;;;;.              '::;:OM0,   .;;;;'     .,;;;;;.  .;,.                                                                                                                                                                                                                                                                 lMMMd                   
llllllllllllllllllllo0MMMo                                                                                                         .oOd;:OKo.                                                                                                                                                                                                                                                                                             lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                           ,::::.                                                                                                                                                                                                                                                                                               lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
OOOOOOOOOOOOOOOOOOOOOXMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
,,,,,,,,,,,,,,,,,,,,;kMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
                     oMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
                     oMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
                     oMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
                     oMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
                     oMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
                     oMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
                     oMMMo                                                                                                                                                                                                                                                                                                                                                                                                                lMMMd                   
                     oMMMk,'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''xMMMd                   
                     oMMMWNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNWMMMd                   
                     cKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKl                   
                      ........................................................................................................................................................................................................................................................................................................................................................................................................................                    
                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                                   Special thanks to:
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMWWWWWWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWNNNNNWWWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWWWWWWMMMMMMMMMMMMM
MMMMMMMMMNKkdoolloodk0XNWMMMMMMMMMMMMMMMMMMMMMMMMMWNKOxdolccccllodk0NWMMMMMMMMMMMMMMMMMMMMMMMMMWNX0kdoolloodkKNMMMMMMMMM
MMMMMMWN0dllooool:;,,,:cok0NWMMMMMMMMMMMMMMMMMMWXOdl;,,,,''''''',,,;cokXWMMMMMMMMMMMMMMMMMMWN0koc:,,,;:loooolld0NMMMMMMM
MMMMMMWX0KXNWWWWNX0xl;,',,;cdOXWMMMMMMMMMMMMMWKxc,',,,,,''''''',,,'',',:d0NMMMMMMMMMMMMMWXOdc;,,',;lx0XNWWWWNXK0XWMMMMMM
MMMMMMMMMMMMMMMMMMMMN0d;,',,,,:oOXWMMMMMMMMW0d:,,,,,,,'''''''',,,,,',,'',;oONMMMMMMMMWXOo:,',,,,:d0NMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMN0l,,,,,,,,cd0NWMMMWXx:,,,,',,,,''''''',,,,,',,''',,,;dKWMMMWN0dc,,''',,,l0WMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMXd;,,'',,,,;cx0NN0l,,,,,''''''',,,,''',,,,,'''''',,,,cONN0xc;,,,'''',;dXMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMNx;,''',,,,',;lo:,',,,''',,,,,,,,,,''''',,,,'''',,,,,;ll;,',,,,,,,,;xNMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMXo,',,',,,,''''''',,,',,,,,,,,,,,,'',''',,,,',,,,,,',,'''',,,,,,,,oXMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMWO:,,,'',,,,,,'',;:;,,,,',,'''''',,,,,''',,',,,,;:;,,,,,'',,,'''':OWMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMKl,'''',,,,','':kKKOxo:,,,,'',''''''''',,,,:oxOKKk:'',,,,,,,,'',oKMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMWk;,,,,,,,''',,oXMMMWWXkl;,,,,,,,,,,',',;lkXWMMMMKl,'','''',,,';kWMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMKl,,,,,,',,,';xNMMMMMMMN0o;,,'',,,,,,;o0NMMMMMMMKl,,''''''',',lKMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMNx;''''',,,,,;xNMMMMMMMMMNOc,,'''',,cONMMMMMMMMW0:',,,,'''''';xNMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMM0c,''',,,''',oXMMMMMMMMMMWKo;,,,';dKWMMMMMMMMWKo,,,,,',,,,',c0MMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMWk;,''',''''';xNMMMMMMMMMMMNx;'';dXNWWWNNXK0kdc,,,,,,',,,,,;kWMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNx;',,,'',,,,;dKWMMMMMMWNKko;'';clooollc::;,,'',,,,'',,,,;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNx;,',,',,,,',cdkOOOkxoc;,,,''''',,''''',,,,',,,,,'',,,;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWO:'',,',,''',,,,,,,,,,',,,''',,''''''',,,,'''''',,,,,':OWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWO:,,,,'',,',,,,,,,'',,,,,,''''''''''',,,,,'',,,'',',,';OWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWO:,,,,'',,,,,,,,,,,''',,,,'''''''''''',,,'',,,,''''',':OWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMXo,',''',,,,,,','',,,,'''''''','''',,,,'',,'','',,,,''cKMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWO:',,,,',,,,,,,,,,,,,,''',,,,,,,,,',,'',,,,,'''',,,,;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNk:,,,,,,'',',,,,,,;;,,,,,,''',,,;;,'',,,,'',',,,,,,oXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNOc,''',,',;;,,,:dOo,',,,,''',,;oOd:,,,,,,',,,,',;dXWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWKd:,',,,ck0kxkKWWk;,,,,'',,,,;kWWKxdxkxc,,''',ckNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMW0d:,,,;coxOKXNWXkxxddddddxxkXWWNXKOxl;,',,ckXWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWKko:,,',,:cloddxxkkOOOOkkxxdol:;,'',;cdOXWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWXOxo:;,,',,''',,,,,,,,,,',',,;cox0XWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWX0Okdolcc::;;;;;:::clodxOKNWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWWNXXXKKKKKXXXNWWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMNo. lNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMNo. lNMMMMMMMMMMMMMMMMWNK0O0KXWMMMMMMMMMWNXXXWMMMMMMMMMWXK000KNWMMMMMMMMMMWNK000KXWMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMM0o:'  .:l0NxcoOkokNMMWXx:'......;llcdXMMWOc,..';xXMMMMW0d;......'clclONMMMXkc'......:dKWMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMO;..   .,kX;  ...;0MNk,  'cudax'    '0MWk. .;:'  oWMMXl. .;lodo:.    cXMWO,  .cudax,. .dNMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMWNKc  cKNWX;  .:OXWNd. .xXWMMMMNx'  '0MWx. 'xKkllOWMXc  ;ONMMMMWKl.  cXWk. .oXWMMMMNk,  lNMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMNo. lNMMX;  cXMMMK; .xWMMMMMMMWk. '0MMNd'...'lONMWk. ,KMMMMMMMMNl  cXX:  oWMMMMMMMWO. 'OMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMNo. lNMMX; .xWMMM0, .kWMMMMMMMWO. '0MMMWXkc;. .lXWx. ;XMMMMMMMMWl  cXX: .dWMMMMMMMMO' .OMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMNo. lNMMX; .kWMMMNo. ;0WMMMMMW0:  '0M0l:xNWWO' .xMK; .lXMMMMMMNk'  cXWd. ,OWMMMMMWKc  :XMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMNo. lNMMX; .kWMMMMXo. .:dkOOdc.   '0Mk. .cddc. 'OMW0;  'lxOOko;.  .oNMNd. .:dkOOxc. .cKMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMNd..oNMMX: 'kWMMMMMNOl'. ..  ..,. ,0MWk;......:OWMMNO:......      'OMMMW0o,.. .  ..ckNMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMWXOOXWMMWKkONMMMMMMMMWX0kdddxOXX0k0NMMMWKOOO0XWMMMMKc...,;;;,..  'kNMMMMMWN0kxddkOXWMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWKd,.  ..  .,oKWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMN0kdoodk0NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
sgd: HylianPablo & Bolu. */

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

	// CUDA stop
	cudaDeviceSynchronize();

	/* 5. Stop global time */
	ttotal = cp_Wtime() - ttotal;

#ifdef DEBUG
	printf("List of cells at the end of the simulation: %d\n\n", num_cells );
	for( i=0; i<num_cells; i++ ) {
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
				cells[i].age );
	}
#endif // DEBUG

	/* 6. Output for leaderboard */
	printf("\n");
	/* 6.1. Total computation time */
	printf("Time: %lf\n", ttotal );

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
	free( culture );
	free( culture_cells );
	free( cells );

	/* 8. End */
	return 0;
}
