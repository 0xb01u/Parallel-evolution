/*
 * Kernel functions for the simplified simulation of life evolution
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2019/2020
 *
 * v1.6
 *
 * (c) 2020, Arturo Gonzalez Escribano
 */
#include "include/taylor_trig.h"
#include "include/glibc_nrand48.h"
#include "include/atomic.h"

#include "include/evolution.h"

/* ================ PROVIDED FUNCTIONS ================ */

/*
 * Function: Choose a new direction of movement for a cell
 * 	This function can be changed and/or optimized by the students
 */
__device__ void cell_new_direction( Cell *cell ) {
	int angle = int_urand48( INT_2PI, cell->random_seq );
	cell->mov_row = taylor_sin( angle );
	cell->mov_col = taylor_cos( angle );
}

/*
 * Function: Mutation of the movement genes on a new cell
 * 	This function can be changed and/or optimized by the students
 */
__device__ void cell_mutation( Cell *cell ) {
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

/* ================ DEVELOPED FUNCTIONS ================ */

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
extern "C" __global__ void initGPU(int *culture_d, short *culture_cells_d, int rows_d, int columns_d, Cell *cells_d1, Cell *cells_d2, int num_cells_d, Statistics *stats)
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
extern "C" __global__ void initCells(unsigned short *random_seqs_d)
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
extern "C" __global__ void step1()
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
			atomicAdd(pos, (short)1);
		}

	} // End cell movements

	// Statistics: Max age of a cell in the simulation history
	reductionMax(cells, num_cells, &sim_stat->history_max_age);
}

/*
 * Function to clean dead cells from the list of cells.
 *
 */
extern "C" __global__ void cleanCells()
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
extern "C" __global__ void placeFood(food_t *food, int num_food)
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
extern "C" __global__ void step2()
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
extern "C" __global__ void recount()
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
extern "C" __global__ void step3()
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
