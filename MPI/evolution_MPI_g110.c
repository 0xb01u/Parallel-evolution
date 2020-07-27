/*
 * Simplified simulation of life evolution
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2019/2020
 *
 * v1.2
 *
 * (c) 2020 Arturo Gonzalez Escribano
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <cputils.h>
#include <mpi.h>
#include <stddef.h>

/* Structure to store data of a cell */
typedef struct
{
    float pos_row, pos_col;       // Position
    float mov_row, mov_col;       // Direction of movement
    float choose_mov[3];          // Genes: Probabilities of 0 turning-left; 1 advance; 2 turning-right
    float storage;                // Food/Energy stored
    int age;                      // Number of steps that the cell has been alive
    unsigned short random_seq[3]; // Status value of its particular random sequence
    bool alive;                   // Flag indicating if the cell is still alive
} Cell;

/* Structure for simulation statistics */
typedef struct
{
    int history_total_cells;     // Accumulated number of cells created
    int history_dead_cells;      // Accumulated number of dead cells
    int history_max_alive_cells; // Maximum number of cells alive in a step
    int history_max_new_cells;   // Maximum number of cells created in a step
    int history_max_dead_cells;  // Maximum number of cells died in a step
    int history_max_age;         // Maximum age achieved by a cell
    float history_max_food;      // Maximum food level in a position of the culture
} Statistics;

/* 
 * Macro function to simplify accessing with two coordinates to a flattened array
 * 	This macro-function can be changed and/or optimized by the students
 *
 */
#define accessMat(arr, exp1, exp2) arr[(int)(exp1)*columns + (int)(exp2)]

/*
 * Function: Choose a new direction of movement for a cell
 * 	This function can be changed and/or optimized by the students
 */
void cell_new_direction(Cell *cell)
{
    // We stan M_TAU.
    float angle = (float)(6.283185307179586 * erand48(cell->random_seq));
    cell->mov_row = sinf(angle);
    cell->mov_col = cosf(angle);
}

/*
 * Function: Mutation of the movement genes on a new cell
 * 	This function can be changed and/or optimized by the students
 */
void cell_mutation(Cell *cell)
{
    /* 1. Select which genes change:
	 	0 Left grows taking part of the Advance part
	 	1 Advance grows taking part of the Left part
	 	2 Advance grows taking part of the Right part
	 	3 Right grows taking part of the Advance part
	*/
    int mutation_type = (int)(4 * erand48(cell->random_seq));
    /* 2. Select the amount of mutation (up to 50%) */
    float mutation_percenTAGe = (float)(0.5 * erand48(cell->random_seq));
    /* 3. Apply the mutation */
    float mutation_value;
    switch (mutation_type)
    {
    case 0:
        mutation_value = cell->choose_mov[1] * mutation_percenTAGe;
        cell->choose_mov[1] -= mutation_value;
        cell->choose_mov[0] += mutation_value;
        break;
    case 1:
        mutation_value = cell->choose_mov[0] * mutation_percenTAGe;
        cell->choose_mov[0] -= mutation_value;
        cell->choose_mov[1] += mutation_value;
        break;
    case 2:
        mutation_value = cell->choose_mov[2] * mutation_percenTAGe;
        cell->choose_mov[2] -= mutation_value;
        cell->choose_mov[1] += mutation_value;
        break;
    case 3:
        mutation_value = cell->choose_mov[1] * mutation_percenTAGe;
        cell->choose_mov[1] -= mutation_value;
        cell->choose_mov[2] += mutation_value;
        break;
    default:
        fprintf(stderr, "Error: Imposible type of mutation\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    /* 4. Correct potential precision problems */
    cell->choose_mov[2] = 1.0f - cell->choose_mov[1] - cell->choose_mov[0];
}

#ifdef DEBUG
/* 
 * Function: Print the current state of the simulation 
 */
void print_status(int iteration, int rows, int columns, float *culture, int num_cells, Cell *cells, int num_cells_alive, Statistics sim_stat)
{
    /* 
	 * You don't need to optimize this function, it is only for pretty printing and debugging purposes.
	 * It is not compiled in the production versions of the program.
	 * Thus, it is never used when measuring times in the leaderboard
	 */
    int i, j;

    printf("Iteration: %d\n", iteration);
    printf("+");
    for (j = 0; j < columns; j++)
        printf("---");
    printf("+\n");
    for (i = 0; i < rows; i++)
    {
        printf("|");
        for (j = 0; j < columns; j++)
        {
            char symbol;
            if (accessMat(culture, i, j) >= 20)
                symbol = '+';
            else if (accessMat(culture, i, j) >= 10)
                symbol = '*';
            else if (accessMat(culture, i, j) >= 5)
                symbol = '.';
            else
                symbol = ' ';

            int t;
            int counter = 0;
            for (t = 0; t < num_cells; t++)
            {
                int row = (int)(cells[t].pos_row);
                int col = (int)(cells[t].pos_col);
                if (cells[t].alive && row == i && col == j)
                {
                    counter++;
                }
            }
            if (counter > 9)
                printf("(M)");
            else if (counter > 0)
                printf("(%1d)", counter);
            else
                printf(" %c ", symbol);
        }
        printf("|\n");
    }
    printf("+");
    for (j = 0; j < columns; j++)
        printf("---");
    printf("+\n");
    printf("Num_cells_alive: %04d\nHistory( Cells: %04d, Dead: %04d, Max.alive: %04d, Max.new: %04d, Max.dead: %04d, Max.age: %04d, Max.food: %6f )\n\n",
           num_cells_alive,
           sim_stat.history_total_cells,
           sim_stat.history_dead_cells,
           sim_stat.history_max_alive_cells,
           sim_stat.history_max_new_cells,
           sim_stat.history_max_dead_cells,
           sim_stat.history_max_age,
           sim_stat.history_max_food);
}
#endif

/*
 * Function: Print usage line in stderr
 */
void show_usage(char *program_name)
{
    fprintf(stderr, "Usage: %s ", program_name);
    fprintf(stderr, "<rows> <columns> <maxIter> <max_food> <food_density> <food_level> <short_rnd1> <short_rnd2> <short_rnd3> <num_cells>\n");
    fprintf(stderr, "\tOptional arguments for special food spot: [ <row> <col> <size_rows> <size_cols> <density> <level> ]\n");
    fprintf(stderr, "\n");
}

/*
 * MAIN PROGRAM
 */
int main(int argc, char *argv[])
{
    int i, j;

    // Simulation data
    int max_iter;         // Maximum number of simulation steps
    int rows, columns;    // Cultivation area sizes
    float *culture;       // Cultivation area values
    short *culture_cells; // Ancillary structure to count the number of cells in a culture space

    float max_food;     // Maximum level of food on any position
    float food_density; // Number of food sources introduced per step
    float food_level;   // Maximum number of food level in a new source

    bool food_spot_active = false;  // Special food spot: Active
    int food_spot_row = 0;          // Special food spot: Initial row
    int food_spot_col = 0;          // Special food spot: Initial row
    int food_spot_size_rows = 0;    // Special food spot: Rows size
    int food_spot_size_cols = 0;    // Special food spot: Cols size
    float food_spot_density = 0.0f; // Special food spot: Food density
    float food_spot_level = 0.0f;   // Special food spot: Food level

    unsigned short init_random_seq[3];      // Status of the init random sequence
    unsigned short food_random_seq[3];      // Status of the food random sequence
    unsigned short food_spot_random_seq[3]; // Status of the special food spot random sequence

    int num_cells; // Number of cells currently stored in the list
    Cell *cells;   // List to store cells information

    // Statistics
    Statistics sim_stat;
    sim_stat.history_total_cells = 0;
    sim_stat.history_dead_cells = 0;
    sim_stat.history_max_alive_cells = 0;
    sim_stat.history_max_new_cells = 0;
    sim_stat.history_max_dead_cells = 0;
    sim_stat.history_max_age = 0;
    sim_stat.history_max_food = 0.0f;

    /* 0. Initialize MPI */
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* 1. Read simulation arguments */
    /* 1.1. Check minimum number of arguments */
    if (argc < 11)
    {
        fprintf(stderr, "-- Error: Not enough arguments when reading configuration from the command line\n\n");
        show_usage(argv[0]);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
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
    for (i = 0; i < 3; i++)
    {
        init_random_seq[i] = (unsigned short)atoi(argv[7 + i]);
    }

    /* 1.5. Read number of cells */
    num_cells = atoi(argv[10]);

    /* 1.6. Read special food spot */
    if (argc > 11)
    {
        if (argc < 17)
        {
            fprintf(stderr, "-- Error in number of special-food-spot arguments in the command line\n\n");
            show_usage(argv[0]);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        else
        {
            food_spot_active = true;
            food_spot_row = atoi(argv[11]);
            food_spot_col = atoi(argv[12]);
            food_spot_size_rows = atoi(argv[13]);
            food_spot_size_cols = atoi(argv[14]);
            food_spot_density = atof(argv[15]);
            food_spot_level = atof(argv[16]);

            // Check non-used trailing arguments
            if (argc > 17)
            {
                fprintf(stderr, "-- Error: too many arguments in the command line\n\n");
                show_usage(argv[0]);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
    }

#ifdef DEBUG
    /* 1.7. Print arguments */
    printf("Arguments, Rows: %d, Columns: %d, max_iter: %d\n", rows, columns, max_iter);
    printf("Arguments, Max.food: %f, Food density: %f, Food level: %f\n", max_food, food_density, food_level);
    printf("Arguments, Init Random Sequence: %hu,%hu,%hu\n", init_random_seq[0], init_random_seq[1], init_random_seq[2]);
    if (food_spot_active)
    {
        printf("Arguments, Food_spot, pos(%d,%d), size(%d,%d), Density: %f, Level: %f\n",
               food_spot_row, food_spot_col, food_spot_size_rows, food_spot_size_cols, food_spot_density, food_spot_level);
    }
    printf("Initial cells: %d\n", num_cells);
#endif // DEBUG

    /* 1.8. Initialize random sequences for food dropping */
    for (i = 0; i < 3; i++)
    {
        food_random_seq[i] = (unsigned short)nrand48(init_random_seq);
        food_spot_random_seq[i] = (unsigned short)nrand48(init_random_seq);
    }

    /* 1.9. Initialize random sequences of cells */
    cells = (Cell *)malloc(sizeof(Cell) * (size_t)num_cells);
    if (cells == NULL)
    {
        fprintf(stderr, "-- Error allocating: %d cells\n", num_cells);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    for (i = 0; i < num_cells; i++)
    {
        // Initialize the cell ramdom sequences
        for (j = 0; j < 3; j++)
            cells[i].random_seq[j] = (unsigned short)nrand48(init_random_seq);
    }

#ifdef DEBUG
    /* 1.10. Print random seed of the initial cells */
    /*
	printf("Initial cells random seeds: %d\n", num_cells );
	for( i=0; i<num_cells; i++ )
		printf("\tCell %d, Random seq: %hu,%hu,%hu\n", i, cells[i].random_seq[0], cells[i].random_seq[1], cells[i].random_seq[2] );
	*/
#endif // DEBUG

    /* 2. Start global timer */
    MPI_Barrier(MPI_COMM_WORLD);
    double ttotal = cp_Wtime();

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */

/*
 * Minimum of two numbers, as a macro function.
 * It doesn't get any simpler than that. 
 *
 */
#define min(a, b) (a < b ? a : b)
    
/*
 * Same as the last macro, but with maximum.
 * Because I feel like it.
 *
 */
#define max(a, b) (a > b ? a : b)

/* 
 * Macro function to simplify accessing with two coordinates to a flattened array.
 * 	This version takes into account matrix division by rank.
 *
 */
#define accessMatSec(arr, exp1, exp2) arr[(int)(exp1)*columns + (int)(exp2)-my_begin]

/*
 * Macro function to get the exact offset in an array for a cell.
 *
 */
#define arrayPos(cell) ((int)(cell.pos_row) * columns + (int)(cell.pos_col))
#define arrayPosCoords(expr1, expr2) ((int)(expr1) * columns + (int)(expr2))

/*
 * Macro function to check if a point in a matrix belongs to this process' section.
 *
 */
#define mine(exp1, exp2) (((int)(exp1)*columns + (int)(exp2)) >= my_begin && ((int)(exp1)*columns + (int)(exp2)) < my_begin + my_size)

/*
 * Macro functions to get the matrix section a cell belongs to.
 *
 */
#define potentialSection(cell) (arrayPos(cell) / fraction)
#define potentialSectionCoords(expr1, expr2) (arrayPosCoords(expr1, expr2) / fraction)
#define section(cell) (potentialSection(cell) - (arrayPos(cell) - potentialSection(cell) * fraction < min(remainder, potentialSection(cell))))
#define sectionCoords(expr1, expr2) (potentialSectionCoords(expr1, expr2) - (arrayPosCoords(expr1, expr2) - potentialSectionCoords(expr1, expr2) * fraction < min(remainder, potentialSectionCoords(expr1, expr2))))

/*
 * Arbitrary minimum number of matrix positions per process.
 *
 */
#define THRESHOLD min(nprocs, rows * columns)

/*
 * Arbitraty minimun number of process dedicated to the simulation.
 *
 */
#define MIN_SIMULATION_PROCESSES 1

/* 
 * Macro function measure execution times for each section, if not in a leaderboard
 * execution.
 *
 */
#ifndef CP_TABLON
#define update_time(timer)          \
	{                               \
		MPI_Barrier(simulators);	\
		timer = MPI_Wtime() - timer;\
	}
#else
#define update_time(timer)
#endif
#ifndef CP_TABLON
    double time3_1 = 0.0;
    double time3_2 = 0.0;
    double time4_1 = 0.0;
    double time4_3 = 0.0;
    double time4_X = 0.0;
    double time4_4 = 0.0;
    double time4_7 = 0.0;
    double time4_8 = 0.0;
    double time4_9 = 0.0;

    double sum_time4_1 = 0.0;
    double sum_time4_3 = 0.0;
    double sum_time4_X = 0.0;
    double sum_time4_4 = 0.0;
    double sum_time4_7 = 0.0;
    double sum_time4_8 = 0.0;
    double sum_time4_9 = 0.0;

    double max_time4_1 = 0.0;
    double max_time4_3 = 0.0;
    double max_time4_X = 0.0;
    double max_time4_4 = 0.0;
    double max_time4_7 = 0.0;
    double max_time4_8 = 0.0;
    double max_time4_9 = 0.0;
#endif

	// Non-statistics printed variables:
	int iter = 0;
	int num_cells_alive, total_cells;
	/*
	 * Choose faster CPU to execute/balance load:
	 * (In case the connection is slow.)
	 *
	 */
	bool chosen = true;
#ifdef CP_TABLON
	int size;
	char name[MPI_MAX_PROCESSOR_NAME];

	MPI_Get_processor_name(name, &size);

	if (!strcmp(name, "heracles")) chosen = false;
#endif

	// Check if no processor was chosen, and choose all processors.
	bool any_chosen;
	MPI_Allreduce(&chosen, &any_chosen, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
	if (!any_chosen)
	{
		chosen = true;
	}

	MPI_Comm all_chosen, universe, simulators;
	MPI_Comm_split(MPI_COMM_WORLD, chosen ? 1 : MPI_UNDEFINED, 0, &all_chosen);

	/*
	 * Terminate non-chosen processes:
	 *
	 */
	if (!chosen)
	{
		// Wait for all the other processes:
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Finalize();
		return 0;
	}

	// Recalculate rank:
	MPI_Comm_rank(all_chosen, &rank);

	/* 
	 * MPI constants initialization:
	 *
	 */
	int nprocs; // Number of processes available.
	MPI_Comm_size(all_chosen, &nprocs);
	#define TAG 1000
	MPI_Request request, food_generator;

	// Use one process to calculate food generation, if available:
	if (nprocs >  MIN_SIMULATION_PROCESSES)
	{
		/*
		 * Matrix division.
		 * Size for each "section" (sub-matrix in each process):
		 *
		 */
		// 3.1
		int fraction = (rows * columns)/(nprocs - 1); // Size of matrix for each process.

		// Check if there are surplus processes:
		if (fraction < THRESHOLD)
		{
			fraction = THRESHOLD;
			nprocs = (rows * columns)/THRESHOLD + 1;
		}
		int remainder = (rows * columns) % (nprocs - 1); // Remaining unassigned positions.

		/*
		 * Create definitive communicator:
		 *
		 */
		MPI_Comm_split(all_chosen, rank < nprocs ? 1 : MPI_UNDEFINED, 0, &universe);

		/*
		 * Terminate unused processes:
		 *
		 */
		if (rank >= nprocs)
		{
			// Wait for all the other processes:
			MPI_Barrier(MPI_COMM_WORLD);

			MPI_Finalize();
			return 0;
		}	

		// Re-recalculate rank:
		MPI_Comm_rank(universe, &rank);

		/*
		 * Split processes in random-calculator and simulators.
		 *
		 */
		MPI_Comm_split(universe, rank < nprocs - 1 ? 1 : MPI_UNDEFINED, rank, &simulators);

		/*
		 * Random-food generator program:
		 *
		 */
		if (rank == nprocs - 1)
		{
			bool simulating = 1;
			int num_new_sources = (int)(rows * columns * food_density);
			int num_new_sources_spot = food_spot_active ? (int)(food_spot_size_rows * food_spot_size_cols * food_spot_density) : 0;
			float **food_spots = (float **)malloc(sizeof(float *) * (size_t)(nprocs - 1));
			for (i = 0; i < (nprocs - 1); i++)
			{
				food_spots[i] = (float *)malloc(sizeof(float) * (size_t)(3 * (num_new_sources + num_new_sources_spot)));
			}
			int number_food_spots[nprocs], food_offsets[nprocs];

			MPI_Irecv(&simulating, 1, MPI_C_BOOL, 0, TAG, universe, &food_generator);
			while (simulating)
			{
				/* 4.1. Spreading new food */
				//update_time(time4_1);
				for (i = 0; i < nprocs; i++)
				{
					number_food_spots[i] = 0;
					food_offsets[i] = 0;
				}
				// Across the whole culture
				for (i = 0; i < num_new_sources; i++)
				{
					int row = (int)(rows * erand48(food_random_seq));
					int col = (int)(columns * erand48(food_random_seq));
					float food = food_level * erand48(food_random_seq);

					int section = sectionCoords(row, col);

					food_spots[section][number_food_spots[section]] = food;
					food_spots[section][number_food_spots[section] + 1] = row;
					food_spots[section][number_food_spots[section] + 2] = col;
					number_food_spots[section] += 3;
				}

				// In the special food spot
				if (food_spot_active)
				{
					for (i = 0; i < num_new_sources_spot; i++)
					{
						int row = (int)(food_spot_row + food_spot_size_rows * erand48(food_spot_random_seq));
						int col = (int)(food_spot_col + food_spot_size_cols * erand48(food_spot_random_seq));
						float food = food_spot_level * erand48(food_spot_random_seq);

						int section = sectionCoords(row, col);

						food_spots[section][number_food_spots[section]] = food;
						food_spots[section][number_food_spots[section] + 1] = row;
						food_spots[section][number_food_spots[section] + 2] = col;
						number_food_spots[section] += 3;
					}
				}

				for (i = 1; i < (nprocs - 1); i++)
				{
					food_offsets[i] = food_offsets[i - 1] + number_food_spots[i - 1];
					memcpy(&food_spots[0][food_offsets[i]], food_spots[i],  sizeof(float) * (size_t)number_food_spots[i]);
				}

				MPI_Wait(&food_generator, MPI_STATUS_IGNORE);
				if (simulating)
				{
					MPI_Scatter(number_food_spots, 1, MPI_INT, NULL, 0, MPI_INT, rank, universe);
					MPI_Scatterv(food_spots[0], number_food_spots, food_offsets, MPI_FLOAT, NULL, 0, MPI_FLOAT, rank, universe);
					MPI_Irecv(&simulating, 1, MPI_C_BOOL, 0, TAG, universe, &food_generator);
				}
				//update_time(time4_1);
			}

			for (i = 0; i < (nprocs - 1); i++)
			{
				free(food_spots[i]);
			}
			free(food_spots);

			/* 5. Stop global time */
			MPI_Barrier(MPI_COMM_WORLD);

			/* 8. End */
			MPI_Finalize();
			return 0;
		}
		/*
		 * Simulation program:
		 *
		 */
		else
		{
			int my_size = fraction + (rank < remainder);

			/*
			 * Beginning for each section:
			 *
			 */
			update_time(time3_1);
			int max_ceil = min(rank, remainder);
			int max_floor = max(rank - remainder, 0);
			int my_begin = max_ceil * (fraction + 1) + max_floor * fraction;

			/*
			 * Create datatype for Cells:
			 *
			 */
			int fields = 9;				   			// Number of field blocks.
			int array_of_blocklengths[] = {			// Number of elements per block
				1, 1, 1, 1, 3, 1, 1, 3, 1};
			MPI_Aint array_of_displacements[] = {	// Block displacements
				offsetof(Cell, pos_row),
				offsetof(Cell, pos_col),
				offsetof(Cell, mov_row),
				offsetof(Cell, mov_col),
				offsetof(Cell, choose_mov),
				offsetof(Cell, storage),
				offsetof(Cell, age),
				offsetof(Cell, random_seq),
				offsetof(Cell, alive)};
			MPI_Datatype array_of_types[] = {		// Block types
				MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT, MPI_UNSIGNED_SHORT, MPI_C_BOOL};
			MPI_Aint lb, extent;
			MPI_Datatype MPI_Cell, MPI_CellExt;
			// Create basic fields structure
			MPI_Type_create_struct(fields, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_Cell);
			// Resize to cover alignment extent
			MPI_Type_get_extent(MPI_Cell, &lb, &extent);
			MPI_Type_create_resized(MPI_Cell, lb, extent, &MPI_CellExt);
			MPI_Type_commit(&MPI_CellExt);

			/* 3. Initialize culture surface and initial cells */
			culture = (float *)calloc(sizeof(float), (size_t)my_size);
			culture_cells = (short *)calloc(sizeof(short), (size_t)my_size);
			update_time(time3_1);

			// Memory errors:
#ifndef CP_TABLON
			if (culture == NULL || culture_cells == NULL)
			{
				fprintf(stderr, "-- Error allocating culture structures for size: %d x %d \n", rows, columns);
				MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			}
#endif	

			// 3.2
			update_time(time3_2);
			total_cells = num_cells; // Total number of cells in the program.

			for (i = 0; i < total_cells; i++)
			{
				// Initial age: Between 1 and 20
				cells[i].age = 1 + (int)(19 * erand48(cells[i].random_seq));
				// Initial storage: Between 10 and 20 units
				cells[i].storage = (float)(10 + 10 * erand48(cells[i].random_seq));
				// Initial position: Anywhere in the culture arena
				cells[i].pos_row = (float)(rows * erand48(cells[i].random_seq));
				cells[i].pos_col = (float)(columns * erand48(cells[i].random_seq));

				// Calculate which section this cells belongs to:
				if (rank == section(cells[i]))
				{
					cells[i].alive = true;
					// Movement direction: Unity vector in a random direction
					cell_new_direction(&cells[i]);
					// Movement genes: Probabilities of advancing or changing direction: The sum should be 1.00
					cells[i].choose_mov[0] = 0.33f;
					cells[i].choose_mov[1] = 0.34f;
					cells[i].choose_mov[2] = 0.33f;
				}
				else
				{
					// Mark cell to be delelted:
					cells[i].alive = false;
				}
			}

			// Delete non-belonging cells (4.6):
			num_cells = 0;
			for (i = 0; i < total_cells; i++)
			{
				if (cells[i].alive)
				{
					if (num_cells != i)
					{
						cells[num_cells] = cells[i];
					}
					num_cells++;
				}
			}
			cells = (Cell *)realloc(cells, sizeof(Cell) * num_cells);
			update_time(time3_2);

			// Space for the list of new cells (maximum number of new cells in 4.4 is num_cells):
			Cell *new_cells = (Cell *)malloc(sizeof(Cell) * num_cells);
			// Number of cells moved to each process each iteration:
			int *cells_moved_to = (int *)calloc(sizeof(int), (size_t)(nprocs - 1));
			int offsets_to[nprocs - 1];
			offsets_to[0] = 0;
			// Number of cells received from each process each iteration:
			int cells_moved_from[nprocs - 1];
			int offsets_from[nprocs - 1];
			offsets_from[0] = 0;
			// Auxiliary index for each process.
			int index[nprocs - 1];

			// Memory errors:
#ifndef CP_TABLON
			if (new_cells == NULL)
			{
				fprintf(stderr, "-- Error allocating new cells structures for: %d cells\n", num_cells);
				MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			}
			if (cells_moved_to == NULL)
			{
				fprintf(stderr, "-- Error allocating send data structures for: %d processes\n", nprocs);
				MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			}
			if (cells_moved_from == NULL)
			{
				fprintf(stderr, "-- Error allocating received data structures for: %d processes\n", nprocs);
				MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			}
#endif
			// Statistics: Initialize total number of cells, and max. alive
			sim_stat.history_total_cells = total_cells;
			sim_stat.history_max_alive_cells = total_cells;

#ifdef DEBUG
			/* Show initial cells data */
			if (rank == 0)
			{
				printf("Initial cells data: %d\n", num_cells);
				for (i = 0; i < num_cells; i++)
				{
					printf("\tCell %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
						   i,
						   cells[i].pos_row,
						   cells[i].pos_col,
						   cells[i].mov_row,
						   cells[i].mov_col,
						   cells[i].choose_mov[0],
						   cells[i].choose_mov[1],
						   cells[i].choose_mov[2],
						   cells[i].storage,
						   cells[i].age);
				}
			}
#endif // DEBUG

			/* 4. Simulation */
			/* 
			 * Simulation variables initialization:
			 *
			 */
			float current_max_food = 0.0f;

			// num_cells helpers:
			int num_max_cells = num_cells;	 // For realloc-ing memory.

			// First 10 iterations are different:
			// (Some conditions are never met.)
			int first_iterations = min(10, max_iter);
			if (total_cells == 0) first_iterations = 0;
			// (There is a leaderboard test with less than 10 iterations lol.)

			bool destiny = first_iterations > 0;
			if (rank == 0)
			{
				MPI_Isend(&destiny, 1, MPI_C_BOOL, nprocs - 1, TAG, universe, &food_generator);
			}
			for (iter = 0; iter < first_iterations && current_max_food <= max_food; iter++)
			{
				/* 4.3. Cell movements */
				for (i = 0; i < (nprocs - 1); i++)
				{
					cells_moved_to[i] = 0;
				}
				// 4.3.0
				for (i = 0; i < num_cells; i++)
				{
					cells[i].age++;
					// Statistics: Max age of a cell in the simulation history
					if (cells[i].age > sim_stat.history_max_age)
						sim_stat.history_max_age = cells[i].age;

					// Consume energy to move
					cells[i].storage -= 1.0f;

					/* 4.3.2. Choose movement direction */
					float prob = (float)erand48(cells[i].random_seq);
					if (prob < cells[i].choose_mov[0])
					{
						// Turn left (90 degrees)
						float tmp = cells[i].mov_col;
						cells[i].mov_col = cells[i].mov_row;
						cells[i].mov_row = -tmp;
					}
					else if (prob >= cells[i].choose_mov[0] + cells[i].choose_mov[1])
					{
						// Turn right (90 degrees)
						float tmp = cells[i].mov_row;
						cells[i].mov_row = cells[i].mov_col;
						cells[i].mov_col = -tmp;
					}
					// else do not change the direction

					/* 4.3.3. Update position moving in the choosen direction */
					cells[i].pos_row += cells[i].mov_row;
					cells[i].pos_col += cells[i].mov_col;

					// Periodic arena: Left/Rigth edges are connected, Top/Bottom edges are connected
					if (cells[i].pos_row < 0)
						cells[i].pos_row += rows;
					if (cells[i].pos_row >= rows) // These can't be elsed.
						cells[i].pos_row -= rows;
					if (cells[i].pos_col < 0)
						cells[i].pos_col += columns;
					if (cells[i].pos_col >= columns) // These can't be elsed.
						cells[i].pos_col -= columns;

					/* 4.3.4. Annotate that there is one more cell in this culture position */
					if (mine(cells[i].pos_row, cells[i].pos_col))
					{
						accessMatSec(culture_cells, cells[i].pos_row, cells[i].pos_col) += 1;
					}
					else
					{
						int cell_section = section(cells[i]);
						cells_moved_to[cell_section]++;
					}
				} // End cell movements

				// Number of cells received from each process:
				MPI_Alltoall(cells_moved_to, 1, MPI_INT, cells_moved_from, 1, MPI_INT, simulators);

				// Fill counts and displacements for MPI_Alltoallv:
				index[0] = 0;
				for (i = 1; i < (nprocs - 1); i++)
				{
					index[i] = 0;
					offsets_to[i] = offsets_to[i - 1] + cells_moved_to[i - 1];
					offsets_from[i] = offsets_from[i - 1] + cells_moved_from[i - 1];
				}
				// Allocate memory for send and receive lists:
				int cells_received = offsets_from[i - 1] + cells_moved_from[i - 1];
				Cell *cells_to_send = (Cell *)malloc(sizeof(Cell) * (size_t)(offsets_to[i - 1] + cells_moved_to[i - 1]));
				Cell *mailbox = (Cell *)malloc(sizeof(Cell) * (size_t)cells_received);

				// Fill cells to send matrix:
				for (i = 0; i < num_cells; i++)
				{
					if (!mine(cells[i].pos_row, cells[i].pos_col))
					{
						int section = section(cells[i]);
						cells_to_send[offsets_to[section] + index[section]++] = cells[i];
						cells[i].alive = false;
					}
				}

				/* 4.6. Clean dead cells from the original list */
				// 4.6.1. Move alive cells to the left to substitute dead cells
				num_cells_alive = 0;
				for (i = 0; i < num_cells; i++)
				{
					if (cells[i].alive)
					{
						if (num_cells_alive != i)
						{
							cells[num_cells_alive] = cells[i];
						}
						num_cells_alive++;
					}
				}

				MPI_Alltoallv(cells_to_send, cells_moved_to, offsets_to, MPI_CellExt, mailbox, cells_moved_from, offsets_from, MPI_CellExt, simulators);
				free(cells_to_send);

				/* 4.7. Join cell lists: Old and new cells list */
				if (cells_received > 0)
				{
					num_cells_alive += cells_received;
					// Reallocate memory, if the list of cells is the biggest one so far:
					if (num_cells_alive > num_max_cells)
					{
						num_max_cells = num_cells_alive;
						cells = (Cell *)realloc(cells, sizeof(Cell) * num_cells_alive);
						new_cells = (Cell *)realloc(new_cells, sizeof(Cell) * num_cells_alive);
					}

					for (i = 0; i < cells_received; i++)
					{
						cells[num_cells_alive - i - 1] = mailbox[i];
						accessMatSec(culture_cells, mailbox[i].pos_row, mailbox[i].pos_col) += 1;
					}
				}
				free(mailbox);

				/* 4.1.X Receive food from food generator */
				int food_amount;
				MPI_Scatter(NULL, 0, MPI_INT, &food_amount, 1, MPI_INT, nprocs - 1, universe);

				float food[food_amount];
				MPI_Scatterv(NULL, NULL, NULL, MPI_FLOAT, food, food_amount, MPI_FLOAT, nprocs - 1, universe);

				for (i = 0; i < food_amount; i += 3)
				{
					culture[(int)food[i + 1]*columns + (int)food[i + 2] - my_begin] += food[i];
				}

				/* 4.4. Cell actions */
				for (i = 0; i < num_cells_alive; i++)
				{
					cells[i].storage += accessMatSec(culture, cells[i].pos_row, cells[i].pos_col) / accessMatSec(culture_cells, cells[i].pos_row, cells[i].pos_col);

				} // End cell actions
				num_cells = num_cells_alive;

				/* 4.8. Decrease non - harvested food */
				/* 4.5.1. Clean the food consumed by the cells in the culture data structure */
				/* 4.5. Clean ancillary data structures */
				for (i = 0; i < my_size; i++)
				{
					culture[i] *= 0.95f; // Reduce 5%
					if (culture_cells[i] > 0)
					{
						culture[i] = 0.0f;
						culture_cells[i] = 0;
					}
					else if (culture[i] > current_max_food)
					{
						current_max_food = culture[i];
					}
				}

				/* 4.9. Statistics */
				// Food reduction:
				float current_max_food_all;
				MPI_Allreduce(&current_max_food, &current_max_food_all, 1, MPI_FLOAT, MPI_MAX, simulators);
				current_max_food = current_max_food_all;
				if (rank == 0)
				{
					// Send destiny to the food generator:
					destiny = (iter + 1) < max_iter && current_max_food <= max_food;
					MPI_Wait(&food_generator, MPI_STATUS_IGNORE);
					MPI_Isend(&destiny, 1, MPI_C_BOOL, nprocs - 1, TAG, universe, &food_generator);

					// Statistics: Max food
					if (current_max_food > sim_stat.history_max_food)
					{
						sim_stat.history_max_food = current_max_food;
					}
				}

#ifdef DEBUG
				/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
				if (rank == 0)
					print_status(iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat);
#endif // DEBUG
			}

			
			for (; iter < max_iter && current_max_food <= max_food && total_cells > 0; iter++)
			{
#ifndef CP_TABLON
				sum_time4_1 += time4_1;
				sum_time4_3 += time4_3;
				sum_time4_X += time4_X;
				sum_time4_4 += time4_4;
				sum_time4_7 += time4_7;
				sum_time4_8 += time4_8;
				sum_time4_9 += time4_9;

				max_time4_1 = max(max_time4_1, time4_1);
				max_time4_3 = max(max_time4_3, time4_3);
				max_time4_X = max(max_time4_X, time4_X);
				max_time4_4 = max(max_time4_4, time4_4);
				max_time4_7 = max(max_time4_7, time4_7);
				max_time4_8 = max(max_time4_8, time4_8);
				max_time4_9 = max(max_time4_9, time4_9);

				time4_1 = 0.0;
				time4_3 = 0.0;
				time4_X = 0.0;
				time4_4 = 0.0;
				time4_7 = 0.0;
				time4_8 = 0.0;
				time4_9 = 0.0;
#endif
				/* 4.3. Cell movements */
				update_time(time4_3);
				for (i = 0; i < (nprocs - 1); i++)
				{
					cells_moved_to[i] = 0;
				}

				int step_dead_cells = 0;
				// 4.3.0
				for (i = 0; i < num_cells; i++)
				{
					cells[i].age++;
					// Statistics: Max age of a cell in the simulation history
					if (cells[i].age > sim_stat.history_max_age)
						sim_stat.history_max_age = cells[i].age;

					/* 4.3.1. Check if the cell has the needed energy to move or keep alive */
					if (cells[i].storage < 0.1f)
					{
						// Cell has died
						cells[i].alive = false;
						step_dead_cells++;
						continue;
					}
					if (cells[i].storage < 1.0f)
					{
						// Almost dying cell, it cannot move, only if enough food is dropped here it will survive
						cells[i].storage -= 0.2f;
					}
					else
					{
						// Consume energy to move
						cells[i].storage -= 1.0f;

						/* 4.3.2. Choose movement direction */
						float prob = (float)erand48(cells[i].random_seq);
						if (prob < cells[i].choose_mov[0])
						{
							// Turn left (90 degrees)
							float tmp = cells[i].mov_col;
							cells[i].mov_col = cells[i].mov_row;
							cells[i].mov_row = -tmp;
						}
						else if (prob >= cells[i].choose_mov[0] + cells[i].choose_mov[1])
						{
							// Turn right (90 degrees)
							float tmp = cells[i].mov_row;
							cells[i].mov_row = cells[i].mov_col;
							cells[i].mov_col = -tmp;
						}
						// else do not change the direction

						/* 4.3.3. Update position moving in the choosen direction */
						cells[i].pos_row += cells[i].mov_row;
						cells[i].pos_col += cells[i].mov_col;

						// Periodic arena: Left/Rigth edges are connected, Top/Bottom edges are connected
						if (cells[i].pos_row < 0)
							cells[i].pos_row += rows;
						if (cells[i].pos_row >= rows) // These can't be elsed.
							cells[i].pos_row -= rows;
						if (cells[i].pos_col < 0)
							cells[i].pos_col += columns;
						if (cells[i].pos_col >= columns) // These can't be elsed.
							cells[i].pos_col -= columns;
					}
					/* 4.3.4. Annotate that there is one more cell in this culture position */
					if (mine(cells[i].pos_row, cells[i].pos_col))
					{
						accessMatSec(culture_cells, cells[i].pos_row, cells[i].pos_col) += 1;
					}
					else
					{
						int cell_section = section(cells[i]);
						cells_moved_to[cell_section]++;
					}
				} // End cell movements
				update_time(time4_3);

				/* 4.X - Cell delivery */
				update_time(time4_X);
				// Number of cells received from each process:
				MPI_Alltoall(cells_moved_to, 1, MPI_INT, cells_moved_from, 1, MPI_INT, simulators);

				// Fill counts and displacements for MPI_Alltoallv:
				index[0] = 0;
				for (i = 1; i < (nprocs - 1); i++)
				{
					index[i] = 0;
					offsets_to[i] = offsets_to[i - 1] + cells_moved_to[i - 1];
					offsets_from[i] = offsets_from[i - 1] + cells_moved_from[i - 1];
				}
				// Allocate memory for send and receive lists:
				int cells_received = offsets_from[i - 1] + cells_moved_from[i - 1];
				Cell *cells_to_send = (Cell *)malloc(sizeof(Cell) * (size_t)(offsets_to[i - 1] + cells_moved_to[i - 1]));
				Cell *mailbox = (Cell *)malloc(sizeof(Cell) * (size_t)cells_received);

				// Fill cells to send matrix:
				for (i = 0; i < num_cells; i++)
				{
					if (!mine(cells[i].pos_row, cells[i].pos_col))
					{
						int section = section(cells[i]);
						cells_to_send[offsets_to[section] + index[section]++] = cells[i];
						cells[i].alive = false;
					}
				}

				/* 4.6. Clean dead cells from the original list */
				// 4.6.1. Move alive cells to the left to substitute dead cells
				num_cells_alive = 0;
				for (i = 0; i < num_cells; i++)
				{
					if (cells[i].alive)
					{
						if (num_cells_alive != i)
						{
							cells[num_cells_alive] = cells[i];
						}
						num_cells_alive++;
					}
				}

				MPI_Alltoallv(cells_to_send, cells_moved_to, offsets_to, MPI_CellExt, mailbox, cells_moved_from, offsets_from, MPI_CellExt, simulators);
				free(cells_to_send);

				/* 4.7. Join cell lists: Old and new cells list */
				if (cells_received > 0)
				{
					num_cells_alive += cells_received;
					// Reallocate memory, if the list of cells is the biggest one so far:
					if (num_cells_alive > num_max_cells)
					{
						num_max_cells = num_cells_alive;
						cells = (Cell *)realloc(cells, sizeof(Cell) * num_cells_alive);
						new_cells = (Cell *)realloc(new_cells, sizeof(Cell) * num_cells_alive);
					}

					for (i = 0; i < cells_received; i++)
					{
						cells[num_cells_alive - i - 1] = mailbox[i];
						accessMatSec(culture_cells, mailbox[i].pos_row, mailbox[i].pos_col) += 1;
					}
				}
				free(mailbox);
				update_time(time4_X);

				/* 4.1.X Receive food from food generator */
				update_time(time4_1);
				int food_amount;
				MPI_Scatter(NULL, 0, MPI_INT, &food_amount, 1, MPI_INT, nprocs - 1, universe);

				float food[food_amount];
				MPI_Scatterv(NULL, NULL, NULL, MPI_FLOAT, food, food_amount, MPI_FLOAT, nprocs - 1, universe);

				for (i = 0; i < food_amount; i += 3)
				{
					culture[(int)food[i + 1]*columns + (int)food[i + 2] - my_begin] += food[i];
				}
				update_time(time4_1);

				/* 4.4. Cell actions */
				update_time(time4_4);
				int step_new_cells = 0;
				for (i = 0; i < num_cells_alive; i++)
				{
					/* 4.4.1. Food harvesting */
					float food = accessMatSec(culture, cells[i].pos_row, cells[i].pos_col);

					short count = accessMatSec(culture_cells, cells[i].pos_row, cells[i].pos_col);

					float my_food = food / count;
					cells[i].storage += my_food;

					/* 4.4.2. Split cell if the conditions are met: Enough maturity and energy */
					if (cells[i].age > 30 && cells[i].storage > 20)
					{
						// Split: Create new cell
						step_new_cells++;

						// Split energy stored and update age in both cells
						cells[i].storage /= 2.0f;
						cells[i].age = 1;

						// New cell is a copy of parent cell
						new_cells[step_new_cells - 1] = cells[i];

						// Random seed for the new cell, obtained using the parent random sequence
						new_cells[step_new_cells - 1].random_seq[0] = (unsigned short)nrand48(cells[i].random_seq);
						new_cells[step_new_cells - 1].random_seq[1] = (unsigned short)nrand48(cells[i].random_seq);
						new_cells[step_new_cells - 1].random_seq[2] = (unsigned short)nrand48(cells[i].random_seq);

						// Both cells start in random directions
						cell_new_direction(&cells[i]);
						cell_new_direction(&new_cells[step_new_cells - 1]);

						// Mutations of the movement genes in both cells
						cell_mutation(&cells[i]);
						cell_mutation(&new_cells[step_new_cells - 1]);
					}
				} // End cell actions
				num_cells = num_cells_alive;
				num_cells_alive += step_new_cells;

				// Cells reductions:
				int sum_stats[2] = { step_new_cells, step_dead_cells };
				int global_sum_stats[2];

				MPI_Iallreduce(sum_stats, global_sum_stats, 2, MPI_INT, MPI_SUM, simulators, &request);
				update_time(time4_4);

				/* 4.7. Join cell lists: Old and new cells list */
				update_time(time4_7);
				if (num_cells_alive > num_max_cells)
				{
					num_max_cells = num_cells_alive;
					cells = (Cell *)realloc(cells, sizeof(Cell) * num_cells_alive);
					new_cells = (Cell *)realloc(new_cells, sizeof(Cell) * num_cells_alive);
				}
				if (step_new_cells > 0)
				{
					for (i = 0; i < step_new_cells; i++)
						cells[num_cells + i] = new_cells[i];
					num_cells += step_new_cells;
				}
				update_time(time4_7);

				/* 4.8. Decrease non - harvested food */
				/* 4.5.1. Clean the food consumed by the cells in the culture data structure */
				/* 4.5. Clean ancillary data structures */
				update_time(time4_8);
				for (i = 0; i < my_size; i++)
				{
					culture[i] *= 0.95f; // Reduce 5%
					if (culture_cells[i] > 0)
					{
						culture[i] = 0.0f;
						culture_cells[i] = 0;
					}
					else if (culture[i] > current_max_food)
					{
						current_max_food = culture[i];
					}
				}
				update_time(time4_8);

				/* 4.9. Statistics */
				update_time(time4_9);
				// Food reduction:
				float current_max_food_all;
				MPI_Allreduce(&current_max_food, &current_max_food_all, 1, MPI_FLOAT, MPI_MAX, simulators);
				current_max_food = current_max_food_all;
				MPI_Wait(&request, MPI_STATUS_IGNORE);
				total_cells += (global_sum_stats[0] - global_sum_stats[1]);
				if (rank == 0)
				{
					// Send destiny to the food generator:
					destiny = (iter + 1) < max_iter && current_max_food <= max_food && total_cells > 0;
					MPI_Wait(&food_generator, MPI_STATUS_IGNORE);
					MPI_Isend(&destiny, 1, MPI_C_BOOL, nprocs - 1, TAG, universe, &food_generator);

					// Statistics: Max food
					if (current_max_food > sim_stat.history_max_food)
						sim_stat.history_max_food = current_max_food;
					// Statistics: Max new cells per step
					sim_stat.history_total_cells += global_sum_stats[0];
					if (global_sum_stats[0] > sim_stat.history_max_new_cells)
						sim_stat.history_max_new_cells = global_sum_stats[0];
					// Statistics: Accumulated dead and Max dead cells per step
					sim_stat.history_dead_cells += global_sum_stats[1];
					if (global_sum_stats[1] > sim_stat.history_max_dead_cells)
						sim_stat.history_max_dead_cells = global_sum_stats[1];
					// Statistics: Max alive cells per step
					if (total_cells > sim_stat.history_max_alive_cells)
						sim_stat.history_max_alive_cells = total_cells;
				}
				update_time(time4_9);

#ifdef DEBUG
				/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
				if (rank == 0)
					print_status(iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat);
#endif // DEBUG
			}

			// Last reduction (age):
			int max_age_root;
			MPI_Reduce(&sim_stat.history_max_age, &max_age_root, 1, MPI_INT, MPI_MAX, 0, simulators);
			sim_stat.history_max_age = max_age_root;

			num_cells_alive = total_cells;

			// Let's not be bad...
			free(new_cells);
			free(cells_moved_to);
		}
	}
	// Not enough processes execute the program:
	else
	{
		/*
		 * Matrix division.
		 * Size for each "section" (sub-matrix in each process):
		 *
		 */
		// 3.1
		int fraction = (rows * columns)/nprocs; // Size of matrix for each process.

		// Check if there are surplus processes:
		if (fraction < THRESHOLD)
		{
			fraction = THRESHOLD;
			nprocs = (rows * columns)/THRESHOLD;
		}

		/*
		 * Create definitive communicator:
		 *
		 */
		MPI_Comm universe;
		MPI_Comm_split(all_chosen, rank < nprocs ? 1 : MPI_UNDEFINED, 0, &universe);

		/*
		 * Terminate unused processes:
		 *
		 */
		if (rank >= nprocs)
		{
			// Wait for all the other processes:
			MPI_Barrier(MPI_COMM_WORLD);

			MPI_Finalize();
			return 0;
		}

		// Re-recalculate rank:
		MPI_Comm_rank(universe, &rank);	

		int remainder = (rows * columns) % nprocs; // Remaining unasigned positions.
		int my_size = fraction + (rank < remainder);

		/*
		 * Beginning for each section:
		 *
		 */
		int max_ceil = min(rank, remainder);
		int max_floor = max(rank - remainder, 0);
		int my_begin = max_ceil * (fraction + 1) + max_floor * fraction;

		total_cells = num_cells; // Total number of cells in the program.

		/*
		 * Create datatype for Cells:
		 *
		 */
		int fields = 9;				   			// Number of field blocks.
		int array_of_blocklengths[] = {			// Number of elements per block
			1, 1, 1, 1, 3, 1, 1, 3, 1};
		MPI_Aint array_of_displacements[] = {	// Block displacements
			offsetof(Cell, pos_row),
			offsetof(Cell, pos_col),
			offsetof(Cell, mov_row),
			offsetof(Cell, mov_col),
			offsetof(Cell, choose_mov),
			offsetof(Cell, storage),
			offsetof(Cell, age),
			offsetof(Cell, random_seq),
			offsetof(Cell, alive)};
		MPI_Datatype array_of_types[] = {		// Block types
			MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT, MPI_UNSIGNED_SHORT, MPI_C_BOOL};
		MPI_Aint lb, extent;
		MPI_Datatype MPI_Cell, MPI_CellExt;
		// Create basic fields structure
		MPI_Type_create_struct(fields, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_Cell);
		// Resize to cover alignment extent
		MPI_Type_get_extent(MPI_Cell, &lb, &extent);
		MPI_Type_create_resized(MPI_Cell, lb, extent, &MPI_CellExt);
		MPI_Type_commit(&MPI_CellExt);

		/* 3. Initialize culture surface and initial cells */
		culture = (float *)calloc(sizeof(float), (size_t)my_size);
		culture_cells = (short *)calloc(sizeof(short), (size_t)my_size);

		// Memory errors:
#ifndef CP_TABLON
		if (culture == NULL || culture_cells == NULL)
		{
			fprintf(stderr, "-- Error allocating culture structures for size: %d x %d \n", rows, columns);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
#endif	

		// 3.2
		for (i = 0; i < total_cells; i++)
		{
			// Initial age: Between 1 and 20
			cells[i].age = 1 + (int)(19 * erand48(cells[i].random_seq));
			// Initial storage: Between 10 and 20 units
			cells[i].storage = (float)(10 + 10 * erand48(cells[i].random_seq));
			// Initial position: Anywhere in the culture arena
			cells[i].pos_row = (float)(rows * erand48(cells[i].random_seq));
			cells[i].pos_col = (float)(columns * erand48(cells[i].random_seq));

			// Calculate which section this cells belongs to:
			if (rank == section(cells[i]))
			{
				cells[i].alive = true;
				// Movement direction: Unity vector in a random direction
				cell_new_direction(&cells[i]);
				// Movement genes: Probabilities of advancing or changing direction: The sum should be 1.00
				cells[i].choose_mov[0] = 0.33f;
				cells[i].choose_mov[1] = 0.34f;
				cells[i].choose_mov[2] = 0.33f;
			}
			else
			{
				// Mark cell to be delelted:
				cells[i].alive = false;
			}
		}

		// Delete non-belonging cells (4.6):
		num_cells = 0;
		for (i = 0; i < total_cells; i++)
		{
			if (cells[i].alive)
			{
				if (num_cells != i)
				{
					cells[num_cells] = cells[i];
				}
				num_cells++;
			}
		}
		cells = (Cell *)realloc(cells, sizeof(Cell) * num_cells);

		// Space for the list of new cells (maximum number of new cells in 4.4 is num_cells):
		Cell *new_cells = (Cell *)malloc(sizeof(Cell) * num_cells);
		// Number of cells moved to each process each iteration:
		int *cells_moved_to = (int *)calloc(sizeof(int), (size_t)nprocs);
		int offsets_to[nprocs];
		offsets_to[0] = 0;
		// Number of cells received from each process each iteration:
		int cells_moved_from[nprocs];
		int offsets_from[nprocs];
		offsets_from[0] = 0;
		// Auxiliary index for each process.
		int index[nprocs];

		// Memory errors:
#ifndef CP_TABLON
		if (new_cells == NULL)
		{
			fprintf(stderr, "-- Error allocating new cells structures for: %d cells\n", num_cells);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
		if (cells_moved_to == NULL)
		{
			fprintf(stderr, "-- Error allocating send data structures for: %d processes\n", nprocs);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
		if (cells_moved_from == NULL)
		{
			fprintf(stderr, "-- Error allocating received data structures for: %d processes\n", nprocs);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
#endif
		// Statistics: Initialize total number of cells, and max. alive
		sim_stat.history_total_cells = total_cells;
		sim_stat.history_max_alive_cells = total_cells;

#ifdef DEBUG
		/* Show initial cells data */
		if (rank == 0)
		{
			printf("Initial cells data: %d\n", num_cells);
			for (i = 0; i < num_cells; i++)
			{
				printf("\tCell %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
					   i,
					   cells[i].pos_row,
					   cells[i].pos_col,
					   cells[i].mov_row,
					   cells[i].mov_col,
					   cells[i].choose_mov[0],
					   cells[i].choose_mov[1],
					   cells[i].choose_mov[2],
					   cells[i].storage,
					   cells[i].age);
			}
		}
#endif // DEBUG

		/* 4. Simulation */
		/* 
		 * Simulation variables initialization:
		 *
		 */
		float current_max_food = 0.0f;

		// For 4.1:
		int num_new_sources = (int)(rows * columns * food_density);
		int num_new_sources_spot = food_spot_active ? (int)(food_spot_size_rows * food_spot_size_cols * food_spot_density) : 0;
		int max_sources = max(num_new_sources, num_new_sources_spot);
		float rand4_1[3 * max_sources];

		// num_cells helpers:
		int num_cells_alive;
		int num_max_cells = num_cells;	 // For realloc-ing memory.

		// First 10 iterations are different:
		// (Some conditions are never met.)
		int first_iterations = min(10, max_iter);
		// (There is a leaderboard test with less than 10 iterations lol.)
		for (iter = 0; iter < first_iterations && current_max_food <= max_food; iter++)
		{
			/* 4.1. Spreading new food */
			// Across the whole culture
			j = 0;
			for (i = 0; i < num_new_sources; i++)
			{
				rand4_1[3 * j] = (int)(rows * erand48(food_random_seq));
				rand4_1[3 * j + 1] = (int)(columns * erand48(food_random_seq));
				rand4_1[3 * j + 2] = food_level * erand48(food_random_seq);
				if (mine(rand4_1[3 * j], rand4_1[3 * j + 1])) j++;
			}
			for (i = 0; i < j; i++)
				accessMatSec(culture, rand4_1[3 * i], rand4_1[3 * i + 1]) += rand4_1[3 * i + 2];

			// In the special food spot
			if (food_spot_active)
			{
				j = 0;
				for (i = 0; i < num_new_sources_spot; i++)
				{
					rand4_1[3 * j] = (int)(food_spot_row + food_spot_size_rows * erand48(food_spot_random_seq));
					rand4_1[3 * j + 1] = (int)(food_spot_col + food_spot_size_cols * erand48(food_spot_random_seq));
					rand4_1[3 * j + 2] = food_spot_level * erand48(food_spot_random_seq);
					if (mine(rand4_1[3 * j], rand4_1[3 * j + 1])) j++;
				}
				for (i = 0; i < j; i++)
					accessMatSec(culture, rand4_1[3 * i], rand4_1[3 * i + 1]) += rand4_1[3 * i + 2];
			}

			/* 4.3. Cell movements */
			for (i = 0; i < nprocs; i++)
			{
				cells_moved_to[i] = 0;
			}
			// 4.3.0
			for (i = 0; i < num_cells; i++)
			{
				cells[i].age++;
				// Statistics: Max age of a cell in the simulation history
				if (cells[i].age > sim_stat.history_max_age)
					sim_stat.history_max_age = cells[i].age;

				// Consume energy to move
				cells[i].storage -= 1.0f;

				/* 4.3.2. Choose movement direction */
				float prob = (float)erand48(cells[i].random_seq);
				if (prob < cells[i].choose_mov[0])
				{
					// Turn left (90 degrees)
					float tmp = cells[i].mov_col;
					cells[i].mov_col = cells[i].mov_row;
					cells[i].mov_row = -tmp;
				}
				else if (prob >= cells[i].choose_mov[0] + cells[i].choose_mov[1])
				{
					// Turn right (90 degrees)
					float tmp = cells[i].mov_row;
					cells[i].mov_row = cells[i].mov_col;
					cells[i].mov_col = -tmp;
				}
				// else do not change the direction

				/* 4.3.3. Update position moving in the choosen direction */
				cells[i].pos_row += cells[i].mov_row;
				cells[i].pos_col += cells[i].mov_col;

				// Periodic arena: Left/Rigth edges are connected, Top/Bottom edges are connected
				if (cells[i].pos_row < 0)
					cells[i].pos_row += rows;
				if (cells[i].pos_row >= rows) // These can't be elsed.
					cells[i].pos_row -= rows;
				if (cells[i].pos_col < 0)
					cells[i].pos_col += columns;
				if (cells[i].pos_col >= columns) // These can't be elsed.
					cells[i].pos_col -= columns;

				/* 4.3.4. Annotate that there is one more cell in this culture position */
				if (mine(cells[i].pos_row, cells[i].pos_col))
				{
					accessMatSec(culture_cells, cells[i].pos_row, cells[i].pos_col) += 1;
				}
				else
				{
					int cell_section = section(cells[i]);
					cells_moved_to[cell_section]++;
				}
			} // End cell movements

			// Number of cells received from each process:
			MPI_Alltoall(cells_moved_to, 1, MPI_INT, cells_moved_from, 1, MPI_INT, universe);

			// Fill counts and displacements for MPI_Alltoallv:
			index[0] = 0;
			for (i = 1; i < nprocs; i++)
			{
				index[i] = 0;
				offsets_to[i] = offsets_to[i - 1] + cells_moved_to[i - 1];
				offsets_from[i] = offsets_from[i - 1] + cells_moved_from[i - 1];
			}
			// Allocate memory for send and receive lists:
			int cells_received = offsets_from[i - 1] + cells_moved_from[i - 1];
			Cell *cells_to_send = (Cell *)malloc(sizeof(Cell) * (size_t)(offsets_to[i - 1] + cells_moved_to[i - 1]));
			Cell *mailbox = (Cell *)malloc(sizeof(Cell) * (size_t)cells_received);

			// Fill cells to send matrix:
			for (i = 0; i < num_cells; i++)
			{
				if (!mine(cells[i].pos_row, cells[i].pos_col))
				{
					int section = section(cells[i]);
					cells_to_send[offsets_to[section] + index[section]++] = cells[i];
					cells[i].alive = false;
				}
			}

			/* 4.6. Clean dead cells from the original list */
			// 4.6.1. Move alive cells to the left to substitute dead cells
			num_cells_alive = 0;
			for (i = 0; i < num_cells; i++)
			{
				if (cells[i].alive)
				{
					if (num_cells_alive != i)
					{
						cells[num_cells_alive] = cells[i];
					}
					num_cells_alive++;
				}
			}

			MPI_Alltoallv(cells_to_send, cells_moved_to, offsets_to, MPI_CellExt, mailbox, cells_moved_from, offsets_from, MPI_CellExt, universe);
			free(cells_to_send);

			/* 4.7. Join cell lists: Old and new cells list */
			if (cells_received > 0)
			{
				num_cells_alive += cells_received;
				// Reallocate memory, if the list of cells is the biggest one so far:
				if (num_cells_alive > num_max_cells)
				{
					num_max_cells = num_cells_alive;
					cells = (Cell *)realloc(cells, sizeof(Cell) * num_cells_alive);
					new_cells = (Cell *)realloc(new_cells, sizeof(Cell) * num_cells_alive);
				}

				for (i = 0; i < cells_received; i++)
				{
					cells[num_cells_alive - i - 1] = mailbox[i];
					accessMatSec(culture_cells, mailbox[i].pos_row, mailbox[i].pos_col) += 1;
				}
			}
			free(mailbox);

			/* 4.4. Cell actions */
			for (i = 0; i < num_cells_alive; i++)
			{
				cells[i].storage += accessMatSec(culture, cells[i].pos_row, cells[i].pos_col) / accessMatSec(culture_cells, cells[i].pos_row, cells[i].pos_col);

			} // End cell actions
			num_cells = num_cells_alive;

			/* 4.8. Decrease non - harvested food */
			/* 4.5.1. Clean the food consumed by the cells in the culture data structure */
			/* 4.5. Clean ancillary data structures */
			for (i = 0; i < my_size; i++)
			{
				culture[i] *= 0.95f; // Reduce 5%
				if (culture_cells[i] > 0)
				{
					culture[i] = 0.0f;
					culture_cells[i] = 0;
				}
				else if (culture[i] > current_max_food)
				{
					current_max_food = culture[i];
				}
			}

			/* 4.9. Statistics */
			// Food reduction:
			float current_max_food_all;
			MPI_Allreduce(&current_max_food, &current_max_food_all, 1, MPI_FLOAT, MPI_MAX, universe);
			current_max_food = current_max_food_all;
			if (rank == 0)
			{
				// Statistics: Max food
				if (current_max_food > sim_stat.history_max_food)
					sim_stat.history_max_food = current_max_food;
			}

#ifdef DEBUG
			/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
			if (rank == 0)
				print_status(iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat);
#endif // DEBUG
		}

		
		for (; iter < max_iter && current_max_food <= max_food && total_cells > 0; iter++)
		{
#ifndef CP_TABLON
			sum_time4_1 += time4_1;
			sum_time4_3 += time4_3;
			sum_time4_X += time4_X;
			sum_time4_4 += time4_4;
			sum_time4_7 += time4_7;
			sum_time4_8 += time4_8;
			sum_time4_9 += time4_9;

			max_time4_1 = max(max_time4_1, time4_1);
			max_time4_3 = max(max_time4_3, time4_3);
			max_time4_X = max(max_time4_X, time4_X);
			max_time4_4 = max(max_time4_4, time4_4);
			max_time4_7 = max(max_time4_7, time4_7);
			max_time4_8 = max(max_time4_8, time4_8);
			max_time4_9 = max(max_time4_9, time4_9);

			time4_1 = 0.0;
			time4_3 = 0.0;
			time4_X = 0.0;
			time4_4 = 0.0;
			time4_7 = 0.0;
			time4_8 = 0.0;
			time4_9 = 0.0;
#endif
			/* 4.1. Spreading new food */
			// Across the whole culture
			j = 0;
			for (i = 0; i < num_new_sources; i++)
			{
				rand4_1[3 * j] = (int)(rows * erand48(food_random_seq));
				rand4_1[3 * j + 1] = (int)(columns * erand48(food_random_seq));
				rand4_1[3 * j + 2] = food_level * erand48(food_random_seq);
				if (mine(rand4_1[3 * j], rand4_1[3 * j + 1])) j++;
			}
			for (i = 0; i < j; i++)
				accessMatSec(culture, rand4_1[3 * i], rand4_1[3 * i + 1]) += rand4_1[3 * i + 2];

			// In the special food spot
			if (food_spot_active)
			{
				j = 0;
				for (i = 0; i < num_new_sources_spot; i++)
				{
					rand4_1[3 * j] = (int)(food_spot_row + food_spot_size_rows * erand48(food_spot_random_seq));
					rand4_1[3 * j + 1] = (int)(food_spot_col + food_spot_size_cols * erand48(food_spot_random_seq));
					rand4_1[3 * j + 2] = food_spot_level * erand48(food_spot_random_seq);
					if (mine(rand4_1[3 * j], rand4_1[3 * j + 1])) j++;
				}
				for (i = 0; i < j; i++)
					accessMatSec(culture, rand4_1[3 * i], rand4_1[3 * i + 1]) += rand4_1[3 * i + 2];
			}

			/* 4.3. Cell movements */
			for (i = 0; i < nprocs; i++)
			{
				cells_moved_to[i] = 0;
			}

			int step_dead_cells = 0;
			// 4.3.0
			for (i = 0; i < num_cells; i++)
			{
				cells[i].age++;
				// Statistics: Max age of a cell in the simulation history
				if (cells[i].age > sim_stat.history_max_age)
					sim_stat.history_max_age = cells[i].age;

				/* 4.3.1. Check if the cell has the needed energy to move or keep alive */
				if (cells[i].storage < 0.1f)
				{
					// Cell has died
					cells[i].alive = false;
					step_dead_cells++;
					continue;
				}
				if (cells[i].storage < 1.0f)
				{
					// Almost dying cell, it cannot move, only if enough food is dropped here it will survive
					cells[i].storage -= 0.2f;
				}
				else
				{
					// Consume energy to move
					cells[i].storage -= 1.0f;

					/* 4.3.2. Choose movement direction */
					float prob = (float)erand48(cells[i].random_seq);
					if (prob < cells[i].choose_mov[0])
					{
						// Turn left (90 degrees)
						float tmp = cells[i].mov_col;
						cells[i].mov_col = cells[i].mov_row;
						cells[i].mov_row = -tmp;
					}
					else if (prob >= cells[i].choose_mov[0] + cells[i].choose_mov[1])
					{
						// Turn right (90 degrees)
						float tmp = cells[i].mov_row;
						cells[i].mov_row = cells[i].mov_col;
						cells[i].mov_col = -tmp;
					}
					// else do not change the direction

					/* 4.3.3. Update position moving in the choosen direction */
					cells[i].pos_row += cells[i].mov_row;
					cells[i].pos_col += cells[i].mov_col;

					// Periodic arena: Left/Rigth edges are connected, Top/Bottom edges are connected
					if (cells[i].pos_row < 0)
						cells[i].pos_row += rows;
					if (cells[i].pos_row >= rows) // These can't be elsed.
						cells[i].pos_row -= rows;
					if (cells[i].pos_col < 0)
						cells[i].pos_col += columns;
					if (cells[i].pos_col >= columns) // These can't be elsed.
						cells[i].pos_col -= columns;
				}
				/* 4.3.4. Annotate that there is one more cell in this culture position */
				if (mine(cells[i].pos_row, cells[i].pos_col))
				{
					accessMatSec(culture_cells, cells[i].pos_row, cells[i].pos_col) += 1;
				}
				else
				{
					int cell_section = section(cells[i]);
					cells_moved_to[cell_section]++;
				}
			} // End cell movements

			/* 4.X - Cell delivery */
			// Number of cells received from each process:
			MPI_Alltoall(cells_moved_to, 1, MPI_INT, cells_moved_from, 1, MPI_INT, universe);

			// Fill counts and displacements for MPI_Alltoallv:
			index[0] = 0;
			for (i = 1; i < nprocs; i++)
			{
				index[i] = 0;
				offsets_to[i] = offsets_to[i - 1] + cells_moved_to[i - 1];
				offsets_from[i] = offsets_from[i - 1] + cells_moved_from[i - 1];
			}
			// Allocate memory for send and receive lists:
			int cells_received = offsets_from[i - 1] + cells_moved_from[i - 1];
			Cell *cells_to_send = (Cell *)malloc(sizeof(Cell) * (size_t)(offsets_to[i - 1] + cells_moved_to[i - 1]));
			Cell *mailbox = (Cell *)malloc(sizeof(Cell) * (size_t)cells_received);

			// Fill cells to send matrix:
			for (i = 0; i < num_cells; i++)
			{
				if (!mine(cells[i].pos_row, cells[i].pos_col))
				{
					int section = section(cells[i]);
					cells_to_send[offsets_to[section] + index[section]++] = cells[i];
					cells[i].alive = false;
				}
			}

			/* 4.6. Clean dead cells from the original list */
			// 4.6.1. Move alive cells to the left to substitute dead cells
			num_cells_alive = 0;
			for (i = 0; i < num_cells; i++)
			{
				if (cells[i].alive)
				{
					if (num_cells_alive != i)
					{
						cells[num_cells_alive] = cells[i];
					}
					num_cells_alive++;
				}
			}


			MPI_Alltoallv(cells_to_send, cells_moved_to, offsets_to, MPI_CellExt, mailbox, cells_moved_from, offsets_from, MPI_CellExt, universe);
			free(cells_to_send);

			/* 4.7. Join cell lists: Old and new cells list */
			if (cells_received > 0)
			{
				num_cells_alive += cells_received;
				// Reallocate memory, if the list of cells is the biggest one so far:
				if (num_cells_alive > num_max_cells)
				{
					num_max_cells = num_cells_alive;
					cells = (Cell *)realloc(cells, sizeof(Cell) * num_cells_alive);
					new_cells = (Cell *)realloc(new_cells, sizeof(Cell) * num_cells_alive);
				}

				for (i = 0; i < cells_received; i++)
				{
					cells[num_cells_alive - i - 1] = mailbox[i];
					accessMatSec(culture_cells, mailbox[i].pos_row, mailbox[i].pos_col) += 1;
				}
			}
			free(mailbox);

			/* 4.4. Cell actions */
			int step_new_cells = 0;
			for (i = 0; i < num_cells_alive; i++)
			{
				/* 4.4.1. Food harvesting */
				float food = accessMatSec(culture, cells[i].pos_row, cells[i].pos_col);

				short count = accessMatSec(culture_cells, cells[i].pos_row, cells[i].pos_col);

				float my_food = food / count;
				cells[i].storage += my_food;

				/* 4.4.2. Split cell if the conditions are met: Enough maturity and energy */
				if (cells[i].age > 30 && cells[i].storage > 20)
				{
					// Split: Create new cell
					step_new_cells++;

					// Split energy stored and update age in both cells
					cells[i].storage /= 2.0f;
					cells[i].age = 1;

					// New cell is a copy of parent cell
					new_cells[step_new_cells - 1] = cells[i];

					// Random seed for the new cell, obtained using the parent random sequence
					new_cells[step_new_cells - 1].random_seq[0] = (unsigned short)nrand48(cells[i].random_seq);
					new_cells[step_new_cells - 1].random_seq[1] = (unsigned short)nrand48(cells[i].random_seq);
					new_cells[step_new_cells - 1].random_seq[2] = (unsigned short)nrand48(cells[i].random_seq);

					// Both cells start in random directions
					cell_new_direction(&cells[i]);
					cell_new_direction(&new_cells[step_new_cells - 1]);

					// Mutations of the movement genes in both cells
					cell_mutation(&cells[i]);
					cell_mutation(&new_cells[step_new_cells - 1]);
				}
			} // End cell actions
			num_cells = num_cells_alive;
			num_cells_alive += step_new_cells;

			// Cells reductions:
			int sum_stats[2] = { step_new_cells, step_dead_cells };
			int global_sum_stats[2];

			MPI_Iallreduce(sum_stats, global_sum_stats, 2, MPI_INT, MPI_SUM, universe, &request);

			/* 4.7. Join cell lists: Old and new cells list */
			if (num_cells_alive > num_max_cells)
			{
				num_max_cells = num_cells_alive;
				cells = (Cell *)realloc(cells, sizeof(Cell) * num_cells_alive);
				new_cells = (Cell *)realloc(new_cells, sizeof(Cell) * num_cells_alive);
			}
			if (step_new_cells > 0)
			{
				for (i = 0; i < step_new_cells; i++)
					cells[num_cells + i] = new_cells[i];
				num_cells += step_new_cells;
			}
			num_cells = num_cells_alive;

			/* 4.8. Decrease non - harvested food */
			/* 4.5.1. Clean the food consumed by the cells in the culture data structure */
			/* 4.5. Clean ancillary data structures */
			for (i = 0; i < my_size; i++)
			{
				culture[i] *= 0.95f; // Reduce 5%
				if (culture_cells[i] > 0)
				{
					culture[i] = 0.0f;
					culture_cells[i] = 0;
				}
				else if (culture[i] > current_max_food)
				{
					current_max_food = culture[i];
				}
			}

			/* 4.9. Statistics */
			// Food reduction:
			float current_max_food_all;
			MPI_Allreduce(&current_max_food, &current_max_food_all, 1, MPI_FLOAT, MPI_MAX, universe);
			current_max_food = current_max_food_all;

			// Wait for the other statistics.
			MPI_Wait(&request, MPI_STATUS_IGNORE);
			total_cells += (global_sum_stats[0] - global_sum_stats[1]);
			if (rank == 0)
			{
				// Statistics: Max food
				if (current_max_food > sim_stat.history_max_food)
					sim_stat.history_max_food = current_max_food;
				// Statistics: Max new cells per step
				sim_stat.history_total_cells += global_sum_stats[0];
				if (global_sum_stats[0] > sim_stat.history_max_new_cells)
					sim_stat.history_max_new_cells = global_sum_stats[0];
				// Statistics: Accumulated dead and Max dead cells per step
				sim_stat.history_dead_cells += global_sum_stats[1];
				if (global_sum_stats[1] > sim_stat.history_max_dead_cells)
					sim_stat.history_max_dead_cells = global_sum_stats[1];
				// Statistics: Max alive cells per step
				if (total_cells > sim_stat.history_max_alive_cells)
					sim_stat.history_max_alive_cells = total_cells;
			}

#ifdef DEBUG
			/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
			if (rank == 0)
				print_status(iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat);
#endif // DEBUG
		}

		// Last reduction (age):
		int max_age_root;
		MPI_Reduce(&sim_stat.history_max_age, &max_age_root, 1, MPI_INT, MPI_MAX, 0, universe);
		sim_stat.history_max_age = max_age_root;

		num_cells_alive = total_cells;

		// Let's not be bad...
		free(new_cells);
		free(cells_moved_to);
	}
#ifndef CP_TABLON
    if (rank == 0)
    {
        printf("Execution times for each subsection:\n");
        printf("\t3.1 - %lf\n", time3_1);
        printf("\t3.2 - %lf\n", time3_2);
        printf("\t4.1 - %lf (max: %lf)\n", sum_time4_1, max_time4_1);
        printf("\t4.3 - %lf (max: %lf)\n", sum_time4_3, max_time4_3);
        printf("\t4.X - %lf (max: %lf)\n", sum_time4_X, max_time4_X);
        printf("\t4.4 - %lf (max: %lf)\n", sum_time4_4, max_time4_4);
        printf("\t4.7 - %lf (max: %lf)\n", sum_time4_7, max_time4_7);
        printf("\t4.8 - %lf (max: %lf)\n", sum_time4_8, max_time4_8);
        printf("\t4.9 - %lf (max: %lf)\n", sum_time4_9, max_time4_9);
    }
#endif
    /*

==================================================================== POWERED BY ======================================================================

MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNOO0NMMMMMMMMMMMMMMMMMMMMMMMMMWNNWWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNOO0NMMMMMMMMMMMMMMMMMMMMMMMMWKOOKNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNOO0NMMMMMMMMMMMMMMMMMMMMMMNK0OOKNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNOO0NMMMMMMMMMMMMMMMMMMMMWX0OOOKWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNOO0NMMMMMMMMMMMMMMMMMMMNKOOO0XWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWNXK0OkkkkkOkkkkkO0KXNWMMMMMMMMWN0OOKXNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWX0kxoc::;,,,,cxOxc,,,;;::loxO0XWMMWXOOOKWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMWXOxl:;,,,,,,,,,,,cxOxc,,,,,,,,,,,;cok00OO0XWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMWKko:;,,,,,,,,,,,,,,,cxOxc,,,,,,,,,,,,,;lxkOKNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMWNOo:;,,,,,,,,,,,,,,,,,,cxOxc,,,,,,,,,,,,:okOOXWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMWXkl;,,,,,,,,,,,,,,,,,,,,,cxOxc,,,,,,,,,,;cdkO0NWMMMNOxkKMMMMMMMW0x0WMMW0x0WMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMXkc;,,,,,,,,,,,,,,,,,,,,,,;lkOxc,,,,,,,,,;lxkOKNMMMMMKl,;kWMMMMMMNo;dNMMXo;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMNOl;,,,,,,,,,,,,,,,,,;:codxkOOOOOOOkxdoc:;:okO0XWMMMMMMKl,;kWMMMMMMNo;dNMMXo;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMXd:,,,,,,,,,,,,,,,;:lxOKNWMMMMNOk0NMMMMWNKkxkO0NMMMMMMMMKl,;kWMMMMMMNo;dNMMXo;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMWKl;,,,,,,,,,,,,,,:okKWMMMMMMMMMNOk0NMMMMMWX0OOKWMMMMMMMMMKl,;kWMMMMMMNo;dNMMXo;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMW0c;,,,,,,,,,,,,,:d0WMMMMMMMMMMMMNOk0NMMMMNKOO0XWMMMMMMMMMMKl,;kWMMMMMMNo;dNMMXo;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMW0c,,,,,,,,,,,,,:dKWMMMMMMMMMMMMMMNOO0NMMMN0OO0NMMMMMMMMMMMMXl,;xWMMMMMMXo;xWMMXo;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMKl,,,,,,,,,,,,;cONMMMMMMMMMMMMMMMMNOO0NMWX0OOKWMMMMMMMMMMMMMWkc:l0WMMMMNk:c0MMMNo;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMNd;,,,,,,,,,,,;lKWMMMMMMMMMMMMMMMMMNOO0NNKOO0XWMMMMMMMMMMMMMMMNKklcoxkkxoco0WMMMNd;xNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
WWWWWWWWWWWWk:,,,,,,,,,,,,l0WWWWWWWWWWWWWWWWWWWXOkOK0OO0XWWWWWWWWWWWWWWWWWWWWX0kxxxxOKNWWWWWN0kKNWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
00000000000kl,,,,,,,,,,,,:dO0000000000000000000OOOOOOOOO0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
OOOOOOOOOOOd:,,,,,,,,,,,;lkOOOOOOOOOOOOOOOOOOOOOOOkOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
NNNNNNNNNNXx;,,,,,,,,,,,cONNNNNNNNNNNNNNNNNNNNX0OkOKXNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
xxxxxxxxxxxl;,,,,,,,,,,,:dxxxxxxxxxxxxxxxxxxxk0OOOkkxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;cdkOOOxc,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
dddddddddddc;,,,,,,,,,,,:lddddddddddddddddooxkOOOOkxdddddddddddddddddddddo:,,,,,,,,,,,,;codddddddddddodddddddddddddddddddddddddddddddddddddddddddddddd
MMMMMMMMMMWO:,,,,,,,,,,,c0WMMMMMMMMMMMMMMNKOOOXXOO0NMMMMMMMMMMMMMMMMMMMMMKl,,,,,,,,,,,,cONNNNWWWNXNNWNXNWMNKNWXXNXNWNXNWWWWNNNWNNNNNNXXKXWNNNWNNXXNNNN
MMMMMMMMMMMKl,,,,,,,,,,,;xNMMMMMMMMMMMMMN0OO0NWNOk0NMMMMMMMMMMMMMMMMMMMMWk:,,,,,,,,,,,,oK00KKNKK00K000K0XMN0XN00000KKKKKKK0OO000XKO0KK0OKX0OOK0KKK00K0
MMMMMMMMMMMNd;,,,,,,,,,,,c0WMMMMMMMMMMWX0OOKWMMNOO0NMMMMMMMMMMMMMMMMMMMMKl,,,,,,,,,,,,;kNXNNKK0K0OKKKKKKNMWXKKKXXNNXKXN00NX0KXXXWNK0KXXKXX000KXNNXKKKK
MMMMMMMMMMMM0c,,,,,,,,,,,;lKMMMMMMMMMNKOO0XWMMMNOO0NMMMMMMMMMMMMMMMMMMMXo;,,,,,,,,,,,,lKMMMMMWWWNNWMMWWMMMMMWWMMMMMMMMMMWMMWWWMMMMWWMMMMWWWWWMMMMMMWWW
MMMMMMMMMMMMWk:,,,,,,,,,,,;lKWMMMMMMN0OO0NMMMMMNOO0NMMMMMMMMMMMMMMMMMMKo;,,,,,,,,,,,,:OWMMWNNWWWMMMWNWWWWWXKNWWWMWWWWWWWWMMNNWNNWNWMWWMWWWMMMMMMMMMMMM
MMMMMMMMMMMMMNd;,,,,,,,,,,,;cOWMMMWXOkOKWMMMMMMNOO0NMMMMMMMMMMMMMMMMW0l;,,,,,,,,,,,,;xNMWKK00K00KWMX0KKKKKOOKKKKKKKK0K0KKXN0OK0OXKKXKNNX00NMMMMMMMMMMM
MMMMMMMMMMMMMMXo;,,,,,,,,,,,,:dKWNKOO0XWMMMMMMMNOO0NMMMMMMMMMMMMMMWXx:,,,,,,,,,,,,,;dNMMN0K000O0KNMNKKKNKKKKX0K00KNXKNKXKKKOk0K0XK00KNX0OkKMMMMMMMMMMM
MMMMMMMMMMMMMMMXd;,,,,,,,,,,,,;cdkkOKNMMMMMMMMMNOO0NMMMMMMMMMMMMWKxc;,,,,,,,,,,,,,:xNMMMMWNNWWNNWMMWWWWWWWWWWNNNWWWWWWWWWWWNNWWNNWWWNNWNNNWMMMMMMMMMMM
MMMMMMMMMMMMMMMMXx:,,,,,,,,,,,,,,;cxKWMMMMMMMMMNOO0NMMMMMMMMMWX0d:;,,,,,,,,,,,,,,:kNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMWOl;,,,,,,,,,,,,,,;cok0XWMMMMMNOO0NMMMMMWNKkdc;,,,,,,,,,,,,,,,;l0WMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMWXx:,,,,,,,,,,,,,,,,,;cldxkO0OxdxO0OOxdoc;;,,,,,,,,,,,,,,,,;ckXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMWKd:,,,,,,,,,,,,,,,,,,,,,;;;;;;;;;,,,,,,,,,,,,,,,,,,,,;;:xKWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMWKd:;,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;cxOKWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMWXkoc:;,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;:oOXWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMWKOO0KKkoc;,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,;cokKWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMN0OO0NWMMWX0koc:;,,,,,,,,,,,,,,,,,,,,,,,;:cox0XWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMWX0OOKNMMMMMMMMWNKOkxolc::;;;;;,;;;::clodkOKNWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMWKOO0XWMMMMMMMMMMMMMMMWWNXKK0kdddO0KKXNWWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMNXKKNMMMMMMMMMMMMMMMMMMMMMMMMNOO0NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMWMMMMMMMMMMMMMMMMMMMMMMMMMMNOO0NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNOO0NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

 */

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

    /* 5. Stop global time */
    MPI_Barrier(MPI_COMM_WORLD);
    ttotal = cp_Wtime() - ttotal;

#ifdef DEBUG
    printf("List of cells at the end of the simulation: %d\n\n", num_cells);
    for (i = 0; i < num_cells; i++)
    {
        printf("Cell %d, Alive: %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
               i,
               cells[i].alive,
               cells[i].pos_row,
               cells[i].pos_col,
               cells[i].mov_row,
               cells[i].mov_col,
               cells[i].choose_mov[0],
               cells[i].choose_mov[1],
               cells[i].choose_mov[2],
               cells[i].storage,
               cells[i].age);
    }
#endif // DEBUG

    /* 6. Output for leaderboard */
    if (rank == 0)
    {
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
               sim_stat.history_max_food);
    }

    /* 7. Free resources */
    free(culture);
    free(culture_cells);
    free(cells);

    /* 8. End */
    MPI_Finalize();
    return 0;
}
