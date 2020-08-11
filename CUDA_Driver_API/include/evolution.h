/*
 * evolution:
 *
 * Constants, custom types, and macro functions used in the
 * evolution simulation program and its kernels.
 */

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
 * Get an uniformly distributed random number between 0 and max
 * It uses glibc_nrand, that returns a number between 0 and 2^31
 */
#ifdef __CUDACC__
#define int_urand48( max, seq )	(int)( (long)(max) * glibc_nrand48( seq ) / 2147483648 )
#else
#define int_urand48( max, seq )	(int)( (long)(max) * nrand48( seq ) / 2147483648 )
#endif // __CUDACC__

/* 
 * Macro function to simplify accessing with two coordinates to a flattened array
 * 	This macro-function can be modified by the students if needed
 *
 */
#define accessMat( arr, exp1, exp2 )	arr[ (int)(exp1) * columns + (int)(exp2) ]

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
