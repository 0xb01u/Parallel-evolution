/*
 * cuda_check:
 *
 * Simple macro functions to check errors on (a) calls to CUDA library, 
 * or (b) after a kernel execution.
 *
 * Examples of use:
 * 	cudaCheckCall( cudaSetDevice(1) );
 *
 * 	kernel<<< grid, block >>>( params );
 * 	cudaCheckLast();
 *
 * The macro functions are ignored and no check is done
 * if the macro CUDA_CHECK is not defined.
 *
 * (c) 2020, Arturo Gonzalez-Escribano
 * To be used in an assignment of the Parallel Computing Course,
 * Computer Engineering degree, Universidad de Valladolid (Spain)
 * Academic year 2019/2020
 *
 * Students can freely modify these macro functions for their purposes 
 * when working in the assignment. 
 */

#ifdef CUDA_CHECK

// If macro CUDA_CHECK is defined do error checks
#define cudaCheckCall( call )	{ \
			cudaError_t ok = call; \
			if ( ok != cudaSuccess ) \
				fprintf(stderr, "-- Error CUDA call in line %d: %s\n", __LINE__, cudaGetErrorString( ok ) ); \
			}

#define cudaCheckLast()	{ \
			cudaError_t ok = cudaGetLastError(); \
			if ( ok != cudaSuccess ) \
				fprintf(stderr, "-- Error CUDA last in line %d: %s\n", __LINE__, cudaGetErrorString( ok ) ); \
			}
#else

// If macro CUDA_CHECK is not defined skip checks
#define cudaCheckCall( call )	call
#define cudaCheckLast()

#endif

