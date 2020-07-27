/*
 * cuda_time:
 *
 * Simple macro functions to time executions of code sections
 * in a CUDA program.
 *
 * The macro functions are ignored and no check is done
 * if the macro DEVELOPMENT is not defined.
 *
 * (c) 2020, Manuel de Castro Caballero
 * To be used in an assignment of the Parallel Computing Course,
 * Computer Engineering degree, Universidad de Valladolid (Spain)
 * Academic year 2019/2020
 */
#ifndef max
#define max(x, y) (x > y ? x : y)
#endif

#ifdef DEVELOPMENT

float time3_1 = .0f;
float time3_2 = .0f;
float time4_1 = .0f;
float time4_2 = .0f;
float time4_3 = .0f;
float time4_4 = .0f;
float time4_5 = .0f;
float time4_6 = .0f;
float time4_7 = .0f;
float time4_8 = .0f;
float time4_9 = .0f;

float sum_time4_1 = .0f;
float sum_time4_2 = .0f;
float sum_time4_3 = .0f;
float sum_time4_4 = .0f;
float sum_time4_5 = .0f;
float sum_time4_6 = .0f;
float sum_time4_7 = .0f;
float sum_time4_8 = .0f;
float sum_time4_9 = .0f;

float max_time4_1 = .0f;
float max_time4_2 = .0f;
float max_time4_3 = .0f;
float max_time4_4 = .0f;
float max_time4_5 = .0f;
float max_time4_6 = .0f;
float max_time4_7 = .0f;
float max_time4_8 = .0f;
float max_time4_9 = .0f;

cudaEvent_t start, stop;
cudaError_t errSync;

#define time_start() { \
	cudaEventCreate(&start); \
	cudaEventCreate(&stop); \
	cudaEventRecord(start, NULL); \
}

#define time_end(timer) { \
	errSync = cudaDeviceSynchronize(); \
	if (errSync != cudaSuccess) printf("Synchronization error: %s\n", cudaGetErrorString(errSync)); \
	cudaEventRecord(stop, NULL); \
	cudaEventSynchronize(stop); \
	cudaEventElapsedTime(&timer, start, stop); \
}

#define update_times() { \
	sum_time4_1 += time4_1; \
	sum_time4_2 += time4_2; \
	sum_time4_3 += time4_3; \
	sum_time4_4 += time4_4; \
	sum_time4_5 += time4_5; \
	sum_time4_6 += time4_6; \
	sum_time4_7 += time4_7; \
	sum_time4_8 += time4_8; \
	sum_time4_9 += time4_9; \
	max_time4_1 = max(max_time4_1, time4_1); \
	max_time4_2 = max(max_time4_2, time4_2); \
	max_time4_3 = max(max_time4_3, time4_3); \
	max_time4_4 = max(max_time4_4, time4_4); \
	max_time4_5 = max(max_time4_5, time4_5); \
	max_time4_6 = max(max_time4_6, time4_6); \
	max_time4_7 = max(max_time4_7, time4_7); \
	max_time4_8 = max(max_time4_8, time4_8); \
	max_time4_9 = max(max_time4_9, time4_9); \
	time4_1 = .0f; \
	time4_2 = .0f; \
	time4_3 = .0f; \
	time4_4 = .0f; \
	time4_5 = .0f; \
	time4_6 = .0f; \
	time4_7 = .0f; \
	time4_8 = .0f; \
	time4_9 = .0f; \
}

#define print_times() { \
	printf("Execution times for each subsection:\n"); \
    printf("\t3.1 - %lf\n", time3_1); \
    printf("\t3.2 - %lf\n", time3_2); \
    printf("\t4.1 - %lf (max: %lf)\n", sum_time4_1, max_time4_1); \
    printf("\t4.2 - %lf (max: %lf)\n", sum_time4_2, max_time4_2); \
    printf("\t4.3 - %lf (max: %lf)\n", sum_time4_3, max_time4_3); \
    printf("\t4.4 - %lf (max: %lf)\n", sum_time4_4, max_time4_4); \
    printf("\t4.5 - %lf (max: %lf)\n", sum_time4_5, max_time4_5); \
    printf("\t4.6 - %lf (max: %lf)\n", sum_time4_6, max_time4_6); \
    printf("\t4.7 - %lf (max: %lf)\n", sum_time4_7, max_time4_7); \
    printf("\t4.8 - %lf (max: %lf)\n", sum_time4_8, max_time4_8); \
    printf("\t4.9 - %lf (max: %lf)\n", sum_time4_9, max_time4_9); \
}
#else
#define time_start()
#define time_end(timer)
#define update_times()
#define print_times()
#endif