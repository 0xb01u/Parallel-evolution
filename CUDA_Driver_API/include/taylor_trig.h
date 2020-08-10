/* 
 * Portable trigonometric functions with fixed point arithmetics
 * Sine approximated by Taylor polynomial of degree 7
 *
 * sin(a) = a - a^3/3! + a^5/5! - a^7/7!
 *
 * Taylor approximation works in range -PI:+PI
 * Input: Angle in radians, valid range -PI:3*PI
 *
 * (c) 2020, Arturo Gonzalez-Escribano
 * Simple implementation for an assignment in a Parallel Computing course, 
 * Computer Engineering degree, Universidad de Valladolid (Spain). 
 * Academic year 2019/2020
 *
 * This approximation is good enough for the purpose of the simulation
 * proposed in the assignment
 *
 */ 
#include "int_float.h"


__host__ __device__ int taylor_sin( int angle ) {
	if ( angle > INT_PI ) angle = -INT_2PI + angle;
	long angle3 = ( (long)angle * (long)angle / PRECISION ) * (long)angle / PRECISION;
	long angle5 = ( angle3 * (long)angle / PRECISION ) * (long)angle / PRECISION;
	long angle7 = ( angle5 * (long)angle / PRECISION ) * (long)angle / PRECISION;
	return (int)(angle - angle3 / 6 + angle5 / 120 - angle7 / 5040 );
}

__host__ __device__ int taylor_cos( int angle ) {
	int turnedAngle = angle + INT_PI/2;
	return taylor_sin( turnedAngle );
}

