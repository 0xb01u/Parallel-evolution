/*
 * int_float:
 *
 * Use 32-bit integers to represent real numbers with fixed point precision.
 * It is presented with examples using as precision a number of decimals 
 * in base 10, which improves human readability of the stored number.
 *
 * (c) 2020, Arturo Gonzalez-Escribano
 * Designed for an assignment in a Parallel Computing course, Computer 
 * Engineering degree, Universidad de Valladolid (Spain). 
 * Academic year 2019/2020
 */
#ifndef __INT_FLOAT
#define __INT_FLOAT

#include<math.h>

/*
 * Change this constat to change precision. 
 *
 * Be carefull with the range of numbers supported: +/- ( 2^31 ) / PRECISION
 * More precision means a narrower range in the integer part.
 * Example: The range of numbers represented with 
 * 	PRECISION = 100000 is [ -21474.83648, 21474.83648 ]
 */
#define	PRECISION	100000

/* 
 * Common transformation operations
 * 	Transform an integer or float number into a fixed precision integer: 
 * 	(int)(num * PRECISION)
 *
 * 	Transform a fixed precision integer into a rounded-down normal integer 
 * 	using integer division directly: 
 * 	(num / PRECISION)
 */
#define float2int( f )	(int)( f * PRECISION )
#define int2float( i )	( i / PRECISION )

/*
 * PI and 2*PI constants as a fixed precision integer
 */
#define INT_PI	(int)(M_PI * PRECISION)
#define INT_2PI	(int)(2 * M_PI * PRECISION)

/*
 * Operators:
 *
 * 1) Additive operators +,- work as usual
 * 2) Multipliction or division by constants work as usual
 * 3) Multiplication of two fixed point integers: 
 * 	Both numbers are stored multiplied by the precision constant
 *	When the numbers are multiplied, the precision constant is duplicated
 *	Example: 2.0 * 3.0 = 6.0
 *		(2.0*PRECISION) * (3.0*PRECISION) = 6.0 * PRECISION * PRECISION
 * 	Eliminate the extra PRECISION constant dividing the result by PRECISION
 * 	To avoid oveflows, use long type for the partial result of the 
 * 	multiplication before eliminating the extra PRECISION constant.
 */
#define intfloatMult( a, b )	(int)((long)(a) * (long)(b) / PRECISION)

#endif // __INT_FLOAT
