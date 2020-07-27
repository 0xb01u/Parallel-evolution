/* Copyright (C) 1995-2020 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@gnu.ai.mit.edu>, August 1995.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

/*
 * Modification to implement only the nrand function and
 * use it outside of glibc.
 *
 * (c) 2020, Arturo Gonzalez-Escribano
 * Designed for an assignment in a Parallel Computing course, 
 * Computer Engineering degree, Universidad de Valladolid (Spain). 
 * Academic year 2019/2020
 */

#include <stdlib.h>

__host__ __device__ void
glibc_rand48_iterate (unsigned short int xsubi[3])
{
  u_int64_t X;
  u_int64_t result;

  /* Do the real work.  We choose a data type which contains at least
     48 bits.  Because we compute the modulus it does not care how
     many bits really are computed.  */

  X = (u_int64_t) xsubi[2] << 32 | (u_int32_t) xsubi[1] << 16 | xsubi[0];

  result = X * 0x5deece66dull + 0xb;

  xsubi[0] = result & 0xffff;
  xsubi[1] = (result >> 16) & 0xffff;
  xsubi[2] = (result >> 32) & 0xffff;
}


__host__ __device__ long int
glibc_nrand48 (unsigned short int xsubi[3])
{
  long int result;

  /* Compute next state.  */
  glibc_rand48_iterate (xsubi);

  /* Store the result.  */
  if (sizeof (unsigned short int) == 2)
    result = xsubi[2] << 15 | xsubi[1] >> 1;
  else
    result = xsubi[2] >> 1;

  return result;
}

