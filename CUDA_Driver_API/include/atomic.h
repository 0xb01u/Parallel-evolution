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
