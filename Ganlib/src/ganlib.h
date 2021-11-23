
/**********************************/
/* C API for Ganlib5 support      */
/* author: A. Hebert (31/05/2009) */
/**********************************/

/*
Copyright (C) 2009 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/

typedef float float_32;
typedef double double_64;

#if __LP64__  || __64BIT__
typedef int int_32;
#else
typedef long int_32;
#endif

void xabort_c(char *);
int_32 * setara_c(int_32);
void rlsara_c(int_32 *);
