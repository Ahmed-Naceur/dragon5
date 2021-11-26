
/**********************************/
/* C API for kdi file support     */
/* author: A. Hebert (08/04/2005) */
/**********************************/

/*
Copyright (C) 2005 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/

#include <stdio.h>
#include "ganlib.h"

typedef struct {
   char nom[73];          /* FILE NAME */
   FILE *fd;              /* FILE POINTER */
} kdi_file;

kdi_file * kdiop_c(char *, int_32);
int_32 kdiput_c(kdi_file *, int_32 *, int_32, int_32);
int_32 kdiget_c(kdi_file *, int_32 *, int_32, int_32);
int_32 kdicl_c(kdi_file *, int_32);
