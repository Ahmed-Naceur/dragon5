/*
   -----------------------------------------------------------------------
   Copyright (C) 2002 Ecole Polytechnique de Montreal

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   DYNAMIC ALLOCATION/LIBERATION OF MEMORY WITH SETARA/RLSARA IN C.

   SETARA PARAMETERS:
   LENGTH : NUMBER OF SINGLE PRECISION WORDS TO BE ALLOCATED BY SETARA.
    IOF : OFFSET OF THE ALLOCATED SPACE.

   RLSARA PARAMETER:
   BASE : FIRST WORD OF THE MEMORY SPACE TO BE DEALLOCATED.

   -----------------------------------------------------------------------
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ganlib.h"

static char AbortString[80];

int_32 * setara_c(int_32 length)
{
   char *nomsub="setara_c";
   int_32  *kaddr;

   kaddr = (int_32 *)malloc(length*sizeof(int_32));
   if (kaddr == 0) {
      sprintf(AbortString,"%s: invalid return code %ld length=%d",nomsub,(long)kaddr,(int)length);
      xabort_c(AbortString);
   }
   return kaddr;
}
/* Function rlsara_c */
void rlsara_c(int_32 *iof)
{
   free(iof);
   return;
}
