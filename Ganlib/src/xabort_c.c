
/*****************************************/
/*               GANLIB API              */
/*     AUTHOR: A. Hebert ; 06/05/09      */
/*****************************************/

/*
Copyright (C) 2009 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/

#include <stdlib.h>
#include <stdio.h>

void xabort_c(char *msg){
  printf("%s\n",msg);
  exit(1);
}
