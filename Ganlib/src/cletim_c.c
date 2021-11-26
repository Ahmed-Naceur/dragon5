/*****************************************/
/*             CLE-2000 API              */
/*     AUTHOR: A. Hebert ; 16/07/10      */
/*****************************************/

#include <stdlib.h>
#include <time.h>
#include "ganlib.h"

void cletim_c(int_32 *sec, int_32 *nsec){
   int_32 value = (int_32) clock();
   *sec = value / CLOCKS_PER_SEC;
   *nsec = (1000000 / CLOCKS_PER_SEC) * (value % CLOCKS_PER_SEC);
}
