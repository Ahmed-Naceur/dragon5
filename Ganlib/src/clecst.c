
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 15/05/09      */
/*****************************************/

#include <string.h>
#include "cle2000.h"
#define index_f(A, B) (strstr(A, B) == NULL ? 0 : strstr(A, B) - A + 1)

int_32 clecst(char *cparm, int_32 *ityp, int_32 *nitma, float_32 *flott, char *text, double_64 *dflot)
{
   static char *ctypes[]  = {"_I", "_R", "_S", "_D", "_L"};
   static char *cinteg[]  = {"$Version", "$XLangLvl", "$c0", "$Date", "$Time", "$True", "$False"};
   static int_32 dinteg[]   = { 2,77,299792458,20000101,0,1,-1 };
   static char *cfloat[]  = {"$Pi", "$E", "$Euler", "$c0", "$Na", "$u", "$eV", "$h"};
   static double_64 dfloat[] = {3.141592653589793, 2.718281828459045, .577215664901533, 299792458.,
                                6.02214199e23, 1.66053873e-27, 1.602176462e-19, 6.62606876e-34};
   static char *cstrin[]  = {"$Code", "$Release", "$XLang", "$Date", "$Time", "$Bang", "$GetIn",
                             "$GetOut"};
   static char *dstrin[]  = {"CLE2000", "3", "Fortran", "20000101", "000000",  "!", ">>",
                             "<<"};

/*     CLE-2000 CONSTANTS: R.ROY (11/1999) */

/*             *CLECST* WILL ATTEMPT TO FIND VALUES FOR */
/*                      A CLE-2000 CONSTANT. */

/*      INPUT: *CPARM * IS THE TENTATIVE CONSTANT NAME */

/*     OUTPUT: *ITYP  * IS THE CONSTANT TYPE (1:INTEGER) */
/*                                           (2:REAL) */
/*                                           (3:CHARACTER STRING) */
/*                                           (4:DOUBLE) */
/*                                           (5:LOGICAL) */
/*             *NITMA * IS AN INTEGER VALUE */
/*                      (= LENGTH OF STRING IF *ITYP* .EQ.3) */
/*                      (= -1 FOR .F. +1 FOR .T. IF *ITYP* .EQ.5) */
/*             *FLOTT * IS AN REAL VALUE */
/*             *TEXT  * IS AN CHARACTER STRING */
/*             *DFLOT * IS AN DOUBLE PRECISION VALUE */

/*       NOTE: *CLECST* = 0, IF WE FOUND THE PARAMETER *CPARM* */

/*             THIS FUNCTION DEPEND ON THE APPLICATION */
/*             THE EXAMPLE GIVEN HERE SHOULD HELP THE DEVELOPER */
/*             TO WRITE ITS OWN APPLICATION-BASED CONSTANT LIST. */

/*             PHYSICAL CONSTANTS GIVEN HERE WERE TAKEN FROM: */
/*             http://physics.nist.gov/cuu/Constants/ */

/*    EXAMPLE: HERE *ITYP* IS IMPOSED AT THE END OF *CPARM* */
/*             (1:INTEGER) => END: _I */
/*             (2:REAL   ) => END: _R */
/*             (3:STRING ) => END: _S */
/*             (4:DOUBLE ) => END: _D */
/*             (5:LOGICAL) => END: _L */
/*             ALL FLOATING (_R, _D ) ARE KEPT IN DOUBLE */
/*             AND CONVERTED INTO THE APPROPRIATE MODE. */

   int_32 iloop1, ret_val = 1;
   char cparin[13];

/* IDENTITY WHICH TYPE: _I ,_R, _D, _S, _L */
   int_32 indlec = 0;
   for (iloop1 = 0; iloop1 < 5; ++iloop1) {
      int_32 idftyp = index_f(cparm, ctypes[iloop1]);
      if (idftyp != 0) {
         indlec = iloop1 + 1;
         strncpy(cparin, cparm, idftyp-1); cparin[idftyp-1] = '\0';
      }
   }

   if (indlec == 1 || indlec == 5) {
/*    LOOK FOR INTEGER VARIABLES */
      for (iloop1 = 0; iloop1 < 7; ++iloop1) {
         if (strcmp(cparin, cinteg[iloop1]) == 0) {
/*          FOUND: RETURN => TYPE=1, INTEGER */
/*                        => TYPE=5, LOGICAL */
            *ityp = indlec;
            *nitma = dinteg[iloop1];
            ret_val = 0;
            goto L666;
         }
      }
   } else if (indlec == 3) {
/*    LOOK FOR STRING VARIABLES */
      for (iloop1 = 0; iloop1 < 8; ++iloop1) {
         if (strcmp(cparin, cstrin[iloop1]) == 0) {
/*          FOUND: RETURN => TYPE=3, STRING AND ITS LENGTH */
            *ityp = 3;
            strcpy(text, dstrin[iloop1]);
            *nitma = strlen(dstrin[iloop1]);
            ret_val = 0;
            goto L666;
         }
      }
   } else if (indlec != 0) {
/*    LOOK FOR FLOATING VARIABLES */
      for (iloop1 = 0; iloop1 < 8; ++iloop1) {
         if (strcmp(cparin, cfloat[iloop1]) == 0) {
            if (indlec == 2) {
/*             FOUND: RETURN => TYPE=2, REAL */
               *flott = (float_32) dfloat[iloop1];
            } else {
/*             FOUND: RETURN => TYPE=4, DOUBLE */
               *dflot = dfloat[iloop1];
            }
            *ityp = indlec;
            ret_val = 0;
            goto L666;
         }
      }
   }

L666:
   return ret_val;
} /* clecst */
