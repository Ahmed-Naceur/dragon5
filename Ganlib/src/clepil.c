
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 15/05/09      */
/*****************************************/

#include <string.h>
#include "cle2000.h"

int_32 clepil(FILE *iredin, FILE *iwrite, kdi_file *iunito,
              int_32 (*dumcst)(char *, int_32 *, int_32 *, float_32 *, char *, double_64 *))
{

/*     CLE-2000 SYSTEM: R.ROY (03/1999), VERSION 2.1 */

/*             *CLEPIL* WILL PERFORM A SYNTACTICAL ANALYSIS */
/*                      AND COMPILE THE INPUT UNIT *IREDIN*. */
/*                      RESULT IS THE OBJECT D.A. UNIT *IUNITO* */
/*                      COMPILER COMMENTS ARE WRITTEN ON UNIT *IWRITE* */
/*                      WORDS ARE SEPARATED AND CLASSIFIED BY TYPES. */
/*                      EVERYTHING IS CHECKED FOR CORRECT EXECUTION. */

/*      INPUT: *IREDIN* IS THE INPUT  UNIT */
/*             *IWRITE* IS THE OUTPUT UNIT FOR COMPILER COMMENTS */
/*             *IUNITO* IS THE DIRECT ACCESS UNIT FOR OBJECT CODE */
/*             *DUMCST* IS THE EXTERNAL FUNCTION FOR *CLE-2000* CONSTANTS */

/*       NOTE: *CLEPIL* = 0 IF NO PROBLEM WAS ENCOUNTERED WHILE COMPILING */

   char *nomsub = "clepil";
   char *clistc[] = {"clelog", "clestk", "clexrf"};
   int_32 iretcd, istepc;
   int_32 ret_val = 0;

/* CONSTRUCT OBJECT FILE AND ANALYSE LOGIC */
   istepc = 0;
   iretcd = clelog(iredin, iwrite, iunito);
   if (iretcd != 0) goto L9002;

/* ADD CLE-2000 VARIABLES */
   istepc = 1;
   iretcd = clestk(iunito, iwrite, dumcst);
   if (iretcd != 0) goto L9002;
   istepc = 3;

/* X-REF CLE-2000 VARIABLES */
   istepc = 2;
   iretcd = clexrf(iunito, iwrite);
   if (iretcd != 0) goto L9002;

L666:
   return ret_val;

L9002:
   printf("! %s: ERROR CODE IN >>%s<< ERROR NUMBER (%d)\n", nomsub, clistc[istepc], (int)iretcd);
   ret_val = iretcd;
   goto L666;

} /* clepil */
