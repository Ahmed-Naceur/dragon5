
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 15/05/09      */
/*****************************************/

#include <string.h>
#include "cle2000.h"

int_32 objpil(kdi_file *iunito, FILE *iwrite, int_32 ldatav)
{

/*     GAN-2000 SYSTEM: R.ROY (12/1999), VERSION 2.1 */

/*             *OBJPIL* WILL COMPLETE SYNTACTICAL ANALYSIS */
/*                      INCLUDING OBJECT CONSISTENCE IN GAN-2000 */
/*                      RESULT IS THE OBJECT D.A. UNIT *IUNITO* */
/*                      COMPILER COMMENTS ARE WRITTEN ON UNIT *IWRITE* */
/*                      EVERYTHING IS CHECKED FOR CORRECT EXECUTION. */

/*      INPUT: *IUNITO* IS THE DIRECT ACCESS UNIT FOR OBJECT CODE */
/*             *IWRITE* IS THE OUTPUT UNIT FOR COMPILER COMMENTS */
/*             *LDATAV* =0/1: PROCEDURE SECTION/DATA SECTION */

/*       NOTE: *OBJPIL* = 0 IF NO PROBLEM WAS ENCOUNTERED WHILE COMPILING */

    char *nomsub = "objpil";
    char *clistc[] = {"objstk", "objxrf"};
    int_32 iretcd, istepc;
    int_32 ret_val = 0;

/*  ADD OBJECTS/MODULES TO CLE-2000 FILES */
    istepc = 0;
    iretcd = objstk(iunito, iwrite, ldatav);
    if (iretcd != 0) goto L9002;

/*  X-REF OBJECTS/MODULES */
    istepc = 1;
    iretcd = objxrf(iunito, iwrite);
    if (iretcd != 0) goto L9002;

L666:
    return ret_val;

L9002:
    printf("! %s: ERROR CODE IN >>%s<< ERROR NUMBER (%d)\n", nomsub, clistc[istepc], (int)iretcd);
    ret_val = iretcd;
    goto L666;

} /* objpil */
