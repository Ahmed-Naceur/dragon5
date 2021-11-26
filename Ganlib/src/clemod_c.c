
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 19/06/09      */
/*****************************************/

/* Call a single module without a CLE-2000 procedure */

#include <string.h>
#include "cle2000.h"
int_32 clemod_c(char *cmodul, FILE *filein, int_32 nentry, char (*hentry)[13], int_32 *ientry,
                int_32 *jentry, lcm **kentry, char (*hparam)[73],
                int_32 (*dummod)(char *, int_32, char (*)[13], int_32 *, int_32 *, lcm **, char (*)[73]))
{
   char *nomsub = "clemod_c";
   int_32 ret_val = 0;
   FILE *jwrite;
   char hsmg[132], filenm[8];
   int_32 iretcd, jrecin;
   kdi_file *iKDI;
   char hwrite[73] = " ";

/* first step, initialize stuff and compile main */
   sprintf(filenm,"_FIL%.3d",0);
   iKDI = kdiop_c(filenm,0);
   if (iKDI == NULL) {
      sprintf(hsmg, "%s: kdiop failure\n", nomsub);
      printf("%s\n", hsmg);
      ret_val = -1;
      goto L10;
   }

/* compile main input into object file */
   iretcd = clepil(filein, stdout, iKDI, clecst);
   if (iretcd != 0) {
      sprintf(hsmg, "%s: COMPILING _MAIN.c2m FILE (ERROR CODE) IRC=%d\n", nomsub,(int)iretcd);
      printf("%s\n", hsmg);
      ret_val = -2;
      goto L10;
   }

/* add objects/modules to object file */
   iretcd = objpil(iKDI, stdout, 1);
   if (iretcd != 0) {
      sprintf(hsmg, "%s: BAD OBJECTS _MAIN.c2m FILE (ERROR CODE) IRC=%d\n", nomsub,(int)iretcd);
      printf("%s\n", hsmg);
      ret_val = -3;
      goto L10;
   }

/* execute a module of the software application */
   redopn_c(iKDI, stdout, hwrite, 0);
   fflush(stdout);
   if (strcmp(cmodul, "END:") == 0) {
      printf("%s: dummy END: module called\n", nomsub);
      ret_val = 0;
   } else {
      iretcd = (*dummod)(cmodul, nentry, hentry, ientry, jentry, kentry, hparam);
      if (iretcd != 0) {
         sprintf(hsmg, "%s: calculation module failure IRC=%d\n", nomsub,(int)iretcd);
         printf("%s\n", hsmg);
         ret_val = -4;
         goto L10;
      }
   }

/* close the REDGET input reader */
   redcls_c(&iKDI, &jwrite, hwrite, &jrecin);
   iretcd = kdicl_c(iKDI, 2);
   if (iretcd != 0) {
      sprintf(hsmg, "%s: kdicl failure IRC=%d\n", nomsub,(int)iretcd);
      printf("%s\n", hsmg);
      ret_val = -5;
   }
L10:
   fflush(stdout);
   return ret_val;
}
