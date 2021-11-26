
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 31/07/10      */
/*****************************************/

#include <stdlib.h>
#include <string.h>
#include "cle2000.h"

void drviox(lifo *my_iptdat, int_32 minput, int_32 *nusec2)
{
   char *nomsub = "drviox";
   static char *ctypes[] = {"_I", "_R", "_S", "_D", "_L"};

/*     GAN-2000 SYSTEM: R.ROY (01/2000), VERSION 2.1 */

/*             *DRVIOX* IS USED TO INPUT/OUTPUT CLE-2000 VALUES */
/*                                 INTO/FROM    DATA STRUCTURE. */

/*      INPUT: *IPTDAT* IS THE DATA STRUCTURE POINTER (ALLOCATED) */
/*             *MINPUT* IS AN INTEGER -1: TO READ DATA INPUT (IN MAIN) */
/*                                     0: TO GET THIS INPUT (IN PROC, AFTER "::") */
/*                                    +1: TO RETURN VALUES   (IN MAIN) */
/*             *NUSEC2* IS THE OFFSET OF NEXT DATA VALUE ENTERED AFTER "::" */

   int_32 ityp, nitma, ntypc2;
   float_32 flott;
   double_64 dflot;
   char text[73], messag[73];
   lifo_node *my_node;

   if (minput == -1) {
      int_32 ndatc2 = 0, nembed = 0;

/*    INPUT CLE-2000 VARIABLES FROM MAIN PROGRAM */
L10:
      redget_c(&ityp, &nitma, &flott, text, &dflot);
      if (ityp == 10) {
         xabort_c("DRVIOX: REDGET HITS THE ; (1).");
      } else if (ityp == 3 && strcmp(text, ";") == 0 && nembed == 0) {
/*       END OF STATEMENT REACHED */
         goto L666;
      } else {
         my_node = (lifo_node *) malloc(sizeof(lifo_node));
         clepush(&my_iptdat, my_node);
         strcpy(my_node->name, "_dummy");
         if (ityp > 0) {
            my_node->type = 10 + ityp;
            my_node->access = 1;
            if (ityp == 3) {
               if (strcmp(text, ":::") == 0) {
                  ++nembed;
               } else if (strcmp(text, ";") == 0) {
                  --nembed;
               } else if (strcmp(text, "::") == 0) {
                  xabort_c("DRVIOX: INPUT DATA MISTAKE (ACT).");
               }
               strcpy(my_node->value.hval, text);
            } else if (ityp == 1 || ityp == 5) {
               my_node->value.ival = nitma;
            } else if (ityp == 2) {
               my_node->value.fval = flott;
            } else if (ityp == 4) {
               my_node->value.dval = dflot;
            } else {
               xabort_c("DRVIOX: INVALID TYPE (ACT).");
            }
         } else {
            my_node->type = -10 + ityp;
            my_node->access = 0;
            strcpy(my_node->name, text);
         }
         ++ndatc2;
      }
      goto L10;
   } else if (minput == 0) {
/*    READ/WRITE CLE-2000 VARIABLES IN THE PROCEDURE CALL */
L20:
      redget_c(&ityp, &nitma, &flott, text, &dflot);
      if (ityp == 10) {
         xabort_c("DRVIOX: REDGET HITS THE ; (2).");
      } else if (ityp == 3 && strcmp(text, ";") == 0) {
         goto L666;
      } else {
         if (*nusec2 + 1 > my_iptdat->nup) {
            printf("%s: INVALID NUMBER OF PARAMETERS (nusec2=%d)\n", nomsub, (int)(*nusec2 + 1));
            sprintf(messag, "%s: PROC WAS CALLED WITH ONLY %d PARAMETERS.\n", nomsub, (int)my_iptdat->nup);
            xabort_c(messag);
         }
         if (ityp < 0) {
            my_node = clepos(&my_iptdat, *nusec2);
            ntypc2 = my_node->type - 10;
            if (-ityp != ntypc2) {
               if (ityp < 0) {
                  printf("%s: DUMMY VARIABLE NAME *%.12s* OF TYPE(%s)\n",
                         nomsub, text, ctypes[-ityp-1]);
               } else {
                  printf("%s: DUMMY VALUE OF TYPE(%s)\n", nomsub, ctypes[ityp-1]);
               }
               if (my_node->type < 0) {
                  strcpy(text, my_node->name);
                  printf("%s: ACTUAL VARIABLE NAME *%.12s* OF TYPE(%s)\n", nomsub, text, ctypes[-ntypc2-1]);
               } else {
                  printf("%s: DUMMY VALUE OF TYPE(%s)\n", nomsub, ctypes[ntypc2-1]);
               }
               xabort_c("DRVIOX: INVALID TYPE (DUMMY) 5.");
            }
            if (strcmp(my_node->name, "_dummy") == 0) strcpy(my_node->name, text);
            nitma = 0;
            flott = 0.f;
            dflot = 0.L;
            strcpy(text, " ");
            if (ityp == -1 || ityp == -5) {
               nitma = my_node->value.ival;
            } else if (ityp == -2) {
               flott = my_node->value.fval;
            } else if (ityp == -3) {
               strcpy(text, my_node->value.hval);
               nitma = strlen(text);
            } else if (ityp == -4) {
               dflot = my_node->value.dval;
            } else {
               xabort_c("DRVIOX: INVALID TYPE (DUMMY) 6.");
            }
            redput_c(&ntypc2, &nitma, &flott, text, &dflot);
         } else {
            my_node = clepos(&my_iptdat, *nusec2);
            ntypc2 = my_node->type + 10;
            if (-ityp != ntypc2) {
               if (ityp < 0) {
                  printf("%s: DUMMY VARIABLE NAME *%.12s* OF TYPE(%s)\n", nomsub, text, ctypes[-ityp-1]);
               } else {
                  printf("%s: DUMMY VALUE OF TYPE(%s)\n", nomsub, ctypes[ityp-1]);
               }
               if (my_node->type < 0) {
                  strcpy(text, my_node->name);
                  printf("%s: ACTUAL VARIABLE NAME *%.12s* OF TYPE(%s)\n", nomsub, text, ctypes[-ntypc2-1]);
               } else {
                  printf("%s: DUMMY VALUE OF TYPE(%s)\n", nomsub, ctypes[ntypc2-1]);
               }
               xabort_c("DRVIOX: INVALID TYPE (DUMMY) 5.");
            }
            if (ityp == 1 || ityp == 5) {
               my_node->value.ival = nitma;
            } else if (ityp == 2) {
               my_node->value.fval = flott;
            } else if (ityp == 3) {
               strcpy(my_node->value.hval, text);
            } else if (ityp == 4) {
               my_node->value.dval = dflot;
            } else {
               xabort_c("DRVIOX: INVALID TYPE (DUMMY) 7.");
            }
            my_node->type = -my_node->type;
         }
         ++(*nusec2);
      }
      goto L20;
   } else if (minput == 1) {
      int_32 iloop;

/*    CONSISTENT RETURN (NOW, *REDPUT* IN THE REVERSE ORDER) */
      for (iloop = my_iptdat->nup - 1; iloop >= 0; --iloop) {
         my_node = clepos(&my_iptdat, iloop);
         if (my_node->type < 10) continue;
         ityp = my_node->type - 10;
         nitma = 0;
         flott = 0.f;
         dflot = 0.L;
         strcpy(text, " ");
         if (my_node->access == 0) {
            if (ityp == 1 || ityp == 5) {
               nitma = my_node->value.ival;
            } else if (ityp == 2) {
               flott = my_node->value.fval;
            } else if (ityp == 3) {
               strcpy(text, my_node->value.hval);
               nitma = strlen(text);
            } else if (ityp == 4) {
               dflot = my_node->value.dval;
            } else {
               xabort_c("DRVIOX: INVALID TYPE (DUMMY) 4.");
            }
            redput_c(&ityp, &nitma, &flott, text, &dflot);
            my_node->access = 1;
         }
      }
   } else {
      xabort_c("DRVIOX: INVALID VALUE FOR *MINPUT* ARG");
   }
L666:
   return;
} /* drviox */
