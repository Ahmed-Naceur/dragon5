
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 31/07/10      */
/*****************************************/

#include <string.h>
#include "cle2000.h"

int_32 kdrdpr(lifo **my_iptrun, int_32 nentry, char (*hentry)[13])
{

/*     GAN-2000 SYSTEM: R.ROY (01/2000), VERSION 2.1 */

/*             *KDRDPR* IS THE MODULE FOR *PROCEDURE   * DECLARATIONS */
/*                      =0 IF NO ERROR */

/*      INPUT: *MY_IPTRUN* IS THE EXEC STRUCTURE POINTER (ALLOCATED) */
/*             *NENTRY* IS THE # OF LINKED LISTS AND FILES USED. */
/*             *HENTRY* NAMES OF EACH OBJECT <- LINKED LIST OR FILE. */
/*                      ( CHARACTER*12 HENTRY(NENTRY) ) */

/*     SYNTAX: */
/*              PROCEDURE *HENTRY(I=1,NENTRY)* ; */
/*                (DEFAULT VALUES, CHECK EXISTENCE) */

/*              PROCEDURE *HENTRY(I=1,NENTRY)* :: *DATA* ; */
/*                (USER DEFINED VALUES, ACCESS TO CLE-2000 COMPILER) */

   char *nomsub = "kdrdpr";
   int_32 ret_val = 0;
   int_32 ityp, nitma, lndata;
   float_32 flott;
   double_64 dflot;
   char text[73], messag[133], filenm[73], filinp[77], filobj[77];
   int_32 iloop1;
   
   if (nentry <= 0) goto L666;

   redget_c(&ityp, &nitma, &flott, text, &dflot);
   lndata = ityp != 10 && (ityp != 3 || strcmp(text, ";") != 0);
   if (lndata) {
      sprintf(messag, "%s: NOT DEVELOPED YET (RR)", nomsub);
      printf("%-132s\n", messag);
      ret_val = -666;
      goto L666;
   }
   for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
      int_32 iparam;
      lifo_node *my_node;
      
      my_node = clenode(my_iptrun, hentry[iloop1]);
      if (my_node == NULL) {
         printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
         ret_val = -665;
         goto L666;
      }
      iparam = my_node->type;
      if (lndata) {
         redget_c(&ityp, &nitma, &flott, text, &dflot);
         if (ityp == 3) {
            strcpy(filenm, text);
         } else {
            goto L8001;
         }
      } else {
         strcpy(filenm, my_node->name);
      }

      sprintf(filinp, "%s.c2m",filenm);
      sprintf(filobj, "%s.o2m",filenm);
      if (iparam == 1) {
/*       ONLY VERIFY IF *filobj* EXISTS */
         FILE *file;
         file = fopen(filobj, "r");
         if (file == NULL) {
            sprintf(messag, "%s: OBJECT FILE *%s* DOES NOT EXIST: MUST BE COMPILED", nomsub, filobj);
            printf("%-132s\n", messag);
            ret_val = -1;
            goto L666;
         } else {
            fclose(file);
         }
      } else {
/*       ONLY VERIFY IF *filinp* EXISTS */
         FILE *file;
         file = fopen(filinp, "r");
         if (file == NULL) {
            sprintf(messag, "%s: INPUT FILE *%s* DOES NOT EXIST",nomsub, filinp);
            printf("%-132s\n", messag);
            ret_val = -1;
            goto L666;
         } else {
            fclose(file);
         }
         my_node->type = -iparam;
         strcpy(my_node->OSname, filobj);
      }
   }

/* CAN WE FOUND *;* AT THE END OF THE SENTENCE ? */
   if (lndata) {
      redget_c(&ityp, &nitma, &flott, text, &dflot);
      if (ityp != 3 || strcmp(text, ";") != 0) goto L8002;
   }
L666:
   return ret_val;

L8001:
   sprintf(messag, "%s: INVALID TYPE IN *PROCEDURE* DATA.", nomsub);
   printf("%-132s\n", messag);
   ret_val = 8001;
   goto L666;
L8002:
   sprintf(messag, "%s: INVALID END  IN *PROCEDURE* DATA.", nomsub);
   printf("%-132s\n", messag);
   ret_val = 8002;
   goto L666;
} /* kdrdpr */
