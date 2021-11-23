
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 31/07/10      */
/*****************************************/

#include <stdlib.h>
#include <string.h>
#include "cle2000.h"
#define mpara2 64
#define ndclkw 8

int_32 kdrprm(lifo **my_iptrun, lifo **my_param, int_32 minput, int_32 nentry, int_32 *jentry, char (*hentry)[13])
{
   char *nomsub = "kdrprm";
   static char *cdclkw[] = {"PROCEDURE", "MODULE", "LINKED_LIST",  "XSM_FILE",
                            "SEQ_BINARY", "SEQ_ASCII", "DIR_ACCESS", "PARAMETER"};

/*     GAN-2000 SYSTEM: R.ROY (01/2000), VERSION 2.1 */

/*             *KDRPRM* IS USED TO PASS DUMMY ARGUMENTS */
/*                              IN CLE-2000 PROCEDURES. */

/*      INPUT: *MY_IPTRUN* IS THE EXEC STRUCTURE POINTER (ALLOCATED) */
/*             *MY_PARAM*  IS THE PARAMETER STRUCTURE POINTER (ALLOCATED) */
/*             *MINPUT* IS AN INTEGER -1: TO READ PARM INPUT (IN MAIN) */
/*                                     0: TO GET  PARM INPUT (IN PROC) */
/*                                    +1: TO RETURN VALUES   (IN MAIN) */
/*             *NENTRY* IS THE # OF LINKED LISTS AND FILES USED. */
/*             *HENTRY* NAMES OF EACH OBJECT <- LINKED LIST OR FILE. */
/*                      ( CHARACTER*12 HENTRY(NENTRY) ) */

   int_32 ret_val = 0;
   int_32 iloop1, iparam, jparam, lparam=0, nusepr;
   lcm *kparam = NULL;
   char hparam[73], messag[133], AbortString[132];
   char hwrite[73] = " ";
   lifo_node *my_node;

   if (minput == -1) {
/*    PUT OBJECTS INTO *IPTRUN* BEFORE CALLING A PROCEDURE */
      cleopn(my_param);
      if (nentry != 0) {
         for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
            my_node = clenode(my_iptrun, hentry[iloop1]);
            if (my_node == NULL) {
               printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
               ret_val = -20;
               goto L666;
            }
            iparam = my_node->type;
            jparam = my_node->access;
            if (iparam == 3) kparam = my_node->value.mylcm;
            if (iparam == 7) lparam = my_node->lparam;
            if (abs(iparam) <= 10) {
               if (strlen(my_node->OSname) > 72) {
                 sprintf(AbortString,"%s: OSname EXCEEDING 72 CHARACTERS.",nomsub);
                 xabort_c(AbortString);
               }
               strcpy(hparam, my_node->OSname);
            } else {
               strcpy(hparam, " ");
            }

            my_node = (lifo_node *) malloc(sizeof(lifo_node));
            strcpy(my_node->name, hentry[iloop1]);
            my_node->type = iparam;
            my_node->access = jentry[iloop1];
            if (iparam == 3) my_node->value.mylcm = kparam;
            if (iparam == 7) my_node->lparam = lparam;
            strcpy(my_node->OSname, hparam);
            clepush(my_param, my_node);
         }
      }

   } else if (minput == 0) {
      int_32 npara2, l2data=0, jrecin, iretcd;
      FILE *jwrite;
      char hpara2[mpara2][13];
      kdi_file *jread;

/*    LINK DUMMY OBJECTS WITH THEIR ACTUAL ARGUMENTS */
      int_32 ityp, nitma, lndata;
      float_32 flott;
      double_64 dflot;
      char text[73], cmodul[13], aparam[13];
      int_32 iprint = 0;

      if (nentry <= 0) {
         sprintf(messag, "%s: *PARAMETER* WITHOUT OBJECTS", nomsub);
         printf("%-132s\n", messag);
         ret_val = -1;
         goto L666;
      }
      redget_c(&ityp, &nitma, &flott, text, &dflot);
      lndata = ityp != 10 && (ityp != 3 || strcmp(text, ";") != 0);
      if (lndata && ityp == 3 && strcmp(text, "EDIT") == 0) {
         redget_c(&ityp, &iprint, &flott, text, &dflot);
         if (ityp != 1 && iprint < 0) {
            sprintf(messag, "%s: AFTER *EDIT*, PUT A POSITIVE INT", nomsub);
            printf("%-132s\n", messag);
            ret_val = -2;
            goto L666;
         }
         redget_c(&ityp, &nitma, &flott, text, &dflot);
         lndata = ityp != 10 && (ityp != 3 || strcmp(text, ";") != 0);
      }
      for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
         my_node = clepos(my_param, iloop1);
         strcpy(my_node->name_daughter, hentry[iloop1]);
         if (strlen(my_node->name) > 12) {
           sprintf(AbortString,"%s: name EXCEEDING 12 CHARACTERS.",nomsub);
           xabort_c(AbortString);
         }
         strcpy(aparam, my_node->name);
         iparam = my_node->type;
         jparam = my_node->access;
         if (iparam == 3) kparam = my_node->value.mylcm;
         if (iparam == 7) lparam = my_node->lparam;
         if (abs(iparam) <= 10) {
            if (strlen(my_node->OSname) > 72) {
              sprintf(AbortString,"%s: OSname EXCEEDING 72 CHARACTERS.",nomsub);
              xabort_c(AbortString);
            }
            strcpy(hparam, my_node->OSname);
         } else {
            strcpy(hparam, " ");
         }

         my_node = clenode(my_iptrun, hentry[iloop1]);
         if (my_node == NULL) {
            printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
            ret_val = -21;
            goto L666;
         }
         my_node->type = iparam;
         my_node->access = jparam;
         if (iparam == 3) my_node->value.mylcm = kparam;
         if (iparam == 7) my_node->lparam = lparam;
         strcpy(my_node->OSname, hparam);
         if (iprint > 0) {
            printf("PARAMETER %s <= %s with name(*%s*)\n", hentry[iloop1], aparam, hparam);
            if (iparam <= 0) {
               printf("          %s UNDEFINED.\n", cdclkw[-iparam-1]);
            } else {
               printf("          %s DEFINED.\n", cdclkw[iparam-1]);
            }
         }
      }

/*    NOW, CALL THE EMBEDDED DECLARATION MODULE IF DATA */
L21:
/*    MAKES IT A *DO WHILE* STRUCTURE */
      if (!lndata) goto L666;
      if (ityp != 3 || strcmp(text, ":::") != 0) {
         sprintf(messag, "%s: INVALID EMBEDDED DECL MODUL DATA", nomsub);
         printf("%-132s\n", messag);
         ret_val = -4;
         goto L666;
      }

/*    GET DECLARATION MODULE NAME */
      redget_c(&ityp, &nitma, &flott, text, &dflot);
      nusepr = 0;
      for (iloop1 = 0; iloop1 < ndclkw; ++iloop1) {
         if (strcmp(text, cdclkw[iloop1]) == 0) nusepr = iloop1 + 1;
      }
      if (nusepr == 0) {
         sprintf(messag, "%s: INVALID DECLARATION KEYWORD", nomsub);
         printf("%-132s\n", messag);
         ret_val = -5;
         goto L666;
      }
      strcpy(cmodul, cdclkw[nusepr-1]);

/*    GET OBJECT LIST */
      npara2 = 0;
      for (iloop1 = 0; iloop1 < nentry + 1; ++iloop1) {
         redget_c(&ityp, &nitma, &flott, text, &dflot);
         if (ityp == 3) {
            if (strcmp(text, ";") == 0) {
               l2data = 0;
               redcls_c(&jread, &jwrite, hwrite, &jrecin);
               break;
            } else if (strcmp(text, "::") == 0) {
               l2data = 1;
               break;
            } else {
/*             READER'S NAME *MUST* BE ONE OF *HENTRY* */
               int_32 jloop2;
               for (jloop2 = 1; jloop2 <= nentry; ++jloop2) {
                  if (strcmp(text, hentry[jloop2-1]) == 0) {
                     ++npara2;
                     strcpy(hpara2[npara2-1], hentry[jloop2-1]);
                     goto L25;
                  }
               }
               if (iloop1 != nentry) {
                  sprintf(messag, "%s: OBJECT *%s* NOT IN THE INPUT LIST", nomsub, text);
                  ret_val = -7;
               } else {
                  sprintf(messag, "%s: TOO MANY OBJECTS IN EMBEDDED MODULE", nomsub);
                  ret_val = -8;
               }
               printf("%-132s\n", messag);
            }
         } else {
            sprintf(messag, "%s: INVALID TYPE IN EMBEDDED MODULE", nomsub);
            printf("%-132s\n", messag);
            ret_val = -6;
         }
         goto L666;
L25:
         continue;
      }
/*    NOW CALL THE SELECTED DECLARATION MODULE */
      if (strcmp(cmodul, "MODULE") == 0) {
/*       *MODULE      * DECLARATION MODULE */
         iretcd = kdrdmd(my_iptrun, npara2, hpara2);
      } else if (strcmp(cmodul, "LINKED_LIST") == 0) {
/*       *LINKED_LIST * DECLARATION MODULE */
         iretcd = kdrdll(my_iptrun, npara2, hpara2);
      } else if (strcmp(cmodul, "XSM_FILE") == 0) {
/*       *XSM_FILE    * DECLARATION MODULE */
         iretcd = kdrdxf(my_iptrun, npara2, hpara2);
      } else if (strcmp(cmodul, "SEQ_BINARY") == 0) {
/*       *SEQ_BINARY  * DECLARATION MODULE */
         iretcd = kdrdsb(my_iptrun, npara2, hpara2);
      } else if (strcmp(cmodul, "SEQ_ASCII") == 0) {
/*       *SEQ_ASCII   * DECLARATION MODULE */
         iretcd = kdrdsa(my_iptrun, npara2, hpara2);
      } else if (strcmp(cmodul, "DIR_ACCESS") == 0) {
/*       *DIR_ACCESS  * DECLARATION MODULE */
         iretcd = kdrdda(my_iptrun, npara2, hpara2);
      } else {
/*       OTHERWISE, DECLARATION MODULE IS NOT AVAILABLE */
         sprintf(messag, "%s: EMBEDDED PARAMETER MODULE *%s* IS NOT AVAILABLE", nomsub, cmodul);
         printf("%-132s\n", messag);
         ret_val = -9;
         goto L666;
      }
      if (iretcd != 0) {
         sprintf(messag, "%s: PROBLEM WITH EMBEDDED MODULE *%s*", nomsub, cmodul);
         printf("%-132s\n", messag);
         ret_val = -666;
         goto L666;
      }

      if (!l2data) redopn_c(jread, jwrite, hwrite, jrecin);
      redget_c(&ityp, &nitma, &flott, text, &dflot);
      lndata = ityp != 10 && (ityp != 3 || strcmp(text, ";") != 0);
      goto L21;

   } else if (minput == 1) {
/*    RETURN OBJECTS CREATED IN THE PROCEDURE BEFORE ENDING */
      if ((*my_param)->nup != 0) {
         char aparam[13];
         for (iloop1 = 0; iloop1 < (*my_param)->nup; ++iloop1) {
            my_node = clepos(my_param, iloop1);
            if ((my_node->type < 0) || (my_node->type >= 10)) continue;
            if (strlen(my_node->name) > 12) {
              sprintf(AbortString,"%s: name EXCEEDING 12 CHARACTERS.",nomsub);
              xabort_c(AbortString);
            }
            strcpy(aparam, my_node->name);
            iparam = my_node->type;
            if (iparam == 3) kparam = my_node->value.mylcm;

            my_node = clenode(my_iptrun, aparam);
            if (my_node == NULL) {
               printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, aparam);
               ret_val = -23;
               goto L666;
            }
            my_node->type = iparam;
            if (iparam == 3) my_node->value.mylcm = kparam;
         }
      }
   } else {
      sprintf(messag, "%s: INVALID VALUE FOR *MINPUT* ARG", nomsub);
      printf("%-132s\n", messag);
      ret_val = -3;
   }
L666:
   return ret_val;
} /* kdrprm */

int_32 kdrdmd(lifo **my_iptrun, int_32 nentry, char (*hentry)[13])
{

/*     GAN-2000 SYSTEM: R.ROY (01/2000), VERSION 2.1 */

/*             *KDRDMD* IS THE MODULE FOR *MODULE      * DECLARATIONS */
/*                      =0 IF NO ERROR */

/*      INPUT: *MY_IPTRUN* IS THE EXEC STRUCTURE POINTER (ALLOCATED) */
/*             *NENTRY* IS THE # OF LINKED LISTS AND FILES USED. */
/*             *HENTRY* NAMES OF EACH OBJECT <- LINKED LIST OR FILE. */
/*                      ( CHARACTER*12 HENTRY(NENTRY) ) */

   char *nomsub = "kdrdmd";
   int_32 ret_val = 0;
   int_32 ityp, nitma, lndata;
   float_32 flott;
   double_64 dflot;
   char text[73], messag[73];
   int_32 iloop1, iparam;

   redget_c(&ityp, &nitma, &flott, text, &dflot);
   lndata = ityp != 10 && (ityp != 3 || strcmp(text, ";") != 0);
   if (lndata) {
      sprintf(messag, "%s: NOT DEVELOPED YET (RR)", nomsub);
      printf("%-132s\n", messag);
      ret_val = -666;
      goto L666;
   }
   for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
      lifo_node *my_node;
      
      my_node = clenode(my_iptrun, hentry[iloop1]);
      if (my_node == NULL) {
         printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
         ret_val = -21;
         goto L666;
      }
      iparam = my_node->type ;
      if (abs(iparam) != 2) goto L8001;
   }

/* CAN WE FOUND *;* AT THE END OF THE SENTENCE ? */
   if (lndata) {
      redget_c(&ityp, &nitma, &flott, text, &dflot);
      if (ityp != 3 || strcmp(text, ";") != 0) goto L8002;
   }
L666:
   return ret_val;

L8001:
   sprintf(messag, "%s: INVALID TYPE (%d) IN *MODULE* DATA.", nomsub, (int)iparam);
   printf("%-132s\n", messag);
   ret_val = 8001;
   goto L666;
L8002:
   sprintf(messag, "%s: INVALID END  IN *MODULE* DATA.", nomsub);
   printf("%-132s\n", messag);
   ret_val = 8002;
   goto L666;
} /* kdrdmd */

int_32 kdrdll(lifo **my_iptrun, int_32 nentry, char (*hentry)[13])
{

/*     GAN-2000 SYSTEM: R.ROY (01/2000), VERSION 2.1 */

/*             *KDRDLL* IS THE MODULE FOR *LINKED_LIST * DECLARATIONS */
/*                      =0 IF NO ERROR */

/*      INPUT: *MY_IPTRUN* IS THE EXEC STRUCTURE POINTER (ALLOCATED) */
/*             *NENTRY* IS THE # OF LINKED LISTS AND FILES USED. */
/*             *HENTRY* NAMES OF EACH OBJECT <- LINKED LIST OR FILE. */
/*                      ( CHARACTER*12 HENTRY(NENTRY) ) */

   char *nomsub = "kdrdll";
   int_32 ret_val = 0;
   int_32 ityp, nitma, lndata;
   float_32 flott;
   double_64 dflot;
   char text[73], messag[73];
   int_32 iloop1, iparam;

   redget_c(&ityp, &nitma, &flott, text, &dflot);
   lndata = ityp != 10 && (ityp != 3 || strcmp(text, ";") != 0);
   if (lndata) {
      sprintf(messag, "%s: NOT DEVELOPED YET (RR)", nomsub);
      printf("%-132s\n", messag);
      ret_val = -666;
      goto L666;
   }
   for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
      lifo_node *my_node;
      
      my_node = clenode(my_iptrun, hentry[iloop1]);
      if (my_node == NULL) {
         printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
         ret_val = -21;
         goto L666;
      }
      iparam = my_node->type ;
      if (abs(iparam) != 3) goto L8001;
   }

/* CAN WE FOUND *;* AT THE END OF THE SENTENCE ? */
   if (lndata) {
      redget_c(&ityp, &nitma, &flott, text, &dflot);
      if (ityp != 3 || strcmp(text, ";") != 0) goto L8002;
   }
L666:
   return ret_val;

L8001:
   sprintf(messag, "%s: INVALID TYPE (%d) IN *LINKED_LIST * DATA.", nomsub, (int)iparam);
   printf("%-132s\n", messag);
   ret_val = 8001;
   goto L666;
L8002:
   sprintf(messag, "%s: INVALID END  IN *LINKED_LIST * DATA.", nomsub);
   printf("%-132s\n", messag);
   ret_val = 8002;
   goto L666;
} /* kdrdll */

int_32 kdrdxf(lifo **my_iptrun, int_32 nentry, char (*hentry)[13])
{

/*     GAN-2000 SYSTEM: R.ROY (01/2000), VERSION 2.1 */

/*             *KDRDXF* IS THE MODULE FOR *XSM_FILE    * DECLARATIONS */
/*                      =0 IF NO ERROR */

/*      INPUT: *MY_IPTRUN* IS THE EXEC STRUCTURE POINTER (ALLOCATED) */
/*             *NENTRY* IS THE # OF LINKED LISTS AND FILES USED. */
/*             *HENTRY* NAMES OF EACH OBJECT <- LINKED LIST OR FILE. */
/*                      ( CHARACTER*12 HENTRY(NENTRY) ) */

   char *nomsub = "kdrdxf";
   int_32 ret_val = 0;
   int_32 ityp, nitma, lndata;
   float_32 flott;
   double_64 dflot;
   char text[73], messag[73];
   int_32 iprint = 0;
   int_32 iloop1, iparam;

   redget_c(&ityp, &nitma, &flott, text, &dflot);
   lndata = ityp != 10 && (ityp != 3 || strcmp(text, ";") != 0);
   if (lndata) {
      if (ityp == 3) {
         if (strcmp(text, "EDIT") == 0) {
            redget_c(&ityp, &iprint, &flott, text, &dflot);
            if (ityp != 1 && iprint < 0) {
               sprintf(messag, "%s: AFTER *EDIT*, PUT A POSITIVE INT", nomsub);
               printf("%-132s\n", messag);
               ret_val = -1;
               goto L666;
            }
            redget_c(&ityp, &nitma, &flott, text, &dflot);
         }
         if (strcmp(text, "FILE") != 0) {
            sprintf(messag, "%s: EXPECTING *FILE* KEYWORD; TEXT=%.12s", nomsub, text);
            printf("%-132s\n", messag);
            ret_val = -2;
            goto L666;
         }
      } else {
         sprintf(messag, "%s: INVALID INPUT", nomsub);
         printf("%-132s\n", messag);
         ret_val = -666;
         goto L666;
      }
   }
   for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
      lifo_node *my_node;
      
      my_node = clenode(my_iptrun, hentry[iloop1]);
      if (my_node == NULL) {
         printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
         ret_val = -21;
         goto L666;
      }
      iparam = my_node->type ;
      if (abs(iparam) != 4) goto L8001;
      if (lndata) {
         FILE *file;
         redget_c(&ityp, &nitma, &flott, text, &dflot);
         if (ityp != 3) {
            sprintf(messag, "%s: INVALID FILE NAME", nomsub);
            printf("%-132s\n", messag);
            ret_val = -666;
            goto L666;
         }
         file = fopen(text, "r");

/*       DEFINE EXISTENCE MODE */
         if (file != NULL) {
            fclose(file);
            if (iprint != 0) printf("OLD/XF: %s\n", text);
            if (iparam < 0) my_node->type = -iparam;
            my_node->access = 2;
         } else {
            if (iprint != 0) printf("NEW/XF: %s\n", text);
            if (iparam > 0) my_node->type = -iparam;
            my_node->access = 0;
         }

/*       REGISTER FILE NAME */
         strcpy(my_node->OSname, text);
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
   sprintf(messag, "%s: INVALID TYPE (%d) IN *XSM_FILE    * DATA.", nomsub, (int)iparam);
   printf("%-132s\n", messag);
   ret_val = 8001;
   goto L666;
L8002:
   sprintf(messag, "%s: INVALID END  IN *XSM_FILE    * DATA.", nomsub);
   printf("%-132s\n", messag);
   ret_val = 8002;
   goto L666;
} /* kdrdxf */

int_32 kdrdsb(lifo **my_iptrun, int_32 nentry, char (*hentry)[13])
{

/*     GAN-2000 SYSTEM: R.ROY (01/2000), VERSION 2.1 */

/*             *KDRDSB* IS THE MODULE FOR *SEQ_BINARY  * DECLARATIONS */
/*                      =0 IF NO ERROR */

/*      INPUT: *MY_IPTRUN* IS THE EXEC STRUCTURE POINTER (ALLOCATED) */
/*             *NENTRY* IS THE # OF LINKED LISTS AND FILES USED. */
/*             *HENTRY* NAMES OF EACH OBJECT <- LINKED LIST OR FILE. */
/*                      ( CHARACTER*12 HENTRY(NENTRY) ) */

   char *nomsub = "kdrdsb";
   int_32 ret_val = 0;
   int_32 ityp, nitma, lndata;
   float_32 flott;
   double_64 dflot;
   char text[73], messag[73];
   int_32 iprint = 0;
   int_32 iloop1, iparam;

   redget_c(&ityp, &nitma, &flott, text, &dflot);
   lndata = ityp != 10 && (ityp != 3 || strcmp(text, ";") != 0);
   if (lndata) {
      if (ityp == 3) {
         if (strcmp(text, "EDIT") == 0) {
            redget_c(&ityp, &iprint, &flott, text, &dflot);
            if (ityp != 1 && iprint < 0) {
               sprintf(messag, "%s: AFTER *EDIT*, PUT A POSITIVE INT", nomsub);
               printf("%-132s\n", messag);
               ret_val = -1;
               goto L666;
            }
            redget_c(&ityp, &nitma, &flott, text, &dflot);
         }
         if (strcmp(text, "FILE") != 0) {
            sprintf(messag, "%s: EXPECTING *FILE* KEYWORD; TEXT=%.12s", nomsub, text);
            printf("%-132s\n", messag);
            ret_val = -2;
            goto L666;
         }
      } else {
         sprintf(messag, "%s: INVALID INPUT", nomsub);
         printf("%-132s\n", messag);
         ret_val = -666;
         goto L666;
      }
   }
   for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
      lifo_node *my_node;
      
      my_node = clenode(my_iptrun, hentry[iloop1]);
      if (my_node == NULL) {
         printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
         ret_val = -21;
         goto L666;
      }
      iparam = my_node->type ;
      if (abs(iparam) != 5) goto L8001;
      if (lndata) {
         FILE *file;
         redget_c(&ityp, &nitma, &flott, text, &dflot);
         if (ityp != 3) {
            sprintf(messag, "%s: INVALID FILE NAME", nomsub);
            printf("%-132s\n", messag);
            ret_val = -666;
            goto L666;
         }
         file = fopen(text, "r");

/*       DEFINE EXISTENCE MODE */
         if (file != NULL) {
            fclose(file);
            if (iprint != 0) printf("OLD/SB: %s\n", text);
            if (iparam < 0) my_node->type = -iparam;
            my_node->access = 2;
         } else {
            if (iprint != 0) printf("NEW/SB: %s\n", text);
            if (iparam > 0) my_node->type = -iparam;
            my_node->access = 0;
         }

/*       REGISTER FILE NAME */
         strcpy(my_node->OSname, text);
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
   sprintf(messag, "%s: INVALID TYPE (%d) IN *SEQ_BINARY  * DATA.", nomsub, (int)iparam);
   printf("%-132s\n", messag);
   ret_val = 8001;
   goto L666;
L8002:
   sprintf(messag, "%s: INVALID END  IN *SEQ_BINARY  * DATA.", nomsub);
   printf("%-132s\n", messag);
   ret_val = 8002;
   goto L666;
} /* kdrdsb */

int_32 kdrdsa(lifo **my_iptrun, int_32 nentry, char (*hentry)[13])
{

/*     GAN-2000 SYSTEM: R.ROY (01/2000), VERSION 2.1 */

/*             *KDRDSA* IS THE MODULE FOR *SEQ_ASCII   * DECLARATIONS */
/*                      =0 IF NO ERROR */

/*      INPUT: *MY_IPTRUN* IS THE EXEC STRUCTURE POINTER (ALLOCATED) */
/*             *NENTRY* IS THE # OF LINKED LISTS AND FILES USED. */
/*             *HENTRY* NAMES OF EACH OBJECT <- LINKED LIST OR FILE. */
/*                      ( CHARACTER*12 HENTRY(NENTRY) ) */

   char *nomsub = "kdrdsa";
   int_32 ret_val = 0;
   int_32 ityp, nitma, lndata;
   float_32 flott;
   double_64 dflot;
   char text[73], messag[73];
   int_32 iprint = 0;
   int_32 iloop1, iparam;

   redget_c(&ityp, &nitma, &flott, text, &dflot);
   lndata = ityp != 10 && (ityp != 3 || strcmp(text, ";") != 0);
   if (lndata) {
      if (ityp == 3) {
         if (strcmp(text, "EDIT") == 0) {
            redget_c(&ityp, &iprint, &flott, text, &dflot);
            if (ityp != 1 && iprint < 0) {
               sprintf(messag, "%s: AFTER *EDIT*, PUT A POSITIVE INT", nomsub);
               printf("%-132s\n", messag);
               ret_val = -1;
               goto L666;
            }
            redget_c(&ityp, &nitma, &flott, text, &dflot);
         }
         if (strcmp(text, "FILE") != 0) {
            sprintf(messag, "%s: EXPECTING *FILE* KEYWORD; TEXT=%.12s", nomsub, text);
            printf("%-132s\n", messag);
            ret_val = -2;
            goto L666;
         }
      } else {
         sprintf(messag, "%s: INVALID INPUT", nomsub);
         printf("%-132s\n", messag);
         ret_val = -666;
         goto L666;
      }
   }
   for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
      lifo_node *my_node;
      
      my_node = clenode(my_iptrun, hentry[iloop1]);
      if (my_node == NULL) {
         printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
         ret_val = -21;
         goto L666;
      }
      iparam = my_node->type ;
      if (abs(iparam) != 6) goto L8001;
      if (lndata) {
         FILE *file;
         redget_c(&ityp, &nitma, &flott,text, &dflot);
         if (ityp != 3) {
            sprintf(messag, "%s: INVALID FILE NAME", nomsub);
            printf("%-132s\n", messag);
            ret_val = -666;
            goto L666;
         }
         file = fopen(text, "r");

/*       DEFINE EXISTENCE MODE */
         if (file != NULL) {
            fclose(file);
            if (iprint != 0) printf("OLD/SA: %s\n", text);
            if (iparam < 0) my_node->type = -iparam;
            my_node->access = 2;
         } else {
            if (iprint != 0) printf("NEW/SA: %s\n", text);
            if (iparam > 0) my_node->type = -iparam;
            my_node->access = 0;
         }

/*       REGISTER FILE NAME */
         strcpy(my_node->OSname, text);
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
   sprintf(messag, "%s: INVALID TYPE (%d) IN *SEQ_ASCII   * DATA.", nomsub, (int)iparam);
   printf("%-132s\n", messag);
   ret_val = 8001;
   goto L666;
L8002:
   sprintf(messag, "%s: INVALID END  IN *SEQ_ASCII   * DATA.", nomsub);
   printf("%-132s\n", messag);
   ret_val = 8002;
   goto L666;
} /* kdrdsa */

int_32 kdrdda(lifo **my_iptrun, int_32 nentry, char (*hentry)[13])
{

/*     GAN-2000 SYSTEM: R.ROY (01/2000), VERSION 2.1 */

/*             *KDRDDA* IS THE MODULE FOR *DIR_ACCESS  * DECLARATIONS */
/*                      =0 IF NO ERROR */

/*      INPUT: *MY_IPTRUN* IS THE EXEC STRUCTURE POINTER (ALLOCATED) */
/*             *NENTRY* IS THE # OF LINKED LISTS AND FILES USED. */
/*             *HENTRY* NAMES OF EACH OBJECT <- LINKED LIST OR FILE. */
/*                      ( CHARACTER*12 HENTRY(NENTRY) ) */

   char *nomsub = "kdrdda";
   int_32 ret_val = 0;
   int_32 ityp, nitma, lnfile=0, lndata;
   float_32 flott;
   double_64 dflot;
   char text[73], messag[73], filenm[73];
   int_32 iprint = 0;
   int_32 iloop1, iparam, lparam;

   redget_c(&ityp, &nitma, &flott, text, &dflot);
   lndata = ityp != 10 && (ityp != 3 || strcmp(text, ";") != 0);
   if (lndata) {
      if (ityp == 3) {
         if (strcmp(text, "EDIT") == 0) {
            redget_c(&ityp, &iprint, &flott, text, &dflot);
            if (ityp != 1 && iprint < 0) {
               sprintf(messag, "%s: AFTER *EDIT*, PUT A POSITIVE INT", nomsub);
               printf("%-132s\n", messag);
               ret_val = -1;
               goto L666;
            }
            redget_c(&ityp, &nitma, &flott, text, &dflot);
         }
         if (strcmp(text, "FILE") == 0) {
            lnfile = 1;
            redget_c(&ityp, &nitma, &flott, text, &dflot);
         } else {
            lnfile = 0;
         }
         if (strcmp(text, "RECL") != 0) {
            sprintf(messag, "%s: KEYWORD *RECL* MUST BE THERE", nomsub);
            printf("%-132s\n", messag);
            ret_val = -2;
            goto L666;
         }
      } else {
         sprintf(messag, "%s: INVALID INPUT", nomsub);
         printf("%-132s\n", messag);
         ret_val = -3;
         goto L666;
      }
   }
   for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
      lifo_node *my_node;
      
      my_node = clenode(my_iptrun, hentry[iloop1]);
      if (my_node == NULL) {
         printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
         ret_val = -21;
         goto L666;
      }
      iparam = my_node->type ;
      if (abs(iparam) != 7) goto L8001;
      if (lndata) {
         lparam = 0;
         strcpy(filenm, ";");
         redget_c(&ityp, &nitma, &flott, text, &dflot);
         if (ityp == 1) lparam = nitma;
         if (ityp == 3) strcpy(filenm, text);
         if (lnfile) {
            redget_c(&ityp, &nitma, &flott, text, &dflot);
            if (ityp == 1) lparam = nitma;
            if (ityp == 3) strcpy(filenm, text);
         }
         if (lparam <= 0) {
            sprintf(messag, "%s: INVALID VALUE FOR *RECL*", nomsub);
            printf("%-132s\n", messag);
            ret_val = -4;
            goto L666;
         } else {
/*          REGISTER RECORD LENGTH OF EACH DA FILES */
            my_node->lparam = lparam;
         }
         if (lnfile) {
            FILE *file;
            if (strcmp(filenm, ";") == 0) {
               sprintf(messag, "%s: INVALID FILE NAME", nomsub);
               printf("%-132s\n", messag);
               ret_val = -5;
               goto L666;
            }
            file = fopen(filenm, "r");

/*          DEFINE EXISTENCE MODE */
            if (file != NULL) {
               fclose(file);
               if (iprint != 0) printf("OLD/DA: %s\n", text);
               if (iparam < 0) my_node->type = -iparam;
               my_node->access = 2;
            } else {
               if (iprint != 0) printf("NEW/DA: %s\n", text);
               if (iparam > 0) my_node->type = -iparam;
               my_node->access = 0;
            }

/*          REGISTER FILE NAME */
            strcpy(my_node->OSname, text);
         }
      } else {
/*       ONLY VERIFY IF *RECL* MAKES SENSE FOR A DEFINED DA */
         if (iparam > 0) {
            lparam = my_node->lparam;
            if (lparam <= 0) {
               sprintf(messag, "%s: INVALID VALUE FOR *RECL*", nomsub);
               printf("%-132s\n", messag);
               ret_val = -6;
               goto L666;
            }
         }
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
   sprintf(messag, "%s: INVALID TYPE (%d) IN *DIR_ACCESS  * DATA.", nomsub, (int)iparam);
   printf("%-132s\n", messag);
   ret_val = 8001;
   goto L666;
L8002:
   sprintf(messag, "%s: INVALID END  IN *DIR_ACCESS  * DATA.", nomsub);
   printf("%-132s\n", messag);
   ret_val = 8002;
   goto L666;
} /* kdrdda */
