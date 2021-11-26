
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 12/05/09      */
/*****************************************/

#include <stdlib.h>
#include <string.h>
#include "cle2000.h"
#include "header.h"
#define index_f(A, B) (strstr(A, B) == NULL ? 0 : strstr(A, B) - A + 1)
#define maxdxt 200 /* maximum number of modules */
#define ndclkw 8
#define nmodst 15
#define nmawrd 36

int_32 objstk(kdi_file *iunito, FILE *iwrite, int_32 ldatav)
{
   char *nomsub="objstk";
   static char cerror[] = "* GAN-2000 VERS 2.1 * ERROR FOUND FOR THIS LINE";
   static char cl2000[] = "CLE2000(V3)";
   static char alphab[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz";
   static char digits[] = "0123456789";
   static char *ckeywd[] = {"INTEGER", "REAL", "STRING", "DOUBLE", "LOGICAL", "EVALUATE", "ECHO", "ELSEIF",
                            "IF", "WHILE", "UNTIL", "ENDWHILE", "REPEAT", "ELSE", "ENDIF", "THEN", "DO",
                            "QUIT", "NOT", "ABS", "CHS", "LN", "SIN", "COS", "TAN", "ARCSIN", "ARCCOS",
                            "ARCTAN", "EXP", "SQRT", "R_TO_I", "D_TO_I", "I_TO_R", "D_TO_R", "I_TO_D",
                            "R_TO_D", "I_TO_S", "I_TO_S4", "_MIN_", "_MAX_", "_TRIM_"};
   static char *cdclkw[] = {"PROCEDURE", "MODULE", "LINKED_LIST", "XSM_FILE", "SEQ_BINARY", "SEQ_ASCII",
                            "DIR_ACCESS", "PARAMETER"};

/*     GAN-2000 SYSTEM: R.ROY (12/1999), VERSION 2.0 */

/*             *OBJSTK* FIRST-PASS COMPILE OF THE D.A. UNIT   *IUNITO* */
/*                                 NOW INCLUDING OBJECTS & MODULES */
/*                      RESULT IS STILL THE OBJECT D.A. UNIT  *IUNITO* */
/*                      COMPILER COMMENTS ARE WRITTEN ON UNIT *IWRITE* */
/*                      STACK IS BUILT AT THE END OF          *IUNITO* */

/*        USE:          MODULE+OBJECT NAMES ARE DEFINED AND ALLOCATED, */
/*                      CONSISTENCE OF OBJECTS IN CALL */
/*                      STATEMENTS ARE ALSO CHECKED. */

/*      INPUT: *IUNITO* IS THE DIRECT ACCESS UNIT FOR OBJECT CODE */
/*             *IWRITE* IS THE OUTPUT UNIT FOR COMPILER COMMENTS */
/*             *LDATAV* =0/1: PROCEDURE SECTION/DATA SECTION */

/*       NOTE: *OBJSTK* = 0 IF NO PROBLEM WAS ENCOUNTERED WHILE COMPILING */

   int_32 ret_val = 0;
   char cmodul[13], cparin[13], cparav[13], myreco[73], cdatin[73];
   int_32 i, iretcd, nrecor, ninput, nstack, idblst, idatin, iofset, iloop1, jloop2, ilines,
          ilevel, indlin, idclin, idefin, iusein;
   float_32 adatin;
   double_64 ddatin;
   int_32 maskck[nmaskc], ipacki[nmaskc];
   int_32 idebwd[nmawrd+1], ifinwd[nmawrd+1], jndlec[nmawrd];
   int_32 ivabeg, ivaend, ilogin, lequal, nobjet, ndatav, imodul, logprv, nembed, nmodul;
   int_32 irecor, krecor=0;

/* READ TOP OF OBJECT FILE */
   iretcd = kdiget_c(iunito, (int_32 *)&header, 0, kdisize(header));
   if (iretcd != 0) goto L9023;
   strcpy(cparin, header.cparin);
   strcpy(myreco, header.cdatin);
   nrecor = header.nrecor;
   ninput = header.ninput;
   nstack = header.nstack;
   idblst = header.idblst;
   if (strcmp(cparin, cl2000) != 0) goto L9025;
   if (idblst > 0) {
      printf("%-72s   LINE\n", cerror);
      printf(" \n");
   }
   ivabeg = ninput + nstack;
   ivaend = ninput + nstack;
   ilogin = 0;
   lequal = 0;
   nobjet = -1;
   ndatav = -1;
   imodul = 0;
   logprv = 1;
   nembed = 0;
   nmodul = 0;

/* ***  MAIN LOOP OVER RECORDS (BEGIN) */
   for (irecor = 2; irecor <= ninput; ++irecor) {
      int_32 iwords = 1;
      int_32 nwords = 1;
      int_32 jbiprv = 0;

/*    READ A NEW RECORD */
      iofset = (irecor - 1) * lrclen;
      iretcd = kdiget_c(iunito, (int_32 *)&record1, iofset, kdisize(record1));
      if (iretcd != 0) goto L9023;
      strcpy(cparin, record1.cparin);
      strcpy(myreco, record1.myreco);
      ilines = record1.ilines;
      ilevel = record1.ilevel;
      for (i = 0; i < nmaskc; i++) maskck[i] = record1.maskck[i];
      for (i = 0; i < nmaskc; i++) ipacki[i] = record1.ipacki[i];

/*    RECORDS  INSIDE CLE-2000, NOTHING TO DO */
      if (ilevel != 0) goto L100;

/*    RECORDS OUTSIDE CLE-2000, SCRUTINIZE... */

/*    BEGIN: MASK RECOVERY */
      for (iloop1 = 1; iloop1 <= 72; ++iloop1) {
         int_32 jbicur;
         jloop2 = (iloop1 + 23) / 24;
         jbicur = maskck[jloop2 - 1] % 2;
         iwords += jbiprv * (1 - jbicur);
         idebwd[nwords - 1] = iloop1;
         ifinwd[iwords - 1] = iloop1;
         nwords += jbicur * (1 - jbiprv);
         jbiprv = jbicur;
         maskck[jloop2 - 1] /= 2;
      }
      --nwords;
/*    END:   MASK RECOVERY */

/*    BEGIN: UNPACK JNDLEC WITH TYPES (ITYP-1) */
      for (iloop1 = 1; iloop1 <= nwords; ++iloop1) {
         jloop2 = ((iloop1 << 1) + 23) / 24;
         jndlec[iloop1 - 1] = ipacki[jloop2 - 1] % 4;
         ipacki[jloop2 - 1] /= 4;
      }
/*    END:   UNPACK JNDLEC WITH TYPES (ITYP-1) */

/*    CHECK ALL DECLARATION STATEMENTS, IF NOT YET FOUND */
      if (logprv) {
         krecor = irecor;
         if (jndlec[0] == 2 && myreco[idebwd[0]-1] != '\'' && ifinwd[0]-idebwd[0] <= 11) {
            strncpy(cparav, &myreco[idebwd[0]-1], ifinwd[0]-idebwd[0]+1);
            cparav[ifinwd[0]-idebwd[0]+1] = '\0';
            for (iloop1 = 1; iloop1 <= ndclkw; ++iloop1) {
               if (strcmp(cparav, cdclkw[iloop1-1]) == 0) ilogin = iloop1;
            }
            if (ilogin != 0) {
               strcpy(cmodul, cdclkw[ilogin-1]);
               imodul = -ilogin;
               ++nmodul;
            }
         }
      }
      logprv = 0;

/*    SCAN OTHER WORDS */
      for (iloop1 = 1; iloop1 <= nwords; ++iloop1) {
         int_32 ileng = ifinwd[iloop1-1] - idebwd[iloop1-1] + 1;

/*       ARE WE IN THE DATA SECTION ? */
         if (ldatav) {
            ++ndatav;
/*          INSIDE  THE DATA SECTION ( ...   :: *HERE* ) */
/*          THEN COUNT EMBEDDED MODULES UNTIL MODULE ENDING */
            if (jndlec[iloop1-1] == 2 && myreco[idebwd[iloop1-1]-1] != '\''
                && ifinwd[iloop1-1]-idebwd[iloop1-1] <= 2) {
               strncpy(cparav, &myreco[idebwd[iloop1-1]-1], ileng);
               cparav[ileng] = '\0';
               if (strcmp(cparav, ":::") == 0) {
                  ++nembed;
               } else if (strcmp(cparav, ";") == 0) {
                  if (iloop1 != nwords) goto L9010;
                  if (nembed == 0) {
/*                          END OF STATEMENT REACHED */
                     logprv = 1;
                  } else {
                     --nembed;
                  }
               } else if (strcmp(cparav, "::") == 0) {
                  goto L5002;
               }
            }
         } else {
            char clisto[17];
/*          OUTSIDE THE DATA SECTION ( *HERE*  ::  ... ) */
/*          NOTE:   EVERY OBJECT/MODULE MUST BE FIRST DECLARED */
            if (jndlec[iloop1-1] != 2 || myreco[idebwd[iloop1-1]-1] == '\''
                || ifinwd[iloop1-1]-idebwd[iloop1 - 1] > 15) goto L5001;
            strncpy(clisto, &myreco[idebwd[iloop1-1]-1], ileng);
            clisto[ileng] = '\0';
            if (ifinwd[iloop1-1]-idebwd[iloop1-1] == 1 && strcmp(clisto, ":=") == 0) {
               lequal = 1;
            } else {
/*             REMAININGS: 1 MODULE & OBJECTS ... */
               if (ifinwd[iloop1-1]-idebwd[iloop1-1] > 11) goto L5001;
               strncpy(cparav, &myreco[idebwd[iloop1-1]-1], ileng);
               cparav[ileng] = '\0';
               if (strcmp(cparav, ";") == 0) {
                  if (iloop1 != nwords) goto L9010;
/*                END OF STATEMENT REACHED " "*/
                  logprv = 1;
               } else if (strcmp(cparav, "::") == 0) {
                  ldatav = 1;
               } else if (strcmp(cparav, ":::") == 0) {
                  goto L5002;
               } else {
/*                USING *CPARAV* VARIABLE, SCAN ALL DECLARED OBJECTS */
                  int_32 ilowrc = ivabeg;
                  int_32 ihigrc = ivaend + 1;
L41:
                  if (ihigrc - ilowrc <= 1) {
                     char cc[2];
/*                   OBJECT/MODULE NOT FOUND */
                     if (ilogin == 0) goto L5004;
                     ++ivaend;

/*                   SHIFT GREATER VALUES */
                     if (ihigrc != ivaend) {
                        for (jloop2 = ivaend - 1; jloop2 >= ihigrc; --jloop2) {
                           iofset = (jloop2 - 1) * lrclen;
                           iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                           if (iretcd != 0) goto L9003;
                           iofset = jloop2 * lrclen;
                           iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                           if (iretcd != 0) goto L9001;
                        }
                     }

/*                   CHECK IF OBJECT/MODULE NAME COMPLIES WITH THE RULES */
                     strncpy(&cc[0], &cparav[0], 1); cc[1] = '\0';
                     if (index_f(alphab, cc) == 0) {
                        printf("%s: CHARACTER *%s* IS NOT ALLOWED\n",nomsub, cc);
                        goto L5001;
                     }
                     for (jloop2 = 2; jloop2 <= strlen(cparav); ++jloop2) {
                        int_32 jin1, jin2, jin3, jin4, jin5;
                        strncpy(&cc[0], &cparav[jloop2-1], 1); cc[1] = '\0';
                        jin1 = index_f(alphab, cc);
                        jin2 = index_f(digits, cc);
                        jin3 = index_f(cc, " ");
                        jin4 = index_f(cc, ".");
                        jin5 = index_f(cc, ":");
                        if (jin1 + jin2 + jin3 + jin4 + jin5 == 0) {
                           printf("%s: CHARACTER *%1s* IN *%s* IS NOT ALLOWED\n",nomsub, cc, cparav);
                           goto L5001;
                        }
                     }

/*                   CHECK IF OBJECT/MODULE NAME IS A KEYWORD */
                     for (jloop2 = 1; jloop2 <= 40; ++jloop2) {
                        if (strcmp(cparav, ckeywd[jloop2-1]) == 0) {
                           printf("%s: OBJECT *%s* IS A CLE-2000 KEYWORD\n", nomsub, cparav);
                           goto L5001;
                        }
                     }

/*                   VALID OBJECT/MODULE NAME, WRITE AT END */
                     for (i = 0; i < 72; i++) cdatin[i] = ' ';
                     cdatin[72] = '\0';
                     indlin = -ilogin;
                     idatin = 0;
                     adatin = 0.f;
                     ddatin = 0.;
                     idclin = ilines;
                     idefin = 0;
                     iusein = 0;
                     for (jloop2 = 0; jloop2 < ndclkw; ++jloop2) {
/*                              TO ACCEPT DECLARATIONS AS DEFINED MODULES */
                        if (strcmp(cparav, cdclkw[jloop2]) == 0) indlin = 2;
                     }

/*                   VALID OBJECT/MODULE NAME, WRITE AT *IHIGRC* */
                     iofset = (ihigrc - 1) * lrclen;
                     strcpy(record2.cparin, cparav);
                     strcpy(record2.cdatin, cdatin);
                     record2.indlin = indlin;
                     record2.idatin = idatin;
                     record2.adatin = adatin;
                     record2.ddatin = ddatin;
                     record2.idclin = idclin;
                     record2.idefin = idefin;
                     record2.iusein = iusein;
                     iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                     if (iretcd != 0) goto L9001;
                     ++nobjet;
                  } else {
                     int_32 imedrc = (ihigrc + ilowrc) / 2;
                     iofset = (imedrc - 1) * lrclen;
                     iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                     if (iretcd != 0) goto L9003;
                     strcpy(cparin, record2.cparin);
                     strcpy(cdatin, record2.cdatin);
                     indlin = record2.indlin;
                     idatin = record2.idatin;
                     adatin = record2.adatin;
                     ddatin = record2.ddatin;
                     idclin = record2.idclin;
                     idefin = record2.idefin;
                     iusein = record2.iusein;
                     if (strcmp(cparin, cparav) == 0) {
/*                      OBJECT/MODULE FOUND */
                        if (ilogin != 0) {
/*                         DECLARATION STATEMENT, VERIFY CONSISTENCY */
/*                         IF IT IS THE DECLARATION MODULE NAME */
                           if (strcmp(cparav, cdclkw[ilogin-1]) == 0) {
                              if (indlin != 2) goto L8000;
                           } else {
                              if (abs(indlin) != ilogin) goto L8000;
                           }
                        } else {
                           if (abs(indlin) == 1 || abs(indlin) == 2) {

/*                            PROC/MODULE NAME IS FOUND */
                              strcpy(cmodul, cparav);
                              imodul = abs(indlin);
                              ++nmodul;
                              lequal = 1;
                           } else {
                              if (lequal && iusein == 0) {
/*                               CHANGE FIRST USED LINE FOR OBJECT */
                                 record2.iusein = ilines;
                                 iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                                 if (iretcd != 0) goto L9001;
                              } else if (!lequal && idefin == 0) {
/*                                          CHANGE FIRST DEFINED LINE FOR OBJECT */
                                 record2.idefin = ilines;
                                 iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                                 if (iretcd != 0) goto L9001;
                              }
                           }
                           ++nobjet;
                        }
                     } else if (strcmp(cparin, cparav) < 0) {
                        ilowrc = imedrc;
                        goto L41;
                     } else {
                        ihigrc = imedrc;
                        goto L41;
                     }
                  }
               }
            }
         }
      }

/*    STATEMENT END WAS REACHED, */
/*    WRITE MODULE NAME IN 1ST RECORD OF THIS STATEMENT */
      if (logprv) {
         char ctcall[9], ctciox[19], ctcobj[19], ctotcl[73];
         if (nmodul == 0) {
/*          NO MODULE FOUND, IMPOSE */
/*                 => CMODUL = 'IOX:' (WHEN NO OBJECTS) */
/*                 => CMODUL = 'EQU:' (OTHERWISE) */
            if (nobjet == -1) {
               strcpy(cmodul, "IOX:");
            } else {
               strcpy(cmodul, "EQU:");
            }
            nmodul = 1;
            imodul = 2;
         } else if (nmodul != 1) {
            goto L5008;
         }

         iofset = (krecor - 1) * lrclen;
         iretcd = kdiget_c(iunito, (int_32 *)&record1, iofset, kdisize(record1));
         if (iretcd != 0) goto L9023;
         strcpy(cparin, record1.cparin);
         strcpy(myreco, record1.myreco);
         ilines = record1.ilines;
         ilevel = record1.ilevel;
         for (i = 0; i < nmaskc; i++) maskck[i] = record1.maskck[i];
         for (i = 0; i < nmaskc; i++) ipacki[i] = record1.ipacki[i];

/*       WITH MODULE NAME, ADD THE NUMBER OF DATA ITEMS (EXCEPT *;*) */
         strcpy(record1.cparin, cmodul);
         record1.irecor = ndatav;
         iretcd = kdiput_c(iunito, (int_32 *)&record1, iofset, kdisize(record1));
         if (iretcd != 0) goto L9023;
         if (idblst > 0) {
            if (imodul > 0) {
               strcpy(ctcall, "CALL");
            } else {
               strcpy(ctcall, "DECLARE ");
               strcpy(cmodul, " ");
            }
            if (nobjet == -1) {
               strcpy(ctcobj, " WITHOUT  OBJ    ,");
            } else {
               sprintf(ctcobj, " WITH %d OBJ VAL", (int)nobjet);
            }
            if (ndatav == -1) {
               strcpy(ctciox, " WITHOUT  I/O");
            } else {
               sprintf(ctciox, " WITH %d I/O VAL", (int)nobjet);
            }
            sprintf(ctotcl, "%8s%12s%12s%18s%18s", ctcall,cdclkw[abs(imodul)-1],cmodul,ctcobj,ctciox);
            fprintf(iwrite, "%-72s %7d\n", ctotcl, (int)ilines);
         }

/*       RESET THINGS BEFORE NEXT MODULE... */
         nmodul = 0;
         lequal = 0;
         ldatav = 0;
         nobjet = -1;
         ndatav = -1;
         imodul = 0;
         ilogin = 0;
      }

L100:
      ;
   }
/* ***  MAIN LOOP OVER RECORDS (END) */
/* ALL VARIABLES ARE NOW SORTED AT THE END OF THE OBJECT FILE */
/* REWRITE TOP OF OBJECT FILE TO UPDATE *NSTACK+/+NRECOR* */
   nobjet = ivaend - ivabeg;
   nrecor = ivaend;
   header.nrecor = nrecor;
   header.nobjet = nobjet;
   iretcd = kdiput_c(iunito, (int_32 *)&header, 0, kdisize(header));
   if (iretcd != 0) goto L9001;
   if (idblst > 0) fprintf(iwrite, " \n");

L666:
   return ret_val;

L5000:
   printf("%-72s   LINE\n", cerror);
   printf("%-72s   %04d\n", myreco, (int)ilines);
   goto L666;
L5001:
   printf("! %s: INVALID OBJECT/MODULE NAME IN RECORD\n", nomsub);
   ret_val = 5001;
   goto L5000;
L5002:
   printf("! %s: INVALID EMBEDDED MODULES, REVIEW SYNTAX\n", nomsub);
   ret_val = 5002;
   goto L5000;
L5004:
   printf("! %s: OBJECT/MODULE NOT YET DECLARED *%s*\n", nomsub, cparav);
   ret_val = 5004;
   goto L5000;
L5008:
   printf("! %s: MORE THAN 1 MODULE FOUND\n", nomsub);
   ret_val = 5008;
   goto L5000;
L8000:
   printf("! %s: *%s* NOW WITH TYPE %s\n", nomsub, cparav, cdclkw[ilogin-1]);
   ret_val = 8000;
   goto L5000;
L9001:
   iretcd = -1;
   printf("! %s: WRITING RETURN CODE =%d\n", nomsub, (int)iretcd);
   ret_val = iretcd;
   goto L666;
L9003:
   iretcd = -1;
   printf("! %s: READING RETURN CODE =%d\n", nomsub, (int)iretcd);
   ret_val = iretcd;
   goto L666;
L9010:
   printf("! %s: UNEXPECTED END OF STATEMENT\n", nomsub);
   ret_val = 9010;
   goto L5000;
L9023:
   iretcd = -1;
   printf("! %s: IOSTAT RETURN CODE =%d\n", nomsub, (int)iretcd);
   printf("! %s: IMPOSSIBLE TO USE THIS *OBJECT* FILE\n", nomsub);
   ret_val = -2;
   goto L666;
L9025:
   printf("! %s: IMPOSSIBLE TO USE OLD  *OBJECT* FILE\n", nomsub);
   ret_val = -3;
   goto L666;
} /* objstk */
