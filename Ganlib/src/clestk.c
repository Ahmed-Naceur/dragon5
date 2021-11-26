
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 11/05/09      */
/*****************************************/

#include <stdlib.h>
#include <string.h>
#include "cle2000.h"
#include "header.h"
#define index_f(A, B) (strstr(A, B) == NULL ? 0 : strstr(A, B) - A + 1)
#define ndimst 128
#define nmawrd 36

int_32 clestk(kdi_file *iunito, FILE *iwrite,
              int_32 (*dumcst)(char *, int_32 *, int_32 *, float_32 *, char *, double_64 *))
{
   char *nomsub = "clestk";
   static char cerror[] = "* CLE-2000 VERS 3.0 * ERROR FOUND FOR THIS LINE *";
   static char cl2000[] = "CLE2000(V3)";
   static char cequal[] = ":=";
   static char alphab[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz";
   static char digits[] = "0123456789";
   static char *ckeywd[] = {"INTEGER", "REAL", "STRING", "DOUBLE", "LOGICAL", "EVALUATE", "ECHO", "ELSEIF",
                            "IF", "WHILE", "UNTIL", "ENDWHILE", "REPEAT", "ELSE", "ENDIF", "THEN", "DO",
                            "QUIT", "NOT", "ABS", "CHS", "LN", "SIN", "COS", "TAN", "ARCSIN", "ARCCOS",
                            "ARCTAN", "EXP", "SQRT", "R_TO_I", "D_TO_I", "I_TO_R", "D_TO_R", "I_TO_D",
                            "R_TO_D", "I_TO_S", "I_TO_S4", "_MIN_", "_MAX_", "_TRIM_"};
   static char *clognd[] = {";", ";", ";", ";", ";", ";", ";", "THEN", "THEN", "DO",
                            ";", ";", "REPEAT", "ELSE", ";"};

/*     CLE-2000 SYSTEM: R.ROY (11/1999), VERSION 3.0 */

/*             *CLESTK* SECOND-PASS COMPILE OF THE D.A. UNIT  *IUNITO* */
/*                      RESULT IS STILL THE OBJECT D.A. UNIT  *IUNITO* */
/*                      COMPILER COMMENTS ARE WRITTEN ON UNIT *IWRITE* */
/*                      STACK IS BUILT AT THE END OF          *IUNITO* */

/*        USE:          VARIABLE NAMES ARE DEFINED AND ALLOCATED, */
/*                      CONSISTENCE OF TYPES IN EVALUATIONS IS CHECKED, */
/*                      <<.>> AND >>.<< STATEMENTS ARE ALSO CHECKED. */

/*      INPUT: *IUNITO* IS THE DIRECT ACCESS UNIT FOR OBJECT CODE */
/*             *IWRITE* IS THE OUTPUT UNIT FOR COMPILER COMMENTS */
/*             *DUMCST* IS THE EXTERNAL FUNCTION FOR *CLE-2000* CONSTANTS */

/*       NOTE: *CLESTK* = 0 IF NO PROBLEM WAS ENCOUNTERED WHILE COMPILING */

   int_32 ret_val = 0;
   char cparin[13], myreco[73], cdatin[73], cparav[13];
   int_32 indrgt[ndimst], indlft[ndimst], irclft[ndimst];
   int_32 idebwd[nmawrd+1], ifinwd[nmawrd+1], jndlec[nmawrd];
   float_32 adatin;
   double_64 ddatin;
   int_32 maskck[nmaskc], ipacki[nmaskc];
   int_32 i, iretcd, nrecor, ninput, maxlvl, nstack, iofset, ilines, ilevel,
          indlin, idclin, idefin, iusein, idatin, iloop1, jloop2;
   int_32 ivabeg, ivaend, ilogin, nstlft, nstrgt;
   int_32 irecor;
   int_32 lequal=0;

/* READ TOP OF OBJECT FILE */
   iretcd = kdiget_c(iunito, (int_32 *)&header, 0, kdisize(header));
   if (iretcd != 0) goto L9023;
   strcpy(cparin, header.cparin);
   strcpy(myreco, header.cdatin);
   nrecor = header.nrecor;
   ninput = header.ninput;
   maxlvl = header.maxlvl;
   nstack = header.nstack;
   if (strcmp(cparin, cl2000) != 0) goto L9025;
   if (nstack != 0) goto L9025;

/* CASE WHERE THERE ARE NO CLE-2000 SENTENCES */
   if (maxlvl == 0) goto L666;
   ivabeg = ninput;
   ivaend = ninput;
   ilogin = 0;
   nstlft = 0;
   nstrgt = 0;

/* ***  MAIN LOOP OVER RECORDS (BEGIN) */
   for (irecor = 2; irecor <= ninput; ++irecor) {

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

/*    RECORDS OUTSIDE CLE-2000, CHECK >>.<< DEFINITIONS */
/*    AND <<.>> USES (FOR DECLARED VARIABLE) */
      if (ilevel == 0) {
L10:
         iloop1 = index_f(myreco, ">>");
         jloop2 = index_f(myreco, "<<");
         if (iloop1 != 0 || jloop2 != 0) {
            int_32 ilowrc = ivabeg;
            int_32 ihigrc = ivaend + 1;

            if (iloop1 == 0) goto L5006;
            if (jloop2 == 0) goto L5006;

/*          RECOVER VARIABLE NAME INSIDE <<.>> OR >>.<< */
            if (jloop2 > iloop1) {
               strncpy(cparav, &myreco[iloop1+1], jloop2-iloop1-2); cparav[jloop2-iloop1-2] = '\0';
               for (i = iloop1-1; i < jloop2+1; i++) myreco[i] = ' ';
            } else {
               strncpy(cparav, &myreco[jloop2+1], iloop1-jloop2-2); cparav[iloop1-jloop2-2] = '\0';
               for (i = jloop2-1; i < iloop1+1; i++) myreco[i] = ' ';
            }
L11:
            if (ihigrc - ilowrc <= 1) {
/*             VARIABLE WAS NOT FOUND (USED BEFORE DECLARED...) */
               goto L5004;
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
                  if (jloop2 > iloop1 && idefin == 0) {
/*                   CHANGE FIRST DEFINED LINE FOR >>.<< */
                     record2.idefin = ilines;
                     iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                     if (iretcd != 0) goto L9001;
                  } else if (iloop1 > jloop2 && iusein == 0) {
/*                   CHANGE FIRST USED    LINE FOR <<.>> */
                     record2.iusein = ilines;
                     iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                     if (iretcd != 0) goto L9001;
                  }
               } else if (strcmp(cparin, cparav) < 0) {
                  ilowrc = imedrc;
                  goto L11;
               } else {
                  ihigrc = imedrc;
                  goto L11;
               }
            }
            goto L10;
         }

/*    RECORDS INSIDE CLE-2000, CHECK ALL DECLARATIONS AND DEFINITIONS */
      } else {
         char chrend[13];
         int_32 logprv = (ilogin == 0);
         int_32 iwords = 1;
         int_32 nwords = 1;
         int_32 jbiprv = 0;
         if (logprv) {
            for (iloop1 = 1; iloop1 <= 11; ++iloop1) {
               if (strcmp(cparin, ckeywd[iloop1-1]) == 0) ilogin = iloop1;
            }
            if (ilogin == 0) goto L100;

/*          KEYWORDS: *ECHO+/+ELSEIF+/+IF+/+WHILE+/+UNTIL* */
/*          ASSUME THAT THERE WAS AN *:=* SIGN */
/*          (AS IF WE WERE THEN ON RIGHT SIDE OF AN EVALUATE) */
            lequal = ilogin > 6;
            if (lequal) {
               nstlft = 1;
               indlft[0] = 5;
            }
            strcpy(chrend, clognd[ilogin-1]);
         }

/*       HERE, WE HAVE FOUND A SENTENCE INCLUDING A STACK... */

/*       BEGIN: MASK RECOVERY */
         for (iloop1 = 1; iloop1 <= 72; ++iloop1) {
            int_32 jbicur;
            jloop2 = (iloop1 + 23) / 24;
            jbicur = maskck[jloop2 - 1] % 2;
            iwords += jbiprv * (1 - jbicur);
            idebwd[nwords-1] = iloop1;
            ifinwd[iwords-1] = iloop1;
            nwords += jbicur * (1 - jbiprv);
            jbiprv = jbicur;
            maskck[jloop2 - 1] /= 2;
         }
         --nwords;
/*       END:   MASK RECOVERY */

         if (logprv) {
/*          THIS IS A NEW *ILOGIN* */
            if (nwords == 1) goto L100;

/*          START AT CURRENT WORD NUMBER 2 */
            iwords = 2;
         } else {
/*          THIS IS NOW THE FIRST WORD, BUT WITH AN OLD *ILOGIN* */
            iwords = 1;
         }

/*       BEGIN: UNPACK JNDLEC WITH TYPES (ITYP-1) */
         for (iloop1 = 1; iloop1 <= nwords; ++iloop1) {
            jloop2 = ((iloop1 << 1) + 23) / 24;
            jndlec[iloop1-1] = ipacki[jloop2-1] % 4;
            ipacki[jloop2-1] /= 4;
         }
/*       END:   UNPACK JNDLEC WITH TYPES (ITYP-1) */

         for (iloop1 = iwords; iloop1 <= nwords; ++iloop1) {
            if (jndlec[iloop1-1] == 2 && myreco[idebwd[iloop1-1]-1] != '\"') {
               strncpy(cparav, &myreco[idebwd[iloop1-1]-1], ifinwd[iloop1-1]-idebwd[iloop1-1]+1);
               cparav[ifinwd[iloop1-1]-idebwd[iloop1-1]+1] = '\0';
               if (strcmp(cparav, chrend) == 0) {
                  if (iloop1 != nwords) goto L9010;
/*                END OF STATEMENT REACHED */

/*                CHECK CONSISTENCY ON NUMBER OF LEFT/RIGHT EVALUATIONS */
                  if (ilogin != 7 && lequal) {
                     if (nstlft != nstrgt) goto L9008;
                     for (jloop2 = 1; jloop2 <= nstlft; ++jloop2) {
                        if (indrgt[jloop2-1] != indlft[jloop2-1]) goto L9009;
                     }
                  }

/*                RESET *ILOGIN* AND LEFT/RIGHT NUMBERS */
                  ilogin = 0;
                  nstlft = 0;
                  nstrgt = 0;
               } else if (lequal) {
/*                RIGHT-SIDE: VALUES, OPERATIONS & VARIABLES */
/*                VALUES ARE NOW EXCLUDED, STILL OPERATIONS & VARIABLES */
/*                CHECK CONVERSION OPERATIONS */
                  if (strcmp(cparav, "R_TO_I") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 2) goto L8001;
                     indrgt[nstrgt-1] = 1;
                  } else if (strcmp(cparav, "D_TO_I") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 4) goto L8001;
                     indrgt[nstrgt-1] = 1;
                  } else if (strcmp(cparav, "I_TO_R") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 1) goto L8001;
                     indrgt[nstrgt-1] = 2;
                  } else if (strcmp(cparav, "D_TO_R") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 4) goto L8001;
                     indrgt[nstrgt-1] = 2;
                  } else if (strcmp(cparav, "I_TO_D") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 1) goto L8001;
                     indrgt[nstrgt-1] = 4;
                  } else if (strcmp(cparav, "R_TO_D") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 2) goto L8001;
                     indrgt[nstrgt-1] = 4;
                  } else if (strcmp(cparav, "I_TO_S") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 1) goto L8001;
                     indrgt[nstrgt-1] = 3;
                  } else if (strcmp(cparav, "I_TO_S4") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 1) goto L8001;
                     indrgt[nstrgt-1] = 3;

/*                CHECK UNARY OPERATIONS */
                  } else if (strcmp(cparav, "NOT") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 5) goto L8002;
                  } else if (strcmp(cparav, "CHS") == 0 || strcmp(cparav, "ABS") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 1 && indrgt[nstrgt-1] != 2 && indrgt[nstrgt-1] != 4) {
                        goto L8002;
                     }
                  } else if (strcmp(cparav, "EXP") == 0 || strcmp(cparav, "LN") == 0
                             || strcmp(cparav, "SIN") == 0 || strcmp(cparav, "COS") == 0
                             || strcmp(cparav, "TAN") == 0 || strcmp(cparav, "ARCSIN") == 0
                             || strcmp(cparav, "ARCCOS") == 0 || strcmp(cparav, "ARCTAN") == 0
                             || strcmp(cparav, "SQRT") == 0) {
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 2 && indrgt[nstrgt-1] != 4) goto L8003;

/*                CHECK BINARY OPERATIONS */
                  } else if (strcmp(cparav, "_MIN_") == 0 || strcmp(cparav, "_MAX_") == 0) {
                     --nstrgt;
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 1 && indrgt[nstrgt-1] != 2 && indrgt[nstrgt-1] != 4) {
                        goto L8006;
                     }
                     if (indrgt[nstrgt-1] != indrgt[nstrgt]) {
                        goto L8006;
                     }
                  } else if (strcmp(cparav, "_TRIM_") == 0) {
                     --nstrgt;
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 2 && indrgt[nstrgt-1] != 4) {
                        goto L8007;
                     }
                     if (indrgt[nstrgt] != 1) {
                        goto L8007;
                     }
                  } else if (strcmp(cparav, "+") == 0 || strcmp(cparav, "-") == 0) {
                     --nstrgt;
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != indrgt[nstrgt]) goto L8004;
                  } else if (strcmp(cparav, "*") == 0 || strcmp(cparav, "/") == 0) {
                     --nstrgt;
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] == 3) goto L8004;
                     if (indrgt[nstrgt-1] != indrgt[nstrgt]) goto L8004;
                  } else if (strcmp(cparav, "%") == 0) {
                     --nstrgt;
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 1) goto L8004;
                     if (indrgt[nstrgt-1] != indrgt[nstrgt]) goto L8004;
                  } else if (strcmp(cparav, "**") == 0) {
                     --nstrgt;
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] != 1 && indrgt[nstrgt-1] != 2 && indrgt[nstrgt-1] != 4) goto L8004;
                     if (indrgt[nstrgt-1] != indrgt[nstrgt]) goto L8004;
                  } else if (strcmp(cparav, "<") == 0 || strcmp(cparav, ">") == 0
                             || strcmp(cparav, "=") == 0 || strcmp(cparav, "<=") == 0
                             || strcmp(cparav, ">=") == 0 || strcmp(cparav, "<>") == 0) {
                     --nstrgt;
                     if (nstrgt <= 0) goto L9006;
                     if (indrgt[nstrgt-1] == 5) goto L8005;
                     if (indrgt[nstrgt-1] != indrgt[nstrgt]) goto L8005;
                     indrgt[nstrgt-1] = 5;
                  } else {
/*                   RIGHT-SIDE VARIABLES (FOR ... := ... ; ) */
/*                   USING *CPARAV* VARIABLE, SCAN ALL DECLARED VARIABLES */
                     int_32 ilowrc = ivabeg;
                     int_32 ihigrc = ivaend + 1;
L41:
                     if (ihigrc - ilowrc <= 1) {
/*                      VARIABLE NOT FOUND */
                        if (cparav[0] != '$') goto L5004;
/*                      ONLY INTERESTED IN PARAMETRIC CONSTANTS */
                        ++ivaend;

/*                      SHIFT GREATER VARIABLES */
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

/*                      NEW PARAMETRIC CONSTANT FOUND ? => CALL *CLECST* */
                        iretcd = (*dumcst)(cparav, &indlin, &idatin, &adatin, cdatin, &ddatin);
                        if (iretcd != 0) goto L5005;

/*                      VALID PARAMETRIC CONSTANT, WRITE AT END */
                        iofset = (ihigrc - 1) * lrclen;
                        strcpy(record2.cparin, cparav);
                        strcpy(record2.cdatin, cdatin);
                        record2.indlin = indlin;
                        record2.idatin = idatin;
                        record2.adatin = adatin;
                        record2.ddatin = ddatin;
                        record2.idclin = ilines;
                        record2.idefin = ilines;
                        record2.iusein = ilines;
                        iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                        if (iretcd != 0) goto L9001;
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
                           if (iusein == 0) {
/*                            CHANGE FIRST USED LINE FOR VARIABLE */
                              record2.iusein = ilines;
                              iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                              if (iretcd != 0) goto L9001;
                           }
                        } else if (strcmp(cparin, cparav) < 0) {
                           ilowrc = imedrc;
                           goto L41;
                        } else {
                           ihigrc = imedrc;
                           goto L41;
                        }
                     }
                     if (nstrgt == ndimst) goto L9005;
                     ++nstrgt;
                     indrgt[nstrgt-1] = abs(indlin);
                  }
               } else {
                  if (strcmp(cparav, cequal) == 0) {
                     lequal = 1;
/*                   FIRST DEFINED LINE FOR LEFT ( := ) VARIABLES FOR */
/*                   KEYWORDS: *INTEGER+/+REAL+/+STRING+/+DOUBLE+/+LOGICAL* */
                     if (ilogin != 6) {
                        for (jloop2 = 1; jloop2 <= nstlft; ++jloop2) {
                           iofset = (irclft[jloop2 - 1] - 1) * lrclen;
                           iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                           if (iretcd != 0) goto L9003;
                           record2.idefin = record2.idclin;
                           iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                           if (iretcd != 0) goto L9001;
                        }
                     }
                  } else {
                     int_32 ilowrc = ivabeg;
                     int_32 ihigrc = ivaend + 1;

/*                   DECLARATION FOR A VARIABLE */
/*                   ( PARAMETRIC CONSTANTS SHOULD NOT BE DECLARED) */
                     if (cparav[0] == '$') goto L5002;
/*                   NEW VARIABLE NAME FOUND ? */
L44:
/*                   BINARY SEARCH (WITH POSSIBLE INSERTION) */
                     if (ihigrc - ilowrc <= 1) {
/*                      NEW VARIABLE NAME FOUND (INVALID *EVALUATE*) */
                        ++ivaend;
                        if (ilogin == 6) goto L5004;

/*                      CHECK IF VARIABLE NAME COMPLIES WITH THE RULES */
                        for (jloop2 = 2; jloop2 <= strlen(cparav); ++jloop2) {
                           char cc[] = {cparav[jloop2-1], '\0'};
                           int_32 jin1 = index_f(alphab, cc);
                           int_32 jin2 = index_f(digits, cc);
                           int_32 jin3 = index_f(cc, " ");
                           if (jin1 + jin2 + jin3 == 0) {
                              printf("%s: CHARACTER *%c* IS NOT ALLOWED\n", nomsub, cc[0]);
                              goto L5001;
                           }
                        }

/*                      CHECK IF VARIABLE NAME IS A KEYWORD */
                        for (jloop2 = 1; jloop2 <= 40; ++jloop2) {
                           if (strcmp(cparav, ckeywd[jloop2-1]) == 0) {
                              printf("%s: VARIABLE *%s* IS A CLE-2000 KEYWORD\n", nomsub, cparav);
                              goto L5001;
                           }
                        }

/*                      SHIFT GREATER VARIABLES */
                        if (ihigrc != ivaend) {
                           for (jloop2 = ivaend - 1; jloop2 >= ihigrc; --jloop2) {
                              iofset = (jloop2 - 1) * lrclen;
                              iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                              if (iretcd != 0) goto L9003;
                              iofset = jloop2 * lrclen;
                              iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                              if (iretcd != 0) goto L9001;
                           }
                           if (nstlft != 0) {
                              for (jloop2 = 0; jloop2 < nstlft; ++jloop2) {
                                 if (irclft[jloop2] >= ihigrc) ++irclft[jloop2];
                              }
                           }
                        }
                        for (i = 0; i < 72; i++) cdatin[i] = ' ';
                        cdatin[72] = '\0';
                        indlin = -ilogin;
                        idatin = 0;
                        adatin = 0.f;
                        ddatin = 0.;
                        idclin = ilines;
                        idefin = 0;
                        iusein = 0;

/*                      VALID VARIABLE NAME, WRITE AT *IHIGRC* */
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

/*                      COUNT LEFT VARIABLES (FOR ... := ... ; ) AND */
/*                      KEEP RECORD NUMBER (IN CASE OF := DEFINITION) */
                        if (nstlft == ndimst) goto L9005;
                        ++nstlft;
                        irclft[nstlft-1] = ihigrc;
                        indlft[nstlft-1] = ilogin;
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
                           if (ilogin == 6) {
                              if (idefin == 0) {
/*                               CHANGE FIRST DEFINED LINE FOR VARIABLE */
                                 record2.idefin = ilines;
                                 iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                                 if (iretcd != 0) goto L9001;
                              } else if (nstlft != 0) {
/*                               VARIABLE DEFINED ONCE IN *EVALUATE* ? */
                                 for (jloop2 = 1; jloop2 <= nstlft; ++jloop2) {
                                    if (imedrc == irclft[jloop2 - 1]) goto L5007;
                                 }
                              }
                              if (nstlft == ndimst) goto L9005;
                              ++nstlft;
                              indlft[nstlft-1] = abs(indlin);
                              irclft[nstlft-1] = imedrc;
                           } else {
/*                            VARIABLE ALREADY DECLARED */
                              goto L5003;
                           }
                        } else if (strcmp(cparin, cparav) < 0) {
                           ilowrc = imedrc;
                           goto L44;
                        } else {
                           ihigrc = imedrc;
                           goto L44;
                        }
                     }
                  }
               }
            } else {
               if (!lequal) goto L5001;
               if (nstrgt == ndimst) goto L9005;
               ++nstrgt;
               indrgt[nstrgt-1] = jndlec[iloop1-1] + 1;
            }
         }
      }

L100:
      ;
   }
/* ***  MAIN LOOP OVER RECORDS (END) */

/* ALL VARIABLES ARE NOW SORTED AT THE END OF THE OBJECT FILE */

/* REWRITE TOP OF OBJECT FILE TO UPDATE *NSTACK+/+NRECOR* */
   nstack = ivaend - ivabeg;
   nrecor = ivaend;
   header.nrecor = nrecor;
   header.nstack = nstack;
   iretcd = kdiput_c(iunito, (int_32 *)&header, 0, kdisize(header));
   if (iretcd != 0) goto L9001;

L666:
   return ret_val;

L5000:
   printf("%-72s   LINE\n", cerror);
   printf("%-72s   %04d\n", myreco, (int)ilines);
   goto L666;
L5001:
   printf("! %s: INVALID VARIABLE NAME IN RECORD\n", nomsub);
   ret_val = 5001;
   goto L5000;
L5002:
   printf("! %s: *%s* CANNOT BE DECLARED\n", nomsub, cparav);
   ret_val = 5002;
   goto L5000;
L5003:
   printf("! %s: VARIABLE DECLARED TWICE *%s*\n", nomsub, cparav);
   ret_val = 5003;
   goto L5000;
L5004:
   printf("! %s: VARIABLE NOT YET DECLARED *%s*\n", nomsub, cparav);
   ret_val = 5004;
   goto L5000;
L5005:
   printf("! %s: INVALID PARAMETER *%s*\n", nomsub, cparav);
   ret_val = 5005;
   goto L5000;
L5006:
   printf("! %s: INVALID VARIABLE FOR >>.<< OR <<.>>\n", nomsub);
   ret_val = 5006;
   goto L5000;
L5007:
   printf("! %s: VARIABLE EVALUATED TWICE *%s*\n", nomsub, cparav);
   ret_val = 5007;
   goto L5000;
L8000:
   printf("! %s: *%s* WITH TYPE %s\n", nomsub, cparav, ckeywd[indrgt[nstrgt-1]-1]);
   goto L5000;
L8001:
   printf("! %s: INVALID TYPE_TO_TYPE CONVERSION\n", nomsub);
   ret_val = 8001;
   goto L8000;
L8002:
   printf("! %s: INVALID *NOT* OR *ABS*\n", nomsub);
   ret_val = 8002;
   goto L8000;
L8003:
   printf("! %s: INVALID TYPE FOR REAL/DOUBLE FUNCTION\n", nomsub);
   ret_val = 8003;
   goto L8000;
L8004:
   printf("! %s: INVALID TYPE FOR +,-,*,/,modulo OR **\n", nomsub);
   ret_val = 8004;
   goto L8000;
L8005:
   printf("! %s: INVALID TYPE FOR <,>,=,<=,>= OR <>\n", nomsub);
   ret_val = 8005;
   goto L8000;
L8006:
   printf("! %s: INVALID TYPE FOR _MIN_ OR _MAX_\n", nomsub);
   ret_val = 8006;
   goto L8000;
L8007:
   printf("! %s: INVALID TYPE FOR _TRIM_\n", nomsub);
   ret_val = 8006;
   goto L8000;
L9001:
   iretcd = -1;
   printf("! %s: WRITING RETURN CODE =%d\n", nomsub, (int)iretcd);
   ret_val = iretcd;
   goto L666;
L9003:
   iretcd = -1;
   printf("! %s: WRITING RETURN CODE =%d\n", nomsub, (int)iretcd);
   ret_val = iretcd;
   goto L666;
L9005:
   printf("! %s: *STACK* MEMORY IS FULL\n", nomsub);
   ret_val = 9005;
   goto L5000;
L9006:
   printf("! %s: *STACK* MEMORY IS EMPTY\n", nomsub);
   ret_val = 9006;
   goto L5000;
L9008:
   printf("! %s: ERROR ON THE NUMBER OF EVALUATIONS\n", nomsub);
   printf("! %s: CLESTK: LEFT=%d VS. RIGHT=%d\n", nomsub, (int)nstlft, (int)nstrgt);
   ret_val = 9008;
   goto L5000;
L9009:
   printf("! %s: ERROR ON THE TYPE OF AN EVALUATION\n", nomsub);
   ret_val = 9009;
   goto L5000;
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
} /* clestk */
