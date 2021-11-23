
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 24/04/09      */
/*****************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cle2000.h"
#include "header.h"
#define sign(A)  (A > 0 ? 1 : (A < 0 ? -1 : 0))
#define index_f(A, B) (strstr(A, B) == NULL ? 0 : strstr(A, B) - A + 1)
#define nlogkw 15
#define ndimst 128
#define ifinal 74
#define nmawrd 36

static kdi_file *iunito = NULL;
static int_32 ivabeg = 0;
static int_32 ivaend = 0;
static int_32 ioulst = 0;
static int_32 ilogin = 0;
static FILE *iwrite = NULL;
static char hwrite[73] = " ";
static int_32 idblst = 0;
static int_32 irecor = 0;
static int_32 ninput = 0;
static int_32 nwords = 0;
static int_32 iwords = 0;
static int_32 nstput = 0;
static int_32 nstlvl = 0;
static int_32 indrgt[ndimst], irclvl[ndimst], irclft[ndimst], ircput[ndimst];
static int_32 idebwd[nmawrd+1], ifinwd[nmawrd+1], jndlec[nmawrd];
static char myreco[73];
static char cl2000[12] = "CLE2000(V3)";
static int_32 text_len = 72; /* CAUTION: text must be declared as text[73] in the calling function*/

static int_32 intstk[ndimst];
static float_32 relstk[ndimst];
static char chrstk[ndimst][73];
static double_64 dblstk[ndimst];
static int_32 logstk[ndimst];

void redget_c(int_32 *ityp, int_32 *nitma, float_32 *flott, char text[73], double_64 *dflot)
{
   static char *clogbg[] = {"INTEGER", "REAL", "STRING", "DOUBLE", "LOGICAL", "EVALUATE",
                            "ECHO", "ELSEIF", "IF", "WHILE", "UNTIL", "ENDWHILE", "REPEAT",
                            "ELSE", "ENDIF"};
   static char *clognd[] = {";", ";", ";", ";", ";", ";", ";", "THEN", "THEN", "DO",
                            ";", ";", "REPEAT", "ELSE", ";"};
   int_32 i, ilines, jlines=0, jlevel=0, jrecor, iretcd, iloop1, lrgtst=0, nstlft=0, imedrc,
           iofset, idefkw=0, nstrgt=0;
   int_32 maskck[nmaskc], ipacki[nmaskc];
   char cparin[17], cparav[13], chrend[13];
   int_32 ilengv, indlec;
   char cdatav[73];
   int_32 idatin, indlin;
   float_32 adatin;
   char cdatin[73];
   double_64 ddatin;
   char *nomsub="redget_c";

/* TAKE ANY LEVEL AS INPUT */
   int_32 ilevel = -1;
/* L01> ***LOOP*** OVER WORDS (BEGIN) */
L10:
   ++iwords;
   if (iwords > nwords) {
      int_32 jbiprv;
/*    L02>    ***LOOP*** OVER RECORDS (BEGIN) */
L20:
      ++irecor;
      if (irecor > ninput) {
         if (ninput == 0) {
/*          REDGET IS CLOSED */
            *ityp = 10;
            if (idblst > 0 && iwrite != NULL) {
               sprintf(myreco,"READER IS CLOSED ON FILE");
               fprintf(iwrite,".|%-72s|.\n",myreco);
            }
         } else {
            int_32 i;
            *ityp = 9;
            if (idblst > 0 && iwrite != NULL) {
               sprintf(myreco,"QUIT \"DEBUG\" ");
               fprintf(iwrite,".|%-72s|.\n",myreco);
            }
            for ( i = 0; i < 72; i++) myreco[i]='-';
            myreco[72]='\0'; fprintf(iwrite,". %s .\n",myreco);
         }
         return;
      }
/*    READ A NEW RECORD */
      iofset = (irecor-1) * lrclen;
      iretcd = kdiget_c(iunito, (int_32 *)&record1, iofset, kdisize(record1));
      if (iretcd != 0) goto L8000;
      strcpy(cparin, record1.cparin);
      strcpy(myreco, record1.myreco);
      ilines = record1.ilines;
      jlevel = record1.ilevel;
      jrecor = record1.irecor;
      for (i = 0; i < nmaskc; i++) maskck[i] = record1.maskck[i];
      for (i = 0; i < nmaskc; i++) ipacki[i] = record1.ipacki[i];
      if (ilevel >= 0 && jlevel != ilevel) goto L20;

/* L02> ***LOOP*** OVER RECORDS (ACCEPTED) */

/*    RECORD IS ACCEPTED, PRINT WHEN REQUESTED */
      if (idblst <= 0 && jlevel == 0) {
         if (ioulst > 0 && iwrite != NULL) fprintf(iwrite,"<|%-72s|<%04d\n",myreco,(int)ilines);
      } else if (idblst > 0 && iwrite != NULL) {
         if (jlevel == 0) {
            fprintf(iwrite,"<|%-72s|<%04d\n",myreco,(int)ilines);
         } else {
            fprintf(iwrite,".|%-72s|.%04d\n",myreco,(int)ilines);
         }
      }

/*    BEGIN: MASK RECOVERY */
      jbiprv = 0;
      iwords=1;
      nwords = 1;
      for (iloop1=0; iloop1<72; ++iloop1) {
         int_32 jloop2 = iloop1/24;
         int_32 jbicur = maskck[jloop2] % 2;
         iwords += jbiprv * (1-jbicur);
         idebwd[nwords-1] = iloop1;
         ifinwd[iwords-1] = iloop1;
         nwords += jbicur * (1-jbiprv);
         jbiprv = jbicur;
         maskck[jloop2] /= 2;
      }
      --nwords;
/*    END:   MASK RECOVERY */

/*    THIS IS NOW THE FIRST WORD */
      iwords = 1;
      if (jlevel != 0) {
/* L03>  THIS IS A NEW *CLE-2000* RECORD (BEGIN) */
         if (strcmp(cparin, " ") != 0) {

/* L04>     TREAT SIGNIFICANT 1ST-WORD OF CLE-2000 STATEMENT (BEGIN) */
            nstlft = 0;
            nstrgt = 0;
            for (iloop1 = 0; iloop1 < nlogkw; ++iloop1) {
               if (strcmp(clogbg[iloop1], cparin) == 0) ilogin = iloop1+1;
            }
            lrgtst = ilogin >= 7 && ilogin <= 11;
            strcpy(chrend, clognd[ilogin-1]);

/*          KEYWORDS: *IF+/+WHILE+/+UNTIL* */
            if (ilogin == 9 || ilogin == 10 || ilogin == 11) {
               ++nstlvl;
               if (nstlvl > ndimst) goto L9002;
               irclvl[nstlvl - 1] = jrecor - 1;

/*          KEYWORDS: *ELSEIF+/+ELSE* */
            } else if (ilogin == 8 || ilogin == 14) {
               if (ilevel == -1) {
                  if (nstlvl == 0) goto L9001;
                  irecor = irclvl[nstlvl-1];
                  nwords = 1;
               } else {
                  ilevel = -1;
               }

/*          KEYWORD:  *ENDWHILE* */
            } else if (ilogin == 12) {
               if (ilevel == -1) {
                  irecor = jrecor - 1;
                  nwords = 1;
               } else {
                  ilevel = -1;
               }

/*          KEYWORD:  *ENDIF* */
            } else if (ilogin == 15) {
               ilevel = -1;
               if (nstlvl == 0) goto L9001;
               --nstlvl;

/*          KEYWORD:  *ECHO* */
            } else if (ilogin == 7) {
               jlines = ilines;
            }

/*          CYCLE ON WORDS WITHOUT UNPACKING (ONE WORD ONLY) */
            if (nwords == iwords) goto L10;
            for (iloop1 = 1; iloop1 <= nwords; ++iloop1) {
/* L05>        UNPACK INDLEC ONLY IF MORE THAN 1 WORD (BEGIN) */
               int_32 jloop2 = ((iloop1 << 1) + 23) / 24;
               jndlec[iloop1-1] = ipacki[jloop2-1] % 4 + 1;
               ipacki[jloop2-1] /= 4;
/* L05>        UNPACK INDLEC ONLY IF MORE THAN 1 WORD (END) */
            }

/* L04>     TREAT SIGNIFICANT 1ST-WORD OF CLE-2000 STATEMENTS (END) */
/*          RETURN TO NEXT WORD */
            goto L10;
         }
/* L03>  THIS IS A NEW *CLE-2000* RECORD (END) */

      } else {
/* L03>  RECORD OUTSIDE *CLE-2000* (BEGIN) */
         ilogin = 0;
/* L03>  RECORD OUTSIDE *CLE-2000* (END) */
      }
      for (iloop1 = 1; iloop1 <= nwords; ++iloop1) {
/* L03>  RECORD OUTSIDE *CLE-2000* OR NEW *CLE-2000* RECORD, */
/*       BUT CONTINUATION STATEMENT (BEGIN) THEN, ALWAYS UNPACK INDLEC */
         int_32 jloop2 = ((iloop1 << 1) + 23) / 24;
         jndlec[iloop1-1] = ipacki[jloop2-1] % 4 + 1;
         ipacki[jloop2-1] /= 4;
/* L03>  RECORD OUTSIDE *CLE-2000* OR NEW *CLE-2000* RECORD, */
/*       BUT CONTINUATION STATEMENT (END) */
      }
/* L02>  ***LOOP*** OVER RECORDS (END) */
   }
/* L01> ***LOOP*** OVER WORDS (ACCEPTED) */

/* DETERMINE NEXT WORD */
   ilengv = ifinwd[iwords-1] - idebwd[iwords-1] + 1;
   indlec = jndlec[iwords-1];

   if (indlec == 3) {
      for ( i = 0; i < ilengv; i++) cdatav[i]=myreco[idebwd[iwords-1]+i];
      cdatav[ilengv] = '\0';
   } else if (indlec == 1) {
      for ( i = 0; i < ilengv; i++) cdatin[i]=myreco[idebwd[iwords-1]+i];
      cdatin[ilengv] = '\0';
      sscanf(cdatin, "%d", (int *)&idatin);
   } else if (indlec == 2) {
      for ( i = 0; i < ilengv; i++) cdatin[i]=myreco[idebwd[iwords-1]+i];
      cdatin[ilengv] = '\0';
      sscanf(cdatin, "%e", &adatin);
   } else {
      int_32 id;
      for ( i = 0; i < ilengv; i++) cdatin[i]=myreco[idebwd[iwords-1]+i];
      cdatin[ilengv] = '\0';
      id = index_f(cdatin, "D");
      if (id > 0) cdatin[id-1] = 'E';
      sscanf(cdatin, "%le", &ddatin);
   }
   if (ilogin == 0) {

/* L02> WORDS OUTSIDE *CLE2000* STATEMENTS: HIT AND RUN... */
      if (indlec == 3) {

/* L03>  STRING, <<.>> OR >>.<< TREATMENT */
         int_32 lrdput = strncmp(cdatav, ">>", 2) == 0;
         if (strncmp(cdatav, "<<", 2) == 0 || lrdput) {
            int_32 ilowrc = ivabeg;
            int_32 ihigrc = ivaend;

/* L04>     CASES <<.>> OR >>.<< */
/*          GET VARIABLE WITH A BINARY SEARCH IN SORTED FILE */
/*          SET UPPER AND LOWER BOUNDS */
            strcpy(cparin, &cdatav[2]);
            strncpy(cparav, cparin, ilengv-4); cparav[ilengv-4] = '\0';
L11:
            imedrc = (ihigrc+ilowrc) / 2;
            iofset = (imedrc-1) * lrclen;
            iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
            if (iretcd != 0) goto L9023;
            strcpy(cparin, record2.cparin);
            strcpy(cdatin, record2.cdatin);
            indlin = record2.indlin;
            idatin = record2.idatin;
            adatin = record2.adatin;
            ddatin = record2.ddatin;
            if (iretcd != 0) goto L9023;
            if (strcmp(cparin, cparav) == 0) {
               if (lrdput) {
                  int_32 ilengt = min(text_len, 12);
/*                KEEP RECORD NUMBER FOR *REDPUT* */
                  ++nstput;
                  if (nstput > ndimst) goto L9004;
                  ircput[nstput-1] = imedrc;
/*                MAKE THE VARIABLE UNDEFINED */
                  indlin = -abs(indlin);
                  record2.indlin = indlin;
                  iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                  if (iretcd != 0) goto L9023;
/*                SEND BACK NEGATIVE TYPE AND VARIABLE NAME TO THE APPLICATION */
                  *ityp = indlin;
                  strncpy(text, cparin, ilengt); text[ilengt]='\0';
               } else {
                  if (indlin <= 0) goto L9008;
/*                SEND BACK TYPE AND DEFINED VALUES */
                  *ityp = indlin;
                  if (*ityp == 1 || *ityp == 5) {
                     *nitma = idatin;
                  } else if (*ityp == 2) {
                     *flott = adatin;
                  } else if (*ityp == 3) {
                     int_32 ilengt = min(text_len, 72);
                     *nitma = min(idatin,ilengt);
                     strncpy(text, cdatin, ilengt); text[ilengt]='\0';
                  } else if (*ityp == 4) {
                     *dflot = ddatin;
                  }
               }
            } else if (strcmp(cparin, cparav) < 0) {
               ilowrc = imedrc;
               goto L11;
            } else {
               ihigrc = imedrc;
               goto L11;
            }
         } else if (cdatav[0] == '\'') {
            int_32 ilengt = min(text_len, 72);
            if (ilengv > 2) {
               strncpy(cdatin, &cdatav[1], ilengv-2); cdatin[ilengv-2]='\0';
            }
            *ityp = 3;
            *nitma = min(ilengv-2, ilengt);
            strncpy(text, cdatin, ilengt); text[ilengt] = '\0';
         } else {
            int_32 ilengt = min(text_len, 72);
            *ityp = 3;
            *nitma = min(ilengv,ilengt);
            strncpy(text, cdatav, ilengt); text[ilengt] = '\0';
         }
      } else {
/* L03>  OTHER THAN STRING TREATMENT */
         *ityp = indlec;
         if (*ityp == 1 || *ityp == 5) {
            *nitma = idatin;
         } else if (*ityp == 2) {
            *flott = adatin;
         } else if (*ityp == 4) {
            *dflot = ddatin;
         }
      }
      return;

/* L02> WORDS OUTSIDE *CLE2000* STATEMENTS: END. */

   } else {

/* L02> PROCESS *CLE2000* STATEMENTS: DRINK, DRIVE (BEGIN) */
/*      WATCH FOR STRINGS... */

/* L03> IF( INDLEC.EQ.3.AND.CDATAV(1:1).EQ.'"' )THEN */
      if (indlec == 3 && cdatav[0] == '"') {
         ++nstrgt;
         indrgt[nstrgt-1] = indlec;
         if (ilengv > 2) {
            strncpy(chrstk[nstrgt-1], &cdatav[1], ilengv-2); chrstk[nstrgt-1][ilengv-2] = '\0';
         }
         intstk[nstrgt-1] = ilengv-2;

/* L03> ELSEIF( LRGTST )THEN */
      } else if (lrgtst) {

/* L04>  IF( INDLEC.EQ.3 )THEN */
         if (indlec == 3) {
            strncpy(cparav, cdatav, 12); cparav[12] = '\0';

/* L05>     IF( CPARAV.EQ.CHREND )THEN */
            if (strcmp(cparav, chrend) == 0) {

/*             TRUEWAY LEFT/RIGHT */
/*             KEYWORDS: *int_32+/+REAL+/+STRING+/+DOUBLE+/+LOGICAL+/+EVALUATE* */
               if (ilogin <= 6) {
/*                PUT VARIABLE VALUES */
L25:
                  indlin = indrgt[nstlft-1];
                  imedrc = irclft[nstlft-1];
                  iofset = (imedrc-1) * lrclen;
                  iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                  if (iretcd != 0) goto L9023;
                  strcpy(cdatin, record2.cdatin);
                  idatin = record2.idatin;
                  adatin = record2.adatin;
                  ddatin = record2.ddatin;
                  if (indlin == 1) {
                     idatin = intstk[nstlft-1];
                  } else if (indlin == 2) {
                     adatin = relstk[nstlft-1];
                  } else if (indlin == 3) {
                     strcpy(cdatin, chrstk[nstlft-1]);
                     idatin = intstk[nstlft-1];
                  } else if (indlin == 4) {
                     ddatin = dblstk[nstlft-1];
                  } else if (indlin == 5) {
                     if (logstk[nstlft - 1]) {
                        idatin = 1;
                     } else {
                        idatin = -1;
                     }
                  }
                  strcpy(record2.cdatin, cdatin);
                  record2.indlin = indlin;
                  record2.idatin = idatin;
                  record2.adatin = adatin;
                  record2.ddatin = ddatin;
                  iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                  if (iretcd != 0) goto L9023;
                  --nstlft;
                  if (nstlft != 0) goto L25;

/*             KEYWORD: *ECHO* (PRINTER UTILITY ) */
               } else if (ilogin == 7) {
                  char cprint[81];
                  int_32 ilgtpr=0, idepar=3;

                  sprintf(cprint, ">|%-72s|>%04d", " ", (int)jlines);
                  for (iloop1 = 0; iloop1 < nstrgt; ++iloop1) {
                     indlin = indrgt[iloop1];
                     if (indlin == 1) {
                        int_32 jloop2, indle2 = abs(intstk[iloop1]);
                        for (jloop2 = 1; jloop2 <= 12; ++jloop2) {
                           indle2 /= 10;
                           if (indle2 == 0) {
                              ilgtpr = jloop2;
                              goto L302;
                           }
                        }
L302:
                        if (intstk[iloop1] < 0) ++ilgtpr;
                        sprintf(cdatin, "%d", (int)intstk[iloop1]);
                     } else if (indlin == 2) {
                        ilgtpr = 13;
                        sprintf(cdatin, "%13.6e", relstk[iloop1]);
                     } else if (indlin == 3) {
                        ilgtpr = intstk[iloop1];
                        strcpy(cdatin, chrstk[iloop1]);
                     } else if (indlin == 4) {
                        ilgtpr = 23;
                        sprintf(cdatin, "%23.15E", dblstk[iloop1]);
                     } else if (indlin == 5) {
                        if (logstk[iloop1]) {
                           ilgtpr = 6;
                           strncpy(cdatin, ".TRUE.", 6); cdatin[6] = '\0';
                        } else {
                           ilgtpr = 7;
                           strncpy(cdatin, ".FALSE.", 7); cdatin[7] = '\0';
                        }
                     } else if (indlin == -1) {
                        ilgtpr = 3;
                        strncpy(cdatin, "?_I", 3); cdatin[3] = '\0';
                     } else if (indlin == -2) {
                        ilgtpr = 3;
                        strncpy(cdatin, "?_R", 3); cdatin[3] = '\0';
                     } else if (indlin == -3) {
                        ilgtpr = 3;
                        strncpy(cdatin, "?_S", 3); cdatin[3] = '\0';
                     } else if (indlin == -4) {
                        ilgtpr = 3;
                        strncpy(cdatin, "?_D", 3); cdatin[3] = '\0';
                     } else if (indlin == -5) {
                        ilgtpr = 3;
                        strncpy(cdatin, "?_L", 3); cdatin[3] = '\0';
                     } else {
                        goto L9007;
                     }
                     if (idepar + ilgtpr >= ifinal) {
                        if (iwrite != NULL) fprintf(iwrite,"%s\n",cprint);
                        sprintf(cprint, ">|%-72s|>%04d", " ", (int)jlines);
                        idepar = 3;
                     }
                     strncpy(&cprint[idepar-1], cdatin, ilgtpr); cprint[80] = '\0';
                     idepar = idepar + ilgtpr + 1;
                     if (idepar >= ifinal) {
                        if (iwrite != NULL) fprintf(iwrite,"%s\n",cprint);
                        sprintf(cprint, ">|%-72s|>%04d", " ", (int)jlines);
                        idepar = 3;
                     }
                  }
                  if (iwrite != NULL && idepar != 3) fprintf(iwrite,"%s\n",cprint);
                  fflush(iwrite);

/*             KEYWORDS: *ELSEIF+/+IF* */
               } else if (ilogin == 8 || ilogin == 9) {
                  if (indrgt[nstrgt-1] != 5) goto L9006;
                  if (logstk[nstrgt-1]) {
                     ilevel = -1;
                  } else {
                     ilevel = jlevel;
                  }

/*             KEYWORD: *UNTIL* */
               } else if (ilogin == 11) {
                  if (nstlvl == 0) goto L9001;
                  if (indrgt[nstrgt-1] != 5) goto L9006;
                  if (!logstk[nstrgt-1]) {
                     irecor = irclvl[nstlvl-1];
                     iwords = nwords;
                  }
                  --nstlvl;

/*                  KEYWORD: *WHILE* */
               } else if (ilogin == 10) {
                  if (nstlvl == 0) goto L9001;
                  if (indrgt[nstrgt-1] != 5) goto L9006;
                  if (!logstk[nstrgt-1]) {
                     ilevel = jlevel;
                     irecor = irclvl[nstlvl-1];
                     iwords = nwords;
                  }
                  --nstlvl;
               }
            } else {
/*             CHECK CONVERSION OPERATIONS */
               if (strcmp(cparav, "R_TO_I") == 0) {
                  indrgt[nstrgt-1] = sign(indrgt[nstrgt-1]);
                  intstk[nstrgt-1] = (int_32) relstk[nstrgt-1];
               } else if (strcmp(cparav, "D_TO_I") == 0) {
                  indrgt[nstrgt-1] = sign(indrgt[nstrgt-1]);
                  intstk[nstrgt-1] = (int_32) dblstk[nstrgt-1];
               } else if (strcmp(cparav, "I_TO_R") == 0) {
                  indrgt[nstrgt-1] = 2 * sign(indrgt[nstrgt-1]);
                  relstk[nstrgt-1] = (float_32) intstk[nstrgt-1];
               } else if (strcmp(cparav, "D_TO_R") == 0) {
                  indrgt[nstrgt-1] = 2 * sign(indrgt[nstrgt-1]);
                  relstk[nstrgt-1] = (float_32) dblstk[nstrgt-1];
               } else if (strcmp(cparav, "I_TO_D") == 0) {
                  indrgt[nstrgt-1] = 4 * sign(indrgt[nstrgt-1]);
                  dblstk[nstrgt-1] = (double_64) intstk[nstrgt-1];
               } else if (strcmp(cparav, "R_TO_D") == 0) {
                  indrgt[nstrgt-1] = 4 * sign(indrgt[nstrgt-1]);
                  dblstk[nstrgt-1] = (double_64) relstk[nstrgt-1];
               } else if (strcmp(cparav, "I_TO_S") == 0) {
                  indrgt[nstrgt-1] = 3 * sign(indrgt[nstrgt-1]);
                  if (intstk[nstrgt-1] > 99999999) goto L9013;
                  if (intstk[nstrgt-1] < -9999999) goto L9014;
                  sprintf(chrstk[nstrgt-1], "%d", (int)intstk[nstrgt-1]);
                  intstk[nstrgt-1] = (int)strlen(chrstk[nstrgt-1]);
               } else if (strcmp(cparav, "I_TO_S4") == 0) {
                  indrgt[nstrgt-1] = 3 * sign(indrgt[nstrgt-1]);
                  if (intstk[nstrgt-1] > 99999999) goto L9013;
                  if (intstk[nstrgt-1] < -9999999) goto L9014;
                  sprintf(chrstk[nstrgt-1], "%04d", (int)intstk[nstrgt-1]);
                  intstk[nstrgt-1] = (int)strlen(chrstk[nstrgt-1]);
/*             CHECK UNARY OPERATIONS */
               } else if (strcmp(cparav, "NOT") == 0) {
                  logstk[nstrgt-1] = !logstk[nstrgt-1];
               } else if (strcmp(cparav, "CHS") == 0) {
                  if (indrgt[nstrgt-1] == 1) {
                     intstk[nstrgt-1] = -intstk[nstrgt-1];
                  } else if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = -relstk[nstrgt-1];
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = -dblstk[nstrgt-1];
                  }
               } else if (strcmp(cparav, "ABS") == 0) {
                  if (indrgt[nstrgt-1] == 1 && intstk[nstrgt-1] < 0) {
                     intstk[nstrgt-1] = -intstk[nstrgt-1];
                  } else if (indrgt[nstrgt-1] == 2 && relstk[nstrgt-1] < 0.f) {
                     relstk[nstrgt-1] = -relstk[nstrgt-1];
                  } else if (indrgt[nstrgt-1] == 4 && dblstk[nstrgt-1] < 0.) {
                     dblstk[nstrgt-1] = -dblstk[nstrgt-1];
                  }
               } else if (strcmp(cparav, "EXP") == 0) {
                  if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = exp(relstk[nstrgt-1]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = exp(dblstk[nstrgt-1]);
                  }
               } else if (strcmp(cparav, "LN") == 0) {
                  if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = log(relstk[nstrgt-1]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = log(dblstk[nstrgt-1]);
                  }
               } else if (strcmp(cparav, "SIN") == 0) {
                  if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = sin(relstk[nstrgt-1]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = sin(dblstk[nstrgt-1]);
                  }
               } else if (strcmp(cparav, "COS") == 0) {
                  if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = cos(relstk[nstrgt-1]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = cos(dblstk[nstrgt-1]);
                  }
               } else if (strcmp(cparav, "TAN") == 0) {
                  if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = tan(relstk[nstrgt-1]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = tan(dblstk[nstrgt-1]);
                  }
               } else if (strcmp(cparav, "ARCSIN") == 0) {
                  if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = asin(relstk[nstrgt-1]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = asin(dblstk[nstrgt-1]);
                  }
               } else if (strcmp(cparav, "ARCCOS") == 0) {
                  if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = acos(relstk[nstrgt-1]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = acos(dblstk[nstrgt-1]);
                  }
               } else if (strcmp(cparav, "ARCTAN") == 0) {
                  if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = atan(relstk[nstrgt-1]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = atan(dblstk[nstrgt-1]);
                  }
               } else if (strcmp(cparav, "SQRT") == 0) {
                  if (indrgt[nstrgt-1] == 2) {
                     if (relstk[nstrgt-1] < 0.f) {
                        idefkw = 2;
                        goto L9009;
                     }
                     relstk[nstrgt-1] = sqrt(relstk[nstrgt-1]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     if (dblstk[nstrgt-1] < 0.) {
                        idefkw = 4;
                        goto L9009;
                     }
                     dblstk[nstrgt-1] = sqrt(dblstk[nstrgt-1]);
                  }

/*             CHECK BINARY OPERATIONS */
               } else if (strcmp(cparav, "_MIN_") == 0) {
                  --nstrgt;
                  indrgt[nstrgt-1] = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  if (indrgt[nstrgt-1] == 1) {
                     intstk[nstrgt-1] = min(intstk[nstrgt-1], intstk[nstrgt]);
                  } else if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = min(relstk[nstrgt-1], relstk[nstrgt]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = min(dblstk[nstrgt-1], dblstk[nstrgt]);
                  }
               } else if (strcmp(cparav, "_MAX_") == 0) {
                  --nstrgt;
                  indrgt[nstrgt-1] = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  if (indrgt[nstrgt-1] == 1) {
                     intstk[nstrgt-1] = max(intstk[nstrgt-1], intstk[nstrgt]);
                  } else if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] = max(relstk[nstrgt-1], relstk[nstrgt]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] = max(dblstk[nstrgt-1], dblstk[nstrgt]);
                  }
               } else if (strcmp(cparav, "_TRIM_") == 0) {
                  --nstrgt;
                  int_32 ndig = intstk[nstrgt];
                  if (indrgt[nstrgt-1] == 2) {
                     int_32 nx = (int_32)floor(log10(relstk[nstrgt-1]))-ndig;
                     float_32 rrr = relstk[nstrgt-1]*pow(10.0,(float_32)(-nx));
                     relstk[nstrgt-1] = floor(rrr)*pow(10.0,(float_32)(nx));
                  } else if (indrgt[nstrgt-1] == 4) {
                     int_32 nx = (int_32)floor(log10(dblstk[nstrgt-1]))-ndig;
                     double_64 ddd = dblstk[nstrgt-1]*pow(10.0,(double_64)(-nx));
                     dblstk[nstrgt-1] = floor(ddd)*pow(10.0,(double_64)(nx));
                  }
               } else if (strcmp(cparav, "+") == 0) {
                  --nstrgt;
                  indrgt[nstrgt-1] = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  if (indrgt[nstrgt-1] == 1) {
                     intstk[nstrgt-1] += intstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] += relstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 3) {
                     char cdata1[73], cdata2[73];
                     int_32 ileng2 = intstk[nstrgt-1];
                     int_32 ileng1 = intstk[nstrgt];
                     strcpy(cdata2, chrstk[nstrgt-1]);
                     strcpy(cdata1, chrstk[nstrgt]);
                     if (ileng1 == 0) {
                        if (ileng2 == 0) {
                           strcpy(chrstk[nstrgt-1], " ");
                        } else {
                           strcpy(chrstk[nstrgt-1], cdata2);
                        }
                     } else if (ileng2 == 0) {
                        strcpy(chrstk[nstrgt-1], cdata1);
                     } else if (ileng1 + ileng2 <= 72) {
                        strcpy(chrstk[nstrgt-1], cdata2);
                        strcat(chrstk[nstrgt-1], cdata1);
                     } else {
                        printf("%s: STRING IS LONGER THAN 72 CHRS", nomsub);
                        goto L9012;
                     }
                     intstk[nstrgt-1] = ileng1 + ileng2;
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] += dblstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 5) {
                     logstk[nstrgt-1] = logstk[nstrgt-1] || logstk[nstrgt];
                  }
               } else if (strcmp(cparav, "-") == 0) {
                  --nstrgt;
                  indrgt[nstrgt-1] = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  if (indrgt[nstrgt-1] == 1) {
                     intstk[nstrgt-1] -= intstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] -= relstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 3) {
                     char cdata1[73], cdata2[73];
                     int_32 ileng2 = intstk[nstrgt-1];
                     int_32 ileng1 = intstk[nstrgt];
                     if (ileng1 != 0) {
                        if (ileng2 < ileng1) {
                           printf("%s: IMPOSSIBLE TO - A SUBSTRING WITH LENGTH LESS THAN STRING", nomsub);
                           goto L9012;
                        }
                        strcpy(cdata2, chrstk[nstrgt-1]);
                        strcpy(cdata1, chrstk[nstrgt]);
                        if (strncmp(&cdata2[ileng2-ileng1], cdata1, ileng1) != 0) {
                           printf("%s: IMPOSSIBLE TO - A SUBSTRING NOT AT THE END OF A STRING", nomsub);
                           goto L9012;
                        } else if (ileng1 == ileng2) {
                           strcpy(chrstk[nstrgt-1], " ");
                        } else {
                           strncpy(chrstk[nstrgt-1], cdata2, ileng2-ileng1);
                           chrstk[nstrgt-1][ileng2-ileng1] = '\0';
                        }
                        intstk[nstrgt-1] = ileng2-ileng1;
                     }
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] -= dblstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 5) {
                     logstk[nstrgt-1] = logstk[nstrgt-1] || !logstk[nstrgt];
                  }
               } else if (strcmp(cparav, "*") == 0) {
                  --nstrgt;
                  indrgt[nstrgt-1] = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  if (indrgt[nstrgt-1] == 1) {
                     intstk[nstrgt-1] *= intstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 2) {
                     relstk[nstrgt-1] *= relstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 4) {
                     dblstk[nstrgt-1] *= dblstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 5) {
                     logstk[nstrgt-1] = logstk[nstrgt-1] && logstk[nstrgt];
                  }
               } else if (strcmp(cparav, "%") == 0) {
                  --nstrgt;
                  indrgt[nstrgt-1] = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  if (indrgt[nstrgt-1] == 1) {
                     if (intstk[nstrgt] == 0) {
                        idefkw = 1;
                        goto L9010;
                     }
                     intstk[nstrgt-1] %= intstk[nstrgt];
                  }
               } else if (strcmp(cparav, "/") == 0) {
                  --nstrgt;
                  indrgt[nstrgt-1] = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  if (indrgt[nstrgt-1] == 1) {
                     if (intstk[nstrgt] == 0) {
                        idefkw = 1;
                        goto L9010;
                     }
                     intstk[nstrgt-1] /= intstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 2) {
                     if (relstk[nstrgt] == 0.f) {
                        idefkw = 2;
                        goto L9010;
                     }
                     relstk[nstrgt-1] /= relstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 4) {
                     if (dblstk[nstrgt] == 0.) {
                        idefkw = 4;
                        goto L9010;
                     }
                     dblstk[nstrgt-1] /= dblstk[nstrgt];
                  } else if (indrgt[nstrgt-1] == 5) {
                     logstk[nstrgt-1] = logstk[nstrgt-1] && !logstk[nstrgt];
                  }
               } else if (strcmp(cparav, "**") == 0) {
                  --nstrgt;
                  indrgt[nstrgt-1] = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  if (indrgt[nstrgt-1] == 1) {
                     if (intstk[nstrgt-1] < 0 && intstk[nstrgt] < 0) {
                        idefkw = 1;
                        goto L9011;
                     }
                     intstk[nstrgt-1] = pow(intstk[nstrgt-1], intstk[nstrgt]);
                  } else if (indrgt[nstrgt-1] == 2) {
                     if (relstk[nstrgt-1] < 0.f && relstk[nstrgt] < 0.f) {
                        idefkw = 2;
                        goto L9011;
                     }
                     relstk[nstrgt-1] = pow(relstk[nstrgt-1], relstk[nstrgt]);
                  } else if (indrgt[nstrgt-1] == 4) {
                     if (dblstk[nstrgt-1] < 0. && dblstk[nstrgt] < 0.) {
                        idefkw = 4;
                        goto L9011;
                     }
                     dblstk[nstrgt-1] = pow(dblstk[nstrgt-1], dblstk[nstrgt]);
                  }
               } else if (strcmp(cparav, "<") == 0) {
                  --nstrgt;
                  indlin = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  indrgt[nstrgt-1] = 5;
                  if (indlin == 1) {
                     logstk[nstrgt-1] = intstk[nstrgt-1] < intstk[nstrgt];
                  } else if (indlin == 2) {
                     logstk[nstrgt-1] = relstk[nstrgt-1] < relstk[nstrgt];
                  } else if (indlin == 3) {
                     logstk[nstrgt-1] = strcmp(chrstk[nstrgt-1], chrstk[nstrgt]) < 0;
                  } else if (indlin == 4) {
                     logstk[nstrgt-1] = dblstk[nstrgt-1] < dblstk[nstrgt];
                  } else {
                     indrgt[nstrgt-1] = -5;
                  }
               } else if (strcmp(cparav, ">") == 0) {
                  --nstrgt;
                  indlin = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  indrgt[nstrgt-1] = 5;
                  if (indlin == 1) {
                     logstk[nstrgt-1] = intstk[nstrgt-1] > intstk[nstrgt];
                  } else if (indlin == 2) {
                     logstk[nstrgt-1] = relstk[nstrgt-1] > relstk[nstrgt];
                  } else if (indlin == 3) {
                     logstk[nstrgt-1] = strcmp(chrstk[nstrgt-1], chrstk[nstrgt]) > 0;
                  } else if (indlin == 4) {
                     logstk[nstrgt-1] = dblstk[nstrgt-1] > dblstk[nstrgt];
                  } else {
                     indrgt[nstrgt-1] = -5;
                  }
               } else if (strcmp(cparav, "=") == 0) {
                  --nstrgt;
                  indlin = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  indrgt[nstrgt-1] = 5;
                  if (indlin == 1) {
                     logstk[nstrgt-1] = intstk[nstrgt-1] == intstk[nstrgt];
                  } else if (indlin == 2) {
                     logstk[nstrgt-1] = relstk[nstrgt-1] == relstk[nstrgt];
                  } else if (indlin == 3) {
                     logstk[nstrgt-1] = strcmp(chrstk[nstrgt-1], chrstk[nstrgt]) == 0;
                  } else if (indlin == 4) {
                     logstk[nstrgt-1] = dblstk[nstrgt-1] == dblstk[nstrgt];
                  } else {
                     indrgt[nstrgt-1] = -5;
                  }
               } else if (strcmp(cparav, "<=") == 0) {
                  --nstrgt;
                  indlin = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  indrgt[nstrgt-1] = 5;
                  if (indlin == 1) {
                     logstk[nstrgt-1] = intstk[nstrgt-1] <= intstk[nstrgt];
                  } else if (indlin == 2) {
                     logstk[nstrgt-1] = relstk[nstrgt-1] <= relstk[nstrgt];
                  } else if (indlin == 3) {
                     logstk[nstrgt-1] = strcmp(chrstk[nstrgt-1], chrstk[nstrgt]) <= 0;
                  } else if (indlin == 4) {
                     logstk[nstrgt-1] = dblstk[nstrgt-1] <= dblstk[nstrgt];
                  } else {
                     indrgt[nstrgt-1] = -5;
                  }
               } else if (strcmp(cparav, ">=") == 0) {
                  --nstrgt;
                  indlin = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  indrgt[nstrgt-1] = 5;
                  if (indlin == 1) {
                     logstk[nstrgt-1] = intstk[nstrgt-1] >= intstk[nstrgt];
                  } else if (indlin == 2) {
                     logstk[nstrgt-1] = relstk[nstrgt-1] >= relstk[nstrgt];
                  } else if (indlin == 3) {
                     logstk[nstrgt-1] = strcmp(chrstk[nstrgt-1], chrstk[nstrgt]) >= 0;
                  } else if (indlin == 4) {
                     logstk[nstrgt-1] = dblstk[nstrgt-1] >= dblstk[nstrgt];
                  } else {
                     indrgt[nstrgt-1] = -5;
                  }
               } else if (strcmp(cparav, "<>") == 0) {
                  --nstrgt;
                  indlin = min(indrgt[nstrgt-1], indrgt[nstrgt]);
                  indrgt[nstrgt-1] = 5;
                  if (indlin == 1) {
                     logstk[nstrgt-1] = intstk[nstrgt-1] != intstk[nstrgt];
                  } else if (indlin == 2) {
                     logstk[nstrgt-1] = relstk[nstrgt-1] != relstk[nstrgt];
                  } else if (indlin == 3) {
                     logstk[nstrgt-1] = strcmp(chrstk[nstrgt-1], chrstk[nstrgt]) != 0;
                  } else if (indlin == 4) {
                     logstk[nstrgt-1] = dblstk[nstrgt-1] != dblstk[nstrgt];
                  } else {
                     indrgt[nstrgt-1] = -5;
                  }
               } else {
/*                NO CHANCE, MAN... TRY IT WITH VARIABLES */
                  int_32 ilowrc = ivabeg;
                  int_32 ihigrc = ivaend;
L50:
                  imedrc = (ihigrc + ilowrc) / 2;
                  iofset = (imedrc-1) * lrclen;
                  iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
                  if (iretcd != 0) goto L9023;
                  strcpy(cparin, record2.cparin);
                  strcpy(cdatin, record2.cdatin);
                  indlin = record2.indlin;
                  idatin = record2.idatin;
                  adatin = record2.adatin;
                  ddatin = record2.ddatin;
                  if (strcmp(cparin, cparav) == 0) {
/*                   STACK VARIABLE VALUE */
                     ++nstrgt;
                     indrgt[nstrgt-1] = indlin;
                     if (indlin == 1) {
                        intstk[nstrgt-1] = idatin;
                     } else if (indlin == 2) {
                        relstk[nstrgt-1] = adatin;
                     } else if (indlin == 3) {
                        strcpy(chrstk[nstrgt-1], cdatin);
                        intstk[nstrgt-1] = idatin;
                     } else if (indlin == 4) {
                        dblstk[nstrgt-1] = ddatin;
                     } else if (indlin == 5) {
                        logstk[nstrgt-1] = idatin == 1;
                     }
                  } else if (strcmp(cparin, cparav) < 0) {
                     ilowrc = imedrc;
                     goto L50;
                  } else {
                     ihigrc = imedrc;
                     goto L50;
                  }
               }
/* L05>        ENDIF( ON CPARAV ) */
            }

/* L04>  ELSEIF( INDLEC.NE.3 )THEN */
         } else {
            ++nstrgt;
            indrgt[nstrgt-1] = indlec;
            if (indlec == 1) {
               intstk[nstrgt-1] = idatin;
            } else if (indlec == 2) {
               relstk[nstrgt-1] = adatin;
            } else if (indlec == 4) {
               dblstk[nstrgt-1] = ddatin;
            }

/* L04>  ENDIF( ON INDLEC ) */
         }

/* L03> ELSEIF( .NOT.LRGTST )THEN */
      } else {
         strcpy(cparav, cdatav);
         if (strcmp(cparav, chrend) == 0) {
            lrgtst = 0;
         } else if (ilogin <= 6 && strcmp(cparav, ":=") == 0) {
            lrgtst = 1;
         } else {
            int_32 ilowrc = ivabeg;
            int_32 ihigrc = ivaend;

            ++nstlft;

/*          FIND RECORD NUMBER FOR VARIABLE *CPARAV* */
L27:
            imedrc = (ihigrc + ilowrc) / 2;
            iofset = (imedrc-1) * lrclen;
            iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
            if (iretcd != 0) goto L9023;
            strcpy(cparin, record2.cparin);
            strcpy(cdatin, record2.cdatin);
            indlin = record2.indlin;
            idatin = record2.idatin;
            adatin = record2.adatin;
            ddatin = record2.ddatin;
            if (strcmp(cparin, cparav) == 0) {
               irclft[nstlft-1] = imedrc;
            } else if (strcmp(cparin, cparav) < 0) {
               ilowrc = imedrc;
               goto L27;
            } else {
               ihigrc = imedrc;
               goto L27;
            }
         }
/* L03>  ENDIF( ON LRGTST ) */
      }
/* L02>  PROCESS *CLE2000* STATEMENTS: AND HIT... (END) */
   }
/* L01> ***LOOP*** OVER WORDS (END) */
   goto L10;

L8000:
   printf("%s: IREC=%d REC=%s\n", nomsub, (int)irecor, myreco);
   xabort_c("REDGET: PROBLEM READING *RECORD*, USE DEBUG");
L9001:
   xabort_c("REDGET: EMBEDDED LOGIC LEVEL .LT.  1 (MIN)");
L9002:
   xabort_c("REDGET: EMBEDDED LOGIC LEVEL .GT. 128 (MAX)");
L9004:
   xabort_c("REDGET: NUMBER OF >>.<< ACCUMULATED.GT. 128 (MAX)");
L9006:
   printf("%s: UNDEFINED LOGICAL IN *%s* OPERATION\n", nomsub, clogbg[ilogin-1]);
   xabort_c("REDGET: IMPOSSIBLE LOGICAL OPERATION");
L9007:
   xabort_c("REDGET: IMPOSSIBLE TO PRINT VALUE");
L9008:
   printf("%s: VARIABLE *%s* HAS STILL NO VALUE\n", nomsub, cparav);
   xabort_c("REDGET: IMPOSSIBLE TO GET VALUE");
L9009:
   printf("%s: *%s* HAS NEGATIVE VALUE\n", nomsub, clogbg[idefkw-1]);
   xabort_c("REDGET: IMPOSSIBLE TO TAKE *SQRT*");
L9010:
   printf("%s: *%s* DIVISION BY 0\n", nomsub, clogbg[idefkw-1]);
   xabort_c("REDGET: IMPOSSIBLE TO DIVIDE");
L9011:
   printf("%s: *%s* .LT. 0 RAISED TO POWER .LT. 0\n", nomsub, clogbg[idefkw-1]);
   xabort_c("REDGET: IMPOSSIBLE TO TAKE POWER");
L9012:
   xabort_c("REDGET: IMPOSSIBLE TO + OR - STRINGS");
L9013:
   xabort_c("REDGET: LONG < 99999999 REQUIRED FOR CONVERSION TO STRING");
L9014:
   xabort_c("REDGET: LONG > -9999999 REQUIRED FOR CONVERSION TO STRING");
   goto L10;
L9023:
   printf("%s: IOSTAT RETURN CODE =%d\n", nomsub, (int)iretcd);
   xabort_c("REDGET: IMPOSSIBLE TO USE THIS *STACK* FILE");
} /* redget_c */

/*          *REDPUT* ENTRY POINT */
/*                   TO INPUT VALUES FOR CLE-2000 VARIABLES */
/*                   USING THE COMMAND >>.<< */
/*     INPUT VARIABLES: */
/*          *ITYP*   TYPE FOR VARIABLE (+1: INT */
/*                                      +2: REAL */
/*                                      +3: STRING */
/*                                      +4: DOUBLE */
/*                                      +5: LOGICAL ) */
/*          *NITMA*  INT VALUE IF *ITYP*.EQ.+1.OR.*ITYP*.EQ.+5 */
/*          *FLOTT*  REAL    VALUE IF *ITYP*.EQ.+2 */
/*          *TEXT*   STRING  VALUE IF *ITYP*.EQ.+3 */
/*          *DFLOT*  DOUBLE  VALUE IF *ITYP*.EQ.+4 */

/*     NOTE: LOGICAL VALUES ARE GIVEN EITHER BY *TRUE* = *NITMA*.EQ.+1 */
/*                                        OR BY *FALSE*= *NITMA*.EQ.-1 */

void redput_c(int_32 *ityp, int_32 *nitma, float_32 *flott, char *text, double_64 *dflot)
{
   char *nomsub="redput_c";
   char cparin[13], cdatin[73];
   int_32 iretcd, ilengt, indlin, idatin;
   float_32 adatin;
   double_64 ddatin;
   int_32 imedrc, iofset;

   if (*ityp == 3) {
      ilengt = strlen(text);
   } else {
      ilengt = 0;
   }
   if (nstput == 0) {
      xabort_c("REDPUT: NOTHING TO PUT");
   } else if (ilengt > 72) {
      xabort_c("REDPUT: STRING LENGTH RESTRICTED TO 72");
   } else if (*ityp <= 0) {
      xabort_c("REDPUT: PLEASE USE *ITYP*.GT.0");
   } else if (irecor == 0 || irecor > ninput) {
      xabort_c("REDPUT: READER IS CLOSED OR FILE END");
   }
   imedrc = ircput[nstput-1];
   iofset = (imedrc-1) * lrclen;
   iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
   if (iretcd != 0) goto L9023;
   strcpy(cparin, record2.cparin);
   strcpy(cdatin, record2.cdatin);
   indlin = record2.indlin;
   idatin = record2.idatin;
   adatin = record2.adatin;
   ddatin = record2.ddatin;
   if (indlin >= 0) xabort_c("REDPUT: CANNOT PUT ON A DEFINED VALUE");
   indlin = abs(indlin);
   if (indlin != *ityp) xabort_c("REDPUT: INCOMPATIBLE TYPE OF THE VARIABLE");
   if (indlin == 1) {
      idatin = *nitma;
   } else if (indlin == 2) {
      adatin = *flott;
   } else if (indlin == 3) {
      if (strncmp(text, " ", ilengt) == 0) {
/*       ALL BLANK STRING IS CONSIDERED AS NULL-STRING => "" */
         idatin = 0;
      } else if (text[0] == '"') {
/*       THIS IS A "...." STRING =>  START "... */
         if (ilengt == 1) {
/*          PROVIDES A WAY FOR APPLICATION TO PUT A '"' STRING */
            idatin = 1;
            cdatin[0] = text[0];
         } else {
/*          LOOK FOR => END ..." */
            idatin = index_f(&text[1], "\"") - 1;
            if (idatin < 0) {
               if (strcmp(&text[1], " ") == 0) {
                  idatin = 1;
                  cdatin[0] = text[0];
               } else {
                  xabort_c("REDPUT: INVALID STRING \" NEVER ENDS)");
               }
            } else if (idatin != 0) {
               if (ilengt == idatin + 2) {
                  strncpy(cdatin, &text[1], idatin); cdatin[idatin] = '\0';
               } else {
                  if (strcmp(&text[idatin+2], " ") != 0) xabort_c("REDPUT: \".\" + OTHER WORDS");
                  strncpy(cdatin, &text[1], idatin); cdatin[idatin] = '\0';
               }
            }
         }
      } else {
         strncpy(cdatin, text, ilengt); cdatin[ilengt] = '\0';
         idatin = ilengt;
      }
   } else if (indlin == 4) {
      ddatin = *dflot;
   } else if (indlin == 5) {
      idatin = *nitma;
      if (idatin != -1 && idatin != 1) xabort_c("REDPUT: LOGICAL IS UNDEFINED");
   } else {
      xabort_c("REDPUT: UNDEFINED TYPE");
   }
   strcpy(record2.cparin, cparin);
   strcpy(record2.cdatin, cdatin);
   record2.indlin = indlin;
   record2.idatin = idatin;
   record2.adatin = adatin;
   record2.ddatin = ddatin;
   iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
   if (iretcd != 0) goto L9023;

/* ONE LESS TO PUT */
   --nstput;
   return;
L9023:
   printf("%s: IOSTAT RETURN CODE =%d\n", nomsub, (int)iretcd);
   xabort_c("REDPUT: IMPOSSIBLE TO USE THIS *STACK* FILE");
} /* redput_c */

/*          *REDOPN* ENTRY POINT */
/*                   TO OPEN THE DA-FILE CONTAINING CLE-2000 DATA */
/*     INPUT VARIABLES: */
/*          *IINP1* IS THE CLE-2000 DA-FILE UNIT */
/*          *IOUT1* IS THE OUTPUT FILE UNIT FOR MESSAGES (NORMALLY STDOUT) */
/*          *FILENAME* IS THE OUTPUT FILE NAME FOR MESSAGES */
/*          *NREC*  IS THE RECORD NUMBER WHERE WE START READING */

void redopn_c(kdi_file *iinp1, FILE *iout1, char *filename, int_32 nrec)
{
   char *nomsub="redopn_c";
   int_32 nstack;
   char cparav[13];
   int_32 iretcd;

   iunito = iinp1;
   iwrite = iout1;
   strcpy(hwrite, filename);

/* READ TOP OF OBJECT FILE */
   irecor = 1;
   iretcd = kdiget_c(iunito, (int_32 *)&header, 0, kdisize(header));
   if (iretcd != 0) goto L9023;
   strcpy(cparav, header.cparin);
   strcpy(myreco, header.cdatin);
   ninput = header.ninput;
   nstack = header.nstack;
   ioulst = header.ioulst;
   idblst = header.idblst;
   if (strcmp(cparav, cl2000) != 0) goto L9025;
   if (nrec != 0) irecor = nrec;
   nwords = 0;
   iwords = 0;
   ivabeg = ninput;
   ivaend = ninput + nstack + 1;
   return;
L9023:
   printf("%s: IOSTAT RETURN CODE =%d\n", nomsub, (int)iretcd);
   xabort_c("REDOPN: IMPOSSIBLE TO USE THIS *STACK* FILE");
L9025:
   xabort_c("REDOPN: UNABLE TO OPEN FILE");
} /* redopn_c */

/*          *REDCLS* ENTRY POINT */
/*                   TO CLOSE THE DA-FILE AT A CURRENT RECORD POSITION */
/*     OUTPUT VARIABLES: */
/*          *IINP1* IS THE CLE-2000 DA-FILE UNIT */
/*          *IOUT1* IS THE OUTPUT FILE UNIT FOR MESSAGES (NORMALLY STDOUT) */
/*          *FILENAME* IS THE OUTPUT FILE NAME FOR MESSAGES */
/*          *NREC*  IS THE RECORD NUMBER WHERE WE STOP READING */

void redcls_c(kdi_file **iinp1, FILE **iout1, char filename[73], int_32 *nrec)
{
   if (nwords != iwords) xabort_c("REDCLS: RECORD NOT FINISHED => CANNOT CLOSE");
   *nrec = irecor;
   *iinp1 = iunito;
   *iout1 = iwrite;
   strcpy(filename, hwrite);

/* WE PUT IRECOR=0 TO CLOSE THE READER (SEE START OF *REDGET*) */
   irecor = 0;
   ninput = 0;
   nwords = 0;
   iwords = 0;
   return;
} /* redcls_c */
