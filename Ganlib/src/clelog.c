
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 09/05/09      */
/*****************************************/

#include <string.h>
#include "cle2000.h"
#include "header.h"
#define min(A,B)  ((A) < (B) ? (A) : (B))
#define max(A,B)  ((A) > (B) ? (A) : (B))
#define index_f(A, B) (strstr(A, B) == NULL ? 0 : strstr(A, B) - A + 1)
#define ndecal 4
#define nmawrd 36
#define ndimst 128
static char AbortString[132];

int_32 clelog(FILE *iredin, FILE *iwrite, kdi_file *iunito)
{
   char *nomsub="clelog";
   static char cl2000[] = "CLE2000(V3)";
   static int_32 lvelbg[] = {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, -1, 0, -1, -1};
   static char *clognd[] = {";", ";", ";", ";", ";", ";", ";", "THEN", "THEN", "DO",
                            ";", ";", "REPEAT", "ELSE", ";"};
   static int_32 lvelnd[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 1, 0, -1};
   static char ctitle[] = "* CLE-2000 VERS 3.0 * R.ROY, EPM COPYRIGHT 1999 *";
   static char *terror[] = {"! CLELOG: UNEXPECTED CHARACTERS TO BE REPLACED WITH BLANKS",
                            "! CLELOG: UNBALANCED OPENING OR CLOSING STRINGS",
                            "! CLELOG: MISPLACED SEMICOLON ...; OR ;... OR ...;...",
                            "! CLELOG: CHARACTERS SUPPRESSED OUTSIDE COLUMN RANGE 1:72",
                            "! CLELOG: << AND >> NOT ALLOWED IN STRINGS (SUPPRESSED)",
                            "! CLELOG: (* ... *) INVALID COMMENTS (USE ! INSTEAD)",
                            "! CLELOG: QUIT \"...\" . SHOULD APPEAR ALONE A SINGLE LINE"};
   static char csemic[] = ";";
   static char digped[] = "0123456789+-.ED";
   static char onelet[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+-*/%<>=;";
   static char *clogbg[] = {"INTEGER", "REAL", "STRING", "DOUBLE", "LOGICAL", "EVALUATE", "ECHO",
                            "ELSEIF", "IF", "WHILE", "UNTIL", "ENDWHILE", "REPEAT", "ELSE", "ENDIF"};

/*     CLE-2000 SYSTEM: R.ROY (11/1999), VERSION 3.0 */

/*             *CLELOG* FIRST-PASS COMPILE OF THE INPUT UNIT  *IREDIN* */
/*                      COMPILER COMMENTS ARE WRITTEN ON UNIT *IWRITE* */
/*                      RESULT IS THE OBJECT D.A. UNIT        *IUNITO* */

/*        USE:          INPUT DATA IS COPIED ON D.A. UNIT, */
/*                      SENTENCES AND LOGIC DESCRIPTORS ARE SEPARATED, */
/*                      LOGICAL LEVELS ARE BUILT AND LOGIC IS CHECKED. */

/*      INPUT: *IREDIN* IS THE INPUT  UNIT */
/*             *IWRITE* IS THE OUTPUT UNIT FOR COMPILER COMMENTS */
/*             *IUNITO* IS THE DIRECT ACCESS UNIT FOR OBJECT CODE */

/*       NOTE: *CLELOG* = 0 IF NO PROBLEM WAS ENCOUNTERED WHILE COMPILING */

   int_32 ret_val = 0, lapos1 = 0, lapos2 = 0;
   int_32 ilines = 0, idblst = -1, ixrlst = -1, ioulst = -1, lnwsen = 1, lrecio = 1, l1lett = 1,
          nequal = 0, nwrsen = -2, nrecio = 0, nlevel = 1, maxlvl = 0, irecor = 0, nstlvl = 0,
          ninput = 1, nstack = 0, nrecor = 1;
   int_32 i, iloop1, jloop2, iretcd, iofset, ilevel, ilogin=0;
   int_32 maskck[nmaskc], ipacki[nmaskc];
   int_32 idebwd[nmawrd], ifinwd[nmawrd], jndlec[nmawrd];
   int_32 irclvl[ndimst], itylvl[ndimst];
   char crecbg[13], chrend[13], cerror[73];
   char recred[133+ndecal], cemask[73], myreco[73], cbla72[72], rwrite[73];
   int_32 jcm0bg, jst2bg, jfndbg;
   int_32 jst1bg, jfndnd, jcm0nd;
   int_32 ilengv=0, jsc0bg;
   int_32 jkthen, jkelse, jkrepe, jkdodo, jquitp;
   int_32 lnbprv, nwords, iapo12=0;
   int_32 iwords;
   char crecnd[13];

   fprintf(iwrite, "%-72s   LINE\n", ctitle);
   for ( i = 0; i < 72; i++) cbla72[i] = ' ';
   strcpy(crecbg, " ");
   strcpy(chrend, csemic);

/* WRITE TOP OF OBJECT FILE */
   strcpy(header.cparin, cl2000);
   strcpy(header.cdatin, " ");
   header.nrecor = nrecor;
   header.ninput = ninput;
   header.maxlvl = maxlvl;
   header.nstack = nstack;
   header.ixrlst = ixrlst;
   header.ioulst = ioulst;
   header.idblst = idblst;
   header.nobjet = 0;
   iretcd = kdiput_c(iunito, (int_32 *)&header, 0, kdisize(header));
   if (iretcd != 0) goto L9001;

/* DUMP INPUT FILE TO OBJECT FILE */
   strncpy(cemask, cbla72, 72); cemask[72] = '\0';
   for ( i = 0; i < ndecal+132; i++) recred[i] = ' ';
   recred[ndecal+132] = '\0';
L10:
   ++ilines;
/* READ A NEW RECORD */
   iretcd = fscanf(iredin, "%133[^\n]\n", &recred[ndecal]);
   if (iretcd == EOF) goto L100;
   if (strlen(&recred[ndecal]) > 132) goto L5000;
   for(i = strlen(&recred[ndecal]); i < 132; i++) recred[ndecal+i] = ' ';
   recred[ndecal+132] = '\0';
   strncpy(myreco, &recred[ndecal], 72); myreco[72] = '\0';

/* OUPUT THE RECORD */
   fprintf(iwrite, "%-72s   %04d\n", myreco, (int)ilines);

/* SUPPRESS RECORD IF ALL COMMENTS OR BLANK */
   if (strncmp(myreco, cbla72, 72) == 0) {
     if (ilines == 1) {
       sprintf(AbortString,"%s: empty line %d",nomsub,ilines);
       xabort_c(AbortString);
     }
     goto L10;
   }

/* SUPPRESS ! EXCLAMATION COMMENTS FROM *RECRED* AND RECORD */
   jcm0bg = index_f(recred, "!");
   if (jcm0bg != 0) {
      for ( i = jcm0bg - 1; i < ndimst; i++) recred[i]=' ';
      if (jcm0bg + ndecal - 1 < 72) {
         for ( i = jcm0bg - ndecal - 1; i < 72; i++) myreco[i]=' ';
      }
   }
   if (myreco[0] == '*') goto L10;

/* ANYTHING OUTSIDE COLUMNS NDECAL+1:NDECAL+72 RANGE IN *RECRED* */
   if (strncmp(&recred[ndecal+72], "       ", 7) != 0) {
      printf("%.132s\n", &recred[ndecal]);
      printf("%72s????????\n", " ");
      printf("%72s\n", terror[3]);
      ++ret_val;
   }
   for (iloop1 = 0; iloop1 < 72; ++iloop1) {
/*    SUPPRESS UNEXPECTED CHARACTERS */
      if (strncmp(&myreco[iloop1], " ", 1) < 0) {
         cemask[iloop1] = '?';
         myreco[iloop1] = ' ';
         recred[iloop1+ndecal] = ' ';
         ++ret_val;
      }
   }

/* SUPPRESS STRINGS OF TYPES: '...' OR "...." FROM *RECRED* */
L25:
   jst1bg = index_f(recred, " \'");
   if (jst1bg == 0) jst1bg = 132;
   jst2bg = index_f(recred, " \"");
   if (jst2bg == 0) jst2bg = 132;
   jfndbg = min(jst1bg,jst2bg);
   if (jfndbg != 132) {
      if (jfndbg == jst1bg) {
         jfndnd = index_f(&recred[jfndbg+1], "\' ") + jfndbg + 1;
      } else {
         jfndnd = index_f(&recred[jfndbg+1], "\" ") + jfndbg + 1;
      }
      if (jfndnd == jfndbg + 1) {
         cemask[jfndbg-ndecal] = '?';
         strcpy(cerror, terror[1]);
         ++ret_val;
         jfndnd = 132;
      }

/*    STRING IS FOUND, CHECK IF << OR >> IS CONTAINED INSIDE */
      strncpy(rwrite, &recred[jfndbg], jfndnd-jfndbg); rwrite[jfndnd-jfndbg] = '\0';
      jcm0bg = index_f(rwrite, "<<") + jfndbg;
      jcm0nd = index_f(rwrite, ">>") + jfndbg;
      if (jcm0bg != jfndbg) {
         memcpy(&cemask[jcm0bg-ndecal-1], "??", 2);
         memcpy(&myreco[jcm0bg-ndecal-1], "  ", 2);
         strcpy(cerror, terror[4]);
         ++ret_val;
      }
      if (jcm0nd != jfndbg) {
         memcpy(&cemask[jcm0nd-ndecal-1], "??", 2);
         memcpy(&myreco[jcm0nd-ndecal-1], "  ", 2);
         strcpy(cerror, terror[4]);
         ++ret_val;
      }
      for ( i = jfndbg; i < jfndnd; i++) recred[i] = ' ';
      goto L25;
   }

/* CONTROL STRANGE FORMS OF TYPES: ...'... OR ..."... */
   jst1bg = index_f(recred, "'");
   if (jst1bg != 0) {
      cemask[jst1bg-ndecal-1] = '?';
      strcpy(cerror, terror[1]);
      ++ret_val;
      recred[jst1bg-1] = ' ';
      myreco[jst1bg-ndecal-1] = ' ';
   }

   if (strncmp(cemask, cbla72, 72) != 0) {
      fprintf(iwrite, "%-72s\n", cemask);
      fprintf(iwrite, "%-72s\n", cerror);
      strncpy(cemask, cbla72, 72);
   }

/* SUPPRESS OLD FORMS OF COMMENTS: (*...  ...*) */
L26:
   jcm0bg = index_f(recred, "(*");
   jcm0nd = index_f(recred, "*)");
   if (jcm0bg != 0) {
      strcpy(cerror, "! WARNING: (* ... *) OBSOLETE COMMENTS (USE ! INSTEAD)");
      if (jcm0nd == 0) {
         strncpy(&myreco[jcm0bg-ndecal-1], cbla72, 72-(jcm0bg-ndecal-1));
         strncpy(&recred[jcm0bg-1], cbla72, 132-(jcm0bg-1));
         strcpy(cerror, terror[5]);
         ++ret_val;
      } else {
         strncpy(&myreco[jcm0bg-ndecal-1], cbla72, jcm0nd + 1 - (jcm0bg - 1));
         strncpy(&recred[jcm0bg-1], cbla72, jcm0nd + 1 - (jcm0bg - 1));
      }
      goto L26;
   } else if (jcm0nd != 0) {
      memcpy(&cemask[jcm0nd-ndecal-1], "??", 2);
      strcpy(cerror, terror[5]);
      ++ret_val;
      strncpy(recred, cbla72, jcm0nd + 1);
      strncpy(myreco, cbla72, jcm0nd - ndecal + 1);
      goto L26;
   }

   if (strncmp(cemask, cbla72, 72) != 0) {
      fprintf(iwrite, "%-72s\n", cemask);
      fprintf(iwrite, "%-72s\n", cerror);
      strncpy(cemask, cbla72, 72);
   }

/* TO SEPARATE LOGIC, PUT RETURNS FOR INPUT LINES ENDING WITH: */
/* *;*, *REPEAT*, *THEN*, *ELSE* OR *DO* */
L30:
   jsc0bg = index_f(recred, " ; ");
   if (jsc0bg == 0) {
/*    CONTROL STRANGE FORMS OF TYPES: ...; OR ;... OR ...;... */
      jsc0bg = index_f(recred, ";");
      if (jsc0bg != 0) {
         cemask[jsc0bg-ndecal-1] = '?';
         strcpy(cerror, terror[2]);
         ++ret_val;
         recred[jsc0bg-1] = ' ';
         myreco[jsc0bg-ndecal-1] = ' ';
      }
      jfndbg = 132;
   } else {
      jfndbg = jsc0bg;
      ilengv = 1;
   }
   jkthen = index_f(recred, " THEN ");
   if (jkthen == 0) jkthen = 132;
   if (jkthen < jfndbg) {
      jfndbg = jkthen;
      ilengv = 4;
   }
   jkelse = index_f(recred, " ELSE ");
   if (jkelse == 0) jkelse = 132;
   if (jkelse < jfndbg) {
      jfndbg = jkelse;
      ilengv = 4;
   }
   jkrepe = index_f(recred, " REPEAT ");
   if (jkrepe == 0) jkrepe = 132;
   if (jkrepe < jfndbg) {
      jfndbg = jkrepe;
      ilengv = 6;
   }
   jkdodo = index_f(recred, " DO ");
   if (jkdodo == 0) jkdodo = 132;
   if (jkdodo < jfndbg) {
      jfndbg = jkdodo;
      ilengv = 2;
   }
   jquitp = index_f(recred, " QUIT ");
   if (jquitp == 0) jquitp = 132;
   if (jquitp < jfndbg) {
      jfndbg = jquitp;
      ilengv = 4;
   }

   if (jfndbg == 132) {
      strncpy(rwrite, myreco, 72);
      strncpy(myreco, cbla72, 72);
   } else {
      if (jfndbg == jquitp) {
         strncpy(myreco, cbla72, jfndbg - ndecal + ilengv);
         goto L200;
      } else {
         strncpy(rwrite, cbla72, 72);
         strncpy(rwrite, myreco, jfndbg - ndecal + ilengv);
         strncpy(myreco, cbla72, jfndbg - ndecal + ilengv);
         strncpy(&recred[jfndbg], cbla72, jfndbg + ilengv - jfndbg);
      }
   }
   if (strncmp(cemask, cbla72, 72) != 0) {
      fprintf(iwrite, "%-72s\n", cemask);
      fprintf(iwrite, "%-72s\n", cerror);
      strncpy(cemask, cbla72, 72);
      goto L30;
   }

/* SUPPRESS RECORD IF ALL IS STILL BLANK */
   if (strncmp(rwrite, cbla72, 72) == 0) goto L10;

/* NEW RECORD FOUND, READY TO PROCESS *RWRITE* */
   for (iloop1 = 0; iloop1 < nmaskc; ++iloop1) {
      maskck[iloop1] = 0;
      ipacki[iloop1] = 0;
   }

/* *** BEWARE **** FROM HERE WORDS ARE IN REVERSE ORDER... */

/* BEGIN: CONSTRUCT MASK NUMBERS */

/* PREVIOUS NON-BLANK CHARACTER (ASSUME BLANK AT START) */
   lnbprv = 0;
   nwords = 0;
   for (iloop1 = 72; iloop1 >= 1; --iloop1) {
      jloop2 = (iloop1 + 23) / 24;
      maskck[jloop2-1] <<= 1;
      if (lapos1) {
/*       ALL CHARACTERS ARE MASKED WHILE *LAPOS1* */
         ++maskck[jloop2-1];
         lapos1 = rwrite[iloop1-1] != '\'';

/*       MAKE AS IF PREVIOUS WAS BLANK */
         lnbprv = 0;
         --idebwd[nwords-1];
      } else if (lapos2) {
/*       ALL CHARACTERS ARE MASKED WHILE *LAPOS1* */
         ++maskck[jloop2-1];
         lapos2 = rwrite[iloop1-1] != '"';

/*       MAKE AS IF PREVIOUS WAS BLANK */
         lnbprv = 0;
         --idebwd[nwords-1];
      } else {
         int_32 lnbcur = rwrite[iloop1-1] != ' ';
         if (lnbcur) {
/*          FIND A NON-BLANK CHARACTER, MASK IT */
            ++maskck[jloop2-1];
            if (lnbprv) {
               --idebwd[nwords-1];
            } else {
/*             PREVIOUS ONE WAS BLANK, LOOK FOR ' OR " */
               lapos1 = rwrite[iloop1-1] == '\'';
               lapos2 = rwrite[iloop1-1] == '\"';
               if (lapos1 || lapos2) iapo12 = iloop1;
/*             BEGIN A NEW WORD (REVERSED ORDER) */
               ++nwords;
               ifinwd[nwords-1] = iloop1;
               idebwd[nwords-1] = iloop1;
            }
         } else if (lnbprv) {
/*          FIND A BLANK CHARACTER, BUT AFTER A NON-BLANK */
/*          THIS COULD BE A MISTAKE IF ' OR " ARE NOT IN USE. */
            if ((!lapos1 && rwrite[iloop1] == '\'') || (!lapos2 && rwrite[iloop1] == '\"')) {
               cemask[iloop1] = '?';
            }
         }
         lnbprv = lnbcur;
      }
   }
   if (lapos1 || lapos2) {
      cemask[iapo12-1] = '?';
      lapos1 = 0;
      lapos2 = 0;
   }
   if (strncmp(cemask, cbla72, 72) != 0) {
      fprintf(iwrite, "%-72s\n", cemask);
      fprintf(iwrite, "%-72s\n", terror[1]);
      strncpy(cemask, cbla72, 72);
   }
/* END:   CONSTRUCT MASK NUMBERS */

/* BEGIN: IDENTITY TYPES AND PACK *JNDLEC* WITH (ITYP-1) DATA */
   for (iwords = 1; iwords <= nwords; ++iwords) {
      char cdatin[73];
      int_32 jindex=0;
      ilengv = ifinwd[iwords-1] - idebwd[iwords-1] + 1;
      strncpy(cdatin, &rwrite[idebwd[iwords-1]-1], ilengv); cdatin[ilengv] = '\0';

/*    DETERMINATION OF TYPE FOR THAT WORD */
      if (cdatin[0] == '\'' || cdatin[0] == '"') {
         jndlec[iwords-1] = 2;
      } else if (ilengv == 1) {
         jindex = index_f(digped, &cdatin[0]);
         if (jindex > 0 && jindex <= 10) {
            jndlec[iwords-1] = 0;
         } else {
            jndlec[iwords-1] = 2;
            l1lett = l1lett && index_f(onelet, &cdatin[0]) != 0;
         }
      } else {
         int_32 ipoint = 0;
         int_32 ifloat = 0;
         int_32 idoubl = 0;
         for (iloop1 = 1; iloop1 <= ilengv; ++iloop1) {
            char cc[] = {cdatin[iloop1-1], '\0'};
            jindex = index_f(digped, cc);
            if (jindex == 0) {
               goto L62;
            } else if (jindex == 11 || jindex == 12) {
               if (iloop1 != 1) {
/*                CHECK SIGN AFTER EXPONENT */
                  if (iloop1 - 1 != ifloat && iloop1 - 1 != idoubl) {
                     jindex = 0;
                     goto L62;
                  }
               }
            } else if (jindex == 13) {
               if (ipoint != 0) {
                  jindex = 0;
                  goto L62;
               }
               ipoint = iloop1;
            } else if (jindex == 14) {
               if (ifloat != 0 || iloop1 == 1) {
                  jindex = 0;
                  goto L62;
               }
               ifloat = iloop1;
            } else if (jindex == 15) {
               if (idoubl != 0 || iloop1 == 1) {
                  jindex = 0;
                  goto L62;
               }
               idoubl = iloop1;
            }
         }
L62:
         if (jindex == 0) {
            jndlec[iwords-1] = 2;
         } else if (idoubl != 0) {
            jndlec[iwords-1] = 3;
         } else if (ifloat != 0 || ipoint != 0) {
            jndlec[iwords-1] = 1;
         } else {
            jndlec[iwords-1] = 0;
         }
      }
      jloop2 = (((nwords - iwords + 1) << 1) + 23) / 24;
      ipacki[jloop2-1] <<= 2;
      ipacki[jloop2-1] += jndlec[iwords-1];

/*    COUNT FOR THE NUMBER OF *:=*, *<<.>>* AND *>>.<<* */
/*    AND THE NUMBER OF *$.* BEFORE *:=* */
      if (jndlec[iwords-1] == 2) {
         if (ilengv == 2 && strcmp(cdatin, ":=") == 0) {
            ++nequal;
         } else if (ilengv >= 2) {
            if (strncmp(cdatin, ">>", 2) == 0) {
               lrecio = ilengv >= 5;
               if (lrecio) {
                  lrecio = strcmp(&cdatin[ilengv-2], "<<") == 0 && cdatin[2] != '$';
               }
               ++nrecio;
            } else if (strncmp(cdatin, "<<", 2) == 0) {
               lrecio = ilengv >= 5;
               if (lrecio) {
                  lrecio = strcmp(&cdatin[ilengv-2], ">>") == 0 && cdatin[2] != '$';
               }
               ++nrecio;
            }
         }
      }
   }
/* END:   IDENTITY TYPES AND PACK *JNDLEC* WITH (ITYP-1) DATA */

/* NOW, LOOK FOR 1-ST WORD OF SENTENCES */
   if (lnwsen && jndlec[nwords-1] == 2 && ifinwd[nwords-1] - idebwd[nwords-1] <= 11) {

/*    RECOVER 1-ST WORD TO CHECK BEGIN OF SENTENCE */
      char cparin[13];
      strncpy(cparin, &rwrite[idebwd[nwords-1]-1], ifinwd[nwords-1] - idebwd[nwords-1] + 1);
      cparin[ifinwd[nwords-1] - idebwd[nwords-1] + 1] = '\0';
      ilogin = 0;
      for (iloop1 = 0; iloop1 < 15; ++iloop1) {
         if (strcmp(cparin, clogbg[iloop1]) == 0) ilogin = iloop1 + 1;
      }

/*    BACKWARD COMPATIBILITY ( *CHARACTER* / *PRINT* ) */
      if (strcmp(cparin, "CHARACTER") == 0) {
         ilogin = 3;
         strcpy(cerror, "! WARNING: *CHARACTER* => *STRING* (REPLACED)");
         fprintf(iwrite, "%s\n", cerror);
      } else if (strcmp(cparin, "PRINT") == 0) {
         ilogin = 7;
         strcpy(cerror, "! WARNING: *PRINT*     => *ECHO*   (REPLACED)");
         fprintf(iwrite, "%s\n", cerror);
      }

      if (ilogin == 0) {
         strcpy(crecbg, " ");
         strcpy(chrend, csemic);
      } else {
         strcpy(crecbg, clogbg[ilogin-1]);
         strcpy(chrend, clognd[ilogin-1]);
/*       KEYWORDS: *IF+/+WHILE+/+REPEAT* */
         if (ilogin == 9 || ilogin == 10 || ilogin == 13) {
            if (nstlvl == ndimst) goto L7000;
            ++nstlvl;
            itylvl[nstlvl-1] = ilogin;
            irclvl[nstlvl-1] = ninput + 1;
/*       KEYWORDS: *ENDWHILE* */
         } else if (ilogin == 12) {
            if (itylvl[nstlvl-1] != 10) goto L7001;
/*       KEYWORDS: *UNTIL* */
         } else if (ilogin == 11) {
            if (itylvl[nstlvl-1] != 13) goto L7001;
/*       KEYWORD:  *ELSEIF+/+ELSE* */
         } else if (ilogin == 8 || ilogin == 14) {
            if (nstlvl == 0) goto L7000;
            if (itylvl[nstlvl-1] != 8 && itylvl[nstlvl-1] != 9) goto L7001;
            itylvl[nstlvl-1] = ilogin;
/*       KEYWORD:  *ENDIF* */
         } else if (ilogin == 15) {
            if (itylvl[nstlvl-1] != 8 && itylvl[nstlvl-1] != 9 && itylvl[nstlvl-1] != 14) goto L7001;
         }

/*       KEYWORDS:  *UNTIL+/+ENDWHILE+/+ENDIF* */
         if (ilogin == 11 || ilogin == 12 || ilogin == 15) {
            int_32 jrecor = ninput + 1;
            if (nstlvl == 0) goto L7002;
            irecor = irclvl[nstlvl-1];
            --nstlvl;

/*          REWRITE OLD RECORD TO KEEP LINK WITH THIS END */
            iofset = (irecor - 1) * lrclen;
            iretcd = kdiget_c(iunito, (int_32 *)&record1, iofset, kdisize(record1));
            if (iretcd != 0) goto L9003;
            record1.irecor = jrecor;
            iretcd = kdiput_c(iunito, (int_32 *)&record1, iofset, kdisize(record1));
            if (iretcd != 0) goto L9001;
         }
      }
   }

   if (ilogin == 0) {
      ilevel = 0;
/*    STATEMENTS OUTSIDE CLE-2000, STRINGS MUST BE '...' */
      for (iwords = 1; iwords <= nwords; ++iwords) {
         if (jndlec[iwords-1] == 2 && rwrite[idebwd[iwords-1]-1] == '\"') {
            cemask[idebwd[iwords-1]-1] = '?';
            cemask[ifinwd[iwords-1]-1] = '?';
            fprintf(iwrite, "%-72s\n", cemask);
            strcpy(cerror, "! WARNING: OUTSIDE CLE-2000, ENCLOSE STRINGS IN '...' (REPLACED)");
            strncpy(cemask, cbla72, 72);
            rwrite[idebwd[iwords-1]-1] = '\'';
            rwrite[ifinwd[iwords-1]-1] = '\'';
         }
      }
   } else {
/*      FOR CLE-2000 STATEMENTS, */
/*      AND INVALID *ILENGV* > 12 FOR STRINGS NOT ENCLOSED BY "..." */
      for (iwords = 1; iwords <= nwords; ++iwords) {
         if (jndlec[iwords-1] == 2 && rwrite[idebwd[iwords-1]-1] != '\"') {
            if (rwrite[idebwd[iwords-1]-1] == '\'') {
               cemask[idebwd[iwords-1]-1] = '?';
               cemask[ifinwd[iwords-1]-1] = '?';
               fprintf(iwrite, "%-72s\n", cemask);
               strcpy(cerror, "! WARNING: INSIDE CLE-2000, ENCLOSE STRINGS IN \"...\" (REPLACED)");
               strncpy(cemask, cbla72, 72);
               rwrite[idebwd[iwords-1]-1] = '\"';
               rwrite[ifinwd[iwords-1]-1] = '\"';
            } else {
               if (ifinwd[iwords-1] - idebwd[iwords-1] > 11) goto L5012;
            }
         }
      }
      ilevel = nlevel + lvelbg[ilogin-1];
      maxlvl = max(maxlvl,ilevel);
      nwrsen += nwords;
   }

/* RECOVER LAST WORD TO CHECK END OF SENTENCE */
   if (jndlec[0] == 2 && ifinwd[0] - idebwd[0] <= 11) {
      strncpy(crecnd, &rwrite[idebwd[0]-1], ifinwd[0] - idebwd[0] + 1);
      crecnd[ifinwd[0] - idebwd[0] + 1] = '\0';
   } else {
      strcpy(crecnd, " ");
   }

/* WRITE RECORD OR PART OF IT */
   ++ninput;
   iofset = (ninput - 1) * lrclen;
   strcpy(record1.cparin, crecbg);
   rwrite[72] = '\0';
   strcpy(record1.myreco, rwrite);
   record1.ilines = ilines;
   record1.ilevel = ilevel;
   record1.irecor = irecor;
   for (i = 0; i < nmaskc; i++) record1.maskck[i] = maskck[i];
   for (i = 0; i < nmaskc; i++) record1.ipacki[i] = ipacki[i];
   iretcd = kdiput_c(iunito, (int_32 *)&record1, iofset, kdisize(record1));
   if (iretcd != 0) goto L9001;
   if (!lrecio) goto L7005;
   irecor = 0;
   strcpy(crecbg, " ");
   lnwsen = (strcmp(crecnd, chrend) == 0);
   if (lnwsen) {
      if (ilogin != 0) {
         nlevel += lvelnd[ilogin-1];
         if (nrecio != 0) goto L7005;
         if (!l1lett) goto L5001;
/*       KEYWORDS: *INTEGER+/+REAL+/+STRING+/+DOUBLE+/+LOGICAL* */
         if (ilogin <= 5) {
            if (nequal == 0) {
               if (nwrsen <= 0) goto L7004;
            } else if (nequal == 1) {
               if (nwrsen <= 2) goto L7004;
            } else {
               goto L7003;
            }
            if (nlevel != 1) goto L7006;
/*       KEYWORD: *EVALUATE+/ */
         } else if (ilogin == 6) {
            if (nequal != 1) goto L7003;
            if (nwrsen <= 2) goto L7004;
/*       KEYWORDS: *ECHO+/+ELSEIF+/+IF+/+WHILE+/+UNTIL* */
         } else if (ilogin >= 7 && ilogin <= 11) {
            if (nequal != 0) goto L7003;
            if (nwrsen <= 0) goto L7004;
/*       KEYWORDS: *REPEAT+/+ELSE* */
         } else if (ilogin == 13 || ilogin == 14) {
            if (nequal != 0) goto L7003;
            if (nwrsen != -1) goto L7004;
/*       KEYWORDS: *ENDWHILE+/+ENDIF* */
         } else if (ilogin == 12 || ilogin == 15) {
            if (nequal != 0) goto L7003;
            if (nwrsen != 0) goto L7004;
         }
      } else {
/*       USE OF <<.>> OR >>.<<, BUT STILL NO CLE-2000 INSTRUCTION */
         if (maxlvl == 0 && nrecio != 0) goto L7005;
      }
/*    RESET NUMBER OF EQUALS, WORDS, <<.>> >>.<<, $. */
      nequal = 0;
      nwrsen = -2;
      nrecio = 0;
      l1lett = 1;
   }
   if (strncmp(myreco, cbla72, 72) != 0) goto L30;
   goto L10;
L100:
   memcpy(myreco, "QUIT .", 6);
   fprintf(iwrite, "%-72sIMPLICIT\n",myreco);
   memcpy(myreco, "     .", 6);
L200:
   if (nlevel != 1 || strcmp(chrend, csemic) != 0) goto L7007;
   fprintf(iwrite, " \n");
   strncpy(rwrite, cbla72, 72);
   jfndbg = index_f(myreco, "\"");
   if (jfndbg != 0) {
      jfndnd = index_f(&myreco[jfndbg], "\"") + jfndbg;
      if (jfndnd == jfndbg) {
         cemask[jfndbg-1] = '?';
         printf("%s\n", cemask);
         printf("%s\n", terror[6]);
         ++ret_val;
         strncpy(cemask, cbla72, 72);
      } else {
         strncpy(rwrite, &myreco[jfndbg], jfndnd - jfndbg - 1);
         strncpy(&myreco[jfndbg-1], cbla72, jfndnd - jfndbg + 1);
      }
   }
   if (index_f(rwrite, "NODEBUG") == 0) {
      if (index_f(rwrite, "DEBUG") != 0) idblst = 1;
   }
   if (index_f(rwrite, "NOXREF") == 0) {
      if (index_f(rwrite, "XREF") != 0) ixrlst = 1;
   }
   if (index_f(rwrite, "NOLIST") == 0) {
      if (index_f(rwrite, "LIST") != 0) ioulst = 1;
   }
   jfndnd = index_f(myreco, ".");
   if (jfndnd == 0) {
      printf("%s\n", terror[6]);
      ++ret_val;
   } else {
      myreco[jfndnd-1] = ' ';
   }
   if (strncmp(myreco, cbla72, 72) != 0) {
      printf("%s\n", terror[6]);
      ++ret_val;
   }

/* REWRITE TOP OF OBJECT FILE TO KEEP THE NUMBER OF RECORDS */
/* AND THE MAXIMUM LEVEL; TRANSMIT LAST STRING AS TITLE */
   nrecor = ninput;
   rwrite[72] = '\0';
   strcpy(header.cdatin, rwrite);
   header.nrecor = nrecor;
   header.ninput = ninput;
   header.maxlvl = maxlvl;
   header.nstack = nstack;
   header.ixrlst = ixrlst;
   header.ioulst = ioulst;
   header.idblst = idblst;
   iretcd = kdiput_c(iunito, (int_32 *)&header, 0, kdisize(header));
   if (iretcd != 0) goto L9001;

L666:
   return ret_val;

L5000:
   printf("! %s: INPUT LINE OVERFLOW\n", nomsub);
   strncpy(myreco, &recred[65], 72); myreco[72] = '\0';
   printf("! -->%s<--\n",myreco);
   ret_val = 5000;
   goto L666;
L5001:
   printf("! %s: INVALID 1-CHARACTER WORD IN CLE-2000\n", nomsub);
   ret_val = 5001;
   goto L666;
L5012:
   printf("! %s: MORE THAN 12-CHARACTER WORD IN CLE-2000\n", nomsub);
   ret_val = 5012;
   goto L666;
L7000:
   printf("! %s: KEYWORD= *%s*, BUT MAXIMUM NUMBER OF LEVELS IS ACHIEVED\n", nomsub, clogbg[ilogin-1]);
   printf("! %s: REVISE YOUR LOGIC\n", nomsub);
   ret_val = 7000;
   goto L666;
L7001:
   printf("! %s: AFTER *%s*, NOT EXPECTING KEYWORD= *%s\n",
          nomsub, clogbg[itylvl[nstlvl-1]-1], clogbg[ilogin-1]);
   printf("! %s: REVISE YOUR LOGIC\n", nomsub);
   ret_val = 7001;
   goto L666;
L7002:
   printf("! %s: KEYWORD= *%s*, BUT NOTHING LEFT FOR THIS LEVEL\n", nomsub, clogbg[ilogin-1]);
   printf("! %s: REVISE YOUR LOGIC\n", nomsub);
   ret_val = 7002;
   goto L666;
L7003:
   printf("! %s: KEYWORD= *%s*, BUT THE NUMBER OF EQUALS *:=* IS %d\n",
          nomsub, clogbg[ilogin-1], (int)nequal);
   ret_val = 7003;
   goto L666;
L7004:
   printf("! %s: KEYWORD= *%s*, BUT THE NUMBER OF WORDS IS %d\n",
          nomsub, clogbg[ilogin-1], (int)nwrsen);
   ret_val = 7003;
   goto L666;
L7005:
   printf("! %s: INVALID <<.>> OR >>.<< INSTRUCTION\n", nomsub);
   ret_val = 7005;
   goto L666;
L7006:
   printf("! %s: DECLARATION AS *%s* MUST APPEAR AT LOGIC LEVEL 1\n", nomsub, clogbg[ilogin-1]);
   ret_val = 7006;
   goto L666;
L7007:
   printf("! %s: INCONSISTENT END-OF-FILE, LOGIC LEVEL IS %d > 1\n", nomsub, (int)nlevel);
   printf("! %s: EXPECTING *%s* AT THE END OF STATEMENT\n", nomsub, crecbg);
   ret_val = 7007;
   goto L666;
L9001:
   iretcd = -1;
   printf("! %s: WRITING RETURN CODE = %d\n", nomsub, (int)iretcd);
   ret_val = iretcd;
   goto L666;
L9003:
   iretcd = -1;
   printf("! %s: READING RETURN CODE = %d\n", nomsub, (int)iretcd);
   ret_val = iretcd;
   goto L666;
} /* clelog */
