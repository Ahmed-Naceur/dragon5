
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 11/05/09      */
/*****************************************/

#include <stdlib.h>
#include <string.h>
#include "cle2000.h"
#include "header.h"
#define ndclkw 8
#define ntotxr 7
#define nmawrd 36

int_32 objxrf(kdi_file *iunito, FILE *iwrite)
{
   char *nomsub = "objxrf";
   static char cl2000[] = "CLE2000(V3)";
   static char ctitdb[] = "* GAN-2000 VERS 2.1 * DEBUG (WARNINGS AND ERRORS)";
   static char ctitxr[] = "* GAN-2000 VERS 2.1 * CROSS REFERENCE LISTING";
   static char *cdclkw[] = {"PROCEDURE", "MODULE", "LINKED_LIST", "XSM_FILE",
                            "SEQ_BINARY", "SEQ_ASCII", "DIR_ACCESS", "PARAMETER"};
   static char *ctypes[] = {"PR", "MD", "LL", "XF", "SB", "SA", "DA", "--"};

/*     CLE-2000 SYSTEM: R.ROY (11/1999), VERSION 2.1 */

/*             *OBJXRF* X-REF FOR OBJECTS ON THE D.A. UNIT  *IUNITO* */
/*                      OUTPUT IS WRITTEN ON UNIT           *IWRITE* */

/*        USE:          DRESS UP A LIST OF OBJECTS AND LINES WHERE USED. */
/*                      IN DEBUG MODE, ATTEMPT TO LIST POSSIBLE ERRORS. */

/*      INPUT: *IUNITO* IS THE DIRECT ACCESS UNIT FOR OBJECT CODE */
/*             *IWRITE* IS THE OUTPUT UNIT */

/*       NOTE: *OBJXRF* = 0 IF NO PROBLEM WAS ENCOUNTERED WHILE COMPILING */

   int_32 ret_val = 0;
   int_32 i, iretcd, ninput, nstack, ixrlst, idblst, nobjet, indlin, idclin, idefin, iusein,
          iofset, ilines, ilevel;
   char cparin[13], myreco[73], cdatin[73], cerror[13], cdefst[21], cusest[21], myparm[13];
   int_32 maskck[nmaskc], ipacki[nmaskc];
   int_32 idebwd[nmawrd+1], ifinwd[nmawrd+1], jndlec[nmawrd];
   int_32 lequal=0, istack, linxrf[ntotxr];
   char my_header[35];

/* READ TOP OF OBJECT FILE */
   iretcd = kdiget_c(iunito, (int_32 *)&header, 0, kdisize(header));
   if (iretcd != 0) goto L9023;
   strcpy(cparin, header.cparin);
   strcpy(myreco, header.cdatin);
   ninput = header.ninput;
   nstack = header.nstack;
   ixrlst = header.ixrlst;
   idblst = header.idblst;
   nobjet = header.nobjet;
   if (strcmp(cparin, cl2000) != 0) goto L9025;

/* CASES WHERE THERE ARE NO OBJECTS/MODULE */
   if (nobjet == 0) goto L666;

/* CASE WHERE DEBUG IS ACTIVE */
   if (idblst > 0) {
      int_32 lfirst = 1;
      int_32 irecor;
      strcpy(cerror, "            ");
      strcpy(cdefst, "                    ");
      strcpy(cusest, "                    ");
      for (irecor = ninput + nstack + 1; irecor <= ninput + nstack + nobjet; ++irecor) {
         iofset = (irecor - 1) * lrclen;
         iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
         if (iretcd != 0) goto L9023;
         strcpy(cparin, record2.cparin);
         strcpy(cdatin, record2.cdatin);
         indlin = record2.indlin;
         idclin = record2.idclin;
         idefin = record2.idefin;
         iusein = record2.iusein;
         if (abs(indlin) == 1 || abs(indlin) == 2 || abs(indlin) == 8) goto L60;
         if (idefin == 0) {
            strcpy(cdefst, "NOT DEFINED");
            if (iusein == 0) {
               strcpy(cusest, "NOT USED");
               strcpy(cerror, "WARNING");
            } else {
               strcpy(cusest, "BEFORE DEFINED");
               strcpy(cerror, "EXTERNAL");
            }
         } else if (iusein == 0) {
            strcpy(cusest, "NOT USED");
            strcpy(cerror, "WARNING");
         } else {
            if (idclin > idefin) {
               strcpy(cdefst, "BEFORE DECLARED");
               strcpy(cerror, "ERROR");
               ++ret_val;
            }
            if (idefin > iusein) {
               strcpy(cusest, "BEFORE DEFINED");
               strcpy(cerror, "ERROR");
               ++ret_val;
            }
         }
         if (lfirst && strcmp(cerror, "            ") != 0) {
            fprintf(iwrite, "\n");
            fprintf(iwrite, "%-72s\n", ctitdb);
            fprintf(iwrite, " REPORT-----/OBJECT------/DEFINED-STATUS------/USED-STATUS---------\n");
            lfirst = 0;
         }
         if (strcmp(cerror, "            ") != 0) {
            fprintf(iwrite, "%-12s/%-12s/%-20s/%-20s\n", cerror, cparin, cdefst, cusest);
         }
L60:
         ;
      }
      if (!lfirst) {
         if (ret_val > 0) {
            fprintf(iwrite, " REPORT-----> NB. OF ERRORS=%7d\n", (int)ret_val);
            fprintf(iwrite, " REPORT-----> MAY STILL EXECUTE WELL...\n");
         }
         fprintf(iwrite, " \n");
      }
   }

/* CASE WHERE NO XREF WAS ASKED */
   if (ixrlst <= 0) goto L666;
   fprintf(iwrite, "%-72s\n", ctitxr);
   fprintf(iwrite, " \n");
   fprintf(iwrite, "     OBJECT      TYPE LIN_DCL    ****   FOUND IN LINES (- MEANS NEW EVALUATION)  ****\n");
   fprintf(iwrite, " \n");

/* ***  MAIN LOOP OVER VARIABLES (BEGIN) */
   for (istack = ninput + nstack + 1; istack <= ninput + nstack + nobjet; ++istack) {
      int_32 ilogin = 0;
      int_32 jlines = 0;
      int_32 iuseln = 0;
      int_32 idefln = 0;
      int_32 nxreft = 0;
      int_32 irecor;

      iofset = (istack - 1) * lrclen;
      iretcd = kdiget_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
      if (iretcd != 0) goto L9023;
      strcpy(myparm, record2.cparin);
      strcpy(cdatin, record2.cdatin);
      indlin = record2.indlin;
      idclin = record2.idclin;
      idefin = record2.idefin;
      iusein = record2.iusein;

/*    IF WE DO NOT WANT TO SEE PROCEDURES */
/*    IF( ABS(INDLIN).EQ.1 ) GO TO 200 */
/*    IF WE DO NOT WANT TO SEE MODULES */
      if (abs(indlin) == 2) goto L200;

/*    PREPARE HEADER FOR VARIABLE *MYPARM* */
      sprintf(my_header, " %-4d %-12s %-2s    %04d_", (int)istack, myparm, ctypes[abs(indlin)-1], (int)idclin);

/*    ***  MAIN LOOP OVER RECORDS   (BEGIN) */
      for (irecor = 2; irecor <= ninput; ++irecor) {
         int_32 iloop1, jloop2;
         char cparav[13];

/*       READ A NEW RECORD */
         iofset = (irecor - 1) * lrclen;
         iretcd = kdiget_c(iunito, (int_32 *)&record1, iofset, kdisize(record1));
         if (iretcd != 0) goto L9023;
         strcpy(cparin, record1.cparin);
         strcpy(myreco, record1.myreco);
         ilines = record1.ilines;
         ilevel = record1.ilevel;
         for (i = 0; i < nmaskc; i++) maskck[i] = record1.maskck[i];
         for (i = 0; i < nmaskc; i++) ipacki[i] = record1.ipacki[i];

/*       RECORDS OUTSIDE CLE-2000 */
         if (ilevel == 0) {
            int_32 iwords = 1;
            int_32 nwords = 1;
            int_32 jbiprv = 0;

/*          SKIP OBJECT/MODULE DECLARATIONS */
            if (strcmp(cparin, " ") != 0) {
               ilogin = 0;
               for (iloop1 = 1; iloop1 <= ndclkw; ++iloop1) {
                  if (strcmp(cparin, cdclkw[iloop1-1]) == 0) ilogin = iloop1;
               }
               lequal = 0;
            }
            if (ilogin != 0) goto L100;

/*          HERE, WE HAVE FOUND A SENTENCE INCLUDING A STACK... */

/*          BEGIN: MASK RECOVERY */
            for (iloop1 = 1; iloop1 <= 72; ++iloop1) {
               int_32 jbicur;
               jloop2 = (iloop1 + 23) / 24;
               jbicur = maskck[jloop2-1] % 2;
               iwords += jbiprv * (1 - jbicur);
               idebwd[nwords-1] = iloop1;
               ifinwd[iwords-1] = iloop1;
               nwords += jbicur * (1 - jbiprv);
               jbiprv = jbicur;
               maskck[jloop2 - 1] /= 2;
            }
            --nwords;
/*          END:   MASK RECOVERY */

/*          BEGIN: UNPACK JNDLEC WITH TYPES (ITYP-1) */
            for (iloop1 = 1; iloop1 <= nwords; ++iloop1) {
               jloop2 = ((iloop1 << 1) + 23) / 24;
               jndlec[iloop1-1] = ipacki[jloop2-1] % 4;
               ipacki[jloop2-1] /= 4;
            }
/*          END:   UNPACK JNDLEC WITH TYPES (ITYP-1) */

            for (iloop1 = 1; iloop1 <= nwords; ++iloop1) {
               if (jndlec[iloop1-1] == 2 && myreco[idebwd[iloop1-1]-1] != '\''
                   && ifinwd[iloop1-1] - idebwd[iloop1-1] <= 11) {
                  strncpy(cparav, &myreco[idebwd[iloop1-1]-1], ifinwd[iloop1-1]-idebwd[iloop1-1]+1);
                  cparav[ifinwd[iloop1-1]-idebwd[iloop1-1]+1] = '\0';
                  if (strcmp(cparav, ";") == 0 || strcmp(cparav, "::") == 0) {
                     if (!lequal && idefln != 0) iuseln = idefln;

/*                   RESET *ILOGIN* AND LEFT/RIGHT NUMBERS */
                     ilogin = -10;
                     idefln = 0;
                     goto L55;
                  } else if (lequal) {
                     if (strcmp(cparav, myparm) == 0) {
/*                              USING THIS VARIABLE */
                        if (iuseln == 0) iuseln = ilines;
                     }
                  } else {
                     if (strcmp(cparav, ":=") == 0) {
                        lequal = 1;
/*                      *:=* SIGN IMPLIES REDEFINITION */
                        if (idefln != 0) iuseln = -idefln;
                        idefln = 0;
                     } else {
/*                      KEEP THE DEFINITION LINE UNTIL *:=* OR END */
                        if (strcmp(cparav, myparm) == 0) idefln = ilines;
                     }
                  }
               }
            }
L55:
            ;
         }

/*       HAVE WE FOUND A NEW XREF LINE ? */
         if (iuseln != 0 && iuseln != jlines) {
            if (nxreft == ntotxr) {
               char xline[81];
               sprintf(&xline[0], "%-24s", my_header);
               for (i = 0; i < ntotxr; i++) sprintf(&xline[24 + 8*i], "    %04d", (int)linxrf[i]);
               fprintf(iwrite, "%-80s\n", xline);
               strcpy(my_header, " ");
               nxreft = 0;
            }
            ++nxreft;
            linxrf[nxreft-1] = iuseln;
            jlines = iuseln;
         }
         iuseln = 0;

L100:
         ;
      }
/*      ***  MAIN LOOP OVER RECORDS   (END) */

/*    POSSIBLE INCOMPLETE LAST LINE... */
      if (nxreft != 0) {
         char xline[81];
         sprintf(&xline[0], "%-24s", my_header);
         for (i = 0; i < nxreft; i++) sprintf(&xline[24 + 8*i], "    %04d", (int)linxrf[i]);
         fprintf(iwrite, "%-80s\n", xline);
      } else if (strcmp(my_header, " ") != 0) {
         fprintf(iwrite, "%-24s <= WARNING: NEVER DEFINED, NEVER USED... POSSIBLE ERROR\n", my_header);
      }
L200:
      ;
   }
/* ***  MAIN LOOP OVER VARIABLES (END) */

   fprintf(iwrite, " \n");
L666:
   return ret_val;

L9023:
   iretcd = -1;
   printf("! %s: IOSTAT RETURN CODE =%d\n", nomsub,(int)iretcd);
   printf("! %s: IMPOSSIBLE TO USE THIS *OBJECT* FILE\n", nomsub);
   ret_val = -2;
   goto L666;
L9025:
   printf("! %s: IMPOSSIBLE TO USE OLD  *OBJECT* FILE\n", nomsub);
   ret_val = -3;
   goto L666;
} /* objxrf */
