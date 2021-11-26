
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 15/05/09      */
/*****************************************/

#include <string.h>
#include "cle2000.h"
#include "header.h"

int_32 clecop(kdi_file *iuniti, kdi_file *iunito)
{
   static char cl2000[] = "CLE2000(V3)";

/*     CLE-2000 SYSTEM: R.ROY (09/1999), VERSION 2.1 */

/*             *CLECOP* WILL COPY ONE OBJECT FILE TO ANOTHER */

/*      INPUT: *IUNITI* IS THE DIRECT ACCESS UNIT FOR OBJECT CODE */
/*             *IUNITO* IS ITS COPY */

/*       NOTE: FACULTATIVE IN THE CLE-2000 SYSTEM. */
/*             THIS UTILITARY FUNCTION CAN BE USED IN THE APPLICATION */
/*             (MUST BE USED IN CASES WHERE THE OBJECT FILE IS RECURSIVE) */

   int_32 iretcd, iofset;
   int_32 ret_val = 0;
   int_32 irecor;

/* READ AND WRITE TOP OF OBJECT FILE */
   iretcd = kdiget_c(iuniti, (int_32 *)&header, 0, kdisize(header));
   if (iretcd != 0) goto L9023;
   iretcd = kdiput_c(iunito, (int_32 *)&header, 0, kdisize(header));
   if (iretcd != 0) goto L9025;

/* VERIFY CONSISTENCY OF VERSION */
   if (strcmp(header.cparin, cl2000) != 0) goto L9002;
   if (header.nrecor <= 0) goto L9003;
   if (header.nrecor != header.ninput + header.nstack + header.nobjet) goto L9003;

/* CASE WHERE THERE ARE NO LOG-FILE (TRIVIAL FILE) */
   if (header.ninput == 1) goto L666;

/* COPY ALL LOG-FILE RECORDS */
   for (irecor = 1; irecor < header.ninput; ++irecor) {
      iofset = irecor * lrclen;
      iretcd = kdiget_c(iuniti, (int_32 *)&record1, iofset, kdisize(record1));
      if (iretcd != 0) goto L9023;
      iretcd = kdiput_c(iunito, (int_32 *)&record1, iofset, kdisize(record1));
      if (iretcd != 0) goto L9025;
   }

/* CASE WHERE THERE ARE NO CLE-2000 VARIABLES */
   if (header.nstack + header.nobjet <= 0) goto L666;

/* COPY ALL VAR-STACK AND OBJECT RECORDS */
   for (irecor = header.ninput; irecor < header.nrecor; ++irecor) {
      iofset = irecor * lrclen;
      iretcd = kdiget_c(iuniti, (int_32 *)&record2, iofset, kdisize(record2));
      if (iretcd != 0) goto L9023;
      iretcd = kdiput_c(iunito, (int_32 *)&record2, iofset, kdisize(record2));
      if (iretcd != 0) goto L9025;
   }

L666:
   return ret_val;
L9002:
   ret_val = 9002;
   goto L666;
L9003:
   ret_val = 9003;
   goto L666;
L9023:
   ret_val = 9023;
   goto L666;
L9025:
   ret_val = 9025;
   goto L666;
} /* clecop */
