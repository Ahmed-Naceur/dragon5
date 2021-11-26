
/**********************************/
/* C API for xsm file support     */
/* author: A. Hebert (30/04/2002) */
/**********************************/

/*
   Copyright (C) 2002 Ecole Polytechnique de Montreal

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
 */

#include <stdlib.h>
#include <string.h>
#include "xsm.h"

#define iprim 3
#define iwrd 3
#define klong 5+iwrd+(3+iwrd)*iofmax
#if !defined(min)
#define min(A,B)  ((A) < (B) ? (A) : (B))
#endif
#define TRUE 1     /* valeur boolenne TRUE */
#define FALSE 0    /* valeur boolenne FALSE */

static char AbortString[132];

/* Table of constant values */

static int_32 c__0 = 0;
static int_32 c__1 = 1;
static int_32 c__2 = 2;
static int_32 c__8 = 8;
static char *bl12="                                  ";

void xsmkep(db1 *ipkeep, int_32 imode, xsm **iplist)
/*
 *-----------------------------------------------------------------------
 *
 * keep the addresses of the open active directories.
 *
 * input parameters:
 *  ipkeep : address of the database handle (always the same).
 *  imode  : =1: add to the database; =2: remove from the database.
 *  iplist : address of an active directory.
 *
 * output parameter:
 *  iplist : last active directory in the database. =0 if the
 *           database is empty.
 *
 * database handle structure:
 *  0 : number of addresses in the database.
 *  1 : maximum slots in the database.
 *  2 : address of the database.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmkep";
   int_32 n = ipkeep->nad;
   if (imode == 1) {
      int_32 i;
      xsm **my_parray;
      if ((*iplist)->header != 200) {
         sprintf(AbortString,"%s: WRONG HEADER(1).",nomsub);
         xabort_c(AbortString);
      } else if (ipkeep->nad + 1 > ipkeep->maxad) {
         ipkeep->maxad += maxit;
         my_parray = (xsm **) malloc((ipkeep->maxad)*sizeof(*my_parray));
         for (i = 0; i < n; ++i) my_parray[i]=ipkeep->idir[i];
         if (n > 0) free(ipkeep->idir);
         ipkeep->idir=my_parray;
      }
      ++ipkeep->nad;
      ipkeep->idir[n] = *iplist;
   } else if (imode == 2) {
      int_32 i, i0=0;
      for (i = n; i >= 1; --i) {
         if (ipkeep->idir[i-1] == *iplist) {
            i0 = i;
            goto L30;
         }
      }
      sprintf(AbortString,"%s: UNABLE TO FIND AN ADDRESS.",nomsub);
      xabort_c(AbortString);
L30:
      for (i = i0; i <= n-1; ++i)
         ipkeep->idir[i-1]=ipkeep->idir[i];
      --ipkeep->nad;
      if (ipkeep->nad == 0) {
         *iplist = NULL;
         free(ipkeep->idir);
         ipkeep->maxad=0;
         ipkeep->idir=NULL;
      } else {
         *iplist = ipkeep->idir[n-1];
         if ((*iplist)->header != 200) {
            sprintf(AbortString,"%s: WRONG HEADER(2).",nomsub);
            xabort_c(AbortString);
         }
      }
   } else {
      sprintf(AbortString,"%s: INVALID VALUE OF IMODE.",nomsub);
      xabort_c(AbortString);
   }
   return;
}

void xsmdir(int_32 *ind, block2 *my_block2)
/*
 *-----------------------------------------------------------------------
 *
 * import or export a directory using the kdi utility.
 *
 * input parameters:
 *  ind   : =1 for import ; =2 for export.
 *  my_block2 : address of memory-resident xsm structure (block 2).
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmdir";
   int_32 i,j,irc,ibuf[8],ipos,iofma2,iivec[6*iofmax];
   ipos = my_block2->idir;
   if (*ind == 1) {
      irc = kdiget_c(my_block2->ifile, ibuf, ipos, c__8);
      if (irc != 0) goto L40;
      if (strncmp((char*)ibuf,"$$$$",4) != 0) goto L30;
      iofma2 = ibuf[1];
      if (iofma2 > iofmax) goto L30;
      my_block2->nmt = ibuf[2];
      my_block2->link = ibuf[3];
      my_block2->iroot = ibuf[4];
      strncpy(my_block2->mynam,(char*)&ibuf[5],12);
      my_block2->mynam[12]='\0';
      for(i=11; i>0; i--) {
         if(my_block2->mynam[i] != ' ') break;
         my_block2->mynam[i]='\0';
      }
      if (my_block2->nmt == 0) return;
      ipos += c__8;
      irc = kdiget_c(my_block2->ifile, iivec, ipos, 6*iofma2);
      if (irc != 0) goto L40;
      for(i=0; i<my_block2->nmt; i++) {
         my_block2->iofs[i] = iivec[i];
         my_block2->jlon[i] = iivec[iofma2+i];
         my_block2->jtyp[i] = iivec[2*iofma2+i];
         strncpy(my_block2->cmt[i],(char*)&iivec[3*(iofma2+i)],12);
         my_block2->cmt[i][12]='\0';
         for(j=11; j>0; j--) {
            if(my_block2->cmt[i][j] != ' ') break;
            my_block2->cmt[i][j]='\0';
         }
      }
   } else if (*ind == 2) {
      my_block2->modif = 0;
      strncpy((char*)ibuf,"$$$$",4);
      ibuf[1] = iofmax;
      ibuf[2] = my_block2->nmt;
      ibuf[3] = my_block2->link;
      ibuf[4] = my_block2->iroot;
      memcpy((char*)&ibuf[5],bl12,12);
      strncpy((char*)&ibuf[5],my_block2->mynam,min(12,strlen(my_block2->mynam)));
      irc = kdiput_c(my_block2->ifile, ibuf, ipos, c__8);
      if (irc != 0) goto L50;
      ipos += c__8;
      memset(iivec, 0, 6*iofmax);
      for(i=0; i<my_block2->nmt; i++) {
         iivec[i] = my_block2->iofs[i];
         iivec[iofmax+i] = my_block2->jlon[i];
         iivec[2*iofmax+i] = my_block2->jtyp[i];
         memcpy((char*)&iivec[3*(iofmax+i)],bl12,12);
         strncpy((char*)&iivec[3*(iofmax+i)],my_block2->cmt[i],min(12,strlen(my_block2->cmt[i])));
      }
      irc = kdiput_c(my_block2->ifile, iivec, ipos, 6*iofmax);
      if (irc != 0) goto L50;
   }
   return;
/* ABORT ON FATAL ERRORS */
L30:
   sprintf(AbortString,"%s: UNABLE TO RECOVER DIRECTORY.",nomsub);
   xabort_c(AbortString);
L40:
   sprintf(AbortString,"%s: kdiget_c ERROR NB. %d.",nomsub,(int)irc);
   xabort_c(AbortString);
L50:
   sprintf(AbortString,"%s: kdiput_c ERROR NB. %d.",nomsub,(int)irc);
   xabort_c(AbortString);
}

void xsmop_c(xsm **iplist, char *namp, int_32 imp, int_32 impx)
/*
 *-----------------------------------------------------------------------
 *
 * open an existing or create a new xsm file.
 *
 * The xsm database results from the juxtaposition of a hierarchical
 * logical structure into a direct access file. The direct access file
 * is ANSI-C or Fortran-77 compatible and is managed by kdiget/put/cl.
 * xsmop/put/get/len/vec/nxt/cl entries provide a set of methods to
 * access a xsm file.
 *
 * The logical structure of a xsm file is made of a root directory fol-
 * lowed by variable-length blocks containing the useful information.
 * Each time a directory is full, an extent is automatically created at
 * the end of the file, so that the total number of blocks in a direc-
 * tory is only limited by the maximum size of the direct access file.
 * Any block can contain a sub-directory in order to create a hierar-
 * chical structure.
 *
 * input parameters:
 *    namp : character name (null terminated) of the xsm file.
 *     imp : type of access.  =0: new file mode;
 *                            =1: modification mode;
 *                            =2: read only mode.
 *    impx : if impx=0, we suppress printing on xsmop.
 *
 * output parameters:
 *  iplist : address of the handle to the xsm file.
 *    namp : character*12 name of the xsm file if imp=1 or imp=2.
 *
 * The active directory is made of two blocks linked together. A block 1
 * is allocated for each scalar directory or vector directory component.
 * Block 2 is unique for a given xsm file; every block 1 is pointing to
 * the same block 2.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmop_c";
   int_32 irc,ibuf;
   char hbuf[5];
   block2 *my_block2;
   db1 *my_db1;
   db2 *my_db2;
   kdi_file *my_file;
   if (imp < 0 || imp > 2) {
      sprintf(AbortString,"%s: INVALID ACTION ( %d ) ON XSM FILE '%s'.",nomsub,(int)imp,namp);
      xabort_c(AbortString);
   } else if (strlen(namp) > 72) {
      sprintf(AbortString,"%s: FINENAME '%s' EXCEEDING 72 CHARACTERS.",nomsub,namp);
      xabort_c(AbortString);
   }
   *iplist = (xsm *) malloc(sizeof(**iplist));
   my_block2 = (block2 *) malloc(sizeof(*my_block2));

   my_db1 = (db1 *) malloc(sizeof(*my_db1));
   my_db1->nad = 0;
   my_db1->maxad = 0;
   my_db1->idir = NULL;
   (*iplist)->icang = my_db1;

   my_db2 = (db2 *) malloc(sizeof(*my_db2));
   my_db2->nad = 0;
   my_db2->maxad = 0;
   my_db2->iref = NULL;
   my_db2->iofset = NULL;
   my_db2->lg = NULL;
   (*iplist)->icang2 = my_db2;

   (*iplist)->header = 200;
   (*iplist)->listlen = -1;
   (*iplist)->impf = imp;
   (*iplist)->ibloc = my_block2;
   (*iplist)->father = NULL;
   strcpy((*iplist)->hname,namp);
   my_file = (kdi_file *) kdiop_c(namp,imp);
   if (my_file == NULL) {
      sprintf(AbortString,"%s: UNABLE TO OPEN XSM FILE '%s'.",nomsub,namp);
      xabort_c(AbortString);
   }
   my_block2->ifile = my_file;
   if (impx > 1)
      printf("%s: KDI FILE OPEN = %ld NAME = '%s' ACTION = %d.\n",nomsub,(long)my_file->fd,namp,(int)imp);
   if (imp >= 1) {
/*    RECOVER THE ROOT DIRECTORY IF THE XSM FILE ALREADY EXISTS. */
      irc=kdiget_c(my_block2->ifile, &ibuf, c__0, c__1);
      if (irc != 0) goto L140;
      strncpy(hbuf,(char*)&ibuf,4);
      hbuf[4]='\0';
      if (strcmp(hbuf,"$XSM") != 0) {
         sprintf(AbortString,"%s: WRONG HEADER ON XSM FILE '%s'.",nomsub,namp);
         xabort_c(AbortString);
      }
      irc=kdiget_c(my_block2->ifile, &(my_block2->ioft), c__1, c__1);
      if (irc != 0) goto L140;
      irc=kdiget_c(my_block2->ifile, &(my_block2->idir), c__2, c__1);
      if (irc != 0) goto L140;
      (*iplist)->idir = my_block2->idir;
      xsmdir(&c__1,my_block2);
      my_block2->modif = 0;
      if (impx > 0) {
         printf("%s: XSM FILE RECOVERY. FILE = '%s'.\n",nomsub,namp);
         printf("%6s HIGHEST ATTAINABLE ADDRESS = %d\n"," ",(int)my_block2->ioft);
         printf("%6s ACTIVE DIRECTORY = %s\n"," ",my_block2->mynam);
      }
   } else {
/*    NEW-FILE MODE. */
      (*iplist)->impf = 1;
      (*iplist)->idir = iprim;
      my_block2->ioft = iprim+klong;
      my_block2->idir = iprim;
      my_block2->iroot = -1;
      my_block2->nmt = 0;
      my_block2->link = iprim;
      sprintf(my_block2->mynam,"/");
      memcpy((char*)&ibuf,"$XSM",sizeof(ibuf));
      irc=kdiput_c(my_block2->ifile, &ibuf, c__0, c__1);
      if (irc != 0) goto L150;
      irc=kdiput_c(my_block2->ifile, &(my_block2->ioft), c__1, c__1);
      if (irc != 0) goto L150;
      irc=kdiput_c(my_block2->ifile, &(my_block2->idir), c__2, c__1);
      if (irc != 0) goto L150;
      xsmdir(&c__2,my_block2);
      my_block2->modif = 1;
   }
   return;
L140:
   sprintf(AbortString,"%s: kdiget_c ERROR NB. %d ON XSM FILE '%s'.",nomsub,(int)irc,namp);
   xabort_c(AbortString);
L150:
   sprintf(AbortString,"%s: kdiput_c ERROR NB. %d ON XSM FILE '%s'.",nomsub,(int)irc,namp);
   xabort_c(AbortString);
}

void xsmrep(const char *namt, int_32 *ind, int_32 *idir, block2 *my_block2, int_32 *iii)
/*
 *-----------------------------------------------------------------------
 *
 * find a block (record or directory) position in the active directory
 * and related extents.
 *
 * input parameters:
 *  namt   : character*12 name of the required block.
 *  ind    : =1 search namt ; =2 search and positionning in an empty
 *           slot of the active directory if namt does not exists.
 *  idir   : offset of active directory on xsm file.
 *  my_block2 : address of memory-resident xsm structure (block 2).
 *
 * output parameter:
 *  iii    : return code. =0 if the block named namt does not exists;
 *           =position in the active directory extent if namt extsts.
 *           =0 or 1 if namt=' '.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmrep";
   int_32 i, ipos, ipos2, irc, irc2, istart;
   char namp[13],nomC[25];

   if (my_block2->idir != *idir) {
/*    SWITCH TO THE CORRECT ACTIVE DIRECTORY (BLOCK 2). */
      if (my_block2->modif == 1) xsmdir(&c__2, my_block2);
      my_block2->idir = *idir;
      xsmdir(&c__1, my_block2);
   }
   if (strcmp(namt,"***HANDLE***") == 0) {
      sprintf(AbortString,"%s: ***HANDLE*** IS A RESERVED KEYWORD.",nomsub);
      xabort_c(AbortString);
   }
   strcpy(namp,namt);
   if (strcmp(namp," ") == 0) strcpy(namp,"***HANDLE***");
   ipos = -1;
   if (my_block2->nmt < iofmax) ipos = my_block2->idir;
   if (my_block2->nmt == 0) goto L50;
   for (i = 1; i <= my_block2->nmt; ++i) {
      if (strcmp(namp,my_block2->cmt[i-1]) == 0) {
/*       THE BLOCK ALREADY EXISTS. */
         *iii = i;
         return;
      }
   }
/* THE BLOCK NAMP DOES NOT EXISTS IN THE ACTIVE DIRECTORY EXTENT. WE
   SEARCH IN OTHER EXTENTS THAT BELONG TO THE ACTIVE DIRECTORY. */
   if (my_block2->idir != my_block2->link) {
/*    RECOVER A NEW DIRECTORY EXTENT. */
      istart = my_block2->link;
      if (my_block2->modif == 1) xsmdir(&c__2, my_block2);
      my_block2->idir = istart;
L30:
      xsmdir(&c__1, my_block2);
      if (my_block2->nmt < iofmax) ipos = my_block2->idir;
      for (i = 1; i <= my_block2->nmt; ++i) {
         if (strcmp(namp,my_block2->cmt[i-1]) == 0) {
/*          THE BLOCK NAMP WAS FOUND IN THE ACTIVE DIRECTORY EXTENT. */
            *iii = i;
            return;
         }
      }
      if (my_block2->link == istart) goto L50;
      my_block2->idir = my_block2->link;
      goto L30;
   }
L50:
   *iii = 0;
   if (*ind == 1) return;
   if (ipos >= 0 && ipos != my_block2->idir) {
/*    AN EXTENT WITH AN EMPTY SLOT WAS FOUND. */
      if (my_block2->modif == 1) xsmdir(&c__2, my_block2);
      my_block2->idir = ipos;
      xsmdir(&c__1, my_block2);
   } else if (ipos == -1) {
/*    THE ACTIVE DIRECTORY IS FULL. CREATE AN EXTENT. */
      ipos = my_block2->link;
      my_block2->link = my_block2->ioft;
      if (my_block2->modif == 1) {
         xsmdir(&c__2, my_block2);
      } else {
         ipos2 = my_block2->idir + 3;
         irc=kdiput_c(my_block2->ifile, &(my_block2->link), ipos2, c__1);
         if (irc != 0) goto L150;
      }
      my_block2->idir = my_block2->link;
      my_block2->link = ipos;
      my_block2->ioft += klong;
      my_block2->nmt = 0;
   }
   ++my_block2->nmt;
   *iii = my_block2->nmt;
   my_block2->modif = 1;
   my_block2->jlon[*iii - 1] = 0;
   my_block2->jtyp[*iii - 1] = 99;
   strcpy(my_block2->cmt[*iii - 1],namp);
   return;

L150:
   irc2=kdicl_c(my_block2->ifile, c__1);
   strcpy(nomC,(my_block2->ifile)->nom);
   if (irc2 != 0) printf("%s: kdicl_c ERROR NB. %d ON XSM FILE '%s'.\n",nomsub,(int)irc2,nomC);
   sprintf(AbortString,"%s: kdiput_c ERROR NB. %d ON XSM FILE '%s'.",nomsub,(int)irc,nomC);
   xabort_c(AbortString);
}

void xsmput_c(xsm **iplist, const char *namp, int_32 ilong, int_32 itype, int_32 *data1)
/*
 *-----------------------------------------------------------------------
 *
 * copy a block from memory into the xsm file.
 *
 * input parameter:
 *  iplist : address of the handle to the xsm file.
 *    namp : character*12 name of the current block.
 *   ilong : number of information elements stored in the current block.
 *   itype : type of information elements stored in the current block.
 *           0: directory                1: integer
 *           2: single precision         3: character*4
 *           4: double precision         5: logical
 *           6: complex
 *   data1 : information elements. dimension data1(ilong)
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmput_c";
   char nomC[13];
   block2 *my_block2;
   int_32 iii,irc;
   if ((*iplist)->impf == 2) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' IS OPEN IN READ-ONLY MODE(1).",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (ilong <= 0) {
      sprintf(AbortString,"%s: INVALID LENGTH (%d) FOR NODE '%s' IN THE XSM FILE '%.60s'.",
              nomsub,(int)ilong,namp,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (itype <= 0 || itype >= 8) {
      sprintf(AbortString,"%s: INVALID TYPE NUMBER (%d) FOR NODE '%s' IN THE XSM FILE '%.60s'.",
              nomsub,(int)itype,namp,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' HAVE THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   my_block2=(*iplist)->ibloc;
   xsmrep(namp, &c__2, &(*iplist)->idir, my_block2, &iii);
   my_block2->modif = 1;
   int_32 jlong = ilong;
   if (itype == 4 || itype == 6) jlong = 2*ilong;
   if (jlong > my_block2->jlon[iii-1]) {
      my_block2->iofs[iii-1] = my_block2->ioft;
      my_block2->ioft += jlong;
   }
   my_block2->jlon[iii-1] = jlong;
   my_block2->jtyp[iii-1] = itype;
   irc=kdiput_c(my_block2->ifile, data1, my_block2->iofs[iii-1], jlong);
   if (irc != 0) {
      strcpy(nomC,(my_block2->ifile)->nom);
      sprintf(AbortString,"%s: kdiput_c ERROR NB. %d ON XSM FILE '%s'.",nomsub,(int)irc,nomC);
      xabort_c(AbortString);
   }
   return;
}
void xsmget_c(xsm **iplist, const char *namp, int_32 *data2)
/*
 *-----------------------------------------------------------------------
 *
 * copy a block from the xsm file into memory.
 *
 * input parameters:
 *  iplist : address of the handle to the xsm file.
 *   namp  : character*12 name of the current block.
 *
 * output parameter:
 *   data2 : information elements. dimension data2(ilong)
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmget_c";
   char nomC[13];
   block2 *my_block2;
   int_32 iii,irc;
   if ((*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' HAVE THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   my_block2=(*iplist)->ibloc;
   xsmrep(namp, &c__1, &(*iplist)->idir, my_block2, &iii);
   if (iii > 0) {
      irc=kdiget_c(my_block2->ifile, data2, my_block2->iofs[iii-1], my_block2->jlon[iii-1]);
      if (irc != 0) {
         strcpy(nomC,(my_block2->ifile)->nom);
         sprintf(AbortString,"%s: kdiget_c ERROR NB. %d ON XSM FILE '%s'.",nomsub,(int)irc,nomC);
         xabort_c(AbortString);
      }
   } else {
      sprintf(AbortString,"%s: UNABLE TO FIND BLOCK '%s' INTO DIRECTORY '%s' IN THE XSM FILE '%.45s'.",
              nomsub,namp,my_block2->mynam,(*iplist)->hname);
      xabort_c(AbortString);
   }
   return;
}

void xsmcl_c(xsm **iplist, int_32 istatu)
/*
 *-----------------------------------------------------------------------
 *
 * close the xsm file.
 *
 * input parameters:
 *  iplist : address of the handle to the xsm file.
 *  istatu : =1 to keep the file at close ; =2 to destroy it.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmcl_c";
   block2 *my_block2;
   int_32 i, irc, iii;
   db1 *ipkeep;
   db2 *ipkep2;
   if ((*iplist)->impf == 2 && istatu == 2) {
      sprintf(AbortString,"%s: CANNOT ERASE THE XSM FILE '%.60s' OPEN IN READ-ONLY MODE.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (istatu < 1 || istatu > 2) {
      sprintf(AbortString,"%s: INVALID ACTION ( %d ) ON XSM FILE '%s'.",
              nomsub,(int)istatu,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' HAVE THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   ipkeep = (*iplist)->icang;
   if (ipkeep->nad > 0) {
      for (i = 1; i <= ipkeep->nad; ++i) free(ipkeep->idir[i-1]);
      free(ipkeep->idir);
   }
   ipkep2 = (*iplist)->icang2;
   if (ipkep2->nad > 0) {
      for (i = 1; i <= ipkep2->nad; ++i) free(ipkep2->iofset[i-1]); /* rlsara_c() */
      free(ipkep2->iref);
      free(ipkep2->iofset);
      free(ipkep2->lg);
   }
   my_block2 = (*iplist)->ibloc;

   if (my_block2->modif == 1) {
      if ((*iplist)->impf == 2) {
         sprintf(AbortString,"%s: THE XSM FILE '%.60s' IS OPEN IN READ-ONLY MODE(2).",
                 nomsub,(*iplist)->hname);
         xabort_c(AbortString);
      }
      xsmdir(&c__2, my_block2);
   }

   if (my_block2->idir != (*iplist)->idir) {
/*    SWITCH TO THE CORRECT ACTIVE DIRECTORY (BLOCK 2). */
      my_block2->idir = (*iplist)->idir;
      xsmdir(&c__1, my_block2);
   }
   if (my_block2->iroot != -1) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' IS NOT ON ROOT DIRECTORY.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   if (my_block2->idir != iprim) {
/*    SWITCH TO THE CORRECT ACTIVE DIRECTORY (BLOCK 2). */
      my_block2->idir = iprim;
      xsmdir(&c__1, my_block2);
   }

   irc = kdiget_c(my_block2->ifile, &iii, c__1, c__1);
   if (irc != 0) goto L140;
   if (my_block2->ioft > iii) {
      if ((*iplist)->impf == 2) {
         sprintf(AbortString,"%s: THE XSM FILE '%.60s' IS OPEN IN READ-ONLY MODE(3).",
                 nomsub,(*iplist)->hname);
         xabort_c(AbortString);
      }
      irc = kdiput_c(my_block2->ifile, &(my_block2->ioft), c__1, c__1);
      if (irc != 0) goto L150;
   }
   irc = kdicl_c(my_block2->ifile, istatu);
   if (irc != 0) goto L160;

/* RELEASE THE XSM FILE HANDLE. */
   free((*iplist)->icang);
   free((*iplist)->icang2);
   my_block2->ifile = NULL;
   free(my_block2);
   (*iplist)->header = 0;
   free(*iplist);
   *iplist = NULL;
   return;

L140:
   sprintf(AbortString,"%s: kdiget_c ERROR NB. %d ON XSM FILE '%s'.",nomsub,(int)irc,(*iplist)->hname);
   xabort_c(AbortString);
L150:
   sprintf(AbortString,"%s: kdiput_c ERROR NB. %d ON XSM FILE '%s'.",nomsub,(int)irc,(*iplist)->hname);
   xabort_c(AbortString);
L160:
   sprintf(AbortString,"%s: kdicl_cS ERROR NB. %d ON XSM FILE '%s'.",nomsub,(int)irc,(*iplist)->hname);
   xabort_c(AbortString);
}

void xsmnxt_c(xsm **iplist, char *namp)
/*
 *-----------------------------------------------------------------------
 *
 * find the name of the next block stored in the active directory.
 *
 * input parameters:
 *  iplist : address of the handle to the xsm file.
 *    namp : character*12 name of a block. if namp=' ' at input, find
 *           any name for any block stored in this directory.
 *
 * output parameters:
 *    namp : character*12 name of the next block. namp=' ' for an empty
 *           directory.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmnxt_c";
   block2 *my_block2;
   int_32 iii;

   if ((*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' HAVE THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   my_block2 = (*iplist)->ibloc;
   if (strcmp(namp," ") == 0) {
      if (my_block2->idir != (*iplist)->idir) {
/*       SWITCH TO THE CORRECT ACTIVE DIRECTORY (BLOCK 2). */
         if (my_block2->modif == 1) xsmdir(&c__2, my_block2);
         my_block2->idir = (*iplist)->idir;
         xsmdir(&c__1, my_block2);
      }
      iii = min(my_block2->nmt,1);
   } else {
      xsmrep(namp, &c__1, &(*iplist)->idir, my_block2, &iii);
   }
   if (iii == 0 && strcmp(namp, " ") == 0) {
/*    EMPTY DIRECTORY */
      sprintf(AbortString,"%s: THE ACTIVE DIRECTORY '%s' OF THE XSM FILE '%.45s' IS EMPTY.",
              nomsub,my_block2->mynam,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iii == 0) {
      sprintf(AbortString,"%s: UNABLE TO FIND BLOCK '%s' INTO DIRECTORY '%s' IN THE XSM FILE '%.45s'.",
              nomsub,namp,my_block2->mynam,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iii + 1 <= my_block2->nmt) {
      strcpy(namp,my_block2->cmt[iii]);
      return;
   }
/* SWITCH TO THE NEXT DIRECTORY. */
   if (my_block2->idir != my_block2->link) {
      if (my_block2->modif == 1) xsmdir(&c__2, my_block2);
      my_block2->idir = my_block2->link;
/*    RECOVER THE NEXT DIRECTORY. */
      xsmdir(&c__1, my_block2);
   }
   strcpy(namp,my_block2->cmt[0]);
   if (strcmp(namp,"***HANDLE***") == 0) strcpy(namp," ");
   return;
}

void xsmlen_c(xsm **iplist, const char *namp, int_32 *ilong, int_32 *itype)
/*
 *-----------------------------------------------------------------------
 *
 * return the length and type of a block. Return 0 if the block does not
 * exists.
 *
 * input parameters:
 *  iplist : address of the handle to the xsm file.
 *   namp  : character*12 name of the current block.
 *   ilong : number of information elements stored in the current block.
 *           ilong=-1 is returned for a scalar directory.
 *           ilong=0 if the block does not exists.
 *   itype : type of information elements stored in the current block.
 *           0: directory                1: integer
 *           2: single precision         3: character*4
 *           4: double precision         5: logical
 *           6: complex                 99: undefined
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmlen_c";
   block2 *my_block2;
   int_32 iii;
   if ((*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' HAVE THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   my_block2=(*iplist)->ibloc;
   xsmrep(namp, &c__1, &(*iplist)->idir, my_block2, &iii);
   if (iii > 0) {
      *ilong = my_block2->jlon[iii-1];
      *itype = my_block2->jtyp[iii-1];
      if (*itype == 4 || *itype == 6) *ilong=*ilong/2;
   } else {
      *ilong = 0;
      *itype = 99;
   }
   return;
}
void xsminf_c(xsm **iplist, char *namxsm, char *nammy, int_32 *empty, int_32 *ilong, int_32 *access)
/*
 *-----------------------------------------------------------------------
 *
 * recover global informations related to an xsm file.
 *
 * input parameters:
 *  iplist : address of the handle to the xsm file.
 *
 * output parameters:
 *  namxsm : character*12 name of the xsm file.
 *   nammy : charecter*12 name of the active directory.
 *   empty : =.true. if the active directory is empty.
 *   ilong : =-1: for a table; >0: number of list items.
 *  access : type of access. =1: object open for modification;
 *           =2: object in read-only mode.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsminf_c";
   block2 *my_block2;
   if ((*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' HAVE THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   my_block2=(*iplist)->ibloc;
   if (my_block2->idir != (*iplist)->idir) {
/*    SWITCH TO THE CORRECT ACTIVE DIRECTORY (BLOCK 2). */
      if (my_block2->modif == 1) xsmdir(&c__2, my_block2);
      my_block2->idir = (*iplist)->idir;
      xsmdir(&c__1, my_block2);
   }
   strcpy(namxsm,(*iplist)->hname);
   strcpy(nammy,my_block2->mynam);

   *empty = (my_block2->nmt == 0);
   *ilong = (*iplist)->listlen;
   *access = (*iplist)->impf;
   return;
}

void xsmsix_c(xsm **iplist, const char *namp, int_32 iact)
/*
 *-----------------------------------------------------------------------
 *
 * move in the scalar hierarchical structure of a xsm file.
 *
 * input parameters:
 *  iplist : address of the father/son table.
 *    namp : character*12 name of the son table if iact=1.
 *           not used if iact=0 or iact=2.
 *    iact : type of movement in the hierarchical structure.
 *           0: move back to the root directory;
 *           1: move to a son vectorial directory;
 *           2: move back to the parent directory.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmsix_c";
   block2 *my_block2;
   int_32 iii, lenold, idir=0, ityold;
   xsm *iofset, *iofpre;
   char nomC[13];

   if (iact < 0 || iact > 2) {
      sprintf(AbortString,"%s: INVALID ACTION (%d) ON THE XSM FILE '%.60s'.",
              nomsub,(int)iact,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' HAVE THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   if (iact == 1) {
/*    MOVE TO A SON DIRECTORY. */
      my_block2=(*iplist)->ibloc;
      xsmrep(namp, &c__2, &(*iplist)->idir, my_block2, &iii);
      lenold = my_block2->jlon[iii-1];
      if (lenold == -1) lenold = 1;
      ityold = my_block2->jtyp[iii-1];
      if (lenold == 0) {
/*       CREATE A NEW SCALAR DIRECTORY EXTENT ON THE XSM FILE. */
         if ((*iplist)->impf == 2) {
            printf("new directory name=%s\n",namp);
            sprintf(AbortString,"%s: THE XSM FILE '%.60s' IS OPEN IN READ-ONLY MODE(4).",
                    nomsub,(*iplist)->hname);
            xabort_c(AbortString);
         }
         my_block2->jlon[iii-1] = -1;
         my_block2->jtyp[iii-1] = 0;
         my_block2->iofs[iii-1] = my_block2->ioft;
         idir = my_block2->iofs[iii-1];
         my_block2->ioft += klong;
         xsmdir(&c__2, my_block2);
         my_block2->iroot = my_block2->idir;
         strcpy(my_block2->mynam,namp);
         my_block2->idir = my_block2->iofs[iii-1];
         my_block2->nmt = 0;
         my_block2->link = my_block2->idir;
         my_block2->modif = 1;
      } else if (lenold == 1 && ityold == 0) {
         idir = my_block2->iofs[iii-1];
      } else if (ityold != 0) {
         sprintf(AbortString,"%s:  BLOCK '%s' IS NOT A DIRECTORY OF THE XSM FILE '%.60s'.",
                 nomsub,namp,(*iplist)->hname);
         xabort_c(AbortString);
      }
      iofset = (xsm *) malloc(sizeof(*iofset));

/*    COPY BLOCK1 */
      iofset->header = (*iplist)->header;
      strcpy(iofset->hname,(*iplist)->hname);
      iofset->listlen = -1;
      iofset->impf = (*iplist)->impf;
      iofset->idir = (*iplist)->idir;
      iofset->ibloc = (*iplist)->ibloc;
      iofset->icang = (*iplist)->icang;
      iofset->icang2 = (*iplist)->icang2;
      iofset->father = (*iplist)->father;

      xsmkep(iofset->icang, c__1, &iofset);
      iofset->idir  = idir;
      iofset->father = *iplist;
      *iplist = iofset;
   } else if (iact == 0 || iact == 2) {
/*    MOVE BACK TO THE ROOT OR PARENT DIRECTORY. */
L50:
      my_block2 = (*iplist)->ibloc;
      if (my_block2->modif == 1) xsmdir(&c__2, my_block2);
      if (my_block2->idir != (*iplist)->idir) {
/*       SWITCH TO THE CORRECT ACTIVE DIRECTORY (BLOCK 2). */
         my_block2->idir = (*iplist)->idir;
         xsmdir(&c__1, my_block2);
      }
      if (my_block2->iroot == -1) {
         if (iact == 0) {
            if ((*iplist)->header != 200) {
               sprintf(AbortString,"%s: WRONG HEADER(1).",nomsub);
               xabort_c(AbortString);
            }
            return;
         }
         sprintf(AbortString,"%s: THE XSM FILE '%.60s' IS ALREADY ON ROOT DIRECTORY.",
                 nomsub,(*iplist)->hname);
         xabort_c(AbortString);
      }
      strcpy(nomC,my_block2->mynam);
      iofset = *iplist;
      *iplist = (*iplist)->father;
      xsmrep(nomC, &c__1, &(*iplist)->idir, my_block2, &iii);
      if (iii == 0) {
         sprintf(AbortString,"%s: UNABLE TO STEP DOWN ON FATHER RECORD '%s' FOR XSM FILE '%.60s'.",
                 nomsub,namp,(*iplist)->hname);
         xabort_c(AbortString);
      }
      iofpre = iofset;
      xsmkep((*iplist)->icang, c__2, &iofpre);
      free(iofset);
      if (iact == 0 && (*iplist)->idir != iprim) {
         goto L50;
      } else if ((*iplist)->header != 200) {
         sprintf(AbortString,"%s: WRONG HEADER(2).",nomsub);
         xabort_c(AbortString);
      }
   }
   return;
}

void xsmdid_c(xsm **iplist, const char *namp, xsm **jplist)
/*
 *-----------------------------------------------------------------------
 *
 * create/access a daughter associative table in a father table.
 *
 * input parameters:
 *  iplist : address of the father table.
 *    namp : character*12 name of the daughter associative table.
 *
 * output parameter:
 *  jplist : address of the daughter associative table.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmdid_c";
   block2 *my_block2;
   int_32 iii, lenold, idir=0, ityold;
   xsm *iofset;

   if ((*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' HAVE THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   my_block2=(*iplist)->ibloc;
   xsmrep(namp, &c__2, &(*iplist)->idir, my_block2, &iii);
   lenold = my_block2->jlon[iii-1];
   ityold = my_block2->jtyp[iii-1];
   if (lenold == 0) {
/*    CREATE A NEW SCALAR DIRECTORY EXTENT ON THE XSM FILE. */
      if ((*iplist)->impf == 2) {
         sprintf(AbortString,"%s: THE XSM FILE '%.60s' IS OPEN IN READ-ONLY MODE(5).",
                 nomsub,(*iplist)->hname);
         xabort_c(AbortString);
      }
      my_block2->jlon[iii-1] = -1;
      my_block2->jtyp[iii-1] = 0;
      my_block2->iofs[iii-1] = my_block2->ioft;
      idir = my_block2->iofs[iii-1];
      my_block2->ioft += klong;
      xsmdir(&c__2, my_block2);
      my_block2->iroot = my_block2->idir;
      strcpy(my_block2->mynam,namp);
      my_block2->idir = my_block2->iofs[iii-1];
      my_block2->nmt = 0;
      my_block2->link = my_block2->idir;
      my_block2->modif = 1;
   } else if (lenold == -1 && ityold == 0) {
      idir = my_block2->iofs[iii-1];
   } else {
      sprintf(AbortString,"%s: BLOCK '%s' IS NOT AN ASSOCIATIVE TABLE OF THE XSM FILE '%.60s'.",
              nomsub,namp,(*iplist)->hname);
      xabort_c(AbortString);
   }
   iofset = (xsm *) malloc(sizeof(*iofset));
   *jplist = iofset;

/*  COPY BLOCK1 */
   (*jplist)->header = (*iplist)->header;
   strcpy((*jplist)->hname,(*iplist)->hname);
   (*jplist)->listlen = -1;
   (*jplist)->impf = (*iplist)->impf;
   (*jplist)->idir  = idir;
   (*jplist)->ibloc = (*iplist)->ibloc;
   (*jplist)->icang = (*iplist)->icang;
   (*jplist)->icang2 = (*iplist)->icang2;
   (*jplist)->father = *iplist;
   xsmkep((*iplist)->icang, c__1, &iofset);
   return;
}

void xsmlid_c(xsm **iplist, const char *namp, int_32 ilong, xsm **jplist)
/*
 *-----------------------------------------------------------------------
 *
 * create/access the hierarchical structure of a list in a xsm file.
 *
 * input parameters:
 *  iplist : address of the father table.
 *    namp : character*12 name of the daughter list.
 *   ilong : dimension of the daughter list.
 *
 * output parameter:
 *  jplist : address of the daughter list.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="xsmlid_c";
   char nomC[13];
   block2 *my_block2;
   int_32 iii, irc, irc2, lenold, idir=0, i, idiold, ityold, iroold, *iivec;
   xsm *iofset;

   if ((*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE XSM FILE '%.60s' HAVE THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (ilong <= 0) {
      sprintf(AbortString,"%s: INVALID LENGTH (%d) FOR NODE '%s' IN THE XSM FILE '%.60s'.",
              nomsub,(int)ilong,namp,(*iplist)->hname);
      xabort_c(AbortString);
   }
   my_block2=(*iplist)->ibloc;
   xsmrep(namp, &c__2, &(*iplist)->idir, my_block2, &iii);
   lenold = my_block2->jlon[iii-1];
   ityold = my_block2->jtyp[iii-1];
   if ((ilong > lenold && ityold == 10) || lenold == 0) {
/*    CREATE ILONG-LENOLD NEW LIST EXTENTS ON THE XSM FILE. */
      if ((*iplist)->impf == 2) {
         sprintf(AbortString,"%s: THE XSM FILE '%.60s' IS OPEN IN READ-ONLY MODE(6).",
                 nomsub,(*iplist)->hname);
         xabort_c(AbortString);
      }
      my_block2->jlon[iii-1] = ilong;
      my_block2->jtyp[iii-1] = 10;
      idiold = my_block2->iofs[iii-1];
      my_block2->iofs[iii-1] = my_block2->ioft;
      idir = my_block2->iofs[iii-1];
      my_block2->ioft += ilong;
      iroold = my_block2->idir;
      xsmdir(&c__2, my_block2);
      iivec = (int_32 *) malloc(ilong * sizeof(*iivec));
      if (lenold > 0) {
         irc = kdiget_c(my_block2->ifile, iivec, idiold, lenold);
         if (irc != 0) goto L110;
      }
      for (i = abs(lenold) + 1; i <= ilong; ++i) {
         iivec[i-1] = my_block2->ioft;
         my_block2->iroot = iroold;
         strcpy(my_block2->mynam,namp);
         my_block2->nmt = 0;
         my_block2->idir = my_block2->ioft;
         my_block2->ioft += klong;
         my_block2->link = my_block2->idir;
         xsmdir(&c__2, my_block2);
      }
      irc = kdiput_c(my_block2->ifile, iivec, idir, ilong);
      if (irc != 0) goto L100;
      free(iivec);
   } else if (ilong <= lenold && ityold == 10) {
      ilong = lenold;
      idir = my_block2->iofs[iii-1];
   } else if (ityold != 10) {
      sprintf(AbortString,"%s: BLOCK '%s' IS NOT A LIST OF THE XSM FILE '%.60s'.",
              nomsub,namp,(*iplist)->hname);
      xabort_c(AbortString);
   }
   iivec = (int_32 *) malloc(ilong * sizeof(*iivec));
   irc = kdiget_c(my_block2->ifile, iivec, idir, ilong);
   if (irc != 0) goto L110;
   *jplist = (xsm *) malloc(ilong * sizeof(**jplist));
   for (i = 0; i < ilong; ++i) {
      iofset = *jplist + i;

/*    COPY BLOCK1 */
      iofset->header = (*iplist)->header;
      strcpy(iofset->hname,(*iplist)->hname);
      iofset->listlen = 0;
      iofset->impf = (*iplist)->impf;
      iofset->idir = iivec[i];
      iofset->ibloc = (*iplist)->ibloc;
      iofset->icang = (*iplist)->icang;
      iofset->icang2 = (*iplist)->icang2;
      iofset->father = *iplist;
   }
   (*jplist)->listlen = ilong;
   xsmkep((*iplist)->icang, c__1, jplist);
   free(iivec);
   return;

L100:
   irc2=kdicl_c(my_block2->ifile, c__1);
   strcpy(nomC,(my_block2->ifile)->nom);
   if (irc2 != 0) printf("%s: kdicl_c ERROR NB. %d ON XSM FILE '%s'.\n",nomsub,(int)irc2,nomC);
   sprintf(AbortString,"%s: kdiput_c ERROR NB. %d ON XSM FILE '%s'.",nomsub,(int)irc,nomC);
   xabort_c(AbortString);
L110:
   irc2=kdicl_c(my_block2->ifile, c__1);
   if (irc2 != 0) printf("%s: kdicl_c ERROR NB. %d ON XSM FILE '%s'.\n",nomsub,(int)irc2,nomC);
   strcpy(nomC,(my_block2->ifile)->nom);
   sprintf(AbortString,"%s: kdiget_c ERROR NB. %d ON XSM FILE '%s'.",nomsub,(int)irc,nomC);
   xabort_c(AbortString);
}

void xsmgpd_c(xsm **iplist, const char *namp, int_32 **iofdum)
/*
 *-----------------------------------------------------------------------
 *
 * get a malloc pointer for an entry in the xsm file.
 *
 * input parameters:
 *  iplist : address of the handle to the xsm file.
 *    namp : character*12 name of the current block.
 *
 * output parameter:
 *  iofdum : malloc pointer to the xsm entry named namp.
 *
 *-----------------------------------------------------------------------
 */
{
   db2 *ipkep2;
   int_32 i, i0, ilong, itylcm, n;
   xsmlen_c(iplist,namp,&ilong,&itylcm);
   if(itylcm == 4 || itylcm == 6) ilong = 2*ilong;
   ipkep2 = (*iplist)->icang2;
   n = ipkep2->nad;
   for (i = n; i >= 1; --i) {
      if (ipkep2->iref[i-1] == iofdum) {
         i0 = i;
         goto L10;
      }
   }
   goto L20;
L10:
   if (ilong == ipkep2->lg[i0-1]) {
      *iofdum = ipkep2->iofset[i0-1];
      xsmget_c(iplist,namp,*iofdum);
      return;
   }
   free(ipkep2->iofset[i0-1]); /* rlsara_c() */
   for (i = i0; i <= n-1; ++i) {
      ipkep2->iref[i-1]=ipkep2->iref[i];
      ipkep2->iofset[i-1]=ipkep2->iofset[i];
      ipkep2->lg[i-1]=ipkep2->lg[i];
   }
   --ipkep2->nad;
   n = ipkep2->nad;
L20:
   *iofdum = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong) */
   xsmget_c(iplist,namp,*iofdum);
   if (n + 1 > ipkep2->maxad) {
      int_32 ***my_iref, **my_iofset, *my_lg;
      ipkep2->maxad += maxit;
      my_iref = (int_32 ***) malloc((ipkep2->maxad)*sizeof(*my_iref));
      my_iofset = (int_32 **) malloc((ipkep2->maxad)*sizeof(*my_iofset));
      my_lg = (int_32 *) malloc((ipkep2->maxad)*sizeof(*my_lg));
      for (i = 0; i < n; ++i) {
         my_iref[i]=ipkep2->iref[i];
         my_iofset[i]=ipkep2->iofset[i];
         my_lg[i]=ipkep2->lg[i];
      }
      if (n > 0) {
         free(ipkep2->iref);
         free(ipkep2->iofset);
         free(ipkep2->lg);
      }
      ipkep2->iref=my_iref;
      ipkep2->iofset=my_iofset;
      ipkep2->lg=my_lg;
   }
   ipkep2->iref[n] = iofdum;
   ipkep2->iofset[n] = *iofdum;
   ipkep2->lg[n] = ilong;
   ++ipkep2->nad;
   return;
}

void xsmppd_c(xsm **iplist, const char *namp, int_32 ilong, int_32 itype, int_32 *iofdum)
/*
 *-----------------------------------------------------------------------
 *
 * add a new malloc pointer entry in the xsm file.
 *
 * input parameter:
 *  iplist : address of the handle to the xsm file.
 *    namp : character*12 name of the current block.
 *   ilong : number of information elements stored in the current block.
 *   itype : type of information elements stored in the current block.
 *           0: directory                1: integer
 *           2: single precision         3: character*4
 *           4: double precision         5: logical
 *           6: complex
 *  iofdum : malloc pointer of the first information element.
 *
 *-----------------------------------------------------------------------
 */
{
   db2 *ipkep2;
   int_32 i, i0, n;
   xsmput_c(iplist,namp,ilong,itype,iofdum);
   ipkep2 = (*iplist)->icang2;
   n = ipkep2->nad;
   for (i = n; i >= 1; --i) {
      if (ipkep2->iofset[i-1] == iofdum) {
         i0 = i;
         goto L10;
      }
   }
   goto L20;
L10:
   for (i = i0; i <= n-1; ++i) {
      ipkep2->iref[i-1]=ipkep2->iref[i];
      ipkep2->iofset[i-1]=ipkep2->iofset[i];
      ipkep2->lg[i-1]=ipkep2->lg[i];
   }
   --ipkep2->nad;
   if (ipkep2->nad == 0) {
      free(ipkep2->iref);
      free(ipkep2->iofset);
      free(ipkep2->lg);
      ipkep2->maxad = 0;
      ipkep2->iref = NULL;
      ipkep2->iofset = NULL;
      ipkep2->lg = NULL;
   }
L20:
   free(iofdum); /* rlsara_c(iofdum) */
   iofdum = NULL;
   return;
}
