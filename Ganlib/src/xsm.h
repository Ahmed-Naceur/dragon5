
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

#define iofmax 30
#define maxit 100
#include "kdi.h"

typedef struct Block1 {    /* active directory resident-memory xsm structure */
   int_32 header;          /* header (=200 for an xsm file) */
   char hname[73];         /* character*72 name of the xsm file */
   int_32 listlen;         /* number of elements in the list */
   int_32 impf;            /* type of access (1:modif or 2:read-only) */
   int_32 idir;            /* offset of active directory on xsm file */
   struct Block2 *ibloc;   /* address of block 2 in memory */
   struct Db1 *icang;      /* address of the database handle */
   struct Block1 *father;  /* address of the father active directory resident-
                              memory xsm structure. =0 for root directory. */
   struct Db2 *icang2;     /* address of the xsmiof database handle */
} xsm ;

typedef struct Block2 {   /* active directory resident-memory xsm structure */
   kdi_file *ifile;       /* xsm (kdi) file handle */
   int_32 idir;           /* offset of active directory on xsm file */
   int_32 modif;          /* =1 if the active directory extent have been modified */
   int_32 ioft;           /* maximum address on xsm file */
   int_32 nmt;            /* exact number of nodes on the active directory extent */
   int_32 link;           /* offset of the next directory extent */
   int_32 iroot;          /* offset of any parent directory extent */
   char mynam[13];        /* character*12 name of the active directory. ='/' for the root level */
   int_32 iofs[iofmax];   /* offset list (position of the first element of each block
                             that belong to the active directory extent) */
   int_32 jlon[iofmax];   /* length of each record (jlong=0 for a directory) that belong
                             to the active directory extent */
   int_32 jtyp[iofmax];   /* type of each block that belong to the active directory extent */
   char cmt[iofmax][13];  /* list of character*12 names of each block (record or
                             directory) that belong to the active directory extent */
} block2 ;

typedef struct Db1{       /* database handle */
   int_32 nad;            /* number of addresses in the database */
   int_32 maxad;          /* maximum slots in the database */
   xsm **idir;            /* address of the array of pointers */
} db1 ;

typedef struct Db2{       /* xsmiof database handle */
   int_32 nad;            /* number of addresses in the database */
   int_32 maxad;          /* maximum slots in the database */
   int_32 ***iref;        /* address of the array of pointers addresses */
   int_32 **iofset;       /* address of the array of pointers */
   int_32 *lg;            /* address of the array of lengths */
} db2 ;

void xsmop_c(xsm **, char *, int_32, int_32);
void xsmput_c(xsm **, const char *, int_32, int_32, int_32 *);
void xsmget_c(xsm **, const char *, int_32 *);
void xsmcl_c(xsm **, int_32);
void xsmnxt_c(xsm **, char *);
void xsmlen_c(xsm **, const char *, int_32 *, int_32 *);
void xsminf_c(xsm **, char *, char *, int_32 *, int_32 *, int_32 *);
void xsmsix_c(xsm **, const char *, int_32 iact);
void xsmdid_c(xsm **, const char *, xsm **);
void xsmlid_c(xsm **, const char *, int_32, xsm **);
void xsmgpd_c(xsm **, const char *, int_32 **);
void xsmppd_c(xsm **, const char *, int_32, int_32, int_32 *);
