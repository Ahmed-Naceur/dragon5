
/**********************************/
/* C API for lcm object support   */
/* author: A. Hebert (30/04/2002) */
/**********************************/

/*
Copyright (C) 2002 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/
#ifndef lcm_H
#define lcm_H

#define maxext 70
#define lhash 251
#include <stdio.h>
#include "xsm.h"

typedef struct Blocka {   /* fixed length resident-memory structure */
   int_32 header;         /* header (=100 remory-resident; =200 xsm file) */
   char hname[73];        /* character*72 name of the lcm object */
   int_32 listlen;        /* number of elements in the list */
   struct Blockb *inext;  /* address of block 2 array in memory */
   struct Blocka *father; /* address of the father lcm object. =null for root directory */
   int_32 ifdir;          /* 0 / record index in father table */
   int_32 imode;          /* 0=closed/1=modification mode/2=read-only mode */
   int_32 imax;           /* maximum number of records in table */
   int_32 inref;          /* exact number of records in table */
   struct Db0 *icang;     /* address of the directory database handle */
   struct Dbref *global;  /* address of the global variable database handle */
   int_32 *hash;          /* hash table address */
} lcm ;

typedef struct Blockb {   /* variable length resident-memory structure */
   int_32 *jdata;         /* data offset */
   int_32 jjlon;          /* record length (in words) */
   int_32 jjtyp;          /* record type */
   int_32 jidat[4];       /* first/last element of record (4 words) */
   int jcmt[3];           /* character*12 name of record (3 words) */
} blockb ;

typedef struct Db0{       /* database handle */
   int_32 nad;            /* number of addresses in the database */
   int_32 maxad;          /* maximum slots in the database */
   lcm **idir;            /* address of the array of pointers */
} db0 ;

typedef struct Dbref{     /* database handle */
   int_32 nad;            /* number of addresses in the database */
   int_32 maxad;          /* maximum slots in the database */
   int_32 **local;        /* address of the array of local references */
} dbref ;

void lcmop_c(lcm **, char *, int_32, int_32, int_32);
void lcmppd_c(lcm **, const char *, int_32, int_32, int_32 *);
void lcmlen_c(lcm **, const char *, int_32 *, int_32 *);
void lcminf_c(lcm **, char *, char *, int_32 *, int_32 *, int_32 *, int_32 *);
void lcmnxt_c(lcm **, char *);
void lcmgpd_c(lcm **, const char *, int_32 **);
void lcmget_c(lcm **, const char *, int_32 *);
void lcmval_c(lcm **, const char *);
void lcmcl_c(lcm **, int_32);
void lcmdel_c(lcm **,const char *);
void lcmsix_c(lcm **, const char *, int_32);
lcm * lcmdid_c(lcm **, const char *);
lcm * lcmlid_c(lcm **, const char *, int_32);
lcm * lcmdil_c(lcm **, int_32);
lcm * lcmlil_c(lcm **, int_32, int_32);
void lcmppl_c(lcm **, int_32, int_32, int_32, int_32 *);
void lcmlel_c(lcm **, int_32, int_32 *, int_32 *);
void lcmgpl_c(lcm **, int_32, int_32 **);
void lcmgdl_c(lcm **, int_32, int_32 *);
lcm * lcmgid_c(lcm **, const char *);
lcm * lcmgil_c(lcm **, int_32);
void lcmequ_c(lcm **,lcm **);
void lcmput_c(lcm **,const char *,int_32,int_32,int_32 *);
void lcmpdl_c(lcm **,int_32,int_32,int_32 itype,int_32 *);
void lcmpcd_c(lcm **,const char *,int_32,char *[]);
void lcmgcd_c(lcm **,const char *,char *[]);
void lcmpcl_c(lcm **,int_32,int_32 iint_32,char *[]);
void lcmgcl_c(lcm **,int_32,char *[]);
void lcmpsd_c(lcm **,const char *,char *);
char * lcmgsd_c(lcm **,const char *);
void lcmpsl_c(lcm **,int_32,char *);
char * lcmgsl_c(lcm **,int_32);
void lcmlib_c(lcm **);
void lcmexp_c(lcm **, int_32, FILE *, int_32, int_32);
void lcmexpv3_c(lcm **, int_32, FILE *, int_32, int_32);
void refpush(lcm **, int_32 *);
int_32 refpop(lcm **, int_32 *);
void strcut_c(char *, char *, int_32);
void strfil_c(char *, char *, int_32);
#endif
