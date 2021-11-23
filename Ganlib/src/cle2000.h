
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. HEBERT ; 31/07/10      */
/*****************************************/

/*
Copyright (C) 2010 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/

#include "lcm.h"
#define nmaskc 3   /* size of maskck and ipacki */
#define kdisize(object) (sizeof(object)+sizeof(int_32)-1)/sizeof(int_32)
#if !defined(max)
#define max(A,B)  ((A) > (B) ? (A) : (B))
#define min(A,B)  ((A) < (B) ? (A) : (B))
#endif
#define lrclen max(kdisize(record1),kdisize(record2))  /* max kdisize of record1 and record2 */

typedef struct LIFO {          /* last-in-first-out (lifo) stack */
   int_32 nup;                 /* number of nodes in stack */
   struct LIFO_NODE *root;     /* address of the first node in stack */
   struct LIFO_NODE *node;     /* address of the last node in stack */
} lifo ;

typedef struct LIFO_NODE {     /* node in last-in-first-out (lifo) stack */
   int_32 type;                /* type of node: 3= lcm object; 4= xsm file; 5= seq binary;
                                  6= seq ascii; 7= da binary; 11= integer value; 12= real value
                                  13= character string; 14= double precision value;
                                  15= logical value */
   int_32 access;              /* 0=creation mode/1=modification mode/2=read-only mode */
   int_32 lparam;              /* record length for DA file objects */
   union {
       int_32 ival;            /* integer or logical value */
       float_32 fval;          /* real value */
       double dval;            /* double precision value */
       lcm *mylcm;             /* handle towards a LCM object */
       char hval[73];          /* character value */
   } value;
   struct LIFO_NODE *daughter; /* address of the daughter node in stack */
   char name[13];              /* name of node in the calling script */
   char name_daughter[13];     /* name of node in the daughter script */
   char OSname[73];            /* physical filename */
} lifo_node ;

void cleopn(lifo **);
lifo_node * clepop(lifo **);
void clepush(lifo **, lifo_node *);
int_32 clecls(lifo **);
lifo_node * clenode(lifo **, const char *);
lifo_node * clepos(lifo **, int_32);
void clelib(lifo **);

void redopn_c(kdi_file *, FILE *, char *, int_32);
void redcls_c(kdi_file **, FILE **, char [73], int_32 *);
void redget_c(int_32 *, int_32 *, float_32 *, char [73], double *);
void redput_c(int_32 *, int_32 *, float_32 *, char *, double *);

int_32 cle2000_c(int_32,
       int_32 (*)(char *, int_32, char (*)[13], int_32 *, int_32 *, lcm **, char (*)[73]),
       char *, int_32, lifo *);
int_32 clemod_c(char *, FILE *, int_32, char (*)[13], int_32 *, int_32 *, lcm **, char (*)[73],
       int_32 (*)(char *, int_32, char (*)[13], int_32 *, int_32 *, lcm **, char (*)[73]));
int_32 kdrcln(lifo *, int_32);

void drviox(lifo *, int_32, int_32 *);
int_32 objpil(kdi_file *, FILE *, int_32);
int_32 objstk(kdi_file *, FILE *, int_32);
int_32 objxrf(kdi_file *, FILE *);
int_32 clepil(FILE *, FILE *, kdi_file *, int_32 (*clecst)(char *, int_32 *, int_32 *, float_32 *, char *, double *));
int_32 clecst(char *, int_32 *, int_32 *, float_32 *, char *, double *);
int_32 clecop(kdi_file *, kdi_file *);
int_32 clexrf(kdi_file *, FILE *);
int_32 clestk(kdi_file *, FILE *,int_32 (*)(char *, int_32 *, int_32 *, float_32 *, char *, double *));
int_32 clelog(FILE *, FILE *, kdi_file *);

int_32 kdrdpr(lifo **, int_32, char (*)[13]);
int_32 kdrprm(lifo **, lifo **, int_32 minput, int_32, int_32 *, char (*)[13]);
int_32 kdrdmd(lifo **, int_32, char (*)[13]);
int_32 kdrdll(lifo **, int_32, char (*)[13]);
int_32 kdrdxf(lifo **, int_32, char (*)[13]);
int_32 kdrdsb(lifo **, int_32, char (*)[13]);
int_32 kdrdsa(lifo **, int_32, char (*)[13]);
int_32 kdrdda(lifo **, int_32, char (*)[13]);
