
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

#include <stdlib.h>
#include <string.h>
#include "lcm.h"

#if !defined(max)
#define min(A,B)  ((A) < (B) ? (A) : (B))
#define max(A,B)  ((A) > (B) ? (A) : (B))
#endif

static int_32 c__1 = 1;
static int_32 c__2 = 2;
static char AbortString[132];

FILE * stdfil_c(char *s)
/*
 *-----------------------------------------------------------------------
 *
 * Return standard file pointers in ANSI C.
 *
 *-----------------------------------------------------------------------
 */
{
   if (strcmp(s, "stdin") == 0) {
      return stdin;
   } else if (strcmp(s, "stdout") == 0) {
      return stdout;
   } else if (strcmp(s, "stderr") == 0) {
      return stderr;
   } else {
      return NULL;
   }
}

void strcut_c(char *s, char *ct, int_32 n)
/*
 *-----------------------------------------------------------------------
 *
 * Copy n characters from string ct to s. Eliminate leading ' ' and '\0'
 * characters in s. Terminate s with a '\0'.
 *
 *-----------------------------------------------------------------------
 */
{
   int i;
   for(i=n-1; i>0; i--) {
      if(ct[i] != ' ' && ct[i] != '\0') break;
   }
   strncpy(s, ct, i+1); s[i+1] = '\0';
}

void strfil_c(char *s, char *ct, int_32 n)
/*
 *-----------------------------------------------------------------------
 *
 * Copy n characters from string ct to s. Eliminate '\0' characters and
 * pack with ' '. Assume that ct is null-terminated.
 *
 *-----------------------------------------------------------------------
 */
{
   int i;
   for(i=0; i<n; i++) s[i] = ' ';
   for(i=min(n,strlen(ct))-1; i>0; i--) {
      if(ct[i] != ' ' && ct[i] != '\0') break;
   }
   strncpy(s, ct, i+1);
}

void refpush(lcm **iplist, int_32 *iplocal)
/*
 *-----------------------------------------------------------------------
 *
 * store a new address in the shared lcm reference database.
 * This database is used to keep track of the elementary array pointers
 * which are shared between LCM and external objects (implemented in a OO
 * language such as C++/Boost, Python or Java). If a pointer is stored
 * in this database (using refpush), the free call on this pointer
 * is replaced by a refpop call.
 *
 * input parameters:
 *  iplist  : address of the object.
 *  iplocal : local reference.
 *
 *-----------------------------------------------------------------------
 */
{
   dbref* ipkeep = (*iplist)->global;
   int_32 n = ipkeep->nad;
   int_32 i;
   for (i = 0; i < n; ++i) {
      if (ipkeep->local[i] == iplocal) return;
   }
   if (ipkeep->nad + 1 > ipkeep->maxad) {
      /* increase the size of the database */
      int_32 **my_local;
      ipkeep->maxad += 50;
      my_local = (int_32 **) malloc((ipkeep->maxad)*sizeof(*my_local));
      for (i = 0; i < n; ++i) my_local[i]=ipkeep->local[i];
      if (n > 0) free(ipkeep->local);
      ipkeep->local=my_local;
   }
   ++ipkeep->nad;
   ipkeep->local[n] = iplocal;
   return;
}

int_32 refpop(lcm **iplist, int_32 *iplocal)
/*
 *-----------------------------------------------------------------------
 *
 * remove one address of the shared lcm reference database.
 *
 * input parameters:
 *  iplist  : address of the object.
 *  iplocal : local reference.
 *
 * output parameter:
 *  refpop : =0:the reference does exists; =1: does not exists.
 *
 *-----------------------------------------------------------------------
 */
{
   dbref* ipkeep = (*iplist)->global;
   int_32 n = ipkeep->nad;
   int_32 i;
   if (n == 0) return 1;
   for (i = 0; i < n; ++i) {
      if (ipkeep->local[i] == iplocal) return 0;
   }
   return 1;
}

void lcmkep(db0 *ipkeep, int_32 imode, lcm **iplist)
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
   char *nomsub="lcmkep";
   int_32 n = ipkeep->nad;
   if (imode == 1) {
      int_32 i;
      lcm **my_parray;
      if ((*iplist)->header != 100) {
         sprintf(AbortString,"%s: WRONG HEADER(1).",nomsub);
         xabort_c(AbortString);
      } else if (ipkeep->nad + 1 > ipkeep->maxad) {
         ipkeep->maxad += 500;
         my_parray = (lcm **) malloc((ipkeep->maxad)*sizeof(*my_parray));
         for (i = 0; i < n; ++i) my_parray[i]=ipkeep->idir[i];
         if (n > 0) free(ipkeep->idir);
         ipkeep->idir=my_parray;
      }
      ++ipkeep->nad;
      ipkeep->idir[n] = *iplist;
   } else if (imode == 2) {
      int_32 i;
      for (i = n; i >= 1; --i) {
         if (ipkeep->idir[i-1] == *iplist) {
            ipkeep->idir[i-1] = NULL;
            return;
         }
      }
      sprintf(AbortString,"%s: UNABLE TO FIND AN ADDRESS.",nomsub);
      xabort_c(AbortString);
   } else {
      sprintf(AbortString,"%s: INVALID VALUE OF IMODE.",nomsub);
      xabort_c(AbortString);
   }
   return;
}

void lcmop_c(lcm **iplist, char *namp, int_32 imp, int_32 medium, int_32 impx)
/*
 *-----------------------------------------------------------------------
 *
 * open an existing or create a new object.
 *
 * input parameters:
 *  iplist : address of the existing object if imp=1 or imp=2.
 *    namp : character name (null terminated) of the object if imp=0.
 *     imp : =0 to create a new object; =1 to modify an existing object;
 *           =2 to access an existing object in read-only mode.
 *  medium : =1 use memory; =2 use an xsm file.
 *    impx : if impx=0, we suppress printing on lcmop.
 *
 * output parameters:
 *  iplist : address of the new object if imp=0.
 *    namp : character name (null terminated) of the object if imp=1 or imp=2.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmop_c";
   char text12[13];
   blockb *my_blockb;
   db0 *my_db0;
   dbref *my_dbref;
   if (medium == 2) {
/*     USE A XSM FILE TO STORE INFORMATION */
      xsmop_c((xsm **)iplist, namp, imp, impx);
      return;
   } else if (medium != 1) {
      sprintf(AbortString,"%s: INVALID MEDIUM (%d).",nomsub,(int)medium);
      xabort_c(AbortString);
   } else if (imp < 0 || imp > 2) {
      sprintf(AbortString,"%s: INVALID ACTION (%d) ON XSM FILE '%s'.",nomsub,(int)imp,namp);
      xabort_c(AbortString);
   } else if (strlen(namp) > 72) {
      sprintf(AbortString,"%s: OBJECT NAME '%s' EXCEEDING 72 CHARACTERS.",nomsub,namp);
      xabort_c(AbortString);
   }
   if (imp == 0) {
      *iplist = (lcm *) malloc(sizeof(**iplist));
      my_db0 = (db0 *) malloc(sizeof(*my_db0));
      my_dbref = (dbref *) malloc(sizeof(*my_dbref));
      (*iplist)->header = 100;
      if (strcmp(namp," ") == 0) {
         strcpy((*iplist)->hname,"*TEMPORARY*");
      } else {
         strcpy((*iplist)->hname,namp);
      }
      (*iplist)->listlen = -1;
      (*iplist)->inext = NULL;
      (*iplist)->father = NULL;
      (*iplist)->ifdir = 0;
      (*iplist)->imode = 1;
      (*iplist)->imax = 0;
      (*iplist)->inref = 0;
      (*iplist)->icang = my_db0;
      (*iplist)->global = my_dbref;
      (*iplist)->hash = NULL;
      my_db0->nad = 0;
      my_db0->maxad = 0;
      my_db0->idir = NULL;
      my_dbref->nad = 0;
      my_dbref->maxad = 0;
      my_dbref->local = NULL;
   } else if ((*iplist)->header != 100) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' WITH ADDRESS =%ld HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname,(long)(*iplist));
      xabort_c(AbortString);
   } else if ((*iplist)->imode != 0) {
      if ((*iplist)->father == NULL) {
         strcpy(text12,"/");
      } else {
         lcm *my_father;
         my_father = (*iplist)->father;
         my_blockb = my_father->inext;
         if (my_blockb == NULL) {
            memcpy(text12,"            ",12);
         } else {
            strncpy(text12,(char*)my_blockb[(*iplist)->ifdir - 1].jcmt,12);
         }
         text12[12]='\0';
      }
      sprintf(AbortString,"%s: DIRECTORY '%s' IN THE OBJECT '%.45s' WITH ADDRESS =%ld IS ALREADY OPEN.",
              nomsub,text12,(*iplist)->hname,(long)(*iplist));
      xabort_c(AbortString);
   } else {
      int_32 i, n;
      db0 *ipkeep;
      strcpy(namp,(*iplist)->hname);
      (*iplist)->imode = imp;
      ipkeep = (*iplist)->icang;
      n = ipkeep->nad;
      if (n > 0) {
         for (i = 0; i < n; ++i) {
           if (ipkeep->idir[i] != NULL) (ipkeep->idir[i])->imode = imp;
         }
      }
   }
   if (impx > 0 && imp == 0) {
      printf("%s: OPEN A NEW OBJECT NAMED '%s' WITH ADDRESS = %ld.\n",
             nomsub,(*iplist)->hname,(long)(*iplist));
   } else if (impx > 0 && imp == 1) {
      printf("%s: MODIFY AN OBJECT NAMED '%s' WITH ADDRESS = %ld.\n",
             nomsub,(*iplist)->hname,(long)(*iplist));
   } else if (impx > 0 && imp == 2) {
      printf("%s: OPEN AN OBJECT NAMED '%s' WITH ADDRESS = %ld IN READ-ONLY MODE.\n",
             nomsub,(*iplist)->hname,(long)(*iplist));
   }
   return;
}

void lcmppd_c(lcm **iplist, const char *namp, int_32 ilong, int_32 itype, int_32 *iofdum)
/*
 *-----------------------------------------------------------------------
 *
 * add a new pointer entry in the table.
 *
 * input parameters:
 *  iplist : address of the object.
 *    namp : character*12 name of the current block.
 *   ilong : number of information elements stored in the current block.
 *   itype : type of information elements stored in the current block.
 *  iofdum : pointer of the first information element.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmppd_c";
   int_32 i, j, ipos, iref, *jofdum;
   int inamt[3];
   blockb *ipnode, *jpnode;
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->listlen >= 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS A LIST.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iofdum == NULL) {
      sprintf(AbortString,"%s: THE MALLOC POINTER OF NODE '%s' IS NOT SET IN THE OBJECT '%.60s'.",
              nomsub,namp,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (ilong <= 0) {
      sprintf(AbortString,"%s: INVALID LENGTH (%d) FOR NODE '%s' IN THE OBJECT '%.60s'.",
              nomsub,(int)ilong,namp,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsmppd_c((xsm **)iplist,namp,ilong,itype,iofdum);
      return;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->imode == 2) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS OPEN IN READ-ONLY MODE.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   /* SCAN THE NODAL TABLE AND INCLUDE THE NEW ENTRY. */
   ipnode = (*iplist)->inext;
   strncpy((char*)inamt,namp,12);
   ipos = abs(inamt[0] + inamt[1] * 3 + inamt[2] * 5) % lhash;
   if (ipnode == NULL) goto L10;
   for (i = (*iplist)->hash[ipos] - 1; i >= 0; --i) {
      if (ipnode[i].jcmt[0] == inamt[0]) {
         if (ipnode[i].jcmt[1] == inamt[1]) {
            if (ipnode[i].jcmt[2] == inamt[2]) {
               iref = i + 1;
               /* REMOVE THE OLD NODE. */
               jofdum = ipnode[i].jdata;
               if (jofdum != iofdum) {
                  if(refpop(iplist,jofdum)) free(jofdum); /* rlsara_c(jofdum); */
               }
               goto L20;
            }
         }
      }
   }
L10:
   iref = (*iplist)->inref + 1;
   if (iref > (*iplist)->imax) {
      /* INCREASE THE SIZE OF THE NODE TABLE. */
      jpnode = (blockb *) malloc(((*iplist)->imax + maxext) * sizeof(*jpnode));
      (*iplist)->inext = jpnode;
      if (ipnode != NULL) {
         for (i = 0; i < (*iplist)->imax; ++i) {
            jpnode[i].jdata = ipnode[i].jdata;
            jpnode[i].jjlon = ipnode[i].jjlon;
            jpnode[i].jjtyp = ipnode[i].jjtyp;
            for (j = 0; j < 4; ++j) jpnode[i].jidat[j] = ipnode[i].jidat[j];
            for (j = 0; j < 3; ++j) jpnode[i].jcmt[j] = ipnode[i].jcmt[j];
         }
         free(ipnode);
      } else {
         (*iplist)->hash = (int_32 *) malloc(lhash*sizeof(int_32));
         for (i = 0; i < lhash; ++i) (*iplist)->hash[i] = 0;
      }
      (*iplist)->imax += maxext;
      ipnode = jpnode;
   }
   (*iplist)->hash[ipos] = max(iref, (*iplist)->hash[ipos]);
   (*iplist)->inref = iref;
   for (j = 0; j < 3; ++j) ipnode[iref-1].jcmt[j] = inamt[j];

   /* STORE THE INFORMATION RELATIVE TO THE NEW INFORMATION ELEMENT. */
L20:
   (ipnode[iref-1]).jdata = iofdum;
   (ipnode[iref-1]).jjlon = ilong;
   (ipnode[iref-1]).jjtyp = itype;

   /* STORE THE FIRST AND LAST ELEMENTS FOR VALIDATION PURPOSE. */
   if (itype == 1 || itype == 2 || itype == 3 || itype == 5) {
      (ipnode[iref-1]).jidat[0] = iofdum[0];
      (ipnode[iref-1]).jidat[1] = iofdum[ilong-1];
   } else if (itype == 4 || itype == 6) {
      (ipnode[iref-1]).jjlon = 2*ilong;
      (ipnode[iref-1]).jidat[0] = iofdum[0];
      (ipnode[iref-1]).jidat[1] = iofdum[1];
      (ipnode[iref-1]).jidat[2] = iofdum[2*ilong-2];
      (ipnode[iref-1]).jidat[3] = iofdum[2*ilong-1];
   }
   return;
}

void lcmlen_c(lcm **iplist, const char *namp, int_32 *ilong, int_32 *itylcm)
/*
 *-----------------------------------------------------------------------
 *
 * return the length and type of a table entry.
 *
 * input parameters:
 *  iplist : address of the object.
 *    namp : character*12 name of the current block.
 *
 * output parameters:
 *  ilong  : number of information elements pointed by the lcm entry.
 *           ilong=0 is returned if the entry does not exists.
 *           ilong=-1 is returned for an associative table.
 *  itylcm : type of information elements pointed by the lcm entry.
 *           0: directory                1: integer
 *           2: single precision         3: character*4
 *           4: double precision         5: logical
 *           6: complex                 99: empty node
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmlen_c";
   blockb *ipnode;
   int_32 i, ipos;
   int inamt[3];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->listlen >= 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS A LIST.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsmlen_c((xsm **)iplist,namp,ilong,itylcm);
      return;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   ipnode = (*iplist)->inext;
   if (ipnode == NULL) goto L10;
   strncpy((char*)inamt,namp,12);
   ipos = abs(inamt[0] + inamt[1] * 3 + inamt[2] * 5) % lhash;
   for (i = (*iplist)->hash[ipos] - 1; i >= 0; --i) {
      if (ipnode[i].jcmt[0] == inamt[0]) {
         if (ipnode[i].jcmt[1] == inamt[1]) {
            if (ipnode[i].jcmt[2] == inamt[2]) {
               *ilong = ipnode[i].jjlon;
               *itylcm = ipnode[i].jjtyp;
               if (*itylcm == 4 || *itylcm == 6) *ilong=*ilong/2;
               return;
            }
         }
      }
   }
L10:
   *ilong = 0;
   *itylcm = 99;
   return;
}

void lcminf_c(lcm **iplist, char *namlcm, char *nammy, int_32 *empty, int_32 *ilong, int_32 *lcml,
     int_32 *access)
/*
 *-----------------------------------------------------------------------
 *
 * find general information about an associative table or list.
 *
 * input parameters:
 *  iplist : address of the object.
 *
 * output parameter:
 *  namlcm : character name (null terminated) of the object.
 *   nammy : character name (null terminated) of the active directory.
 *   empty : =.true. if the active directory is empty.
 *   ilong : =-1: for a table; >0: number of list items.
 *    lcml : =.true.: memory used; =.false.: xsm file used.
 *  access : type of access. =0: object closed (lcm only); =1: object
 *           open for modification; =2: object in read-only mode.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcminf_c";
   blockb *my_blockb;
   if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      *lcml = 0;
      xsminf_c((xsm **)iplist,namlcm,nammy,empty,ilong,access);
      return;
   } else if ((*iplist)->header != 100) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   *lcml = 1;
   *access = (*iplist)->imode;
   if ((*iplist)->father == NULL) {
      strcpy(nammy,"/");
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      my_blockb = my_father->inext;
      if (my_blockb == NULL) {
         strncpy(nammy,"            ",12);
      } else {
         strncpy(nammy,(char*)my_blockb[(*iplist)->ifdir - 1].jcmt,12);
      }
      nammy[12]='\0';
   }
   strcpy(namlcm,(*iplist)->hname);
   *ilong = (*iplist)->listlen;
   my_blockb = (*iplist)->inext;
   *empty = (my_blockb == NULL);
   if ((!*empty) && (*ilong == -1)) {
      char namp[13];
      int_32 iref, i;
      iref = 0;
L10:
      ++iref;
      if (iref > (*iplist)->inref) {
         *empty = 1;
         return;
      }
      strncpy(namp,(char*)my_blockb[iref-1].jcmt,12);
      namp[12]=' ';
      for(i=12; i>0; i--) {
         if (namp[i] != ' ') break;
         namp[i]='\0';
      }
      if (strcmp(namp," ") == 0) goto L10;
   }
   return;
}

void lcmnxt_c(lcm **iplist,char *namp)
/*
 *-----------------------------------------------------------------------
 *
 * input parameters:
 *  iplist : address of the object.
 *    namp : character*12 name of a block. if namp=' ' at input, find
 *           any name for any block stored in this directory.
 *
 * output parameter:
 *    namp : character*12 name of the next block.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmnxt_c";
   char nammy[13];
   blockb *ipnode;
   int_32 i, ipos, iref, lcheck;
   int inamt[3];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->listlen >= 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS A LIST.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsmnxt_c((xsm **)iplist,namp);
      return;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   ipnode = (*iplist)->inext;
   if (ipnode == NULL) {
      if ((*iplist)->father == NULL) {
         strcpy(nammy,"/");
      } else {
         lcm *my_father;
         my_father = (*iplist)->father;
         ipnode = my_father->inext;
         if (ipnode == NULL) {
            memcpy(nammy,"            ",12);
         } else {
            strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
         }
         nammy[12]='\0';
      }
      sprintf(AbortString,"%s: EMPTY DIRECTORY '%s' IN THE OBJECT '%.60s' (1).",
              nomsub,nammy,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (strcmp(namp," ") == 0) {
      iref = 0;
L10:
      ++iref;
      if (iref > (*iplist)->inref) {
         if ((*iplist)->father == NULL) {
            strcpy(nammy,"/");
         } else {
            lcm *my_father;
            my_father = (*iplist)->father;
            ipnode = my_father->inext;
            if (ipnode == NULL) {
               memcpy(nammy,"            ",12);
            } else {
               strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
            }
            nammy[12]='\0';
         }
         sprintf(AbortString,"%s: EMPTY DIRECTORY '%s' IN THE OBJECT '%.60s' (2).",
                 nomsub,nammy,(*iplist)->hname);
         xabort_c(AbortString);
      }
      strncpy(namp,(char*)ipnode[iref-1].jcmt,12);
      namp[12]=' ';
      for(i=12; i>0; i--) {
         if (namp[i] != ' ') break;
         namp[i]='\0';
      }
      if (strcmp(namp," ") == 0) goto L10;
      return;
   }
   strncpy((char*)inamt,namp,12);
   ipos = abs(inamt[0] + inamt[1] * 3 + inamt[2] * 5) % lhash;
   for (i = (*iplist)->hash[ipos] - 1; i >= 0; --i) {
      iref = i+1;
      if (ipnode[i].jcmt[0] == inamt[0]) {
         if (ipnode[i].jcmt[1] == inamt[1]) {
            if (ipnode[i].jcmt[2] == inamt[2]) {
               lcheck = 0;
L20:
               if (iref < (*iplist)->inref) {
                  ++iref;
               } else {
                  if (lcheck == 1) {
                     sprintf(AbortString,"%s: INFINITE LOOPING.",nomsub);
                     xabort_c(AbortString);
                  }
                  iref = 1;
                  lcheck = 1;
               }
               strncpy(namp,(char*)ipnode[iref-1].jcmt,12);
               namp[12]=' ';
               for(i=12; i>0; i--) {
                  if (namp[i] != ' ') break;
                  namp[i]='\0';
               }
               if (strcmp(namp," ") == 0) goto L20;
               return;
            }
         }
      }
   }
   if ((*iplist)->father == NULL) {
      strcpy(nammy,"/");
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode == NULL) {
         memcpy(nammy,"            ",12);
      } else {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
      }
      nammy[12]='\0';
   }
   sprintf(AbortString,"%s: UNABLE TO FIND BLOCK '%s' INTO DIRECTORY '%s' IN THE OBJECT '%.50s'.",
           nomsub,namp,nammy,(*iplist)->hname);
   xabort_c(AbortString);
}

void lcmgpd_c(lcm **iplist, const char *namp, int_32 **iofdum)
/*
 *-----------------------------------------------------------------------
 *
 * get a malloc pointer for an entry in the table.
 *
 * input parameters:
 *  iplist : address of the object.
 *    namp : character*12 name of the current block.
 *
 * output parameter:
 *  iofdum : malloc pointer to the lcm entry named namp.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmgpd_c";
   char nammy[13];
   blockb *ipnode;
   int_32 i, ipos;
   int inamt[3];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->listlen >= 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS A LIST.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsmgpd_c((xsm **)iplist,namp,iofdum);
      return;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   ipnode = (*iplist)->inext;
   if (ipnode == NULL) goto L10;
   strncpy((char*)inamt,namp,12);
   ipos = abs(inamt[0] + inamt[1] * 3 + inamt[2] * 5) % lhash;
   for (i = (*iplist)->hash[ipos] - 1; i >= 0; --i) {
      if (ipnode[i].jcmt[0] == inamt[0]) {
         if (ipnode[i].jcmt[1] == inamt[1]) {
            if (ipnode[i].jcmt[2] == inamt[2]) {
               *iofdum = ipnode[i].jdata;
               return;
            }
         }
      }
   }
L10:
   if ((*iplist)->father == NULL) {
      strcpy(nammy,"/");
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode == NULL) {
         memcpy(nammy,"            ",12);
      } else {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
      }
      nammy[12]='\0';
   }
   sprintf(AbortString,"%s: UNABLE TO FIND BLOCK '%s' INTO DIRECTORY '%s' IN THE OBJECT '%.50s'.",
           nomsub,namp,nammy,(*iplist)->hname);
   xabort_c(AbortString);
}

void lcmget_c(lcm **iplist, const char *namp, int_32 *idata)
/*
 *-----------------------------------------------------------------------
 *
 * copy a block of data from a table into memory.
 *
 * input parameters:
 *  iplist : address of the object.
 *    namp : character*12 name of the current block.
 *
 * output parameter:
 *   idata : information elements. dimension idata1(ilong) where ilong
 *           is the number of information elements.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmget_c";
   char nammy[13];
   blockb *ipnode;
   int_32 i, j, ipos;
   int inamt[3];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->listlen >= 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS A LIST.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsmget_c((xsm **)iplist,namp,idata);
      return;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   ipnode = (*iplist)->inext;
   if (ipnode == NULL) goto L10;
   strncpy((char*)inamt,namp,12);
   ipos = abs(inamt[0] + inamt[1] * 3 + inamt[2] * 5) % lhash;
   for (i = (*iplist)->hash[ipos] - 1; i >= 0; --i) {
      if (ipnode[i].jcmt[0] == inamt[0]) {
         if (ipnode[i].jcmt[1] == inamt[1]) {
            if (ipnode[i].jcmt[2] == inamt[2]) {
               for (j = 0; j < ipnode[i].jjlon; ++j) idata[j] = ipnode[i].jdata[j];
               return;
            }
         }
      }
   }
L10:
   if ((*iplist)->father == NULL) {
      strcpy(nammy,"/");
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode == NULL) {
         memcpy(nammy,"            ",12);
      } else {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
      }
      nammy[12]='\0';
   }
   sprintf(AbortString,"%s: UNABLE TO FIND BLOCK '%s' INTO DIRECTORY '%s' IN THE OBJECT '%.50s'.",
           nomsub,namp,nammy,(*iplist)->hname);
   xabort_c(AbortString);
}

void lcmval_part2(int_32 ilong,lcm *iplist);

void lcmval_part1(lcm *iplist,const char *namp)
/* ASSOCIATIVE TABLE VALIDATION. */
{
   char *nomsub="lcmval_part1";
   char namt[13];
   int_32 iref, ilong, itylcm, lerr1, lerr2;
   blockb *inode;
   inode = iplist->inext;
   for (iref = 0; iref < iplist->inref; ++iref) {
      strncpy(namt,(char*)inode[iref].jcmt,12);
      namt[12]='\0';
      if ( ((strcmp(namp," ") == 0 || strcmp(namp,namt) == 0)) &&
           (strcmp(namt," ") != 0) ) {
         ilong = inode[iref].jjlon;
         itylcm = inode[iref].jjtyp;
         lerr1 = 0;
         lerr2 = 0;
         if (itylcm == 0 && ilong == -1) {
            /* ASSOCIATIVE TABLE. */
            lcmval_part1((lcm *)inode[iref].jdata," ");
         } else if (itylcm == 10) {
            /* LIST. */
            lcmval_part2(ilong,(lcm *)inode[iref].jdata);
         } else if (itylcm == 1 || itylcm == 2 || itylcm == 3 || itylcm == 5) {
            lerr1 = inode[iref].jidat[0] != inode[iref].jdata[0];
            lerr2 = inode[iref].jidat[1] != inode[iref].jdata[ilong - 1];
         } else if (itylcm == 4 || itylcm == 6) {
            lerr1 = (inode[iref].jidat[0] != inode[iref].jdata[0]) ||
                    (inode[iref].jidat[1] != inode[iref].jdata[1]);
            lerr2 = (inode[iref].jidat[2] != inode[iref].jdata[ilong - 2]) ||
                    (inode[iref].jidat[3] != inode[iref].jdata[ilong - 1]);
         }
         if (lerr1) {
            sprintf(AbortString,"%s: BLOCK '%s' OF THE OBJECT '%.50s' HAS BEEN OVERWRITTEN (1).",
                    nomsub,namt,iplist->hname);
            xabort_c(AbortString);
         } else if (lerr2) {
            sprintf(AbortString,"%s: BLOCK '%s' OF THE OBJECT '%.50s' HAS BEEN OVERWRITTEN (2).",
                    nomsub,namt,iplist->hname);
            xabort_c(AbortString);
         }
      }
   }
   return;
}

void lcmval_part2(int_32 kjlon,lcm *iplist)
/* LIST VALIDATION. */
{
   char *nomsub="lcmval_part2";
   int_32 ilong, itylcm, lerr1, lerr2, ivec;
   blockb *knode;
   for (ivec = 0; ivec < kjlon; ++ivec) {
      knode = iplist[ivec].inext;
      if (knode) {
         ilong = knode[0].jjlon;
         itylcm = knode[0].jjtyp;
         lerr1 = 0;
         lerr2 = 0;
         if (itylcm == 0 && ilong == -1) {
            /* ASSOCIATIVE TABLE. */
            lcmval_part1((lcm *)knode[0].jdata," ");
         } else if (itylcm == 10) {
            /* LIST. */
            lcmval_part2(ilong,(lcm *)knode[0].jdata);
         } else if (itylcm == 1 || itylcm == 2 || itylcm == 3 || itylcm == 5) {
            lerr1 = knode[0].jidat[0] != knode[0].jdata[0];
            lerr2 = knode[0].jidat[1] != knode[0].jdata[ilong - 1];
         } else if (itylcm == 4 || itylcm == 6) {
            lerr1 = (knode[0].jidat[0] != knode[0].jdata[0]) ||
                    (knode[0].jidat[1] != knode[0].jdata[1]);
            lerr2 = (knode[0].jidat[2] != knode[0].jdata[ilong - 2]) ||
                    (knode[0].jidat[3] != knode[0].jdata[ilong - 1]);
         }
         if (lerr1) {
            sprintf(AbortString,"%s: LIST ELEMENT %d OF THE OBJECT '%.50s' HAS BEEN OVERWRITTEN (1).",
                    nomsub,(int)ivec,iplist->hname);
            xabort_c(AbortString);
         } else if (lerr2) {
            sprintf(AbortString,"%s: LIST ELEMENT %d OF THE OBJECT '%.50s' HAS BEEN OVERWRITTEN (2).",
                    nomsub,(int)ivec,iplist->hname);
            xabort_c(AbortString);
         }
      }
   }
   return;
}

void lcmval_c(lcm **iplist,const char *namp)
{
   char *nomsub="lcmval_c";
   if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      return;
   } else if ((*iplist)->header != 100) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   if ((*iplist)->listlen == -1) {
      lcmval_part1(*iplist,namp);
   } else {
      lcmval_part2((*iplist)->listlen,*iplist);
   }
   return;
}

void lcmcl_part2(int_32 kjlon, lcm *iplist);

void lcmcl_part1(lcm *iplist)
/* ASSOCIATIVE TABLE DESTRUCTION. */
{
   int_32 iref, ilong, itylcm;
   blockb *inode;
   lcm *kplist;
   inode = iplist->inext;
   for (iref = 0; iref < iplist->inref; ++iref) {
      ilong = inode[iref].jjlon;
      itylcm = inode[iref].jjtyp;
      if (itylcm == 0 && ilong == -1) {
         /* ASSOCIATIVE TABLE. */
         kplist = (lcm *)inode[iref].jdata;
         lcmcl_part1(kplist);
         lcmkep(iplist->icang,c__2,&kplist);
         free(inode[iref].jdata);
      } else if (itylcm == 10) {
         /* LIST. */
         kplist = (lcm *)inode[iref].jdata;
         lcmcl_part2(ilong, kplist);
         lcmkep(iplist->icang,c__2,&kplist);
         free(inode[iref].jdata);
      } else if (itylcm != 99) {
         if(refpop(&iplist,inode[iref].jdata)) free(inode[iref].jdata); /* rlsara_c() */
      }
   }
   if (inode != NULL) free(inode);
   if (iplist->hash != NULL) free(iplist->hash);
   return;
}

void lcmcl_part2(int_32 kjlon, lcm *iplist)
/* LIST DESTRUCTION. */
{
   int_32 ilong, itylcm, ivec;
   lcm *kplist;
   blockb *knode;
   for (ivec = 0; ivec < kjlon; ++ivec) {
      knode = iplist[ivec].inext;
      if (knode) {
         ilong = knode[0].jjlon;
         itylcm = knode[0].jjtyp;
         if (itylcm == 0 && ilong == -1) {
            /* ASSOCIATIVE TABLE. */
            kplist = (lcm *)knode[0].jdata;
            lcmcl_part1(kplist);
            lcmkep(iplist[ivec].icang,c__2,&kplist);
            free(knode[0].jdata);
         } else if (itylcm == 10) {
            /* LIST. */
            kplist = (lcm *)knode[0].jdata;
            lcmcl_part2(ilong,kplist);
            lcmkep(iplist[ivec].icang,c__2,&kplist);
            free(knode[0].jdata);
         } else if (itylcm != 99) {
            if(refpop(&iplist,knode[0].jdata)) free(knode[0].jdata); /* rlsara_c() */
         }
         free(knode);
      }
   }
   return;
}

void lcmcl_c(lcm **iplist,int_32 iact)
/*
 *-----------------------------------------------------------------------
 *
 * close, destroy or erase an object with validation.
 *
 * input parameters:
 *  iplist : address of the existing object.
 *    iact : =1 to close; =2 to destroy; =3 to erase and close.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmcl_c";
   db0 *ipkeep;
   dbref *ipkref;
   if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      if (iact == 3) {
         sprintf(AbortString,"%s: THE XSM FILE '%s' CANNOT BE ERASED.",
                 nomsub,(*iplist)->hname);
         xabort_c(AbortString);
      }
      xsmcl_c((xsm **)iplist, iact);
      return;
   } else if ((*iplist)->header != 100) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iact != 1 && iact != 2 && iact != 3) {
      sprintf(AbortString,"%s: INVALID ACTION (%d) ON THE OBJECT '%.60s'.",
              nomsub,(int)iact,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS ALREADY CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->father != NULL) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS NOT ROOT.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   lcmval_c(iplist, " ");
   if (iact == 2 || iact == 3) {
      /* DESTROY OR ERASE THE OBJECT. */
      if ((*iplist)->imode == 2) {
         sprintf(AbortString,"%s: CANNOT DESTROY OR ERASE THE OBJECT '%.60s' OPEN IN READ-ONLY MODE.",
                 nomsub,(*iplist)->hname);
         xabort_c(AbortString);
      }
      /* RECURSIVE DESTRUCTION OF THE OBJECT CONTENT WITH ADDRESS IPLIST. */
      if ((*iplist)->listlen == -1) {
         lcmcl_part1(*iplist);
      } else {
         lcmcl_part2((*iplist)->listlen,*iplist);
      }
   } else {
      int_32 i, n;
      ipkeep = (*iplist)->icang;
      n = ipkeep->nad;
      if (n > 0) {
         for (i = 0; i < n; ++i) {
           if (ipkeep->idir[i] != NULL) (ipkeep->idir[i])->imode = 0;
         }
      }
   }
   if ((*iplist)->father == NULL && iact >= 2) {
      /*  DESTROY THE TABLE. */
      int_32 i, n;
      ipkeep = (*iplist)->icang;
      n = ipkeep->nad;
      if (n > 0) {
         for (i = 0; i < n; ++i) {
           if (ipkeep->idir[i] != NULL) free(ipkeep->idir[i]);
         }
      }
      free(ipkeep);
      ipkref = (*iplist)->global;
      n = ipkref->nad;
      if (n > 0) free(ipkref->local);
      free(ipkref);
      if (iact == 2) {
         free(*iplist);
         *iplist = NULL;
      } else if (iact == 3) {
         /* ERASE THE TABLE. */
         db0 *my_db0;
         dbref *my_dbref;
         my_db0 = (db0 *) malloc(sizeof(*my_db0));
         my_dbref = (dbref *) malloc(sizeof(*my_dbref));
         (*iplist)->inext = NULL;
         (*iplist)->imode = 0;
         (*iplist)->imax = 0;
         (*iplist)->inref = 0;
         (*iplist)->icang = my_db0;
         (*iplist)->global = my_dbref;
         (*iplist)->hash = NULL;
         my_db0->nad = 0;
         my_db0->maxad=0;
         my_db0->idir=NULL;
         my_dbref->nad = 0;
         my_dbref->maxad=0;
         my_dbref->local=NULL;
      }
   } else {
      (*iplist)->imode = 0;
   }
   return;
}

void lcmsix_c(lcm **iplist,const char *namp,int_32 iact)
/*
 *-----------------------------------------------------------------------
 *
 * move in the hierarchical structure of a node table.
 *
 * input parameters:
 *  iplist : address of the table.
 *    namp : character*12 name of the directory if iact=1. not used if
 *           iact.ne.1.
 *    iact : type of movement in the hierarchical structure.
 *           0: return to the root directory;
 *           1: move to a son directory;
 *           2: move back to the parent directory.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmsix_c";
   blockb *ipnode, *jpnode;
   int_32 i, j, ipos, mode, iref;
   int inamt[3];
   lcm *jplist;
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iact < 0 || iact > 2) {
      sprintf(AbortString,"%s: INVALID ACTION (%d) ON THE OBJECT '%.60s'.",
              nomsub,(int)iact,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->listlen >= 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS A LIST.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsmsix_c((xsm **)iplist,namp,iact);
      return;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   if (iact == 0 && (*iplist)->father == NULL) return;
   ipnode = (*iplist)->inext;
   mode = (*iplist)->imode;
   if (iact == 1) {
      /* MOVE TO A SON DIRECTORY. */
      /* CHECK IF THE DIRECTORY EXISTS IN THE NODE TABLE. */
      strncpy((char*)inamt,namp,12);
      ipos = abs(inamt[0] + inamt[1] * 3 + inamt[2] * 5) % lhash;
      if (ipnode == NULL) goto L10;
      for (i = (*iplist)->hash[ipos] - 1; i >= 0; --i) {
         if (ipnode[i].jcmt[0] == inamt[0]) {
            if (ipnode[i].jcmt[1] == inamt[1]) {
               if (ipnode[i].jcmt[2] == inamt[2]) {
                  iref = i+1;
                  goto L20;
               }
            }
         }
      }
      /* THE DIRECTORY DOES NOT EXISTS IN THE NODE TABLE. */
L10:
      if ((*iplist)->imode == 2) {
         char nammy[13];
         if ((*iplist)->father == NULL) {
            strcpy(nammy,"/");
         } else {
            lcm *my_father;
            my_father = (*iplist)->father;
            ipnode = my_father->inext;
            if (ipnode == NULL) {
               memcpy(nammy,"            ",12);
            } else {
               strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
            }
            nammy[12]='\0';
         }
         sprintf(AbortString,"%s: UNABLE TO CREATE DIRECTORY '%s' FROM DIRECTORY '%s' IN READ-ONLY OBJECT '%.35s'.",
                 nomsub,namp,nammy,(*iplist)->hname);
         xabort_c(AbortString);
      }
      /* CREATE A NEW NODE TABLE. */
      iref = (*iplist)->inref + 1;
      if (iref > (*iplist)->imax) {
         /* INCREASE THE SIZE OF THE NODE TABLE. */
         jpnode = (blockb *) malloc(((*iplist)->imax + maxext) * sizeof(*jpnode));
         (*iplist)->inext = jpnode;
         if (ipnode != NULL) {
            for (i = 0; i < (*iplist)->imax; ++i) {
               jpnode[i].jdata = ipnode[i].jdata;
               jpnode[i].jjlon = ipnode[i].jjlon;
               jpnode[i].jjtyp = ipnode[i].jjtyp;
               for (j = 0; j < 4; ++j) jpnode[i].jidat[j] = ipnode[i].jidat[j];
               for (j = 0; j < 3; ++j) jpnode[i].jcmt[j] = ipnode[i].jcmt[j];
            }
            free(ipnode);
         } else {
            (*iplist)->hash = (int_32 *) malloc(lhash*sizeof(int_32));
            for (i = 0; i < lhash; ++i) (*iplist)->hash[i] = 0;
         }
         (*iplist)->imax += maxext;
         ipnode = jpnode;
      }
      (*iplist)->hash[ipos] = max(iref, (*iplist)->hash[ipos]);
      (*iplist)->inref = iref;
      for (j = 0; j < 3; ++j) ipnode[iref-1].jcmt[j]=inamt[j];
      (ipnode[iref-1]).jjlon = -1;
      (ipnode[iref-1]).jjtyp = 0;

      jplist = (lcm *) malloc(sizeof(*jplist));
      jplist->header = (*iplist)->header;
      strcpy(jplist->hname, (*iplist)->hname);
      jplist->listlen = -1;
      jplist->inext = NULL;
      jplist->father = *iplist;
      jplist->ifdir = iref;
      jplist->imode = 0;
      jplist->imax = 0;
      jplist->inref = 0;
      jplist->icang = (*iplist)->icang;
      jplist->global = (*iplist)->global;
      jplist->hash = NULL;
      lcmkep(jplist->icang, c__1, &jplist);
      (ipnode[iref-1]).jdata = (int_32 *)jplist;

      /* SWITCH TO THE SON DIRECTORY. */
L20:
      if ((ipnode[iref-1]).jjlon != -1 || (ipnode[iref-1]).jjtyp != 0) {
         sprintf(AbortString,"%s: '%s' IS AN INVALID DIRECTORY TYPE. OBJECT='%s'.",
                 nomsub,namp,(*iplist)->hname);
         xabort_c(AbortString);
      }
      *iplist = (lcm *)(ipnode[iref-1]).jdata;
      (*iplist)->imode = mode;
   } else if (iact == 0 || iact == 2) {
      /* MOVE BACK TO THE ROOT OR PARENT DIRECTORY. */
      if ((*iplist)->father == NULL) {
         sprintf(AbortString,"%s: THE OBJECT '%.60s' IS ON ROOT DIRECTORY.",
                 nomsub,(*iplist)->hname);
         xabort_c(AbortString);
      }
L30:
      jplist = *iplist;
      *iplist = jplist->father;
      if (iact == 0 && (*iplist)->father != NULL) goto L30;
   }
   return;
}

lcm * lcmdid_c(lcm **iplist,const char *namp)
/*
 *-----------------------------------------------------------------------
 *
 * create/access a son table in a father table.
 *
 * input parameters:
 *  iplist : address of the father table.
 *    namp : character*12 name of the son table.
 *
 * output parameter:
 *  lcmdid_c : address of the son table.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmdid_c";
   blockb *ipnode, *jpnode;
   int_32 i, j, ipos, mode, iref;
   int inamt[3];
   lcm *jplist;
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->listlen >= 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS A LIST.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *kplist;
      xsmdid_c((xsm **)iplist,namp,&kplist);
      return (lcm*)kplist;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->imode == 2) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS OPEN IN READ-ONLY MODE.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   ipnode = (*iplist)->inext;
   mode = (*iplist)->imode;
   strncpy((char*)inamt,namp,12);
   ipos = abs(inamt[0] + inamt[1] * 3 + inamt[2] * 5) % lhash;
   if (ipnode == NULL) goto L10;
   for (i = (*iplist)->hash[ipos] - 1; i >= 0; --i) {
      if (ipnode[i].jcmt[0] == inamt[0]) {
         if (ipnode[i].jcmt[1] == inamt[1]) {
            if (ipnode[i].jcmt[2] == inamt[2]) {
               iref = i+1;
               goto L20;
            }
         }
      }
   }
   /* THE DIRECTORY DOES NOT EXISTS IN THE NODE TABLE. */
L10:
   /* CREATE A NEW NODE TABLE. */
   iref = (*iplist)->inref + 1;
   if (iref > (*iplist)->imax) {
      /* INCREASE THE SIZE OF THE NODE TABLE. */
      jpnode = (blockb *) malloc(((*iplist)->imax + maxext) * sizeof(*jpnode));
      (*iplist)->inext = jpnode;
      if (ipnode != NULL) {
         for (i = 0; i < (*iplist)->imax; ++i) {
            jpnode[i].jdata = ipnode[i].jdata;
            jpnode[i].jjlon = ipnode[i].jjlon;
            jpnode[i].jjtyp = ipnode[i].jjtyp;
            for (j = 0; j < 4; ++j) jpnode[i].jidat[j] = ipnode[i].jidat[j];
            for (j = 0; j < 3; ++j) jpnode[i].jcmt[j] = ipnode[i].jcmt[j];
         }
         free(ipnode);
      } else {
         (*iplist)->hash = (int_32 *) malloc(lhash*sizeof(int_32));
         for (i = 0; i < lhash; ++i) (*iplist)->hash[i] = 0;
      }
      (*iplist)->imax += maxext;
      ipnode = jpnode;
   }
   (*iplist)->hash[ipos] = max(iref, (*iplist)->hash[ipos]);
   (*iplist)->inref = iref;
   for (j = 0; j < 3; ++j) ipnode[iref-1].jcmt[j]=inamt[j];
   (ipnode[iref-1]).jjlon = -1;
   (ipnode[iref-1]).jjtyp = 0;

   jplist = (lcm *) malloc(sizeof(*jplist));
   jplist->header = (*iplist)->header;
   strcpy(jplist->hname, (*iplist)->hname);
   jplist->listlen = -1;
   jplist->inext = NULL;
   jplist->father = *iplist;
   jplist->ifdir = iref;
   jplist->imode = 0;
   jplist->imax = 0;
   jplist->inref = 0;
   jplist->icang = (*iplist)->icang;
   jplist->global = (*iplist)->global;
   jplist->hash = NULL;
   lcmkep(jplist->icang, c__1, &jplist);
   (ipnode[iref-1]).jdata = (int_32 *)jplist;

   /* SWITCH TO THE SON DIRECTORY. */
L20:
   if ((ipnode[iref-1]).jjlon != -1 || (ipnode[iref-1]).jjtyp != 0) {
      sprintf(AbortString,"%s: '%s' IS AN INVALID DIRECTORY TYPE. OBJECT='%s'.",
              nomsub,namp,(*iplist)->hname);
      xabort_c(AbortString);
   }
   jplist = (lcm *)(ipnode[iref-1]).jdata;
   jplist->imode = mode;
   return jplist;
}

lcm * lcmlid_c(lcm **iplist,const char *namp,int_32 ilong)
/*
 *-----------------------------------------------------------------------
 *
 * create/access the hierarchical structure of a list in a father table.
 *
 * input parameters:
 *  iplist : address of the father table.
 *    namp : character*12 name of the list.
 *   ilong : dimension of the list.
 *
 * output parameter:
 *  lcmlid_c : address of the list.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmlid_c";
   blockb *ipnode, *jpnode;
   int_32 i, j, ipos, lenold, ityold, mode, iref=0;
   int inamt[3];
   lcm *jplist;
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->listlen >= 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS A LIST.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsmlid_c((xsm **)iplist,namp,ilong,(xsm **)(&jplist));
      return jplist;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->imode == 2) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS OPEN IN READ-ONLY MODE.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (ilong <= 0) {
      sprintf(AbortString,"%s: INVALID LENGTH (%d) FOR NODE '%s' IN THE OBJECT '%.60s'.",
              nomsub,(int)ilong,namp,(*iplist)->hname);
      xabort_c(AbortString);
   }
   ipnode = (*iplist)->inext;
   mode = (*iplist)->imode;
   lenold = 0;
   ityold = 10;
   strncpy((char*)inamt,namp,12);
   ipos = abs(inamt[0] + inamt[1] * 3 + inamt[2] * 5) % lhash;
   if (ipnode == NULL) goto L10;
   for (i = (*iplist)->hash[ipos] - 1; i >= 0; --i) {
      if (ipnode[i].jcmt[0] == inamt[0]) {
         if (ipnode[i].jcmt[1] == inamt[1]) {
            if (ipnode[i].jcmt[2] == inamt[2]) {
               iref = i+1;
               lenold = (ipnode[i]).jjlon;
               ityold = (ipnode[i]).jjtyp;
               goto L10;
            }
         }
      }
   }
L10:
   if (ityold != 10) {
      sprintf(AbortString,"%s: '%s' IS AN INVALID LIST TYPE. OBJECT='%s'.",
              nomsub,namp,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (lenold != 0 && lenold > ilong) {
      ilong = lenold;
   }
   if (lenold == 0) {
      /* CREATE A NEW NODE TABLE. */
      iref = (*iplist)->inref + 1;
      if (iref > (*iplist)->imax) {
         /* INCREASE THE SIZE OF THE NODE TABLE. */
         jpnode = (blockb *) malloc(((*iplist)->imax + maxext) * sizeof(*jpnode));
         (*iplist)->inext = jpnode;
         if (ipnode != NULL) {
            for (i = 0; i < (*iplist)->imax; ++i) {
               jpnode[i].jdata = ipnode[i].jdata;
               jpnode[i].jjlon = ipnode[i].jjlon;
               jpnode[i].jjtyp = ipnode[i].jjtyp;
               for (j = 0; j < 4; ++j) jpnode[i].jidat[j] = ipnode[i].jidat[j];
               for (j = 0; j < 3; ++j) jpnode[i].jcmt[j] = ipnode[i].jcmt[j];
            }
            free(ipnode);
         } else {
            (*iplist)->hash = (int_32 *) malloc(lhash*sizeof(int_32));
            for (i = 0; i < lhash; ++i) (*iplist)->hash[i] = 0;
         }
         (*iplist)->imax += maxext;
         ipnode = jpnode;
      }
      (*iplist)->hash[ipos] = max(iref, (*iplist)->hash[ipos]);
      (*iplist)->inref = iref;
      for (j = 0; j < 3; ++j) ipnode[iref-1].jcmt[j] = inamt[j];
      (ipnode[iref-1]).jjtyp = 10;
   }
   if (ilong != lenold) {
      lcm *iofset, *iofold;
      (ipnode[iref-1]).jjlon = ilong;
      iofset = (lcm *) malloc(ilong*sizeof(*iofset));
      for (i = 0; i < ilong; ++i) {
         if (i < lenold) {
            iofold = (lcm *)(ipnode[iref-1]).jdata;
            iofset[i].header = iofold[i].header;
            strcpy(iofset[i].hname, iofold[i].hname);
            iofset[i].listlen = 0;
            iofset[i].inext = iofold[i].inext;
            iofset[i].father = iofold[i].father;
            iofset[i].ifdir = iofold[i].ifdir;
            iofset[i].imode = iofold[i].imode;
            iofset[i].imax = iofold[i].imax;
            iofset[i].inref = iofold[i].inref;
            iofset[i].icang = iofold[i].icang;
            iofset[i].global = iofold[i].global;
            iofset[i].hash = iofold[i].hash;
            /* PUT THE OLD OBJECT IN READ-ONLY MODE */
            iofold[i].imode = 2;
         } else {
            iofset[i].header = (*iplist)->header;
            strcpy(iofset[i].hname, (*iplist)->hname);
            iofset[i].listlen = 0;
            iofset[i].inext = NULL;
            iofset[i].father = *iplist;
            iofset[i].ifdir = iref;
            iofset[i].imode = 0;
            iofset[i].imax = 0;
            iofset[i].inref = 0;
            iofset[i].icang = (*iplist)->icang;
            iofset[i].global = (*iplist)->global;
            iofset[i].hash = NULL;
         }
      }
      iofset[0].listlen = ilong;
      lcmkep(iofset->icang, c__1, &iofset);
      (ipnode[iref-1]).jdata = (int_32 *)iofset;
   }
   /* SWITCH TO THE SON LIST. */
   jplist = (lcm *)(ipnode[iref-1]).jdata;
   for (i = 0; i < ilong; ++i) jplist[i].imode = mode;
   return jplist;
}

lcm * lcmgid_c(lcm **iplist, const char *namp)
/*
 *-----------------------------------------------------------------------
 *
 * get the address of a table or of a list located in a father table.
 *
 * input parameters:
 *  iplist : address of the father table.
 *    namp : character*12 name of the son table or list.
 *
 * output parameter:
 *  lcmgid_c : address of the table or of the list named namp.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmgid_c";
   char nammy[13];
   blockb *ipnode;
   int_32 i, ipos;
   lcm *jplist;
   int inamt[3];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->listlen >= 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS A LIST.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *jplist;
      int_32 ilong, itylcm;
      xsmlen_c((xsm**)iplist, namp, &ilong, &itylcm);
      if (ilong == -1) {
         xsmdid_c((xsm **)iplist,namp,&jplist);
      } else {
         xsmlid_c((xsm **)iplist,namp,ilong,&jplist);
      }
      return (lcm*)jplist;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   ipnode = (*iplist)->inext;
   if (ipnode == NULL) goto L10;
   strncpy((char*)inamt,namp,12);
   ipos = abs(inamt[0] + inamt[1] * 3 + inamt[2] * 5) % lhash;
   for (i = (*iplist)->hash[ipos] - 1; i >= 0; --i) {
      if (ipnode[i].jcmt[0] == inamt[0]) {
         if (ipnode[i].jcmt[1] == inamt[1]) {
            if (ipnode[i].jcmt[2] == inamt[2]) {
               if ((ipnode[i].jjtyp != 0) && (ipnode[i].jjtyp != 10)) {
                  sprintf(AbortString,"%s: BLOCK '%s' IN OBJECT '%s' IS NOT A TABLE/LIST.",
                          nomsub,namp,(*iplist)->hname);
                  xabort_c(AbortString);
               }
               jplist = (lcm *)(ipnode[i]).jdata;
               return jplist;
            }
         }
      }
   }
L10:
   if ((*iplist)->father == NULL) {
      strcpy(nammy,"/");
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode == NULL) {
         memcpy(nammy,"            ",12);
      } else {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
      }
      nammy[12]='\0';
   }
   sprintf(AbortString,"%s: UNABLE TO FIND BLOCK '%s' INTO DIRECTORY '%s' IN THE OBJECT '%.50s'.",
           nomsub,namp,nammy,(*iplist)->hname);
   xabort_c(AbortString);
   return NULL;
}

void lcmdel_c(lcm **iplist,const char *namp)
/*
 *-----------------------------------------------------------------------
 *
 * delete an entry in the table.
 *
 * input parameters:
 *  iplist : address of the table.
 *    namp : character*12 name of the block to delete.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmdel_c";
   char nammy[13];
   blockb *ipnode;
   int_32 i, ipos;
   int inamt[3];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->listlen >= 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS A LIST.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      sprintf(AbortString,"%s: UNABLE TO DELETE RECORD '%s' FROM AN XSM FILE.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   }
   ipnode = (*iplist)->inext;
   if (ipnode == NULL) goto L10;
   strncpy((char*)inamt,namp,12);
   ipos = abs(inamt[0] + inamt[1] * 3 + inamt[2] * 5) % lhash;
   for (i = (*iplist)->hash[ipos] - 1; i >= 0; --i) {
      if (ipnode[i].jcmt[0] == inamt[0]) {
         if (ipnode[i].jcmt[1] == inamt[1]) {
            if (ipnode[i].jcmt[2] == inamt[2]) {
               if (ipnode[i].jjtyp == 0) {
                  /* DELETE AN ASSOCIATIVE TABLE. */
                  lcm *kplist;
                  kplist = (lcm *)ipnode[i].jdata;
                  lcmcl_part1(kplist);
                  lcmkep((*iplist)->icang,c__2,&kplist);
                  free(ipnode[i].jdata);
               } else if (ipnode[i].jjtyp == 10) {
                  /* DELETE A LIST. */
                  lcm *kplist;
                  kplist = (lcm *)ipnode[i].jdata;
                  lcmcl_part2(kplist->listlen,kplist);
                  lcmkep((*iplist)->icang,c__2,&kplist);
                  free(ipnode[i].jdata);
               } else if (ipnode[i].jjtyp == 99) {
                  sprintf(AbortString,"%s: BLOCK '%s' IN THE OBJECT '%.60s' IS ARLEADY DELETED.",
                          nomsub,nammy,(*iplist)->hname);
                  xabort_c(AbortString);
               } else {
                  /* DELETE A NODE. */
                  if(refpop(iplist,ipnode[i].jdata)) free(ipnode[i].jdata); /* rlsara_c */
               }
               ipnode[i].jdata = NULL;
               ipnode[i].jjlon = 0;
               ipnode[i].jjtyp = 99;
               memcpy((char*)ipnode[i].jcmt,"            ",12);
               if (i+1 == (*iplist)->inref) --(*iplist)->inref;
               return;
            }
         }
      }
   }
L10:
   if ((*iplist)->father == NULL) {
      strcpy(nammy,"/");
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode == NULL) {
         memcpy(nammy,"            ",12);
      } else {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
      }
      nammy[12]='\0';
   }
   sprintf(AbortString,"%s: UNABLE TO FIND BLOCK '%s' INTO DIRECTORY '%s' IN THE OBJECT '%.50s'.",
           nomsub,namp,nammy,(*iplist)->hname);
   xabort_c(AbortString);
}

lcm * lcmdil_c(lcm **iplist,int_32 iset)
/*
 *-----------------------------------------------------------------------
 *
 * create/access the hierarchical structure of a node table located in
 * a list.
 *
 * input parameters:
 *  iplist : address of the list.
 *    iset : position in the father list of the son table.
 *
 * output parameter:
 *  lcmdil_c : address of the son table.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmdil_c";
   lcm *jplist;
   blockb *ipnode;
   int_32 mode, lenold, ityold;
   char nammy[13];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iset < 0 || iset >= (*iplist)->listlen) {
      sprintf(AbortString,"%s: LIST INDEX %d OUT OF BOUNDS (%d,%d) IN OBJECT '%.60s'.",
              nomsub,(int)iset,0,(int)((*iplist)->listlen-1),(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *ipxsm;
      ipxsm = (xsm *)*iplist + iset;
      xsmdid_c(&ipxsm," ",(xsm **)(&jplist));
      return jplist;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->imode == 2) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS OPEN IN READ-ONLY MODE.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->father == NULL) {
      sprintf(AbortString,"%s: THE FATHER OBJECT ('/') IS NOT A LIST.",
              nomsub);
      xabort_c(AbortString);
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode[(*iplist)->ifdir - 1].jjtyp != 10) {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
         nammy[12]='\0';
         sprintf(AbortString,"%s: THE FATHER OBJECT ('%s') IS NOT A LIST.",
                 nomsub,nammy);
         xabort_c(AbortString);
      }
   }
   ipnode = ((*iplist)[iset]).inext;
   mode = ((*iplist)[iset]).imode;
   if (ipnode == NULL) {
      lenold = 0;
      ityold = 0;
   } else {
      lenold = (ipnode[0]).jjlon;
      ityold = (ipnode[0]).jjtyp;
   }

   if (ityold != 0) {
      sprintf(AbortString,"%s: LIST ELEMENT %d IS AN INVALID DIRECTORY TYPE. OBJECT='%.60s'.",
              nomsub,(int)iset,(*iplist)[iset].hname);
      xabort_c(AbortString);
   } else if (lenold != 0 && lenold != -1) {
      sprintf(AbortString,"%s: LIST ELEMENT %d OF THE OBJECT '%.60s' HAS AN INVALID LENGTH ( %d ).",
              nomsub,(int)iset,(*iplist)[iset].hname,(int)lenold);
      xabort_c(AbortString);
   } else if (lenold == 0) {
      /* CREATE A NEW NODE TABLE. */
      ipnode = (blockb *) malloc(sizeof(*ipnode));
      (*iplist)[iset].inext = ipnode;
      (*iplist)[iset].imax = 1;
      (*iplist)[iset].inref = 1;
      memcpy((char*)ipnode[0].jcmt,"            ",12);
      (ipnode[0]).jjlon = -1;
      (ipnode[0]).jjtyp = 0;
      jplist = (lcm *) malloc(sizeof(*jplist));
      jplist->header = (*iplist)->header;
      strcpy(jplist->hname, (*iplist)->hname);
      jplist->listlen = -1;
      jplist->inext = NULL;
      jplist->father = *iplist;
      jplist->ifdir = 1;
      jplist->imax = 0;
      jplist->inref = 0;
      jplist->icang = (*iplist)->icang;
      jplist->global = (*iplist)->global;
      jplist->hash = NULL;
      lcmkep(jplist->icang, c__1, &jplist);
      (ipnode[0]).jdata = (int_32 *)jplist;
   }
   /* SWITCH TO THE SON DIRECTORY. */
   jplist = (lcm *)(ipnode[0]).jdata;
   jplist->imode = mode;
   return jplist;
}

lcm * lcmlil_c(lcm **iplist,int_32 iset,int_32 ilong)
/*
 *-----------------------------------------------------------------------
 *
 * create/access the hierarchical structure of a list embedded in
 * another list.
 *
 * input parameters:
 *  iplist : address of the father list.
 *    iset : position of the embedded list in the father list.
 *   ilong : dimension of the embedded list.
 *
 * output parameter:
 *  lcmlil_c : address of the embedded list.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmlil_c";
   lcm *jplist, *iofset, *iofold;
   blockb *ipnode;
   int_32 i, mode, lenold, ityold;
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iset < 0 || iset >= (*iplist)->listlen) {
      sprintf(AbortString,"%s: LIST INDEX %d OUT OF BOUNDS (%d,%d) IN OBJECT '%.60s'.",
              nomsub,(int)iset,0,(int)((*iplist)->listlen-1),(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *ipxsm;
      ipxsm = (xsm *)*iplist + iset;
      xsmlid_c(&ipxsm," ",ilong,(xsm **)(&jplist));
      return jplist;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->imode == 2) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS OPEN IN READ-ONLY MODE.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (ilong <= 0) {
      sprintf(AbortString,"%s: INVALID LENGTH (%d) FOR LIST ELEMENT %d IN THE OBJECT '%.45s'.",
              nomsub,(int)ilong,(int)iset,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->father == NULL) {
      sprintf(AbortString,"%s: THE FATHER OBJECT ('/') IS NOT A LIST.",
              nomsub);
      xabort_c(AbortString);
   }
   ipnode = ((*iplist)[iset]).inext;
   mode = ((*iplist)[iset]).imode;
   if (ipnode == NULL) {
      lenold = 0;
      ityold = 10;
   } else {
      lenold = (ipnode[0]).jjlon;
      ityold = (ipnode[0]).jjtyp;
   }

   if (ityold != 10) {
      sprintf(AbortString,"%s: LIST ELEMENT %d IS AN INVALID LIST TYPE. TYPE=%d OBJECT='%.60s'.",
              nomsub,(int)iset,(int)ityold,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (lenold != 0 && lenold > ilong) {
      ilong = lenold;
   }

   if (lenold == 0) {
      /* CREATE A NEW NODE TABLE. */
      ipnode = (blockb *) malloc(sizeof(*ipnode));
      (*iplist)[iset].inext = ipnode;
      (*iplist)[iset].imax = 1;
      (*iplist)[iset].inref = 1;
      memcpy((char*)ipnode[0].jcmt,"            ",12);
      (ipnode[0]).jjtyp = 10;
   }

   if (ilong != lenold) {
      (ipnode[0]).jjlon = ilong;
      iofset = (lcm *) malloc(ilong*sizeof(*iofset));
      for (i = 0; i < ilong; ++i) {
         if (i < lenold) {
            iofold = (lcm *)(ipnode[0]).jdata;
            iofset[i].header = iofold[i].header;
            strcpy(iofset[i].hname, iofold[i].hname);
            iofset[i].listlen = 0;
            iofset[i].inext = iofold[i].inext;
            iofset[i].father = iofold[i].father;
            iofset[i].ifdir = iofold[i].ifdir;
            iofset[i].imode = iofold[i].imode;
            iofset[i].imax = iofold[i].imax;
            iofset[i].inref = iofold[i].inref;
            iofset[i].icang = iofold[i].icang;
            iofset[i].global = iofold[i].global;
            iofset[i].hash = iofold[i].hash;
            /* PUT THE OLD TABLE IN READ-ONLY MODE */
            iofold[i].imode = 2;
         } else {
            iofset[i].header = (*iplist)->header;
            strcpy(iofset[i].hname, (*iplist)->hname);
            iofset[i].listlen = 0;
            iofset[i].inext = NULL;
            iofset[i].father = *iplist + iset;
            iofset[i].ifdir = 1;
            iofset[i].imode = 0;
            iofset[i].imax = 0;
            iofset[i].inref = 0;
            iofset[i].icang = (*iplist)->icang;
            iofset[i].global = (*iplist)->global;
            iofset[i].hash = NULL;
         }
      }
      iofset[0].listlen = ilong;
      lcmkep(iofset->icang, c__1, &iofset);
      (ipnode[0]).jdata = (int_32 *)iofset;
   }
   /* SWITCH TO THE SON LIST. */
   jplist = (lcm *)(ipnode[0]).jdata;
   for (i = 0; i < ilong; ++i) jplist[i].imode = mode;
   return jplist;
}

void lcmppl_c(lcm **iplist,int_32 iset,int_32 ilong,int_32 itype,int_32 *iofdum)
/*
 *-----------------------------------------------------------------------
 *
 * add a new malloc pointer entry in the list.
 *
 * input parameters:
 *  iplist : address of the list.
 *    iset : position of the specific element.
 *   ilong : number of information elements stored in the current block.
 *   itype : type of information elements stored in the current block.
 *  iofdum : malloc pointer of the first information element.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmppl_c";
   blockb *ipnode;
   int_32 *jofdum;
   char nammy[13];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iset < 0 || iset >= (*iplist)->listlen) {
      sprintf(AbortString,"%s: LIST INDEX %d OUT OF BOUNDS (%d,%d) IN OBJECT '%.60s'.",
              nomsub,(int)iset,0,(int)((*iplist)->listlen-1),(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iofdum == NULL) {
      sprintf(AbortString,"%s: THE MALLOC POINTER OF LIST ELEMENT %d IS NOT SET IN THE OBJECT '%.45s'.",
              nomsub,(int)iset,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (ilong <= 0) {
      sprintf(AbortString,"%s: INVALID LENGTH ( %d ) FOR LIST ELEMENT %d IN THE OBJECT '%.45s'.",
              nomsub,(int)ilong,(int)iset,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *ipxsm;
      ipxsm = (xsm *)*iplist + iset;
      xsmppd_c(&ipxsm," ",ilong,itype,iofdum);
      return;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->imode == 2) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS OPEN IN READ-ONLY MODE.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->father == NULL) {
      sprintf(AbortString,"%s: THE FATHER OBJECT ('/') IS NOT A LIST.",
              nomsub);
      xabort_c(AbortString);
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode[(*iplist)->ifdir - 1].jjtyp != 10) {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
         nammy[12]='\0';
         sprintf(AbortString,"%s: THE FATHER OBJECT ('%s') IS NOT A LIST.",
                 nomsub,nammy);
         xabort_c(AbortString);
      }
   }
   ipnode = (*iplist)[iset].inext;
   if (ipnode == NULL) {
      ipnode = (blockb *) malloc(sizeof(*ipnode));
      (*iplist)[iset].inext = ipnode;
      (*iplist)[iset].imax = 1;
      (*iplist)[iset].inref = 1;
      memcpy((char*)ipnode[0].jcmt,"            ",12);
   } else {
      jofdum = (ipnode[0]).jdata;
      if (jofdum != iofdum) {
         if(refpop(iplist,jofdum)) free(jofdum); /* rlsara_c(jofdum); */
      }
   }

   /* STORE THE INFORMATION RELATIVE TO THE NEW INFORMATION ELEMENT. */
   (ipnode[0]).jdata = iofdum;
   (ipnode[0]).jjlon = ilong;
   (ipnode[0]).jjtyp = itype;

   /* STORE THE FIRST AND LAST ELEMENTS FOR VALIDATION PURPOSE. */
   if (itype == 1 || itype == 2 || itype == 3 || itype == 5) {
      (ipnode[0]).jidat[0] = iofdum[0];
      (ipnode[0]).jidat[1] = iofdum[ilong-1];
   } else if (itype == 4 || itype == 6) {
      (ipnode[0]).jjlon = 2*ilong;
      (ipnode[0]).jidat[0] = iofdum[0];
      (ipnode[0]).jidat[1] = iofdum[1];
      (ipnode[0]).jidat[2] = iofdum[2*ilong-2];
      (ipnode[0]).jidat[3] = iofdum[2*ilong-1];
   }
   return;
}

void lcmlel_c(lcm **iplist, int_32 iset, int_32 *ilong, int_32 *itylcm)
/*
 *-----------------------------------------------------------------------
 *
 * return the length and type of a list entry.
 *
 * input parameters:
 *  iplist : address of the list.
 *    iset : position of the specific element.
 *
 * output parameters:
 *  ilong  : number of information elements pointed by the lcm entry.
 *           ilong=0 is returned if the entry does not exists.
 *  itylcm : type of information elements pointed by the lcm entry.
 *           0: directory                1: integer
 *           2: single precision         3: character*4
 *           4: double precision         5: logical
 *           6: complex                 99: empty node
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmlel_c";
   blockb *ipnode;
   char nammy[13];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iset < 0 || iset >= (*iplist)->listlen) {
      sprintf(AbortString,"%s: LIST INDEX %d OUT OF BOUNDS (%d,%d) IN OBJECT '%.60s'.",
              nomsub,(int)iset,0,(int)((*iplist)->listlen-1),(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *ipxsm;
      ipxsm = (xsm *)*iplist + iset;
      xsmlen_c(&ipxsm," ",ilong,itylcm);
      return;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->father == NULL) {
      sprintf(AbortString,"%s: THE FATHER OBJECT ('/') IS NOT A LIST.",
              nomsub);
      xabort_c(AbortString);
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode[(*iplist)->ifdir - 1].jjtyp != 10) {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
         nammy[12]='\0';
         sprintf(AbortString,"%s: THE FATHER OBJECT ('%s') IS NOT A LIST.",
                 nomsub,nammy);
         xabort_c(AbortString);
      }
   }
   ipnode = (*iplist)[iset].inext;
   if (ipnode == NULL) {
      *ilong = 0;
      *itylcm = 99;
   } else {
      *ilong = ipnode[0].jjlon;
      *itylcm = ipnode[0].jjtyp;
      if (*itylcm == 4 || *itylcm == 6) *ilong=*ilong/2;
   }
   return;
}

void lcmgpl_c(lcm **iplist, int_32 iset, int_32 **iofdum)
/*
 *-----------------------------------------------------------------------
 *
 * get a malloc pointer for a list entry.
 *
 * input parameters:
 *  iplist : address of the list.
 *    iset : position of the specific element.
 *
 * output parameter:
 *  iofdum : malloc pointer to the iset-th list entry.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmgpl_c";
   blockb *ipnode;
   char nammy[13];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iset < 0 || iset >= (*iplist)->listlen) {
      sprintf(AbortString,"%s: LIST INDEX %d OUT OF BOUNDS (%d,%d) IN OBJECT '%.60s'.",
              nomsub,(int)iset,0,(int)((*iplist)->listlen-1),(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *ipxsm;
      ipxsm = (xsm *)*iplist + iset;
      xsmgpd_c(&ipxsm," ",iofdum);
      return;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->father == NULL) {
      sprintf(AbortString,"%s: THE FATHER OBJECT ('/') IS NOT A LIST.",
              nomsub);
      xabort_c(AbortString);
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode[(*iplist)->ifdir - 1].jjtyp != 10) {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
         nammy[12]='\0';
         sprintf(AbortString,"%s: THE FATHER OBJECT ('%s') IS NOT A LIST.",
                 nomsub,nammy);
         xabort_c(AbortString);
      }
   }
   ipnode = (*iplist)[iset].inext;
   if (ipnode != NULL) {
      *iofdum = ipnode[0].jdata;
      return;
   }
   ipnode = (*iplist)->father->inext;
   strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
   nammy[12]='\0';
   sprintf(AbortString,"%s: UNABLE TO FIND LIST ELEMENT %d INTO DIRECTORY '%s' IN OBJECT '%.45s'.",
           nomsub,(int)iset,nammy,(*iplist)->hname);
   xabort_c(AbortString);
}

void lcmgdl_c(lcm **iplist, int_32 iset, int_32 *idata)
/*
 *-----------------------------------------------------------------------
 *
 * copy a block of data from a list into memory.
 *
 * input parameters:
 *  iplist : address of the list.
 *    iset : position of the specific element.
 *
 * output parameter:
 *   idata : information elements. dimension idata1(ilong) where ilong
 *           is the number of information elements.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmgdl_c";
   blockb *ipnode;
   int_32 j;
   char nammy[13];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iset < 0 || iset >= (*iplist)->listlen) {
      sprintf(AbortString,"%s: LIST INDEX %d OUT OF BOUNDS (%d,%d) IN OBJECT '%.60s'.",
              nomsub,(int)iset,0,(int)((*iplist)->listlen-1),(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *ipxsm;
      ipxsm = (xsm *)*iplist + iset;
      xsmget_c(&ipxsm," ",idata);
      return;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->father == NULL) {
      sprintf(AbortString,"%s: THE FATHER OBJECT ('/') IS NOT A LIST.",
              nomsub);
      xabort_c(AbortString);
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode[(*iplist)->ifdir - 1].jjtyp != 10) {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
         nammy[12]='\0';
         sprintf(AbortString,"%s: THE FATHER OBJECT ('%s') IS NOT A LIST.",
                 nomsub,nammy);
         xabort_c(AbortString);
      }
   }
   ipnode = (*iplist)[iset].inext;
   if (ipnode != NULL) {
      for (j = 0; j < ipnode[0].jjlon; ++j) idata[j] = ipnode[0].jdata[j];
      return;
   }
   ipnode = (*iplist)->father->inext;
   strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
   nammy[12]='\0';
   sprintf(AbortString,"%s: UNABLE TO FIND LIST ELEMENT %d INTO DIRECTORY '%s' IN OBJECT '%.45s'.",
           nomsub,(int)iset,nammy,(*iplist)->hname);
   xabort_c(AbortString);
}

lcm * lcmgil_c(lcm **iplist, int_32 iset)
/*
 *-----------------------------------------------------------------------
 *
 * get the address of a table or of a list located in a father list.
 *
 * input parameters:
 *  iplist : address of the father list.
 *    iset : position of the specific element.
 *
 * output parameter:
 *  lcmgil_c : address of the table or of the list named namp.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmgil_c";
   blockb *ipnode;
   char nammy[13];
   if ((*iplist)->header != 100 && (*iplist)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if (iset < 0 || iset >= (*iplist)->listlen) {
      sprintf(AbortString,"%s: LIST INDEX %d OUT OF BOUNDS (%d,%d) IN OBJECT '%.60s'.",
              nomsub,(int)iset,0,(int)((*iplist)->listlen-1),(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *ipxsm, *jpxsm;
      int_32 iilong, itylcm;
      ipxsm = (xsm *)*iplist + iset;
      xsmlen_c(&ipxsm," ",&iilong,&itylcm);
      if (itylcm == 0) {
         xsmdid_c(&ipxsm," ",&jpxsm);
      } else {
         xsmlid_c(&ipxsm," ",iilong,&jpxsm);
      }
      return (lcm *)jpxsm;
   } else if ((*iplist)->imode == 0) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' IS CLOSED.",
              nomsub,(*iplist)->hname);
      xabort_c(AbortString);
   } else if ((*iplist)->father == NULL) {
      sprintf(AbortString,"%s: THE FATHER OBJECT ('/') IS NOT A LIST.",
              nomsub);
      xabort_c(AbortString);
   } else {
      lcm *my_father;
      my_father = (*iplist)->father;
      ipnode = my_father->inext;
      if (ipnode[(*iplist)->ifdir - 1].jjtyp != 10) {
         strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
         nammy[12]='\0';
         sprintf(AbortString,"%s: THE FATHER OBJECT ('%s') IS NOT A LIST.",
                 nomsub,nammy);
         xabort_c(AbortString);
      }
   }
   ipnode = (*iplist)[iset].inext;
   if (ipnode != NULL) {
      lcm *jplist;
      if ((ipnode[0].jjtyp != 0) && (ipnode[0].jjtyp != 10)) {
         sprintf(AbortString,"%s: LIST ELEMENT %d IN LIST '%s' IS NOT A TABLE/LIST.",
                 nomsub,(int)iset,(*iplist)->hname);
         xabort_c(AbortString);
      }
      jplist = (lcm *)(ipnode[0]).jdata;
      return jplist;
   }
   ipnode = (*iplist)->father->inext;
   strncpy(nammy,(char*)ipnode[(*iplist)->ifdir - 1].jcmt,12);
   nammy[12]='\0';
   sprintf(AbortString,"%s: UNABLE TO FIND LIST ELEMENT %d INTO DIRECTORY '%s' IN OBJECT '%.45s'.",
           nomsub,(int)iset,nammy,(*iplist)->hname);
   xabort_c(AbortString);
   return NULL;
}

void lcmequ_part2(int_32 ilong, lcm *iplis1, lcm *iplis2);

void lcmequ_part1(lcm *iplis1, lcm *iplis2)
/* FAST COPY OF AN ASSOCIATIVE TABLE */
{
   int_32 i, iref, ilong, itylcm;
   int inamt[3];
   blockb *inode;
   lcm *kdata2;
   char namt[13];
   inode = iplis1->inext;
   for (iref = 0; iref < iplis1->inref; ++iref) {
      memcpy((char*)inamt,"            ",12);
      if ( (inode[iref].jcmt[0] != inamt[0]) ||
           (inode[iref].jcmt[1] != inamt[1]) ||
           (inode[iref].jcmt[2] != inamt[2]) ) {
         ilong = inode[iref].jjlon;
         itylcm = inode[iref].jjtyp;
         strncpy(namt,(char*)inode[iref].jcmt,12);
         namt[12]='\0';
         if (itylcm == 0 && ilong == -1) {
            /* ASSOCIATIVE TABLE. */
            kdata2 = lcmdid_c(&iplis2, namt);
            lcmequ_part1((lcm *)inode[iref].jdata, kdata2);
         } else if (itylcm == 10) {
            /* LIST. */
            kdata2 = lcmlid_c(&iplis2, namt, ilong);
            lcmequ_part2(ilong, (lcm *)inode[iref].jdata, kdata2);
         } else {
            int_32 *iass;
            int_32 jlong = ilong;
            if(itylcm == 4 || itylcm == 6) jlong = ilong/2;
            iass = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
            for (i = 0; i < ilong; ++i) iass[i] = inode[iref].jdata[i];
            lcmppd_c(&iplis2, namt, jlong, itylcm, iass);
         }
      }
   }
   return;
}
void lcmequ_part2(int_32 ilong, lcm *iplis1, lcm *iplis2)
/* FAST COPY OF A LIST */
{
   int_32 i, ivec, kjlon, itylcm;
   blockb *knode;
   lcm *kdata2;
   for (ivec = 0; ivec < ilong; ++ivec) {
      knode = iplis1[ivec].inext;
      if (knode) {
         kjlon = knode[0].jjlon;
         itylcm = knode[0].jjtyp;
         if (itylcm == 0 && kjlon == -1) {
            /* ASSOCIATIVE TABLE. */
            kdata2 = lcmdil_c(&iplis2, ivec);
            lcmequ_part1((lcm *)knode[0].jdata, kdata2);
         } else if (itylcm == 10) {
            /* LIST. */
            kdata2=lcmlil_c(&iplis2, ivec, kjlon);
            lcmequ_part2(kjlon, (lcm *)knode[0].jdata, kdata2);
         } else {
            int_32 *iass;
            int_32 jlong = kjlon;
            if(itylcm == 4 || itylcm == 6) jlong = kjlon/2;
            iass = (int_32 *)malloc(kjlon*sizeof(int_32)); /* setara_c(kjlon); */
            for (i = 0; i < kjlon; ++i) iass[i] = knode[0].jdata[i];
            lcmppl_c(&iplis2, ivec, jlong, itylcm, iass);
         }
      }
   }
   return;
}

void lcmequ_part4(int_32 ilong, lcm *iplis1, lcm *iplis2);

void lcmequ_part3(lcm *iplis1, lcm *iplis2)
/* GENERAL COPY OF AN ASSOCIATIVE TABLE */
{
   char namlcm[73], myname[13], namt[13], first[13];
   int_32 empty, ilong, lcml, access, itylcm;
   lcm *kdata1, *kdata2;
   lcminf_c(&iplis1, namlcm, myname, &empty, &ilong, &lcml, &access);
   if (empty) return;
   strcpy(namt," ");
   lcmnxt_c(&iplis1,namt);
   strcpy(first,namt);
L10:
   lcmlen_c(&iplis1, namt, &ilong, &itylcm);
   if (ilong != 0 && itylcm == 0) {
      /* ASSOCIATIVE TABLE. */
      kdata1 = lcmgid_c(&iplis1, namt);
      kdata2 = lcmdid_c(&iplis2, namt);
      lcmequ_part3(kdata1, kdata2);
   } else if (ilong != 0 && itylcm == 10) {
      /* LIST. */
      kdata1 = lcmgid_c(&iplis1, namt);
      kdata2 = lcmlid_c(&iplis2, namt, ilong);
      lcmequ_part4(ilong, kdata1, kdata2);
   } else if (ilong != 0 && itylcm <= 6) {
      int_32 *iass;
      int_32 jlong = ilong;
      if (itylcm == 4 || itylcm == 6) jlong = 2*ilong;
      iass = (int_32 *)malloc(jlong*sizeof(int_32)); /* setara_c(jlong); */
      lcmget_c(&iplis1, namt, iass);
      lcmppd_c(&iplis2, namt, ilong, itylcm, iass);
   }
   lcmnxt_c(&iplis1, namt);
   if (strcmp(namt,first) != 0) goto L10;
   return;
}

void lcmequ_part4(int_32 ilong, lcm *iplis1, lcm *iplis2)
/* GENERAL COPY OF A LIST */
{
   int_32 ivec, kjlon, itylcm;
   lcm *kdata1, *kdata2;
   for (ivec = 0; ivec < ilong; ++ivec) {
      lcmlel_c(&iplis1, ivec, &kjlon, &itylcm);
      if (kjlon != 0 && itylcm == 0) {
         /* ASSOCIATIVE TABLE. */
         kdata1 = lcmgil_c(&iplis1, ivec);
         kdata2 = lcmdil_c(&iplis2, ivec);
         lcmequ_part3(kdata1, kdata2);
      } else if (kjlon != 0 && itylcm == 10) {
         /* LIST. */
         kdata1=lcmgil_c(&iplis1, ivec);
         kdata2=lcmlil_c(&iplis2, ivec, kjlon);
         lcmequ_part4(kjlon, kdata1, kdata2);
      } else if (kjlon != 0 && itylcm <= 6) {
         int_32 *iass;
         int_32 jlong = kjlon;
         if (itylcm == 4 || itylcm == 6) jlong = 2*kjlon;
         iass = (int_32 *)malloc(jlong*sizeof(int_32)); /* setara_c(jlong); */
         lcmgdl_c(&iplis1, ivec, iass);
         lcmppl_c(&iplis2, ivec, kjlon, itylcm, iass);
      }
   }
   return;
}

void lcmequ_c(lcm **iplis1,lcm **iplis2)
/*
 *-----------------------------------------------------------------------
 *
 * copy the information contained in the active directory of the memory
 * or xsm file object pointed by iplis1 into the table or xsm file
 * pointed by iplis2. iplis2 is not created by lcmequ.
 *
 * input parameters:
 *  iplis1 : address of the existing object.
 *
 * output parameter:
 *  iplis2 : address of the object where the copy is performed.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmequ_c";
   if ((*iplis1)->header != 100 && (*iplis1)->header !=  200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER(1).",
              nomsub,(*iplis1)->hname);
      xabort_c(AbortString);
   } else if ((*iplis2)->header != 100 && (*iplis2)->header != 200) {
      sprintf(AbortString,"%s: THE OBJECT '%.60s' HAS THE WRONG HEADER(2).",
              nomsub,(*iplis1)->hname);
      xabort_c(AbortString);
   }
   if ((*iplis1)->header == 100 && (*iplis2)->header == 100) {
      /* USE A FAST COPY ALGORITHM. */
      if ((*iplis1)->listlen == -1) {
         lcmequ_part1(*iplis1,*iplis2);
      } else {
         lcmequ_part2((*iplis1)->listlen,*iplis1,*iplis2);
      }
   } else {
      /* USE A GENERAL COPY ALGORITHM. */
      if ((*iplis1)->listlen == -1) {
         lcmequ_part3(*iplis1,*iplis2);
      } else {
         lcmequ_part4((*iplis1)->listlen,*iplis1,*iplis2);
      }
   }
   return;
}

typedef char String8[9];
typedef char String10[11];

void Ote_blanc (char *chaine)
/*
 *----------------------------------------------------------------------
 *
 * Remove lagging blank characters from a C string.
 *
 *----------------------------------------------------------------------
 */
{
   int len, i;

   len = strlen(chaine);
   for (i = len-1; i > -1; i--) {
      if(chaine[i] == '\n') {
         chaine[i] = '\0';
         break;
      }
      else if(chaine[i] != ' '  && chaine[i] != '\0') {
         chaine[i+1] = '\0';
         break;
      }
   }
}

void lcmlib_c(lcm **iplist)
/*
 *----------------------------------------------------------------------
 *
 * list the lcm entries contained in memory or a xsm file.
 *
 * input parameters:
 *  iplist : address of the object or handle to the xsm file.
 *
 *----------------------------------------------------------------------
 */
{
   char *nomsub = "LCMLIB";
   char *nomlist = "LIST";
   char *nomtable = "TABLE";
   char namlcm[73], myname[13], namt[13], first[13], isign[13];
   int_32 empty, ilong, lcm, access, imed, itylcm, ilon, iset;
   int_32 itot, inmt;
   char* ctype[]={"DIRECTORY","INTEGER","REAL","CHARACTER","DOUBLE PRECISION",
                  "LOGICAL","COMPLEX","UNDEFINED"," "," ","LIST"};
   char* cmediu[]={"TABLE","XSM FILE"};

   lcminf_c(iplist, namlcm, myname, &empty, &ilong, &lcm, &access);
   if(lcm == 0) {
      imed=2;
   }
   else{
      imed=1;
   }
   if(ilong > 0) {
      printf("\n\n %s: CONTENT OF ACTIVE %s NAMED '%s' IN THE %8s-LOCATED LCM OBJECT '%.50s':\n",
             nomsub, nomlist, myname, cmediu[imed-1], namlcm);
      itot=0;
      printf(" LIST ITEM ---    LENGTH    TYPE\n");
      for ( iset = 0; iset < ilong; iset++) {
         lcmlel_c(iplist, iset, &ilon, &itylcm);
         if(itylcm == 0 || itylcm ==10) {
            printf(" %13d%10d    %-16s\n", (int)(iset+1), (int)ilon, ctype[itylcm]);
         }
         else if(itylcm >= 1 && itylcm <= 6) {
            printf(" %13d%10d    %-16s\n", (int)(iset+1), (int)ilon, ctype[itylcm]);
            itot=itot+ilon;
         }
         else{
            printf(" %13d%10d    %-16s\n", (int)(iset+1), (int)ilon, ctype[7]);
         }
      }
      printf(" TOTAL NUMBER OF WORDS ON LIST =%10d\n", (int)itot);
   }
   else{
      printf("\n\n %s: CONTENT OF ACTIVE %5s NAMED '%s' IN THE %8s-LOCATED LCM OBJECT '%.50s':\n",
             nomsub, nomtable, myname, cmediu[imed-1], namlcm);
      if(empty == 1) {
         printf(" %s: EMPTY TABLE.\n", nomsub);
         return;
      }
      strcpy(namt, " ");
      lcmnxt_c(iplist,namt);
      strcpy(first,namt);
      printf(" BLOCK NAME------------    LENGTH    TYPE\n");
      itot=0;
      inmt=0;
      while (strcmp(namt, first) != 0 || inmt == 0) {
         inmt++;
         lcmlen_c(iplist, namt, &ilong, &itylcm);
         if(itylcm == 0 || itylcm ==10) {
            printf(" %6d  '%-12s'%10d    %-16s\n",(int)inmt,namt,(int)ilong,ctype[itylcm]);
         }
         else if(itylcm >= 1 && itylcm <= 6) {
            if((strcmp(namt,"SIGNATURE") == 0) && itylcm == 3) {
               int_32 i, ndata[13];
               lcmget_c(iplist,namt,ndata);
               for (i=0; i<3; i++) strncpy ((isign+4*i),(char *) &ndata[i], 4);
               isign[12] = '\0';
               printf(" %6d  '%-12s'%10d    %-16s='%-12s'\n",
                      (int)inmt,namt,(int)ilong,ctype[itylcm],isign);
            }
            else{
               printf(" %6d  '%-12s'%10d    %-16s\n",(int)inmt,namt,(int)ilong,ctype[itylcm]);
            }
            itot=itot+ilong;
         }
         else{
            printf(" %6d  '%-12s'%10d    %-16s\n",(int)inmt,namt,(int)ilong,ctype[7]);
         }
         lcmnxt_c(iplist, namt);
      }
      printf("\n\n TOTAL NUMBER OF WORDS IN TABLE =%10d\n", (int)itot);
   }
   fflush(stdout);
   return;
}

/****************************************/
/* C API for lcm export/import support  */
/****************************************/

void lcmnod_c(FILE *file, int_32 imode, int_32 idir, int_32 jlong,
              int_32 itylcm, int_32 *iass)
{
   char *ccc = NULL;
   int_32 *iii;
   float_32 *rrr;
   double_64 *ddd;
   int_32 *lll;
   int_32 i, j, nb_ligne, reste;
   int_32 lendat = 4;
   String10 typelogic[8];

   if (idir == 1) {
      /*     EXPORT A NODE.*/
      if(itylcm == 1) {
         /* INTEGER DATA */
         iii = (int_32*)iass;
         if( file != NULL && imode == 1) {
            fwrite(iii, sizeof(int_32), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 8 DATA BY LINE */
            nb_ligne = (int)jlong/8;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 8; j++) fprintf(file, "%10d", (int)iii[i*8+j]);
               fprintf(file, "\n");
            }
            reste = jlong%8;
            for (j = 0; j < reste; j++) fprintf(file, "%10d", (int)iii[nb_ligne*8+j]);
            if(reste != 0) fprintf(file, "\n");
         }
      } else if(itylcm == 2) {
         /* SINGLE PRECISION DATA */
         rrr = (float_32*)iass;
         if( file != NULL && imode == 1) {
            fwrite(rrr, sizeof(float_32), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 5 DATA BY LINE */
            nb_ligne = (int)jlong/5;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 5; j++) fprintf(file, "%16.8E", rrr[i*5+j]);
               fprintf(file, "\n");
            }
            reste = jlong%5;
            for (j = 0; j < reste; j++) fprintf(file, "%16.8E", rrr[nb_ligne*5+j]);
            if(reste != 0) fprintf(file, "\n");
         }
      } else if(itylcm == 3) {
         /* CHARACTER*4 DATA */
         int i;
         ccc = (char *) malloc ((int)jlong*lendat + 1); /* +1 for \0 */
         for (i=0; i<jlong; i++) strncpy ((ccc+lendat*i),(char *) (iass + i), (int)lendat);
         ccc[(int)jlong*lendat] = '\0';
         if( file != NULL && imode == 1) {
            for (i = 0; i < jlong; i++) fwrite(&lendat, sizeof(int), 1, file);
            fwrite(ccc, sizeof(char), (int)jlong*lendat, file);
         }
         else if( file != NULL && imode == 2) {
            /* 8 DATA BY LINE */
            nb_ligne = (int)jlong/8;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 8; j++) fprintf(file, "%10d", (int)lendat);
               fprintf(file, "\n");
            }
            reste = jlong%8;
            for (j = 0; j < reste; j++) fprintf(file, "%10d", (int)lendat);
            if(reste != 0) fprintf(file, "\n");
            /* 20 DATA BY LINE */
            nb_ligne = (int)jlong/20;
            reste = jlong%20;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 20*lendat; j++)
                  fprintf(file, "%1c",ccc[i*20*lendat+j]);
               fprintf(file, "\n");
            }
            for (j = 0; j < reste*lendat; j++)
               fprintf(file, "%1c",ccc[nb_ligne*20*lendat+j]);
            if(reste != 0) fprintf(file, "\n");
         }
         if(ccc != NULL) {
            free(ccc);
            ccc = NULL;
         }
      } else if(itylcm == 4) {
         /* DOUBLE PRECISION DATA */
         ddd = (double_64*)iass;
         if( file != NULL && imode == 1) {
            fwrite(ddd, sizeof(double_64), (int)jlong, file);
         }
         else if( file != NULL && imode == 2) {
            /* 4 DATA BY LINE */
            nb_ligne = (int)jlong/4;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 4; j++) fprintf(file, "%20.12E", ddd[i*4+j]);
               fprintf(file, "\n");
            }
            reste = jlong%4;
            for (j = 0; j < reste; j++)
               fprintf(file, "%20.12E", ddd[nb_ligne*4+j]);
            if(reste != 0) fprintf(file, "\n");
         }
      } else if(itylcm == 5) {
         /* LOGICAL DATA */
         lll = (int_32*)iass;
         if( file != NULL && imode == 1) {
            fwrite(lll, sizeof(int_32), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 8 DATA BY LINE */
            nb_ligne = (int)jlong/8;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 8; j++) {
                  if (lll[i*8+j] == 0) {
                     fprintf(file, "         F");
                  } else {
                     fprintf(file, "         T");
                  }
               }
               fprintf(file, "\n");
            }
            reste = jlong%8;
            for (j = 0; j < reste; j++) {
               if (lll[nb_ligne*8+j] == 0) {
                  fprintf(file, "         F");
               } else {
                  fprintf(file, "         T");
               }
            }
            fprintf(file, "\n");
         }
      } else if(itylcm == 6) {
         /* COMPLEX  DATA */
         rrr = (float_32*)iass;
         if( file != NULL && imode == 1) {
            fwrite(rrr, sizeof(float_32), 2*(int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 5 DATA BY LINE */
            nb_ligne = (int)(2*jlong/5);
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 5; j++) fprintf(file, "%16.8E", rrr[i*5+j]);
               fprintf(file, "\n");
            }
            reste = 2*jlong%5;
            for (j = 0; j < reste; j++) fprintf(file, "%16.8E", rrr[nb_ligne*5+j]);
            if(reste != 0) fprintf(file, "\n");
         }
      }
      free(iass); /* rlsara_c(iass); */
   } else if (idir == 2) {
      /* IMPORT A NODE. */
      if (itylcm == 1) {
         /* INTEGER  DATA */
         iii = (int_32*)iass;
         if( file != NULL && imode == 1) {
            fread(iii, sizeof(int_32), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 8 DATA BY LINE */
            nb_ligne = (int)jlong/8;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 8; j++) fscanf(file, "%10d", (int *)&iii[i*8+j]);
               fscanf(file, "\n");
            }
            reste = jlong%8;
            for (j = 0; j < reste; j++) fscanf(file, "%10d", (int *)&iii[nb_ligne*8+j]);
            if(reste != 0) fscanf(file, "\n");
         }
      } else if (itylcm == 2) {
         /* SINGLE PRECISION DATA */
         rrr = (float_32*)iass;
         if( file != NULL && imode == 1) {
            fread(rrr, sizeof(float_32), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 5 DATA BY LINE */
            nb_ligne = (int)jlong/5;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 5; j++) fscanf(file, "%e", &rrr[i*5+j]);
               fscanf(file, "\n");
            }
            reste = jlong%5;
            for (j = 0; j < reste; j++) fscanf(file, "%e", &rrr[nb_ligne*5+j]);
            if(reste != 0) fscanf(file, "\n");
         }
      } else if(itylcm == 3) {
         /* CHARACTER*4 DATA */
         ccc= (char *)malloc((int)jlong*lendat + 1);
         if( file != NULL && imode == 1) {
            for (i = 0; i < jlong; i++) fread(&lendat, sizeof(int), 1, file);
            fread(ccc, sizeof(char), (int)jlong*lendat, file);
         } else if( file != NULL && imode == 2) {
            /* 8 DATA BY LINE */
            nb_ligne = (int)jlong/8;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 8; j++) fscanf(file, "%10d", (int *)&lendat);
               fgetc(file);
            }
            reste = jlong%8;
            for (j = 0; j < reste; j++) fscanf(file, "%10d", (int *)&lendat);
            if(reste != 0) fgetc(file);
            /* 20 DATA BY LINE */
            nb_ligne = (int)jlong/20;
            reste = jlong%20;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 20*lendat; j++)
                  ccc[i*20*lendat+j] = fgetc(file);
               fgetc(file);
            }
            if(reste != 0) {
               for (j = 0; j < reste*lendat; j++)
                  ccc[nb_ligne*20*lendat+j] = fgetc(file);
               fgetc(file);
            }
         }
         ccc[(int)jlong*lendat] = '\0';
         strncpy((char*)iass, ccc, (int)jlong*lendat);
         if(ccc != NULL) {
            free(ccc);
            ccc = NULL;
         }
      }
      else if(itylcm == 4) {
         /* DOUBLE PRECISION DATA */
         ddd = (double_64*)iass;
         if( file != NULL && imode == 1) {
            fread(ddd, sizeof(double_64), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 4 DATA BY LINE */
            nb_ligne = (int)jlong/4;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 4; j++) fscanf(file, "%lE", &ddd[i*4+j]);
               fscanf(file, "\n");
            }
            reste = jlong%4;
            for (j = 0; j < reste; j++) fscanf(file, "%lE", &ddd[nb_ligne*4+j]);
            if(reste != 0) fscanf(file, "\n");
         }
      } else if(itylcm == 5) {
         /* LOGICAL DATA */
         lll = (int_32*)iass;
         if( file != NULL && imode == 1) {
            fread(lll, sizeof(int), (int)jlong, file);
         }
         else if( file != NULL && imode == 2) {
            /* 8 DATA BY LINE */
            nb_ligne = (int)jlong/8;
            for (i = 0; i < nb_ligne; i++) {
               fscanf(file, "%s %s %s %s %s %s %s %s\n",
                      typelogic[0], typelogic[1], typelogic[2], typelogic[3],
                      typelogic[4], typelogic[5], typelogic[6], typelogic[7]);
               for (j = 0; j < 8; j++) {
                  if(strcmp(typelogic[j],"F") == 0) lll[i*8+j] = 0;
                  else lll[i*8+j] = 1;
               }
            }
            reste = jlong%8;
            switch(reste) {
            case 1:
               fscanf(file, "%s\n",
                      typelogic[0]);
               break;
            case 2:
               fscanf(file, "%s %s\n",
                      typelogic[0], typelogic[1]);
               break;
            case 3:
               fscanf(file, "%s %s %s\n",
                      typelogic[0], typelogic[1], typelogic[2]);
               break;
            case 4:
               fscanf(file, "%s %s %s %s\n",
                      typelogic[0], typelogic[1], typelogic[2], typelogic[3]);
               break;
            case 5:
               fscanf(file, "%s %s %s %s %s\n",
                      typelogic[0], typelogic[1], typelogic[2],
                      typelogic[3], typelogic[4]);
               break;
            case 6:
               fscanf(file, "%s %s %s %s %s %s\n",
                      typelogic[0], typelogic[1], typelogic[2],
                      typelogic[3], typelogic[4], typelogic[5]);
               break;
            case 7:
               fscanf(file, "%s %s %s %s %s %s %s\n",
                      typelogic[0], typelogic[1], typelogic[2],
                      typelogic[3], typelogic[4], typelogic[5], typelogic[6]);
               break;
            }
            for (j = 0; j < reste; j++) {
               if(strcmp(typelogic[j],"F") == 0) lll[nb_ligne*8+j] = 0;
               else lll[nb_ligne*8+j] = 1;
            }
         }
      }
      else if(itylcm == 6) {
         /* COMPLEX  DATA */
         rrr = (float_32*)iass;
         if( file != NULL && imode == 1) {
            fread(rrr, sizeof(float_32), 2*(int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 5 DATA BY LINE */
            nb_ligne = (int)(2*jlong/5);
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 5; j++) fscanf(file, "%e", &rrr[i*5+j]);
               fscanf(file, "\n");
            }
            reste = 2*jlong%5;
            for (j = 0; j < reste; j++) fscanf(file, "%e", &rrr[nb_ligne*5+j]);
            if(reste != 0) fscanf(file, "\n");
         }
      }
   }
   return;
}

void lcmexp_part2(int_32 ilong, int_32 impx, int_32 imode, lcm *iplist, int_32 idir,
                  FILE *file,int_32 imed, int_32 *ilev, int_32 *itot, int_32 *lennam,
                  char namlcm[]);

void lcmexp_part1(lcm *iplist, int_32 impx, int_32 imode, int_32 idir, FILE *file,
                  int_32 imed, int_32 *ilev, int_32 *itot, int_32 *lennam)
/* GENERAL EXPORT OF AN ASSOCIATIVE TABLE */
{
   char *nomsub = "lcmexp_part1";
   String8 cmediu[2];
   char namlcm[73], myname[13], namt[13], first[13];
   int_32 empty, ilong, licm, access, itylcm, jlong;
   int zero = 0;
   lcm *kdata1;

   strcpy(cmediu[0], "TABLE");
   strcpy(cmediu[1], "XSM FILE");
   /*     FILE EXPORT.*/
   /*     ASSOCIATIVE TABLE.*/
   lcminf_c(&iplist, namlcm, myname, &empty, &ilong, &licm, &access);
   if(empty == 1) {
      if( file != NULL && imode == 1) {
         int_32 negilev = -(*ilev);
         fwrite(&negilev, sizeof(int_32), 1, file);
         fwrite(&zero, sizeof(int_32), 1, file);
         fwrite(&zero, sizeof(int_32), 1, file);
         fwrite(&zero, sizeof(int_32), 1, file);
      } else if( file != NULL && imode == 2) {
         fprintf(file, "->%8d%8d%8d%8d%32s <-   \n",(int)(-(*ilev)),zero,zero,zero," ");
      }
      return;
   }
   strcpy(namt, " ");
   lcmnxt_c(&iplist, namt);
   *lennam = 12;
   if(strcmp(namt, " ") == 0) *lennam = 0;
   strcpy(first,namt);
L10:
   lcmlen_c(&iplist, namt, &jlong, &itylcm);
   if (jlong != 0 ) {
      if (impx > 0) {
         printf(" %5d  '%-12s'%8d%8d\n", (int)(*ilev), namt, (int)itylcm, (int)jlong);
         fflush(stdout);
      }
      if( file != NULL && imode == 1) {
         fwrite(ilev, sizeof(int_32), 1, file);
         fwrite(lennam, sizeof(int_32), 1, file);
         fwrite(&itylcm, sizeof(int_32), 1, file);
         fwrite(&jlong, sizeof(int_32), 1, file);
         if ( *lennam > 0) fwrite(namt, sizeof(char), *lennam, file);
      } else if( file != NULL && imode == 2) {
         fprintf(file, "->%8d%8d%8d%8d%32s <-   \n",(int)(*ilev),(int)(*lennam),(int)itylcm,(int)jlong," ");
         if(*lennam > 0) fprintf(file, "%-80s\n", namt);
      }
      if(itylcm == 0 ) {
         /* EXPORT ASSOCIATIVE TABLE DATA.*/
         *ilev = *ilev + 1;
         kdata1 = lcmgid_c(&iplist, namt);
         lcmexp_part1(kdata1, impx, imode, idir, file, imed, ilev, itot, lennam);
         *ilev =*ilev - 1;
      } else if(itylcm ==10) {
         /* EXPORT LIST DATA.*/
         *ilev = *ilev + 1;
         kdata1 = lcmgid_c(&iplist, namt);
         lcmexp_part2(jlong, impx, imode, kdata1, idir, file, imed,
                      ilev, itot, lennam, namlcm);
         *ilev =*ilev - 1;
      } else if(itylcm <= 6) {
         int_32 *iass;
         *itot = *itot + jlong;
         ilong = jlong;
         if(itylcm == 4 || itylcm == 6) ilong = 2*jlong;
         iass = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
         lcmget_c(&iplist, namt, iass);

         /*--------------- EXPORT A NODE ---------------*/
         lcmnod_c(file, imode, idir, jlong, itylcm, iass);
         /*---------------------------------------------*/
      } else {
         sprintf(AbortString,"%s: TRY TO EXPORT UNKNOWN TYPE RECORD %d ON THE "
                 "%8s NAMED %.45s.",nomsub, (int)itylcm, cmediu[imed-1], namlcm);
         xabort_c(AbortString);
      }
   }
   lcmnxt_c(&iplist, namt);
   if (strcmp(namt,first) != 0) goto L10;
   if( file != NULL && imode == 1) {
      int_32 negilev = -(*ilev);
      fwrite(&negilev, sizeof(int_32), 1, file);
      fwrite(&zero, sizeof(int_32), 1, file);
      fwrite(&zero, sizeof(int_32), 1, file);
      fwrite(&zero, sizeof(int_32), 1, file);
   } else if( file != NULL && imode == 2) {
      fprintf(file, "->%8d%8d%8d%8d%32s <-   \n",(int)(-(*ilev)),zero,zero,zero," ");
   }
   return;
}

void lcmexp_part2(int_32 ilong, int_32 impx, int_32 imode, lcm *iplist, int_32 idir,
                  FILE *file, int_32 imed, int_32 *ilev, int_32 *itot, int_32 *lennam,
                  char namlcm[])
/* GENERAL COPY OF A LIST */
{
   char *nomsub = "lcmexp_part2";
   String8 cmediu[2];
   int_32 ivec;
   lcm *kdata1;

   strcpy(cmediu[0], "TABLE");
   strcpy(cmediu[1], "XSM FILE");
   for (ivec = 0; ivec < ilong; ++ivec) {
      int_32 jlong, itylcm;
      lcmlel_c(&iplist, ivec, &jlong, &itylcm);
      if (impx > 0) {
         printf(" %5d  '%-12s'%8d%8d\n", (int)*ilev, " ", (int)itylcm, (int)jlong);
         fflush(stdout);
      }
      if( file != NULL && imode == 1) {
         int_32 zero = 0;
         fwrite(ilev, sizeof(int_32), 1, file);
         fwrite(&zero, sizeof(int_32), 1, file);
         fwrite(&itylcm, sizeof(int_32), 1, file);
         fwrite(&jlong, sizeof(int_32), 1, file);
      } else if( file != NULL && imode == 2) {
         fprintf(file, "->%8d%8d%8d%8d%32s <-   %08d\n",(int)(*ilev),0,(int)itylcm,(int)jlong," ",(int)(ivec+1));
      }
      if (jlong != 0 && itylcm == 0) {
         /* EXPORT ASSOCIATIVE TABLE DATA. */
         *ilev = *ilev +1;
         kdata1 = lcmgil_c(&iplist, ivec);
         lcmexp_part1(kdata1, impx, imode, idir, file, imed, ilev, itot, lennam);
         *ilev =*ilev - 1;
      } else if (jlong != 0 && itylcm == 10) {
         /* EXPORT LIST DATA. */
         *ilev = *ilev +1;
         kdata1=lcmgil_c(&iplist, ivec);
         lcmexp_part2(jlong, impx, imode, kdata1, idir, file, imed,
                      ilev, itot, lennam, namlcm);
         *ilev =*ilev - 1;
      } else if (jlong != 0 && itylcm <= 6) {
         int_32 *iass, kjlon;
         *itot = *itot + jlong;
         kjlon = jlong;
         if(itylcm == 4 || itylcm == 6) kjlon = 2*jlong;
         iass = (int_32 *)malloc(kjlon*sizeof(int_32)); /* setara_c(kjlon); */
         lcmgdl_c(&iplist, ivec, iass);
         /*--------------- EXPORT A NODE ---------------*/
         lcmnod_c(file, imode, idir, jlong, itylcm, iass);
         /*---------------------------------------------*/
      } else if (jlong != 0) {
         sprintf(AbortString, "%s: TRY TO IMPORT BAD TYPE RECORD %d ON THE %8s "
                 "NAMED %.50s.",nomsub, (int)itylcm, cmediu[imed-1], namlcm);
         xabort_c(AbortString);
      }
   }
   return;
}

void lcmexp_part4(int_32 ilong, int_32 impx, int_32 imode, lcm *iplist, int_32 idir,
                  FILE *file, int_32 imed, int_32 *ilev, int_32 *itot, int_32 *lennam,
                  char namlcm[]);

void lcmexp_part3(lcm *iplist, int_32 impx, int_32 imode, int_32 idir, FILE *file,
                  int_32 imed, int_32 *ilev, int_32 *itot, int_32 *lennam, char namlcm[])
/* GENERAL IMPORT OF AN ASSOCIATIVE TABLE */
{
   char *nomsub = "lcmexp_part3";
   String8 cmediu[2];
   char namt[13];
   int_32 ilong, itylcm, jlong;
   int jtylcm;
   lcm *kdata1;
   int_32 jlev;

   strcpy(cmediu[0], "TABLE");
   strcpy(cmediu[1], "XSM FILE");
L10:
   if( file != NULL && imode == 1) {
      if(fread(&jlev, sizeof(int_32), 1, file) == EOF) return;
      fread(lennam, sizeof(int_32), 1, file);
      fread(&itylcm, sizeof(int_32), 1, file);
      fread(&ilong, sizeof(int_32), 1, file);
   } else if( file != NULL && imode == 2) {
      if(fscanf(file, "->%d%d%d%d%*s \n", (int *)(&jlev),(int *)lennam,(int *)(&itylcm),(int *)(&ilong)) == EOF) {
         return;
      }
   }
   jtylcm = itylcm;
   if (jlev == *ilev) {
      namt[12] = '\0';
      if ( *lennam == 0) strcpy(namt, " ");
      else if( file != NULL && imode == 1) {
         if ( *lennam > 0) fread(namt, sizeof(char), *lennam, file);
      } else if( file != NULL && imode == 2) {
         if(*lennam > 0) {
            fscanf(file, "%c%c%c%c%c%c%c%c%c%c%c%c\n",
                   &namt[0], &namt[1], &namt[2], &namt[3], &namt[4],
                   &namt[5], &namt[6], &namt[7], &namt[8], &namt[9],
                   &namt[10], &namt[11]);
            Ote_blanc(namt);
         }
      }
      if(impx > 0) {
         printf("\n %5d  '%-12s'%8d%8d", (int)jlev, namt, (int)itylcm, (int)ilong);
         fflush(stdout);
      }
      if(jtylcm == 0 ) {
         /* IMPORT ASSOCIATIVE TABLE DATA.*/
         *ilev = *ilev + 1;
         kdata1 = lcmdid_c(&iplist, namt);
         lcmexp_part3(kdata1, impx, imode, idir, file, imed, ilev,
                      itot, lennam, namlcm);
         *ilev =*ilev - 1;
      } else if (jtylcm == 10) {
         /* IMPORT LIST DATA.*/
         *ilev = *ilev + 1;
         kdata1 = lcmlid_c(&iplist, namt, ilong);
         lcmexp_part4(ilong, impx, imode, kdata1, idir, file, imed,
                      ilev, itot, lennam, namlcm);
         *ilev =*ilev - 1;
      } else if (jtylcm <= 6) {
         int_32 *iass;
         jlong = ilong;
         if(jtylcm == 4 || jtylcm == 6) jlong = 2*ilong;
         iass = (int_32 *)malloc(jlong*sizeof(int_32)); /* setara_c(jlong); */
         /*--------------- IMPORT A NODE ---------------*/
         lcmnod_c(file, imode, idir, ilong, itylcm, iass);
         /*---------------------------------------------*/
         lcmppd_c(&iplist, namt, ilong, itylcm, iass);
         *itot = *itot + jlong;
      } else {
         sprintf(AbortString, "%s: TRY TO IMPORT UNKNOWN TYPE RECORD %d ON "
                 "THE %8s NAMED %.50s.",nomsub, (int)itylcm, cmediu[imed-1], namlcm);
         xabort_c(AbortString);
      }
      goto L10;
   } else if (jlev == -(*ilev)) {
      return;
   } else {
      sprintf(AbortString, "%s: UNABLE TO IMPORT '%8s' NAMED '%.50s'.",
              nomsub, cmediu[imed-1], namlcm);
      xabort_c(AbortString);
   }
}

void lcmexp_part4(int_32 jlong, int_32 impx, int_32 imode, lcm *iplist, int_32 idir,
                  FILE *file, int_32 imed, int_32 *ilev, int_32 *itot, int_32 *lennam,
                  char namlcm[])
/* GENERAL IMPORT OF A LIST */
{
   char *nomsub = "lcmexp_part4";
   String8 cmediu[2];
   int_32 ivec, ilong, itylcm;
   int jtylcm;
   lcm *kdata1;
   int_32 jlev;

   strcpy(cmediu[0], "TABLE");
   strcpy(cmediu[1], "XSM FILE");
   for (ivec = 0; ivec < jlong; ++ivec) {
      if( file != NULL && imode == 1) {
         if(fread(&jlev, sizeof(int_32), 1, file) == EOF) return;
         fread(lennam, sizeof(int_32), 1, file);
         fread(&itylcm, sizeof(int_32), 1, file);
         fread(&ilong, sizeof(int_32), 1, file);
      } else if (file != NULL && imode == 2) {
         if(fscanf(file, "->%d%d%d%d%*s%*d \n",
                   (int *)(&jlev),(int *)lennam,(int *)(&itylcm),(int *)(&ilong)) == EOF) {
            return;
         }
      }
      jtylcm = itylcm;
      if (jlev != *ilev) {
         sprintf(AbortString,"%s: INVALID LIST LEVEL ON THE '%8s' NAMED '%.50s'.",
                 nomsub, cmediu[imed-1], namlcm);
         xabort_c(AbortString);
      }
      if(impx > 0) {
         printf("\n %5d  '%-12s'%8d%8d", (int)jlev, " ", (int)itylcm,(int)ilong);
         fflush(stdout);
      }
      if (ilong != 0 && jtylcm == 0) {
         /* IMPORT ASSOCIATIVE TABLE DATA. */
         *ilev = *ilev + 1;
         kdata1 = lcmdil_c(&iplist, ivec);
         lcmexp_part3(kdata1, impx, imode, idir, file, imed, ilev, itot,
                      lennam, namlcm);
         *ilev =*ilev - 1;
      } else if (ilong != 0 && jtylcm == 10) {
         /* IMPORT LIST DATA. */
         *ilev = *ilev + 1;
         kdata1=lcmlil_c(&iplist, ivec, ilong);
         lcmexp_part4(ilong, impx, imode, kdata1, idir, file, imed,
                      ilev, itot, lennam, namlcm);
         *ilev =*ilev - 1;
      } else if (ilong != 0 && jtylcm <= 6) {
         int_32 *iass, kjlon;
         kjlon = ilong;
         if(jtylcm == 4 || jtylcm == 6) kjlon = 2*ilong;
         iass = (int_32 *)malloc(kjlon*sizeof(int_32)); /* setara_c(kjlon); */
         /*--------------- IMPORT A NODE ---------------*/
         lcmnod_c(file, imode, idir, ilong, itylcm, iass);
         /*---------------------------------------------*/
         lcmppl_c(&iplist, ivec, ilong, itylcm, iass);
         *itot = *itot + jlong;
      } else if(ilong != 0) {
         sprintf(AbortString, "%s: TRY TO IMPORT UNKNOWN TYPE RECORD %d ON "
                 "THE %8s NAMED %.50s.",nomsub, (int)itylcm, cmediu[imed-1], namlcm);
         xabort_c(AbortString);
      }
   }
   return;
}

void lcmexp_c(lcm **iplist, int_32 impx, FILE *file, int_32 imode, int_32 idir)
/*
 *-----------------------------------------------------------------------
 *
 * export/import the content of a table or xsm file using the contour
 * method. Export start from the active directory.
 *
 *  iplist : address of the table or handle to the xsm file.
 *  impx   : equal to zero for no print.
 *  nunit  : file unit number where the export/import is performed.
 *  imode  : type of export/import file:
 *           =1 sequential unformatted; =2 sequential formatted (ascii).
 *  idir   : =1 to export ; =2 to import.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub = "lcmexp_c";
   String8 cmediu[2];
   char namlcm[73], myname[13];
   int_32 empty, ilong, lcm, access, imed;
   int_32 itot, ilev, lennam;
   FILE *fileout = NULL;

   strcpy(cmediu[0], "TABLE");
   strcpy(cmediu[1], "XSM FILE");

   lcminf_c(iplist, namlcm, myname, &empty, &ilong, &lcm, &access);
   imed=2;
   if(lcm == 1) imed=1;
   if(imode < 1 || imode > 2) {
      sprintf(AbortString, "%s: INVALID FILE TYPE ON THE %8s NAMED '%.50s'.",
              nomsub, cmediu[imed-1], namlcm);
      xabort_c(AbortString);
   } else if(idir != 1 && idir != 2) {
      sprintf(AbortString, "%s: INVALID ACTION ON THE %8s NAMED '%.50s'.",
              nomsub, cmediu[imed-1], namlcm);
      xabort_c(AbortString);
   }
   if (file == NULL) {
      sprintf(AbortString, "%s: NULL IMPORT/EXPORT FILE.", nomsub);
      xabort_c(AbortString);
   } else {
      fileout = file;
   }
   itot = 0;
   ilev = 1;
   if (idir == 1) {
      /* FILE EXPORT. ALGORITHM. */
      if ((*iplist)->listlen == -1) {
         lcmval_c(iplist," ");
         lcmexp_part1(*iplist, impx, imode, idir, fileout, imed, &ilev, &itot,
                      &lennam);
      } else {
         lcmexp_part2((*iplist)->listlen, impx, imode, *iplist, idir, fileout, imed,
                      &ilev, &itot, &lennam, namlcm);
      }
      if(impx > 0) printf("\n TOTAL NUMBER OF WORDS EXPORTED =%10d\n",(int)itot);
   } else {
      /* FILE IMPORT. ALGORITHM. */
      if(impx > 0) {
         printf("\n\n%s: %6s %8s NAMED '%-12s' FROM ACTIVE DIRECTORY '%.50s' :"
                "\n\n LEVEL  BLOCK NAME---     TYPE  LENGTH\n",
                nomsub, "IMPORT",cmediu[imed-1],namlcm,myname);
         fflush(stdout);
      }
      if ((*iplist)->listlen == -1) {
         lcmexp_part3(*iplist, impx, imode, idir, fileout, imed, &ilev, &itot,
                      &lennam, namlcm);
      } else {
         lcmexp_part4((*iplist)->listlen, impx, imode, *iplist, idir, fileout, imed,
                      &ilev, &itot, &lennam, namlcm);
      }
      if(impx > 0) printf("\n TOTAL NUMBER OF WORDS IMPORTED =%10d\n",(int)itot);
   }
   fflush(stdout);
   return;
}

void lcmnodv3_c(FILE *file, int_32 imode, int_32 idir, int_32 jlong,
              int_32 itylcm, int_32 *iass)
{
   char *ccc = NULL;
   int_32 *iii;
   float_32 *rrr;
   double_64 *ddd;
   int_32 *lll;
   int_32 i, j, nb_ligne, reste;
   int_32 lendat = 4;
   String10 typelogic[8];
   
   if (idir == 1) {
      /*     EXPORT A NODE.*/
      if(itylcm == 1) {
         /* INTEGER DATA */
         iii = (int_32*)iass;
         if( file != NULL && imode == 1) {
            fwrite(iii, sizeof(int_32), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 8 DATA BY LINE */
            nb_ligne = (int)jlong/8;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 8; j++) fprintf(file, "%10d", (int)iii[i*8+j]);
               fprintf(file, "\n");
            }
            reste = jlong%8;
            for (j = 0; j < reste; j++) fprintf(file, "%10d", (int)iii[nb_ligne*8+j]);
            if(reste != 0) fprintf(file, "\n");
         }
      } else if(itylcm == 2) {
         /* SINGLE PRECISION DATA */
         rrr = (float_32*)iass;
         if( file != NULL && imode == 1) {
            fwrite(rrr, sizeof(float_32), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 5 DATA BY LINE */
            nb_ligne = (int)jlong/5;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 5; j++) fprintf(file, "%16.8E", rrr[i*5+j]);
               fprintf(file, "\n");
            }
            reste = jlong%5;
            for (j = 0; j < reste; j++) fprintf(file, "%16.8E", rrr[nb_ligne*5+j]);
            if(reste != 0) fprintf(file, "\n");
         }
      } else if(itylcm == 3) {
         /* CHARACTER*4 DATA */
         int i;
         ccc = (char *) malloc ((int)jlong*lendat + 1); /* +1 for \0 */
         for (i=0; i<jlong; i++) strncpy ((ccc+lendat*i),(char *) (iass + i), (int)lendat);
         ccc[(int)jlong*lendat] = '\0';
         if( file != NULL && imode == 1) {
            fwrite(ccc, sizeof(char), (int)jlong*lendat, file);
         }
         else if( file != NULL && imode == 2) {
            /* 20 DATA BY LINE */
            nb_ligne = (int)jlong/20;
            reste = jlong%20;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 20*lendat; j++)
                  fprintf(file, "%1c",ccc[i*20*lendat+j]);
               fprintf(file, "\n");
            }
            for (j = 0; j < reste*lendat; j++)
               fprintf(file, "%1c",ccc[nb_ligne*20*lendat+j]);
            if(reste != 0) fprintf(file, "\n");
         }
         if(ccc != NULL) {
            free(ccc);
            ccc = NULL;
         }
      } else if(itylcm == 4) {
         /* DOUBLE PRECISION DATA */
         ddd = (double_64*)iass;
         if( file != NULL && imode == 1) {
            fwrite(ddd, sizeof(double_64), (int)jlong, file);
         }
         else if( file != NULL && imode == 2) {
            /* 4 DATA BY LINE */
            nb_ligne = (int)jlong/4;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 4; j++) fprintf(file, "%20.12E", ddd[i*4+j]);
               fprintf(file, "\n");
            }
            reste = jlong%4;
            for (j = 0; j < reste; j++)
               fprintf(file, "%20.12E", ddd[nb_ligne*4+j]);
            if(reste != 0) fprintf(file, "\n");
         }
      } else if(itylcm == 5) {
         /* LOGICAL DATA */
         lll = (int_32*)iass;
         if( file != NULL && imode == 1) {
            fwrite(lll, sizeof(int_32), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 8 DATA BY LINE */
            nb_ligne = (int)jlong/8;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 8; j++) {
                  if (lll[i*8+j] == 0) {
                     fprintf(file, "         F");
                  } else {
                     fprintf(file, "         T");
                  }
               }
               fprintf(file, "\n");
            }
            reste = jlong%8;
            for (j = 0; j < reste; j++) {
               if (lll[nb_ligne*8+j] == 0) {
                  fprintf(file, "         F");
               } else {
                  fprintf(file, "         T");
               }
            }
            fprintf(file, "\n");
         }
      } else if(itylcm == 6) {
         /* COMPLEX  DATA */
         rrr = (float_32*)iass;
         if( file != NULL && imode == 1) {
            fwrite(rrr, sizeof(float_32), 2*(int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 5 DATA BY LINE */
            nb_ligne = (int)(2*jlong/5);
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 5; j++) fprintf(file, "%16.8E", rrr[i*5+j]);
               fprintf(file, "\n");
            }
            reste = 2*jlong%5;
            for (j = 0; j < reste; j++) fprintf(file, "%16.8E", rrr[nb_ligne*5+j]);
            if(reste != 0) fprintf(file, "\n");
         }
      }
      free(iass); /* rlsara_c(iass); */
   } else if (idir == 2) {
      /* IMPORT A NODE. */
      if (itylcm == 1) {
         /* INTEGER  DATA */
         iii = (int_32*)iass;
         if( file != NULL && imode == 1) {
            fread(iii, sizeof(int_32), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 8 DATA BY LINE */
            nb_ligne = (int)jlong/8;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 8; j++) fscanf(file, "%10d", (int *)&iii[i*8+j]);
               fscanf(file, "\n");
            }
            reste = jlong%8;
            for (j = 0; j < reste; j++) fscanf(file, "%10d", (int *)&iii[nb_ligne*8+j]);
            if(reste != 0) fscanf(file, "\n");
         }
      } else if (itylcm == 2) {
         /* SINGLE PRECISION DATA */
         rrr = (float_32*)iass;
         if( file != NULL && imode == 1) {
            fread(rrr, sizeof(float_32), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 5 DATA BY LINE */
            nb_ligne = (int)jlong/5;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 5; j++) fscanf(file, "%e", &rrr[i*5+j]);
               getc(file);
            }
            reste = jlong%5;
            for (j = 0; j < reste; j++) fscanf(file, "%e", &rrr[nb_ligne*5+j]);
            if(reste != 0) getc(file);
         }
      } else if(itylcm == 3) {
         /* CHARACTER*4 DATA */
         ccc= (char *)malloc((int)jlong*lendat + 1);
         if( file != NULL && imode == 1) {
            fread(ccc, sizeof(char), (int)jlong*lendat, file);
         } else if( file != NULL && imode == 2) {
            /* 20 DATA BY LINE */
            nb_ligne = (int)jlong/20;
            reste = jlong%20;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 20*lendat; j++) {
                  ccc[i*20*lendat+j] = getc(file);
               }
               getc(file);
            }
            if(reste != 0) {
               for (j = 0; j < reste*lendat; j++) {
                  ccc[nb_ligne*20*lendat+j] = getc(file);
               }
               getc(file);
            }
         }
         ccc[(int)jlong*lendat] = '\0';
         strncpy((char*)iass, ccc, (int)jlong*lendat);
         if(ccc != NULL) {
            free(ccc);
            ccc = NULL;
         }
      }
      else if(itylcm == 4) {
         /* DOUBLE PRECISION DATA */
         ddd = (double_64*)iass;
         if( file != NULL && imode == 1) {
            fread(ddd, sizeof(double_64), (int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 4 DATA BY LINE */
            nb_ligne = (int)jlong/4;
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 4; j++) fscanf(file, "%lE", &ddd[i*4+j]);
               fscanf(file, "\n");
            }
            reste = jlong%4;
            for (j = 0; j < reste; j++) fscanf(file, "%lE", &ddd[nb_ligne*4+j]);
            if(reste != 0) fscanf(file, "\n");
         }
      } else if(itylcm == 5) {
         /* LOGICAL DATA */
         lll = (int_32*)iass;
         if( file != NULL && imode == 1) {
            fread(lll, sizeof(int), (int)jlong, file);
         }
         else if( file != NULL && imode == 2) {
            /* 8 DATA BY LINE */
            nb_ligne = (int)jlong/8;
            for (i = 0; i < nb_ligne; i++) {
               fscanf(file, "%s %s %s %s %s %s %s %s\n",
                     typelogic[0], typelogic[1], typelogic[2], typelogic[3],
                     typelogic[4], typelogic[5], typelogic[6], typelogic[7]);
               for (j = 0; j < 8; j++) {
                  if(strcmp(typelogic[j],"F") == 0) lll[i*8+j] = 0;
                  else lll[i*8+j] = 1;
               }
            }
            reste = jlong%8;
            switch(reste) {
               case 1:
                  fscanf(file, "%s\n",
                        typelogic[0]);
                  break;
               case 2:
                  fscanf(file, "%s %s\n",
                        typelogic[0], typelogic[1]);
                  break;
               case 3:
                  fscanf(file, "%s %s %s\n",
                        typelogic[0], typelogic[1], typelogic[2]);
                  break;
               case 4:
                  fscanf(file, "%s %s %s %s\n",
                        typelogic[0], typelogic[1], typelogic[2], typelogic[3]);
                  break;
               case 5:
                  fscanf(file, "%s %s %s %s %s\n",
                        typelogic[0], typelogic[1], typelogic[2],
                        typelogic[3], typelogic[4]);
                  break;
               case 6:
                  fscanf(file, "%s %s %s %s %s %s\n",
                        typelogic[0], typelogic[1], typelogic[2],
                        typelogic[3], typelogic[4], typelogic[5]);
                  break;
               case 7:
                  fscanf(file, "%s %s %s %s %s %s %s\n",
                        typelogic[0], typelogic[1], typelogic[2],
                        typelogic[3], typelogic[4], typelogic[5], typelogic[6]);
                  break;
            }
            for (j = 0; j < reste; j++) {
               if(strcmp(typelogic[j],"F") == 0) lll[nb_ligne*8+j] = 0;
               else lll[nb_ligne*8+j] = 1;
            }
         }
      }
      else if(itylcm == 6) {
         /* COMPLEX  DATA */
         rrr = (float_32*)iass;
         if( file != NULL && imode == 1) {
            fread(rrr, sizeof(float_32), 2*(int)jlong, file);
         } else if( file != NULL && imode == 2) {
            /* 5 DATA BY LINE */
            nb_ligne = (int)(2*jlong/5);
            for (i = 0; i < nb_ligne; i++) {
               for (j = 0; j < 5; j++) fscanf(file, "%e", &rrr[i*5+j]);
               fscanf(file, "\n");
            }
            reste = 2*jlong%5;
            for (j = 0; j < reste; j++) fscanf(file, "%e", &rrr[nb_ligne*5+j]);
            if(reste != 0) fscanf(file, "\n");
         }
      }
   }
   return;
}

void lcmexpv3_part1(lcm *iplist, int_32 impx, int_32 imode, int_32 idir, FILE *file,
                  int_32 imed, int_32 *ilev, int_32 *itot, int_32 *lennam)
/* GENERAL EXPORT OF AN ASSOCIATIVE TABLE */
{
   char *nomsub = "lcmexpv3_part1";
   String8 cmediu[2];
   char namlcm[73], myname[13], namt[13], first[13];
   int_32 empty, ilong, licm, access, itylcm, jlong;
   lcm *kdata1;

   strcpy(cmediu[0], "TABLE");
   strcpy(cmediu[1], "XSM FILE");
   /*     FILE EXPORT.*/
   /*     ASSOCIATIVE TABLE.*/
   lcminf_c(&iplist, namlcm, myname, &empty, &ilong, &licm, &access);
   if(empty == 1) return;
   strcpy(namt, " ");
   lcmnxt_c(&iplist, namt);
   *lennam = 12;
   if(strcmp(namt, " ") == 0) *lennam = 0;
   strcpy(first,namt);
L10:
   lcmlen_c(&iplist, namt, &jlong, &itylcm);
   if (jlong != 0) {
      if(itylcm == 3) jlong = jlong*4;
      if (impx > 0) {
         printf(" %5d  '%-12s'%8d%8d\n", (int)(*ilev), namt, (int)itylcm, (int)jlong);
         fflush(stdout);
      }
      if( file != NULL && imode == 1) {
         fwrite(ilev, sizeof(int_32), 1, file);
         fwrite(namt, sizeof(char), *lennam, file);
         fwrite(&itylcm, sizeof(int_32), 1, file);
         fwrite(&jlong, sizeof(int_32), 1, file);
      } else if( file != NULL && imode == 2) {
         fprintf(file," %5d  '%-12s'%8d%8d\n", (int)(*ilev), namt, (int)itylcm, (int)jlong);
      }
      if(itylcm == 3) jlong = jlong/4;
      if(itylcm == 0) {
         /* EXPORT ASSOCIATIVE TABLE DATA.*/
         *ilev = *ilev + 1;
         kdata1 = lcmgid_c(&iplist, namt);
         lcmexpv3_part1(kdata1, impx, imode, idir, file, imed, ilev, itot, lennam);
         *ilev =*ilev - 1;
      } else if(itylcm <= 6) {
         int_32 *iass;
         *itot = *itot + jlong;
         ilong = jlong;
         if(itylcm == 4 || itylcm == 6) ilong = 2*jlong;
         iass = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
         lcmget_c(&iplist, namt, iass);

         /*---------------- EXPORT A NODE ----------------*/
         lcmnodv3_c(file, imode, idir, jlong, itylcm, iass);
         /*-----------------------------------------------*/
      } else {
         sprintf(AbortString,"%s: TRY TO EXPORT UNKNOWN TYPE RECORD %d ON THE "
                 "%8s NAMED %.45s.",nomsub, (int)itylcm, cmediu[imed-1], namlcm);
         xabort_c(AbortString);
      }
   }
   lcmnxt_c(&iplist, namt);
   if (strcmp(namt,first) != 0) goto L10;
   return;
}

void lcmexpv3_part3(lcm *iplist, int_32 impx, int_32 imode, int_32 idir, FILE *file,
                  int_32 imed, int_32 *ilev, int_32 *itot, int_32 *lennam, char namlcm[])
/* GENERAL IMPORT OF AN ASSOCIATIVE TABLE */
{
   char *nomsub = "lcmexpv3_part3";
   String8 cmediu[2];
   char namt[13];
   int_32 ilong, itylcm;
   int jtylcm;
   lcm *kdata1[100];
   int_32 jlev;
   kdata1[0]=iplist;

   strcpy(cmediu[0], "TABLE");
   strcpy(cmediu[1], "XSM FILE");
L10:
   if( file != NULL && imode == 1) {
      if(fread(&jlev, sizeof(int_32), 1, file) == EOF) return;
      fread(namt, sizeof(char), 12, file);
      fread(&itylcm, sizeof(int_32), 1, file);
      fread(&ilong, sizeof(int_32), 1, file);
   } else if( file != NULL && imode == 2) {
      if(fscanf(file, " %5d  '%c%c%c%c%c%c%c%c%c%c%c%c'%8d%8d", (int *)(&jlev),
              &namt[0], &namt[1], &namt[2], &namt[3], &namt[4],
              &namt[5], &namt[6], &namt[7], &namt[8], &namt[9],
              &namt[10], &namt[11],
                (int *)(&itylcm),(int *)(&ilong)) == EOF) {
         return;
      } else {
         getc(file);
      }
   }
   if(itylcm == 3 ) ilong = ilong/4;
   jtylcm = itylcm;
   if (jlev <= *ilev) {
      *ilev=jlev;
      namt[12] = '\0';
      Ote_blanc(namt);
      if(jtylcm == 0 ) {
         /* IMPORT ASSOCIATIVE TABLE DATA.*/
         *ilev = *ilev + 1;
         kdata1[*ilev-1] = lcmdid_c((&kdata1[*ilev-2]), namt);
      } else if (jtylcm <= 6) {
         int_32 *iass;
         int_32 jlong = ilong;
         if(jtylcm == 4 || jtylcm == 6) jlong = 2*ilong;
         iass = (int_32 *)malloc(jlong*sizeof(int_32)); /* setara_c(jlong); */
         /*---------------- IMPORT A NODE ----------------*/
         lcmnodv3_c(file, imode, idir, ilong, itylcm, iass);
         /*-----------------------------------------------*/
         lcmppd_c(&kdata1[*ilev-1], namt, ilong, itylcm, iass);
         *itot = *itot + jlong;
      } else {
         sprintf(AbortString, "%s: TRY TO IMPORT UNKNOWN TYPE RECORD %d ON "
               "THE %8s NAMED %.50s.",nomsub, (int)itylcm, cmediu[imed-1], namlcm);
         xabort_c(AbortString);
      }
      goto L10;
   } else {
      printf("\n%8d<>%8d\n", (int)jlev, (int ) *ilev);
      sprintf(AbortString, "%s: UNABLE TO IMPORT '%8s' NAMED '%.50s'.",
            nomsub, cmediu[imed-1], namlcm);
      xabort_c(AbortString);
   }
}

void lcmexpv3_c(lcm **iplist, int_32 impx, FILE *file, int_32 imode, int_32 idir)
/*
 *-----------------------------------------------------------------------
 *
 * export/import the content of a table or xsm file using the contour
 * method for version 3. Export start from the active directory.  
 *
 *  iplist : address of the table or handle to the xsm file.
 *  impx   : equal to zero for no print.
 *  nunit  : file unit number where the export/import is performed.
 *  imode  : type of export/import file:
 *           =1 sequential unformatted; =2 sequential formatted (ascii).
 *  idir   : =1 to export ; =2 to import.
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub = "lcmexpv3_c";
   String8 cmediu[2];
   char namlcm[73], myname[13];
   int_32 empty, ilong, lcm, access, imed;
   int_32 itot, ilev, lennam;
   FILE *fileout = NULL;
   
   strcpy(cmediu[0], "TABLE");
   strcpy(cmediu[1], "XSM FILE");
   
   lcminf_c(iplist, namlcm, myname, &empty, &ilong, &lcm, &access);
   imed=2;
   if(lcm == 1) imed=1;
   if(imode < 1 || imode > 2) {
      sprintf(AbortString, "%s: INVALID FILE TYPE ON THE %8s NAMED '%.50s'.",
            nomsub, cmediu[imed-1], namlcm);
      xabort_c(AbortString);
   } else if(idir != 1 && idir != 2) {
      sprintf(AbortString, "%s: INVALID ACTION ON THE %8s NAMED '%.50s'.",
            nomsub, cmediu[imed-1], namlcm);
      xabort_c(AbortString);
   }
   if (file == NULL) {
      sprintf(AbortString, "%s: NULL IMPORT/EXPORT FILE.", nomsub);
      xabort_c(AbortString);
   } else {
      fileout = file;
   }
   itot = 0;
   ilev = 1;
   if (idir == 1) {
      /* FILE EXPORT. ALGORITHM. */
      lcmval_c(iplist," ");
      lcmexpv3_part1(*iplist, impx, imode, idir, fileout, imed, &ilev, &itot, &lennam);
      if(impx > 0) printf("\n TOTAL NUMBER OF WORDS EXPORTED =%10d\n",(int)itot);
   } else {
      /* FILE IMPORT. ALGORITHM. */
      if(impx > 0) {
         printf("\n\n%s: %6s %8s NAMED '%-12s' FROM ACTIVE DIRECTORY '%-12s' :"
               "\n\n LEVEL  BLOCK NAME---     TYPE  LENGTH\n",
               nomsub, "IMPORT",cmediu[imed-1],namlcm,myname);
         fflush(stdout);
      }
         lcmexpv3_part3(*iplist, impx, imode, idir, fileout, imed, &ilev, &itot,
                   &lennam, namlcm);
      if(impx > 0) printf("\n TOTAL NUMBER OF WORDS IMPORTED =%10d\n",(int)itot);
   }
   fflush(stdout);
   return;
}

long lcmcast_c(lcm **iplist)
/* cast a LCM pointer into an integer (not 64-bit clean) */
{
   long ret_val = (long) *iplist;
   return ret_val;
}
