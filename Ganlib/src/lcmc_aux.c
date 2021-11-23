
/**********************************/
/* C API for lcm object support   */
/* (auxiliary functions)          */
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

static char AbortString[132];

void lcmput_c(lcm **iplist,const char *nom,int_32 ilong,int_32 itype,int_32 *idata)
/*
 *----------------------------------------------------------------------
 *
 * COPY A BLOCK OF DATA FROM MEMORY INTO A TABLE.
 *
 * INPUT PARAMETERS:
 *  IPLIST : ADDRESS OF THE TABLE.
 *    NAMP : CHARACTER*12 NAME OF THE CURRENT BLOCK.
 *   ILONG : NUMBER OF INFORMATION ELEMENTS STORED IN THE CURRENT BLOCK.
 *   ITYPE : TYPE OF INFORMATION ELEMENTS STORED IN THE CURRENT BLOCK.
 *           0: DIRECTORY                1: INTEGER
 *           2: SINGLE PRECISION         3: CHARACTER*4
 *           4: DOUBLE PRECISION         5: LOGICAL
 *           6: COMPLEX                 99: UNDEFINED
 *   IDATA : INFORMATION ELEMENTS. DIMENSION IDATA(ILONG)
 *
 *----------------------------------------------------------------------
 */
{
   if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsmput_c((xsm **)iplist,nom,ilong,itype,idata);
   } else {
      int_32 i, *iofdat;
      int_32 jlong = ilong;
      if (itype == 4 || itype == 6) jlong = 2*ilong;
      iofdat = (int_32 *)malloc(jlong*sizeof(int_32)); /* setara_c(jlong); */
      for (i = 0; i < jlong; ++i) iofdat[i] = idata[i];
      lcmppd_c(iplist,nom,ilong,itype,iofdat);
   }
}

void lcmpdl_c(lcm **iplist,int_32 iset,int_32 ilong,int_32 itype,int_32 *idata)
/*
 *----------------------------------------------------------------------
 *
 * COPY A BLOCK OF DATA FROM MEMORY INTO A LIST.
 *
 * INPUT PARAMETERS:
 *  IPLIST : ADDRESS OF THE LIST.
 *    ISET : POSITION OF THE SPECIFIC ELEMENT.
 *   ILONG : NUMBER OF INFORMATION ELEMENTS STORED IN THE CURRENT BLOCK.
 *   ITYPE : TYPE OF INFORMATION ELEMENTS STORED IN THE CURRENT BLOCK.
 *           0: DIRECTORY                1: INTEGER
 *           2: SINGLE PRECISION         3: CHARACTER*4
 *           4: DOUBLE PRECISION         5: LOGICAL
 *           6: COMPLEX                 99: UNDEFINED
 *   IDATA : INFORMATION ELEMENTS. DIMENSION IDATA(ILONG)
 *
 *----------------------------------------------------------------------
 */
{
   if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *ipxsm = (xsm *)*iplist + iset;
      xsmput_c(&ipxsm," ",ilong,itype,idata);
   } else {
      int_32 i, *iofdat;
      int_32 jlong = ilong;
      if (itype == 4 || itype == 6) jlong = 2*ilong;
      iofdat = (int_32 *)malloc(jlong*sizeof(int_32)); /* setara_c(jlong); */
      for (i = 0; i < jlong; ++i) iofdat[i] = idata[i];
      lcmppl_c(iplist,iset,ilong,itype,iofdat);
   }
}

void lcmpcd_c(lcm **iplist,const char *namp,int_32 ilong,char *hdata[])
/*
 *----------------------------------------------------------------------
 *
 * COPY AN ARRAY OF C STRING VARIABLES FROM MEMORY INTO A TABLE.
 *
 * INPUT PARAMETERS:
 *  IPLIST : ADDRESS OF THE TABLE.
 *    NAMP : CHARACTER*12 NAME OF THE BLOCK.
 *   ILONG : DIMENSION OF THE STRING ARRAY.
 *   HDATA : ARRAY OF ILONG STRINGS.
 *
 *----------------------------------------------------------------------
 */
{
   int_32 iset;
   lcm *jplist;
   jplist = lcmlid_c(iplist, namp, ilong);
   for (iset=0; iset<ilong; iset++) {
      int_32 i, ilen, *iofset;
      ilen = (strlen(hdata[iset]) + 4 ) / 4;
      iofset = (int_32 *)malloc(ilen*sizeof(int_32)); /* setara_c(ilen); */
      for (i=0; i<ilen; i++) strncpy ((char *)(iofset+i), hdata[iset]+4*i, 4);
      lcmppl_c(&jplist, iset, ilen, 3, iofset);
   }
}

void lcmgcd_c(lcm **iplist,const char *namp,char *hdata[])
/*
 *-----------------------------------------------------------------------
 *
 * COPY AN ARRAY OF C STRING VARIABLES FROM A TABLE INTO MEMORY.
 *
 * INPUT PARAMETERS:
 *  IPLIST : ADDRESS OF THE TABLE.
 *    NAMP : CHARACTER*12 NAME OF THE EXISTING BLOCK.
 *
 * OUTPUT PARAMETER:
 *   HDATA : ARRAY OF ILONG STRINGS (ALLOCATED BY LCMGCD).
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmgcd_c";
   int_32 iset, ilong, itylcm;
   lcm *jplist;
   lcmlen_c(iplist, namp, &ilong, &itylcm);
   if (itylcm != 10) {
      sprintf(AbortString,"%s: LIST EXPECTED.",nomsub);
      xabort_c(AbortString);
   }
   jplist = lcmgid_c(iplist, namp);
   for (iset=0; iset<ilong; iset++) {
      int_32 j, ilcmg, *iass;
      lcmlel_c(&jplist, iset, &ilcmg, &itylcm);
      iass = (int_32 *)malloc(ilcmg*sizeof(int_32)); /* setara_c(ilcmg); */
      lcmgdl_c(&jplist, iset, iass);
      hdata[iset] = (char *)malloc((int)4*ilcmg+1);
      for (j=0; j<ilcmg; j++) strncpy ((hdata[iset]+4*j),(char *) (iass + j), 4);
      hdata[iset][4*ilcmg]=' ';
      free(iass); /* rlsara_c(iass); */
      for(j=4*ilcmg; j>0; j--) {
         if (hdata[iset][j] != ' ') break;
         hdata[iset][j]='\0';
      }
   }
}

void lcmpcl_c(lcm **iplist,int_32 iset,int_32 ilong,char *hdata[])
/*
 *----------------------------------------------------------------------
 *
 * COPY AN ARRAY OF C STRING VARIABLES FROM MEMORY INTO A LIST.
 *
 * INPUT PARAMETERS:
 *  IPLIST : ADDRESS OF THE TABLE.
 *    ISET : POSITION OF THE BLOCK IN THE LIST.
 *   ILONG : DIMENSION OF THE CHARACTER VARIABLE.
 *   HDATA : ARRAY OF ILONG STRINGS.
 *
 *----------------------------------------------------------------------
 */
{
   int_32 jset;
   lcm *jplist;
   jplist=lcmlil_c(iplist, iset, ilong);
   for (jset=0; jset<ilong; jset++) {
      int_32 i, ilen, *iofset;
      ilen = (strlen(hdata[jset]) + 4 ) / 4;
      iofset = (int_32 *)malloc(ilen*sizeof(int_32)); /* setara_c(ilen); */
      for (i=0; i<ilen; i++) strncpy ((char *)(iofset+i), hdata[jset]+4*i, 4);
      lcmppl_c(&jplist, jset, ilen, 3, iofset);
   }
}

void lcmgcl_c(lcm **iplist,int_32 iset,char *hdata[])
/*
 *-----------------------------------------------------------------------
 *
 * COPY AN ARRAY OF C STRING VARIABLES FROM A LIST INTO MEMORY.
 *
 * INPUT PARAMETERS:
 *  IPLIST : ADDRESS OF THE TABLE.
 *    ISET : POSITION OF THE BLOCK IN THE LIST.
 *
 * OUTPUT PARAMETER:
 *   HDATA : ARRAY OF ILONG STRINGS (ALLOCATED BY LCMGCL).
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmgcl_c";
   int_32 jset, ilong, itylcm;
   lcm *jplist;
   lcmlel_c(iplist, iset, &ilong, &itylcm);
   if (itylcm != 10) {
      sprintf(AbortString,"%s: LIST EXPECTED.",nomsub);
      xabort_c(AbortString);
   }
   jplist = lcmgil_c(iplist, iset);
   for (jset=0; jset<ilong; jset++) {
      int_32 j, ilcmg, *iass;
      lcmlel_c(&jplist, jset, &ilcmg, &itylcm);
      iass = (int_32 *)malloc(ilcmg*sizeof(int_32)); /* setara_c(ilcmg); */
      lcmgdl_c(&jplist, jset, iass);
      hdata[jset] = (char *)malloc((int)4*ilcmg+1);
      for (j=0; j<ilcmg; j++) strncpy ((hdata[jset]+4*j),(char *) (iass + j), 4);
      hdata[jset][4*ilcmg]=' ';
      free(iass); /* rlsara_c(iass); */
      for(j=4*ilcmg; j>0; j--) {
         if (hdata[jset][j] != ' ') break;
         hdata[jset][j]='\0';
      }
   }
}

void lcmpsd_c(lcm **iplist,const char *namp,char *hdata)
/*
 *----------------------------------------------------------------------
 *
 * COPY A SINGLE C STRING VARIABLE FROM MEMORY INTO A TABLE.
 *
 * INPUT PARAMETERS:
 *  IPLIST : ADDRESS OF THE TABLE.
 *    NAMP : CHARACTER*12 NAME OF THE BLOCK.
 *   HDATA : C STRING.
 *
 *----------------------------------------------------------------------
 */
{
   int_32 i, ilong, *iofset;
   ilong = (strlen(hdata) + 4 ) / 4;
   iofset = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
   for (i=0; i<ilong; i++) strncpy ((char *)(iofset+i), hdata+4*i, 4);
   lcmppd_c(iplist, namp, ilong, 3, iofset);
}

char * lcmgsd_c(lcm **iplist,const char *namp)
/*
 *-----------------------------------------------------------------------
 *
 * COPY A SINGLE C STRING VARIABLE FROM A TABLE INTO MEMORY.
 *
 * INPUT PARAMETERS:
 *  IPLIST : ADDRESS OF THE TABLE.
 *    NAMP : CHARACTER*12 NAME OF THE EXISTING BLOCK.
 *
 * OUTPUT PARAMETER:
 *   lcmgsd_c : C STRING.
 *
 *-----------------------------------------------------------------------
 */
{
   static char nomstatic[133];
   char *nomsub="lcmgsd_c";
   int_32 i, ilong, itylcm, *iass;
   lcmlen_c(iplist, namp, &ilong, &itylcm);
   if (itylcm != 3) {
      sprintf(AbortString,"%s: CHARACTER DATA EXPECTED.",nomsub);
      xabort_c(AbortString);
   } else if (ilong*4 > 132) {
      sprintf(AbortString,"%s: CHARACTER DATA OVERFLOW.",nomsub);
      xabort_c(AbortString);
   }

   iass = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
   lcmget_c(iplist, namp, iass);
   for (i=0; i<ilong; i++) strncpy ((nomstatic+4*i),(char *) (iass+i), 4);
   nomstatic[ilong*4] = ' ';
   free(iass); /* rlsara_c(iass); */
   for(i=ilong*4; i>0; i--) {
      if(nomstatic[i] != ' ') break;
      nomstatic[i] = '\0';
   }
   return nomstatic;
}

void lcmpsl_c(lcm **iplist,int_32 iset,char *hdata)
/*
 *----------------------------------------------------------------------
 *
 * COPY A SINGLE C STRING VARIABLE FROM MEMORY INTO A LIST.
 *
 * INPUT PARAMETERS:
 *  IPLIST : ADDRESS OF THE TABLE.
 *    ISET : POSITION OF THE BLOCK IN THE LIST.
 *   HDATA : C STRING.
 *
 *----------------------------------------------------------------------
 */
{
   int_32 i, ilong, *iofset;
   ilong = (strlen(hdata) + 4 ) / 4;
   iofset = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
   for (i=0; i<ilong; i++) strncpy ((char *)(iofset+i), hdata+4*i, 4);
   lcmppl_c(iplist, iset, ilong, 3, iofset);
}

char * lcmgsl_c(lcm **iplist,int_32 iset)
/*
 *-----------------------------------------------------------------------
 *
 * COPY A SINGLE C STRING VARIABLE FROM A LIST INTO MEMORY.
 *
 * INPUT PARAMETERS:
 *  IPLIST : ADDRESS OF THE TABLE.
 *    ISET : POSITION OF THE BLOCK IN THE LIST.
 *
 * OUTPUT PARAMETER:
 *   lcmgsd_c : C STRING.
 *
 *-----------------------------------------------------------------------
 */
{
   static char nomstatic[133];
   char *nomsub="lcmgsl_c";
   int_32 i, ilong, itylcm, *iass;
   lcmlel_c(iplist, iset, &ilong, &itylcm);
   if (itylcm != 3) {
      sprintf(AbortString,"%s: CHARACTER DATA EXPECTED.",nomsub);
      xabort_c(AbortString);
   } else if (ilong*4 > 132) {
      sprintf(AbortString,"%s: CHARACTER DATA OVERFLOW.",nomsub);
      xabort_c(AbortString);
   }

   iass = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
   lcmgdl_c(iplist, iset, iass);
   for (i=0; i<ilong; i++) strncpy ((nomstatic+4*i),(char *) (iass+i), 4);
   nomstatic[ilong*4] = ' ';
   free(iass); /* rlsara_c(iass); */
   for(i=ilong*4; i>0; i--) {
      if(nomstatic[i] != ' ') break;
      nomstatic[i] = '\0';
   }
   return nomstatic;
}
