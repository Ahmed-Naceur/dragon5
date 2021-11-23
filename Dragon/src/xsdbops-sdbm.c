/*
 *     $Id: xsdbops-sdbm.c,v 1.4 1995/12/12 22:21:27 laughton Exp laughton $
 *
 *---------------------------------------------------------------------------
 *
 *      Peter J. Laughton
 *      AECL
 *      Chalk River Laboratories
 *      Chalk River, Ontario
 *      CANADA  K0J 1J0
 *
 *      Phone:  (613) 584-8811, extension 4267
 *      FAX:    (613) 584-1108
 *
 *      Internet:  laughtonp@crl.aecl.ca
 *
 *---------------------------------------------------------------------------
 *
 *      Revision history (as of 1995 October 25):
 *
 *      $Log:	xsdbops-sdbm.c,v $
 * Revision 1.5  96/06/28  13:04:34  13:04:34  laughton
 * support for NDAS file conversion started
 * 
 *      Revision 1.4  1995/12/12 22:21:27  laughton
 *      more GDBM updates
 *
 * Revision 1.3  1995/12/04  14:46:02  laughton
 * continuing development
 *
 * Revision 1.2  1995/11/30  20:26:36  laughton
 * XSDB files now opened in this file
 *
 * Revision 1.1  1995/11/07  16:36:31  laughton
 * Initial revision
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xsdbops.h"
#include "sdbm.h"

#define KEYBUFSIZE 200
static char keybuf[KEYBUFSIZE];

static datum key;

void xsdbReadInit(char *flist, int *status) {
     key.dptr=keybuf;		/* use the same buffer for all keys */
     *status=initRead(flist);
}

int xsdbCountRecs(int *status, int dbFileIndex) {
     int n;
     *status=countRecs(&n,dbFileIndex);
     return n;
}

void xsdbWriteInit(int *status) {
     key.dptr=keybuf;		/* use the same buffer for all keys */
     *status=initWrite();
}

void xsdbCloseWrite() {
     closeSDBWrite();
}

void formItemKey(char *itemName) {
     strcpy(key.dptr,itemName);
     key.dsize=strlen(key.dptr)+1;
}

void formNuclideItemKey(char *nuclideName, char *itemName) {
     strcpy(key.dptr,nuclideName);
     strcat(key.dptr,"/");
     strcat(key.dptr,itemName);
     key.dsize=strlen(key.dptr)+1;
}

void formNuclideTempItemKey(char *nuclideName, char *itemName, float temp) {
     int endOffset;
     strcpy(key.dptr,nuclideName);
     strcat(key.dptr,"/");
     strcat(key.dptr,itemName);
     strcat(key.dptr,"/");
     endOffset=strlen(key.dptr);
     sprintf(key.dptr+endOffset,"%.2f",temp);
     key.dsize=strlen(key.dptr)+1;
}

void formNuclideGroupItemKey(char *nuclideName, char *itemName, int g) {
     int endOffset;
     strcpy(key.dptr,nuclideName);
     strcat(key.dptr,"/");
     strcat(key.dptr,itemName);
     strcat(key.dptr,"/");
     endOffset=strlen(key.dptr);
     sprintf(key.dptr+endOffset,"G%d",g);
     key.dsize=strlen(key.dptr)+1;
}



void xsdbStoreItem(char *itemName, void *from, 
		   int nbytes, int dataType, int verbose) {

     int failure;

     formItemKey(itemName);
     failure=writeRecord(key.dptr, from, nbytes, dataType);
     if (failure) {
	  fprintf(stderr,"error in xsdbStoreItem()\n");
	  exit(1);
     }
     if (verbose)
	  printf("Wrote SDBM record with length %6d bytes, type code %3d -- %s\n",
		 nbytes, dataType, key.dptr);
}

void xsdbStoreNuclideItem(char *nuclideName, 
			  char *itemName, void *from, int nbytes, 
			  int dataType, int verbose) {

     formNuclideItemKey(nuclideName,itemName);
     xsdbStoreItem(keybuf,from,nbytes,dataType,verbose);
}

void xsdbStoreNuclideTempItem(char *nuclideName, 
			      char *itemName, float temp, void *from, 
			      int nbytes, int dataType, int verbose) {

     formNuclideTempItemKey(nuclideName,itemName,temp);
     xsdbStoreItem(keybuf,from,nbytes,dataType,verbose);
}

void xsdbStoreNuclideGroupItem(char *nuclideName, 
			       char *itemName, int g, void *from, int nbytes, 
			       int dataType, int verbose) {

     formNuclideGroupItemKey(nuclideName,itemName,g);
     xsdbStoreItem(keybuf,from,nbytes,dataType,verbose);
}

/* ~~~~~~~~~~~~~~~~~~ */
/* Retrieval Routines */
/* ~~~~~~~~~~~~~~~~~~ */

/* The following routines make no check as to whether a record with
   the specified key data exists, and run-time failure will result if
   no such record is found. */

void *xsdbRetrieveItem(char *itemName, int *nbytes, int verbose) {

     datum record;

     formItemKey(itemName);
     record=readRecord(keybuf);

     *nbytes=record.dsize;
     return record.dptr;
}

void *xsdbRetrieveNuclideItem(char *nuclideName, 
			      char *itemName, int *nbytes, int verbose) {

     formNuclideItemKey(nuclideName,itemName);
     return xsdbRetrieveItem(keybuf,nbytes,verbose);
}

void *xsdbRetrieveNuclideTempItem(char *nuclideName, 
				  char *itemName, float temp, int *nbytes, int verbose) {

     formNuclideTempItemKey(nuclideName,itemName,temp);
     return xsdbRetrieveItem(keybuf,nbytes,verbose);
}

void *xsdbRetrieveNuclideGroupItem(char *nuclideName, 
				   char *itemName, int g, int *nbytes, int verbose) {

     formNuclideGroupItemKey(nuclideName,itemName,g);
     return xsdbRetrieveItem(keybuf,nbytes,verbose);
}
