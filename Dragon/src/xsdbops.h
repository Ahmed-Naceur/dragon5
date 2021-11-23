/*
 *     $Id: xsdbops.h,v 1.6 1995/12/12 22:21:27 laughton Exp laughton $
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
 *      Revision history (as of 1995 August 9):
 *
 *      $Log:	xsdbops.h,v $
 * Revision 1.7  96/06/28  13:10:19  13:10:19  laughton
 * support for NDAS file conversion started.
 * 
 *      Revision 1.6  1995/12/12 22:21:27  laughton
 *      more GDBM updates
 *
 * Revision 1.5  1995/12/12  19:59:38  laughton
 * updated the GDBM routines
 *
 * Revision 1.4  1995/12/04  14:46:36  laughton
 * continuing development
 *
 * Revision 1.3  1995/11/30  20:31:14  laughton
 * XSDB files now opened in xsdbops files
 *
 * Revision 1.2  1995/10/10  18:26:19  laughton
 * continuing development
 *
 * Revision 1.1  1995/08/24  20:27:28  laughton
 * Initial revision
 *
 */

void xsdbReadInit(char *nomC, int *status);

int xsdbCountRecs(int *status, int dbFileIndex);

void xsdbWriteInit(int *status);

void xsdbCloseWrite();

void xsdbStoreItem(char *itemName, void *from, 
		   int nbytes, int dataType, int verbose);

void xsdbStoreNuclideItem(char *nuclideName, 
			  char *itemName, void *from, int nbytes, int dataType, int verbose);

void xsdbStoreNuclideTempItem(char *nuclideName, 
			      char *itemName, float temp, void *from, int nbytes, int dataType, int verbose);

void xsdbStoreNuclideGroupItem(char *nuclideName, 
			       char *itemName, int g, void *from, int nbytes, int dataType, int verbose);


/* ~~~~~~~~~~~~~~~~~~ */
/* Retrieval Routines */
/* ~~~~~~~~~~~~~~~~~~ */

/* The following routines make no check as to whether a record with
   the specified key data exists, and run-time failure will result if
   no such record is found. */


void *xsdbRetrieveItem(char *itemName, int *nbytes, int verbose);

void *xsdbRetrieveNuclideItem(char *nuclideName, 
			      char *itemName, int *nbytes, int verbose);

void *xsdbRetrieveNuclideTempItem(char *nuclideName, 
				  char *itemName, float temp, int *nbytes, int verbose);

void *xsdbRetrieveNuclideGroupItem(char *nuclideName, 
				   char *itemName, int g, int *nbytes, int verbose);
