/*
 *
 *     $Id: sdbm.h,v 1.5 1995/12/12 22:21:27 laughton Exp laughton $
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
 *      $Log:	sdbm.h,v $
 * Revision 1.6  96/06/28  13:10:19  13:10:19  laughton
 * support for NDAS file conversion started.
 * 
 *      Revision 1.5  1995/12/12 22:21:27  laughton
 *      more GDBM updates
 *
 * Revision 1.4  1995/12/12  19:59:38  laughton
 * updated the GDBM routines
 *
 * Revision 1.3  1995/12/04  14:46:47  laughton
 * continuing development
 *
 * Revision 1.2  1995/11/30  20:48:03  laughton
 * merged with main trunk
 *
 * Revision 1.1.1.2  1995/11/30  20:31:40  laughton
 * continuing development
 *
 * Revision 1.1.1.1  1995/11/07  16:41:01  laughton
 * side-branch for experiment
 *
 * Revision 1.1  1995/11/07  16:32:06  laughton
 * Initial revision
 *
 */

#define FAIL (-1)
#define FileFull 1 
#define InvalidKey 2
#define OpenFailure 3
#define BadFile 4
#define MemoryAllocFailure 5
#define FixedLimitExceeded 6

typedef struct {
     char *dptr;
     int   dsize;
} datum;

extern void closeSDBWrite();

extern void closeSDBRead();

extern int writeRecord(char *recordName, 
		       void *data, int nbytes, int dataType);

extern int initRead(char *flist);

extern int countRecs(int *nRecs, int dbFile);

extern int initWrite();

extern datum readRecord(char *recordName);

extern datum readIndexedRecord(int keyIndex, int dbFileIndex,
			       char *recordKey, int *dataType);

extern void emitTitles();
