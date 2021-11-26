/*
 *
 *     $Id: sdbm.c,v 1.10 1997/05/21 20:47:59 laughton Exp laughton $
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
 *      $Log:	sdbm.c,v $
 * Revision 1.11  99/03/11  17:03:56  17:03:56  laughton
 * proper type-dependent calloc calls now used
 * 
 * Revision 1.10  1997/05/21  20:47:59  laughton
 * maximum number of records in one file doubled
 *
 *      Revision 1.9  1996/12/19 19:34:56  laughton
 *      Intel-x86/Windows95 support added (the list sep is a ';', not a ':')
 *
 *      Revision 1.8  1996/06/28 13:04:34  laughton
 *      support for NDAS file conversion started
 *
 *      Revision 1.7  1995/12/12 22:21:27  laughton
 *      more GDBM updates
 *
 * Revision 1.6  1995/12/07  20:10:45  laughton
 * minor SGI compiler warnings rectified
 *
 * Revision 1.5  1995/12/05  21:08:47  laughton
 * emitTitle bug fixed.
 *
 * Revision 1.4  1995/12/04  14:46:12  laughton
 * continuing development
 *
 * Revision 1.3  1995/12/01  21:27:19  laughton
 * more error handling and identification of library titles
 *
 * Revision 1.2  1995/11/30  20:22:37  laughton
 * continuing development
 *
 * Revision 1.1  1995/11/07  16:32:06  laughton
 * Initial revision
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "sdbm.h"

#define MAXKEYLENGTH 64
#define StringTypeC 8
#define MAKE_ARRAY(thing,number) \
               ((thing *) calloc((unsigned) (number),sizeof(thing)))

static int magicNumber=0x198802ab; /* file-type marker */

typedef struct {
  char key[MAXKEYLENGTH+1];	/* string of printable characters */
  int dataType;		/* numeric flag indicating real, integer, etc */
  int size;			/* in bytes */
  int dataOffset;	/* from the beginning of the file */
} HeaderElement;

typedef struct {
  char *filename;
  FILE *filePointer;
  int nValidRecords;		/* > 0 means active */
  HeaderElement *header;
} SDBFile;

#define size_HE (sizeof(HeaderElement))

#define MAXDBFILES 10
#define MAXFNLENGTH 100
static SDBFile dbFile[MAXDBFILES];
static int dbfCount=0;

static int initReadF=0;

/* =============================================== */
/* Prototypes of static functions appearing below: */
/* =============================================== */

static int installDBF(char *filename);
static int compareHeaderElements(const void *h1, const void *h2);
static int searchHeader(HeaderElement *header, char *searchKey, int n);

/* =============================================== */

int initRead(char *flist) {
  int i;
  char nbuf[MAXFNLENGTH+1];
  int fc, lc;
  int retval;
  char flistSep;

#ifdef PC_TargetCPU
  flistSep=';';
#else
  flistSep=':';
#endif

  initReadF=1;

  for (i=0; i<MAXDBFILES; i++) {
    dbFile[i].nValidRecords=0;
  }
     
  fc=0;
  lc=0;
  while (flist[fc]!=0) {		/* end of string reached? */
    while (flist[fc]==flistSep) fc++;
    if (flist[fc]=='\0')	/* at the end of the list? */
      break;

    lc=fc;
    while(flist[lc]!=flistSep 
	  && flist[lc] != 0) lc++;

    if (lc-fc > MAXFNLENGTH) 
      return FixedLimitExceeded;
	       
    strncpy(nbuf, flist+fc, lc-fc);
    nbuf[lc-fc]='\0';
    retval=installDBF(nbuf);

    if (retval) 
      return retval;	/* non-zero is an error code */

    fc=lc;
  }

  emitTitles();
  return 0;			/* success */
}
     
int countRecs(int *nRecs, int dbFileIndex) {
  *nRecs=0;
  if (dbFileIndex < 0 || dbFileIndex > MAXDBFILES-1)
    return BadFile;

  *nRecs=dbFile[dbFileIndex].nValidRecords;
  return 0;			/* success */
}

static int installDBF(char *filename) {
  int mn=0;
  int i;

  /* Later, add a check to make sure this is not a new file. */

  if (dbfCount >= MAXDBFILES) 
    return FixedLimitExceeded;

  dbFile[dbfCount].filename=strdup(filename);
  if (!dbFile[dbfCount].filename)
    return MemoryAllocFailure;

  dbFile[dbfCount].filePointer=fopen(filename,"rb");
  if (!dbFile[dbfCount].filePointer)
    return OpenFailure;
	  
  /* Check that the magic number is present. */
  fread(&mn,sizeof(int),1,dbFile[dbfCount].filePointer);
  if (mn!=magicNumber)
    return BadFile;

  fread(&dbFile[dbfCount].nValidRecords, sizeof(int), 1,
	dbFile[dbfCount].filePointer);

  dbFile[dbfCount].header=MAKE_ARRAY(HeaderElement,
				     dbFile[dbfCount].nValidRecords);
  if (!dbFile[dbfCount].header)
    return MemoryAllocFailure;

  /* Load the header into memory (for now) */
  for (i=0; i < dbFile[dbfCount].nValidRecords; i++)
    fread(&dbFile[dbfCount].header[i], size_HE, 1,
	  dbFile[dbfCount].filePointer);

  dbfCount++;
  return 0;			/* success */
}

static int compareHeaderElements(const void *h1, const void *h2) {
  return strcmp(((HeaderElement *)h1)->key, ((HeaderElement *)h2)->key);
}

void closeSDBRead() {
  int i,irc;

  for (i=0; i<dbfCount; i++) {
    irc = fclose(dbFile[i].filePointer);
    if (irc != 0) perror ("close error for the ndas binary file: ");
  }
  dbfCount = 0;
}

#ifdef SequentialSearch
/* For now, make a brute-force sequential search. */

static int searchHeader(char *searchKey) {
  int i;

  for (i=0; i<nValidRecords; i++)
    if (!strcmp (header[i].key, searchKey))
      return i;

  return FAIL;
}
#endif

/* Binary search -- see Knuth, volume 3, page 407. */

static int searchHeader(HeaderElement *header, char *searchKey, int n) {
  int l, u, midPoint;
  int sign;

  l=0; 
  u=n-1;

  while (u >= l) {
    midPoint=(l+u)/2;	/* floor results automatically */
    sign=strcmp(searchKey,header[midPoint].key);

    if (sign < 0)
      u=midPoint-1;
    else if (sign > 0)
      l=midPoint+1;
    else
      return midPoint;
  }

  return FAIL;
}

datum readRecord(char *recordName) {

  int keyIndex;
  datum record;
  int nbytes;
  int i;
  int found=0;

  /* Failure is indicated if these values remain intact on exit. */
  record.dptr=0;
  record.dsize=0;

  if (!initReadF)
    return record;

  for (i=0; i<dbfCount; i++) {
    keyIndex=searchHeader(dbFile[i].header, recordName, dbFile[i].nValidRecords);
    if (keyIndex!=FAIL) {
      found=1;
      break;
    }
  }

  if (!found) {
    printf("\n XSDB error: record not found -- %s\n",recordName);
    fflush(stdout);
    exit(1);
  }

  nbytes=dbFile[i].header[keyIndex].size;

  if (dbFile[i].header[keyIndex].dataType == StringTypeC)
       record.dptr=calloc(nbytes, 1);
  else
       record.dptr=calloc(nbytes >> 2, 4);

  if (!record.dptr)
    return record;

  record.dsize=nbytes;

  /* Jump to the location in the file and read the record into memory. */
  fseek(dbFile[i].filePointer, dbFile[i].header[keyIndex].dataOffset, SEEK_SET);
  fread((char*)record.dptr, 1, nbytes, dbFile[i].filePointer);

  return record;
}

datum readIndexedRecord(int keyIndex, int dbFileIndex, char *recordKey, int *dataType) {

  datum record;
  int nbytes;

  /* Failure is indicated if these values remain intact on exit. */
  record.dptr=0;
  record.dsize=0;

  if (dbFileIndex < 0 || dbFileIndex > MAXDBFILES-1)
    return record;

  if (!initReadF)
    return record;

  if (keyIndex >= dbFile[dbFileIndex].nValidRecords)
    return record;

  strcpy(recordKey, dbFile[dbFileIndex].header[keyIndex].key);
  *dataType=dbFile[dbFileIndex].header[keyIndex].dataType;

  nbytes=dbFile[dbFileIndex].header[keyIndex].size;

  if (*dataType == StringTypeC)
       record.dptr=calloc(nbytes, 1);
  else
       record.dptr=calloc(nbytes >> 2, 4);

  if (!record.dptr)
    return record;

  record.dsize=nbytes;

  /* Jump to the location in the file and read the record into memory. */
  fseek(dbFile[dbFileIndex].filePointer, 
	dbFile[dbFileIndex].header[keyIndex].dataOffset, 
	SEEK_SET);

  fread((char*)record.dptr, 1, nbytes, dbFile[dbFileIndex].filePointer);

  return record;
}

void emitTitles() {
  int keyIndex;
  datum record;
  int i;

  for (i=0; i<dbfCount; i++) {
    keyIndex=searchHeader(dbFile[i].header, "XSDBTitle", dbFile[i].nValidRecords);
    if (keyIndex!=FAIL) {
      record.dsize=dbFile[i].header[keyIndex].size;
      record.dptr=calloc(record.dsize,1);
      /* The title must include a trailing zero. */

      if (record.dptr) {
	fseek(dbFile[i].filePointer, 
	      dbFile[i].header[keyIndex].dataOffset, SEEK_SET);
	fread((char*)record.dptr, 1, record.dsize, dbFile[i].filePointer);
        printf("\n Library %d -- %s\n",i+1,record.dptr);
        fflush(stdout);
	free(record.dptr);
      } else {
	fprintf(stderr,"memory allocaiton failure\n");
	exit(1);
      }
    } else {
      printf(" Library %d -- ***Untitled***\n",i+1);
      fflush(stdout);
    }
  }
}

/* Writing: */

static int nValidRecords=0;

#define MAXRECORDS 10000
HeaderElement header[MAXRECORDS];

static int initWriteF=0;
static long int nextDataOffset;
static int writeRecordIndex;
static FILE *sdbFile=0;

int initWrite() {
  int i, j;
  long int headerSize;
  char *filename="new.sdb";

  writeRecordIndex=0;
  headerSize=MAXRECORDS*size_HE+2*sizeof(int);

  /* Later add a check to make sure this is a new file. */

  sdbFile=fopen(filename,"wb");
  if (!sdbFile)
    return OpenFailure;
	  
  fseek(sdbFile,headerSize,SEEK_END); /* extend the file */
  nextDataOffset=headerSize;

  /* Preset all keys to 0. */
  for (i=0; i < MAXRECORDS; i++)
    for (j=0; j < MAXKEYLENGTH+1; j++)
      header[i].key[j]='\0';

  return 0;			/* success */
}

void closeSDBWrite() {
  long int offset;
  int i,irc;

  /* Sort the header elements into some order.  (Alphabetical for now.) */

  nValidRecords=writeRecordIndex;
  qsort(header,nValidRecords,size_HE,compareHeaderElements);

  /* Write the file header and then close the file. */
 
  offset=0;
  fseek(sdbFile,offset,SEEK_SET);
  fwrite(&magicNumber,sizeof(int),1,sdbFile);
  fwrite(&nValidRecords,sizeof(int),1,sdbFile); /* stores the record count */

  for (i=0; i<MAXRECORDS; i++) {
    fwrite(&header[i],size_HE,1,sdbFile);
  }

  irc = fclose(sdbFile);
  if (irc != 0) perror ("close error for the ndas binary file: ");
}

int writeRecord(char *recordName, 
		void *data, int nbytes, int dataType) {

  long int offset;
  int retval;

  if (!initWriteF) {
    initWriteF=1;
    retval=initWrite();
    if (retval!=0)
      return retval;
  }

  if (writeRecordIndex==MAXRECORDS) {
    return FileFull;
  }
     
  if (strlen(recordName) > MAXKEYLENGTH)
    return InvalidKey;

  strcat(header[writeRecordIndex].key,recordName); 
  header[writeRecordIndex].dataType=dataType;
  header[writeRecordIndex].size=nbytes;

  header[writeRecordIndex].dataOffset=nextDataOffset;
  offset=nbytes;
  fseek(sdbFile,offset,SEEK_END); /* extend the file */
  fseek(sdbFile,nextDataOffset,SEEK_SET); /* position for write */
  fwrite((char*)data, 1, nbytes, sdbFile);

  writeRecordIndex++;
  nextDataOffset+=offset;

  return 0;
}
