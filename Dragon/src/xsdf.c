
/***************************************/
/* Fortran bindings for the NDAS C API */
/* Copyright: Peter J. Laughton, AECL  */
/***************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xsdbops.h"
#include "sdbm.h"
#include "xsdb-defs.h"

static char AbortString[132];
static int packIndex;
typedef enum {Int, Float, Char8} PackType;
static char *packedBurnupData=0;

typedef struct {
  int ib1, ib2;
  float rb1, rb2;
} BurnQuad;

typedef struct {
  /* From subinx.inc: NBURN, ISOID, AW, IAN, NFISS, NTEMP, NR,
     NSUBNK, NNA, NP1, NFSPEC, IENDFB */
     /* nomenclature of crnl-2866, page 3 */

  int nburn;
  int numericId;		
  float aw;
  int iz;
  int nf;
  int nt;
  int nr;
  int ndat2;			/* not in xs block on unf-seq file */
  int ndat3;                 /* not in xs block on unf-seq file */
  int np1;
  int ns;
  int iendfb;                /* not in xs block on unf-seq file */
  char name[9];
  int installed;             /* flagged true when all of this block has been loaded */
  BurnQuad *burnQuad;
} Nuclide;

Nuclide *nuclide;		/* array of all nuclides */

static int nLoadedNuclides=0;	/* number actually available */
static int lastNuclideAccessed=0; /* library index (1,...,nel) of last nuclide accessed */

typedef struct {
  int code;
  char *name;
} CodeNamePair;

static CodeNamePair itemMap[]={
  {BickleyFunctionTablesKi3, "BickleyFunctionTablesKi3"},
  {BickleyFunctionTablesKi35, "BickleyFunctionTablesKi35"},
  {BurnCount, "BurnCount"},
  {BurnInteger, "BurnInteger"},
  {BurnReal, "BurnReal"},
  {Absorption, "Absorption"},
  {Transport, "Transport"},
  {Fission, "Fission"},
  {NuFission, "NuFission"},
  {N2n, "N2n"},
  {FissionSpectrum, "FissionSpectrum"},
  {PotScatSlowingDown, "PotScatSlowingDown"},
  {LengthsThermalP0, "LengthsThermalP0"},
  {GCLambda, "GCLambda"},
  {Header, "Header"},
  {ResHeader, "ResHeader"},
  {LengthsScatP1, "LengthsScatP1"},
  {ScatP0, "ScatP0"},
  {ScatP1, "ScatP1"},
  {ThermalXSTemp, "ThermalXSTemp"},
  {ThermalP1Temp, "ThermalP1Temp"},
  {GroupBoundaries, "GroupBoundaries"},
  {PotScat, "PotScat"},
  {NJOYFlux, "NJOYFlux"},
  {Hequivalence, "Hequivalence"},
  {HeqHeader, "HeqHeader"},
  {TransportCorrection, "TransportCorrection"}
};

static char *unknownMessage="--unknown--";

#define NITEMS (sizeof(itemMap)/sizeof(CodeNamePair))

char *itemName(int code) {
  int i;
  for (i=0; i<NITEMS; i++) {
    if (itemMap[i].code == code)
      return itemMap[i].name;
  }
  return unknownMessage;
}

#ifdef VerboseXSDB
static int verbose=1;
#else
static int verbose=0;
#endif

typedef int DBItem;		/* fortran integer type */

static FILE *logFile=0;
static FILE *indexFile;

/* constant parameters local to this file */

static int nel;
static int ng, ng0, ng1, ng2, ng3; 
static int nnfpd;
static int p1NuclideCount;
static int fissileNuclideCount;
static int n1rc, m1rc, n1m1rc, nresmc, lsctfl;
static int jp0max, jp1max;

/* Prototypes of functions that appears below: */

int nuclideId(char *name);
int nuclideIndex(char *name);
int bcodeNumber(char *bname);
void unpackBD(void *ptr, PackType type);

#define CAPTURE  1
#define DECAY    2
#define FISSPROD 3
#define FISSION  4
#define N2N      5
#define min_macro(A,B)  ((A) < (B) ? (A) : (B))
#define MAKE_ARRAY(thing,number) \
               ((thing *) calloc((unsigned) (number),sizeof(thing)))

void xsdopn_c(char *nomC, int *status)
{
  int i;
  int endFlag;
  int nbytes;
  char *flist, *idxfn;
  logFile=fopen("xsdb.log","w");
  if (verbose) {
    fprintf(logFile,"In XSDOPN open file -->%s<--\n",nomC);
    fflush(logFile);
  }
  flist=strchr(nomC,':');
  if(!flist) {
    fprintf(logFile,"index file missing; namfil=%s\n",nomC);
    *status=OPEN_FAILURE;
    return;
  } 
  flist++;
  xsdbReadInit(flist,status);
  if (*status) 		/* non-zero means something is wrong */
    return;

  idxfn=strtok(nomC,":\n");

  indexFile=fopen(idxfn,"r");
  if (!indexFile) {
    perror(idxfn);
    fprintf(logFile,"open failure for index file %s\n",idxfn);
    *status=OPEN_FAILURE;
    return;
  }

  fscanf(indexFile,"%d",&nel);
  fscanf(indexFile,"%d",&ng);
  fscanf(indexFile,"%d",&ng0);
  fscanf(indexFile,"%d",&ng1);
  fscanf(indexFile,"%d",&ng2);
  fscanf(indexFile,"%d",&ng3);
  fscanf(indexFile,"%d",&fissileNuclideCount);
  fscanf(indexFile,"%d",&nnfpd);
  fscanf(indexFile,"%d",&p1NuclideCount);
  fscanf(indexFile,"%d",&nresmc);
  fscanf(indexFile,"%d",&n1rc);
  fscanf(indexFile,"%d",&m1rc);
  fscanf(indexFile,"%d",&n1m1rc);
  fscanf(indexFile,"%d",&lsctfl);
  fscanf(indexFile,"%d",&jp0max);
  fscanf(indexFile,"%d",&jp1max);

  nuclide=MAKE_ARRAY(Nuclide,nel);

  if (!nuclide) {
    fprintf(logFile,"error:  memory allocation failure\n");
    exit(1);
  }

  for (i=0; i<nel; i++) {
    fscanf(indexFile,"%s %d",nuclide[i].name,&nuclide[i].numericId);
    nuclide[i].installed=0;
    nuclide[i].nburn=0;  /* default length of burnup data record */
  }

  fclose(indexFile);

  packedBurnupData=xsdbRetrieveItem("BurnupData", &nbytes, verbose);

  packIndex=0;
  endFlag=0;
  do {
    char thisNuclideName[9];
    char type[9];
    int thisNburn;
    int sourceIndex;
    int bcode;

    if (packIndex==nbytes) {
      endFlag=1;
      continue;
    }

    unpackBD(thisNuclideName,Char8);
    unpackBD(&thisNburn,Int);

    if (thisNburn==0)	/* no burnup data follows */
      continue;
    sourceIndex=nuclideIndex(thisNuclideName);

    if (sourceIndex == -1) {
      *status=UnknownRequest;
      return;
    }

    nuclide[sourceIndex].nburn=thisNburn;
    nuclide[sourceIndex].burnQuad=MAKE_ARRAY(BurnQuad,thisNburn);

    if (!nuclide[sourceIndex].burnQuad) {
      fprintf(logFile,"error:  memory allocation failure\n");
      exit(1);
    }

    for (i=0; i<thisNburn; i++) {
      unpackBD(type,Char8);
      bcode=bcodeNumber(type);

      if (bcode== -1) {
	*status=UnknownRequest;
	return;
      }

      nuclide[sourceIndex].burnQuad[i].ib2=bcode;

      if (bcode==FISSION) {
	unpackBD(&nuclide[sourceIndex].burnQuad[i].ib1,Int);
	unpackBD(&nuclide[sourceIndex].burnQuad[i].rb1,Float);
	unpackBD(&nuclide[sourceIndex].burnQuad[i].rb2,Float);
      } else {
	unpackBD(thisNuclideName,Char8);
	unpackBD(&nuclide[sourceIndex].burnQuad[i].rb1,Float);
	unpackBD(&nuclide[sourceIndex].burnQuad[i].rb2,Float);

	nuclide[sourceIndex].burnQuad[i].ib1=nuclideId(thisNuclideName);

	if (nuclide[sourceIndex].burnQuad[i].ib1== -1) {
	  *status=UnknownNuclide;
	  return;
	}
      }
    }
  } while (!endFlag);

  free(packedBurnupData);

  *status=NORMAL;
}

void xsdcl_c()
{
  int irc;
  closeSDBRead();
  irc = fclose(logFile);
  if (irc != 0) perror ("close error of xsdb.log ");  
  }

void xsdnam_c(int *iset, int *numericId, char *isonam, int *status)
{
  char *nomsub="XSDNAM";
  char *cp;
  int j, len;

  if (verbose) {
    fprintf(logFile,"In xsdbnam...\n");
    fflush(logFile);
  }

  if (*iset-1 > nel) {
    sprintf(AbortString,"%s: Insufficent allocation to hold isotope names",nomsub);
    *status=RECORD_INDEX_OVERFLOW;
    return;
  }

  if (verbose) {
    fprintf(logFile,"%10s %d\n",nuclide[*iset-1].name,nuclide[*iset-1].numericId);
    fflush(logFile);
  }
  *numericId=nuclide[*iset-1].numericId;
  
  len=strlen(nuclide[*iset-1].name);
  cp=isonam;
  for (j=0; j<len; j++) {
    *cp=nuclide[*iset-1].name[j];
    cp++;
  }
  *cp='\0';

  *status=NORMAL;
}

void xsdbld_c(DBItem *item, int *where, int *status) {
  float *fp, *fpTarget;
  int *ip, *ipTarget;
  int nbytes;
  int i;

  if (verbose) {
    fprintf(logFile,"In xsdbl with item %25s and location %016llX\n",
	    itemName(*item),(long long)where);
    fflush(logFile);
  }

  switch (*item) {

  case FissionSpectrum:	/* library default fission spectrum */
    fp=xsdbRetrieveItem("FissionSpectrum", &nbytes, verbose);
    if (nbytes/4 != ng0) {
      fprintf(logFile,"warning:  fission-spectrum vector length mismatch\n");
    }
    fpTarget=(float *) where;
    for (i=0, fpTarget=(float *)where; i<ng0; i++, fpTarget++) {
      *fpTarget=fp[i];
    }
    free(fp);
    break;

  case Header:
    ipTarget=(int *) where;
    ipTarget[0]=nel;
    ipTarget[1]=ng;
    ipTarget[2]=ng0;
    ipTarget[3]=ng1;
    ipTarget[4]=ng2;
    ipTarget[5]=ng3;
    ipTarget[6]=fissileNuclideCount;
    ipTarget[7]=nnfpd;
    ipTarget[8]=p1NuclideCount;
    ipTarget[9]=nresmc;
    ipTarget[10]=n1rc;
    ipTarget[11]=m1rc;
    ipTarget[12]=n1m1rc;
    ipTarget[13]=lsctfl;
    ipTarget[14]=jp0max;
    ipTarget[15]=jp1max;
    break;

  case ResHeader:
    ip=xsdbRetrieveItem("ResHeader", &nbytes, verbose);
    ipTarget=(int *) where;
    for (i=0; i<nbytes/4; i++) {
      ipTarget[i]=ip[i];
    }
    free(ip);
    break;

  case GroupBoundaries:
    fp=xsdbRetrieveItem("GroupBoundaries", &nbytes, verbose);
    if (nbytes/4 != ng+1) {
      fprintf(logFile,"warning:  group-boundary vector length mismatch\n");
    }
    fpTarget=(float *) where;
    for (i=0, fpTarget=(float *)where; i<ng+1; i++, fpTarget++) {
      *fpTarget=fp[i];
    }
    free(fp);
    break;

  default:
    fprintf(logFile,"Serious error: unknown xsdb request %d!\n",*item);
    *status=UnknownRequest;
    return;
  }
  *status=NORMAL;
}

void xsdiso_c(int *groupRange, DBItem *item, int *nuclideIndex, int *where, int *status)
{
  float *fp, *fpTarget;
  int *ip, *ipTarget;
  int nbytes;
  int i;
  int j;
  fpTarget=(float *)where;
  ipTarget=(int *)where;
  j=*nuclideIndex;

  if (j==LastNuclideAccessed)
    j=lastNuclideAccessed;
  else
    lastNuclideAccessed=j;
  j--;			/* offset now */

  if (verbose) {
    if (nuclide[j].installed)
      fprintf(logFile,"In xsdbloadnuclide with item %25s, nuclide %10s (%5d) and location %016llX\n",
	      itemName(*item),nuclide[j].name,j+1,(long long)where);
    else
      fprintf(logFile,"In xsdbloadnuclide with item %25s, nuclide %5d and location %016llX\n",
	      itemName(*item),j+1,(long long)where);
    fflush(logFile);
  }

  switch (*item) {

  case BurnCount:
    ipTarget[0]=nuclide[j].nburn;
    break;

  case BurnInteger:
    if (verbose)
      fprintf(logFile,"Nuclide index %3d (%s)\n",j,nuclide[j].name);
    for (i=0; i<nuclide[j].nburn; i++) {
      ipTarget[i*2]=nuclide[j].burnQuad[i].ib1;
      ipTarget[i*2+1]=nuclide[j].burnQuad[i].ib2;
    }
    break;

  case BurnReal:
    for (i=0; i<nuclide[j].nburn; i++) {
      fpTarget[i*2]=nuclide[j].burnQuad[i].rb1;
      fpTarget[i*2+1]=nuclide[j].burnQuad[i].rb2;
    }
    break;

  case GCLambda:
    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "GCLambda/Res",&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i]; 
    free(fp);
    break;

  case Transport:

    switch(*groupRange) {

    case FastRes:
      fp=xsdbRetrieveNuclideItem(nuclide[j].name,
				 "Transport/Epithermal",&nbytes,verbose);
      for (i=0; i< nbytes/4; i++) 
	fpTarget[i]=fp[i];
      break;

    default:
      fprintf(logFile,"error:  unknown groupRange request\n");
      exit(1);
    }
    free(fp);
    break;

  case Absorption:

    switch(*groupRange) {

    case FastRes:
      fp=xsdbRetrieveNuclideItem(nuclide[j].name,
				 "ParticleBalance/Epithermal",&nbytes,verbose);
      for (i=0; i< nbytes/4; i++) 
	fpTarget[i]=fp[i];
      break;

    default:
      fprintf(logFile,"error:  unknown groupRange request\n");
      exit(1);
    }
    free(fp);
    break;

  case TransportCorrection:

    switch(*groupRange) {

    case FastRes:
      fp=xsdbRetrieveNuclideItem(nuclide[j].name,
				 "TransportCorrection/Epithermal",&nbytes,verbose);
      for (i=0; i< nbytes/4; i++) 
	fpTarget[i]=fp[i];
      break;

    default:
      fprintf(logFile,"error:  unknown groupRange request\n");
      exit(1);
    }
    free(fp);
    break;

  case Fission:
    switch(*groupRange) {

    case FastRes:
      fp=xsdbRetrieveNuclideItem(nuclide[j].name,
				 "Fission/Epithermal",&nbytes,verbose);
      for (i=0; i< nbytes/4; i++) 
	fpTarget[i]=fp[i];
      break;

    default:
      fprintf(logFile,"error:  unknown groupRange request\n");
      exit(1);
    }
    free(fp);
    break;

  case NuFission:
    switch(*groupRange) {

    case FastRes:
      fp=xsdbRetrieveNuclideItem(nuclide[j].name,
				 "NuFission/Epithermal",&nbytes,verbose);
      for (i=0; i< nbytes/4; i++) 
	fpTarget[i]=fp[i];
      break;

    default:
      fprintf(logFile,"error:  unknown groupRange request\n");
      exit(1);
    }
    free(fp);
    break;

  case N2n:

    switch(*groupRange) {

    case Fast:
      fp=xsdbRetrieveNuclideItem(nuclide[j].name,
				 "n2n/Fast",&nbytes,verbose);
      for (i=0; i< nbytes/4; i++) 
	fpTarget[i]=fp[i];
      break;

    default:
      fprintf(logFile,"error:  unknown groupRange request\n");
      exit(1);
    }
    free(fp);
    break;

  case FissionSpectrum:	/* read NGFISS chi values */

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "FissionSpectrum",&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);
    break;

  case PotScatSlowingDown:	

    /* potential scattering and slowing-down power back-to-back.  */
    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "PotScat/Res",&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "SDP/Res",&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      fpTarget[ng2+i]=fp[i];
    free(fp);
    break;

  case PotScat:

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "PotScat/Res",&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);
    break;

  case LengthsThermalP0:     /* includes only the thermal-scattering lengths */

    ip=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "LThermalP0",&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      ipTarget[i]=ip[i];
    free(ip);
    break;

  case LengthsScatP1:        /* includes the epithermal length, too, at the end */

    ip=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "LScatP1",&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      ipTarget[i]=ip[i];
    free(ip);
    break;

  case ScatP0:

    switch(*groupRange) {

    case FastRes:
      fp=xsdbRetrieveNuclideItem(nuclide[j].name,
				 "ScatP0/Epithermal",&nbytes,verbose);
      for (i=1; i< nbytes/4; i++) /* first element holds ndat; skip it */
	fpTarget[i-1]=fp[i];
      break;

    default:
      fprintf(logFile,"error:  unknown groupRange request\n");
      exit(1);
    }
    free(fp);
    break;

  case ScatP1:

    switch(*groupRange) {

    case FastRes:
      fp=xsdbRetrieveNuclideItem(nuclide[j].name,
				 "ScatP1/Epithermal",&nbytes,verbose);
      for (i=1; i< nbytes/4; i++) /* first element holds ndat; skip it */
	fpTarget[i-1]=fp[i];
      break;

    default:
      fprintf(logFile,"error:  unknown groupRange request\n");
      exit(1);
    }
    free(fp);
    break;

  case ThermalXSTemp:

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "ThermalXSTemp",&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);
    break;

  case ThermalP1Temp:

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "ThermalP1Temp",&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);
    break;

  case HeqHeader:

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "HeqHeader",&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);
    break;

  case Header:
    if ( !nuclide[j].installed) {
      ip=xsdbRetrieveNuclideItem(nuclide[j].name,"Header",&nbytes,verbose);
      fp=(float *)ip;
      nuclide[j].numericId=ip[0];
      nuclide[j].aw=fp[1];
      nuclide[j].iz=ip[2];
      nuclide[j].nf=ip[3];
      nuclide[j].nt=ip[4];
      nuclide[j].nr=ip[5];
      nuclide[j].ndat2=0; /* temporary -- ndat2, length of epithermal P0 scat rec */
      nuclide[j].ndat3=0; /* temporary -- ndat3, length of nt'th thermal P0 scat rec */
      nuclide[j].np1=ip[6];
      nuclide[j].ns=ip[7];
      nuclide[j].iendfb=ip[8];
      nuclide[j].installed=1;
      nLoadedNuclides++;
      free(ip);
    }
    ipTarget[0]=nuclide[j].nburn;	
    ipTarget[1]=nuclide[j].numericId;
    fpTarget[2]=nuclide[j].aw;
    ipTarget[3]=nuclide[j].iz;
    ipTarget[4]=nuclide[j].nf;
    ipTarget[5]=nuclide[j].nt;
    ipTarget[6]=nuclide[j].nr;
    ipTarget[7]=nuclide[j].ndat2;
    ipTarget[8]=nuclide[j].ndat3;
    ipTarget[9]=nuclide[j].np1;
    ipTarget[10]=nuclide[j].ns;
    ipTarget[11]=nuclide[j].iendfb;
    break;

  default:
    fprintf(logFile,"Serious error: unknown xsdb request %d!\n",*item);
    *status=UnknownRequest;
    return;
  }
  *status=NORMAL;
}

#ifdef NotNow

for (i=0; i<9; i++) {
  fprintf(logFile,"%4d: %d\n",i,((int *)where)[i]);
}
for (i=9; i<2510; i++) {
  fprintf(logFile,"%4d: %f\n",i,((float *)where)[i]);
}

#endif

void xsdthe_c(int *groupRange, DBItem *item, int *nuclideIndex, int *index, int *where, int *status)
{
  float *fp, *fpTarget;
  int nbytes;
  int i,j;
  float thisT;
  fpTarget=(float *)where;
  j=*nuclideIndex;

  if (j==LastNuclideAccessed)
    j=lastNuclideAccessed;
  else
    lastNuclideAccessed=j;
  j--;			/* offset now */

  if (verbose) {
    if (nuclide[j].installed)
      fprintf(logFile,"In xsdbloadnuclidetset with item %25s, nuclide %10s (%5d), index %5d, and location %016llX\n",
	      itemName(*item),nuclide[j].name,j+1,*index,(long long)where);
    else
      fprintf(logFile,"In xsdbloadnuclidetset with item %25s, nuclide %5d index %5d, and location %016llX\n",
	      itemName(*item),j+1,*index,(long long)where);
    fflush(logFile);
  }

  switch (*item) {

  case ScatP0:
    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "ThermalXSTemp",&nbytes,verbose);
    thisT=fp[*index-1];
    free(fp);

    fp=xsdbRetrieveNuclideTempItem(nuclide[j].name,
				   "ScatP0/Thermal",thisT,&nbytes,verbose);

    for (i=1; i< nbytes/4; i++) 
      fpTarget[i-1]=fp[i];
    free(fp);
    break;

  case ScatP1:

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "ThermalP1Temp",&nbytes,verbose);
    thisT=fp[*index-1];
    free(fp);

    if (verbose) {
      fprintf(logFile,"read thermal P1 -- %d, %e\n",*index,thisT);
      fflush(logFile);
    }

    fp=xsdbRetrieveNuclideTempItem(nuclide[j].name,
				   "ScatP1/Thermal",thisT,&nbytes,verbose);

    for (i=1; i< nbytes/4; i++) 
      fpTarget[i-1]=fp[i];
    free(fp);
    break;

  case Absorption:

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "ThermalXSTemp",&nbytes,verbose);
    thisT=fp[*index-1];
    free(fp);

    fp=xsdbRetrieveNuclideTempItem(nuclide[j].name,
				   "ParticleBalance/Thermal",thisT,&nbytes,verbose);

    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);
    break;

  case NuFission:

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "ThermalXSTemp",&nbytes,verbose);
    thisT=fp[*index-1];
    free(fp);

    fp=xsdbRetrieveNuclideTempItem(nuclide[j].name,
				   "NuFission/Thermal",thisT,&nbytes,verbose);

    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);
    break;

  case Fission:

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "ThermalXSTemp",&nbytes,verbose);
    thisT=fp[*index-1];
    free(fp);

    fp=xsdbRetrieveNuclideTempItem(nuclide[j].name,
				   "Fission/Thermal",thisT,&nbytes,verbose);

    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);
    break;

  case Transport:

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "ThermalXSTemp",&nbytes,verbose);
    thisT=fp[*index-1];
    free(fp);

    fp=xsdbRetrieveNuclideTempItem(nuclide[j].name,
				   "Transport/Thermal",thisT,&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);
    break;

  case TransportCorrection:

    fp=xsdbRetrieveNuclideItem(nuclide[j].name,
			       "ThermalXSTemp",&nbytes,verbose);
    thisT=fp[*index-1];
    free(fp);

    fp=xsdbRetrieveNuclideTempItem(nuclide[j].name,
				   "TransportCorrection/Thermal",thisT,&nbytes,verbose);
    for (i=0; i< nbytes/4; i++) 
      fpTarget[i]=fp[i];
    free(fp);
    break;

  default:

    fprintf(logFile,"Serious error: unknown xsdb request %d!\n",*item);
    *status=UnknownRequest;
    return;
  }
  *status=NORMAL;
}

void xsdres_c(int *nuclideIndex, int *where, int *status)
{
  int i, j;
  int nbytes;
  float *fp, *fpTarget;
  int *ip, *ipTarget;
  fpTarget=(float *)where+2;
  ipTarget=(int *)where;
  j=*nuclideIndex;

  if (j==LastNuclideAccessed)
    j=lastNuclideAccessed;
  else
    lastNuclideAccessed=j;
  j--;			/* offset now */

  if (verbose) 
    fprintf(logFile,"In XSDRES with nuclide %5d and location %016llX\n",j+1,(long long)where);

  ip=xsdbRetrieveNuclideItem(nuclide[j].name,"ResHeader",&nbytes,verbose);
  fp=(float *)ip+2;
  ipTarget[0]=ip[0]; /* m1 */
  ipTarget[1]=ip[1]; /* m2 */

  for (i=0; i<ip[0]+ip[1]; i++)
	fpTarget[i]=fp[i]; /* sigma-0 and T values */
  free(ip);

  *status=NORMAL;
  return;
}

void xsdtab_c(DBItem *item, int *nuclideIndex, int *resGroup, int *where, int *status)
{
  int i, j;
  int nbytes;
  float *fp, *fpTarget;
  fpTarget=(float *)where;
  j=*nuclideIndex;

  if (j==LastNuclideAccessed)
    j=lastNuclideAccessed;
  else
    lastNuclideAccessed=j;
  j--;			/* offset now */


  if (verbose)
    fprintf(logFile,"In XSDTAB with item %25s, nuclide %5d, resGroup %5d, and location %016llX\n",
	    itemName(*item),j+1,*resGroup,(long long)where);

  *status=UnknownNuclide;

  switch (*item) {

  case Absorption:
    fp=xsdbRetrieveNuclideGroupItem(nuclide[j].name,
				    "ResTable/ParticleBalance",*resGroup,
				    &nbytes,verbose);
    break;

  case NuFission:
    fp=xsdbRetrieveNuclideGroupItem(nuclide[j].name,
				    "ResTable/NuFission",*resGroup,
				    &nbytes,verbose);
    break;

  case ScatP0:
    fp=xsdbRetrieveNuclideGroupItem(nuclide[j].name,
				    "ResTable/ScatP0",*resGroup,
				    &nbytes,verbose);
    break;

  case NJOYFlux:
    fp=xsdbRetrieveNuclideGroupItem(nuclide[j].name,
				    "ResTable/NJOYFlux",*resGroup,
				    &nbytes,verbose);
    break;

  case Hequivalence:
    fp=xsdbRetrieveNuclideGroupItem(nuclide[j].name,
				    "Hequivalence",*resGroup,
				    &nbytes,verbose);
    break;

  default:
    fprintf(logFile,"Serious error: unknown xsdb request %d!\n",*item);
    *status=UnknownRequest;
    return;
  }
  for (i=0; i<nbytes/4; i++) {
    fpTarget[i]=fp[i];
  }
  free(fp);
  *status=NORMAL;
  return;
}

static char *nullName="NULL";

int nuclideId(char *name) {
  int i;
  if (!strcmp(name,nullName))
    return 0;
  for (i=0; i<nel; i++)
    if (!strcmp(nuclide[i].name,name))
      return nuclide[i].numericId;
  return -1;
}

int nuclideIndex(char *name) {
  int i;
  for (i=0; i<nel; i++)
    if (!strcmp(nuclide[i].name,name))
      return i;
  return -1;
}

static CodeNamePair btyp[]=
{
  {FISSION, "FISSION"},
  {FISSPROD,"FISSPROD"},
  {CAPTURE, "CAPTURE"},
  {DECAY,   "DECAY"},
  {N2N,     "N2N"},
};

int bcodeNumber(char *bname) {
  int i;
  for (i=0; i<5; i++)
    if (!strcmp(btyp[i].name,bname))
      return btyp[i].code;
  return -1;
}

/* unpackBD unravels fragments in the burnup-data structure */

void unpackBD(void *ptr, PackType type) {
  char *cptr=ptr;
  int i;
  switch (type) {
  case Char8:
    cptr[8]='\0';
    for (i=0; i<8; i++) {
      if (packedBurnupData[packIndex+i]==' ')
	cptr[i]='\0';
      else
	cptr[i]=packedBurnupData[packIndex+i];
    }
    packIndex+=8;
    break;
  case Int:
  case Float:
    for (i=0; i<4; i++) {
      cptr[i]=packedBurnupData[packIndex+i];
    }
    packIndex+=4;
    break;
  }
}
