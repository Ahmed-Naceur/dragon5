*DECK LZC
      SUBROUTINE LZC(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read specification for the liquid zone controllers; add the new data
* to the existing device object.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* HENTRY  name of the data structures.
* IENTRY  data structure type where:
*         IENTRY=1 for LCM memory object;
*         IENTRY=2 for XSM file;
*         IENTRY=3 for sequential binary file;
*         IENTRY=4 for sequential ASCII file.
* JENTRY  access permission for the data structure where:
*         JENTRY=0 for a data structure in creation mode;
*         JENTRY=1 for a data structure in modifications mode;
*         JENTRY=2 for a data structure in read-only mode.
* KENTRY  data structure pointer.
*
*Comments:
* The LZC: module specification is:
* DEVICE MATEX := LZC: [ DEVICE ] MATEX :: (desclzc) ;
* where
*   DEVICE : name of the \emph{device} object.
* Note, if the rod-type devices are not present in the reactor core, then 
* DEVICE object must appear only on the LHS (i.e. in create mode), it will 
* contain the information only with respect to the liquid zone controllers. 
* However, if the rod-type devices are present in the reactor core, then they 
* must be specified first (i.e. before the liquid controllers) using the DEVINI:
* module. In the last case, the DEVICE object must also appear on the RHS 
* (i.e. in modification mode), it will contain the additional and separate 
* information with respect to the liquid zone controllers.
*   MATEX  : name of the \emph{matex} object 
* that will be updated by the module. The lzc-devices material mixtures are
* appended to the previous material index and the lzc-devices indices are
* also modified, accordingly.
*   (desclzc) : structure describing the input data to the LZC: module.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      CHARACTER HSIGN*12,TEXT12*12
      INTEGER ISTATE(NSTATE)
      REAL LIMIT(6)
      LOGICAL LNEW
      TYPE(C_PTR) IPDEV,IPMTX
      REAL, ALLOCATABLE, DIMENSION(:) :: XXX,YYY,ZZZ
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.NE.2)CALL XABORT('@LZC: TWO PARAMETERS EXPECTED')
      TEXT12=HENTRY(1)
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@LZ'
     1 //'C: LCM OBJECT EXPECTED FOR L_DEVICE ('//TEXT12//').')
      IF(JENTRY(1).EQ.1)THEN
        CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_DEVICE')CALL XABORT('@LZC: MISSING L_DEV'
     1   //'ICE OBJECT.')
        LNEW=.FALSE.
      ELSEIF(JENTRY(1).EQ.0)THEN
        HSIGN='L_DEVICE'
        CALL LCMPTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
        LNEW=.TRUE.
      ELSE
       CALL XABORT('@LZC: ONLY CREATE OR MODIFICATION MODE EXPEC'
     1 //'TED FOR L_DEVICE OBJECT.')
      ENDIF
      IPDEV=KENTRY(1)
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2))CALL XABORT('@LZ'
     1 //'C: LCM OBJECT EXPECTED FOR L_MATEX.')
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MATEX')CALL XABORT('@LZC: MISSING L_MATEX.')
      IF(JENTRY(2).NE.1)CALL XABORT('@LZC: MODIFICATION MODE EX'
     1 //'PECTED FOR L_MATEX.')
      IPMTX=KENTRY(2)
*----
*  RECOVER INFORMATION
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMTX,'STATE-VECTOR',ISTATE)
      IGEO=ISTATE(6)
      IF(IGEO.NE.7)CALL XABORT('@LZC: ONLY'
     1 //' 3D-CARTESIAN GEOMETRY ALLOWED.')
      NMIX=ISTATE(2)
      NTOT=ISTATE(5)
      LX=ISTATE(8)
      LY=ISTATE(9)
      LZ=ISTATE(10)
*     LIMITS ALONG X-AXIS
      ALLOCATE(XXX(LX+1))
      CALL XDRSET(XXX,LX+1,0.)
      CALL LCMGET(IPMTX,'MESHX',XXX)
      LIMIT(1)=XXX(1)
      LIMIT(2)=XXX(LX+1)
      DEALLOCATE(XXX)
*     LIMITS ALONG Y-AXIS
      ALLOCATE(YYY(LY+1))
      CALL XDRSET(YYY,LY+1,0.)
      CALL LCMGET(IPMTX,'MESHY',YYY)
      LIMIT(3)=YYY(1)
      LIMIT(4)=YYY(LY+1)
      DEALLOCATE(YYY)
*     LIMITS ALONG Z-AXIS
      ALLOCATE(ZZZ(LZ+1))
      CALL XDRSET(ZZZ,LZ+1,0.)
      CALL LCMGET(IPMTX,'MESHZ',ZZZ)
      LIMIT(5)=ZZZ(1)
      LIMIT(6)=ZZZ(LZ+1)
      DEALLOCATE(ZZZ)
*     READ LZC INPUT DATA
      CALL LZCDRV(IPDEV,IPMTX,IGEO,NMIX,NTOT,LIMIT,LNEW)
      RETURN
      END
