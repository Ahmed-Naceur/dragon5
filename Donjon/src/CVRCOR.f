*DECK CVRCOR
      SUBROUTINE CVRCOR(IPMAP,NCH,NB,NFUEL,NX,NY,NZ,IVOID,NVOID,NPARM,
     1 PNAME,PVALUE,VCOOL,LCOOL,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Modify channels data according to the specified core-voiding pattern.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input/output
* IPMAP   pointer to the perturbed fuel-map.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* NFUEL   number of fuel types.
* NX      number of elements along x-axis in fuel map.
* NY      number of elements along y-axis in fuel map.
* NZ      number of elements along z-axis in fuel map.
* IVOID   index associated with the core-voiding pattern:
*          =1 full-core; =2 half-core; =3 quarter-core;
*          =4 checkerboard-full; =5 checkerboard-half;
*          =6 checkerboard-quarter.
* NVOID   total number of voided channels.
* NPARM   total number of recorded parameters.
* PNAME   recorded parameter name for the coolant density.
* PVALUE  structure containing the coolant density values
*         throughout the reactor core.
* VCOOL   coolant density value for voided channels.
* LCOOL   flag with respect to the coolant densities:
*          =.true. to modify these values;
*          =.false. coolant densities not provided.
* IMPX    printing index (=0 for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NCH,NB,NFUEL,NX,NY,NZ,IVOID,NVOID,NPARM,IMPX
      REAL PVALUE(NCH,NB),VCOOL
      CHARACTER PNAME*12
      LOGICAL LCOOL
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER FLMIX(NCH,NB),NSCH(NCH),INAME(3)
      CHARACTER TEXT*20,TEXT12*12
      LOGICAL LCHK
      TYPE(C_PTR) JPMAP,KPMAP
*----
*  RECOVER INFORMATION
*----
      LCHK=.FALSE.
      CALL XDISET(FLMIX,NCH*NB,0)
      CALL LCMGET(IPMAP,'FLMIX',FLMIX)
      CALL LCMLEN(IPMAP,'REF-SCHEME',LENGT,ITYP)
      IF(LENGT.EQ.0)CALL XABORT('@CVRCOR: MISSI'
     1 //'NG REF-SCHEME DATA IN FUEL-MAP.')
      CALL XDISET(NSCH,NCH,0)
      CALL LCMGET(IPMAP,'REF-SCHEME',NSCH)
      IF(IVOID.EQ.1)THEN
        TEXT='FULL-CORE'
      ELSEIF(IVOID.EQ.2)THEN
        TEXT='HALF-CORE'
      ELSEIF(IVOID.EQ.3)THEN
        TEXT='QUARTER-CORE'
      ELSEIF(IVOID.EQ.4)THEN
        TEXT='CHECKERBOARD-FULL'
        LCHK=.TRUE.
      ELSEIF(IVOID.EQ.5)THEN
        TEXT='CHECKERBOARD-HALF'
        LCHK=.TRUE.
      ELSEIF(IVOID.EQ.6)THEN
        TEXT='CHECKERBOARD-QUARTER'
        LCHK=.TRUE.
      ENDIF
      IF(IMPX.GT.0)WRITE(IOUT,1000)TEXT,NVOID
*----
*  MODIFY CHANNEL DATA
*----
      ITOT=0
      JPMAP=LCMGID(IPMAP,'FUEL')
      DO IFUEL=1,NFUEL
        KPMAP=LCMGIL(JPMAP,IFUEL)
        CALL LCMGET(KPMAP,'MIX',MIXF)
        CALL LCMGET(KPMAP,'MIX-VOID',MIXV)
        DO 20 ICH=1,NVOID
        IF(LCHK)THEN
          IF(NSCH(ICH).LT.0)GOTO 20
*         POSITIVE DIRECTION ONLY
        ENDIF
        DO 10 IB=1,NB
        IF(FLMIX(ICH,IB).NE.MIXF)GOTO 10
        FLMIX(ICH,IB)=MIXV
        IF(LCOOL) PVALUE(ICH,IB)=VCOOL
        ITOT=ITOT+1
   10   CONTINUE
   20   CONTINUE
      ENDDO
      IF(IMPX.GT.0)WRITE(IOUT,1001)ITOT
      IF(IMPX.LT.2)GOTO 30
*     PRINTING
      CALL CVRPRN(IPMAP,NCH,NB,NX,NY,NZ,FLMIX,PVALUE,LCOOL,IMPX)
*     STORE NEW DATA
   30 CALL LCMPUT(IPMAP,'FLMIX',NCH*NB,1,FLMIX)
      IF(.NOT.LCOOL)GOTO 40
      JPMAP=LCMGID(IPMAP,'PARAM')
      DO IPAR=1,NPARM
        KPMAP=LCMGIL(JPMAP,IPAR)
        CALL LCMGET(KPMAP,'P-NAME',INAME)
        WRITE(TEXT12,'(3A4)') (INAME(I),I=1,3)
        IF(PNAME.EQ.TEXT12)THEN
          CALL LCMPUT(KPMAP,'P-VALUE',NCH*NB,2,PVALUE)
          GOTO 40
        ENDIF
      ENDDO
   40 RETURN
*
 1000 FORMAT(/2X,'SELECTED VOIDING PATTERN',2X,'=>',2X,A20
     1     //2X,'TOTAL NUMBER OF VOIDED CHANNELS =',1X,I3/)
 1001 FORMAT(2X,'TOTAL NUMBER OF MODIFIED VALUES :',1X,I4/)
      END
