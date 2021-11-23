*DECK CVRUSR
      SUBROUTINE CVRUSR(IPMAP,NCH,NB,NFUEL,NX,NY,NZ,NVOID,NAMXV,NAMYV,
     1 NPARM,PNAME,PVALUE,VCOOL,LCOOL,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Modify channels data according to the user-defined voiding pattern.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input/output
* IPMAP   pointer to fuel-map information.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* NFUEL   number of fuel types.
* NX      number of elements along x-axis in fuel map.
* NY      number of elements along y-axis in fuel map.
* NZ      number of elements along z-axis in fuel map.
* NVOID   total number of voided channels.
* NAMXV   names of voided channels along x-axis.
* NAMYV   names of voided channels along y-axis.
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
      INTEGER NCH,NB,NFUEL,NX,NY,NZ,NVOID,NPARM,IMPX,NAMXV(NVOID),
     1 NAMYV(NVOID)
      REAL PVALUE(NCH,NB),VCOOL
      CHARACTER PNAME*12
      LOGICAL LCOOL
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER MIX(NX*NY*NZ),FLMIX(NCH,NB),NAMX(NX),NAMY(NY),INAME(3)
      CHARACTER TEXT*12,CHANX*2,CHANY*2
      TYPE(C_PTR) JPMAP,KPMAP
      INTEGER, ALLOCATABLE, DIMENSION(:) :: CNANV
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(CNANV(NVOID))
*----
*  RECOVER INFORMATION
*----
      CALL XDISET(MIX,NX*NY*NZ,0)
      CALL LCMGET(IPMAP,'BMIX',MIX)
      CALL XDISET(FLMIX,NCH*NB,0)
      CALL LCMGET(IPMAP,'FLMIX',FLMIX)
*     CHANNEL NAMES
      CALL XDISET(NAMX,NX,0)
      CALL LCMGET(IPMAP,'XNAME',NAMX)
      CALL XDISET(NAMY,NY,0)
      CALL LCMGET(IPMAP,'YNAME',NAMY)
      TEXT='USER-DEFINED'
      IF(IMPX.GT.0)WRITE(IOUT,1000)TEXT,NVOID
*----
*  CHECK VOIDED CHANNELS
*----
      DO 20 IVD=1,NVOID
      IEL=0
      ICH=0
      DO 15 J=1,NY
      DO 10 I=1,NX
        IEL=IEL+1
        IF(MIX(IEL).EQ.0)GOTO 10
        ICH=ICH+1
        IF(NAMXV(IVD).NE.NAMX(I))GOTO 10
        IF(NAMYV(IVD).NE.NAMY(J))GOTO 10
        CNANV(IVD)=ICH
        GOTO 20
   10 CONTINUE
   15 CONTINUE
      WRITE(CHANX,'(A2)') (NAMXV(IVD))
      WRITE(CHANY,'(A2)') (NAMYV(IVD))
      WRITE(IOUT,1001)CHANY,CHANX
      CALL XABORT('@CVRUSR: INVALID INPUT DATA.')
   20 CONTINUE
*----
*  MODIFY CHANNEL DATA
*----
      ITOT=0
      JPMAP=LCMGID(IPMAP,'FUEL')
      DO IFUEL=1,NFUEL
        KPMAP=LCMGIL(JPMAP,IFUEL)
        CALL LCMGET(KPMAP,'MIX',MIXF)
        CALL LCMGET(KPMAP,'MIX-VOID',MIXV)
        DO IVD=1,NVOID
          ICH=CNANV(IVD)
          DO 30 IB=1,NB
          IF(FLMIX(ICH,IB).NE.MIXF)GOTO 30
          FLMIX(ICH,IB)=MIXV
          IF(LCOOL) PVALUE(ICH,IB)=VCOOL
          ITOT=ITOT+1
   30     CONTINUE
        ENDDO
      ENDDO
      IF(IMPX.GT.0)WRITE(IOUT,1002)ITOT
      IF(IMPX.LT.2)GOTO 40
*     PRINTING
      CALL CVRPRN(IPMAP,NCH,NB,NX,NY,NZ,FLMIX,PVALUE,LCOOL,IMPX)
*     STORE NEW DATA
   40 CALL LCMPUT(IPMAP,'FLMIX',NCH*NB,1,FLMIX)
      IF(.NOT.LCOOL)GOTO 50
      JPMAP=LCMGID(IPMAP,'PARAM')
      DO IPAR=1,NPARM
        KPMAP=LCMGIL(JPMAP,IPAR)
        CALL LCMGET(KPMAP,'P-NAME',INAME)
        WRITE(TEXT,'(3A4)') (INAME(I),I=1,3)
        IF(PNAME.EQ.TEXT)THEN
          CALL LCMPUT(KPMAP,'P-VALUE',NCH*NB,2,PVALUE)
          GOTO 50
        ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
   50 DEALLOCATE(CNANV)
      RETURN
*
 1000 FORMAT(/2X,'SELECTED VOIDING PATTERN',2X,'=>',2X,A20
     1     //2X,'TOTAL NUMBER OF VOIDED CHANNELS =',1X,I3/)
 1001 FORMAT(/1X,'@CVRUSR: UNABLE TO FIND THE CHANN',
     1 'EL NAME:',1X,A2,A2)
 1002 FORMAT(2X,'TOTAL NUMBER OF MODIFIED VALUES :',1X,I4/)
      END
