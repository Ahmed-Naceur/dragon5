*DECK CVRPRN
      SUBROUTINE CVRPRN(IPMAP,NCH,NB,NX,NY,NZ,MIXNEW,PVALUE,LCOOL,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print the fuel-type indices per bundle for each reactor channel.
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
* NX      number of elements along x-axis in fuel map.
* NY      number of elements along y-axis in fuel map.
* NZ      number of elements along z-axis in fuel map.
* MIXNEW  updated fuel-type index.
* PVALUE  structure containing the modified coolant density
*         values throughout the reactor core.
* LCOOL   flag with respect to the coolant densities:
*          =.true. coolant densities were modified;
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
      INTEGER NCH,NB,NX,NY,NZ,MIXNEW(NCH,NB),IMPX
      REAL PVALUE(NCH,NB)
      LOGICAL LCOOL
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER MIX(NX*NY*NZ),MIXOLD(NCH,NB),
     1 NSCH(NCH),NAMX(NX),NAMY(NY),FMIX(NX,NY,NB)
      CHARACTER TEXT*12,CHANX*2,CHANY*2,FORM1*14,FORM2*14
*----
*  RECOVER INFORMATION
*----
      CALL XDISET(MIX,NX*NY*NZ,0)
      CALL LCMGET(IPMAP,'BMIX',MIX)
      CALL XDISET(MIXOLD,NCH*NB,0)
      CALL LCMGET(IPMAP,'FLMIX',MIXOLD)
*     CHANNEL NAMES
      CALL XDISET(NAMX,NX,0)
      CALL LCMGET(IPMAP,'XNAME',NAMX)
      CALL XDISET(NAMY,NY,0)
      CALL LCMGET(IPMAP,'YNAME',NAMY)
*     REFUELLING SCHEME
      CALL LCMLEN(IPMAP,'REF-SCHEME',LENGT,ITYP)
      IF(LENGT.EQ.0)CALL XABORT('@CVRPRN: MISSI'
     1 //'NG REF-SCHEME DATA IN FUEL-MAP.')
      CALL XDISET(NSCH,NCH,0)
      CALL LCMGET(IPMAP,'REF-SCHEME',NSCH)
      IF((IMPX.NE.2).AND.(IMPX.LT.4))GOTO 20
*----
*  PRINTING OVER EACH CHANNEL
*----
      IEL=0
      ICH=0
      DO 15 J=1,NY
      DO 10 I=1,NX
      IEL=IEL+1
      IF(MIX(IEL).EQ.0)GOTO 10
      ICH=ICH+1
      WRITE(TEXT,'(A9,I3.3)')'CHANNEL #',ICH
      WRITE(CHANX,'(A2)') (NAMX(I))
      WRITE(CHANY,'(A2)') (NAMY(J))
      WRITE(IOUT,1001)TEXT,CHANY,CHANX,NSCH(ICH)
      WRITE(IOUT,1002) (MIXOLD(ICH,IB),IB=1,NB)
      WRITE(IOUT,1003) (MIXNEW(ICH,IB),IB=1,NB)
      IF(LCOOL)WRITE(IOUT,1004)(PVALUE(ICH,IB),IB=1,NB)
   10 CONTINUE
   15 CONTINUE
   20 IF((IMPX.NE.3).AND.(IMPX.LT.4))GOTO 50
*----
*  PRINTING PER RADIAL PLANE
*----
      WRITE(FORM1,'(A4,I2,A8)')'(A4,',NX,'(A3,1X))'
      WRITE(FORM2,'(A4,I2,A8)')'(A2,',NX,'(I3,1X))'
      CALL XDISET(FMIX,NX*NY*NB,0)
      WRITE(IOUT,1005)
      IEL=0
      DO IB=1,NB
        ICH=0
        DO 35 J=1,NY
        DO 30 I=1,NX
        IEL=IEL+1
        IF(MIX(IEL).EQ.0)GOTO 30
        ICH=ICH+1
        FMIX(I,J,IB)=MIXNEW(ICH,IB)
   30   CONTINUE
   35   CONTINUE
      ENDDO
      DO IB=1,NB
        WRITE(IOUT,1006)IB
        WRITE(IOUT,FORM1)' ',(NAMX(I),I=1,NX)
        WRITE(IOUT,*)' '
        DO 40 J=1,NY
        WRITE(CHANY,'(A2)') (NAMY(J))
        IF(INDEX(CHANY,'-').EQ.1)GOTO 40
        WRITE(IOUT,FORM2)CHANY,(FMIX(I,J,IB),I=1,NX)
   40   CONTINUE
      ENDDO
   50 RETURN
*
 1001 FORMAT(/10X,'* ',A12,' *',5X,'NAME:',1X,A2,A2,
     1 5X,'REF-SCHEME:',1X,I2)
 1002 FORMAT(2X,'OLD FUEL INDICES :',2X,12(I2,2X))
 1003 FORMAT(2X,'NEW FUEL INDICES :',2X,12(I2,2X)/)
 1004 FORMAT(2X,'COOLANT DENSITIES:',1X,12(F4.2,1X)/)
 1005 FORMAT(//20X,'** FUEL-TYPE INDICES PER RADIAL PLANE **')
 1006 FORMAT(//2X,'RADIAL PLANE',1X,'#',I2.2/)
      END
