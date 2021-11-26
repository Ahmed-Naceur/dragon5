*DECK RESBRN
      SUBROUTINE RESBRN(IPMAP,NCH,NB,NCOMB,NX,NY,NZ,LRSCH,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Initialize the axial shape and compute the first burnup limits per
* bundle for every channel (used for the time-average model).
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
*
*Author(s): 
*  D. Sekki
*
*Parameters: input
* IPMAP   pointer to fuel-map information.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* NCOMB   number of combustion zones.
* NX      number of elements along x-axis in fuel map.
* NY      number of elements along y-axis in fuel map.
* NZ      number of elements along z-axis in fuel map.
* LRSCH   flag for the refuelling scheme of channels:
*          =.true. it was read from the input file;
*          =.false. otherwise.
* IMPX    printing index (=0 for no print).
*
*----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NCH,NB,NCOMB,NX,NY,NZ,IMPX
      LOGICAL LRSCH
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER IVECT(NCOMB,NB),NSCH(NCH),IZONE(NCH),MIX(NX*NY*NZ),
     1 NAMX(NX),NAMY(NY),RSCH(NX,NY)
      REAL BVAL(NCOMB),DELT(NB),B0(NB),B1(NB),SHAP(NCH,NB),
     1 BURN0(NCH,NB),BURN1(NCH,NB)
      CHARACTER TEXT*12,CHANY*2,FORM1*14,FORM2*14
*----
*  RECOVER INFORMATION
*----
      CALL LCMLEN(IPMAP,'REF-SCHEME',LENG1,ITYP)
      CALL LCMLEN(IPMAP,'BURN-AVG',LENG2,ITYP)
      IF((LENG1.EQ.0).OR.(LENG2.EQ.0))GOTO 100
      CALL LCMLEN(IPMAP,'AX-SHAPE',LENG3,ITYP)
      IF(LENG3.EQ.0) THEN
*       INITIAL FLAT AXIAL-SHAPE
        IF(IMPX.GT.0)WRITE(IOUT,1000)
        CALL XDRSET(SHAP,NCH*NB,1./NB)
        CALL LCMPUT(IPMAP,'AX-SHAPE',NCH*NB,2,SHAP)
      ELSE
        CALL LCMGET(IPMAP,'AX-SHAPE',SHAP)
      ENDIF
      CALL LCMGET(IPMAP,'REF-VECTOR',IVECT)
      CALL LCMGET(IPMAP,'REF-SCHEME',NSCH)
      CALL LCMGET(IPMAP,'BURN-AVG',BVAL)
      CALL LCMGET(IPMAP,'B-ZONE',IZONE)
      CALL LCMGET(IPMAP,'BMIX',MIX)
      CALL XDRSET(BURN0,NCH*NB,0.)
      CALL XDRSET(BURN1,NCH*NB,0.)
      IF(IMPX.GT.2)WRITE(IOUT,1004)
*----
*  COMPUTE FIRST BURNUP LIMITS
*----
      ICH=0
      DO 60 IEL=1,NX*NY
      IF(MIX(IEL).EQ.0)GOTO 60
      ICH=ICH+1
      IBSH=ABS(NSCH(ICH))
      DO IB=1,NB
        DELT(IB)=IBSH*BVAL(IZONE(ICH))*SHAP(ICH,IB)
        B0(IB)=0.
        B1(IB)=0.
      ENDDO
      IF(NSCH(ICH).LT.0)THEN
        DO 40 IB=1,NB
        KK=NB-IB+1
        KA=NB-IVECT(IZONE(ICH),IB)+1
        IF(IVECT(IZONE(ICH),IB).LE.0)THEN
          B0(KK)=0.
        ELSE
          B0(KK)=B1(KA)
        ENDIF
        B1(KK)=B0(KK)+DELT(KK)
   40   CONTINUE
      ELSE
        DO 50 IB=1,NB
        IF(IVECT(IZONE(ICH),IB).LE.0)THEN
          B0(IB)=0.
        ELSE
          B0(IB)=B1(IVECT(IZONE(ICH),IB))
        ENDIF
        B1(IB)=B0(IB)+DELT(IB)
   50   CONTINUE
      ENDIF
      DO IB=1,NB
        BURN0(ICH,IB)=B0(IB)
        BURN1(ICH,IB)=B1(IB)
      ENDDO
      IF(IMPX.GE.3) THEN
*       CHECK BURNUP LIMITS
        WRITE(TEXT,'(A9,I3.3)')'CHANNEL #',ICH
        WRITE(IOUT,1001)TEXT,NSCH(ICH),IZONE(ICH)
        WRITE(IOUT,1002)'B0',(B0(IB),IB=1,NB)
        WRITE(IOUT,1002)'B1',(B1(IB),IB=1,NB)
      ENDIF
   60 CONTINUE
      CALL LCMPUT(IPMAP,'BURN-BEG',NB*NCH,2,BURN0)
      CALL LCMPUT(IPMAP,'BURN-END',NB*NCH,2,BURN1)
      IF((.NOT.LRSCH).OR.(IMPX.LT.2))GOTO 100
*----
*  PRINT CHANNELS REFUELLING SCHEMES
*----
      WRITE(FORM1,'(A4,I2,A8)')'(A4,',NX,'(A3,1X))'
      WRITE(FORM2,'(A4,I2,A8)')'(A2,',NX,'(I3,1X))'
      CALL LCMGET(IPMAP,'XNAME',NAMX)
      CALL LCMGET(IPMAP,'YNAME',NAMY)
      CALL XDISET(RSCH,NX*NY,0)
      WRITE(IOUT,1003)
      IEL=0
      ICH=0
      DO 75 J=1,NY
      DO 70 I=1,NX
      IEL=IEL+1
      IF(MIX(IEL).EQ.0)GOTO 70
      ICH=ICH+1
      RSCH(I,J)=NSCH(ICH)
   70 CONTINUE
   75 CONTINUE
      WRITE(IOUT,FORM1)' ',(NAMX(I),I=1,NX)
      WRITE(IOUT,*)' '
      DO 80 J=1,NY
      WRITE(CHANY,'(A2)') (NAMY(J))
      IF(INDEX(CHANY,'-').EQ.1)GOTO 80
      WRITE(IOUT,FORM2)CHANY,(RSCH(I,J),I=1,NX)
   80 CONTINUE
  100 RETURN
*
 1000 FORMAT(/1X,'INITIALIZING THE FLAT AXIAL POWER-SHAPE'/
     1 1X,'COMPUTING THE FIRST BURNUP LIMITS PER EACH CHANNEL'/)
 1001 FORMAT(/10X,
     1 A12,10X,'REFUELLING SCHEME:',I3,10X,'ZONE-INDEX:',I3)
 1002 FORMAT(A3,12(F8.1,1X))
 1003 FORMAT(//20X,'** CHANNELS REFUELLING SCHEMES **'/)
 1004 FORMAT(/20X,'** FIRST BURNUP LIMITS PER EACH CHANNEL **'/)
      END
