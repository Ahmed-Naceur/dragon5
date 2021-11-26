*DECK LIBDI9
      SUBROUTINE LIBDI9 (MAXDIL,NGROUP,NAMFIL,HNISOR,HSHI,NDIL,DILUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the dilutions corresponding to a resonant isotope within a
* library in WIMS-E format.
*
*Copyright:
* Copyright (C) 2016 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXDIL  maximum number of dilutions.
* NGROUP  number of energy groups.
* NAMFIL  name of the WIMS-E format file.
* HNISOR  library name of the isotope.
* HSHI    library name of the self-shielding data.
*
*Parameters: output
* NDIL    number of finite dilutions.
* DILUT   dilutions.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER  NAMFIL*(*),HNISOR*12,HSHI*12
      INTEGER    MAXDIL,NGROUP,NDIL
      REAL       DILUT(MAXDIL)
*----
*  LOCAL VARIABLES
*----
      INTEGER    IUTYPE,LRIND,IACTO,IACTC
      PARAMETER (IUTYPE=2,LRIND=0,IACTO=2,IACTC=1,LPZ=8,MAXISO=246)
      CHARACTER  FMT*6,HSMG*131,CWISO(MAXISO)*8
      INTEGER    NPZ(LPZ),IWISO(2*MAXISO)
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR,SCR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GAR(MAXDIL+1))
*----
*  OPEN WIMS-E LIBRARY AND READ GENERAL DIMENSIONING
*----
      IDRES=INDEX(HSHI,'.')
      IUNIT=KDROPN(NAMFIL,IACTO,IUTYPE,LRIND)
      IF(IUNIT.LE.0) THEN
         WRITE (HSMG,'(35HLIBDI9: UNABLE TO OPEN LIBRARY FILE,1X,A16,
     1   8H. IUNIT=,I4,1H.)') NAMFIL,IUNIT
         CALL XABORT(HSMG)
      ENDIF
      READ(IUNIT) (NPZ(II),II=1,LPZ)
      IF(NPZ(2).NE.NGROUP) THEN
        CALL XABORT('LIBDI9: INVALID NUMBER OF GROUPS')
      ENDIF
      NEL=NPZ(1)
      NGR=NPZ(5)
      NW=4*NGR+2*NPZ(3)+NPZ(4)
      ALLOCATE(SCR(NW))
      CALL XDRSET(SCR,NW,0.0)
*----
*  READ ISOTOPE ID NUMBER AND CREATE EQUIVALENT ISOTOPE NAME
*  SCAN TO ASSOCIATE WIMS ISOTOPE NUMBER WITH DRAGON ISOTOPE NUMBER
*----
      IRISO=0
      IRIND=0
      READ(IUNIT) (IWISO(ITC),ITC=1,NEL)
      DO 10 IEL=1,NEL
        CWISO(IEL)='        '
        IF     (IWISO(IEL).LT.10) THEN
          WRITE(CWISO(IEL),'(I1)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.100) THEN
          WRITE(CWISO(IEL),'(I2)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.1000) THEN
          WRITE(CWISO(IEL),'(I3)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.10000) THEN
          WRITE(CWISO(IEL),'(I4)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.100000) THEN
          WRITE(CWISO(IEL),'(I5)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.1000000) THEN
          WRITE(CWISO(IEL),'(I6)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.10000000) THEN
          WRITE(CWISO(IEL),'(I7)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.100000000) THEN
          WRITE(CWISO(IEL),'(I8)') IWISO(IEL)
        ENDIF
        IF(CWISO(IEL).EQ.HNISOR(1:8)) THEN
          IRISO=IEL
          IF(IDRES.EQ.0) THEN
            IRIND=IWISO(IRISO)
          ENDIF
          GO TO 20
        ENDIF
   10 CONTINUE
      CALL XABORT('LIBDI9: ISOTOPE NOT FOUND ON LIBRARY')
   20 CONTINUE
*----
*  READ GROUP STRUCTURE
*----
      READ(IUNIT) (DUMMY,ITC=1,NGROUP)
*----
*  READ DEPLETION CHAIN
*----
      DO 30 IEL=1,NEL
        READ(IUNIT) DUMMY
   30 CONTINUE
*----
* READ FILE FOR TEMPERATURE DEPENDENT XS
*----
      NRTOT=0
      DO 50 IELRT=1,NEL
        READ(IUNIT) IDIEL,AWR,IZ,NFIEL,NTMP,NRIEL,ISOF,IP1OPT
        IF(NRIEL.GT.0) NRTOT=NRTOT+NRIEL
        IF(IELRT.EQ.IRISO) THEN
          READ(IUNIT) (SCR(I),I=1,NW)
        ELSE
          READ(IUNIT) DUMMY
        ENDIF
        IF(NFIEL.GT.1) READ(IUNIT) DUMMY
        READ(IUNIT) NSCT
        IF(NTMP.GT.0) THEN
          READ(IUNIT) DUMMY
          DO 40 IT=1,NTMP
            READ(IUNIT) DUMMY
            IF(NFIEL.GT.1) READ(IUNIT) DUMMY
            READ(IUNIT) NSCT
   40     CONTINUE
        ENDIF
        IF(ISOF.NE.0) READ(IUNIT) DUMMY
        IF(IP1OPT.NE.1) THEN
          DO 45 IT=1,NTMP
            READ(IUNIT) DUMMY
   45     CONTINUE
        ENDIF
   50 CONTINUE
*----
*  IDENTIFY RESONANT ISOTOPE
*----
        IF(IDRES.GT.0) THEN
          WRITE(FMT,'(2H(F,I1,3H.1))') IDRES+1
          READ(HSHI,FMT) RIND
        ELSE
          RIND=FLOAT(IWISO(IRISO))
        ENDIF
*----
*  SCAN OVER RESONANCE SETS+1 AND READ RESONANCE INFO
*----
      DO 90 IGR=1,NGR
        DO 70 IRS=1,3*NRTOT+1
          READ(IUNIT) XIDR,M1,M2,(DUMMY,IT=1,M1),(GAR(ID),ID=1,M2),
     1    ((DUMMY,ID=1,M2),IT=1,M1)
         IF(M2.GT.MAXDIL) CALL XABORT('LIBDI9: MAXDIL OVERFLOW.')
          IF(IDRES.EQ.0) THEN
            XRS1=FLOAT(INT((XIDR+0.01)*10.)-INT(XIDR+0.01)*10)/10.
            XRS1=ABS(XIDR-XRS1-RIND)
          ELSE
            XRS1=ABS(XIDR-RIND)
          ENDIF
          IF(XRS1.LE.0.01) THEN
            NDIL=M2-1
            DO 60 ID=1,NDIL
            DSIGPL=SCR(IGR)*SCR(NW-NGR+IGR)
            IF(GAR(ID)-DSIGPL.GT.0.0) THEN
              DILUT(ID)=GAR(ID)-DSIGPL
            ELSE
              DILUT(ID)=0.0
            ENDIF
   60       CONTINUE
            DILUT(M2)=MIN(GAR(M2),1.0E10)
            GO TO 100
          ENDIF
          IF(XIDR.EQ.0.0) GO TO 80
   70   CONTINUE
   80   CONTINUE
        READ(IUNIT) ENDR
   90 CONTINUE
      NDIL=0
      DILUT(1)=1.0E10
  100 IERR=KDRCLS(IUNIT,IACTC)
      IF(IERR.LT.0) THEN
         HNISOR=NAMFIL
         CALL XABORT('LIBDI9: WIMS-E LIBRARY '//HNISOR//
     1   ' CANNOT BE CLOSED.')
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SCR,GAR)
      RETURN
      END
