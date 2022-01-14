*DECK READMT
      SUBROUTINE READMT (MAXPTS,IPGEOM,IR,MAT,VOL,ILK,ISTAT,NSUPCE,IPAS,
     1 NMC,RAYRE,PROCEL,POURCE,SURFA,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the input data for the description of an arbitrary multicell
* using the 'do-it-yourself' approach.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* MAXPTS  allocated storage for arrays of dimension NSUPCE or IPAS.
* IPGEOM  pointer to the geometry L_GEOM signature).
* IR      number of mixtures.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* ILK     leakage flag (ILK=.true. if neutron leakage through external
*         boundary is present).
* ISTAT   statistical approximation flag (set with ISTAT=1).
* NSUPCE  number of cells.
* IPAS    number of volumes.
* NMC     offset of the first volume in each cell.
* RAYRE   radius of the tubes in each cell.
* PROCEL  user supplied geometrical matrix.
* POURCE  weight assign to each cell.
* SURFA   surface assign to each cell.
* IMPX    print flag (equal to 0 for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM
      INTEGER MAXPTS,IR,ISTAT,NSUPCE,IPAS,MAT(IPAS),NMC(NSUPCE+1),IMPX
      REAL VOL(IPAS),RAYRE(NSUPCE+IPAS),PROCEL(NSUPCE**2),
     1 POURCE(NSUPCE),SURFA(NSUPCE)
      LOGICAL ILK
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PI=3.141592654,EPS1=1.0E-5,EPS2=1.0E-4,NSTATE=40)
      LOGICAL EMPTY,LCM
      CHARACTER GEONAM*12,TEXT12*12,HSMG*131
      DOUBLE PRECISION RRINT,RRZON,PAS,RJ,RJ1,RJN,RJN1
      INTEGER ISTATE(NSTATE),ISTAT2(NSTATE)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBREG,MILIEU,CELL,IGEN
      REAL, ALLOCATABLE, DIMENSION(:) :: RAYZON
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NBREG(MAXPTS),MILIEU(MAXPTS),CELL(3*MAXPTS),IGEN(MAXPTS))
      ALLOCATE(RAYZON(MAXPTS))
*
      CALL LCMINF(IPGEOM,GEONAM,TEXT12,EMPTY,ILONG,LCM)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      NSUPCE=ISTATE(3)
      IF(NSUPCE.GT.MAXPTS) CALL XABORT('READMT: INSUFFICIENT STORAGE.')
      CALL LCMLEN(IPGEOM,'PROCEL',ILEN,ITYLCM)
      IF(ILEN.EQ.0) THEN
         ISTAT=1
      ELSE
         ISTAT=0
         CALL LCMGET(IPGEOM,'PROCEL',PROCEL)
      ENDIF
      CALL LCMGET(IPGEOM,'POURCE',POURCE)
*
      SURFAT=0.0
      IPAS=0
      NMC(1)=0
      PTOT=1.0
      CALL LCMLEN(IPGEOM,'MIX',ILON,ITYLCM)
      IF(ILON.NE.NSUPCE) CALL XABORT('READMT: INVALID NUMBER OF SUB'
     1 //' GEOMETRIES.')
      CALL LCMGET(IPGEOM,'CELL',CELL)
      CALL LCMGET(IPGEOM,'MIX',IGEN)
      DO 50 IKK=1,NSUPCE
      PTOT=PTOT-POURCE(IKK)
      IKG=-IGEN(IKK)
      IF(IKG.LT.0) CALL XABORT('READMT: GENERATING CELLS EXPECTED.')
      WRITE (TEXT12(:4),'(A4)') CELL(3*(IKG-1)+1)
      WRITE (TEXT12(5:8),'(A4)') CELL(3*(IKG-1)+2)
      WRITE (TEXT12(9:),'(A4)') CELL(3*(IKG-1)+3)
      CALL LCMLEN(IPGEOM,TEXT12,ILEN,ITYLCM)
      IF((ILEN.EQ.0).OR.(ITYLCM.NE.0)) CALL XABORT('READMT: SUB GEOME'
     1 //'TRY '//TEXT12//' IS MISSING FROM GEOMETRY '//GEONAM//'.')
      CALL LCMSIX(IPGEOM,TEXT12,1)
         CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTAT2)
         IF(ISTAT2(1).NE.3) CALL XABORT('READMT: A SUB GEOMETRY IS N'
     1   //'OT A TUBE.')
         NZONE=ISTAT2(2)
         CALL LCMLEN(IPGEOM,'SPLITR',ILEN,ITYLCM)
         IF(ILEN.GT.0) THEN
            CALL LCMGET(IPGEOM,'SPLITR',NBREG)
         ELSE
            DO 10 I=1,NZONE
            NBREG(I)=1
   10       CONTINUE
         ENDIF
         CALL LCMGET(IPGEOM,'MIX',MILIEU)
         CALL LCMGET(IPGEOM,'RADIUS',RAYZON)
      CALL LCMSIX(IPGEOM,' ',2)
      RAYRE(IPAS+IKK)=0.0
      RINT=0.0
      RRINT=0.0D0
      RJN=0.0D0
      RJ=0.0D0
      DO 40 N=1,NZONE
      RRZON=RAYZON(N+1)
      PAS=(RRZON-RRINT)/DBLE(ABS(NBREG(N)))
      IF(NBREG(N).LT.0) PAS=PAS*(RRZON+RRINT)
      DO 30 I=1,ABS(NBREG(N))
      J=I+IPAS
      IF(J.GT.MAXPTS) THEN
         WRITE (HSMG,200) J,MAXPTS
         CALL XABORT(HSMG)
      ENDIF
      IF(NBREG(N).GT.0) THEN
         RJ1=RJ+PAS
         RJN1=RJ1**2
      ELSE
         RJN1=RJ*RJ+PAS
         RJ1=SQRT(RJN1)
      ENDIF
      VOL(J)=PI*REAL(RJN1-RJN)*POURCE(IKK)
      RAYRE(J+IKK)=REAL(RJ1)
      RJ=RJ1
      RJN=RJN1
      MAT(J)=MILIEU(N)
   30 CONTINUE
      IPAS=IPAS+ABS(NBREG(N))
      RRINT=RRZON
      RINT=RAYZON(N+1)
   40 CONTINUE
      NMC(IKK+1)=IPAS
      RINT=RAYRE(IPAS+IKK)
      SURFA(IKK)=2.0*PI*RINT*POURCE(IKK)
      SURFAT=SURFAT+SURFA(IKK)
   50 CONTINUE
*
      IF(ABS(PTOT).GT.EPS1) CALL XABORT('READMT: INVALID CELL PROPORT'
     1 //'IONS.')
      IF(ISTAT.EQ.1) THEN
*        STATISTICAL OPTION.
         DO 70 IKK=1,NSUPCE
         DO 60 JKK=1,NSUPCE
         PROCEL((JKK-1)*NSUPCE+IKK)=SURFA(JKK)/SURFAT
   60    CONTINUE
   70    CONTINUE
      ENDIF
      DO 85 IKK=1,NSUPCE
      DO 80 JKK=1,IKK-1
      X1=SURFA(IKK)*PROCEL((JKK-1)*NSUPCE+IKK)
      X2=SURFA(JKK)*PROCEL((IKK-1)*NSUPCE+JKK)
      IF(ABS(X1-X2).GT.EPS2) THEN
         WRITE (HSMG,210) ABS(X1-X2),IKK,JKK
         CALL XABORT(HSMG)
      ENDIF
   80 CONTINUE
   85 CONTINUE
*
      IR=0
      DO 90 I=1,IPAS
      IR=MAX(IR,MAT(I))
   90 CONTINUE
      ILK=.FALSE.
*
      IF(IMPX.GT.0) THEN
         WRITE (6,'(/46H DO-IT-YOURSELF MULTICELL OPTION BASED ON GEOM,
     1   5HETRY ,A12,1H./)') GEONAM
         WRITE (6,'(/34H USER SUPPLIED GEOMETRICAL MATRIX:/)')
         DO 100 I=1,NSUPCE
         WRITE (6,'(1X,1P,10E12.4)') (PROCEL((J-1)*NSUPCE+I),J=1,NSUPCE)
  100    CONTINUE
         I1=0
         DO 120 IKK=1,NSUPCE
         I2=NMC(IKK+1)-NMC(IKK)
         WRITE (6,'(//)')
         WRITE (6,230) IKK
         J1=I1+1
         DO 110 I=1,1+(I2-1)/8
         J2=MIN0(J1+7,I1+I2)
         WRITE (6,240) (J,J=J1,J2)
         WRITE (6,250) (RAYRE(J+IKK),J=J1,J2)
         J1=J1+8
  110    CONTINUE
         I1=I1+I2
  120    CONTINUE
         WRITE (6,'(//)')
         WRITE (6,260) IPAS,MAXPTS,IR
         WRITE (6,'(17H INFINITE DOMAIN./)')
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(RAYZON)
      DEALLOCATE(IGEN,CELL,MILIEU,NBREG)
      RETURN
*
  200 FORMAT (24HREADMT: A CELL NUMBER (=,I4,18H) EXCEED MAXPTS (=,
     1 I4,1H))
  210 FORMAT (53HREADMT: THE USER SUPPLIED GEOMETRICAL MATRIX DO NOT M,
     1 31HEET THE RECIPROCITY CONDITION (,1P,E10.1,7H); CELL,I5,4H <->,
     2 5H CELL,I5,1H.)
  230 FORMAT (1X,11(1H*)/7H * CELL,I3,2H */1X,11(1H*))
  240 FORMAT (/11H VOLUME    ,8(4X,I4,6X,1HI))
  250 FORMAT ( 11H ABSCISSA  ,8(F12.6,2X,1HI))
  260 FORMAT (/20H NUMBER OF VOLUMES =,I4,5X,23HAVAILABLE STORAGE: MAXB,
     1 4HLK =,I4/21H NUMBER OF MIXTURES =,I4/)
      END
