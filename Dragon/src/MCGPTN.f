*DECK MCGPTN
      SUBROUTINE MCGPTN(IMPX,NREG,NSOU,NANGL,NMU,VOLSUR,VOLNUM,SURNUM,
     1                  DENSTY,WZMU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute errors and normalization factors for 3D prismatic extended
* tracking (adapted from NXTTLS.f).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IMPX    print flag.
* NREG    number of regions.
* NSOU    number of external surfaces.
* NANGL   number of plan tracking angles.
* NMU     number of polar angles.
* SURNUM  numerical surfaces.
* DENSTY  plan tracking track density per angle.
* WZMU    polar quadrature weights.
*
*Parameters: input/output
* VOLSUR  analytical volume and surfaces
*         / analytical volumes and numerical surfaces.
* VOLNUM  numerical volumes per angle / normalization factors per angle.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,NREG,NSOU,NANGL,NMU
      REAL VOLSUR(NREG+NSOU),WZMU(NMU)
      DOUBLE PRECISION DENSTY(NANGL),VOLNUM(NREG,NANGL,NMU,2),
     >                 SURNUM(NSOU)
*----
*  LOCAL VARIABLES
*----
      INTEGER IOUT
      CHARACTER NAMSBR*6,CSGMU(2)*1
      PARAMETER (IOUT=6,NAMSBR='MCGPTN')
      INTEGER IR,IS,IANGL,IMU,ISGMU,NBVERG,NBVERA,NBSERR,NBV0,ITDIR
      REAL RTEMP
      DOUBLE PRECISION DCERR,DSVERG,DMVERG,DAVERG,DSVERA,DMVERA,DAVERA,
     1 DSSERR,DMSERR,DASERR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: VGLNUM
      DATA CSGMU / '+','-' /
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(VGLNUM(NREG))
      CALL XDDSET(VGLNUM,NREG,0.D0)
*
      ITDIR=0
      DO 10 IANGL=1,NANGL
      IF(DENSTY(IANGL).EQ.0.0D0) GOTO 10
      ITDIR=ITDIR+1
      DO 20 IMU=1,NMU
      DO 30 ISGMU=1,2
         NBVERA=0
         DSVERA=0.D0
         DMVERA=0.D0
         DAVERA=0.D0
         NBV0=0
         DO IR=1,NREG
            IF(VOLNUM(IR,ITDIR,IMU,ISGMU).EQ.0.D0) THEN
               IF(IMPX.GE.10) WRITE(IOUT,300) NAMSBR,IR,ITDIR,IMU,ISGMU
               VOLNUM(IR,ITDIR,IMU,ISGMU)=1.D0
               NBV0=NBV0+1
            ELSE
               VGLNUM(IR)=VGLNUM(IR)
     1                   +2.0*WZMU(IMU)*VOLNUM(IR,ITDIR,IMU,ISGMU)
               VOLNUM(IR,ITDIR,IMU,ISGMU)=VOLNUM(IR,ITDIR,IMU,ISGMU)
     1                                  *DENSTY(IANGL)
               VOLNUM(IR,ITDIR,IMU,ISGMU)=DBLE(VOLSUR(IR))
     1                             /VOLNUM(IR,ITDIR,IMU,ISGMU)
               NBVERA=NBVERA+1
            ENDIF
            DCERR=100.0D0*(1.D0-VOLNUM(IR,ITDIR,IMU,ISGMU))
            DMVERA=MAX(DMVERA,ABS(DCERR))
            DSVERA=DSVERA+DCERR*DCERR
            DAVERA=DAVERA+DCERR           
         ENDDO
         DSVERA=SQRT(DSVERA/DBLE(NBVERA))
         DAVERA=DAVERA/DBLE(NBVERA)
         IF(NBV0.GT.0) THEN
            WRITE(IOUT,500) NAMSBR,NBV0,ITDIR,IMU,CSGMU(ISGMU)
         ELSE
            IF(IMPX.GE.2) WRITE(IOUT,100) DSVERA,DMVERA,DAVERA,ITDIR,
     1                                     IMU,CSGMU(ISGMU)
         ENDIF
 30   CONTINUE
 20   CONTINUE
 10   CONTINUE
*
      NBVERG=0
      DSVERG=0.D0
      DMVERG=0.D0
      DAVERG=0.D0
      NBV0=0
      DO IR=1,NREG
         IF(VGLNUM(IR).EQ.0.D0) THEN
            WRITE(IOUT,300) NAMSBR,IR
            VGLNUM(IR)=1.D0
            NBV0=NBV0+1
         ELSE
            VGLNUM(IR)=DBLE(VOLSUR(IR))/VGLNUM(IR)
            NBVERG=NBVERG+1
         ENDIF
         DCERR=100.0D0*(1.D0-VGLNUM(IR))
         DMVERG=MAX(DMVERG,ABS(DCERR))
         DSVERG=DSVERG+DCERR*DCERR
         DAVERG=DAVERG+DCERR           
      ENDDO
      IF(NBV0.GT.0) THEN
         WRITE(IOUT,500) NAMSBR,NBV0
      ENDIF
      DSVERG=SQRT(DSVERG/DBLE(NBVERG))
      DAVERG=DAVERG/DBLE(NBVERG)
      IF(IMPX.GE.1) WRITE(IOUT,150) DSVERG,DMVERG,DAVERG
*
      NBSERR=0
      DSSERR=0.D0
      DMSERR=0.D0
      DASERR=0.D0
      DO IR=1,NSOU
         IS=NREG+IR
         IF(SURNUM(IR).EQ.0.D0) THEN
            WRITE(IOUT,400) NAMSBR,-IR
            SURNUM(IR)=1.D0
         ELSE
            RTEMP=VOLSUR(IS)
            VOLSUR(IS)=REAL(SURNUM(IR))
            SURNUM(IR)=RTEMP/VOLSUR(IS)
            NBSERR=NBSERR+1
         ENDIF
         DCERR=100.0D0*(1.D0-SURNUM(IR))
         DMSERR=MAX(DMSERR,ABS(DCERR))
         DSSERR=DSSERR+DCERR*DCERR
         DASERR=DASERR+DCERR
      ENDDO
      DSSERR=SQRT(DSSERR/DBLE(NBSERR))
      DASERR=DASERR/DBLE(NBSERR)
      IF(IMPX.GE.1) WRITE(IOUT,200) DSSERR,DMSERR,DASERR
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VGLNUM)
      RETURN
*
 100  FORMAT(' Angular RMS, maximum and average errors (%) ',
     1       'on region volumes :',3(2X,F10.5),
     2       ' for plan tracking angle',I3,' for polar angle',I3,
     3       '(',A1,')')
 150  FORMAT(' Global RMS, maximum and average errors (%) ',
     1       'on region volumes :',3(2X,F10.5))
 200  FORMAT(' Global RMS, maximum and average errors (%) ',
     1       'on surface areas  :',3(2X,F10.5))
 300  FORMAT(1X,'***** Warning in ',A6,'*****'/
     1       7X,'For region ',I8,
     2       1X,'no crossing by angle ',I8,I8,'(',I1,')')
 400  FORMAT(1X,'***** Warning in ',A6,'*****'/
     1       7X,'For surface ',I8,
     2       1X,'no crossing by any angle ')
 500  FORMAT(1X,'***** Warning in ',A6,'*****'/
     1       7X,I8,' regions not tracked for direction  ',I8,I8,
     2       '(',A1,')')
      END
