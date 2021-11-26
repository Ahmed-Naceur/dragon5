*DECK EPCRMV
      SUBROUTINE EPCRMV(IPEPC,IPCOV,IPRINT,IFMT,NGR,NIS,NXS,NCV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Extract variances and covariances from database and store on
* EPC data structure.
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPEPC   pointer to EPC deat structure.
* IPCOV   pointer to vaqriance and cavariance file.
* IPRINT  print level.
* IFMT    format of covariance file:
*         = 1 for ASCII file;
*         =-1 for BINARY file.
* NGR     number of groups.
* NIS     number of isotopes.
* NXS     number of cross section types per.
* NCV     maximum dimension of symmetrized covariance matrix.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPEPC
      INTEGER          IPCOV,IPRINT,IFMT,NGR,NIS,NXS,NCV
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='EPCRMV')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
*----
*  Local variables
*----
      INTEGER          IPRTL,ISO,NTYPE,ITYPE,NXSR,IXSR,IFCV,IPOC,
     >                 ILCV,ICMG,IGR,JGR
      CHARACTER        ISONAM*12,UNAME*8,RECNAM*12,FNAME*50,XSN*8
      INTEGER          ITC,NEL
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIS,ICOV
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NAMISO,IDXS
      REAL, ALLOCATABLE, DIMENSION(:) :: VAR,COV
*----
*  Scratch storage allocation
*   IDIS    array containing the isotope ID.
*   NAMISO  array containing the isotope names.
*   IDXS    array containing the cross section types (names).
*   VAR     array to store the variances.
*   ICOV    array to store indices to reconstructe full covariance
*           matrix from compressed covariance matrix.
*   COV     array to store compressed covariance matrix.
*----
      ALLOCATE(IDIS(NIS),NAMISO(3,NIS),IDXS(2,NXS),ICOV(NGR))
      ALLOCATE(VAR(NGR),COV(NCV))
*----
*  Scan over isotopes
*----
      IPRTL=IPRINT
      IF(IPRTL .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      IFCV=NGR+1
      NXSR=0
      CALL LCMSIX(IPEPC,'XSVariances ',ILCMUP)
      DO ISO=1,NIS
*----
*  Get isotope ID
*----
        IF(IFMT .GT. 0) THEN
          READ(IPCOV,1000) IDIS(ISO),UNAME,NTYPE,FNAME
        ELSE
          READ(IPCOV) IDIS(ISO),UNAME,NTYPE,FNAME
        ENDIF
        IF(IDIS(ISO) .GT. 999) THEN
          WRITE(ISONAM,'(I4,8X)') IDIS(ISO)
        ELSE IF(IDIS(ISO) .GT. 99) THEN
          WRITE(ISONAM,'(I3,9X)') IDIS(ISO)
        ELSE IF(IDIS(ISO) .GT. 9) THEN
          WRITE(ISONAM,'(I2,10X)') IDIS(ISO)
        ELSE
          WRITE(ISONAM,'(I1,11X)') IDIS(ISO)
        ENDIF
        READ(ISONAM,'(3A4)') (NAMISO(ITC,ISO),ITC=1,3)
        CALL LCMSIX(IPEPC,ISONAM,ILCMUP)
        DO ITYPE=1,NTYPE
*----
*  Get xs name and verify if in the list
*----
          IF(IFMT .GT. 0) THEN
            READ(IPCOV,1001) UNAME
          ELSE
            READ(IPCOV) UNAME
          ENDIF
          DO IXSR=1,NXSR
            WRITE(XSN,'(2A4)') IDXS(1,IXSR),IDXS(2,IXSR)
            IF(XSN .EQ. UNAME) GO TO 100
          ENDDO
          NXSR=NXSR+1
          IF(NXSR .GT. NXS) CALL XABORT(NAMSBR//
     >': number of cross section types insufficient')
          READ(UNAME,'(2A4)') IDXS(1,NXSR),IDXS(2,NXSR)
 100      CONTINUE
*----
*  Get variances and covariances
*----
          IF(IFMT .GT. 0) THEN
            READ(IPCOV,*) (VAR(IGR),IGR=1,NGR)
            READ(IPCOV,*) (COV(IGR),IGR=IFCV,NCV)
          ELSE
            READ(IPCOV) (VAR(IGR),IGR=1,NGR)
            READ(IPCOV) (COV(IGR),IGR=IFCV,NCV)
          ENDIF
*----
*  Compress variance and covariance matrix
*----
          IPOC=1
          ILCV=IFCV-1
          DO IGR=1,NGR
*----
*  Store variance for next element
*----
            COV(IPOC)=0.01*VAR(IGR)
            ICMG=0
*----
*  Scan covariance and remove trailing 0.0
*  Start at the end of COV for group IGR
*----
            DO JGR=NGR-IGR,1,-1
              IF(ICMG .EQ. 0) THEN
                IF(COV(ILCV+JGR) .NE. 0.0) THEN
*----
* First non 0.0 elements
* Add at the correct position in COV
*----
                  ICMG=ICMG+1
                  COV(IPOC+JGR)=COV(ILCV+JGR)
                ENDIF
              ELSE
*----
* Other elements including 0.0
* Add at the correct position in COV
*----
                ICMG=ICMG+1
                COV(IPOC+JGR)=COV(ILCV+JGR)
              ENDIF
            ENDDO
            ILCV=ILCV+NGR-IGR
            IPOC=IPOC+ICMG+1
            ICOV(IGR)=ICMG+1
          ENDDO
          NEL=IPOC-1
          RECNAM=UNAME//'    '
          CALL LCMPUT(IPEPC,RECNAM,NEL,2,COV)
          RECNAM='INDX'//UNAME
          CALL LCMPUT(IPEPC,RECNAM,NGR,1,ICOV)
        ENDDO
        CALL LCMSIX(IPEPC,ISONAM,ILCMDN)
      ENDDO
      CALL LCMPUT(IPEPC,'NAMEXS      ',2*NXSR,3,IDXS)
      CALL LCMPUT(IPEPC,'NAMEISO     ',3*NIS,3,NAMISO)
      CALL LCMSIX(IPEPC,'XSVariances ',ILCMDN)
      IF(IPRTL .GE. 10) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(COV,VAR)
      DEALLOCATE(ICOV,IDXS,NAMISO,IDIS)
      RETURN
*----
*  Formats
*----
 1000 FORMAT(I8,5X,A8,5X,I8,5X,A50)
 1001 FORMAT(A8)
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
