*DECK EPCRMU
      SUBROUTINE EPCRMU(IPEPC,IPMIC,IPRINT,NGR,NIS,NXS,NCV,
     >                  NBISO,NMIXT,NIFISS,ITOTL,ISCAT,NAMDXS,
     >                  NAMISO,NISOU,ISOMIX,IDVF,IDMF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Update MACROLIB from random error distribution.
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
* IPEPC   pointer to EPC data structure.
* IPMIC   pointer to MICROLIB.
* IPRINT  print level.
* NGR     number of groups.
* NIS     number of isotopes on EPC.
* NXS     number of cross section types on EPC.
* NCV     maximum dimension of symmetrized covariance matrix.
* NBISO   number of isotopes on MICROLIB.
* NMIXT   number of mixtures on MICROLIB.
* NIFISS  number of fissiles isotopes on MICROLIB.
* ITOTL   position of total XS in NAMDXS.
* ISCAT   position of scattering XS in NAMDXS.
* NAMDXS  names of XS.
* NAMISO  array containing the isotope names.
* NISOU   MICROLIB isotopes used.
* ISOMIX  MICROLIB isotopes mixtures.
* IDVF    variance isotopes to analyze and fission id.
* IDMF    MICROLIB isotopes to analyze and fission id.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPEPC,IPMIC
      INTEGER          IPRINT,NGR,NIS,NXS,NCV,
     >                 NBISO,NMIXT,NIFISS,ITOTL,ISCAT
      CHARACTER*6      NAMDXS(NXS)
      INTEGER          NAMISO(3,NIS),NISOU(3,NBISO),ISOMIX(NBISO),
     >                 IDVF(2,NIS),IDMF(2,NBISO)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='EPCRMU')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
*----
*  Local variables
*----
      INTEGER          ILCMLN,ILCMTY,IPRTL
      INTEGER          ISO,IXS,JSO,IGR,ISOF
      CHARACTER        ISONAM*12,RECNAM*12,NAMMIC*12
      INTEGER          ITC,IMIX
      REAL             DENSI
      INTEGER          ISEED,IGS,IORD,MINLEG,MAXLEG
      SAVE             ISEED,IGS,IORD,MINLEG,MAXLEG
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ICOV
      REAL, ALLOCATABLE, DIMENSION(:) :: DENS,XSMAC,XSREC,VAR,VAROLD,RST
      REAL, ALLOCATABLE, DIMENSION(:,:) :: COV
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DRVAR
*----
*  Data statement
*----
      DATA             ISEED,IGS,IORD,MINLEG,MAXLEG
     >                /0,    -1, 1,   0,     0/
*----
*  Scratch storage allocation
*   ICOV    array to store indices to reconstructe full covariance
*           matrix from compressed covariance matrix.
*   COV     array to store compressed covariance matrix.
*   DENS    MICROLIB isotopes densities
*   DRVAR   Random variance distribution (width=1)
*----
      ALLOCATE(ICOV(NGR,NXS))
      ALLOCATE(DENS(NBISO),COV(NCV,NXS))
      ALLOCATE(DRVAR(NGR,NXS))
*----
*  Write header
*----
      IPRTL=IPRINT
      IF(IPRTL .GE. 2) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Isotope densities
*----
      ALLOCATE(XSMAC(NGR*NXS*NMIXT*NIFISS))
      CALL XDRSET(XSMAC,NGR*NXS*NMIXT*NIFISS,0.0)
      ALLOCATE(XSREC(NGR*(NXS+1)),VAR(NGR),VAROLD(NGR),RST(NGR))
      CALL LCMGET(IPMIC,'ISOTOPESDENS',DENS)
      DO ISO=1,NIS
*----
*  Test if isotope used in Microlib
*----
        IF(IDVF(1,ISO) .GT .0) THEN
*----
*  Isotope is used
*  read covariance matrices
*----
          ICOV(:NGR,:NXS)=0
          COV(:NCV,:NXS)=0.0
          WRITE(ISONAM,'(3A4)') (NAMISO(ITC,ISO),ITC=1,3)
          IF(IPRTL .GE. 2) THEN
            WRITE(IOUT,*) ISONAM
          ENDIF
          CALL LCMSIX(IPEPC,ISONAM,ILCMUP)
          DO IXS=1,NXS
*----
*  Get covariance matrices
*----
            RECNAM='INDX'//NAMDXS(IXS)//'  '
            CALL LCMLEN(IPEPC,RECNAM,ILCMLN,ILCMTY)
            IF(ILCMLN .EQ. NGR) THEN
              CALL LCMGET(IPEPC,RECNAM,ICOV(1,IXS))
              RECNAM=NAMDXS(IXS)//'      '
              CALL LCMGET(IPEPC,RECNAM,COV(1,IXS))
*----
*  Generate random numbers from normal distribution
*----
              CALL RANDDN(ISEED,NGR,DRVAR(1,IXS))
              IF(IPRTL .GE. 5) THEN
                WRITE(IOUT,*) NAMDXS(IXS),'DRVAR=['
                WRITE(IOUT,6010) (DRVAR(IGR,IXS),IGR=1,NGR)
                WRITE(IOUT,*) '];'
              ENDIF
            ENDIF
          ENDDO
          CALL LCMSIX(IPEPC,ISONAM,ILCMDN)
*----
*  Scan over Microlib isotopes associated with this variance isotope
*----
          DO JSO=1,NBISO
            IF(IDMF(1,JSO) .EQ. ISO) THEN
              ISOF=IDMF(2,JSO)
*----
*  Read microlib for isotope
*----
              WRITE(NAMMIC,'(3A4)') (NISOU(ITC,JSO),ITC=1,3)
              CALL LCMSIX(IPMIC,NAMMIC,ILCMUP)
*----
*  Get microscopic xs
*----
              CALL XDRLXS(IPMIC,IGS,IPRINT,NXS,NAMDXS,IORD,NGR,XSREC)
*----
*  Add contribution to macrolib
*----
              DENSI=DENS(JSO)
              IMIX=ISOMIX(JSO)
              CALL LCMSIX(IPMIC,'VARIANCES   ',ILCMUP)
              CALL EPCRMA(IPMIC,IPRINT,NGR,NXS,NCV,NMIXT,NIFISS,
     >                    IMIX,ISOF,ITOTL,ISCAT,NAMDXS,DENSI,ICOV,
     >                    COV,DRVAR,XSREC,XSMAC,VAR,VAROLD,RST)
              CALL LCMSIX(IPMIC,'VARIANCES   ',ILCMDN)
              CALL LCMSIX(IPMIC,NAMMIC,ILCMDN)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IF(IPRTL .GE. 2) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      DEALLOCATE(RST,VAROLD,VAR,XSREC)
*----
*  Update macrolib
*----
      CALL LCMSIX(IPMIC,'MACROLIB    ',ILCMUP)
      CALL EPCRMS(IPMIC,IPRINT,NGR,NXS,NMIXT,NIFISS,NAMDXS,XSMAC)
      CALL LCMSIX(IPMIC,'MACROLIB    ',ILCMDN)
      DEALLOCATE(XSMAC)
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(DRVAR)
      DEALLOCATE(COV,DENS)
      DEALLOCATE(ICOV)
      RETURN
*----
*  Formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(10F10.4,:,' ...')
      END
