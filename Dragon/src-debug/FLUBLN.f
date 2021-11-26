*DECK FLUBLN
      SUBROUTINE FLUBLN(IPMACR,IPRINT,NGROUP,NBMIX,NREGIO,NUNKNO,
     >                  NIFISS,MATCOD,VOLUME,KEYFLX,FUNKNO,IHETL,
     >                  REFKEF,B2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of a directional buckling from the critical neutron
* balance.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): I. Petrovic and G. Marleau
*
*Parameters: input
* IPMACR  pointer to the macrolib LCM object.
* IPRINT  print selection for flux modules.
* NGROUP  number of groups.
* NBMIX   number of mixtures.
* NREGIO  number of regions.
* NUNKNO  number of unknowns in the system.
* NIFISS  number of fissile isotopes.
* MATCOD  material code in regions.
* IHETL   type of buckling calculation:
*         = 1 x-direction search;
*         = 2 y-direction search;
*         = 3 z-direction search;
*         = 4 r-direction search (X=Y);
*         = 5 global-direction search (X=Y=Z).
* VOLUME  volume of regions.
* KEYFLX  flux elements in unknown system.
* FUNKNO  flux and directional currents.
* REFKEF  target K-effective for type B or L.
*
*Parameters: output
* B2      directional buckling (X, Y, Z, hom).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR
      INTEGER   IPRINT,NGROUP,NBMIX,NREGIO,NUNKNO,NIFISS,MATCOD(NREGIO),
     >          KEYFLX(NREGIO),IHETL
      REAL      VOLUME(NREGIO),FUNKNO(NUNKNO,NGROUP),B2(4)
      DOUBLE PRECISION REFKEF
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IUNOUT=6)
      TYPE(C_PTR) JPMACR,KPMACR
      DOUBLE PRECISION BIL1,SUM(0:3)
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGT0,SIGS0
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGFIS,QTOTL
*----
*  COMPUTE THE TOTAL NEUTRON PRODUCTION
*----
      ALLOCATE(SIGFIS(NBMIX,NIFISS),QTOTL(NREGIO,NIFISS))
      NUN4=NUNKNO/4
      CALL XDRSET(QTOTL,NREGIO*NIFISS,0.0)
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO 30 IGR=1,NGROUP
        KPMACR=LCMGIL(JPMACR,IGR)
        CALL LCMGET(KPMACR,'NUSIGF',SIGFIS)
        DO 20 IFIS=1,NIFISS
          DO 10 IREG=1,NREGIO
            IBM=MATCOD(IREG)
            IF(IBM.GT.0) QTOTL(IREG,IFIS)=QTOTL(IREG,IFIS)
     >      +FUNKNO(KEYFLX(IREG),IGR)*SIGFIS(IBM,IFIS)
 10       CONTINUE
 20     CONTINUE
 30   CONTINUE
      BIL1=0.0D0
      DO 60 IGR=1,NGROUP
        KPMACR=LCMGIL(JPMACR,IGR)
        CALL LCMGET(KPMACR,'CHI',SIGFIS)
        DO 50 IFIS=1,NIFISS
          DO 40 IREG=1,NREGIO
            IBM=MATCOD(IREG)
            IF(IBM.GT.0) BIL1=BIL1+DBLE(VOLUME(IREG)*QTOTL(IREG,IFIS)*
     >      SIGFIS(IBM,IFIS))
 40       CONTINUE
 50     CONTINUE
 60   CONTINUE
      DEALLOCATE(QTOTL,SIGFIS)
*----
*  COMPUTE FISSION SOURCE AND EVALUATE NEUTRON BALANCE
*----
      ALLOCATE(SIGT0(0:NBMIX),SIGS0(0:NBMIX))
      SUM(0)=BIL1/REFKEF
      SUM(1)=0.0D0
      SUM(2)=0.0D0
      SUM(3)=0.0D0
      SIGT0(0)=0.0
      SIGS0(0)=0.0
      DO 80 IGR=1,NGROUP
        KPMACR=LCMGIL(JPMACR,IGR)
        CALL LCMGET(KPMACR,'NTOT0',SIGT0(1))
        CALL LCMGET(KPMACR,'SIGS00',SIGS0(1))
        DO 70 IREG=1,NREGIO
          IBM=MATCOD(IREG)
          IND=KEYFLX(IREG)
          SUM(0)=SUM(0)+(SIGS0(IBM)-SIGT0(IBM))*
     >              VOLUME(IREG)*FUNKNO(IND,IGR)
          SUM(1)=SUM(1)+VOLUME(IREG)*FUNKNO(NUN4+IND,IGR)
          SUM(2)=SUM(2)+VOLUME(IREG)*FUNKNO(2*NUN4+IND,IGR)
          SUM(3)=SUM(3)+VOLUME(IREG)*FUNKNO(3*NUN4+IND,IGR)
 70     CONTINUE
 80   CONTINUE
      IF(IHETL.EQ.1)THEN
        B2(1)=REAL((SUM(0)-B2(2)*SUM(2)-B2(3)*SUM(3))/SUM(1))
      ELSEIF(IHETL.EQ.2)THEN
        B2(2)=REAL((SUM(0)-B2(1)*SUM(1)-B2(3)*SUM(3))/SUM(2))
      ELSEIF(IHETL.EQ.3)THEN
        B2(3)=REAL((SUM(0)-B2(1)*SUM(1)-B2(2)*SUM(2))/SUM(3))
      ELSEIF(IHETL.EQ.4)THEN
        B2(1)=REAL((SUM(0)-B2(3)*SUM(3))/(SUM(1)+SUM(2)))
        B2(2)=B2(1)
      ELSEIF(IHETL.EQ.5)THEN
        B2(1)=REAL(SUM(0)/(SUM(1)+SUM(2)+SUM(3)))
        B2(2)=B2(1)
        B2(3)=B2(1)
      ELSE
        CALL XABORT('FLUBLN: WHICH DIRECTIONAL BUCKLING '//
     >              'WOULD YOU LIKE TO CALCULATE ? ')
      ENDIF
      B2(4)=B2(1)+B2(2)+B2(3)
      IF(IPRINT.GE.10) WRITE(IUNOUT,6000) (B2(IDIR),IDIR=1,3)
      DEALLOCATE(SIGS0,SIGT0)
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT(1X,'FLUBLN OUTPUT'/1X,'HETEROGENEOUS B2 = ',1P,3E15.7)
      END
