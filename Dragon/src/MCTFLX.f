*DECK MCTFLX
      SUBROUTINE MCTFLX(IPTRK,IPOUT,IPRINT,NMIX,NGRP,NL,NFM,NDEL,NED,
     <                  NAMEAD,XSTOT,XSS,XSSNN,XSNUSI,XSCHI,XSN2N,
     <                  XSN3N,XSEDI,NSRCK,IKZ,KCT,ISEED,XYZL,NBSCO,
     <                  NMERGE,NGCOND,KEFF,REKEFF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Power iteration with the Monte Carlo method in 1D/2D/3D Cartesian
* geometry.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): B. Arsenault
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* IPOUT   pointer to the MC data structure.
* IPRINT  print flag.
* NMIX    number of mixtures in the geometry.
* NGRP    number of energy groups.
* NL      number of Legendre orders required in the estimations
*         (NL=1 or higher).
* NFM     number of fissile isotopes.
* NDEL    number of delayed precursor groups.
* NED     number of extra edit vectors.
* NAMEAD  names of these extra edits.
* XSTOT   total macroscopic cross sections for each mixture and energy
*         group.
* XSS     total scattering cross sections for each mixture and energy
*         group.
* XSSNN   in-group and out-of-group macroscopic transfert cross sections
*         for each mixture.
* XSNUSI  the values of Nu time the fission cross sections for each
*         isotope per mixture and energy group.
* XSCHI   the values of fission spectrum per isotope per mixture for
*         each energy group.
* XSN2N   N2N macroscopic cross sections for each mixture and energy
*         group.
* XSN3N   N3N macroscopic cross sections for each mixture and energy
*         group.
* XSEDI   extra edit cross sections for each mixture and energy group.
* NSRCK   number of neutrons generated per cycle.
* IKZ     number of inactive cycles.
* KCT     number of active cycles.
* ISEED   the seed for the generation of random numbers.
* XYZL    Cartesian boundary coordinates.
* NBSCO   number of macrolib-related scores.
* NMERGE  number of homogenized regions.
* NGCOND  number of condensed energy groups.
*
*Parameters: output
* KEFF    effective multiplication factor.
* REKEFF  standard deviation on the effective multiplication factor.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPOUT
      INTEGER IPRINT,NMIX,NGRP,NL,NFM,NDEL,NED,NAMEAD(2,NED),NSRCK,IKZ,
     <        KCT,ISEED,NBSCO,NMERGE,NGCOND
      REAL    XSTOT(NMIX,NGRP),XSS(NMIX,NGRP,NL),XSN2N(NMIX,NGRP),
     <        XSN3N(NMIX,NGRP),XSSNN(NGRP,NGRP,NMIX,NL),
     <        XSCHI(NMIX,NFM,NGRP,1+NDEL),XSNUSI(NMIX,NFM,NGRP,1+NDEL),
     <        XSEDI(NMIX,NGRP,NED)
      DOUBLE PRECISION XYZL(2,3),KEFF,REKEFF
      CHARACTER NAMREC*12
*----
*  LOCAL VARIABLES 
*----
      INTEGER NSTATE
      PARAMETER(NSTATE=40)
      INTEGER ICYCLE,ILOOP,NLOOP,IDIR,IFIRST,T1,T2,MIX,ISONBR,NFREG,
     <        NFSUR,ANGBC,NTRK,JJ,ITYPBC,ITRK,IDIRG,NBOCEL,NBUCEL,
     <        IDIAG,ISAXIS(3),NOCELL(3),NUCELL(3),ICODE(6),MXMSH,
     <        MAXREG,NBTCLS,MAXPIN,MAXMSP,MAXRSP,MXGSUR,MXGREG,NUNK,
     <        MAXMSH,NDIM,ESTATE(NSTATE),GSTATE(NSTATE),ITALLY,IND,IOF,
     <        IGR,ILON1,ILON2,ITYLCM,IBANK1,IBANK2
      REAL    ALBEDO(6),RAND,NUCALL,NULIMIT,SCORE1(3),FACT1,FACT2
      LOGICAL LKEEP
      DOUBLE PRECISION ABSC(3,2),POS(3),KCYCLE,WEIGHT,NU,SUM1,SUM2,
     <        ASCORE1(3),BSCORE1(3)
      CHARACTER HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INMIX,IBCRT,IUNFLD,INDEX,
     1 IDREG,ITPIN,INDGRP,IMERGE,IGCR
      REAL, ALLOCATABLE, DIMENSION(:) :: NUCYCLE
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCORE2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DGMESH,DCMESH,
     1 DRAPIN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ASCORE2,BSCORE2
      INTEGER, POINTER, DIMENSION(:,:) :: INGEN1,INGEN2,INGAR
      DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: DNGEN1,DNGEN2,DNGAR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SCORE2(NBSCO,NMERGE,NGCOND))
      ALLOCATE(ASCORE2(NBSCO,NMERGE,NGCOND),
     <         BSCORE2(NBSCO,NMERGE,NGCOND))
*----
*  SET THE RANDOM NUMBER GENERATOR
*----
      IFIRST=1
      IF(ISEED.EQ.0) THEN
         CALL CLETIM(T1,T2)
         ISEED=T1+T2
         DO JJ=0,MOD(ISEED,10)
           CALL RANDF(ISEED,IFIRST,RAND)
         ENDDO
      ENDIF
*----
*  RECOVER SOME BASIC NXT GEOMETRY ANALYSIS INFO AND ALLOCATE RELATED
*  MEMORY
*----     
      CALL XDISET(GSTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'STATE-VECTOR',GSTATE)
      NFREG    =GSTATE( 1)
      NMIX     =GSTATE( 4)
      NFSUR    =GSTATE( 5)
      IF(GSTATE(7).NE.4)
     1 CALL XABORT('MCTFLX: ONLY NXT: GEOMETRY ANALYSIS IS PERMITTED')
      ANGBC    =ABS(GSTATE(10))
*----
*  READ THE MATERIAL NUMBER ASSOCIATED TO EACH REGION NUMBER 
*----
      ALLOCATE(INMIX(NFREG),IBCRT(NFSUR))
      CALL LCMGET(IPTRK,'MATCOD',INMIX)
      CALL LCMGET(IPTRK,'BC-REFL+TRAN',IBCRT)
      CALL LCMGET(IPTRK,'ALBEDO',ALBEDO)
      CALL LCMGET(IPTRK,'ICODE',ICODE)
      CALL LCMSIX(IPTRK,'NXTRecords',1)
      CALL XDISET(ESTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'G00000001DIM',ESTATE)
      NDIM     =ESTATE( 1)
      ITYPBC   =ESTATE( 2)
      IDIRG    =ESTATE( 3)
      NBOCEL   =ESTATE( 4)
      NBUCEL   =ESTATE( 5)
      IDIAG    =ESTATE( 6)
      ISAXIS(1)=ESTATE( 7)
      ISAXIS(2)=ESTATE( 8)
      ISAXIS(3)=ESTATE( 9)
      NOCELL(1)=ESTATE(10)
      NOCELL(2)=ESTATE(11)
      NOCELL(3)=ESTATE(12)
      NUCELL(1)=ESTATE(13)
      NUCELL(2)=ESTATE(14)
      NUCELL(3)=ESTATE(15)
      MXMSH    =ESTATE(16)
      MAXREG   =ESTATE(17)
      NBTCLS   =ESTATE(18)
      MAXPIN   =ESTATE(19)
      MAXMSP   =ESTATE(20)
      MAXRSP   =ESTATE(21)
      IF(NFSUR.NE.ESTATE(22))
     1 CALL XABORT('MCTFLX: INCONSISTENT NUMBER OF OUTER SURFACES')
      IF(NFREG.NE.ESTATE(23))
     1 CALL XABORT('MCTFLX: INCONSISTENT NUMBER OF REGIONS')
      MXGSUR   =ESTATE(24)
      MXGREG   =ESTATE(25)
      NUNK=NFSUR+NFREG+1
      MAXMSH=MAX(MXMSH,MAXMSP,MAXREG)
*     cell index and orientation for the cells filling the geometry
      ALLOCATE(IUNFLD(2*NBUCEL))
      NAMREC='G00000001CUF'
      CALL LCMGET(IPTRK,NAMREC,IUNFLD)
*     global mesh for geometry
      ALLOCATE(DGMESH((MAXMSH+2)*4))
      CALL XDDSET(DGMESH,(MAXMSH+2)*4,0.D0)
*
*     An offset of 2 has been used to be compatible with the definition
*     of the pin geometries where the first two values give the offset
*     of the pin.  With a Cartesian geometry the array doesn't contain
*     the offests. 
      CALL NXTXYZ(IPTRK,IPRINT,NDIM,ITYPBC,MAXMSH,NUCELL,ABSC,DGMESH)
      DO IDIR=1,NDIM
        XYZL(1,IDIR)=ABSC(IDIR,2)-ABSC(IDIR,1)
        XYZL(2,IDIR)=ABSC(IDIR,2)
      ENDDO
      ALLOCATE(INDEX(2*5*(MXGSUR+MXGREG+1)),IDREG(2*MXGREG),
     1 ITPIN(3*(MAXPIN+1)))
      ALLOCATE(DCMESH(2*4*(MAXMSH+2)),DRAPIN(6*(MAXPIN+1)))
*----
*  TALLY INITIALIZATION
*----
      CALL LCMGET(IPOUT,'STATE-VECTOR',GSTATE)
      ITALLY=GSTATE(6)
      IF(ITALLY.EQ.2) THEN
        ALLOCATE(INDGRP(NGRP))
        CALL XDISET(INDGRP,NGRP,0)
        CALL LCMLEN(IPOUT,'REF:IMERGE',ILON1,ITYLCM)
        CALL LCMLEN(IPOUT,'REF:IGCOND',ILON2,ITYLCM)
        ALLOCATE(IMERGE(ILON1),IGCR(ILON2))
        CALL LCMGET(IPOUT,'REF:IMERGE',IMERGE)
        CALL LCMGET(IPOUT,'REF:IGCOND',IGCR)
        IOF=1
        JJ=IGCR(1)
        DO IND=1,NGRP
          IF(IND.GT.JJ) THEN
            IOF=IOF+1
            IF(IOF.GT.NGCOND) CALL XABORT('MCTFLX: NGCOND OVERFLOW.')
            JJ=IGCR(IOF)
          ENDIF
          INDGRP(IND)=IOF
        ENDDO
      ENDIF
*----
*  MEMORY ALLOCATION FOR THE POWER ITERATION
*----
      ALLOCATE(INGEN1(2,2*NSRCK),DNGEN1(4,2*NSRCK))
      ALLOCATE(INGEN2(2,2*NSRCK),DNGEN2(4,2*NSRCK))
      ALLOCATE(NUCYCLE(2*NSRCK))
      CALL XDISET(INGEN1,2*NSRCK*2,0)
      CALL XDISET(INGEN2,2*NSRCK*2,0)
      CALL XDDSET(DNGEN1,2*NSRCK*4,0.0D0)
      CALL XDDSET(DNGEN2,2*NSRCK*4,0.0D0)
*----
*  UNIFORM INITIAL ESTIMATE OF SOURCE NEUTRONS
*----
      DO ILOOP=1,NSRCK
        DO IDIR=1,NDIM
          CALL RANDF(ISEED,IFIRST,RAND)
          POS(IDIR)=RAND*(XYZL(2,IDIR)-XYZL(1,IDIR))+XYZL(1,IDIR)
        ENDDO
        INGEN1(1,ILOOP) = -1
        INGEN1(2,ILOOP) = -1
        DNGEN1(1,ILOOP) = 1.0D0
        DNGEN1(2,ILOOP) = POS(1)
        DNGEN1(3,ILOOP) = POS(2)
        DNGEN1(4,ILOOP) = POS(3)
      ENDDO
      IBANK1=NSRCK
*----
*  POWER ITERATION
*----
      KCYCLE=1.D0
      SUM1=0.0D0
      SUM2=0.0D0
      IF(ITALLY.GT.0) THEN
        CALL XDDSET(ASCORE1,3,0.0D0)
        CALL XDDSET(BSCORE1,3,0.0D0)
        IF(ITALLY.EQ.2) THEN
          CALL XDDSET(ASCORE2,NBSCO*NMERGE*NGCOND,0.0D0)
          CALL XDDSET(BSCORE2,NBSCO*NMERGE*NGCOND,0.0D0)
        ENDIF
      ENDIF
      DO ICYCLE=1,KCT
        IBANK2=0
        CALL XDRSET(SCORE1,3,0.0)
        CALL XDRSET(SCORE2,NBSCO*NMERGE*NGCOND,0.0)
        WEIGHT=KCYCLE
        KCYCLE=0.D0
        DO NTRK=1,IBANK1
          NU = DNGEN1(1,NTRK)
          NLOOP=MAX(1,INT(NU/WEIGHT))
          DNGEN1(1,NTRK)=NU/(DBLE(NLOOP)*WEIGHT)
          DO ILOOP=1,NLOOP
*----
*  TRACK EACH NEUTRON INDIVIDUALLY
*----
            MIX    = INGEN1(1,NTRK)
            ISONBR = INGEN1(2,NTRK)
            NUCALL = REAL(DNGEN1(1,NTRK))
            POS(1) = DNGEN1(2,NTRK)
            POS(2) = DNGEN1(3,NTRK)
            POS(3) = DNGEN1(4,NTRK)
            CALL MCTRK(IPTRK,IPRINT,NFREG,NFSUR,NDIM,NMIX,ANGBC,ITYPBC,
     1      MAXMSH,NUCELL,MXGSUR,MXGREG,MAXPIN,IBCRT,ICODE,ALBEDO,
     2      IUNFLD,DGMESH,XYZL,INDEX,IDREG,DCMESH,ITPIN,DRAPIN,ISEED,
     3      NGRP,NL,NFM,NDEL,NED,INMIX,XSTOT,XSS,XSN2N,XSN3N,XSSNN,
     4      XSNUSI,XSCHI,XSEDI,MIX,ISONBR,NUCALL,POS,ITALLY,NBSCO,
     5      NMERGE,NGCOND,IMERGE,INDGRP,SCORE1,SCORE2)
            IF(ISONBR.GT.0) THEN
              IBANK2=IBANK2+1
              IF(IBANK2.GT.2*NSRCK) THEN
                CALL XABORT('MCTFLX: TOO MANY NEUTRON TRACKS BEING'//
     1          ' BANKED.')
              ENDIF
              INGEN2(1,IBANK2) = MIX
              INGEN2(2,IBANK2) = ISONBR
              DNGEN2(1,IBANK2) = NUCALL
              DNGEN2(2,IBANK2) = POS(1)
              DNGEN2(3,IBANK2) = POS(2)
              DNGEN2(4,IBANK2) = POS(3)
              NUCYCLE(IBANK2) = NUCALL
              KCYCLE=KCYCLE+NUCALL
            ENDIF
          ENDDO
*         END OF THE NTRK CYCLE
        ENDDO
        KCYCLE=KCYCLE/DBLE(NSRCK)
*----
*  RUSSIAN ROULETTE
*----
        IF(IBANK2.GT.NSRCK) THEN
          CALL SORTRE(IBANK2,NUCYCLE)
          NULIMIT=NUCYCLE(IBANK2-NSRCK+1)
          IBANK1=0
          DO ITRK=1,IBANK2
            MIX    = INGEN2(1,ITRK)
            ISONBR = INGEN2(2,ITRK)
            NUCALL = REAL(DNGEN2(1,ITRK))
            POS(1) = DNGEN2(2,ITRK)
            POS(2) = DNGEN2(3,ITRK)
            POS(3) = DNGEN2(4,ITRK)
            LKEEP=(NUCALL.GE.NULIMIT)
            IF(.NOT.LKEEP) THEN
              CALL RANDF(ISEED,IFIRST,RAND)
              LKEEP=RAND.LE.(NUCALL/NULIMIT)
              NUCALL=NULIMIT
            ENDIF
            IF(LKEEP) THEN
              IBANK1=IBANK1+1
              INGEN1(1,IBANK1) = MIX
              INGEN1(2,IBANK1) = ISONBR
              DNGEN1(1,IBANK1) = NUCALL
              DNGEN1(2,IBANK1) = POS(1)
              DNGEN1(3,IBANK1) = POS(2)
              DNGEN1(4,IBANK1) = POS(3)
            ENDIF
          ENDDO
        ELSE
          IBANK1 = IBANK2
          INGAR  => INGEN2
          INGEN2 => INGEN1
          INGEN1 => INGAR
          DNGAR  => DNGEN2
          DNGEN2 => DNGEN1
          DNGEN1 => DNGAR
        ENDIF
*----
*  K-EFFECTIVE OF THE PROBLEM WITH THE RELATIVE ERROR
*----
        IF(ICYCLE.GT.IKZ) THEN
          SUM1=SUM1+KCYCLE
          SUM2=SUM2+KCYCLE**2
          IF(ITALLY.GT.0) THEN
            DO IND=1,3
              ASCORE1(IND)=ASCORE1(IND)+SCORE1(IND)
              BSCORE1(IND)=BSCORE1(IND)+SCORE1(IND)**2
            ENDDO
            IF(ITALLY.EQ.2) THEN
              DO IND=1,NBSCO
                DO MIX=1,NMERGE
                  DO IGR=1,NGCOND
                    ASCORE2(IND,MIX,IGR)=ASCORE2(IND,MIX,IGR)+
     1              SCORE2(IND,MIX,IGR)
                    BSCORE2(IND,MIX,IGR)=BSCORE2(IND,MIX,IGR)+
     1              SCORE2(IND,MIX,IGR)**2
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        IF(IPRINT.GT.0) WRITE(6,1000) ICYCLE,KCYCLE
      ENDDO
      FACT1=REAL(KCT-IKZ)
      FACT2=REAL(KCT-IKZ-1)
      KEFF=SUM1/FACT1
      REKEFF=SQRT((SUM2-FACT1*KEFF**2)/FACT1/FACT2)
      WRITE(6,2000) KCT,KCYCLE,KEFF,REKEFF
      IF(ITALLY.GT.0) THEN
        DO IND=1,3
          ASCORE1(IND)=ASCORE1(IND)/FACT1
          BSCORE1(IND)=SQRT((BSCORE1(IND)-FACT1*ASCORE1(IND)**2)/FACT1/
     1    FACT2)
        ENDDO
        KEFF=ASCORE1(2)/ASCORE1(3)
        REKEFF=KEFF*(BSCORE1(2)/ASCORE1(2)-BSCORE1(3)/ASCORE1(3))
        WRITE(6,3000) KEFF,REKEFF
      ENDIF
      IF(ITALLY.EQ.2) THEN
        DO IND=1,NBSCO
          DO MIX=1,NMERGE
            DO IGR=1,NGCOND
*---- 
*     CHECK FIRST FOR 0-TALLIES IN TOTAL CROSS-SECTIONS PER MIXTURE
*---- 
              IF(ASCORE2(1,MIX,IGR).EQ.0.0D0) THEN
                WRITE(HSMG,'(28HMCPTFLX: ZERO TALLY FOR MIX ,I5,
     1          10H IN GROUP ,I5,28H. INCREASE KCODE PARAMETERS.)')
     2          MIX,IGR
                CALL XABORT(HSMG)
              ENDIF
              ASCORE2(IND,MIX,IGR)=ASCORE2(IND,MIX,IGR)/FACT1
              BSCORE2(IND,MIX,IGR)=SQRT((BSCORE2(IND,MIX,IGR)-FACT1*
     1        ASCORE2(IND,MIX,IGR)**2)/FACT1/FACT2)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
*----
*  RECONSTRUCT THE TALLY-GENERATED MACROLIB
*----
      IF(ITALLY.EQ.2) THEN
        CALL MCTOUT(IPOUT,NL,NFM,NDEL,NED,NAMEAD,NBSCO,NMERGE,NGCOND,
     1  ASCORE1,ASCORE2)
        DEALLOCATE(IGCR,IMERGE)
      ENDIF
*----
*  DEALLOCATE MEMORY
*----
*     POWER ITERATION RELATED
      DEALLOCATE(NUCYCLE,DNGEN2,INGEN2,DNGEN1,INGEN1)
      IF(ITALLY.EQ.2) DEALLOCATE(INDGRP)
*
*     TRACKING RELATED
      CALL LCMSIX(IPTRK,' ',2)
      DEALLOCATE(DRAPIN,ITPIN,DCMESH,IDREG,INDEX,DGMESH,IUNFLD,IBCRT,
     1 INMIX)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(BSCORE2,ASCORE2)
      DEALLOCATE(SCORE2)
      RETURN
*
 1000 FORMAT(' CYCLE NUMBER: ',I5,' K-EFFECTIVE CYCLE: ',F8.6)
 2000 FORMAT(/' CYCLE NUMBER: ',I5,' K-EFFECTIVE CYCLE: ',F8.6,
     <        ' K-EFFECTIVE AVERAGE: ',F8.6,
     <        ' SIGMA: ',F8.6)
 3000 FORMAT(/' VIRTUAL COLLISION ESTIMATION:  K-EFFECTIVE AVERAGE: ',
     <        F8.6,' SIGMA: ',F8.6)
      END
