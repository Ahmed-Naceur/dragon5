*DECK MCTRK
      SUBROUTINE MCTRK(IPTRK,IPRINT,NFREG,NFSUR,NDIM,NMIX,ANGBC,ITYPBC,
     1                 MAXMSH,NUCELL,MXGSUR,MXGREG,MAXPIN,BCRT,ICODE,
     2                 ALBEDO,IUNFLD,DGMESH,XYZL,INDEX,IDREG,DCMESH,
     3                 ITPIN,DRAPIN,ISEED,NGRP,NL,NFM,NDEL,NED,MATCOD,
     4                 XST,XSS,XSN2N,XSN3N,XSSNN,XSNUSI,XSCHI,XSEDI,
     5                 MIX,ISONBR,NU,POS,ITALLY,NBSCO,NMERGE,NGCOND,
     6                 IMERGE,INDGRP,SCORE1,SCORE2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Simulation of a single particle from source to death.
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
* IPRINT  print level.
* NFREG   number of regions.
* NFSUR   number of surfaces.
* NDIM    problem dimensions.
* NMIX    number of mixtures in the geometry.
* ANGBC   angular treatment for boundary conditions (=0 isotropic;
*         =1 specular).
* ITYPBC  type of boundary.
* MAXMSH  maximum number of elements in MESH array.
* NUCELL  number of cell after unfolding in $X$, $Y$ and $Z$ directions.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* MAXPIN  maximum number of pins in a cell.
* BCRT    reflection/translation array.
* ICODE   albedo index array.
* ALBEDO  albedo array.
* IUNFLD  description of unfolded geometry.
* DGMESH  meshing vector for global geometry.
* XYZL    Cartesian boundary coordinates.
* ISEED   the seed for the generation of random numbers.
* NGRP    number of energy groups.
* NL      number of Legendre orders required in the estimations
*         (NL=1 or higher).
* NFM     number of fissile isotopes.
* NDEL    number of delayed precursor groups.
* NED     number of extra edit vectors.
* XST     total macroscopic cross sections for each mixture and energy
*         group.
* XSS     total scattering cross sections for each mixture and energy
*         group.
* XSN2N   N2N macroscopic cross sections for each mixture and energy
*         group.
* XSN3N   N3N macroscopic cross sections for each mixture and energy
*         group.
* MATCOD  region material.
* XSSNN   in-group and out-of-group macroscopic transfert cross sections
*         for each mixture.
* XSNUSI  the values of Nu time the fission cross sections for each
*         isotope per mixture and energy group.
* XSCHI   the values of fission spectrum per isotope per mixture for
*         each energy group.
* XSEDI   extra edit cross sections for each mixture and energy group.
* MIX     the mixture number where the fission occurs.
* ISONBR  the isotopic number where the fission occurs.
* NU      the value of the particle weight.
* POS     location of the particle in the x, y, and z directions.
* ITALLY  type of tally (=0 no tally; =1 score effective
*         multiplication factor; =2 also score macrolib information).
* NBSCO   number of macrolib related scores.
* NMERGE  number of homogenized regions.
* NGCOND  number of condensed energy groups.
* IMERGE  homogenized regions indices.
* INDGRP  condensed groups indices.
* INDEX   undefined.
* IDREG   undefined.
* DCMESH  undefined.
* ITPIN   undefined.
* DRAPIN  undefined.
*
*Parameters: input/output
* SCORE1  score for total flux and effective multiplication factor.
* SCORE2  macrolib score matrix.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER IPRINT,NFREG,NFSUR,NDIM,NMIX,ANGBC,ITYPBC,MAXMSH,
     1 NUCELL(3),MXGSUR,MXGREG,MAXPIN,MIX,ISONBR,
     2 BCRT(NFSUR),ICODE(6),IUNFLD(2,NUCELL(1),NUCELL(2),NUCELL(3)),
     3 INDEX(5,-MXGSUR:MXGREG,2),IDREG(MXGREG,2),ITPIN(3,MAXPIN),
     4 ISEED,NGRP,NL,NFM,NDEL,NED,MATCOD(NFREG),ITALLY,NBSCO,NMERGE,
     5 NGCOND,IMERGE(NFREG),INDGRP(NGRP)
      REAL ALBEDO(6),XST(NMIX,NGRP),XSS(NMIX,NGRP,NL),XSN2N(NMIX,NGRP),
     1 XSN3N(NMIX,NGRP),XSSNN(NMIX,NGRP,NGRP,NL),
     2 XSNUSI(NMIX,NFM,NGRP,1+NDEL),XSCHI(NMIX,NFM,NGRP,1+NDEL),
     3 XSEDI(NMIX,NGRP,NED),NU,SCORE1(3),SCORE2(NBSCO,NMERGE,NGCOND)
      DOUBLE PRECISION DGMESH(-1:MAXMSH,4),XYZL(2,NDIM),
     1 DCMESH(-1:MAXMSH,4,2),DRAPIN(-1:4,MAXPIN),POS(3)
*----
*  LOCAL VARIABLES
*----
      INTEGER NBIND
      PARAMETER (NBIND=7)
      DOUBLE PRECISION XDRCST,PI,TWOPI
      INTEGER INDX(NBIND,0:2),ODIR(3),IREG,IDIRC(2),MESHC(4,2),NSURC(2),
     1 NREGC(2),NTPIN,IDIR,ILEV,IDS,IDSO,IPSP,IEV,II,IGROUP,IGR,
     2 ISUR,ISUR2,IBM,IDIM,IFIRST
      DOUBLE PRECISION CELLPO(3,2),PINCEN(3),LENGTH,MU,MUS,VIS,PHI,
     2 VDIR(3),OMEGA(3),ALB,NUSIGF0
      LOGICAL TRACK,ACTIVE,KILLED,VIRTUAL
      PARAMETER (ACTIVE=.TRUE.,KILLED=.FALSE.) 
      REAL    RAND,XSM
      DOUBLE PRECISION PSCAT,PN2N,PN3N,XPROB,XX
*----
*  DATA VALUES
*----
      INTEGER INDOS(2,3),INDIC(-3:3)
      DATA INDOS / 2,3,
     1             3,1,
     2             1,2 /
      DATA INDIC / 5,3,1,0,2,4,6/
*----
*  LOAD CONSTANTS
*----
      PI=XDRCST('Pi',' ')
      TWOPI=2.D0*PI
      IFIRST=1
*----
*  TRACK THE NEUTRON PATH
*----
      DO ILEV=0,2
         DO IDIR=1,3
           INDX(IDIR,ILEV)=1
         ENDDO
         DO IDIR=4,NBIND
           INDX(IDIR,ILEV)=0
         ENDDO
      ENDDO
      CALL XDISET(ODIR,3,1)

      TRACK=ACTIVE
      LENGTH=0.D0
      IDSO=0
      VIRTUAL=.FALSE.
     
      DO WHILE(TRACK.EQV.ACTIVE)
        IF(IDSO.EQ.0) THEN
          IF(LENGTH.EQ.0.D0) THEN
            IF(MIX.EQ.-1) THEN
              IGR=1
            ELSE
              CALL RANDF(ISEED,IFIRST,RAND)
              XPROB=0.D0
              IGR=0
              DO WHILE((RAND.GE.XPROB).AND.(IGR.LT.NGRP))
                IGR=IGR+1
                XPROB=XPROB+XSCHI(MIX,ISONBR,IGR,1)
              ENDDO
            ENDIF
*----
*  SET THE TRACK DIRECTION FOR THE FIRST FLIGHT
*----
            XSM=0.0
            DO II=1,NMIX
              XSM=MAX(XSM,XST(II,IGR))
            ENDDO
            CALL RANDF(ISEED,IFIRST,RAND)
            MU=2.D0*RAND-1.D0
            MUS=DSQRT(1.D0-MU*MU)
            CALL RANDF(ISEED,IFIRST,RAND)
            PHI=TWOPI*RAND
            VDIR(1)=MUS*COS(PHI)
            VDIR(2)=MUS*SIN(PHI)
            VDIR(3)=MU
          ELSE
*----
*  IF IT IS A VIRTUAL COLLISION, KEEP THE SAME DIRECTION,
*  IF IT IS NOT A VIRTUAL COLLISION IT CONSISTS OF AN ISOTROPIC
*  SCATTERING REACTION.
*----
            IF(.NOT.VIRTUAL) THEN
              CALL RANDF(ISEED,IFIRST,RAND)
              PHI=TWOPI*RAND
              CALL RANDF(ISEED,IFIRST,RAND)
              MU=COS(PI*(2.D0*RAND-1.D0))
              MUS=DSQRT(1.D0-MU*MU)
              VIS=DSQRT(1.D0-VDIR(1)*VDIR(1))
              OMEGA(1)=VDIR(1)*MU-VIS*MUS*COS(PHI)
              OMEGA(2)=VDIR(2)*MU+MUS*(VDIR(1)*VDIR(2)*
     1                 COS(PHI)-VDIR(3)*SIN(PHI))/VIS
              OMEGA(3)=VDIR(3)*MU+MUS*(VDIR(1)*VDIR(3)*
     1                 COS(PHI)+VDIR(2)*SIN(PHI))/VIS
              DO II=1,3
                VDIR(II)=OMEGA(II)
              ENDDO 
            ENDIF
          ENDIF
*----
*  SAMPLE THE FREE PATH DISTANCE
*----
          CALL RANDF(ISEED,IFIRST,RAND)
          LENGTH=-LOG(RAND)/XSM
        ENDIF
*----
*  LOCATE THE NEUTRON IN THE GEOMETRY
*----
        CALL MCTPTR(IPTRK,IPRINT,NDIM,MAXMSH,ITYPBC,NUCELL,MXGSUR,
     1     MXGREG,MAXPIN,IUNFLD,DGMESH,XYZL,NBIND,POS,LENGTH,VDIR,
     2     ODIR,IDS,IDSO,IREG,INDX,IDIRC,MESHC,NSURC,NREGC,NTPIN,
     3     CELLPO,PINCEN,INDEX,IDREG,DCMESH,ITPIN,DRAPIN)
*---
*  TALLY PROCESSING
*---
        IF(ITALLY.GT.0) THEN
          CALL MCTALLY(ITALLY,NFREG,NMIX,NGRP,NL,NFM,NDEL,NED,NBSCO,
     1    NMERGE,NGCOND,IREG,IGR,NU,MATCOD,IMERGE,INDGRP,XSM,XST,XSS,
     2    XSN2N,XSN3N,XSSNN,XSNUSI,XSCHI,XSEDI,SCORE1,SCORE2)
        ENDIF
*----
*  AN INTERACTION HAS BEEN DETECTED, DETECT IF IT IS A VIRTUAL OR A 
*  REAL COLLISION
*----
        IF(IREG.GT.0) THEN
          IBM=MATCOD(IREG)
          CALL RANDF(ISEED,IFIRST,RAND)
          VIRTUAL=(RAND.LE.((XSM-XST(IBM,IGR))/XSM))
*---
*  DETERMINE THE TYPE OF REACTION
*---
          IF(.NOT.VIRTUAL) THEN
            CALL RANDF(ISEED,IFIRST,RAND)
            PSCAT=XSS(IBM,IGR,1)/XST(IBM,IGR)
            PN2N=XSN2N(IBM,IGR)/XST(IBM,IGR)
            PN3N=XSN3N(IBM,IGR)/XST(IBM,IGR)
            IF(RAND.LE.PSCAT+PN2N+PN3N) THEN
*---
*  ISOTROPIC SCATTERING OR NxN EVENT
*---
              IF(RAND.GT.PSCAT+PN2N) THEN
                NU=3.0*NU
              ELSE IF(RAND.GT.PSCAT) THEN
                NU=2.0*NU
              ENDIF
              XX = 0.D0
              CALL RANDF(ISEED,IFIRST,RAND)
              IGROUP=NGRP
              DO WHILE((RAND.GT.XX).AND.(IGROUP.GE.1))
                XX=XX+XSSNN(IBM,IGROUP,IGR,1)/XSS(IBM,IGR,1)
                IGROUP=IGROUP-1
              ENDDO
              IGR=IGROUP+1
              TRACK=ACTIVE
              XSM=0.0
              DO II=1,NMIX
                XSM=MAX(XSM,XST(II,IGR))
              ENDDO
            ELSE
*---
*  CAPTURE OR FISSION EVENT
*---
              TRACK = KILLED
              MIX   = IBM
              NUSIGF0=0.0D0
              DO II=1,NFM
                NUSIGF0=NUSIGF0+XSNUSI(IBM,II,IGR,1)
              ENDDO
              IF((NUSIGF0.GT.0.0).AND.(NFM.EQ.1)) THEN
                NU=NU*XSNUSI(IBM,1,IGR,1)/(XST(IBM,IGR)-XSS(IBM,IGR,1)-
     1          XSN2N(IBM,IGR)-XSN3N(IBM,IGR))
                ISONBR=1
              ELSE IF(NUSIGF0.GT.0.0) THEN
                CALL RANDF(ISEED,IFIRST,RAND)
                XPROB=0.D0
                ISONBR=0
                DO WHILE((RAND.GE.XPROB).AND.(ISONBR.LT.NFM))
                  ISONBR=ISONBR+1
                  XPROB=XPROB+XSNUSI(IBM,ISONBR,IGR,1)/NUSIGF0
                ENDDO
                NU=NU*XSNUSI(IBM,ISONBR,IGR,1)/(XST(IBM,IGR)-
     1          XSS(IBM,IGR,1)-2.0*XSN2N(IBM,IGR)-3.0*XSN3N(IBM,IGR))
              ELSE
                NU=0.0
                ISONBR=0
              ENDIF
            ENDIF
          ENDIF
        ELSE
*----
*  A BOUNDARY CONDITION HAS BEEN ENCOUNTERED
*----
          IF(ITYPBC.EQ.0) THEN
*----
*  CARTESIAN BOUNDARY
*----
            ISUR=-IREG
            ISUR2=BCRT(ISUR)
            IREG=-ISUR2
            IF((ISUR2.EQ.ISUR).AND.(ANGBC.EQ.1)) THEN
*             SPECULAR REFLECTIVE BOUNDARY CONDITION
              IF(IPRINT.GT.4) WRITE(6,*) 'SPECULAR REFLECTION ON ',ISUR
              ALB=ALBEDO(-ICODE(INDIC(IDSO)))
              IF(ALB.EQ.0.0) THEN
*               no leakage
                ISONBR = 0
                TRACK = KILLED
              ELSE IF(ALB.EQ.1.0) THEN
                VDIR(IDS)=-VDIR(IDS)
              ELSE
                CALL RANDF(ISEED,IFIRST,RAND)
                TRACK=(RAND.LE.ALB)
                VDIR(IDS)=-VDIR(IDS)
              ENDIF
            ELSE IF(ISUR2.NE.ISUR) THEN
*             PERIODIC BOUNDARY CONDITION
              IF(IPRINT.GT.4) WRITE(6,*) 'BC TRANSLATION FROM ',ISUR,
     1        ' TO',ISUR2
              IPSP=-ISUR
              IEV=-INDIC(IDSO)
              IF(IPRINT.GT.99) CALL MCTPSP(IPTRK,POS,IPSP,IEV)
              IF(IDSO.GT.0) THEN
                POS(IDS)=POS(IDS)-(XYZL(2,IDS)-XYZL(1,IDS))
              ELSE
                POS(IDS)=POS(IDS)+(XYZL(2,IDS)-XYZL(1,IDS))
              ENDIF
              IDSO=-IDSO
            ELSE IF((ISUR2.EQ.ISUR).AND.(ANGBC.EQ.0)) THEN
*             WHITE BOUNDARY CONDITION
              ALB=ALBEDO(-ICODE(INDIC(IDSO)))
              IF(ALB.EQ.0.0) THEN
*               no leakage
                ISONBR = 0
                TRACK = KILLED
              ELSE
*               otherwise, choose randomly according to albedo
                CALL RANDF(ISEED,IFIRST,RAND)
                TRACK=(RAND.LE.ALB)
                IF(TRACK.EQV.ACTIVE) THEN
*                 FOR ISOTROPIC BC, CHOOSE RANDOMLY THE REENTERING
*                 DIRECTION
                  CALL RANDF(ISEED,IFIRST,RAND)
                  MU=2.D0*RAND-1.D0
                  MUS=DSQRT(1.D0-MU*MU)
                  CALL RANDF(ISEED,IFIRST,RAND)
                  PHI=PI*RAND
                  VDIR(IDS)=-SIGN(1,IDSO)*MUS*SIN(PHI)
                  VDIR(INDOS(1,IDS))=MUS*COS(PHI)
                  VDIR(INDOS(2,IDS))=MU
                  DO IDIM=1,NDIM
                    IF(IDIM.NE.IDS) THEN
                      CALL RANDF(ISEED,IFIRST,RAND)
                      POS(IDIM)=XYZL(1,IDIM)+RAND*(XYZL(2,IDIM)-
     1                XYZL(1,IDIM))
                    ENDIF
                  ENDDO
                  CALL RANDF(ISEED,IFIRST,RAND)
                  LENGTH=-LOG(RAND)/XSM
                ELSE
                  ISONBR = 0
                ENDIF
              ENDIF
            ELSE
              CALL XABORT('MCTRK: INVALID TYPE OF BOUNDARY CONDITION.')
            ENDIF
          ELSE
*----
*  CYLINDRICAL BOUNDARY
*----
            CALL XABORT('MCTRK: CYLINDRICAL BOUNDARY NOT IMPLEMENTED.')
          ENDIF
          IF(IPRINT.GT.99) THEN
             IPSP=-IREG
             IEV=-INDIC(IDSO)
             CALL MCTPSP(IPTRK,POS,IPSP,IEV)
          ENDIF
        ENDIF
*----
*  SAVE NEUTRON PATHS IN NXT TABLE FOR PSP DISPLAY
*----
        IF(IPRINT.GT.99) THEN
           IF(TRACK.EQV.KILLED) THEN
             IPSP=-ABS(IREG)
           ELSE
             IPSP=ABS(IREG)
           ENDIF
           CALL MCTPSP(IPTRK,POS,IPSP,1)
        ENDIF
      ENDDO
      RETURN
      END
