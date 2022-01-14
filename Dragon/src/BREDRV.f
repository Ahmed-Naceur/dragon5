*DECK BREDRV
      SUBROUTINE BREDRV(NC,IPGEO1,IPMAC1,IPGEO2,IPEDI2,NG,LX1,NMIX1,
     > NMIX2,ITRIAL,IMIX1,IGAP,HMREFL,ISPH,LALB,NGET,ADFREF,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the 1D reflector calculation.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NC      number of sn macrolibs (=1: DF-NEM; =2:Lefebvre-Lebigot or
*         Koebke method; >=2 ERM-NEM).
* IPGEO1  nodal geometry.
* IPMAC1  nodal macrolib.
* IPGEO2  sn geometry.
* IPEDI2  sn edition.
* NG      number of energy groups.
* LX1     number of nodes in the nodal calculation.
* NMIX1   number of mixtures in the nodal calculation.
* NMIX2   number of mixtures in the sn calculation after edition.
* ITRIAL  type of expansion functions in the nodal calculation.
*         (=1: polynomial; =2: hyperbolic).
* IMIX1   mix index of node (equal to zero if the node is not used).
* IGAP    mix index of the right gap where the surface flux is
*         recovered (equal to zero if no gap is defined).
* HMREFL  type of reflector model.
* ISPH    SPH flag (=0: use discontinuity factors; =1: use SPH factors).
* LALB    albedo flag (=.TRUE.: compute an equivalent albedo with FD-NEM
*         and ERM-NEM methods).
* NGET    type of NGET normalization if discontinuity factors
*         (=0: simple; =1: imposed ADF on fuel assembly; =2: recover
*         fuel assembly ADF from input macrolib).
* ADFREF  imposed ADF values on fuel assembly side.
* IPRINT  edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NC
      TYPE(C_PTR) IPGEO1,IPMAC1,IPGEO2,IPEDI2(NC)
      INTEGER NG,LX1,NMIX1,NMIX2,ITRIAL(NG),IMIX1(LX1),IGAP(LX1),ISPH,
     1 NGET,IPRINT
      REAL ADFREF(NG)
      CHARACTER HMREFL*12
      LOGICAL LALB
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE),NCODE(6),ICODE(6)
      REAL ZCODE(6)
      CHARACTER HSMG*131,HCASE*12
      LOGICAL LREFL
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IHOM,IMIX,IMIXS,ISPLTX,
     1 ISTOP
      REAL, ALLOCATABLE, DIMENSION(:) :: XXX,XXXS,XXX1,ENER,ZKEFF,B2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: VOL1
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: FLX1,DC1,TOT1,CHI1,SIGF1,
     1 JXM,JXP,FHETXM,FHETXP,ADF1
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SCAT1
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPMAC2
*----
*  RECOVER SN MACROLIBS
*----
      ALLOCATE(IPMAC2(NC))
      ILEAKS=0
      IDF=0
      DO IC=1,NC
        CALL LCMGET(IPEDI2(IC),'STATE-VECTOR',ISTATE)
        IF(IC.EQ.1) THEN
          LXMS=ISTATE(17)
          ALLOCATE(IHOM(LXMS),ENER(NG+1))
        ELSE
          IF(ISTATE(17).NE.LXMS) CALL XABORT('BREDRV: INVALID LXMS.')
        ENDIF
        CALL LCMGET(IPEDI2(IC),'REF:MATCOD',IHOM)
        NMIX_SN=MAXVAL(IHOM)
        CALL LCMGET(IPEDI2(IC),'REF:IMERGE',IHOM)
        CALL LCMGTC(IPEDI2(IC),'LAST-EDIT',12,1,HCASE)
        IPMAC2(IC)=LCMGID(IPEDI2(IC),HCASE)
        IPMAC2(IC)=LCMGID(IPMAC2(IC),'MACROLIB')
        CALL LCMGET(IPMAC2(IC),'STATE-VECTOR',ISTATE)
        IF(ISTATE(1).NE.NG) THEN
          CALL XABORT('BREDRV: INVALID NUMBER OF CONDENSED GROUPS.')
        ELSE IF(ISTATE(2).NE.NMIX2) THEN
          CALL XABORT('BREDRV: INVALID NUMBER OF SN MIXTURES.')
        ELSE IF(ISTATE(4).NE.1) THEN
          CALL XABORT('BREDRV: ONE FISSILE ISOTOPE EXPECTED.')
        ENDIF
        IF(IC.EQ.1) THEN
          ILEAKS=ISTATE(9)
          IDF=ISTATE(12)
          CALL LCMGET(IPMAC2(1),'ENERGY',ENER)
        ENDIF
      ENDDO
      IF((NGET.EQ.2).AND.(IDF.NE.2)) THEN
        CALL XABORT('BREDRV: MISSING ADF INFO IN INPUT MACROLIB.')
      ENDIF
*----
*  RECOVER AND PROCESS GEOMETRY INFORMATION
*----
      CALL LCMGET(IPGEO2,'STATE-VECTOR',ISTATE)
      CALL LCMGET(IPGEO2,'NCODE',NCODE)
      CALL LCMGET(IPGEO2,'ICODE',ICODE)
      CALL LCMGET(IPGEO2,'ZCODE',ZCODE)
      LREFL=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.6).AND.(ZCODE(1).EQ.1.0))
      IF(ISTATE(1).NE.2) THEN
        CALL XABORT('BREDRV: 1D SLAB GEOMETRY EXPECTED.')
      ELSE IF(ISTATE(8).NE.0) THEN
        CALL XABORT('BREDRV: CELL OPTION IS FORBIDDEN.')
      ELSE IF(ISTATE(7).NE.NMIX_SN) THEN
        WRITE(HSMG,'(40H BREDRV: INVALID NUMBER OF SN MIX (GEOM=,I5,
     1  7H MACRO=,I5,2H).)') ISTATE(7),NMIX_SN
        CALL XABORT(HSMG)
      ELSE IF(.NOT.LREFL) THEN
        CALL XABORT('BREDRV: REFLEXION MANDATORY AT LEFT BOUNDARY.')
      ENDIF
      LX=ISTATE(3)
      ALLOCATE(IMIX(LX),XXX(LX+1),ISPLTX(LX))
      CALL LCMGET(IPGEO2,'MIX',IMIX)
      CALL LCMGET(IPGEO2,'MESHX',XXX)
      CALL LCMLEN(IPGEO2,'SPLITX',ILEN2,ITYLCM)
      IF(ILEN2.GT.0) THEN
        CALL LCMGET(IPGEO2,'SPLITX',ISPLTX)
      ELSE
        ISPLTX(:LX)=1
      ENDIF
      LXS=0
      DO I=1,LX
        LXS=LXS+ABS(ISPLTX(I))
      ENDDO
      IF(LXS.NE.LXMS) THEN
        WRITE(HSMG,'(41H BREDRV: INVALID NUMBER OF REGIONS (GEOM=,I5,
     1  9H EDITION=,I5,2H).)') LXS,LXMS
        CALL XABORT(HSMG)
      ENDIF
      ALLOCATE(IMIXS(LXS),XXXS(LXS+1))
      IF(NCODE(2).EQ.5) THEN
        DEL=XXX(LX+1)-XXX(LX)
        IF(MOD(ISPLTX(LX),2).EQ.0) THEN
          ISPLTX(LX)=ISPLTX(LX)/2
          NCODE(2)=2
          XXX(LX+1)=XXX(LX)+REAL(0.5*DEL)
        ELSE
          IGAR=ISPLTX(LX)
          ISPLTX(LX)=(ISPLTX(LX)+1)/2
          XXX(LX+1)=XXX(LX)+REAL(DEL*(DBLE(ISPLTX(LX))/DBLE(IGAR)))
        ENDIF
      ENDIF
      K=LXS+1
      GAR=XXX(LX+1)
      DO IOLD=LX,1,-1
        ISP=ISPLTX(IOLD)
        DEL=(GAR-XXX(IOLD))/REAL(ISP)
        GAR=XXX(IOLD)
        DO I=ABS(ISP),1,-1
          XXXS(K)=REAL(GAR+DEL*DBLE(I))
          K=K-1
          IMIXS(K)=IMIX(IOLD)
        ENDDO
      ENDDO
      XXXS(1)=XXX(1)
      DEALLOCATE(ISPLTX,XXX,IMIX)
      ALLOCATE(XXX1(LX1+1),ISTOP(LX1))
      ISTOP(:LX1)=0
      DO I=1,LXS
        DO J=1,LX1
          IF(IHOM(I).EQ.IGAP(J)) THEN
            IF((ISTOP(J).NE.0).AND.(IGAP(J).NE.0)) THEN
              WRITE(HSMG,'(23H BREDRV: GAP WITH INDEX,I5,10H IS DEFINE,
     1        8HD TWICE.)') IGAP(J)
              CALL XABORT(HSMG)
            ENDIF
            ISTOP(J)=I
          ENDIF
        ENDDO
      ENDDO
      IF(IPRINT.GE.0) THEN
        WRITE(6,'(/20H BREDRV: SN GEOMETRY)')
        WRITE(6,'(1P,10E12.4)') XXXS(:LXS+1)
      ENDIF
      XXX1(1)=XXXS(1)
      XXX1(LX1+1)=XXXS(LXS+1)
      IOF=0
      DO J=1,LX1
        IOF=IOF+1
        IF((IMIX1(J).NE.0) .AND.(IMIX1(J).NE.IOF)) THEN
          CALL XABORT('BREDRV: INCONSISTENT MIX VALUE.')
        ENDIF
        IF(ISTOP(J).GT.0) THEN
          IOF=IOF+1
          XXX1(J+1)=XXXS(ISTOP(J))
        ENDIF
      ENDDO
      NCODE(1)=2
      ICODE(:2)=0
      IF(LALB.AND.(IGAP(LX1).NE.0)) THEN
        NCODE(2)=6
        ICODE(2)=1
      ENDIF
      DEALLOCATE(ISTOP)
      ISTATE(:)=0
      ISTATE(1)=2
      ISTATE(3)=LX1
      ISTATE(6)=LX1
      ISTATE(7)=NMIX1
      CALL LCMPUT(IPGEO1,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPGEO1,'NCODE',6,1,NCODE)
      CALL LCMPUT(IPGEO1,'ICODE',6,1,ICODE)
      CALL LCMPUT(IPGEO1,'ZCODE',6,2,ZCODE)
      ALLOCATE(IMIX(LX1))
      IMIX(:)=0
      IOF=0
      DO J=1,LX1
        IF(IMIX1(J).NE.0) THEN
          IOF=IOF+1
          IMIX(J)=IOF
        ENDIF
      ENDDO
      CALL LCMPUT(IPGEO1,'MIX',LX1,1,IMIX)
      CALL LCMPUT(IPGEO1,'MESHX',LX1+1,2,XXX1)
      IF(IPRINT.GE.0) THEN
        WRITE(6,'(/23H BREDRV: NODAL GEOMETRY)')
        WRITE(6,'(1P,10E12.4)') XXX1(:LX1+1)
      ENDIF
*----
*  COMPUTE MACROSCOPIC CROSS SECTIONS AND SURFACIC DATA
*----
      ALLOCATE(VOL1(NMIX1,NC),FLX1(NMIX1,NG,NC),DC1(NMIX1,NG,NC),
     1 TOT1(NMIX1,NG,NC),CHI1(NMIX1,NG,NC),SIGF1(NMIX1,NG,NC),
     2 SCAT1(NMIX1,NG,NG,NC),JXM(NMIX1,NG,NC),JXP(NMIX1,NG,NC),
     3 FHETXM(NMIX1,NG,NC),FHETXP(NMIX1,NG,NC),ADF1(NMIX1,NG,NC),
     4 ZKEFF(NC),B2(NC))
      CALL BREMAC(NC,IPMAC2,NG,LX1,NMIX1,NMIX2,IMIX,IMIX1,IGAP,ILEAKS,
     1 IDF,IPRINT,ZKEFF,B2,VOL1,FLX1,DC1,TOT1,CHI1,SIGF1,SCAT1,JXM,JXP,
     2 FHETXM,FHETXP,ADF1)
*----
*  SELECT A REFLECTOR MODEL
*----
      IF(HMREFL.EQ."DF-NEM") THEN
        IF(NC.NE.1) CALL XABORT('BREDRV: NC=1 EXPECTED.')
        CALL BRENEM(IPMAC1,NG,LX1,NMIX1,ITRIAL,IMIX,ICODE,ISPH,ZKEFF,
     1  B2,ENER,VOL1,FLX1,DC1,TOT1,CHI1,SIGF1,SCAT1,JXM,JXP,FHETXM,
     2  FHETXP,ADF1,NGET,ADFREF,IPRINT)
      ELSE IF(HMREFL.EQ."ERM-NEM") THEN
        CALL BREERM(IPMAC1,NC,NG,LX1,NMIX1,ITRIAL,IMIX,ICODE,ISPH,ZKEFF,
     1  B2,ENER,VOL1,FLX1,DC1,TOT1,CHI1,SIGF1,SCAT1,JXM,JXP,FHETXM,
     2  FHETXP,ADF1,NGET,ADFREF,IPRINT)
      ELSE IF(HMREFL.EQ."LEFEBVRE-LEB") THEN
        IF(NC.NE.2) CALL XABORT('BREDRV: NC=2 EXPECTED.')
        CALL BRELLB(IPMAC1,NC,NG,NMIX1,ENER,JXM,FHETXM,IPRINT)
      ELSE IF(HMREFL.EQ."KOEBKE") THEN
        IF(NC.NE.2) CALL XABORT('BREDRV: NC=2 EXPECTED.')
        CALL BREKOE(IPMAC1,NC,NG,NMIX1,ISPH,B2,ENER,DC1,TOT1,SCAT1,JXM,
     1  FHETXM,IPRINT)
      ELSE
        WRITE(HSMG,'(25H BREDRV: REFLECTOR MODEL ,A,12H IS UNKNOWN.)')
     1  HMREFL
        CALL XABORT(HSMG)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IMIX,B2,ZKEFF,ADF1,FHETXP,FHETXM,JXP,JXM,SCAT1,SIGF1,
     1 CHI1,TOT1,DC1,FLX1,VOL1)
      DEALLOCATE(XXX1,XXXS,IMIXS,ENER,IHOM)
      DEALLOCATE(IPMAC2)
      RETURN
      END
