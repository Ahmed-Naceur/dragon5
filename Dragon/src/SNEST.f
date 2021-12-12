*DECK SNEST
      SUBROUTINE SNEST (IPTRK,IMPX,NREG,NUN,MAT,IG,KEYFLX,KEYSPN,
     1   FUNKNO) 
*
*-----------------------------------------------------------------------
*
*Purpose:
* Rearrange SPn flux in the Sn order so that SPn can be used to 
* initialise Sn calculation. Use SPn flux to obtain rough estimate
* of boundary fluxes.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* IPTRK   pointer to the tracking (L_TRACK signature).
* IMPX    print level.
* NREG    total number of regions for which specific values of the
*         neutron flux and reactions rates are required.
* NUN     total number of unknowns in vectors SUNKNO and FUNKNO.
* MAT     index-number of the mixture type assigned to each volume.
* IG      group number.
* KEYFLX  position of averaged flux elements in FUNKNO vector.
* KEYSPN  position of SPn unknowns in FUNKNO vector.
*
*Parameters: input/output
* FUNKNO  SPn (in) / SN (out) unknown vector.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER     IMPX,NREG,NUN,MAT(NREG),IG,KEYFLX(NREG),KEYSPN(NREG) 
      REAL        FUNKNO(NUN)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IUNOUT=6,NSTATE=40)
      INTEGER ISTRK(NSTATE),NLOZH,SPLTL,SBMSH,REM,ISPLH
      REAL ZCODE(6)
*
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMPKEY,ORIKEY
      REAL, ALLOCATABLE, DIMENSION(:) :: FUNSPN
*
*----
*  RECOVER TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'ZCODE',ZCODE)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTRK)
      ITYPE=ISTRK(6)
      NSCT=ISTRK(7)
      IELEM=ISTRK(8)
      NDIM=ISTRK(9)
      LL4=ISTRK(11)
      LX=ISTRK(12)
      LY=ISTRK(13)
      LZ=ISTRK(14)
      ISPLH=1
      NHEX=1
      IF((ITYPE.EQ.8).OR.(ITYPE.EQ.9)) THEN
         ISPLH=ISTRK(26)
         NHEX =LX/(3*ISPLH**2)
      ENDIF
*----
*  PARAMETER VALIDATION
*---- 
      IF((KEYSPN(NREG)).GT.LL4)THEN
         CALL XABORT('SNEST: MORE SPN UNKNOWNS THAN SN UNKNOWNS. '
     1   //'CANNOT GUARANTEE INTEGRITY OF IMPORTED SOLUTION. '
     2   //'CONSIDER INCREASING SPATIAL ORDER FOR SN OR DECREASING '
     3   //'SPATIAL ORDER FOR SPN/DIFF.')
      ENDIF
      IF((ITYPE.NE.2).AND.(ITYPE.NE.5).AND.(ITYPE.NE.7).AND.
     1   (ITYPE.NE.8).AND.(ITYPE.NE.9))CALL XABORT('SNEST: TYPE '
     2   //'OF DISCRETIZATION NOT IMPLEMENTED.')
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(FUNSPN(NUN))
*----
*  PRINT IMPORTED FLUX
*----
      IF(IMPX.GT.5) THEN
         WRITE(IUNOUT,'(//33H I M P O R T E D    F L U X E S (,I5,
     1   3H ):)') IG
         WRITE(IUNOUT,'(1P,4(5X,E15.7))') (FUNKNO(I),I=1,NUN)
      ENDIF
*
*----
*  REBUILD KEYFLX FOR 2D HEXAGONAL CASE
*----
      ! NLOZH - num. of loz. per hexagon
      ! SBMSH - num. of submeshes per lozenge (integer)
      ! SPLTL - split of the lozenge (ISPLH)
      IF((ITYPE.EQ.8).OR.(ITYPE.EQ.9))THEN
         ALLOCATE(TMPKEY(NREG),ORIKEY(NREG))
         TMPKEY(:) = 0
         ORIKEY(1:NREG) = KEYFLX(1:NREG)
         IND = 0
         JND = 0
         NLOZH  = 3*ISPLH**2
         SBMSH  = NLOZH/3
         SPLTL  = ISPLH
         DO IZ=1,LZ
         DO IH=1,NHEX
            DO IM=1,SBMSH
               REM=MOD(IM-1,SPLTL)
               IF((REM.EQ.0).AND.(SBMSH.NE.1))THEN
                  JND = (IH-1)*NLOZH + SBMSH - (IM/SPLTL)
                  JND = JND + (IZ-1)*LX
               ELSEIF((REM.NE.0).AND.(SBMSH.NE.1))THEN
                  JND = JND - (SBMSH*3) - SPLTL
               ENDIF
               DO ILZ=1,3
                  IND = (IZ-1)*LX +(IH-1)*NLOZH +(IM-1)*3 +(ILZ-1) +1
                  IF(SBMSH.EQ.1) JND = IND
                  TMPKEY(IND) = KEYFLX(JND)
                  JND = JND + SBMSH
               ENDDO
            ENDDO
         ENDDO
         ENDDO
         KEYFLX(:) = TMPKEY(:)
         DEALLOCATE(TMPKEY)
      ENDIF
*
*----
*  REARRANGE SPN FLUX IN SN ORDER, FOR P0 ISOTROPIC FLUX ONLY
*----
      CALL XDRSET(FUNSPN,NUN,0.0)
      FUNSPN(1:NUN) = FUNKNO(1:NUN)
      CALL XDRSET(FUNKNO,NUN,0.0)

      DO 100 IR=1,NREG
         IF(MAT(IR).LE.0) GO TO 100
         INDSN=KEYFLX(IR)
         INDPN=KEYSPN(IR)
         FUNKNO(INDSN)=FUNSPN(INDPN)
  100 CONTINUE
*
*----
*  RECUPERATE ORIGINAL KEYFLX FOR HEXAGONAL CASES
*----
      IF((ITYPE.EQ.8).OR.(ITYPE.EQ.9))THEN
         KEYFLX(1:NREG) = ORIKEY(1:NREG)
         DEALLOCATE(ORIKEY)
      ENDIF
*----
*  PRINT REARRANGED FLUX
*----
      IF(IMPX.GT.3) THEN
         WRITE(IUNOUT,'(//37H R E A R R A N G E D    F L U X E S (,I5,
     1   3H ):)') IG
         WRITE(IUNOUT,'(1P,4(5X,E15.7))') (FUNKNO(I),I=1,NUN)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FUNSPN)
      RETURN
      END
