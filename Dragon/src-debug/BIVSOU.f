*DECK BIVSOU
      SUBROUTINE BIVSOU(MAX1,IG,IPTRK,KPMACR,NANIS,NREG,NMAT,NUNKNO,
     > NGRP,MATCOD,VOL,FLUX,SOURCE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the source for the solution of diffusion or PN equations.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAX1    first dimension of FLUX and SOURCE arrays.
* IG      secondary group.
* IPTRK   pointer to the tracking LCM object.
* KPMACR  pointer to the secondary-group related macrolib information.
* NANIS   maximum cross section Legendre order.
* NREG    number of regions.
* NMAT    number of mixtures.
* NUNKNO  number of unknowns per energy group including spherical
*         harmonic terms, interface currents and fundamental
*         currents.
* NGRP    number of energy groups.
* MATCOD  mixture indices.
* VOL     volumes.
* FLUX    fluxes.
*
*Parameters: output
* SOURCE  sources.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,KPMACR
      INTEGER MAX1,IG,NANIS,NREG,NMAT,NUNKNO,NGRP,MATCOD(NREG)
      REAL VOL(NREG),FLUX(MAX1,NGRP),SOURCE(MAX1,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,PI4=12.5663706144)
      INTEGER JPAR(NSTATE)
      REAL R1DP(5,5)
      CHARACTER CAN(0:9)*2
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS,KN
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT,RR
*----
*  DATA STATEMENTS
*----
      DATA CAN /'00','01','02','03','04','05','06','07','08','09'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(0:NMAT),NJJ(0:NMAT),IPOS(0:NMAT))
      ALLOCATE(XSCAT(0:NMAT*NGRP))
*----
*  RECOVER BIVAC SPECIFIC PARAMETERS.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      IF(JPAR(1).NE.NREG) CALL XABORT('BIVSOU: INCONSISTENT NREG.')
      IF(JPAR(2).NE.NUNKNO) CALL XABORT('BIVSOU: INCONSISTENT NUNKNO.')
      ITYPE=JPAR(6)
      IELEM=JPAR(8)
      ICOL=JPAR(9)
      ISPLH=JPAR(10)
      L4=JPAR(11)
      LX=JPAR(12)
      LY=JPAR(13)
      NLF=JPAR(14)
      ISPN=JPAR(15)
      ISCAT=JPAR(16)
      IF(IELEM.NE.1) CALL XABORT('BIVSOU: DUAL AND LINEAR FLUX APPROXI'
     1 //'MATION IN SPACE EXPECTED.')
*----
*  SPHERICAL HARMONICS (PN) CASE.
*----
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      ALLOCATE(KN(MAXKN))
      CALL LCMGET(IPTRK,'KN',KN)
*----
*  RECOVER THE FINITE ELEMENT UNIT STIFFNESS MATRIX.
*----
      IF(NLF.GT.0) THEN
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(RR(LC*LC))
         CALL LCMGET(IPTRK,'R',RR)
         CALL LCMSIX(IPTRK,' ',2)
         DO 15 IC=1,LC
         DO 10 JC=1,LC
         R1DP(IC,JC)=RR((JC-1)*LC+IC)
   10    CONTINUE
   15    CONTINUE
         DEALLOCATE(RR)
      ENDIF
*----
*  COMPUTE THE SOURCE
*----
      IF(NLF.EQ.0) THEN
*
*  ++++ DIFFUSION THEORY ++++
*
         CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
         CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
         CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
         CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
         IF((ITYPE.EQ.2).OR.(ITYPE.EQ.5)) THEN
*----
*  CARTESIAN 2D DUAL (THOMAS-RAVIART) CASE.
*----
            NUM1=0
            DO 30 IR=1,NREG
            IBM=MATCOD(IR)
            IF(IBM.LE.0) GO TO 30
            IND=KN(NUM1+1)
            JG=IJJ(IBM)
            DO 20 JND=1,NJJ(IBM)
            IF(JG.NE.IG) THEN
               SOURCE(IND,IG)=SOURCE(IND,IG)+FLUX(IND,JG)*
     >         XSCAT(IPOS(IBM)+JND-1)
            ENDIF
            JG=JG-1
   20       CONTINUE
            NUM1=NUM1+5
   30       CONTINUE
         ELSE
            CALL XABORT('BIVSOU: DISCRETIZATION NOT IMPLEMENTED(1).')
         ENDIF
      ELSE
*
*  ++++ PN OR SPN THEORY ++++
*
         DO 310 IL=0,MIN(ABS(ISCAT)-1,NANIS)
         FACT=REAL(2*IL+1)
         CALL LCMGET(KPMACR,'NJJS'//CAN(IL),NJJ(1))
         CALL LCMGET(KPMACR,'IJJS'//CAN(IL),IJJ(1))
         CALL LCMGET(KPMACR,'IPOS'//CAN(IL),IPOS(1))
         CALL LCMGET(KPMACR,'SCAT'//CAN(IL),XSCAT(1))
         IF((ITYPE.EQ.2).OR.((ITYPE.EQ.5).AND.(ISPN.EQ.1))) THEN
*----
*  CARTESIAN 2D SPN CASE.
*----
            NUM1=0
            DO 300 IR=1,NREG
            IBM=MATCOD(IR)
            IF(IBM.LE.0) GO TO 300
            IF(MOD(IL,2).EQ.0) THEN
               IND=(IL/2)*L4+KN(NUM1+1)
               JG=IJJ(IBM)
               DO 250 JND=1,NJJ(IBM)
               IF(JG.NE.IG) THEN
                  SOURCE(IND,IG)=SOURCE(IND,IG)+FACT*FLUX(IND,JG)*
     >            XSCAT(IPOS(IBM)+JND-1)
               ENDIF
               JG=JG-1
  250          CONTINUE
            ELSE
               DO 275 IC=1,2
               IIC=1+(IC-1)*IELEM
               IND1=(IL/2)*L4+ABS(KN(NUM1+1+IC))
               S1=REAL(SIGN(1,KN(NUM1+1+IC)))
               DO 270 JC=1,2
               JJC=1+(JC-1)*IELEM
               IND2=(IL/2)*L4+ABS(KN(NUM1+1+JC))
               IF((KN(NUM1+1+IC).NE.0).AND.
     1            (KN(NUM1+1+JC).NE.0)) THEN
                  S2=REAL(SIGN(1,KN(NUM1+1+JC)))
                  AUXX=S1*S2*FACT*R1DP(IIC,JJC)*VOL(IR)
                  JG=IJJ(IBM)
                  DO 260 JND=1,NJJ(IBM)
                  IF(JG.NE.IG) THEN
                     SOURCE(IND1,IG)=SOURCE(IND1,IG)-AUXX*FLUX(IND2,JG)*
     1               XSCAT(IPOS(IBM)+JND-1)
                  ENDIF
                  JG=JG-1
  260             CONTINUE
               ENDIF
  270          CONTINUE
  275          CONTINUE
               DO 295 IC=3,4
               IIC=1+(IC-3)*IELEM
               IND1=(IL/2)*L4+ABS(KN(NUM1+1+IC))
               S1=REAL(SIGN(1,KN(NUM1+1+IC)))
               DO 290 JC=3,4
               JJC=1+(JC-3)*IELEM
               IND2=(IL/2)*L4+ABS(KN(NUM1+1+JC))
               IF((KN(NUM1+1+IC).NE.0).AND.
     1            (KN(NUM1+1+JC).NE.0)) THEN
                  S2=REAL(SIGN(1,KN(NUM1+1+JC)))
                  AUXX=S1*S2*FACT*R1DP(IIC,JJC)*VOL(IR)
                  JG=IJJ(IBM)
                  DO 280 JND=1,NJJ(IBM)
                  IF(JG.NE.IG) THEN
                     SOURCE(IND1,IG)=SOURCE(IND1,IG)-AUXX*FLUX(IND2,JG)*
     1               XSCAT(IPOS(IBM)+JND-1)
                  ENDIF
                  JG=JG-1
  280             CONTINUE
               ENDIF
  290          CONTINUE
  295          CONTINUE
            ENDIF
            NUM1=NUM1+5
  300       CONTINUE
         ELSE
            CALL XABORT('BIVSOU: DISCRETIZATION NOT IMPLEMENTED(2).')
         ENDIF
  310    CONTINUE
      ENDIF
      DEALLOCATE(KN)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XSCAT)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
