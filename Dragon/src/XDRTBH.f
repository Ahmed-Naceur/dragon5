*DECK XDRTBH
      SUBROUTINE XDRTBH(IPGEOM,IPTRK,IQUA10,IBIHET,IMPX,FRTM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the double-heterogeneity (Bihet) data from the geometry
* object IPGEOM and update the tracking object IPTRK.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK   pointer to the excell tracking (L_TRACK).
* IPGEOM  pointer to the geometry (L_GEOM).
* IQUA10  quadrature parameter for the double heterogeneity option.
* IBIHET  type of double-heterogeneity method: =1 Sanchez-Pomraning
*         model; =2 Hebert model; =3 She-Liu-Shi model (no shadow);
*         =4 She-Liu-Shi model (with shadow).
* IMPX    tracking print level.
* FRTM    minimum volume fraction of the grain in the representative 
*         volume for She-Liu-Shi model.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM,IPTRK
      INTEGER IQUA10,IBIHET,IMPX
      REAL FRTM
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTRAK(NSTATE),ISTATE(NSTATE),IPARAM(8)
      CHARACTER CDOOR*12
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NS,IBI,MAT,IDIL,MIXGR,KEYF1,
     1 KEYF2
      REAL, ALLOCATABLE, DIMENSION(:) :: RS,FRACT,VOLK,VOL
*
      IF(IQUA10.EQ.0) CALL XABORT('XDRTBH: INVALID IQUA10.')
      IF(IBIHET.EQ.0) CALL XABORT('XDRTBH: INVALID IBIHET.')
*
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTRAK)
      NREG2=ISTRAK(1)
      NUN2=ISTRAK(2)
      IR2=ISTRAK(4)
      CALL LCMSIX(IPGEOM,'BIHET',1)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      NG=ISTATE(1)
      NSMAX=ISTATE(2)-1
      ALLOCATE(NS(NG))
      CALL LCMGET(IPGEOM,'NS',NS)
      NSMAX=0
      DO 10 I=1,NG
      NSMAX=MAX(NSMAX,NS(I))
   10 CONTINUE
*
      ALLOCATE(IBI(NREG2))
      ALLOCATE(RS(NG*(1+NSMAX)),FRACT(NG*IR2),VOLK(NG*NSMAX))
*
      MAXPTS=NREG2*(NSMAX+1)*NG
      ALLOCATE(MAT(MAXPTS),VOL(MAXPTS))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMSIX(IPTRK,'BIHET',1)
      CALL LCMPUT(IPTRK,'IBI',NREG2,1,MAT)
      CALL LCMPUT(IPTRK,'VOLUME',NREG2,2,VOL)
      CALL LCMSIX(IPTRK,' ',2)
*
      DO 20 I=1,NREG2
      IBI(I)=MAT(I)
   20 CONTINUE
*----
*  RECOVER DOUBLE-HETEROGENEITY INFORMATION FROM GEOMETRY OBJECT
*----
      ALLOCATE(IDIL(IR2),MIXGR(NG*NSMAX*IR2))
      CALL READBH(MAXPTS,IPGEOM,IR1,IR2,NREG,NREG2,MAT,VOL,NG,NSMAX,
     1 MICRO,NS,IBI,RS,FRACT,VOLK,IMPX,IDIL,MIXGR)
      DEALLOCATE(IBI)
      IF(IMPX.GE.1) THEN
         WRITE (6,'(/" QUADRATURE PARAMETER FOR THE MICRO STRUC",
     1   "TURES =",I2/)') IQUA10
         WRITE (6,'(" TYPE OF DOUBLE HETEROGENEITY MODEL (1/2: ",
     1   "SANCHEZ-POMRANING/HEBERT)=",I2/)') IBIHET
      ENDIF
      CALL LCMSIX(IPGEOM,' ',2)
*----
*  RESET STATE-VECTOR INFORMATION
*----
      IPARAM(1)=IR1
      IPARAM(2)=IR2
      IPARAM(3)=NREG2
      IPARAM(4)=NG
      IPARAM(5)=NSMAX
      IPARAM(6)=IBIHET
      IPARAM(7)=MICRO
      IPARAM(8)=IQUA10
      CALL LCMSIX(IPTRK,'BIHET',1)
      CALL LCMPUT(IPTRK,'PARAM',8,1,IPARAM)
      CALL LCMPUT(IPTRK,'NS',NG,1,NS)
      CALL LCMPUT(IPTRK,'RS',NG*(1+NSMAX),2,RS)
      CALL LCMPUT(IPTRK,'FRACT',NG*IR2,2,FRACT)
      CALL LCMPUT(IPTRK,'VOLK',NG*NSMAX,2,VOLK)
      CALL LCMPUT(IPTRK,'IDIL',IR2-IR1,1,IDIL)
      CALL LCMPUT(IPTRK,'MIXGR',NG*NSMAX*(IR2-IR1),1,MIXGR)
      CALL LCMPUT(IPTRK,'FRTM',1,2,FRTM)
      DEALLOCATE(MIXGR,IDIL,NS)
      DEALLOCATE(VOLK,FRACT,RS)
      CALL LCMSIX(IPTRK,' ',2)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTRAK)
      ISTRAK(1)=NREG
      ISTRAK(2)=NUN2+(NREG-NREG2)
      ISTRAK(4)=IR1
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTRAK)
      CALL LCMPUT(IPTRK,'MATCOD',NREG,1,MAT)
      CALL LCMPUT(IPTRK,'VOLUME',NREG,2,VOL)
      DEALLOCATE(VOL,MAT)
*----
*  RESET KEYFLX AND KEYFLX$ANIS
*----
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,CDOOR)
      IF((CDOOR.EQ.'MCCG').OR.(CDOOR.EQ.'SN')) THEN
         CALL LCMLEN(IPTRK,'KEYFLX$ANIS',LKFL,ITYLCM)
         NFUNL=LKFL/NREG2
         ALLOCATE(KEYF1(NREG*NFUNL),KEYF2(NREG2*NFUNL))
         CALL LCMGET(IPTRK,'KEYFLX',KEYF2)
         CALL XDISET(KEYF1,NREG*NFUNL,0)
         DO 35 INF=1,NFUNL
         DO 30 I=1,NREG2
         IOF1=(INF-1)*NREG+I
         IOF2=(INF-1)*NREG2+I
         KEYF1(IOF1)=KEYF2(IOF2)
   30    CONTINUE
   35    CONTINUE
         IUNK=NUN2
         DO 40 I=NREG2+1,NREG
         IUNK=IUNK+1
         KEYF1(I)=IUNK
   40    CONTINUE
         CALL LCMPUT(IPTRK,'KEYFLX',NREG,1,KEYF1(:NREG))
         CALL LCMPUT(IPTRK,'KEYFLX$ANIS',NREG*NFUNL,1,KEYF1)
         DEALLOCATE(KEYF2,KEYF1)
      ELSE
        ALLOCATE(KEYF1(NREG))
        CALL XDISET(KEYF1,NREG,0)
        CALL LCMGET(IPTRK,'KEYFLX',KEYF1(:NREG2))
        IUNK=NUN2
        DO 50 I=NREG2+1,NREG
        IUNK=IUNK+1
        KEYF1(I)=IUNK
   50    CONTINUE
        CALL LCMPUT(IPTRK,'KEYFLX',NREG,1,KEYF1)
        DEALLOCATE(KEYF1)
      ENDIF
      RETURN
      END
