*DECK SOURMONO
      SUBROUTINE SOURMONO(IPTRK,IPGEOM,NREG,LX,LY,LZ,NG,NUNS,NDIM,
     1 NSOUR,ISOUR,ISISOUR,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,XXX,YYY,ZZZ,
     2 MESH_LEN,BSINFO,BS,MAXL,NORM,DIR,MONOP,NBS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute moments of fixed isotropic sources.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Bienvenue 
*
*Parameters: input
*
*Parameters: output
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPGEOM
      INTEGER NREG,NBS,LX,LY,LZ,NG,NUNS,NDIM,NSOUR,MESH_LEN(3),MONOP,
     1 ISISOUR(NG),BSINFO(3,NBS)
      REAL XMIN(NSOUR),XMAX(NSOUR),YMIN(NSOUR),YMAX(NSOUR),ZMIN(NSOUR),
     1 ZMAX(NSOUR),ISOUR(NG),XXX(MESH_LEN(1)),
     2 YYY(MESH_LEN(2)),ZZZ(MESH_LEN(3)),NORM,DIR(3),BS(MAXL,NBS)
*----
*  LOCAL VARIABLES
*----

      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE),SPLIT_LEN(3),XP(NREG),YP(NREG),ZP(NREG),
     1 IPN
      REAL X(LX),Y(LY),Z(LZ),IPVAL
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) U_PTR,DU_PTR,DE_PTR,DZ_PTR,W_PTR
      INTEGER,ALLOCATABLE,DIMENSION(:) :: SPLITX,SPLITY,SPLITZ,MN
      REAL, POINTER, DIMENSION(:) :: U,DU,DE,DZ,W
      REAL,ALLOCATABLE,DIMENSION(:) :: P
*----
*  RECOVER TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NLF=ISTATE(15)

      IF(NDIM.EQ.1) THEN
      CALL LCMGPD(IPTRK,'U',U_PTR)
      CALL LCMGPD(IPTRK,'W',W_PTR)
      CALL C_F_POINTER(U_PTR,U,(/ NLF /))
      CALL C_F_POINTER(W_PTR,W,(/ NLF /))
      ELSE IF(NDIM.EQ.2) THEN
      CALL LCMLEN(IPTRK,'DU',NPQ,ITYLCM)
      CALL LCMGPD(IPTRK,'DU',DU_PTR)
      CALL LCMGPD(IPTRK,'DE',DE_PTR)
      CALL LCMGPD(IPTRK,'W',W_PTR)
      CALL C_F_POINTER(DU_PTR,DU,(/ NPQ /))
      CALL C_F_POINTER(DE_PTR,DE,(/ NPQ /))
      CALL C_F_POINTER(W_PTR,W,(/ NLF /))
      ELSE IF(NDIM.EQ.3) THEN
      CALL LCMLEN(IPTRK,'DU',NPQ,ITYLCM)
      CALL LCMGPD(IPTRK,'DU',DU_PTR)
      CALL LCMGPD(IPTRK,'DE',DE_PTR)
      CALL LCMGPD(IPTRK,'DZ',DZ_PTR)
      CALL LCMGPD(IPTRK,'W',W_PTR)
      CALL C_F_POINTER(DU_PTR,DU,(/ NPQ /))
      CALL C_F_POINTER(DE_PTR,DE,(/ NPQ /))
      CALL C_F_POINTER(DZ_PTR,DZ,(/ NPQ /))
      CALL C_F_POINTER(W_PTR,W,(/ NLF /))
      ENDIF

*----
*  RECOVER GEOMETRY INFORMATION
*----
      IF(NDIM.EQ.2) THEN
      IF(ABS(MONOP).NE.1) THEN     
        CALL LCMLEN(IPGEOM,'SPLITX',SPLIT_LEN(1),ITYLCM)
        ALLOCATE(SPLITX(SPLIT_LEN(1)))
        CALL LCMGET(IPGEOM,'SPLITX',SPLITX)
      ELSE IF(ABS(MONOP).NE.2) THEN
        CALL LCMLEN(IPGEOM,'SPLITY',SPLIT_LEN(2),ITYLCM)
        ALLOCATE(SPLITY(SPLIT_LEN(2)))
        CALL LCMGET(IPGEOM,'SPLITY',SPLITY)
      ENDIF
      ELSE IF(NDIM.EQ.3) THEN
      IF(ABS(MONOP).NE.1.) THEN     
        CALL LCMLEN(IPGEOM,'SPLITX',SPLIT_LEN(1),ITYLCM)
        ALLOCATE(SPLITX(SPLIT_LEN(1)))
        CALL LCMGET(IPGEOM,'SPLITX',SPLITX)
      ELSE IF(ABS(MONOP).NE.2) THEN
        CALL LCMLEN(IPGEOM,'SPLITY',SPLIT_LEN(2),ITYLCM)
        ALLOCATE(SPLITY(SPLIT_LEN(2)))
        CALL LCMGET(IPGEOM,'SPLITY',SPLITY)
      ELSE IF(ABS(MONOP).NE.3) THEN
        CALL LCMLEN(IPGEOM,'SPLITZ',SPLIT_LEN(3),ITYLCM)
        ALLOCATE(SPLITZ(SPLIT_LEN(3)))
        CALL LCMGET(IPGEOM,'SPLITZ',SPLITZ)
      ENDIF 
      ENDIF

      ! CALCULATE X- OR Y- COORDINATES OF EACH VOXELS
      IF(ABS(MONOP).NE.1) THEN 
      K=1
      DO I=1,SPLIT_LEN(1)
      DO J=1,SPLITX(I)
      XP(K)=I
      K=K+1
      ENDDO
      ENDDO
      DO 10 IX=1,LX
      STEPX=(XXX(XP(IX)+1)-XXX(XP(IX)))/SPLITX(XP(IX))
      IF(XP(IX).EQ.1) THEN
        X(IX)=XXX(XP(IX))+0.5*STEPX+STEPX*(IX-1)
      ELSE
        X(IX)=XXX(XP(IX))+0.5*STEPX+STEPX*(IX-SUM(SPLITX(1:XP(IX)-1))-1)
      ENDIF
   10 CONTINUE
      ELSEIF(ABS(MONOP).NE.2.AND.NDIM.GE.2) THEN
      K=1
      DO I=1,SPLIT_LEN(2)
      DO J=1,SPLITY(I)
      YP(K)=I
      K=K+1
      ENDDO
      ENDDO
      DO 110 IY=1,LY
      STEPY=(YYY(YP(IY)+1)-YYY(YP(IY)))/SPLITY(YP(IY))
      IF(YP(IY).EQ.1) THEN
        Y(IY)=YYY(YP(IY))+0.5*STEPY+STEPY*(IY-1)
      ELSE
        Y(IY)=YYY(YP(IY))+0.5*STEPY+STEPY*(IY-SUM(SPLITY(1:YP(IY)-1))-1)
      ENDIF
  110 CONTINUE
      ELSEIF(ABS(MONOP).NE.3.AND.NDIM.GE.3) THEN
      K=1
      DO I=1,SPLIT_LEN(3)
      DO J=1,SPLITZ(I)
      ZP(K)=I
      K=K+1
      ENDDO
      ENDDO
      DO 220 IZ=1,LZ   
      STEPZ=(ZZZ(ZP(IZ)+1)-ZZZ(ZP(IZ)))/SPLITZ(ZP(IZ))
      IF(ZP(IZ).EQ.1) THEN
        Z(IZ)=ZZZ(ZP(IZ))+0.5*STEPZ+STEPZ*(IZ-1)
      ELSE
        Z(IZ)=ZZZ(ZP(IZ))+0.5*STEPZ+STEPZ*(IZ-SUM(SPLITZ(1:ZP(IZ)-1))-1)
      ENDIF        
  220 CONTINUE
      ENDIF


*----
*  1D CARTESIAN CASE
*----

      IF(NDIM.EQ.1) THEN

      ! Define the corresponding discrete ordinate  
      IPN=1
      IPVAL=ABS(U(1)-DIR(1))
      DO IP=2,NLF
        IF(ABS(U(IP)-DIR(1)).LT.IPVAL) THEN
        IPN=IP
        IPVAL=(U(IP)-DIR(1))**2
        ENDIF
      ENDDO

      ! Normalization
      NORM=SUM(ISOUR)
 
      ! Save source information
      IND=1
      DO IG=1,NG
        IF(ISISOUR(IG).EQ.1) THEN
        BSINFO(1,IND)=IG
        BSINFO(2,IND)=IPN
        BSINFO(3,IND)=MONOP
        BS(1,IND)=ISOUR(IG)/W(IPN)
        IND=IND+1
        ENDIF
      ENDDO 
      
*----
*  2D CARTESIAN CASE
*----

      ELSE IF(NDIM.EQ.2) THEN

      ! Define the corresponding discrete ordinate  
      ALLOCATE(MN(2),P(2))
      MN(1)=1
      UVAL1=ABS(DU(1)-DIR(1))
      EVAL1=ABS(DE(1)-DIR(2))
      MN(2)=1
      UVAL2=UVAL1
      EVAL2=EVAL1
      DO M=2,NPQ
        IF((ABS(DU(M)-DIR(1))+ABS(DE(M)-DIR(2))).LT.UVAL2+EVAL2) THEN
        IF((ABS(DU(M)-DIR(1))+ABS(DE(M)-DIR(2))).LT.UVAL1+EVAL1) THEN
        MN(2)=MN(1)
        UVAL2=UVAL1
        EVAL2=EVAL1
        MN(1)=M
        UVAL1=ABS(DU(M)-DIR(1))
        EVAL1=ABS(DE(M)-DIR(2))
        ELSE
        MN(2)=M
        UVAL2=ABS(DU(M)-DIR(1))
        EVAL2=ABS(DE(M)-DIR(2))
        ENDIF
        ELSE IF(M.EQ.2) THEN
        MN(2)=M
        UVAL2=ABS(DU(M)-DIR(1))
        EVAL2=ABS(DE(M)-DIR(2))
        ENDIF
      ENDDO

      VALTOT=(UVAL1+EVAL1)**2+(UVAL2+EVAL2)**2
      P(1)=1.0-(UVAL1+EVAL1)**2/VALTOT
      P(2)=1.0-(UVAL2+EVAL2)**2/VALTOT
      PRINT *, P(1),P(2),MN(1),MN(2),'TRALALA'

      ! Normalization
      NORM=SUM(ISOUR)

      ! Save source information
      IND=1
      DO IDIR=1,2
      DO IG=1,NG
        IF(ISISOUR(IG).EQ.1) THEN
        BSINFO(1,IND)=IG
        BSINFO(2,IND)=MN(IDIR)
        BSINFO(3,IND)=MONOP
        IF(ABS(MONOP).NE.1) THEN
          DO IX=1,LX
          DO NS=1,NSOUR
          IF(XMIN(NS).LE.X(IX).AND.XMAX(NS).GE.X(IX)) THEN
            BS(IX,IND)=ISOUR(IG)/W(MN(IDIR))*P(IDIR)
          ELSE
            BS(IX,IND)=0.0
          ENDIF
          ENDDO
          ENDDO
        ELSEIF(ABS(MONOP).NE.2) THEN
          DO IY=1,LY
          DO NS=1,NSOUR
          IF(YMIN(NS).LE.Y(IY).AND.YMAX(NS).GE.Y(IY)) THEN
            BS(IY,IND)=ISOUR(IG)/W(MN(IDIR))*P(IDIR)
          ELSE
            BS(IY,IND)=0.0
          ENDIF
          ENDDO
          ENDDO
        ENDIF
        IND=IND+1
        ENDIF
      ENDDO 
      ENDDO  

*----
*  3D CARTESIAN CASE
*----

      ELSE IF(NDIM.EQ.3) THEN

      CALL XABORT('SOUR: 3D CARTESIAN BOUNDARY SOURCES NOT IMPLEMENTED.'
     1 )

      ENDIF

      RETURN
      END
