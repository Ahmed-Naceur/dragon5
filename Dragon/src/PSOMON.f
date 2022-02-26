*DECK PSOMON
      SUBROUTINE PSOMON(IPTRK,IPGEOM,NREG,LX,LY,LZ,NG,NDIM,
     1 NSOUR,ISOUR,ISISOUR,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,XXX,YYY,ZZZ,
     2 MESHL,BSINFO,BS,MAXL,NORM,DIR,MONOP,NBST,NBS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute moments of fixed boundary monodirectionnal sources.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Bienvenue 
*
*Parameters: input
* IPTRK   pointer to the tracking LCM object.
* IPGEOM  porinter to the geometry LCM object.
* NREG    number of regions.
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* LZ      number of meshes along Z axis.
* NG      number of energy groups.
* NDIM    geometry dimension.
* NSOUR   number of sources defined.
* ISOUR   intensity of the sources. 
* ISISOUR array with 0/1 values to indicate if energy group contain
*         sources
* XMIN    lower boundaries of the source along X axis.
* XMAX    upper boundaries of the source along X axis.
* YMIN    lower boundaries of the source along Y axis.
* YMAX    upper boundaries of the source along Y axis.
* ZMIN    lower boundaries of the source along Z axis.
* ZMAX    upper boundaries of the source along Z axis.
* XXX     regions boundaries along X axis. 
* YYY     regions boundaries along Y axis.
* ZZZ     regions boundaries along Z axis.
* MESHL   number of regions along X-, Y- and Z-axis
* MAXL    number of intensity values needed to fully described each 
*         boundary source.
* DIR     direction of the source particles.
* MONOP   value describing the boundary source location.
* NBST    total number of sources.
*
*Parameters: output
* NORM    normalization factor.
* BS      intensity values for each sources.
* BSINFO  energy group, discrete ordinate and location corresponding to
*         each sources.
* NBS     number of sources in each energy group.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPGEOM
      INTEGER NREG,NBST,LX,LY,LZ,NG,NDIM,NSOUR,MESHL(3),MONOP,
     1 ISISOUR(NG),BSINFO(2,NG,NBST),NBS(NG)
      REAL XMIN(NSOUR),XMAX(NSOUR),YMIN(NSOUR),YMAX(NSOUR),ZMIN(NSOUR),
     1 ZMAX(NSOUR),ISOUR(NG),XXX(MESHL(1)),
     2 YYY(MESHL(2)),ZZZ(MESHL(3)),NORM,DIR(3),BS(MAXL,NG,NBST)
*----
*  LOCAL VARIABLES
*----
      INTEGER SPLIT_LEN(3),XP(NREG),YP(NREG),ZP(NREG),
     1 IPN,MONOP2
      REAL X(LX),Y(LY),Z(LZ),IPVAL,S
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) U_PTR,DU_PTR,DE_PTR,DZ_PTR,W_PTR
      INTEGER,ALLOCATABLE,DIMENSION(:) :: SPLITX,SPLITY,SPLITZ,MN
      REAL, POINTER, DIMENSION(:) :: U,DU,DE,DZ,W
      REAL,ALLOCATABLE,DIMENSION(:) :: P
*----
*  RECOVER GEOMETRY INFORMATION
*----
      IF(ABS(MONOP).NE.1) THEN     
        CALL LCMLEN(IPGEOM,'SPLITX',SPLIT_LEN(1),ITYLCM)
        ALLOCATE(SPLITX(SPLIT_LEN(1)))
        CALL LCMGET(IPGEOM,'SPLITX',SPLITX)
      ENDIF
      IF(ABS(MONOP).NE.2.AND.NDIM.GE.2) THEN
        CALL LCMLEN(IPGEOM,'SPLITY',SPLIT_LEN(2),ITYLCM)
        ALLOCATE(SPLITY(SPLIT_LEN(2)))
        CALL LCMGET(IPGEOM,'SPLITY',SPLITY)
      ENDIF
      IF(ABS(MONOP).NE.3.AND.NDIM.GE.3) THEN
        CALL LCMLEN(IPGEOM,'SPLITZ',SPLIT_LEN(3),ITYLCM)
        ALLOCATE(SPLITZ(SPLIT_LEN(3)))
        CALL LCMGET(IPGEOM,'SPLITZ',SPLITZ)
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
      ENDIF
      IF(ABS(MONOP).NE.2.AND.NDIM.GE.2) THEN
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
      ENDIF
      IF(ABS(MONOP).NE.3.AND.NDIM.GE.3) THEN
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

      IF(MONOP.EQ.-1) THEN
        MONOP2=1
      ELSE IF(MONOP.EQ.1) THEN
        MONOP2=2
      ELSE IF(MONOP.EQ.-2) THEN
        MONOP2=3
      ELSE IF(MONOP.EQ.2) THEN
        MONOP2=4
      ELSE IF(MONOP.EQ.-3) THEN
        MONOP2=5
      ELSE IF(MONOP.EQ.3) THEN
        MONOP2=6
      ELSE
        MONOP2=0
      ENDIF

*----
*  1D CARTESIAN CASE
*----

      IF(NDIM.EQ.1) THEN

      ! Define the corresponding discrete ordinate  

      CALL LCMLEN(IPTRK,'U',NLF,ITYLCM)
      CALL LCMGPD(IPTRK,'U',U_PTR)
      CALL LCMGPD(IPTRK,'W',W_PTR)
      CALL C_F_POINTER(U_PTR,U,(/ NLF /))
      CALL C_F_POINTER(W_PTR,W,(/ NLF /))
      
      IPN=1
      IPVAL=(U(1)-DIR(1))**2
      DO IP=2,NLF
        IF((U(IP)-DIR(1))**2.LT.IPVAL) THEN
        IPN=IP
        IPVAL=(U(IP)-DIR(1))**2
        ENDIF
      ENDDO

      IF(U(IPN).LT.0.0.AND.MONOP.EQ.-1) THEN
        CALL XABORT('PSOMON: X- AND BACKWARD ORIENTED SOURCE.')
      ELSEIF(U(IPN).GT.0.0.AND.MONOP.EQ.1) THEN
        CALL XABORT('PSOMON: X+ AND FOWARD ORIENTED SOURCE.')
      ENDIF

      ! Normalization
      NORM=SUM(ISOUR)
 
      ! Save source information
      DO IG=1,NG
        IND=1
        IF(ISISOUR(IG).EQ.1) THEN
        BSINFO(1,IG,IND)=MONOP2
        BSINFO(2,IG,IND)=IPN
        BS(1,IG,IND)=ISOUR(IG)/W(IPN)
        IND=IND+1
        ENDIF
        NBS(IG)=IND-1
      ENDDO 
      
*----
*  2D CARTESIAN CASE
*----

      ELSE IF(NDIM.EQ.2) THEN

      ! Define the corresponding discrete ordinate
      
      CALL LCMLEN(IPTRK,'DU',NPQ,ITYLCM)
      CALL LCMGPD(IPTRK,'DU',DU_PTR)
      CALL LCMGPD(IPTRK,'DE',DE_PTR)
      CALL LCMGPD(IPTRK,'W',W_PTR)
      CALL C_F_POINTER(DU_PTR,DU,(/ NPQ /))
      CALL C_F_POINTER(DE_PTR,DE,(/ NPQ /))
      CALL C_F_POINTER(W_PTR,W,(/ NPQ /))

      ALLOCATE(MN(2),P(2))
      INIT1=0
      VAL1=8.0
      VAL2=8.0
      DO 300 M=1,NPQ
        IF(W(M).EQ.0.0) GO TO 300
        IF(INIT1.EQ.0) THEN
          MN(1)=M
          VAL1=(DU(M)-DIR(1))**2+(DE(M)-DIR(2))**2
          INIT1=1
        ELSE
          VAL=(DU(M)-DIR(1))**2+(DE(M)-DIR(2))**2
          IF(VAL1.GT.VAL) THEN
            MN(2)=MN(1)
            MN(1)=M
            VAL2=VAL1
            VAL1=VAL
          ELSEIF(VAL2.GT.VAL) THEN
            MN(2)=M
            VAL2=VAL
          ENDIF
        ENDIF
  300 CONTINUE

      VALTOT=VAL1+VAL2
      P(1)=VAL1/VALTOT
      P(2)=VAL2/VALTOT

      IF(DU(MN(1)).LT.0.0.AND.DU(MN(2)).LT.0.0.AND.MONOP.EQ.-1) THEN
        CALL XABORT('PSOMON: X- AND BACKWARD ORIENTED SOURCE.')
      ELSEIF(DU(MN(1)).GT.0.0.AND.DU(MN(2)).GT.0.0.AND.MONOP.EQ.1) THEN
        CALL XABORT('PSOMON: X+ AND FOWARD ORIENTED SOURCE.')
      ELSEIF(DE(MN(1)).LT.0.0.AND.DE(MN(2)).LT.0.0.AND.MONOP.EQ.-2) THEN
        CALL XABORT('PSOMON: Y- AND BACKWARD ORIENTED SOURCE.')
      ELSEIF(DE(MN(1)).GT.0.0.AND.DE(MN(2)).GT.0.0.AND.MONOP.EQ.2) THEN
        CALL XABORT('PSOMON: Y+ AND FOWARD ORIENTED SOURCE.')
      ENDIF

      ! Normalization
      S=0.0
      IF(ABS(MONOP).EQ.1) THEN
        DO NS=1,NSOUR
          S=S+(YMAX(NS)-YMIN(NS))
        ENDDO
      ENDIF
      IF(ABS(MONOP).EQ.2) THEN
        DO NS=1,NSOUR
          S=S+(XMAX(NS)-XMIN(NS))
        ENDDO
      ENDIF
      NORM=SUM(ISOUR)*S

      ! Save source information
      DO IG=1,NG
      IND=1
      DO IDIR=1,2
        IF(ISISOUR(IG).EQ.1) THEN
        BSINFO(1,IG,IND)=MONOP2
        BSINFO(2,IG,IND)=MN(IDIR)
        IF(ABS(MONOP).EQ.1) THEN
          DO IY=1,LY
          DO NS=1,NSOUR
          IF(YMIN(NS).LE.Y(IY).AND.YMAX(NS).GE.Y(IY)) THEN
            BS(IY,IG,IND)=ISOUR(IG)/W(MN(IDIR))*P(IDIR)
          ELSE
            BS(IY,IG,IND)=0.0
          ENDIF
          ENDDO
          ENDDO
        ELSEIF(ABS(MONOP).EQ.2) THEN
          DO IX=1,LX
          DO NS=1,NSOUR
          IF(XMIN(NS).LE.X(IX).AND.XMAX(NS).GE.X(IX)) THEN
            BS(IX,IG,IND)=ISOUR(IG)/W(MN(IDIR))*P(IDIR)
          ELSE
            BS(IX,IG,IND)=0.0
          ENDIF
          ENDDO
          ENDDO
        ENDIF
        IND=IND+1
        ENDIF
      ENDDO 
      NBS(IG)=IND-1
      ENDDO

      DEALLOCATE(MN,P)
      IF(ALLOCATED(SPLITX)) DEALLOCATE(SPLITX)
      IF(ALLOCATED(SPLITY)) DEALLOCATE(SPLITY)

*----
*  3D CARTESIAN CASE
*----

      ELSE IF(NDIM.EQ.3) THEN

      ! Define the corresponding discrete ordinate
     
      CALL LCMLEN(IPTRK,'DU',NPQ,ITYLCM)
      CALL LCMGPD(IPTRK,'DU',DU_PTR)
      CALL LCMGPD(IPTRK,'DE',DE_PTR)
      CALL LCMGPD(IPTRK,'DZ',DZ_PTR)
      CALL LCMGPD(IPTRK,'W',W_PTR)
      CALL C_F_POINTER(DU_PTR,DU,(/ NPQ /))
      CALL C_F_POINTER(DE_PTR,DE,(/ NPQ /))
      CALL C_F_POINTER(DZ_PTR,DZ,(/ NPQ /))
      CALL C_F_POINTER(W_PTR,W,(/ NPQ /))

      ALLOCATE(MN(4),P(4))
      INIT1=0
      INIT2=0
      INIT3=0
      VAL1=12.0
      VAL2=12.0
      VAL3=12.0
      VAL4=12.0
      DO 400 M=1,NPQ
        IF(W(M).EQ.0.0) GO TO 400
        IF(INIT1.EQ.0) THEN
          MN(1)=M
          VAL1=(DU(M)-DIR(1))**2+(DE(M)-DIR(2))**2+(DZ(M)-DIR(3))**2
          INIT1=1
        ELSEIF(INIT2.EQ.0) THEN
          VAL=(DU(M)-DIR(1))**2+(DE(M)-DIR(2))**2+(DZ(M)-DIR(3))**2
          IF(VAL1.GT.VAL) THEN
            MN(2)=MN(1)
            MN(1)=M
            VAL2=VAL1
            VAL1=VAL
          ELSE
            MN(2)=M
            VAL2=VAL
          ENDIF
          INIT2=1
        ELSEIF(INIT3.EQ.0) THEN
          VAL=(DU(M)-DIR(1))**2+(DE(M)-DIR(2))**2+(DZ(M)-DIR(3))**2
          IF(VAL1.GT.VAL) THEN
            MN(3)=MN(2)
            MN(2)=MN(1)
            MN(1)=M
            VAL3=VAL2
            VAL2=VAL1
            VAL1=VAL
          ELSEIF(VAL2.GT.VAL) THEN
            MN(3)=MN(2)
            MN(2)=M
            VAL3=VAL2
            VAL2=VAL
          ELSE
            MN(3)=M
            VAL3=VAL
          ENDIF
          INIT3=1
        ELSE
          VAL=(DU(M)-DIR(1))**2+(DE(M)-DIR(2))**2+(DZ(M)-DIR(3))**2
          IF(VAL1.GT.VAL) THEN
            MN(4)=MN(3)
            MN(3)=MN(2)
            MN(2)=MN(1)
            MN(1)=M
            VAL4=VAL3
            VAL3=VAL2
            VAL2=VAL1
            VAL1=VAL
          ELSEIF(VAL2.GT.VAL) THEN
            MN(4)=MN(3)
            MN(3)=MN(2)
            MN(2)=M
            VAL4=VAL3
            VAL3=VAL2
            VAL2=VAL
          ELSEIF(VAL3.GT.VAL) THEN
            MN(4)=MN(3)
            MN(3)=M
            VAL4=VAL3
            VAL3=VAL
          ELSEIF(VAL4.GT.VAL) THEN
            MN(4)=M
            VAL4=VAL
          ENDIF
        ENDIF
  400 CONTINUE

      VALTOT=VAL1+VAL2+VAL3+VAL4
      P(1)=VAL1/VALTOT
      P(2)=VAL2/VALTOT
      P(3)=VAL3/VALTOT
      P(4)=VAL4/VALTOT

      IF(DU(MN(1)).LT.0.0.AND.DU(MN(2)).LT.0.0.AND.
     1   DU(MN(3)).LT.0.0.AND.DU(MN(4)).LT.0.0.AND.MONOP.EQ.-1) THEN
        CALL XABORT('PSOMON: X- AND BACKWARD X-ORIENTED SOURCE.')
      ELSEIF(DU(MN(1)).GT.0.0.AND.DU(MN(2)).GT.0.0.AND.
     1       DU(MN(3)).GT.0.0.AND.DU(MN(4)).GT.0.0.AND.MONOP.EQ.1) THEN
        CALL XABORT('PSOMON: X+ AND FOWARD X-ORIENTED SOURCE.')
      ELSEIF(DE(MN(1)).LT.0.0.AND.DE(MN(2)).LT.0.0.AND.
     1       DE(MN(3)).LT.0.0.AND.DE(MN(4)).LT.0.0.AND.MONOP.EQ.-2) THEN
        CALL XABORT('PSOMON: Y- AND BACKWARD Y-ORIENTED SOURCE.')
      ELSEIF(DE(MN(1)).GT.0.0.AND.DE(MN(2)).GT.0.0.AND.
     1       DE(MN(3)).GT.0.0.AND.DE(MN(4)).GT.0.0.AND.MONOP.EQ.2) THEN
        CALL XABORT('PSOMON: Y+ AND FOWARD Y-ORIENTED SOURCE.')
      ELSEIF(DZ(MN(1)).LT.0.0.AND.DZ(MN(2)).LT.0.0.AND.
     1       DZ(MN(3)).LT.0.0.AND.DZ(MN(4)).LT.0.0.AND.MONOP.EQ.-3) THEN
        CALL XABORT('PSOMON: Z- AND BACKWARD Z-ORIENTED SOURCE.')
      ELSEIF(DZ(MN(1)).GT.0.0.AND.DZ(MN(2)).GT.0.0.AND.
     1       DZ(MN(3)).GT.0.0.AND.DZ(MN(4)).GT.0.0.AND.MONOP.EQ.3) THEN
        CALL XABORT('PSOMON: Z+ AND FOWARD Z-ORIENTED SOURCE.')
      ENDIF

      ! Normalization
      S=0.0
      IF(ABS(MONOP).EQ.1) THEN
        DO NS=1,NSOUR
          S=S+(YMAX(NS)-YMIN(NS))*(ZMAX(NS)-ZMIN(NS))
        ENDDO
      ENDIF
      IF(ABS(MONOP).EQ.2) THEN
        DO NS=1,NSOUR
          S=S+(XMAX(NS)-XMIN(NS))*(ZMAX(NS)-ZMIN(NS))
        ENDDO
      ENDIF
      IF(ABS(MONOP).EQ.3) THEN
        DO NS=1,NSOUR
          S=S+(XMAX(NS)-XMIN(NS))*(YMAX(NS)-YMIN(NS))
        ENDDO
      ENDIF
      NORM=SUM(ISOUR)*S

      ! Save source information
      DO IG=1,NG
      IND=1
      DO IDIR=1,4
        IF(ISISOUR(IG).EQ.1) THEN
        BSINFO(1,IG,IND)=MONOP2
        BSINFO(2,IG,IND)=MN(IDIR)
        IF(ABS(MONOP).EQ.1) THEN
          DO IY=1,LY
          DO IZ=1,LZ
            IPOS=(IY-1)*LZ+IZ
            DO NS=1,NSOUR
            IF(YMIN(NS).LE.Y(IY).AND.YMAX(NS).GE.Y(IY).AND.
     1      ZMIN(NS).LE.Z(IZ).AND.ZMAX(NS).GE.Z(IZ)) THEN
              BS(IPOS,IG,IND)=ISOUR(IG)/W(MN(IDIR))*P(IDIR)
            ELSE
              BS(IPOS,IG,IND)=0.0
            ENDIF
            ENDDO
          ENDDO
          ENDDO
        ELSEIF(ABS(MONOP).EQ.2) THEN
          DO IX=1,LX
          DO IZ=1,LZ
            IPOS=(IX-1)*LZ+IZ
            DO NS=1,NSOUR
            IF(XMIN(NS).LE.X(IX).AND.XMAX(NS).GE.X(IX).AND.
     1      ZMIN(NS).LE.Z(IZ).AND.ZMAX(NS).GE.Z(IZ)) THEN
              BS(IPOS,IG,IND)=ISOUR(IG)/W(MN(IDIR))*P(IDIR)
            ELSE
              BS(IPOS,IG,IND)=0.0
            ENDIF
            ENDDO
          ENDDO
          ENDDO
        ELSEIF(ABS(MONOP).EQ.3) THEN
          DO IX=1,LX
          DO IY=1,LY
            IPOS=(IX-1)*LY+IY
            DO NS=1,NSOUR
            IF(XMIN(NS).LE.X(IX).AND.XMAX(NS).GE.X(IX).AND.
     1      YMIN(NS).LE.Y(IY).AND.YMAX(NS).GE.Y(IY)) THEN
              BS(IPOS,IG,IND)=ISOUR(IG)/W(MN(IDIR))*P(IDIR)
            ELSE
              BS(IPOS,IG,IND)=0.0
            ENDIF
            ENDDO
          ENDDO
          ENDDO
        ENDIF
        IND=IND+1
        ENDIF
      ENDDO 
      NBS(IG)=IND-1
      ENDDO  

      DEALLOCATE(MN,P)
      IF(ALLOCATED(SPLITX)) DEALLOCATE(SPLITX)
      IF(ALLOCATED(SPLITY)) DEALLOCATE(SPLITY)
      IF(ALLOCATED(SPLITZ)) DEALLOCATE(SPLITZ)

      ELSE
      CALL XABORT('SOUR: INVALID GEOMETRY, ONLY 1D, 2D AND 3D CARTESIAN'
     1 //' GEOMETRY ARE ACTUALLY IMPLEMENTED.')
      ENDIF
      RETURN
      END
