*DECK NXTRTS
      SUBROUTINE NXTRTS(IPRINT,ITYPG ,MAXMSH,NREG  ,ITRN  ,ITST  ,
     >                  ITSYM ,NM    ,MIX   ,ISPLT ,DAMESH,
     >                  NMS   ,MIXR  ,ISPLTR,DAMESR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Rotate heagon with triangles according to reference turn and test,
* if required, in such a way that it satisfies intrinsic symmetries.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPRINT  intermediate printing level for output.
* ITYPG   geometry type.
* MAXMSH  maximum number of elements in MESH array.
* NREG    number of elements in MIX array.
* ITRN    geometry original turn number.
* ITST    flag for testing symmetry.
* ITSYM   flag for symmetries to test.
*
*Parameters: input/output
* NM      mesh size in all directions ($X$, $Y$, $Z$ and $R$).
* MIX     final mixture description for geometry (including HMIX).
* ISPLT   final split desctiption for geometry.
* DAMESH  final mesh description for geometry.
* NMS     mesh size after splitting.
*
*Parameters: temporary storage
* MIXR    mixture description for rotated geometry (including HMIX).
* ISPLTR  split desctiption for rotated geometry.
* DAMESR  mesh description for rotated geometry.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,ITYPG,MAXMSH,NREG,ITRN,ITST,
     >                 NM(4),ITSYM(4)
      INTEGER          ISPLT(0:MAXMSH-1,4),MIX(0:NREG-1,2)
      DOUBLE PRECISION DAMESH(-1:MAXMSH,4)
      INTEGER          NMS(4),ISPLTR(0:MAXMSH-1,4,2),
     >                 MIXR(0:NREG-1,2,2)
      DOUBLE PRECISION DAMESR(-1:MAXMSH,4,2)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTRTS')
      DOUBLE PRECISION DCUTOF,DZERO,DONE,DTWO
      PARAMETER       (DCUTOF=1.0D-6,DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Functions
*----
      INTEGER          NXTHRS
*----
*  Local variables
*----
      INTEGER          NX,NZ,NRP,NPG,IPG,IG,ICT,ITG,IKT,IX,IZ,IREG,
     >                 ITMI,IDMI,NNZ,NMT,NMTS
      INTEGER          ITM(3,2)
      DOUBLE PRECISION DDD
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Turn reference geometry (IPG=1)
*  and symmetric geometries (IPG=2,3,4,5)
*----
      NX=NM(1)
      NZ=MAX(NM(3),1)
      NRP=6*NX*NX
      NPG=1
      ITM(1,1)=1
      ITM(1,2)=1
      ITM(2,1)=2
      ITM(2,2)=2
      ITM(3,1)=3
      ITM(3,2)=3
      ICT=0
      NMT=0
      IF(ITST .EQ. 1) NPG=4
      DO IPG=1,NPG
        IF(IPG .EQ. 1) THEN
          IG=1
          ICT=ITRN
          DO IX=0,2*NX
            DAMESR(IX,1,IG)=DAMESH(IX,1)
          ENDDO
          ISPLTR(1,1,IG)=ISPLT(1,1)
          DAMESR(-1,1,IG)=DAMESH(-1,1)
        ELSE
          IG=2
          ITG=IPG-1
          IF(ABS(ITSYM(ITG)) .GE. 1) THEN
*----
*  Symmetry is valid
*  Determine final turn after applying symmetry on
*  current turn
*----
            IF(ITG .EQ. 1) THEN
*----
*  Hexagonal symetry
*----
              ICT=NXTHRS(ITRN,1)
            ELSE IF(ITG .EQ. 3) THEN
*----
*  Symmetry in Z
*----
              ICT=NXTHRS(ITRN,-1)
            ENDIF
          ELSE
*----
*  No need to test the geometry for this
*  intrinsic symmetry.
*----
            GO TO 1005
          ENDIF
        ENDIF
        IF(ICT .GT. 12 ) THEN
          IKT=12-ICT
        ELSE
          IKT=ICT
        ENDIF
        IF(IKT .LT. 0) THEN
          DAMESR(-1,3,IG)=-DAMESH(-1,3)
        ELSE
          DAMESR(-1,3,IG)=DAMESH(-1,3)
        ENDIF
        DO IX=0,2*NX
          DAMESR(IX,1,IG)=DAMESH(IX,1)
        ENDDO
        ISPLTR(1,1,IG)=ISPLT(1,1)
        DAMESR(-1,1,IG)=DAMESH(-1,1)
        IF     (ABS(IKT) .EQ. 1) THEN
          IF(IKT .LT. 0) THEN
            DO IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(NZ-IZ,3)-DAMESH(NZ-IZ-1,3)
              ISPLTR(IZ,3,IG)=ISPLT(NZ-IZ-1,3)
              ITMI=IZ*NRP
              IDMI=(NZ-IZ-1)*NRP
              DO IX=0,NX-1
                MIXR(ITMI,IG,1)=MIX(IDMI,1)
                MIXR(ITMI,IG,2)=MIX(IDMI,2)
                ITMI=ITMI+1
                IDMI=IDMI+1
              ENDDO
            ENDDO
          ELSE
            DO IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(IZ+1,3)-DAMESH(IZ,3)
              ISPLTR(IZ,3,IG)=ISPLT(IZ,3)
              ITMI=IZ*NRP
              IDMI=ITMI
              DO IX=0,NRP
                MIXR(ITMI,IG,1)=MIX(IDMI,1)
                MIXR(ITMI,IG,2)=MIX(IDMI,2)
                ITMI=ITMI+1
                IDMI=IDMI+1
             ENDDO
            ENDDO
          ENDIF
        ELSE
          CALL XABORT(NAMSBR//
     >      ': Symmetry not yet programmed')
        ENDIF
        IF(IPRINT .GE. 100) THEN
*----
*  Print turned mesh if required
*----
          WRITE(IOUT,6010) NX,NZ,NREG
          WRITE(IOUT,6011) 'MESHH ='
          WRITE(IOUT,6012) (DAMESR(IX,1,IG),IX=-1,2*NX)
          WRITE(IOUT,6011) 'SPLTH ='
          WRITE(IOUT,6013) ISPLTR(1,1,IG)
          IF(ITYPG .EQ. 13) THEN
            WRITE(IOUT,6011) 'MESHZ ='
            WRITE(IOUT,6012) (DAMESR(IZ,3,IG),IZ=-1,NZ)
            WRITE(IOUT,6011) 'SPLTZ ='
            WRITE(IOUT,6013) (ISPLTR(IZ,3,IG),IZ=0,NZ)
          ENDIF
          WRITE(IOUT,6011) 'MIX   ='
          WRITE(IOUT,6013) (MIXR(IREG,IG,1),IREG=0,NREG-1)
          WRITE(IOUT,6011) 'HMIX  ='
          WRITE(IOUT,6013) (MIXR(IREG,IG,2),IREG=0,NREG-1)
        ENDIF
        IF(IPG .GT. 1) THEN
*----
*  COMPARE GEOMETRY
*  1- MESH AND SPLIT IN Z
*  2- MIXTURES
*  3- OFFCENTER
*----
          NNZ=NM(ITM(3,1))
          IF(NNZ .NE. NM(ITM(3,2))) CALL XABORT(NAMSBR//
     >      ': Symmetry invalid with this mesh')
          DO IZ=0,NNZ-1
            DDD=ABS(DAMESR(IZ,3,1)-DAMESR(IZ,3,2))
            IF(DDD .GT. DCUTOF) CALL XABORT(NAMSBR//
     >      ': Symmetry invalid with this mesh')
            IF(ISPLTR(IZ,3,1) .NE. ISPLTR(IZ,3,2) )
     >      CALL XABORT(NAMSBR//
     >      ': Symmetry invalid with this split')
          ENDDO
          DO IREG=0,NREG-1
            IF(MIXR(IREG,1,1) .NE. MIXR(IREG,2,1)) CALL XABORT(NAMSBR//
     >      ': Symmetry invalid with this mixtures')
            IF(MIXR(IREG,1,2) .NE. MIXR(IREG,2,2)) CALL XABORT(NAMSBR//
     >      ': Symmetry invalid with this merging mixtures')
          ENDDO
          IF(DAMESR(-1,3,1) .NE. DAMESR(-1,3,2)) CALL XABORT(NAMSBR//
     >    ': Symmetry invalid with this off center')
        ELSE
*----
*  Reset reference geometry for turn
*----
          DO IX=0,2*NX
            DAMESH(IX,1)=DAMESR(IX,1,IG)
          ENDDO
          ISPLT(1,1)=ISPLTR(1,1,IG)
          DAMESH(-1,1)=DAMESR(-1,1,IG)
          NNZ=NM(ITM(3,1))
          NMT=NNZ
          NMTS=0
          DO IZ=0,NNZ-1
            NMTS=NMTS+ABS(ISPLTR(IZ,3,1))
          ENDDO
          IF(NMTS .NE. NMS(ITM(3,1))) CALL XABORT(NAMSBR//
     >       ': Global symmetry invalid with this split')
        ENDIF
 1005   CONTINUE
      ENDDO
*----
*  Reset final mesh (center+original turn)
*----
      NZ=NMT
      CALL XDDSET(DAMESH(-1,3),NM(3)+2,DZERO)
      NM(3)=NZ
      DDD=DZERO
      DO IZ=0,NZ-1
        DDD=DDD+DAMESR(IZ,3,1)
      ENDDO
      DDD=DDD/DTWO
      DAMESH(-1,3)=DAMESR(-1,3,1)
      DAMESH(0,3)=-DDD
      DO IZ=1,NZ
        DAMESH(IZ,3)=DAMESH(IZ-1,3)+DAMESR(IZ-1,3,1)
      ENDDO
      DO IZ=0,NZ
        ISPLT(IZ,3)=ISPLTR(IZ,3,1)
      ENDDO
      NMS(3)=0
      DO IZ=0,NZ-1
        NMS(3)=NMS(3)+ABS(ISPLT(IZ,3))
      ENDDO
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6010) NX,NZ,NREG
        WRITE(IOUT,6011) 'MESHH ='
        WRITE(IOUT,6012) (DAMESH(IX,1),IX=-1,2*NX)
        WRITE(IOUT,6011) 'SPLTH ='
        WRITE(IOUT,6013) ISPLT(1,1)
        IF(ITYPG .EQ. 13) THEN
          WRITE(IOUT,6011) 'MESHZ ='
          WRITE(IOUT,6012) (DAMESH(IZ,3),IZ=-1,NZ)
          WRITE(IOUT,6011) 'SPLTZ ='
          WRITE(IOUT,6013) (ISPLT(IZ,3),IZ=0,NZ)
        ENDIF
        WRITE(IOUT,6011) 'MIX   ='
        WRITE(IOUT,6013) (MIX(IREG,1),IREG=0,NREG-1)
        WRITE(IOUT,6011) 'HMIX  ='
        WRITE(IOUT,6013) (MIX(IREG,2),IREG=0,NREG-1)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  FORMATS
*----
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(1X,' DIMENSIONS =',5I10/1X,' ORIGINAL MESH ')
 6011 FORMAT(1X,A7)
 6012 FORMAT(5F15.9)
 6013 FORMAT(5I15)
      END
