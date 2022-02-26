*DECK NXTLHA
      FUNCTION NXTLHA(IPRINT,ITST  ,NDIM  ,MXMESH,LINMAX,
     >                MESH  ,ORITRK,DIRTRK,DCMESH,
     >                NBCOR ,NBSINT,ISINT ,TRKLSI)
*
*----------
*
*Purpose:
* To track an hexagonal assembly in 2-D or 3-D geometry
* using the NXT tracking procedure.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* ITST    type of tracking, where:
*         =-1   only the exact geometry
*               is considered taking into account the
*               submesh in each direction;
*         = 0   only the global geometry
*               is considered without taking into account the
*               submesh in each direction;
*         = 1   both the global
*               geometry (as a first step) and the exact geometry
*               are considered taking into account the
*               submesh in each direction.
* NDIM    dimension of problem.
* MXMESH  maximum number of spatial subdivision in
*         $X$, $Y$ or $Z$.
* LINMAX  maximum number of segments in a track.
* MESH    effective number of spatial subdivision in
*         each direction ($X$, $Y$ and $Z$).
* ORITRK  a point on the track (origin).
* DIRTRK  the track direction (director cosines).
* DCMESH  spatial description of the assembly.
*
*Parameters: output
* NXTLHA  number of side intersections.
* NBCOR   number of corner found for each external faces.
* NBSINT  number of surface crossed by track.
* ISINT   direction of plane intersected and
*         the surfaces crossed by the track.
* TRKLSI  the surface intersection distance.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,ITST,NDIM,MXMESH,LINMAX
      INTEGER          MESH(NDIM)
      DOUBLE PRECISION ORITRK(NDIM),
     >                 DIRTRK(NDIM),
     >                 DCMESH(-1:MXMESH,5)
      INTEGER          NBCOR(2),NBSINT
      INTEGER          ISINT(0:5,LINMAX)
      DOUBLE PRECISION TRKLSI(LINMAX)
      INTEGER          NXTLHA
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTLHA')
      DOUBLE PRECISION DCUTOF,DZERO,DONE,DTWO
      PARAMETER       (DCUTOF=1.0D-9,DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          IC1,IH,IDIR,IZ,ITF,ILF,ICELL
      DOUBLE PRECISION SQ3,SO2,SQ3S,SQ3O2S,HO2,OT(3)
      DOUBLE PRECISION PUX,PUY,PVX,PVY,PWX,PWY,PZ,
     >                 SLPUX,SLPUY,SLPVX,SLPVY,SLPWX,SLPWY,SLPZ
      DOUBLE PRECISION TINT(2,4),XY(4,12)
      INTEGER          IS,NS,IKS(8),IBL,IEL,ITL,IIS
      DOUBLE PRECISION DK(8),DDK
      INTEGER          NSETL,IB,IC,NTTP,NMOVE
*----
*  Verify ITST option and reset to default value if invalid
*----
      IF(ITST .LT. -1 .OR. ITST .GT. 1) THEN
*----
*  Reset ITST=1 (complete analysis) if the value of ITST is invalid
*----
        ITST=1
      ENDIF
      NBCOR(1)=1
      NBCOR(2)=1
      NBSINT=0
*----
*  Initialise output vectors
*----
      CALL XDISET(ISINT ,6*LINMAX,0)
      CALL XDDSET(TRKLSI,LINMAX,DZERO)
*----
*  Print header if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6011) 'HexagonSIDE={         '
        WRITE(IOUT,6012) DCMESH(0,1)
        WRITE(IOUT,6013)
        WRITE(IOUT,6011) 'HexagonX={         '
        WRITE(IOUT,6012) (DCMESH(IC1,1),IC1=1,MESH(1))
        WRITE(IOUT,6013)
        WRITE(IOUT,6011) 'HexagonY={         '
        WRITE(IOUT,6012) (DCMESH(IC1,2),IC1=1,MESH(1))
        WRITE(IOUT,6013)
        IF(NDIM .EQ. 3) THEN
          WRITE(IOUT,6011) 'HexagonZ={          '
          WRITE(IOUT,6012) (DCMESH(IC1,3),IC1=0,MESH(3))
          WRITE(IOUT,6013)
        ENDIF
        WRITE(IOUT,6011) 'trackorigin={       '
        WRITE(IOUT,6012) (ORITRK(IC1),IC1=1,NDIM)
        WRITE(IOUT,6013)
        WRITE(IOUT,6011) 'trackdirection={    '
        WRITE(IOUT,6012) (DIRTRK(IC1),IC1=1,NDIM)
        WRITE(IOUT,6013)
      ENDIF
      SQ3=SQRT(3.0D0)
      SO2=DCMESH(0,1)/DTWO
      SQ3S=SQ3*DCMESH(0,1)
      SQ3O2S=SQ3*SO2
      IEL=0
*----
*  Scan over each hexagon in assembly
*----
      DO IH=1,MESH(1)
*----
*  Move origin to the (x,y) cell center
*----
        NS=0
        DO IDIR=1,2
          OT(IDIR)=ORITRK(IDIR)-DCMESH(IH,IDIR)
        ENDDO
*----
*  Find intersection distance of U faces
*----
        PUX=OT(1)
        PUY=OT(2)
        SLPUX=DIRTRK(1)
        SLPUY=DIRTRK(2)
        TINT(1,1)=-(SQ3O2S+PUX)/SLPUX
        TINT(2,1)= (SQ3O2S-PUX)/SLPUX
*----
*  Find intersection distance of V faces
*----
        PVX=(OT(1)-SQ3*OT(2))/DTWO
        PVY=(OT(2)+SQ3*OT(1))/DTWO
        SLPVX=(DIRTRK(1)-SQ3*DIRTRK(2))/DTWO
        SLPVY=(DIRTRK(2)+SQ3*DIRTRK(1))/DTWO
        TINT(1,2)=-(SQ3O2S+PVX)/SLPVX
        TINT(2,2)= (SQ3O2S-PVX)/SLPVX
*----
*  Find intersection distance of W faces
*----
        PWX=(OT(1)+SQ3*OT(2))/DTWO
        PWY=(OT(2)-SQ3*OT(1))/DTWO
        SLPWX=(DIRTRK(1)+SQ3*DIRTRK(2))/DTWO
        SLPWY=(DIRTRK(2)-SQ3*DIRTRK(1))/DTWO
        TINT(1,3)=-(SQ3O2S+PWX)/SLPWX
        TINT(2,3)= (SQ3O2S-PWX)/SLPWX
        IF(NDIM .EQ. 2) THEN
*----
*  Test for U faces
*----
          DO ITF=1,2
            XY(ITF,1)=PUY+SLPUY*TINT(ITF,1)
            IF(ABS(XY(ITF,1)) .LE. SO2) THEN
              NS=NS+1
              DK(NS)=TINT(ITF,1)
              IKS(NS)=3*(2-ITF)+1
            ENDIF
          ENDDO
*----
*  Test for V faces
*----
          DO ITF=1,2
            XY(ITF,2)=PVY+SLPVY*TINT(ITF,2)
            IF(ABS(XY(ITF,2)) .LE. SO2) THEN
              NS=NS+1
              DK(NS)=TINT(ITF,2)
              IKS(NS)=3+3*(ITF-1)
            ENDIF
          ENDDO
*----
*  Test for W faces
*----
          DO ITF=1,2
            XY(ITF,3)=PWY+SLPWY*TINT(ITF,3)
            IF(ABS(XY(ITF,3)) .LE. SO2) THEN
              NS=NS+1
              DK(NS)=TINT(ITF,3)
              IKS(NS)=3*(2-ITF)+2
            ENDIF
          ENDDO
          IF( NS .EQ. 2 ) THEN
*----
*  Save cell crossing info
*----
            IF(DK(1) .GT. DK(2)) THEN
              ITF=IKS(1)
              IKS(1)=IKS(2)
              IKS(2)=ITF
              DDK=DK(1)
              DK(1)=DK(2)
              DK(2)=DDK
            ENDIF
*----
*  Combine segments in ISINT and TRKLSI
*----
            IF(IEL .EQ. 0) THEN
*----
*  First segment
*----
              DO IS=1,NS
                IEL=IEL+1
*----
*  HEX-FACE
*----
                ISINT(0,IEL)=IH
                ISINT(1,IEL)=IH
                ISINT(2,IEL)=-IKS(IS)
                ISINT(3,IEL)=1
                ISINT(4,IEL)=0
                TRKLSI(IEL)=DK(IS)
              ENDDO
            ELSE
*----
*  Remaining segments
*----
              IBL=IEL
              DO ITL=1,IEL
                DDK=TRKLSI(ITL)-DK(2)
*                write(6,*) IEL,ITL,IBL,TRKLSI(ITL),ISINT(0,ITL),
*     >                     DK(2),IH,DDK
                IF(DDK .GT. -DCUTOF) THEN
                  DO IBL=IEL,ITL,-1
                    TRKLSI(IBL+2)=TRKLSI(IBL)
                    ISINT(0,IBL+2)=ISINT(0,IBL)
                    ISINT(1,IBL+2)=ISINT(1,IBL)
                    ISINT(2,IBL+2)=ISINT(2,IBL)
                    ISINT(3,IBL+2)=ISINT(3,IBL)
                    ISINT(4,IBL+2)=ISINT(4,IBL)
                  ENDDO
                  IBL=ITL-1
                  GO TO 100
                ENDIF
              ENDDO
 100          CONTINUE
              IEL=IEL+2
              DO IS=1,NS
                IBL=IBL+1
*----
*  HEX-FACE
*----
                ISINT(0,IBL)=IH
                ISINT(1,IBL)=IH
                ISINT(2,IBL)=-IKS(IS)
                ISINT(3,IBL)=1
                ISINT(4,IBL)=0
                TRKLSI(IBL)=DK(IS)
              ENDDO
            ENDIF
*            WRITE(IOUT,*) '2-D tracking for cell= ',IH
*            WRITE(IOUT,'(2I10,F20.15)')
*     >      (ITF,IKS(ITF),DK(ITF),ITF=1,NS)
          ELSE IF( NS .GE. 1) THEN
            WRITE(IOUT,9000) IH,NS
            WRITE(IOUT,9001) (ITF,IKS(ITF),IH,DK(ITF),ITF=1,NS)
            WRITE(IOUT,9002) (XY(ITF,1),XY(ITF,2),XY(ITF,3),ITF=1,2)
            WRITE(IOUT,9003) IH
            WRITE(IOUT,9004)
     >      'U  ',PUX,SLPUX,PUY,SLPUY,TINT(1,1),TINT(2,1)
            WRITE(IOUT,9004)
     >      'V  ',PVX,SLPVX,PVY,SLPVY,TINT(1,2),TINT(2,2)
            WRITE(IOUT,9004)
     >      'W  ',PWX,SLPWX,PWY,SLPWY,TINT(1,3),TINT(2,3)
            CALL XABORT(NAMSBR//
     >      ': Problem with 2-D tracking -> '//
     >      ' line can only cross 0 or 2 surfaces in a cell')
          ENDIF
        ELSE
*----
*  scan over Z planes
*----
          DO IZ=1,MESH(3)
            ICELL=IH+(IZ-1)*MESH(1)
*----
*  Move origin to center of the plane in Z
*----
            OT(3)=ORITRK(3)-(DCMESH(IZ,IDIR)+DCMESH(IZ-1,IDIR))/DTWO
            HO2=(DCMESH(IZ,IDIR)-DCMESH(IZ-1,IDIR))/DTWO
*----
*  Find intersection distance of Z faces
*----
            PZ=OT(3)
            SLPZ=DIRTRK(3)
            TINT(1,4)=-(HO2+PZ)/SLPZ
            TINT(2,4)=(HO2-PZ)/SLPZ
*----
*  Test for U faces
*----
            DO ITF=1,2
              XY(ITF,1)=PUY+SLPUY*TINT(ITF,1)
              XY(ITF,2)=PZ+SLPZ*TINT(ITF,1)
              IF( (ABS(XY(ITF,1)) .LE. SO2) .AND.
     >            (ABS(XY(ITF,2)) .LE. HO2   ) ) THEN
                NS=NS+1
                DK(NS)=TINT(ITF,1)
                IKS(NS)=3*(2-ITF)+1
              ENDIF
            ENDDO
*----
*  Test for V faces
*----
            DO ITF=1,2
              XY(ITF,3)=PVY+SLPVY*TINT(ITF,2)
              XY(ITF,4)=PZ+SLPZ*TINT(ITF,2)
              IF( (ABS(XY(ITF,3)) .LE. SO2) .AND.
     >            (ABS(XY(ITF,4)) .LE. HO2   ) ) THEN
                NS=NS+1
                DK(NS)=TINT(ITF,2)
                IKS(NS)=3+3*(ITF-1)
              ENDIF
            ENDDO
*----
*  Test for W faces
*----
            DO ITF=1,2
              XY(ITF,5)=PWY+SLPWY*TINT(ITF,3)
              XY(ITF,6)=PZ+SLPZ*TINT(ITF,3)
              IF( (ABS(XY(ITF,5)) .LE. SO2) .AND.
     >            (ABS(XY(ITF,6)) .LE. HO2   ) ) THEN
                NS=NS+1
                DK(NS)=TINT(ITF,3)
                IKS(NS)=3*(2-ITF)+2
              ENDIF
            ENDDO
*----
*  Test for Z faces
*----
            DO ITF=1,2
              XY(ITF,7)=PUX+SLPUX*TINT(ITF,4)
              XY(ITF,8)=PUY+SLPUY*TINT(ITF,4)
              XY(ITF,9)=PVX+SLPVX*TINT(ITF,4)
              XY(ITF,10)=PVY+SLPVY*TINT(ITF,4)
              XY(ITF,11)=PWX+SLPWX*TINT(ITF,4)
              XY(ITF,12)=PWY+SLPWY*TINT(ITF,4)
              IF( ((ABS(XY(ITF,7)) .LE. SQ3O2S) .AND.
     >             (ABS(XY(ITF,8)) .LE. SO2   ) ) .OR.
     >            ((ABS(XY(ITF,9)) .LE. SQ3O2S) .AND.
     >             (ABS(XY(ITF,10)) .LE. SO2   ) ) .OR.
     >            ((ABS(XY(ITF,11)) .LE. SQ3O2S) .AND.
     >             (ABS(XY(ITF,12)) .LE. SO2   ) ) ) THEN
                NS=NS+1
                DK(NS)=TINT(ITF,4)
                IKS(NS)=6+ITF
              ENDIF
            ENDDO
****
*  remove ENDDO and put at the end of loop.
****
*          ENDDO
          IF( NS .EQ. 2 ) THEN
*----
*  Save cell crossing info
*----
            IF(DK(1) .GT. DK(2)) THEN
              ITF=IKS(1)
              IKS(1)=IKS(2)
              IKS(2)=ITF
              DDK=DK(1)
              DK(1)=DK(2)
              DK(2)=DDK
            ENDIF
*----
*  Combine segments in ISINT and TRKLSI
*----
            IF(IEL .EQ. 0) THEN
*----
*  First segment
*----
              DO IS=1,NS
                IEL=IEL+1
                ISINT(0,IEL)=(IZ-1)*MESH(1)+IH
                IF(IKS(IS) .GT. 6) THEN
*----
*  Z-FACE
*----
                  ISINT(1,IEL)=IH
                  ISINT(2,IEL)=0
                  ISINT(3,IEL)=IZ
                  ISINT(4,IEL)=-IKS(IS)
                ELSE
*----
*  HEX-FACE
*----
                  ISINT(1,IEL)=IH
                  ISINT(2,IEL)=-IKS(IS)
                  ISINT(3,IEL)=IZ
                  ISINT(4,IEL)=0
                ENDIF
                TRKLSI(IEL)=DK(IS)
              ENDDO
            ELSE
*----
*  Remaining segments
*----
              IBL=IEL
              DO ITL=1,IEL
                IF(DK(2) .LE. TRKLSI(ITL) ) THEN
                  DO IBL=IEL,ITL,-1
                    TRKLSI(IBL+2)=TRKLSI(IBL)
                    ISINT(0,IBL+2)=ISINT(0,IBL)
                    ISINT(1,IBL+2)=ISINT(1,IBL)
                    ISINT(2,IBL+2)=ISINT(2,IBL)
                    ISINT(3,IBL+2)=ISINT(3,IBL)
                    ISINT(4,IBL+2)=ISINT(4,IBL)
                  ENDDO
                  IBL=ITL-1
                  GO TO 110
                ENDIF
              ENDDO
 110          CONTINUE
              IEL=IEL+2
              DO IS=1,NS
                IBL=IBL+1
                ISINT(0,IBL)=(IZ-1)*MESH(1)+IH
                IF(IKS(IS) .GT. 6) THEN
*----
*  Z-FACE
*----
                  ISINT(1,IBL)=IH
                  ISINT(2,IBL)=0
                  ISINT(3,IBL)=IZ
                  ISINT(4,IBL)=-IKS(IS)
                ELSE
*----
*  HEX-FACE
*----
                  ISINT(1,IBL)=IH
                  ISINT(2,IBL)=-IKS(IS)
                  ISINT(3,IBL)=IZ
                  ISINT(4,IBL)=0
                ENDIF
                TRKLSI(IBL)=DK(IS)
              ENDDO
            ENDIF
*            write(IOUT,*) '3-D tracking results'
*            write(IOUT,'(2I10/(2I10,F20.15))') IH,IZ,
*     >      (ITF,IKS(ITF),DK(ITF),ITF=1,2)
          ELSE IF( NS .GE. 1) THEN
            WRITE(IOUT,9010) IH,IZ,NS
            WRITE(IOUT,9011) (ITF,IKS(ITF),IH,IZ,DK(ITF),ITF=1,NS)
            WRITE(IOUT,9002) ((XY(ITF,ILF),ITF=1,2),ILF=1,12)
            WRITE(IOUT,9003) IH
            WRITE(IOUT,9004)
     >      'U  ',PUX,SLPUX,PUY,SLPUY,TINT(1,1),TINT(2,1)
            WRITE(IOUT,9004)
     >      'V  ',PVX,SLPVX,PVY,SLPVY,TINT(1,2),TINT(2,2)
            WRITE(IOUT,9004)
     >      'W  ',PWX,SLPWX,PWY,SLPWY,TINT(1,3),TINT(2,3)
            WRITE(IOUT,9004)
     >      'Z  ',PZ,SLPZ,TINT(1,4),TINT(2,4)
            CALL XABORT(NAMSBR//
     >      ': Problem with 3-D tracking -> '//
     >      ' line can only cross 0 or 2 surfaces in a cell')
          ENDIF
****
*  ENDDO inserted here.
****
          ENDDO
        ENDIF
      ENDDO
      NXTLHA=IEL
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6030)
        DO IS=1,NXTLHA
          IF(ISINT(0,IS) .EQ. 0) THEN
            WRITE(IOUT,6020) IS,
     >      (ISINT(IDIR,IS),IDIR=0,5),TRKLSI(IS)
          ELSE
            WRITE(IOUT,6022) IS,
     >      (ISINT(IDIR,IS),IDIR=0,5),TRKLSI(IS)
          ENDIF
        ENDDO
      ENDIF
*----
*  Test if all faces are consecutive otherwise
*  Add intermediate region containing a cell with
*  id=0
*----
      NSETL=0
      ITL=IEL-1
      DO IS=IEL/2,2,-1
        DDK=TRKLSI(ITL)-TRKLSI(ITL-1)
*        DCUT=DCUTOF*ABS(TRKLSI(ITL)+TRKLSI(ITL-1))/DTWO
        IF(ABS(DDK) .GT. DCUTOF) THEN
          NMOVE=NXTLHA+NSETL
          NSETL=NSETL+1
          DO IBL=NMOVE,ITL,-1
            TRKLSI(IBL+1)=TRKLSI(IBL)
            DO IIS=0,4
              ISINT(IIS,IBL+1)=ISINT(IIS,IBL)
            ENDDO
          ENDDO
          DO IIS=0,4
            ISINT(IIS,ITL)=0
          ENDDO
          TRKLSI(ITL)=DDK
        ENDIF
        ITL=ITL-2
      ENDDO
      NXTLHA=NXTLHA+NSETL
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6032)
        DO IS=1,NXTLHA
          IF(ISINT(0,IS) .EQ. 0) THEN
            WRITE(IOUT,6020) IS,
     >      (ISINT(IDIR,IS),IDIR=0,5),TRKLSI(IS)
          ELSE
            WRITE(IOUT,6022) IS,
     >      (ISINT(IDIR,IS),IDIR=0,5),TRKLSI(IS)
          ENDIF
        ENDDO
      ENDIF
      NBSINT=NXTLHA+NSETL+1
      ITL=NBSINT
      IBL=0
      NTTP=0
      DO IS=NXTLHA,2,-1
        IC=ISINT(0,IS)
        IF(IC .EQ. 0) THEN
*----
*  Block outside cell
*----
*          write(6,*) IS,ITL,'Block outside cell'
          TRKLSI(ITL)=TRKLSI(IS)
          DO IIS=0,4
            ISINT(IIS,ITL)=ISINT(IIS,IS)
          ENDDO
          ITL=ITL-1
          IBL=0
          NTTP=NTTP+1
        ELSE
          IF(IBL .EQ. 0) THEN
*----
* Block last face
*----
*          write(6,*) IS,ITL,'Block last face'
            TRKLSI(ITL)=TRKLSI(IS)
            ISINT(0,ITL)=-ISINT(0,IS)
            DO IIS=1,4
              ISINT(IIS,ITL)=ISINT(IIS,IS)
            ENDDO
            ITL=ITL-1
*----
* Block last region
*----
*          write(6,*) IS,ITL,'Block last region'
            TRKLSI(ITL)=TRKLSI(IS)-TRKLSI(IS-1)
            DO IIS=0,4
              ISINT(IIS,ITL)=ISINT(IIS,IS-1)
            ENDDO          
            ISINT(2,ITL)=ISINT(1,IS-1)
            ITL=ITL-1
            IBL=1
            NTTP=NTTP+2
          ELSE
            IB=ISINT(0,IS-1)
            IF(IB .EQ. 0) THEN
*----
*  Block initial face
*----
*          write(6,*) IS,ITL,'Block initial face'
              TRKLSI(ITL)=TRKLSI(IS)
              ISINT(0,ITL)=-ISINT(0,IS)
              DO IIS=1,4
                ISINT(IIS,ITL)=ISINT(IIS,IS)
              ENDDO
              ITL=ITL-1
              NTTP=NTTP+1
            ELSE
*----
* Block region
*----
*          write(6,*) IS,ITL,'Block region'
              TRKLSI(ITL)=TRKLSI(IS)-TRKLSI(IS-1)
              DO IIS=0,4
                ISINT(IIS,ITL)=ISINT(IIS,IS-1)
              ENDDO          
              ISINT(2,ITL)=ISINT(1,IS-1)
              ITL=ITL-1
              NTTP=NTTP+1
            ENDIF
          ENDIF
        ENDIF
      ENDDO
*----
*  Initial face
*----
      IS=1
*      write(6,*) IS,ITL,'Line initial face'
      TRKLSI(ITL)=TRKLSI(IS)
      ISINT(0,ITL)=-ISINT(0,IS)
      DO IIS=1,4
        ISINT(IIS,ITL)=ISINT(IIS,IS)
      ENDDO
      NBSINT=NTTP
*----
*  Compress file for successive regions
*----
      IIS=0
      DO IS=1,NBSINT+1
        IF(ISINT(0,IS) .GT. 0) THEN
          DDK=ABS(TRKLSI(IS))
          IF(DDK .GT. DCUTOF) THEN
            IIS=IIS+1
            TRKLSI(IIS)=TRKLSI(IS)
            DO IDIR=0,4
              ISINT(IDIR,IIS)=ISINT(IDIR,IS)
            ENDDO
          ENDIF
        ELSE
          IIS=IIS+1
          TRKLSI(IIS)=TRKLSI(IS)
          DO IDIR=0,4
            ISINT(IDIR,IIS)=ISINT(IDIR,IS)
          ENDDO
        ENDIF  
      ENDDO
      NBSINT=IIS-1
*----
*  Print final track information
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6031) NBSINT+1
        DO IS=1,NBSINT+1
          IF(ISINT(0,IS) .EQ. 0) THEN
            WRITE(IOUT,6020) IS,
     >      (ISINT(IDIR,IS),IDIR=0,5),TRKLSI(IS)
          ELSE IF(ISINT(0,IS) .GT. 0) THEN
            WRITE(IOUT,6021) IS,
     >      (ISINT(IDIR,IS),IDIR=0,5),TRKLSI(IS)
          ELSE
            WRITE(IOUT,6022) IS,
     >      (ISINT(IDIR,IS),IDIR=0,5),TRKLSI(IS)
          ENDIF
        ENDDO
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6011 FORMAT(A20)
 6012 FORMAT(6(1X,F25.16,:,','))
 6013 FORMAT('};')
 6020 FORMAT('Tracks point  =',I10,' is outside cell   :',6I10,F25.16)
 6021 FORMAT('Tracks segment=',I10,' is in region      :',6I10,F25.16)
 6022 FORMAT('Tracks point  =',I10,' is on cell surface:',6I10,F25.16)
 6030 FORMAT('Intersection point with hexagonal faces')
 6031 FORMAT('Final ',I10,' track segments ')
 6032 FORMAT('Intersection point plus outside cells')
 9000 FORMAT(' Problem in 2-D tracking for cell= ',2I10)
 9001 FORMAT(3I10,F25.16)
 9002 FORMAT(2F25.16)
 9003 FORMAT('Cell =  ',I10)
 9004 FORMAT(A3,1X,6F25.16)
 9010 FORMAT(' Problem in 3-D tracking for cell= ',3I10)
 9011 FORMAT(4I10,F25.16)
      END
