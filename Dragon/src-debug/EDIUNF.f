*DECK EDIUNF
      SUBROUTINE EDIUNF(IPGEO1,IPGEO2,HSYM)
*-----------------------------------------------------------------------
*
*Purpose:
* Unfold a geometry.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPGEO1  pointer to the original geometry (L_GEOM signature).
* IPGEO2  pointer to the unfolded geometry (L_GEOM signature).
* HSYM    type of symmetry: 
*         'DIAG' for diagonal symmetry; 
*         'SYMX' for symmetry relative to X axis;
*         'SYMY' for symmetry relative to Y axis.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)  IPGEO1,IPGEO2
      CHARACTER HSYM*4
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE),NCODE(6),ICODE(6),ITDI(8),ITSX(8),ITSY(8)
      REAL ZCODE(6)
      CHARACTER*12 HSIGN
      LOGICAL LDIAG,LSYMX,LSYMY,LMESHX,LMESHY,LHMIX
      SAVE ITDI,ITSX,ITSY
      DATA ITDI/6,5,8,7,2,1,4,3/
      DATA ITSX/7,6,5,8,3,2,1,4/
      DATA ITSY/5,8,7,6,1,4,3,2/
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX,MERGE,ITURN,IHMIX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MIX2,MERGE2,ITURN2,IHMIX2
      REAL, ALLOCATABLE, DIMENSION(:) :: MESHX,MESHY,MESHX2,MESHY2
*
      CALL LCMGTC(IPGEO1,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN .NE. 'L_GEOM') THEN
        CALL XABORT('EDIUNF: SIGNATURE OF INPUT LCM OBJECT IS '//HSIGN//
     >  '. L_GEOM EXPECTED.')
      ENDIF
      CALL LCMGET(IPGEO1,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.5) CALL XABORT('EDIUNF: CARTESIAN GEOMETRY EXPEC'
     > //'TED.')
      IF(ISTATE(8).NE.1) CALL XABORT('EDIUNF: EMBEDDED CELLS EXPECTED.')
      NX=ISTATE(3)
      NY=ISTATE(4)
      NK=ISTATE(6)
      CALL LCMGET(IPGEO1,'NCODE',NCODE)
      CALL LCMGET(IPGEO1,'ICODE',ICODE)
      CALL LCMGET(IPGEO1,'ZCODE',ZCODE)
      LDIAG=(NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)
      LSYMX=(NCODE(3).EQ.5).AND.(NCODE(1).NE.3)
      LSYMY=(NCODE(1).EQ.5).AND.(NCODE(4).NE.3)
*----
*  Recover original geometry
*----
      ALLOCATE(MIX(NK),MERGE(NK),ITURN(NK),IHMIX(NK))
      ALLOCATE(MESHX(NX+1),MESHY(NY+1))
      CALL LCMLEN(IPGEO1,'MESHX',ILONG,ITYLCM)
      LMESHX=ILONG.NE.0
      IF(LMESHX) CALL LCMGET(IPGEO1,'MESHX',MESHX)
      CALL LCMLEN(IPGEO1,'MESHY',ILONG,ITYLCM)
      LMESHY=ILONG.NE.0
      IF(LMESHY) CALL LCMGET(IPGEO1,'MESHY',MESHY)
      CALL LCMGET(IPGEO1,'MIX',MIX)
      CALL LCMLEN(IPGEO1,'MERGE',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
        CALL LCMGET(IPGEO1,'MERGE',MERGE)
      ELSE
        DO I=1,NK
          MERGE(I)=I
        ENDDO
      ENDIF
      CALL LCMLEN(IPGEO1,'TURN',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
        CALL LCMGET(IPGEO1,'TURN',ITURN)
      ELSE
        CALL XDISET(ITURN,NK,1)
      ENDIF
      CALL LCMLEN(IPGEO1,'HMIX',ILONG,ITYLCM)
      LHMIX=ILONG.NE.0
      IF(LHMIX) CALL LCMGET(IPGEO1,'HMIX',IHMIX)
*----
*  Build unfolded geometry geometry
*----
      NX2=NX
      NY2=NY
      IF(HSYM.EQ.'DIAG') THEN
        IF(.NOT.LDIAG) CALL XABORT('EDIUNF: DIAG UNFOLDING FAILURE.')
      ELSE IF(HSYM.EQ.'SYMX') THEN
        IF(.NOT.LSYMX) CALL XABORT('EDIUNF: SYMX UNFOLDING FAILURE.')
        NY2=2*NY-1
      ELSE IF(HSYM.EQ.'SYMY') THEN
        IF(.NOT.LSYMY) CALL XABORT('EDIUNF: SYMY UNFOLDING FAILURE.')
        NX2=2*NX-1
      ELSE
        CALL XABORT('EDIUNF: INVALID TYPE OF SYMMETRY AXIS.')
      ENDIF
      ALLOCATE(MIX2(NX2,NY2),MERGE2(NX2,NY2),ITURN2(NX2,NY2),
     > IHMIX2(NX2,NY2))
      ALLOCATE(MESHX2(NX2+1),MESHY2(NY2+1))
      IF(HSYM.EQ.'DIAG') THEN
        NCODE(1)=NCODE(3)
        NCODE(4)=NCODE(2)
        ICODE(1)=ICODE(3)
        ICODE(4)=ICODE(2)
        ZCODE(1)=ZCODE(3)
        ZCODE(4)=ZCODE(2)
        IOF=0
        DO IY=1,NY
          DO IX=IY,NX
            IOF=IOF+1
            IF(IOF.GT.NK) CALL XABORT('EDIUNF: NK OVERFLOW(1).')
            MIX2(IX,IY)=MIX(IOF)
            MIX2(IY,IX)=MIX(IOF)
            MERGE2(IX,IY)=MERGE(IOF)
            MERGE2(IY,IX)=MERGE(IOF)
            ITURN2(IX,IY)=ITURN(IOF)
            ITURN2(IY,IX)=ITDI(ITURN(IOF))
            IF(LHMIX) THEN
              IHMIX2(IX,IY)=IHMIX(IOF)
              IHMIX2(IY,IX)=IHMIX(IOF)
            ENDIF
          ENDDO
        ENDDO
        IF(LMESHX) THEN
          DO IX=1,NX+1
            MESHX2(IX)=MESHX(IX)
          ENDDO
        ENDIF
        IF(LMESHY) THEN
          DO IY=1,NY+1
            MESHY2(IY)=MESHY(IY)
          ENDDO
        ENDIF
      ELSE IF(HSYM.EQ.'SYMX') THEN
        NCODE(3)=NCODE(4)
        ICODE(3)=ICODE(4)
        ZCODE(3)=ZCODE(4)
        DO IY=1,NY
          DO IX=1,NX
            IOF=(IY-1)*NX+IX
            IF(IOF.GT.NK) CALL XABORT('EDIUNF: NK OVERFLOW(2).')
            MIX2(IX,NY+IY-1)=MIX(IOF)
            MIX2(IX,NY-IY+1)=MIX(IOF)
            MERGE2(IX,NY+IY-1)=MERGE(IOF)
            MERGE2(IX,NY-IY+1)=MERGE(IOF)
            ITURN2(IX,NY+IY-1)=ITURN(IOF)
            ITURN2(IX,NY-IY+1)=ITSX(ITURN(IOF))
            IF(LHMIX) THEN
              IHMIX2(IX,NY+IY-1)=IHMIX(IOF)
              IHMIX2(IX,NY-IY+1)=IHMIX(IOF)
            ENDIF
          ENDDO
        ENDDO
        IF(LMESHX) THEN
          DO IX=1,NX+1
            MESHX2(IX)=MESHX(IX)
          ENDDO
        ENDIF
        IF(LMESHY) THEN
          DO IY=1,NY+1
            MESHY2(NY+IY-1)=MESHY(IY)
          ENDDO
          DO IY=3,NY+1
            MESHY2(NY-IY+2)=MESHY2(NY-IY+3)-(MESHY(IY)-MESHY(IY-1))
          ENDDO
        ENDIF
      ELSE IF(HSYM.EQ.'SYMY') THEN
        NCODE(1)=NCODE(2)
        ICODE(1)=ICODE(2)
        ZCODE(1)=ZCODE(2)
        DO IY=1,NY
          DO IX=1,NX
            IOF=(IY-1)*NX+IX
            IF(IOF.GT.NK) CALL XABORT('EDIUNF: NK OVERFLOW(3).')
            MIX2(NX+IX-1,IY)=MIX(IOF)
            MIX2(NX-IX+1,IY)=MIX(IOF)
            MERGE2(NX+IX-1,IY)=MERGE(IOF)
            MERGE2(NX-IX+1,IY)=MERGE(IOF)
            ITURN2(NX+IX-1,IY)=ITURN(IOF)
            ITURN2(NX-IX+1,IY)=ITSX(ITURN(IOF))
            IF(LHMIX) THEN
              IHMIX2(NX+IX-1,IY)=IHMIX(IOF)
              IHMIX2(NX-IX+1,IY)=IHMIX(IOF)
            ENDIF
          ENDDO
        ENDDO
        IF(LMESHX) THEN
          DO IX=1,NX+1
            MESHX2(NX+IX-1)=MESHX(IX)
          ENDDO
          DO IX=3,NX+1
            MESHX2(NX-IX+2)=MESHX2(NX-IX+3)-(MESHX(IX)-MESHX(IX-1))
          ENDDO
        ENDIF
        IF(LMESHY) THEN
          DO IY=1,NY+1
            MESHY2(IY)=MESHY(IY)
          ENDDO
        ENDIF
      ENDIF
      CALL LCMEQU(IPGEO1,IPGEO2)
      CALL LCMPUT(IPGEO2,'NCODE',6,1,NCODE)
      CALL LCMPUT(IPGEO2,'ICODE',6,1,ICODE)
      CALL LCMPUT(IPGEO2,'ZCODE',6,2,ZCODE)
      IF(LMESHX) CALL LCMPUT(IPGEO2,'MESHX',NX2+1,2,MESHX)
      IF(LMESHY) CALL LCMPUT(IPGEO2,'MESHY',NY2+1,2,MESHY)
      CALL LCMPUT(IPGEO2,'MIX',NX2*NY2,1,MIX2)
      CALL LCMPUT(IPGEO2,'MERGE',NX2*NY2,1,MERGE2)
      CALL LCMPUT(IPGEO2,'TURN',NX2*NY2,1,ITURN2)
      IF(LHMIX) CALL LCMPUT(IPGEO2,'HMIX',NX2*NY2,1,IHMIX2)
      ISTATE(3)=NX2
      ISTATE(4)=NY2
      ISTATE(6)=NX2*NY2
      CALL LCMPUT(IPGEO2,'STATE-VECTOR',NSTATE,1,ISTATE)
*----
*  Release memory
*----
      DEALLOCATE(MESHY2,MESHX2,MESHY,MESHX)
      DEALLOCATE(IHMIX2,ITURN2,MERGE2,MIX2,IHMIX,ITURN,MERGE,MIX)
      RETURN
      END
