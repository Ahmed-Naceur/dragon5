*DECK RESIND
      SUBROUTINE RESIND(IPMAP,IPMTX,NX,NY,NZ,LX,LY,LZ,MIX,NFUEL,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Update material index, it will store the negative fuel mixtures.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin, D. Sekki
*
*Parameters: input/output
* IPMAP  pointer to fuel-map information.
* IPMTX  pointer to matex information.
* NX     number of elements along x-axis in fuel map.
* NY     number of elements along y-axis in fuel map.
* NZ     number of elements along z-axis in fuel map.
* LX     number of elements along x-axis in geometry.
* LY     number of elements along y-axis in geometry.
* LZ     number of elements along z-axis in geometry.
* MIX    renumbered index over the fuel-map geometry.
* NFUEL  number of fuel types.
* IMPX   printing index (=0 for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP,IPMTX
      INTEGER NX,NY,NZ,LX,LY,LZ,MIX(NX*NY*NZ),NFUEL,IMPX
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER INDX(LX*LY*LZ),FTOT(NFUEL)
      REAL    MAPXX(NX+1),MAPYY(NY+1),MAPZZ(NZ+1),
     1        GEOXX(LX+1),GEOYY(LY+1),GEOZZ(LZ+1)
      TYPE(C_PTR) JPMAP
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMAT
*----
*  UPDATE MATERIAL INDEX
*----
      CALL LCMGET(IPMTX,'MESHX',GEOXX)
      CALL LCMGET(IPMTX,'MESHY',GEOYY)
      CALL LCMGET(IPMTX,'MESHZ',GEOZZ)
      JPMAP=LCMGID(IPMAP,'GEOMAP')
      CALL LCMGET(JPMAP,'MESHX',MAPXX)
      CALL LCMGET(JPMAP,'MESHY',MAPYY)
      CALL LCMGET(JPMAP,'MESHZ',MAPZZ)
      CALL XDISET(INDX,LX*LY*LZ,0)
      IF(IMPX.GT.2)WRITE(IOUT,*)'UPDATING MATERIAL INDEX'
      I1=0
      I2=0
      J1=0
      J2=0
      K1=0
      K2=0
      DO 52 KM=1,NZ
      DO 51 JM=1,NY
      DO 50 IM=1,NX
      DO IG=1,LX
      IF(MAPXX(IM).EQ.GEOXX(IG)) I1=IG
      IF(MAPXX(IM+1).EQ.GEOXX(IG+1))THEN
        I2=IG
        GOTO 10
      ENDIF
      ENDDO
   10 DO JG=1,LY
      IF(MAPYY(JM).EQ.GEOYY(JG)) J1=JG
      IF(MAPYY(JM+1).EQ.GEOYY(JG+1))THEN
        J2=JG
        GOTO 20
      ENDIF
      ENDDO
   20 DO KG=1,LZ
      IF(MAPZZ(KM).EQ.GEOZZ(KG)) K1=KG
      IF(MAPZZ(KM+1).EQ.GEOZZ(KG+1))THEN
        K2=KG
        GOTO 30
      ENDIF
      ENDDO
   30 IELM=(KM-1)*NX*NY+(JM-1)*NX +IM
      DO 42 KG=K1,K2
      DO 41 JG=J1,J2
      DO 40 IG=I1,I2
      IELG=(KG-1)*LX*LY+(JG-1)*LX+IG
      INDX(IELG)=MIX(IELM)
   40 CONTINUE
   41 CONTINUE
   42 CONTINUE
   50 CONTINUE
   51 CONTINUE
   52 CONTINUE
*     CHECK TOTAL NUMBER
      ITOT=0
      DO 60 IEL=1,LX*LY*LZ
      IF(INDX(IEL).NE.0)ITOT=ITOT+1
   60 CONTINUE
      NTOT=0
      CALL LCMGET(IPMTX,'FTOT',FTOT)
      DO 70 IFUEL=1,NFUEL
      NTOT=NTOT+FTOT(IFUEL)
   70 CONTINUE
      IF(ITOT.NE.NTOT) THEN
         WRITE(IOUT,'(/15H @RESIND: ITOT=,I8,6H NTOT=,I8)') ITOT,NTOT
         CALL XABORT('@RESIND: FOUND DIFFERENT TOTAL NUMBER OF FUEL MI'
     1   //'XTURES IN FUEL-MAP AND MATEX.')
      ENDIF
*     STORE NEGATIVE FUEL MIXTURES
      CALL LCMLEN(IPMTX,'MAT',LENGT,ITYP)
      ALLOCATE(IMAT(LENGT))
      CALL XDISET(IMAT,LENGT,0)
      CALL LCMGET(IPMTX,'MAT',IMAT)
      DO 100 IEL=1,LX*LY*LZ
      IF(INDX(IEL).NE.0)IMAT(IEL)=-INDX(IEL)
  100 CONTINUE
      CALL LCMPUT(IPMTX,'MAT',LENGT,1,IMAT)
      DEALLOCATE(IMAT)
      RETURN
      END
