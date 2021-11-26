!
!-----------------------------------------------------------------------
!
!Purpose:
! Perform an homogenization based on a surfacic file.
!
!Copyright:
! Copyright (C) 2017 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s): A. Hebert
!
!-----------------------------------------------------------------------
!
MODULE EDIG2S_MOD
  
  USE PRECISION_AND_KINDS, ONLY : PDB

CONTAINS
  !
  FUNCTION EDIBAR(NODE,CX,CY) RESULT(TLAMB)
    !----
    !  Compute the barycentric coordinates of point (CX,CY) in a triangle
    !----
    REAL(PDB) :: NODE(6),CX,CY,TLAMB(3)
    !
    TLAMB(1) = ((NODE(4) - NODE(6))*(CX - NODE(5)) + (NODE(5) - NODE(3))*(CY - NODE(6))) / &
            ((NODE(4) - NODE(6))*(NODE(1) - NODE(5)) + (NODE(5) - NODE(3))*(NODE(2) - NODE(6)))
    TLAMB(2) = ((NODE(6) - NODE(2))*(CX - NODE(5)) + (NODE(1) - NODE(5))*(CY - NODE(6))) / &
           ((NODE(4) - NODE(6))*(NODE(1) - NODE(5)) + (NODE(5) - NODE(3))*(NODE(2) - NODE(6)))
    TLAMB(3) = 1.0D0 - TLAMB(1) - TLAMB(2)
  END FUNCTION EDIBAR
  !
  SUBROUTINE EDIG2S(IPRINT,IFGEO,NREG,NMERGE,IMERGE)
    !----
    !  Process RECT and TRIA data options
    !
    !Parameters: input
    ! IPRINT  print flag.
    ! IFGEO   unit file number of the surfacic file.
    ! NREG    number of regions.
    !
    !Parameters: input
    ! NMERGE  number of merged indices in array IMERGE.
    ! IMERGE  merged regions position.
    !
    !----
    USE SALGET_FUNS_MOD
    !----
    !  Subroutine arguments
    !----
    INTEGER IPRINT,IFGEO,NREG,NMERGE,IMERGE(NREG)
    !----
    !  Local variables
    !----
    INTEGER PREC,DATAIN(25),IPAR(5)
    REAL DATARE(25)
    REAL(PDB) CX,CY,DX,DY,SAA,SAB,ANGL,RPAR(5),TLAMB1(3),TLAMB2(3)
    REAL(PDB) NODX1,NODX2,NODY1,NODY2
    REAL(PDB), PARAMETER :: CONV=3.141592654_PDB/180.0_PDB
    PARAMETER(IFOUT0=0)
    CHARACTER NAME_GEOM*12,CARLIR*8,HSMG*131
    DOUBLE PRECISION DBLLIR
    !----
    !  Allocatable arrays
    !----
    INTEGER, DIMENSION(:), ALLOCATABLE :: NUM_MERGE,IFLUX,ITNODE
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICOUNT
    REAL(PDB), DIMENSION(:,:), ALLOCATABLE :: NODE
    !----
    !  Read homogeneous node definitions
    !----
    CALL REDGET(ITYPLU,NMERGE,REALIR,CARLIR,DBLLIR)
    IF(ITYPLU.NE.1) CALL XABORT('EDIG2S: INTEGER VARIABLE EXPECTED.')
    IF(NMERGE.LE.0) CALL XABORT('EDIG2S: INVALID VALUE OF NMERGE.')
    ALLOCATE(NODE(6,NMERGE),ITNODE(NMERGE))
    DO IM=1,NMERGE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU.NE.3) CALL XABORT('EDIG2S: CHARACTER VARIABLE EXPECTED.')
      IF(CARLIR.EQ.'RECT') THEN
        ITNODE(IM)=1
        DO I=1,4
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.2) CALL XABORT('EDIG2S: REAL VARIABLE EXPECTED(1).')
          NODE(I,IM)=REALIR
        ENDDO
        NODE(5:6,IM)=0.0D0
      ELSE IF(CARLIR.EQ.'TRIA') THEN
        ITNODE(IM)=2
        DO I=1,6
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.2) CALL XABORT('EDIG2S: REAL VARIABLE EXPECTED(2).')
          NODE(I,IM)=REALIR
        ENDDO
      ELSE
        CALL XABORT('EDIG2S: *RECT* OR *TRIA* KEYWORD EXPECTED.')
      ENDIF
    ENDDO
    !----
    !  Determine homogenization indices
    !----
    IF(IFGEO.EQ.0) CALL XABORT('EDIG2S: surfacic file not defined.')
    CALL SALGET(DATAIN,6,IFGEO,IFOUT0,'dimensions for geometry')
    NBNODE=DATAIN(3)
    NBELEM=DATAIN(4)
    NBFLUX=DATAIN(6)
    CALL SALGET(DATAIN,3,IFGEO,IFOUT0,'index kndex prec')
    INDEX=DATAIN(1)
    KNDEX=DATAIN(2)
    PREC=DATAIN(3)
    CALL SALGET(DATARE,1,IFGEO,IFOUT0,'eps')
    EPS=DATARE(1)
    ALLOCATE(NUM_MERGE(NBNODE))
    CALL SALGET(NUM_MERGE,NBNODE,IFGEO,IFOUT0,'FLUX INDEX PER NODE')
    IF(MAXVAL(NUM_MERGE).NE.NBFLUX) CALL XABORT('EDIG2S: inconsistent NBFLUX.')
    CALL SALGET(NAME_GEOM,IFGEO,IFOUT0,'NAMES OF MACROS')
    ALLOCATE(IFLUX(NBFLUX))
    CALL SALGET(IFLUX,NBFLUX,IFGEO,IFOUT0,'macro order number per flux region.')
    DEALLOCATE(IFLUX)
    ALLOCATE(ICOUNT(NBNODE,NMERGE))
    ICOUNT(:NBNODE,:NMERGE)=0
    DO IELEM=1,NBELEM
      IPAR(:)=0
      RPAR(:)=0.0
      CALL SALGET(IPAR,3,IFGEO,IFOUT0,'integer descriptors')
      ITYPE=IPAR(1)
      SELECT CASE (ITYPE)
        CASE (1)
        NBER=4
        CASE (2)
        NBER=3
        CASE (3)
        NBER=5
      CASE DEFAULT
        WRITE(6,'(1X,''==> SAL126: unknown type '',I3)') ITYPE
        CALL XABORT('EDIG2S: unknown element type.')
      END SELECT
      CALL SALGET(RPAR,NBER,IFGEO,IFOUT0,PREC,'real descriptors')
      IF(ITYPE.EQ.1) THEN
        CX=RPAR(1) ; CY=RPAR(2)
        DX=CX+RPAR(3) ; DY=CY+RPAR(4)
        DO IM=1,NMERGE
          IF(ITNODE(IM).EQ.1) THEN
            NODX1=NODE(1,IM) ; NODX2=NODE(2,IM)
            NODY1=NODE(3,IM) ; NODY2=NODE(4,IM)
            IF((CX.GE.NODX1-EPS).AND.(DX.LE.NODX2+EPS).AND. &
               (CY.GE.NODY1-EPS).AND.(DY.LE.NODY2+EPS)) THEN
              IF((ABS(CX-DX).LE.EPS).AND.(ABS(CX-NODX1).LE.EPS)) THEN ! left vertical side
                IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              ELSE IF((ABS(CX-DX).LE.EPS).AND.(ABS(CX-NODX2).LE.EPS)) THEN ! right vertical side
                IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
              ELSE IF((ABS(CY-DY).LE.EPS).AND.(ABS(CY-NODY1).LE.EPS)) THEN ! lower horizontal side
                IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
              ELSE IF((ABS(CY-DY).LE.EPS).AND.(ABS(CY-NODY2).LE.EPS)) THEN ! upper horizontal side
                IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              ELSE IF((ABS(CX-DX).LE.EPS).OR.(ABS(CY-DY).LE.EPS)) THEN
                IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
                IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
              ENDIF
            ENDIF
          ELSE IF(ITNODE(IM).EQ.2) THEN
            TLAMB1=EDIBAR(NODE(1,IM),CX,CY)
            TLAMB2=EDIBAR(NODE(1,IM),DX,DY)
            IF((TLAMB1(1).GE.-EPS).AND.(TLAMB1(2).GE.-EPS).AND.(TLAMB1(3).GE.-EPS).AND. &
               (TLAMB2(1).GE.-EPS).AND.(TLAMB2(2).GE.-EPS).AND.(TLAMB2(3).GE.-EPS)) THEN
              IF((ABS(TLAMB1(1)).LE.EPS).AND.(ABS(TLAMB2(1)).LE.EPS)) THEN
                IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
              ELSE IF((ABS(TLAMB1(2)).LE.EPS).AND.(ABS(TLAMB2(2)).LE.EPS)) THEN
                IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
              ELSE IF((ABS(TLAMB1(3)).LE.EPS).AND.(ABS(TLAMB2(3)).LE.EPS)) THEN
                IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
              ELSE
                IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
                IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ELSE IF(ITYPE.EQ.2) THEN
        CX=RPAR(1) ; CY=RPAR(2)
        DO IM=1,NMERGE
          IF(ITNODE(IM).EQ.1) THEN
            NODX1=NODE(1,IM) ; NODX2=NODE(2,IM)
            NODY1=NODE(3,IM) ; NODY2=NODE(4,IM)
            IF((CX.GE.NODX1-EPS).AND.(CX.LE.NODX2+EPS).AND. &
               (CY.GE.NODY1-EPS).AND.(CY.LE.NODY2+EPS)) THEN
              IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
            ENDIF
          ELSE IF(ITNODE(IM).EQ.2) THEN
            TLAMB1=EDIBAR(NODE(1,IM),CX,CY)
            IF((TLAMB1(1).GE.-EPS).AND.(TLAMB1(2).GE.-EPS).AND.(TLAMB1(3).GE.-EPS)) THEN
              IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
            ENDIF
          ENDIF
        ENDDO
      ELSE IF(ITYPE.EQ.3) THEN
        SAA=RPAR(4) ; SAB=SAA+RPAR(5)
        IF(SAB>SAA) THEN
          ANGL=(SAB+SAA)*0.5
        ELSE
          ANGL=(SAB+SAA)*0.5+180.0
        ENDIF
        CX=RPAR(1)+COS(ANGL*CONV)*RPAR(3) ; CY=RPAR(2)+SIN(ANGL*CONV)*RPAR(3)
        DO IM=1,NMERGE
          IF(ITNODE(IM).EQ.1) THEN
            NODX1=NODE(1,IM) ; NODX2=NODE(2,IM)
            NODY1=NODE(3,IM) ; NODY2=NODE(4,IM)
            IF((CX.GE.NODX1-EPS).AND.(CX.LE.NODX2+EPS).AND. &
               (CY.GE.NODY1-EPS).AND.(CY.LE.NODY2+EPS)) THEN
              IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
            ENDIF
          ELSE IF(ITNODE(IM).EQ.2) THEN
            TLAMB1=EDIBAR(NODE(1,IM),CX,CY)
            IF((TLAMB1(1).GE.-EPS).AND.(TLAMB1(2).GE.-EPS).AND.(TLAMB1(3).GE.-EPS)) THEN
              IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    IMERGE(:NREG)=0
    ITEST=0
    DO IM=1,NMERGE
      DO INODE=1,NBNODE
        IF(ICOUNT(INODE,IM).GT.0) THEN
          IF(IMERGE(NUM_MERGE(INODE)).NE.0) THEN
            WRITE(HSMG,'(46HEDIG2S: inconsistent homogenization in mixture,I8, &
            & 11h, g2s node=,I8,1h.)') IM,INODE
            CALL XABORT(HSMG)
          ENDIF
          IMERGE(NUM_MERGE(INODE))=IM
          ITEST=ITEST+1
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(NUM_MERGE,ICOUNT,ITNODE,NODE)
    IF(IPRINT.GT.0) THEN
      WRITE(6,'(53H EDIG2S: NUMBER OF NODES PROCESSED BY HOMOGENIZATION=,I8/ &
      & 9X,32HNUMBER OF NODES IN THE GEOMETRY=,12X,I8/ &
      & 9X,31HNUMBER OF HOMOGENEOUS MIXTURES=,13X,I8)') ITEST,NBNODE,NMERGE
    ENDIF
    RETURN
  END SUBROUTINE EDIG2S
END MODULE EDIG2S_MOD
