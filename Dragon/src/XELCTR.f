*DECK XELCTR
      SUBROUTINE XELCTR(IFOLD,IFTRK,MXSUBO,MXSEGO,CUTOFX,ALBEDO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* EXCELL prismatic tracking.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IFOLD   unnormalized tracking file number (at input).     
* IFTRK   normalized tracking file number (at output).    
* MXSUBO  undefined.                                   
* MXSEGO  undefined.                                   
* CUTOFX  cutoff factor.                               
* ALBEDO  geometric albedos on external faces.         
*
*-----------------------------------------------------------------------
*


      IMPLICIT NONE

      INTEGER IFOLD,IFTRK,MXSUBO,MXSEGO
      REAL CUTOFX,ALBEDO(6)

      INTEGER NCOMNT,NSCRP,NDIM,ISPEC,NREG,NSOUT,NALBG,NCOR,NANGL,NRS,
     1 ICODE(6),II,JJ,NBTRK,MXSUB,MXSEG,NSUB,LINE,ITRAK,NOLDS,NNEWS,
     2 NCSEG
      REAL VOLMIN,ASCRP
      DOUBLE PRECISION WEIGHT,RCUT,DASCRP
      CHARACTER CTRK*4,COMENT*80
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATALB,NRSEG,KANGL
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLSUR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ANGLE,DENSTY
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SEGLEN
*---
*  Read Old Tracking File
*---
      READ (IFOLD) CTRK,NCOMNT,NSCRP,NSCRP
      DO II=1,NCOMNT
         READ(IFOLD) COMENT
      ENDDO
      READ (IFOLD) NDIM,ISPEC,NREG,NSOUT,NALBG,NCOR,NANGL,NSCRP,NSCRP
      IF(NALBG.LE.0.OR.NALBG.GT.6)THEN
         CALL XABORT('XELCTR: NALBG.GT.6.OR.NALBG.LE.0'//
     1               ' ON TRACKING FILE')
      ENDIF
      NRS=NREG+NSOUT+1
      ALLOCATE(MATALB(NRS),NRSEG(MXSEGO),KANGL(MXSUBO))
      ALLOCATE(VOLSUR(NRS),ANGLE(NDIM*NANGL),DENSTY(NANGL),
     1 SEGLEN(MXSEGO))
      READ (IFOLD) (VOLSUR(II),II=1,NRS)
      READ (IFOLD) (MATALB(II),II=1,NRS)
      READ (IFOLD) (ICODE(II),II=1,NALBG)
      READ (IFOLD) (ALBEDO(II),II=1,NALBG)
      READ (IFOLD) ((ANGLE((JJ-1)*NDIM+II),II=1,NDIM),JJ=1,NANGL)
      READ (IFOLD) (DENSTY(II),II=1,NANGL)
      VOLMIN=VOLSUR(NSOUT+2)
      DO II= NSOUT+2,NSOUT+NREG
        VOLMIN=MIN(VOLMIN,VOLSUR(II+1))
      ENDDO
      RCUT=VOLMIN*CUTOFX
      NBTRK= 0
      MXSUB= 0
      MXSEG= 0
 20   CONTINUE
         READ(IFOLD,END=40) NSUB,LINE,WEIGHT,(KANGL(II),II=1,NSUB),
     1                      (NRSEG(II),II=1,LINE),(SEGLEN(II),II=1,LINE)
         MXSUB=MAX(MXSUB,NSUB)
         MXSEG=MAX(MXSEG,LINE)
         NBTRK=NBTRK+1
      GOTO 20
   40 CONTINUE
*---
*  Construct New Tracking File
*---
      REWIND IFOLD
      READ (IFOLD) CTRK,NSCRP,NSCRP,NSCRP
      WRITE(IFTRK) CTRK,NCOMNT,NBTRK,0
      DO II=1,NCOMNT
         READ (IFOLD) COMENT
         WRITE(IFTRK) COMENT
      ENDDO
      READ (IFOLD) (NSCRP,II=1,8)
      WRITE(IFTRK) NDIM,ISPEC,NREG,NSOUT,NALBG,NCOR,NANGL,MXSUB,MXSEG
      READ (IFOLD) (ASCRP,II=-NSOUT,NREG)
      WRITE(IFTRK) (VOLSUR(II),II=1,NRS)
      READ (IFOLD) (NSCRP,II=-NSOUT,NREG)
      WRITE(IFTRK) (MATALB(II),II=1,NRS)
      READ (IFOLD) (NSCRP,II=1,NALBG)
      WRITE(IFTRK) (ICODE(II),II=1,NALBG)
      READ (IFOLD) (ASCRP,II=1,NALBG)
      WRITE(IFTRK) (ALBEDO(II),II=1,NALBG)
      READ (IFOLD) ((DASCRP,II=1,NDIM),JJ=1,NANGL)
      WRITE(IFTRK) ((ANGLE((JJ-1)*NDIM+II),II=1,NDIM),JJ=1,NANGL)
      READ (IFOLD) (DASCRP,II=1,NANGL)
      WRITE(IFTRK) (DENSTY(II),II=1,NANGL)
      DO ITRAK=1, NBTRK
         READ(IFOLD) NSUB,LINE,WEIGHT,(KANGL(II),II=1,NSUB),
     1                   (NRSEG(II),II=1,LINE),(SEGLEN(II),II=1,LINE)
         IF (RCUT.GT.0.0)THEN
            II=0
   23       CONTINUE
               IF (II.EQ.LINE) GO TO 25
               II=II+1
               IF (SEGLEN(II).LT.RCUT) THEN
                  IF (II.NE.LINE) THEN
                     DO JJ= II+1, LINE
                        NRSEG(JJ-1)=NRSEG(JJ)
                        SEGLEN(JJ-1)=SEGLEN(JJ)
                     ENDDO
                  ELSE
                     LINE=LINE-1
                     GOTO 25
                  ENDIF
                  LINE=LINE-1
                  II=II-1
               ENDIF
               GOTO 23
   25       CONTINUE
         ENDIF
         NOLDS=NRSEG(1)
         NCSEG=1
         DO II=2,LINE
            NNEWS=NRSEG(II)
            IF ((NNEWS.LT.0).OR.(NNEWS.NE.NOLDS)) THEN
               NOLDS=NNEWS
               NCSEG=NCSEG+1
               NRSEG(NCSEG)=NRSEG(II)
               SEGLEN(NCSEG)=SEGLEN(II)
            ELSEIF (NNEWS.EQ.NOLDS) THEN
               SEGLEN(NCSEG)=SEGLEN(NCSEG)+SEGLEN(II)
            ENDIF
         ENDDO
         WRITE(IFTRK) NSUB,NCSEG,WEIGHT,(KANGL(II),II=1,NSUB),
     1                (NRSEG(II),II=1,NCSEG),(SEGLEN(II),II=1,NCSEG)
      ENDDO
      DEALLOCATE(SEGLEN,DENSTY,ANGLE,VOLSUR)
      DEALLOCATE(KANGL,NRSEG,MATALB)
*
      RETURN
      END
