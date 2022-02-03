*DECK NXTMCB
      SUBROUTINE NXTMCB(IPTRK,NUCELL,MXGSUR,MXGREG,MAXPIN,IUNFLD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Add MC: specific geometry analysis info to NXTRecords.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): Romain Le Tellier
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* NUCELL  number of cell after unfolding in 
*         $X$, $Y$ and $Z$ directions.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* MAXPIN  maximum number of pins in a cell.
* IUNFLD  description of unfolded geometry.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER NUCELL(3),MXGSUR,MXGREG,MAXPIN,IUNFLD(2,NUCELL(1),
     > NUCELL(2),NUCELL(3))
*----
*  LOCAL VARIABLES
*----
      INTEGER NSTATE,IOUT
      PARAMETER(NSTATE=40,IOUT=6)
      INTEGER ESTATE(NSTATE)
      INTEGER K,J,I,ICEL,NREGC,NSURC,NTPIN,NREGF,NSURF,JJ,IPINO,IPIN
      CHARACTER NAMCEL*9,NAMREC*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IDSUR,IDREG,ITPIN
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: INDEX
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDEX(5,-MXGSUR:MXGREG,2),IDSUR(MXGSUR,2),
     1 IDREG(MXGREG,2),ITPIN(3,MAXPIN))
*----
*  CREATE COMPRESSED INDEX FOR ALL THE CELLS/PINS
*----
      DO 12 K=1,MAX(NUCELL(3),1)
      DO 11 J=1,NUCELL(2)
      DO 10 I=1,NUCELL(1)
      IF (IUNFLD(2,I,J,K).NE.1) GOTO 10
*     CELL LEVEL (1)
         ICEL=IUNFLD(1,I,J,K)
         WRITE(NAMCEL,'(A1,I8.8)') 'C',ICEL
         NAMREC=NAMCEL//'DIM'
         CALL XDISET(ESTATE,NSTATE,0)
         CALL LCMGET(IPTRK,NAMREC,ESTATE)
         NREGC=ESTATE(8)
         NSURC=ESTATE(9)
         NTPIN=ESTATE(18)
         NAMREC=NAMCEL//'VSI'
         CALL LCMGET(IPTRK,NAMREC,INDEX)
         NAMREC=NAMCEL//'SID'
         CALL LCMGET(IPTRK,NAMREC,IDSUR)
         NAMREC=NAMCEL//'RID'
         CALL LCMGET(IPTRK,NAMREC,IDREG)
         CALL NXTMCC(IPTRK,NAMCEL,NREGC,NSURC,NREGF,NSURF,INDEX,IDSUR,
     1        IDREG)
         IF(NTPIN.GT.0) THEN
*        PIN LEVEL (2)
            NAMREC=NAMCEL//'PNT'
            CALL LCMGET(IPTRK,NAMREC,ITPIN)
            IPINO=0
            DO 20 JJ=1,NTPIN
               IPIN=ITPIN(2,JJ)
               IF (IPIN.EQ.IPINO) GOTO 20
               WRITE(NAMCEL,'(A1,I8.8)') 'P',IPIN
               NAMREC=NAMCEL//'DIM'
               CALL XDISET(ESTATE,NSTATE,0)
               CALL LCMGET(IPTRK,NAMREC,ESTATE)
               NREGC=ESTATE(8)
               NSURC=ESTATE(9)
               NAMREC=NAMCEL//'VSI'
               CALL LCMGET(IPTRK,NAMREC,INDEX)
               NAMREC=NAMCEL//'SID'
               CALL LCMGET(IPTRK,NAMREC,IDSUR)
               NAMREC=NAMCEL//'RID'
               CALL LCMGET(IPTRK,NAMREC,IDREG)
               CALL NXTMCC(IPTRK,NAMCEL,NREGC,NSURC,NREGF,NSURF,INDEX,
     1              IDSUR,IDREG)
               IPINO=IPIN
 20         CONTINUE
         ENDIF
 10   CONTINUE
 11   CONTINUE
 12   CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ITPIN,IDREG,IDSUR,INDEX)
      RETURN
      END
