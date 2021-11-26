*DECK XELCOR
      SUBROUTINE XELCOR(IFILE1,IFILE2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Produce an equivalent tracking with NCOR=1.
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IFILE1  input tracking file.
* IFILE2  output tracking file.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER            IFILE1,IFILE2
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION   WEIGHT,WEIGHT2
      INTEGER            NCOMNT,NTRK,IFMT,IREC,IC,IR,NDIM,ISPEC,NV,NS,
     >                   NALBG,NCOR,NANGL,MXSUB,MXSEG,NSUB, LINE,NUNKNO
      CHARACTER          CTRK*4, COMENT*80
      INTEGER            IOUT
      PARAMETER        ( IOUT=6 )
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATALB,ICODE,NRSEG,KANGL
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLSUR,ALBEDO
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ANGLES,DENSTY
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SEGLEN
*----
*  SET NCOMNT, NUNKNO, NDIM, NALBG, NCOR ,NANGL AND MXSEG
*----
      READ (IFILE1,ERR=991) CTRK,NCOMNT,NTRK,IFMT
      DO 10 IC= 1, NCOMNT
         READ (IFILE1,ERR=991)
   10 CONTINUE
      READ (IFILE1,ERR=991) NDIM,ISPEC,NV,NS,NALBG,NCOR,NANGL,MXSUB,
     > MXSEG
      DO 20 IC= 1, 6
         READ (IFILE1,ERR=991)
   20 CONTINUE
*----
*  ALLOCATE SPACE TO COPY SUBSEQUENT RECORDS
*----
      NUNKNO= NV+NS+1
      ALLOCATE(MATALB(NUNKNO),ICODE(NALBG),NRSEG(MXSEG),KANGL(MXSUB))
      ALLOCATE(VOLSUR(NUNKNO),ALBEDO(NALBG),ANGLES(NDIM*NANGL),
     > DENSTY(NANGL),SEGLEN(MXSEG))
*----
*  COMPUTE THE NUMBER OF TRACKS
*----
      NTRK2=0
   30 CONTINUE
         READ (IFILE1,END=40,ERR=991) NSUB,LINE,WEIGHT,
     >                      (KANGL(IR),IR=1,NSUB),
     >                      (NRSEG(IR+1),IR=0,LINE-1),
     >                      (SEGLEN(IR+1),IR=0,LINE-1)
         IF(NSUB.GT.MXSUB) CALL XABORT('XELCOR: MXSUB OVERFLOW.')
         IF(NCOR.EQ.1) THEN
            NTRK2=NTRK2+1
         ELSE
           I1=1
           DO IR=2,NCOR
             IF(NRSEG(IR).NE.NRSEG(1)) I1=NCOR
           ENDDO
           I2=1
           DO IR=2,NCOR
             IF(NRSEG(LINE-NCOR+IR).NE.NRSEG(LINE-NCOR+1)) I2=NCOR
           ENDDO
           NTRK2=NTRK2+I1*I2
         ENDIF
      GO TO 30
   40 CONTINUE
*----
*  READ AND COPY FIRST RECORDS (HEADER, COMMENTS)
*----
      REWIND IFILE1
      IREC= 1
      READ (IFILE1,ERR=991) CTRK,NCOMNT,NTRK,IFMT
      WRITE(IFILE2,ERR=992) CTRK,NCOMNT,NTRK2,IFMT
      DO 50 IC= 1, NCOMNT
         IREC= IREC+1
         READ (IFILE1,ERR=991) COMENT
         WRITE(IFILE2,ERR=992) COMENT
   50 CONTINUE
*----
*  READ AND COPY MAIN RECORD AND GET USEFUL DIMENSIONS
*----
      IREC= IREC+1
      READ (IFILE1,ERR=991) NDIM,ISPEC,NV,NS,NALBG,NCOR,NANGL,MXSUB,
     > MXSEG
      WRITE(IFILE2,ERR=992) NDIM,ISPEC,NV,NS,NALBG,1,NANGL,MXSUB,MXSEG
      NUNKNO= NV+NS+1
*----
*  COPY ALL RECORDS BEFORE TRACKS
*----
      IREC= IREC+1
      READ (IFILE1,ERR=991) (VOLSUR(IR),IR=1,NUNKNO)
      WRITE(IFILE2,ERR=992) (VOLSUR(IR),IR=1,NUNKNO)
      IREC= IREC+1
      READ (IFILE1,ERR=991) (MATALB(IR),IR=1,NUNKNO)
      WRITE(IFILE2,ERR=992) (MATALB(IR),IR=1,NUNKNO)
      IREC= IREC+1
      READ (IFILE1,ERR=991) (ICODE(IR),IR=1,NALBG)
      WRITE(IFILE2,ERR=992) (ICODE(IR),IR=1,NALBG)
      IREC= IREC+1
      READ (IFILE1,ERR=991) (ALBEDO(IR),IR=1,NALBG)
      WRITE(IFILE2,ERR=992) (ALBEDO(IR),IR=1,NALBG)
      IREC= IREC+1
      READ (IFILE1,ERR=991) (ANGLES(IR),IR=1,NDIM*NANGL)
      WRITE(IFILE2,ERR=992) (ANGLES(IR),IR=1,NDIM*NANGL)
      IREC= IREC+1
      READ (IFILE1,ERR=991) (DENSTY(IR),IR=1,NANGL)
      WRITE(IFILE2,ERR=992) (DENSTY(IR),IR=1,NANGL)
*----
*  NOW, COPY ALL TRACKS
*----
   60 CONTINUE
         IREC= IREC + 1
         READ (IFILE1,END=70,ERR=991) NSUB,LINE,WEIGHT,
     >                      (KANGL(IR),IR=1,NSUB),
     >                      (NRSEG(IR),IR=1,LINE),(SEGLEN(IR),IR=1,LINE)
         IF(NCOR.EQ.1) THEN
            WRITE(IFILE2,ERR=992) NSUB,LINE,WEIGHT,
     >                  (KANGL(IR),IR=1,NSUB),
     >                  (NRSEG(IR),IR=1,LINE),(SEGLEN(IR),IR=1,LINE)
         ELSE
           I1=1
           DO IR=2,NCOR
             IF(NRSEG(IR).NE.NRSEG(1)) I1=NCOR
           ENDDO
           I2=1
           DO IR=2,NCOR
             IF(NRSEG(LINE-NCOR+IR).NE.NRSEG(LINE-NCOR+1)) I2=NCOR
           ENDDO
           DO IS=1,I1
             DO JS=1,I2
               WEIGHT2=WEIGHT
               IF(I1.GT.1) WEIGHT2=SEGLEN(IS)*WEIGHT2
               IF(I2.GT.1) WEIGHT2=SEGLEN(LINE-NCOR+JS)*WEIGHT2
               ISURF=NRSEG(IS)
               JSURF=NRSEG(LINE-NCOR+JS)
               WRITE(IFILE2,ERR=992) NSUB,LINE-2*NCOR+2,WEIGHT2,
     >                  (KANGL(IR),IR=1,NSUB),
     >                  ISURF,(NRSEG(IR+1),IR=NCOR,LINE-NCOR-1),JSURF,
     >                  1.0D0,(SEGLEN(IR+1),IR=NCOR,LINE-NCOR-1),1.0D0
             ENDDO
           ENDDO
         ENDIF
      GO TO 60
   70 CONTINUE
*----
*  RELEASE TEMPORARY SPACE AND REWIND BOTH FILES
*----
      DEALLOCATE(KANGL,SEGLEN,DENSTY,ANGLES,ALBEDO,VOLSUR)
      DEALLOCATE(NRSEG,ICODE,MATALB)
      REWIND IFILE1
      REWIND IFILE2
      RETURN
*
  991 WRITE(IOUT,'(30H ERROR= RECORD DESTROYED...    )')
      WRITE(IOUT,'(31H ERROR= UNABLE TO READ  RECORD ,I10)') IREC
      WRITE(IOUT,'(31H ERROR=              ON FILE FT,I2.2)') IFILE1
      CALL XABORT( 'XELCOR: --- READ  TRACKING FILE FAILED' )
  992 WRITE(IOUT,'(30H ERROR= NOT ENOUGH SPACE...    )')
      WRITE(IOUT,'(31H ERROR= UNABLE TO WRITE RECORD ,I10)') IREC
      WRITE(IOUT,'(31H ERROR=              ON FILE FT,I2.2)') IFILE1
      CALL XABORT( 'XELCOR: --- WRITE TRACKING FILE FAILED' )
      END
