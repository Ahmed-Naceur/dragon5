*DECK XELMRG
      SUBROUTINE XELMRG (  IPRT,   NSUR,   NVOL,   NSBC, NTOTCL, INDEX,
     >                   MINDIM, MAXDIM, LCLSYM, LCLTRA,    LL1,   LL2,
     >                   MRGCEL, MATALB, KEYMRG, INCELL, MATRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct keymrg according to implicit merging imposed by the 
* boundary conditions.
*
*Copyright:
* Copyright (C) 1990 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPRT    printing level.                               
* NSUR    number of surfaces.                               
* NVOL    number of zones.                                  
* NSBC    number of surfaces with independent BC.           
* NTOTCL  number of cylindres+3.                            
* INDEX   numbering of surfaces and zones.                  
* MINDIM  minimum index values for all axes (rect/cyl). 
* MAXDIM  maximum index values for all axes (rect/cyl).
* LCLSYM  symmetry flags (0: no; -1/+1: syme; -2/+2: ssym). 
* LCLTRA  translation flags (0: no; -1/+1: tra).           
* LL1     diagonal symmetry (2,3).                     
* LL2     diagonal symmetry (1,4).                      
* MRGCEL  merging cell numbering.                           
* MATALB  material types.                              
*
*Parameters: input/output
* KEYMRG  initial numbering at input, merged at output.     
* INCELL  block   numbering at input, merged at output.     
*
*Parameters: output
* MATRT   reflection/transmission vector.
*
*-----------------------------------------------------------------------
*
      IMPLICIT     NONE
*
      INTEGER      IPRT,   NSUR,   NVOL,   NSBC, NTOTCL
      INTEGER      LCLSYM(3), LCLTRA(3), IORD(4), INDEX(4,*), KEYMRG(*),
     >             MATALB(*), MINDIM(NTOTCL), MAXDIM(NTOTCL), INCELL(*),
     >             MRGCEL(*), MATRT(-NSUR,2)
      LOGICAL      LL1, LL2
*
      LOGICAL      SWOK, SWSUR, SWSTOP
      INTEGER      NUM, I, ISUR, ITRA, IVS1, IVS2, ISYM, IORD4, ICC1,
     >             INDEX4, INCR, NO1, NO2, IB1, IB2, NZSU, NZVO, NZABS,
     >             NMBLK, IBLK, NZBLK, IMRG, MINV, MAXV, NVOLM, NSURM,
     >             ICMP1, ICMP2, ITRAC1, NSURC, IP, IR, NVOLC, NMVO
      CHARACTER*4  CORIEN(-6:0)
      INTEGER      IOUT
      PARAMETER  ( IOUT=6 )
*
      DATA              CORIEN
     >        / ' Z+ ',' Z- ',' Y+ ',' Y- ',' X+ ',' X- ','    ' /
*
      NUM(I)= I + 1 - NSUR
*
*     INITIALIZE MATRT TO REFLECTION FOR ORIGINAL SURFACES
      NVOLM=0
      DO 300 ISUR=1, -NSUR
        MATRT(ISUR,1)=0
        MATRT(ISUR,2)=ISUR
 300  CONTINUE
*
*  0) TREAT TRANSLATION SYMMETRIES *************************************
      DO 310 ITRA =1,3
        IF( LCLTRA(ITRA) .EQ. 1) THEN
          DO 320 IVS1 = NSUR, -1
            IF( (KEYMRG(NUM(IVS1)) .NE.  0   ) .AND.
     >          (MATRT(-IVS1,2)    .EQ. -IVS1)       ) THEN
*
*             LOCATE SURFACE IN X, Y, Z AND R
              IORD(1)= INDEX(1,NUM(IVS1))
              IORD(2)= INDEX(2,NUM(IVS1))
              IORD(3)= INDEX(3,NUM(IVS1))
              IORD(4)= INDEX(4,NUM(IVS1))
*
*             LOCATE TRANSLATED SURFACE IN X, Y, Z AND R
              IF( (IORD(ITRA) .GE. MINDIM(ITRA)) .AND.
     >            (IORD(ITRA) .LT. MAXDIM(ITRA)) ) GO TO 345
              IORD(ITRA)= (MAXDIM(ITRA)+MINDIM(ITRA))-(IORD(ITRA)+1)
*              INDEX(1,NUM(0))= IORD(1)
*              INDEX(2,NUM(0))= IORD(2)
*              INDEX(3,NUM(0))= IORD(3)
*              INDEX(4,NUM(0))= IORD(4)
*
*             FOR CYLINDERS, *IORD4* IS ABSOLUTE.
              IORD4 = IORD(4)
              IF(IORD(4) .NE. 0 )THEN
                DO 330 ICC1= NTOTCL, 4, -1
                  IF( IORD(4) .LT. MAXDIM(ICC1) )THEN
                    IORD4 = IORD(4)-MINDIM(ICC1)
                  ENDIF
  330           CONTINUE
              ENDIF
              DO 340 IVS2=NSUR,-1
                IF( IORD(1) .EQ. INDEX(1,NUM(IVS2)).AND.
     >              IORD(2) .EQ. INDEX(2,NUM(IVS2)).AND.
     >              IORD(3) .EQ. INDEX(3,NUM(IVS2))      ) THEN
                  INDEX4= INDEX(4,NUM(IVS2))
                  IF( IORD(4).NE.0 )THEN
                    DO 350 ICC1= NTOTCL, 4, -1
                      IF( INDEX(4,NUM(IVS2)).LT.MAXDIM(ICC1) )THEN
                        INDEX4= INDEX(4,NUM(IVS2))-MINDIM(ICC1)
                      ENDIF
  350               CONTINUE
                  ENDIF
*
*                 SYMMETRIC SURFACE LOCATED FOR TRANSMISSION BC
*                 STORE SURFACES IDENTIFIER IN MATRT AND
*                 EXIT TO 345
                  IF( INDEX4.EQ.IORD4) THEN
                    MATRT(-IVS1,2)=-IVS2
                    MATRT(-IVS2,2)=-IVS1
                    GO TO 345
                  ENDIF
                ENDIF
 340          CONTINUE
              CALL XABORT( 'XELMRG: TRANSLATED SURFACE NO FOUND.' )
            ENDIF
 345        CONTINUE
 320      CONTINUE
        ENDIF
 310  CONTINUE
*
*  1) TREAT AXIAL SYMMETRIES ******************************************
      DO 20 ISYM= 1, 3
         IF( LCLSYM(ISYM).NE.0 )THEN
            DO 10 IVS1= NSUR, NVOL
*
*              FOR REGIONS ABSENT FROM FINAL CELL
*              DO NOT BOTHER TO SYMMETRIZE
               IF( IVS1 .EQ. 0 .OR. KEYMRG(NUM(IVS1)) .EQ. 0) GO TO 10
               IORD(1)= INDEX(1,NUM(IVS1))
               IORD(2)= INDEX(2,NUM(IVS1))
               IORD(3)= INDEX(3,NUM(IVS1))
               IORD(4)= INDEX(4,NUM(IVS1))
*
*  1.1)        RECOMPOSE *ISYM* VALUE TO  GET THE SYMMETRIC COORDINATE
               IORD(ISYM)= (MAXDIM(ISYM)+MINDIM(ISYM))-(IORD(ISYM)+1)
               IF( IVS1.GT.0 )THEN
                  IVS2= NVOL
                  INCR= -1
               ELSE
                  IVS2= NSUR
                  INCR= +1
               ENDIF
*               INDEX(1,NUM(0))= IORD(1)
*               INDEX(2,NUM(0))= IORD(2)
*               INDEX(3,NUM(0))= IORD(3)
*               INDEX(4,NUM(0))= IORD(4)
*
*  1.2)        TO SEARCH FOR THE GOOD CYLINDER, *IORD4* IS ABSOLUTE.
               IORD4 = IORD(4)
               IF( IORD(4).NE.0 )THEN
                  DO 110 ICC1= NTOTCL, 4, -1
                     IF( IORD(4).LT.MAXDIM(ICC1) )THEN
                        IORD4 = IORD(4)-MINDIM(ICC1)
                     ENDIF
  110             CONTINUE
               ENDIF
   11          CONTINUE
                  IF( IORD(1).EQ.INDEX(1,NUM(IVS2)).AND.
     >                IORD(2).EQ.INDEX(2,NUM(IVS2)).AND.
     >                IORD(3).EQ.INDEX(3,NUM(IVS2)) )THEN
                     INDEX4= INDEX(4,NUM(IVS2))
                     IF( IORD(4).NE.0 )THEN
                        DO 112 ICC1= NTOTCL, 4, -1
                           IF( INDEX(4,NUM(IVS2)).LT.MAXDIM(ICC1) )THEN
                              INDEX4= INDEX(4,NUM(IVS2))-MINDIM(ICC1)
                           ENDIF
  112                   CONTINUE
                     ENDIF
                     IF( INDEX4.EQ.IORD4) GO TO 12
                  ENDIF
                  IVS2= IVS2 + INCR
               GO TO 11
   12          IF( IVS2.EQ.0 )THEN
                  CALL XABORT( 'XELMRG: RARE AXIAL SYMMETRY PROBLEM.' )
               ENDIF
               NO1= KEYMRG(NUM(IVS1))
               NO2= KEYMRG(NUM(IVS2))
               IB1= INCELL(NUM(IVS1))
               IB2= INCELL(NUM(IVS2))
*
*  1.3)        SELECT THE MAX OR MIN VALUE TO CORRECTLY # ZONES
               IF( LCLSYM(ISYM).GT.0 )THEN
                  KEYMRG(NUM(IVS1))= MIN(NO1,NO2)
                  KEYMRG(NUM(IVS2))= MIN(NO1,NO2)
                  INCELL(NUM(IVS1))= MIN(IB1,IB2)
                  INCELL(NUM(IVS2))= MIN(IB1,IB2)
               ELSE
                  KEYMRG(NUM(IVS1))= MAX(NO1,NO2)
                  KEYMRG(NUM(IVS2))= MAX(NO1,NO2)
                  INCELL(NUM(IVS1))= MAX(IB1,IB2)
                  INCELL(NUM(IVS2))= MAX(IB1,IB2)
               ENDIF
   10       CONTINUE
         ENDIF
   20 CONTINUE
*
*  2) TREAT DIAGONAL SYMMETRIES ***************************************
*                       (SIDE #3)
*               (SIDE #1) GEOM  (SIDE #2)
*                       (SIDE #4)
*
      IF( LL1.OR.LL2 )THEN
         DO 30 IVS1= NSUR, NVOL
*
*           FOR REGIONS ABSENT FROM FINAL CELL
*           DO NOT BOTHER TO SYMMETRIZE
            IF( IVS1 .EQ. 0 .OR. KEYMRG(NUM(IVS1)) .EQ. 0 ) GO TO 30
*
*  2.1)        FOR (SIDE #1).EQ.(SIDE #4)
*              AND (SIDE #2).EQ.(SIDE #3) *** DIAGONAL SYMMETRY (\) ***
*           NOTE: ***NOT*** ACCEPTED IN DRAGON.
***         IORD(1)= (MAXDIM(2)+MINDIM(1)) - (INDEX(2,NUM(IVS1))+1)
***         IORD(2)= (MAXDIM(1)+MINDIM(2)) - (INDEX(1,NUM(IVS1))+1)
*  2.2)        FOR (SIDE #2).EQ.(SIDE #4)
*              AND (SIDE #1).EQ.(SIDE #3) *** DIAGONAL SYMMETRY (/) ***
            IORD(1)= INDEX(2,NUM(IVS1)) + MINDIM(1) - MINDIM(2)
            IORD(2)= INDEX(1,NUM(IVS1)) + MINDIM(2) - MINDIM(1)
            IORD(3)= INDEX(3,NUM(IVS1))
            IORD(4)= INDEX(4,NUM(IVS1))
            IF( IVS1.GT.0 )THEN
               IVS2= NVOL
               INCR= -1
            ELSE
               IVS2= NSUR
               INCR= +1
            ENDIF
*            INDEX(1,NUM(0))= IORD(1)
*            INDEX(2,NUM(0))= IORD(2)
*            INDEX(3,NUM(0))= IORD(3)
*            INDEX(4,NUM(0))= IORD(4)
            IORD4 = IORD(4)
            IF( IORD(4).NE.0 )THEN
               DO 33 ICC1= NTOTCL, 4, -1
                  IF( IORD(4).LT.MAXDIM(ICC1) )THEN
                     IORD4 = IORD(4)-MINDIM(ICC1)
                  ENDIF
   33          CONTINUE
            ENDIF
   31       CONTINUE
               IF( IORD(1).EQ.INDEX(1,NUM(IVS2)).AND.
     >             IORD(2).EQ.INDEX(2,NUM(IVS2)).AND.
     >             IORD(3).EQ.INDEX(3,NUM(IVS2)) )THEN
                  INDEX4= INDEX(4,NUM(IVS2))
                  IF( IORD(4).NE.0 )THEN
                     DO 34 ICC1= NTOTCL, 4, -1
                        IF( INDEX(4,NUM(IVS2)).LT.MAXDIM(ICC1) )THEN
                           INDEX4= INDEX(4,NUM(IVS2))-MINDIM(ICC1)
                        ENDIF
   34                CONTINUE
                  ENDIF
                  IF( INDEX4.EQ.IORD4) GO TO 32
               ENDIF
               IVS2= IVS2 + INCR
            GO TO 31
   32       IF( IVS2.EQ.0 )THEN
               CALL XABORT( 'XELMRG: RARE DIAGONAL SYMMETRY PROBLEM.' )
            ENDIF
            NO1= KEYMRG(NUM(IVS1))
            NO2= KEYMRG(NUM(IVS2))
            IB1= INCELL(NUM(IVS1))
            IB2= INCELL(NUM(IVS2))
*
*  2.3)        SELECT THE MAX OR MIN VALUE TO CORRECTLY # ZONES
            IF( LL2 )THEN
               KEYMRG(NUM(IVS1))= MIN(NO1,NO2)
               KEYMRG(NUM(IVS2))= MIN(NO1,NO2)
               INCELL(NUM(IVS1))= MIN(IB1,IB2)
               INCELL(NUM(IVS2))= MIN(IB1,IB2)
            ELSE
               KEYMRG(NUM(IVS1))= MAX(NO1,NO2)
               KEYMRG(NUM(IVS2))= MAX(NO1,NO2)
               INCELL(NUM(IVS1))= MAX(IB1,IB2)
               INCELL(NUM(IVS2))= MAX(IB1,IB2)
            ENDIF
   30    CONTINUE
      ENDIF
*
*  3) NOW, STOCK NEW INCREASING VALUES IN *KEYMRG* AND *INCELL* ********
      NZSU= 0
      DO 40 IVS1= -1, NSUR,-1
        DO 41 IVS2= -1, NSUR, -1
*
*  3.1.1) COUNT THE NUMBER OF SURFACES.
          IF( KEYMRG(NUM(IVS2)).EQ.IVS1 )THEN
            NZSU= NZSU-1
            GO TO 40
          ENDIF
   41   CONTINUE
   40 CONTINUE
      NZVO=0
      DO 42 IVS1= 1, NVOL
        DO 43 IVS2= 1, NVOL
*
*  3.1.2) COUNT THE NUMBER OF VOLUMES.
          IF( KEYMRG(NUM(IVS2)).EQ.IVS1 )THEN
            NZVO= NZVO+1
            GO TO 42
          ENDIF
   43   CONTINUE
   42 CONTINUE
      NZABS= -1
      DO 50 IVS1= -1, NSUR, -1
        SWOK= .FALSE.
        DO 51 IVS2= -1, NSUR, -1
*
*  3.2.1) RENUMBER SURFACES.
          IF( KEYMRG(NUM(IVS2)).EQ.IVS1 )THEN
            SWOK= .TRUE.
            KEYMRG(NUM(IVS2))= NZABS
          ENDIF
   51   CONTINUE
        IF( SWOK )THEN
          NZABS= NZABS - 1
        ENDIF
   50 CONTINUE
      IF( NZABS.NE.NZSU-1 )THEN
         CALL XABORT( 'XELMRG: PROBLEMS TO MERGE SURFACES' )
      ENDIF
      KEYMRG(NUM(0))= 0
      NZABS=  1
      DO 52 IVS1= 1, NVOL
        SWOK= .FALSE.
        DO 53 IVS2= 1, NVOL
*
*  3.2.2) RENUMBER VOLUMES.
          IF( KEYMRG(NUM(IVS2)).EQ.IVS1 )THEN
            SWOK= .TRUE.
            KEYMRG(NUM(IVS2))= NZABS
          ENDIF
   53   CONTINUE
        IF( SWOK )THEN
          NZABS= NZABS + 1
        ENDIF
   52 CONTINUE
      IF( NZABS.NE.NZVO+1 )THEN
         CALL XABORT( 'XELMRG: PROBLEMS TO MERGE VOLUMES' )
      ENDIF
      NMBLK= 0
      DO 60 IVS2= NSUR, NVOL
*
*  3.3) COUNT NUMBER OF BLOCKS.
        IF( KEYMRG(NUM(IVS2)).NE.0 )THEN
          IBLK=INCELL(NUM(IVS2))
          IF( IBLK.NE.0 )THEN
            NMBLK=MAX(NMBLK,IBLK)
          ENDIF
        ENDIF
   60 CONTINUE
      NZBLK= 1
      DO 70 IBLK= 1, NMBLK
        SWOK= .FALSE.
        DO 71 IVS2= NSUR, NVOL
*
*  3.4)   RENUMBER BLOCKS.
          IF( KEYMRG(NUM(IVS2)).NE.0 )THEN
            IF( INCELL(NUM(IVS2)).EQ.IBLK )THEN
              SWOK= .TRUE.
              INCELL(NUM(IVS2))= NZBLK
            ENDIF
          ENDIF
   71   CONTINUE
        IF( SWOK )THEN
          NZBLK= NZBLK + 1
        ENDIF
   70 CONTINUE
      NZBLK= NZBLK-1
      IF( NZBLK .LE. 0 .OR. NZBLK .GT. NMBLK)THEN
         CALL XABORT( 'XELMRG: PROBLEMS TO MERGE BLOCKS' )
      ENDIF
*
*  3.5) RENUMBER CELL BLOCKS ACCORDING TO THE MERGE INDEX *MRGCEL*
*     *** THIS WILL RENUMBER VOLUMES, BUT NOT SURFACES. ***
      NMVO = 0
      NMBLK= 1
      SWSTOP= .FALSE.
      DO 290 IMRG= 1, NZBLK
        SWOK= .FALSE.
        DO 280 IBLK=1, NZBLK
          IF( MRGCEL(IBLK).EQ.IMRG )THEN
            IF( IMRG.NE.NMBLK )
     >         CALL XABORT('XELMRG: INCREASING MERGE #ING REQUIRED')
            MINV= +100000000
            MAXV= 0
            SWSUR=.FALSE.
            DO 210 IVS1= 1, NVOL
              IF(KEYMRG(NUM(IVS1)).GT.0) THEN
                IF( INCELL(NUM(IVS1)).EQ.IBLK )THEN
                  MINV= MIN(MINV,KEYMRG(NUM(IVS1)))
                  MAXV= MAX(MAXV,KEYMRG(NUM(IVS1)))
                ENDIF
              ENDIF
  210       CONTINUE
            IF( SWOK )THEN
              SWSTOP= SWSTOP.OR.(NVOLM.NE.MAXV+1-MINV)
            ELSE
              NVOLM = MAXV+1-MINV
            ENDIF
            MINV= MINV-NMVO
            DO 220 IVS1= 1, NVOL
              IF(KEYMRG(NUM(IVS1)).GT.0) THEN
                IF( INCELL(NUM(IVS1)).EQ.IBLK )THEN
                  KEYMRG(NUM(IVS1))= KEYMRG(NUM(IVS1))-(MINV-1)
                ENDIF
              ENDIF
  220       CONTINUE
            SWOK= .TRUE.
          ENDIF
  280   CONTINUE
        IF( SWOK )THEN
          NMVO= NMVO+NVOLM
          NMBLK= NMBLK+1
        ENDIF
  290 CONTINUE
      NMBLK= NMBLK-1
      NSBC=-NZSU
*
*  4) RESET *MATRT* FOR MERGED SURFACES INSTEAD OF ORIGINAL SURFACES ***
      NZSU=0
      DO 360 IVS1=-1,NSUR,-1
        ICMP1=KEYMRG(NUM(IVS1))
        IVS2=-MATRT(-IVS1,2)
        ICMP2=KEYMRG(NUM(IVS2))
        IF( (ICMP1 .LT. 0) .AND. (ICMP2 .LT. 0) ) THEN
          ITRAC1=MATRT(-ICMP1,1)
          IF(ITRAC1 .EQ. 0) THEN
            MATRT(-ICMP1,1)=-ICMP2
            MATRT(-ICMP2,1)=-ICMP1
          ENDIF
        ENDIF
 360  CONTINUE
*
*  5) PRINTING *********************************************************
      IF( IPRT.GT.2 )THEN
         NSURC = -1
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(/40H  SURFACE #ING ( BEFORE CELL MERGE )     )')
         DO 580 IP  = 1, (9 - NSUR) / 10
            NSURM= MAX( NSUR, NSURC-9 )
            WRITE(IOUT,'(10X,10(A5,I7))')
     >           (' SUR ',-IR,IR= NSURC, NSURM, -1)
            WRITE(IOUT,'(8H ORIENT ,2X,10A12)')
     >           (CORIEN(MATALB(NUM(IR))),IR=NSURC,NSURM,-1)
            WRITE(IOUT,'(8H CELL # ,2X,10I12)')
     >           (INCELL(NUM(IR)),IR=NSURC,NSURM,-1)
            WRITE(IOUT,'(9H MERGE TO ,1X,10(A5,I7))')
     >           (' SUR ',-KEYMRG(NUM(IR)),IR=NSURC,NSURM,-1)
            WRITE(IOUT,'(1H )')
            NSURC = NSURC - 10
  580    CONTINUE
         NVOLC= 1
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'( 40H  VOLUME  #ING ( BEFORE CELL MERGE )      )')
         DO 590 IP  = 1, (9 + NVOL) / 10
            NVOLM= MIN( NVOL, NVOLC+9 )
            WRITE(IOUT,'(10X,10(A5,I7))')
     >               (' VOL ',IR,IR=NVOLC,NVOLM, 1)
            WRITE(IOUT,'(8H CELL # ,2X,10I12)')
     >           (INCELL(NUM(IR)),IR=NVOLC,NVOLM, 1)
            WRITE(IOUT,'(9H MERGE TO ,1X,10(A5,I7))')
     >           (' VOL ', KEYMRG(NUM(IR)),IR=NVOLC,NVOLM, 1)
            WRITE(IOUT,'(8H MIX    ,2X,10I12)')
     >           (MATALB(NUM(IR)),IR=NVOLC,NVOLM, 1)
            WRITE(IOUT,'(9H          ,1X,10(A5,I7))')
     >           (' CELL',MRGCEL(INCELL(NUM(IR))),IR=NVOLC,NVOLM, 1)
            WRITE(IOUT,'(1H )')
            NVOLC = NVOLC + 10
  590    CONTINUE
         WRITE(IOUT,'( 40H  BC MATRIX (BEFORE MERGE)                )')
         WRITE(IOUT,'(8(5X,I10,I10))') (IR,MATRT(IR,2),IR=1,-NSUR)
         WRITE(IOUT,'( 40H  BC MATRIX (AFTER MERGE)                 )')
         WRITE(IOUT,'(8(5X,I10,I10))') (IR,MATRT(IR,1),IR=1,NSBC)
      ELSEIF( IPRT.GT.1 )THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(32H # OF SURFACES AFTER SYMMETRIES: ,I8)') -NZSU
         WRITE(IOUT,'(32H # OF ZONES    AFTER SYMMETRIES: ,I8)')  NZVO
         WRITE(IOUT,'(32H # OF CELLS    AFTER SYMMETRIES: ,I8)')  NZBLK
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(32H # OF ZONES    AFTER MERGING   : ,I8)')  NMVO
         WRITE(IOUT,'(32H # OF CELLS    AFTER MERGING   : ,I8)')  NMBLK
         WRITE(IOUT,'(1H )')
      ENDIF
      IF( SWSTOP )THEN
         CALL XABORT('XELMRG: MERGE CELL ONLY WITH SAME # OF ZONES')
      ENDIF
*
      RETURN
      END
