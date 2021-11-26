*DECK LELCHK
      LOGICAL FUNCTION LELCHK(  NSOLD,  NVOLD, VOLOLD, MATOLD,  ICOLD,
     >                          NSNEW,  NVNEW, VOLNEW, MATNEW,  ICNEW,
     >                           IPRT )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Check compatibility between data in the old tracking file and
* in the new geometry.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NSOLD   number of surfaces in tracking file.
* NVOLD   number of zones in tracking file.
* VOLOLD  volumes and surfaces in tracking file.
* MATOLD  numbering of surfaces and zones in tracking file.
* ICOLD   index of B.C. in tracking file.
* NSNEW   number of surfaces in new geometry. 
* NVNEW   number of zones in new geometry. 
* VOLNEW  volumes & surfaces in new geometry. 
* MATNEW  numbering of surfaces and zones in new geometry. 
* ICNEW   index of B.C. in new geometry. 
* IPRT    printing level ( 0: no print)              
*
*Parameters: output
* LELCHK  checking flag: =.true. if everything was compatible
*         =.false. if incompatibility were detected.
*
*-----------------------------------------------------------------------
*
      IMPLICIT    NONE
*
      INTEGER     NSOLD,NVOLD,MATOLD(-NSOLD:NVOLD),ICOLD(6),IPRT,IOUT,
     >            NSNEW,NVNEW,MATNEW(-NSOLD:NVOLD),ICNEW(6),IR,NERROC
      REAL        VOLOLD(-NSOLD:NVOLD),VOLNEW(-NSNEW:NVNEW),
     >            ZERO,HUND,EMAX
      PARAMETER ( IOUT=6, ZERO=0.0, HUND=100.0, EMAX=1.E-5 )
      LELCHK= .TRUE.
*
*1.1) CHECK # OF ZONES ------------------------------------------------
      IF( NVOLD.NE.NVNEW )THEN
         IF( IPRT.GT.0 )THEN
            WRITE(IOUT,'(40H *** INCONSISTENT # OF ZONES            )')
         ENDIF
         LELCHK=.FALSE.
         GO TO 999
      ENDIF
*
*1.2) CHECK # OF FACES ------------------------------------------------
      IF( NSOLD.NE.NSNEW )THEN
         IF( IPRT.GT.0 )THEN
            WRITE(IOUT,'(40H *** INCONSISTENT # OF FACES            )')
         ENDIF
         LELCHK=.FALSE.
         GO TO 999
      ENDIF
*
*1.3) CHECK CONSISTENCY OF INDEX *ICODE* ------------------------------
      DO 10 IR= 1, 6
         IF( ICOLD(IR).NE.ICNEW(IR) )THEN
            IF( IPRT.GT.0 )THEN
               WRITE(IOUT,'(9H   ICODE(,I1,3H)= ,I6,5H(WAS ,I6,1H))')
     >                                  IR,      ICNEW(IR), ICOLD(IR)
            ENDIF
            IF( ICOLD(IR).LE.0.OR.ICNEW(IR).LE.0 )THEN
               LELCHK=.FALSE.
               GO TO 999
            ENDIF
         ENDIF
   10 CONTINUE
*
*1.4) CHECK IF SOME FACES HAVE ICODE=0 --------------------------------
      DO 20 IR= -NSOLD, -1
         IF( ICNEW(-MATNEW(IR)).EQ.0 )THEN
            IF( IPRT.GT.0 )THEN
               WRITE(IOUT,'(9H    FACE(,I1,3H)= ,I6,12H HAS ICODE=0 )')
     >                                 -IR,      MATNEW(IR)
            ENDIF
            LELCHK=.FALSE.
            GO TO 999
         ENDIF
   20 CONTINUE
*
*2)   CHECK CONSISTENCY OF VECTORS *VOLSUR* AND *MATALB* --------------
      NERROC= 0
      DO 30 IR= -NSOLD, NVOLD
         IF( VOLOLD(IR)-VOLNEW(IR).GT.ZERO )THEN
            NERROC= NERROC+1
            IF( IR.EQ.0 ) GO TO 30
            LELCHK= LELCHK.AND.
     >              ABS((VOLNEW(IR)-VOLOLD(IR))/VOLOLD(IR)).LE.EMAX
         ENDIF
         IF( MATOLD(IR).NE.MATNEW(IR) )THEN
            NERROC= NERROC+1
            IF( IR.LE.0 ) LELCHK= .FALSE.
         ENDIF
   30 CONTINUE
      IF( IPRT.GT.0 )THEN
         WRITE(IOUT,'(1H )')
         IF( NERROC.EQ.0 )THEN
            WRITE(IOUT,'(60H ECHO = >>> CONSISTENCY BETWEEN '//
     >                 'TRACKING FILE AND GEOMETRY                 /)')
         ELSE
            WRITE(IOUT,'(60H ECHO = >>> WARNING: INCONSISTENT '//
     >                 'TRACKING FILE                              /)')
            DO 40 IR= -NSOLD, NVOLD
               IF( IR.EQ.0 ) GO TO 40
               IF( VOLOLD(IR)-VOLNEW(IR).GT.ZERO )THEN
               IF( IR.LE.0 )THEN
                  WRITE(IOUT,'(15H ERROR ON FACE(,I4,3H)= ,F10.7,1H%)')
     >                      -IR,HUND*(VOLNEW(IR)-VOLOLD(IR))/VOLOLD(IR)
               ELSE
                  WRITE(IOUT,'(15H ERROR ON ZONE(,I4,3H)= ,F10.7,1H%)')
     >                       IR,HUND*(VOLNEW(IR)-VOLOLD(IR))/VOLOLD(IR)
               ENDIF
               ENDIF
               IF( MATOLD(IR).NE.MATNEW(IR) )THEN
               IF( IR.LE.0 )THEN
                  WRITE(IOUT,'(9H    FACE(,I1,3H)= ,I6,5H(WAS ,I6,1H))')
     >                                    -IR,     MATNEW(IR),MATOLD(IR)
               ELSE
                  WRITE(IOUT,'(9H MIXTURE(,I1,3H)= ,I6,5H(WAS ,I6,1H))')
     >                                     IR,     MATNEW(IR),MATOLD(IR)
               ENDIF
               ENDIF
   40       CONTINUE
         ENDIF
      ENDIF
*
  999 RETURN
      END
