*DECK TRIZNR
      SUBROUTINE TRIZNR(IMPX,ICOTE,CENTER,CELEM,IAXIS,NR0,RR0,XR0,ANG,
     1 QFR,QTR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculates the correcting factor for cylinder outside elements.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
*  IMPX   print parameter (equal to zero for no print).
*  ICOTE  face number (1=X-, 2=X+, 3=Y-, 4=Y+, 5=Z-, 6=Z+).
*  CENTER coordinates for center of cylinder.
*  CELEM  coordinates for center of the element.
*  IAXIS  principal axis for cylinder.
*  NR0    number of radii.
*  RR0    radii.
*  XR0    coordinates on principal axis.
*  ANG    angles for applying circular correction.
*
*Parameters: output
*  QFR    used to compute transmission factor ( K0/COST ).
*  QTR    used to compute transmission factor ( K0*(R0-RELEM) ).
*
*-----------------------------------------------------------------------
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,ICOTE,IAXIS,NR0
      REAL CENTER(3),CELEM(3),RR0(NR0),XR0(NR0),ANG(NR0),QFR,QTR
*----
*  LOCAL VARIABLES
*----
      CHARACTER*4 AXE(6)
      PARAMETER ( PI= 3.1415926535, PIO2= 0.5*PI, EPSERR=0.05  )
      DATA AXE/ '1=X-', '2=X+', '3=Y-', '4=Y+', '5=Z-', '6=Z+'/
      R0 = RR0(1)
      TET0= 0.0
      DO 5 IR = 1, NR0
         IF(CELEM(IAXIS).GT.XR0(IR)) THEN
            R0 = RR0(IR)
            TET0= ANG(IR)
         ENDIF
    5 CONTINUE
      PITET0 = PIO2 - TET0
*
      IC = (ICOTE+1)/2
      IF(IC.EQ.IAXIS) CALL XABORT('TRIZNR: NOT POSSIBLE TO PROJECT CYL'
     1 //'INDERS ON THAT AXIS.')
*
*     FIND THE ANGLE OF THE ELEMENT
      IX= MOD(IAXIS  ,3) + 1
      IY= MOD(IAXIS+1,3) + 1
      THETA = ABS(ATAN2( CELEM(IX)-CENTER(IX), CELEM(IY)-CENTER(IY) ))
      IF( THETA.LE.PITET0 )THEN
*        NO  CORRECTION
         QFR  = 1.0
         QTR  = 0.0
      ELSE
*        CIRCULAR BOUNDARY CONDITION IS APPLIED
         JC1 = 0
         CCFR0  = 0.0
         RELEM  = 0.0
         DO 10 JC= 1, 3
            IF( JC.NE.IAXIS ) THEN
*           CALCULATE THE RADIUS OF THE ELEMENT
               RELEM= RELEM + (CENTER(JC)-CELEM(JC))**2
*           CALCULATE THE DISTANCE BETWEEN THE "JC" COORDINATES OF THE
*           CENTER OF THE CYLINDER AND OF ACTUAL CYLINDRICAL BOUNDARY
*           IN THE IC DIRECTION
               IF( JC.NE.IC ) THEN
                  JC1 = 2*JC
                  CCFR0= (CENTER(JC) - CELEM(JC))
               ENDIF
            ENDIF
   10    CONTINUE
         RELEM= SQRT(RELEM)
         IF((IMPX.GT.0).AND.(ABS((RELEM-R0)/R0).GT.EPSERR)) THEN
            WRITE(6,1001) CELEM, THETA, RELEM, R0
         ENDIF
*
*        THEN, CALCULATE
*         THE DISTANCE BETWEEN THE CENTER OF THE BOUNDARY
*         ELEMENT AND THE ACTUAL BOUNDARY IN THE IC DIRECTION (DELT)
*        AND
*         THE DIRECTION COSINE OF THE OUTWARD DIRECTED
*         NORMAL AT THE ACTUAL BOUNDARY IN THE IC DIRECTION (COST)
*
         DELT = (R0*R0-CCFR0*CCFR0)
         IF( DELT.LT.0.0)THEN
            JC = JC1/2
            WRITE(6,'(7H ICOTE=,I4,7H IAXIS=,I4)') ICOTE,IAXIS
            WRITE(6,2001) AXE(JC1), CELEM(JC), AXE(JC1), CELEM(IC),
     >                 AXE(ICOTE), R0, DELT, AXE(ICOTE), CELEM(IAXIS)
            WRITE(6,2002)
            DO 20 IR=1, NR0
               WRITE(6,2003) IR, XR0(IR), RR0(IR)
   20       CONTINUE
            CALL XABORT('TRIZNR: ALGORITHM FAILURE.')
         ENDIF
         DELT = SQRT(DELT)
         COST = DELT / R0
         DELT = DELT - ABS( CELEM(IC)-CENTER(IC) )
*
         QFR  = COST
         QTR  = DELT*COST
      ENDIF
      RETURN
*
 1001 FORMAT(  1X,' SURFACE POINT:', 3E11.4,' ANGLE: ', F6.4,
     >            ' RAYON ELEMENT:',  E11.4,' CYLINDRE: ', E11.4 )
 2001 FORMAT( /1X,'*** ERREUR / REACTEUR CYLINDRIQUE ***'/ 5X,
     >' LA COTE SUR L AXE ',A4,' DE L ELEMENT SITUE A',E15.6,
     >' (AXE ',A4,') ET',E15.6,' (AXE ',A4,')'/5X,' EST ',
     >'SUPERIEURE AU RAYON DU CYLINDRE (R0 = ',E15.6,')'/
     > 5X,' DISTANCE (DELT) :',E15.6,' A LA FRONTIERE SUR L AXE ',A4/
     > 5X,' VALEUR EN ALTITUDE:',E15.6/
     >1X,'*** IMPOSSIBLE - ARRET DE L EXECUTION ***')
 2002 FORMAT( /1X,'*** ON DONNE LES ALTITUDES ET LES RAYONS'/
     >'  NREG         Z(NREG)      R(NREG)'/)
 2003 FORMAT(1X,I4,2X,E15.6,2X,E15.6)
      END
