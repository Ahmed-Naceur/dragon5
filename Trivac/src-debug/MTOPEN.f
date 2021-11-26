*DECK MTOPEN
      SUBROUTINE MTOPEN(IMPX,IPTRK,LL4)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Examine and print information related to the automatic matrix
* processor (MTLDLS and MTLDLM).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IMPX    print parameter (equal to zero for no print).
* IPTRK   L_TRACK pointer to the tracking information.
* LL4     order of the coefficient matrix.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER IMPX,LL4
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER HSMG*90,CMODUL*12
      LOGICAL LMU,LMUW,LMUX,LMUY,LMUZ
      INTEGER ITP(NSTATE)
*
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,CMODUL)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ITP)
      ICHX=0
      NLF=0
      ISEG=0
      IF(CMODUL.EQ.'BIVAC') THEN
         NLF=ITP(14)
         ISEG=ITP(17)
      ELSE IF(CMODUL.EQ.'TRIVAC') THEN
         ICHX=ITP(12)
         ISEG=ITP(17)
         NLF=ITP(30)
      ENDIF
      IF((IMPX.GT.0).AND.(ISEG.GT.0)) THEN
         IMPV=ITP(18)
         LTSW=ITP(19)
         WRITE(6,'(9X,36HSUPERVECTORIZATION OPTION ON. ISEG =,I4,
     1   8H  IMPV =,I3,8H  LTSW =,I3)') ISEG,IMPV,LTSW
      ENDIF
      CALL LCMLEN(IPTRK,'MU',LL40,ITYLCM)
      LMU=(LL40.NE.0).AND.(ITYLCM.EQ.1)
      CALL LCMLEN(IPTRK,'MUW',LL4W,ITYLCM)
      LMUW=(LL4W.NE.0).AND.(ITYLCM.EQ.1)
      CALL LCMLEN(IPTRK,'MUX',LL4X,ITYLCM)
      LMUX=(LL4X.NE.0).AND.(ITYLCM.EQ.1)
      CALL LCMLEN(IPTRK,'MUY',LL4Y,ITYLCM)
      LMUY=(LL4Y.NE.0).AND.(ITYLCM.EQ.1)
      CALL LCMLEN(IPTRK,'MUZ',LL4Z,ITYLCM)
      LMUZ=(LL4Z.NE.0).AND.(ITYLCM.EQ.1)
      IDIM=1
      IF(LMU) THEN
         LL4TST=LL40
         HSMG='INVERSE POWER METHOD.'
      ELSE IF(LMUW) THEN
         IDIM=2
         IF((.NOT.LMUX).OR.(.NOT.LMUY)) CALL XABORT('MTOPEN: X- OR Y-C'
     1   //'OMPONENT MISSING IN HEXAGONAL GEOMETRY CASE.')
         IF(LMUZ) IDIM=3
         CALL LCMLEN(IPTRK,'IPW',LL4TST,ITYLCM)
         IF(ICHX.EQ.2) LL4TST=ITP(25)+LL4W+LL4X+LL4Y+LL4Z
         WRITE(HSMG,'(I1,33H-AXIS HEXAGONAL ADI POWER METHOD.)') IDIM+1
      ELSE IF(LMUX) THEN
         IF(LMUY) IDIM=2
         IF(LMUZ) IDIM=3
         CALL LCMLEN(IPTRK,'IPX',LL4TST,ITYLCM)
         IF(ICHX.EQ.2) LL4TST=ITP(25)+LL4W+LL4X+LL4Y+LL4Z
         WRITE(HSMG,'(I1,33H-AXIS CARTESIAN ADI POWER METHOD.)') IDIM
      ELSE IF(LMUY) THEN
         IDIM=2
         IF(LMUZ) IDIM=3
         CALL LCMLEN(IPTRK,'IPY',LL4TST,ITYLCM)
         IF(ICHX.EQ.2) LL4TST=ITP(25)+LL4W+2*LL4Y+LL4Z
         WRITE(HSMG,'(I1,42H-AXIS CARTESIAN ADI POWER METHOD (DIAGONAL,
     1   10H SYMMETRY))') IDIM
      ELSE
         CALL XABORT('MTOPEN: MISSING MU INFO ON LCM.')
      ENDIF
*
      IF(NLF.GT.0) LL4TST=LL4TST*NLF/2
      IF(LL4TST.LE.0) CALL XABORT('MTOPEN: UNABLE TO FIND THE NUMBER O'
     1 //'F UNKNOWNS.')
      IF(IMPX.GT.0) WRITE(6,'(/29H MTOPEN: NUMBER OF UNKNOWNS =,I8,
     1 2H. ,A90)') LL4TST,HSMG
      IF(LL4TST.NE.LL4) CALL XABORT('MTOPEN: INVALID NB OF UNKNOWNS.')
      RETURN
      END
