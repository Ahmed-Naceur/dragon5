*DECK SPHGAP
      SUBROUTINE SPHGAP(IPTRK2,IPRINT,NREG,NUN,MAT,KEY,FUNKNO,COUGAP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the average flux at the boundary 
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modIFy it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Chambon
*
*Parameters: input
* IPTRK2  pointer to the TRIVAC tracking of the macro-geometry
*         (L_TRACK signature).
* IPRINT  print flag (equal to 0 for no print).
* NREG    number of macro-regions (in the macro calculation).
* NUN     number of unknowns in the macro-calculation.
* MAT     mixture index per macro-region.
* KEY     position of the flux components associated with each volume.
* FUNKNO  neutron flux.
*
*Parameters: output
* COUGAP  boundary average flux.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK2
      INTEGER NREG,NUN,MAT(NREG),KEY(NREG)
      REAL FUNKNO(NUN),COUGAP
*----
*  LOCAL VARIABLES
*----
      INTEGER     NSTATE,NP
      PARAMETER  (NSTATE=40,NP=3)
      INTEGER     IPAR(NSTATE),IELEM,NCODE(6),LX,LY,LP,ITYPE,IPRINT,
     +            ICHX,IDIM,LC,L4,MAXKN,MKN,ITYLCM
      INTEGER     I,J
      REAL        XM(1),XP(1),YM(1),YP(1),SXM,SXP,SYM,SYP,DG,LTOT,FACT
      REAL        E(25)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KN
      REAL, ALLOCATABLE, DIMENSION(:) :: X,Y,XX,YY,XXX,YYY,AXYZ
*----
*  RECOVER TRIVAC SPECIFIC TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK2,'STATE-VECTOR',IPAR)
      ITYPE=IPAR(6)
      IF(ITYPE.NE.5) CALL XABORT('SPHGAP: 2D Cartesian geometry '
     1  //'expected')
      IELEM=IPAR(9)
      L4=IPAR(11)
      ICHX=IPAR(12)
      LX=IPAR(14)
      LY=IPAR(15)
      LP=MAX(LX,LY)*NP
      IDIM=2
      ALLOCATE(XX(LX*LY),YY(LX*LY),XXX(LX+1),YYY(LY+1))
      ALLOCATE(X(3*LX),Y(3*LY),AXYZ(LP))
      CALL LCMGET(IPTRK2,'XX',XX)
      CALL LCMGET(IPTRK2,'YY',YY)
      CALL LCMGET(IPTRK2,'NCODE',NCODE)
*----
*  Compute the coordinate of the point on the boundary
*----
      XXX(1)=0.0
      DO 10 I=1,LX
      XXX(I+1)=XXX(I)+XX(I)
   10 CONTINUE
      YYY(1)=0.0
      DO 20 I=1,LY
      YYY(I+1)=YYY(I)+YY((I-1)*LX+1)
   20 CONTINUE
      IF(NCODE(1).EQ.5)THEN
        XM(1)=(XXX(1)+XXX(2))/2.0
      ELSE
        XM(1)=XXX(1)
      ENDIF
      IF(NCODE(2).EQ.5)THEN
        XP(1)=(XXX(LX+1)+XXX(LX))/2.0
      ELSE
        XP(1)=XXX(LX+1)
      ENDIF
      IF(NCODE(3).EQ.5)THEN
        YM(1)=(YYY(1)+YYY(2))/2.0
      ELSE
        YM(1)=YYY(1)
      ENDIF
      IF(NCODE(4).EQ.5)THEN
        YP(1)=(YYY(LY+1)+YYY(LY))/2.0
      ELSE
        YP(1)=YYY(LY+1)
      ENDIF
      DO I=1,NP
        FACT=REAL(2*I-1)/REAL(2*NP)
        X(I)=(XXX(2)-XM(1))*FACT
        X(NP*(LX-1)+I)=XXX(LX)+(XP(1)-XXX(LX))*FACT
        DO J=1,LX-2
          X(NP*J+I)=XXX(J+1)+(XXX(J+2)-XXX(J+1))*FACT
        ENDDO
        Y(I)=(YYY(2)-YM(1))*FACT
        Y(NP*(LY-1)+I)=YYY(LY)+(YP(1)-YYY(LY))*FACT
        DO J=1,LY-2
          Y(NP*J+I)=YYY(J+1)+(YYY(J+2)-YYY(J+1))*FACT
        ENDDO
      ENDDO
      IF(IPRINT.GE.100) then
        WRITE(6,*)'FUNKNO: ='
        do I=1,LY
          WRITE(6,*) I,'#',(FUNKNO(KEY(J+(I-1)*LX)),J=1,LX)
        enddo  
        WRITE(6,*)'NCODE: =',(NCODE(I),I=1,6)  
        WRITE(6,*)'XXX: =',(XXX(I),I=1,LX+1)  
        WRITE(6,*)'YYY: =',(YYY(I),I=1,LY+1)  
        WRITE(6,*)'X: =',(X(I),I=1,3*LX)  
        WRITE(6,*)'Y: =',(Y(I),I=1,3*LY)  
      endIF
*----
*  Interpolate the flux
*----
      COUGAP=0.0
      LTOT=0.0
      IF(NCODE(1).EQ.5)THEN
        SXM=-1.0
      ELSE
        IF(ICHX.EQ.1) THEN
*         Variational collocation method
          CALL LCMLEN(IPTRK2,'KN',MAXKN,ITYLCM)
          MKN=MAXKN/(LX*LY)
          ALLOCATE(KN(MAXKN))
          CALL LCMGET(IPTRK2,'KN',KN)
          CALL LCMSIX(IPTRK2,'BIVCOL',1)
          CALL LCMLEN(IPTRK2,'T',LC,ITYLCM)
          CALL LCMGET(IPTRK2,'E',E)
          CALL LCMSIX(IPTRK2,' ',2)
          CALL VALU2B(LC,MKN,LX,LY,L4,XM,Y,XXX,YYY,FUNKNO,MAT,KN,
     +    1,3*LY,E,AXYZ)
        ELSE IF(ICHX.EQ.2) THEN
*         Raviart-Thomas finite element method
          CALL VALU4B(IELEM,NUN,LX,LY,XM,Y,XXX,YYY,FUNKNO,MAT,KEY,
     +    1,3*LY,AXYZ)
        ELSE IF(ICHX.EQ.3) THEN
*         Nodal collocation method (MCFD)
          CALL VALU1B(IDIM,LX,LY,L4,XM,Y,XXX,YYY,FUNKNO,MAT,IELEM,
     +    1,3*LY,AXYZ)
        ELSE
          CALL XABORT('SPHGAP: INTERPOLATION NOT IMPLEMENTED(1).')
        ENDIF
        IF(IPRINT.GE.100) WRITE(6,*)'SPHGAP: AXYZ =',(AXYZ(I),I=1,3*LY)
        SXM=0.0
        DO J=1,LY
          DG=(MIN(YP(1),YYY(J+1))-MAX(YM(1),YYY(J)))/REAL(NP)
          DO I=1,NP
            SXM=SXM+AXYZ((J-1)*NP+I)*DG
          ENDDO
        ENDDO
        COUGAP=COUGAP+SXM
        LTOT=LTOT+YP(1)-YM(1)
      ENDIF
      IF(NCODE(2).EQ.5)THEN
        SXP=-1.0
      ELSE
        IF(ICHX.EQ.1) THEN
*         Variational collocation method
          CALL LCMLEN(IPTRK2,'KN',MAXKN,ITYLCM)
          MKN=MAXKN/(LX*LY)
          ALLOCATE(KN(MAXKN))
          CALL LCMGET(IPTRK2,'KN',KN)
          CALL LCMSIX(IPTRK2,'BIVCOL',1)
          CALL LCMLEN(IPTRK2,'T',LC,ITYLCM)
          CALL LCMGET(IPTRK2,'E',E)
          CALL LCMSIX(IPTRK2,' ',2)
          CALL VALU2B(LC,MKN,LX,LY,L4,XP,Y,XXX,YYY,FUNKNO,MAT,KN,
     +    1,3*LY,E,AXYZ)
        ELSE IF(ICHX.EQ.2) THEN
*         Raviart-Thomas finite element method
          CALL VALU4B(IELEM,NUN,LX,LY,XP,Y,XXX,YYY,FUNKNO,MAT,KEY,
     +    1,3*LY,AXYZ)
        ELSE IF(ICHX.EQ.3) THEN
*         Nodal collocation method (MCFD)
          CALL VALU1B(IDIM,LX,LY,L4,XP,Y,XXX,YYY,FUNKNO,MAT,IELEM,
     +    1,3*LY,AXYZ)
        ELSE
          CALL XABORT('SPHGAP: INTERPOLATION NOT IMPLEMENTED(2).')
        ENDIF
        IF(IPRINT.GE.100) WRITE(6,*)'SPHGAP: AXYZ =',(AXYZ(I),I=1,3*LY)
        SXP=0.0
        DO J=1,LY
          DG=(MIN(YP(1),YYY(J+1))-MAX(YM(1),YYY(J)))/REAL(NP)
          DO I=1,NP
            SXP=SXP+AXYZ((J-1)*NP+I)*DG
          ENDDO
        ENDDO
        COUGAP=COUGAP+SXP
        LTOT=LTOT+YP(1)-YM(1)
      ENDIF
      IF(NCODE(3).EQ.5)THEN
        SYM=-1.0
      ELSE
        IF(ICHX.EQ.1) THEN
*         Variational collocation method
          CALL LCMLEN(IPTRK2,'KN',MAXKN,ITYLCM)
          MKN=MAXKN/(LX*LY)
          ALLOCATE(KN(MAXKN))
          CALL LCMGET(IPTRK2,'KN',KN)
          CALL LCMSIX(IPTRK2,'BIVCOL',1)
          CALL LCMLEN(IPTRK2,'T',LC,ITYLCM)
          CALL LCMGET(IPTRK2,'E',E)
          CALL LCMSIX(IPTRK2,' ',2)
          CALL VALU2B(LC,MKN,LX,LY,L4,X,YM,XXX,YYY,FUNKNO,MAT,KN,
     +    3*LX,1,E,AXYZ)
        ELSE IF(ICHX.EQ.2) THEN
*         Raviart-Thomas finite element method
          CALL VALU4B(IELEM,NUN,LX,LY,X,YM,XXX,YYY,FUNKNO,MAT,KEY,
     +    3*LX,1,AXYZ)
        ELSE IF(ICHX.EQ.3) THEN
*         Nodal collocation method (MCFD)
          CALL VALU1B(IDIM,LX,LY,L4,X,YM,XXX,YYY,FUNKNO,MAT,IELEM,
     +    3*LX,1,AXYZ)
        ELSE
          CALL XABORT('SPHGAP: INTERPOLATION NOT IMPLEMENTED(3).')
        ENDIF
        IF(IPRINT.GE.100) WRITE(6,*)'SPHGAP: AXYZ =',(AXYZ(I),I=1,3*LX)
        SYM=0.0
        DO J=1,LX
          DG=(MIN(XP(1),XXX(J+1))-MAX(XM(1),XXX(J)))/REAL(NP)
          DO I=1,NP
            SXM=SXM+AXYZ((J-1)*NP+I)*DG
          ENDDO
        ENDDO
        COUGAP=COUGAP+SYM
        LTOT=LTOT+XP(1)-XM(1)
      ENDIF
      IF(NCODE(4).EQ.5)THEN
        SYP=-1.0
      ELSE
        IF(ICHX.EQ.1) THEN
*         Variational collocation method
          CALL LCMLEN(IPTRK2,'KN',MAXKN,ITYLCM)
          MKN=MAXKN/(LX*LY)
          ALLOCATE(KN(MAXKN))
          CALL LCMGET(IPTRK2,'KN',KN)
          CALL LCMSIX(IPTRK2,'BIVCOL',1)
          CALL LCMLEN(IPTRK2,'T',LC,ITYLCM)
          CALL LCMGET(IPTRK2,'E',E)
          CALL LCMSIX(IPTRK2,' ',2)
          CALL VALU2B(LC,MKN,LX,LY,L4,X,YP,XXX,YYY,FUNKNO,MAT,KN,
     +    3*LX,1,E,AXYZ)
        ELSE IF(ICHX.EQ.2) THEN
*         Raviart-Thomas finite element method
          CALL VALU4B(IELEM,NUN,LX,LY,X,YP,XXX,YYY,FUNKNO,MAT,KEY,
     +    3*LX,1,AXYZ)
        ELSE IF(ICHX.EQ.3) THEN
*         Nodal collocation method (MCFD)
          CALL VALU1B(IDIM,LX,LY,L4,X,YP,XXX,YYY,FUNKNO,MAT,IELEM,
     +    3*LX,1,AXYZ)
        ELSE
          CALL XABORT('SPHGAP: INTERPOLATION NOT IMPLEMENTED(4).')
        ENDIF
        IF(IPRINT.GE.100) WRITE(6,*)'SPHGAP: AXYZ =',(AXYZ(I),I=1,3*LX)
        SYP=0.0
        DO J=1,LX
          DG=(MIN(XP(1),XXX(J+1))-MAX(XM(1),XXX(J)))/REAL(NP)
          DO I=1,NP
            SXP=SXP+AXYZ((J-1)*NP+I)*DG
          ENDDO
        ENDDO
        COUGAP=COUGAP+SYP
        LTOT=LTOT+XP(1)-XM(1)
      ENDIF
      IF(IPRINT.GE.100) WRITE(6,*)'SPHGAP: S-XY-PM =',SXM,SXP,SYM,SYP

*  Compute the average flux
      IF(LTOT.EQ.0.0) CALL XABORT('SPHGAP: Error boundary flux = 0.0')
      IF(IPRINT.GE.100) WRITE(6,*)'SPHGAP: COUGAP =',COUGAP,' LTOT=',
     1  LTOT,'Before normalization'
      COUGAP=COUGAP/LTOT
      IF(IPRINT.GE.5) WRITE(6,*)'SPHGAP: COUGAP =',COUGAP,' LTOT=',LTOT
*----
*  DEALLOCATE
*----
      DEALLOCATE(AXYZ,Y,X)
      DEALLOCATE(YYY,XXX,YY,XX)
      RETURN
      END
