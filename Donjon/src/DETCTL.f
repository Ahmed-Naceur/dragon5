*DECK DETCTL
      SUBROUTINE DETCTL(NX,NY,NZ,NEL,VECT,RESP,NDET,XCT,YCT,ZCT,COR,
     1           KEYF,IPRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Call the subroutines that perform the parabolic interpolation.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin, M. Guyot
*
*Parameters:
* NX     number of x mesh-splitted elements 
* NY     number of y mesh-splitted elements 
* NZ     number of z mesh-splitted elements
* NEL    number of finite elements
* VECT
* RESP   flux reads by the detector
* NDET   number of detectors
* XCT    center coordinates of each mesh-splitted elements for x
* YCT    center coordinates of each mesh-splitted elements for y
* ZCT    center coordinates of each mesh-splitted elements for z
* COR    coordinates of the center of the detector
* KEYF   keyflux recover from L_TRACK object
* IPRT   printing index
*
*-----------------------------------------------------------------------                                                                    *
 
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER    NX,NY,NZ,NEL,NDET,KEYF(NEL),IPRT
      REAL       VECT(*),COR(*),XCT,YCT,ZCT,RESP(NDET)
*----
*  LOCAL VARIABLES
*----
      REAL      D1,D2,D3,X1,X2,X3,XX1,XX2,XX3,Y1,Y2,
     1          Y3,YY1,YY2,YY3,Z1,Z2,Z3,ZZ1,ZZ2,ZZ3,PD1,PD2,PD3,PPD1,
     2          PPD2,PPD3,CE,BE,AH
      INTEGER   I,III,NIJK,I1,I2,I3,IP1,IP2,IP3,J1,J2,J3,JP1,JP2,JP3,
     1          K1,K2,K3,KP1,KP2,KP3,K0,JJJ

      IF(IPRT.GT.4) WRITE(6,1000)

      IF (NDET.LE.0) RETURN
      NIJK = NX*NY

      DO 10 III=1,NDET
         I = (III-1)*3
         D1 = COR(I+1)
         D2 = COR(I+2)
         D3 = COR(I+3)
*----
*  DETERMINE CENTER OF INTERPOLATE RANGE
*----
         CALL DETRTR(D1,XCT,NX,XX1,XX2,XX3,IP1,IP2,IP3)
         X1 = XX1
         X2 = XX2
         X3 = XX3
         I1 = IP1
         I2 = IP2
         I3 = IP3

         CALL DETRTR(D2,YCT,NY,YY1,YY2,YY3,JP1,JP2,JP3)
         Y1 = YY1
         Y2 = YY2
         Y3 = YY3
         J1 = JP1
         J2 = JP2
         J3 = JP3

         CALL DETRTR(D3,ZCT,NZ,ZZ1,ZZ2,ZZ3,KP1,KP2,KP3)
         Z1 = ZZ1
         Z2 = ZZ2
         Z3 = ZZ3
         K1 = KP1
         K2 = KP2
         K3 = KP3

         IF (IPRT.GT.4) THEN
           IF (MOD(III,25).EQ.0) WRITE(6,1000)
         ENDIF

         IF(IPRT.GT.4) THEN
           WRITE(6,2000) III,D1,X1,X2,X3,D2,Y1,Y2,Y3,D3,Z1,Z2,Z3,
     >                          I1,I2,I3,   J1,J2,J3,   K1,K2,K3
         ENDIF
*----
*  INTERPOLATION IN X AT PLANE Z=K1
*----
         K0   = (K1-1)*NIJK
*----
*  INTERPOLATION IN X AT PLANE Y=J1,Z=K1
*----
         JJJ  = NX*(J1-1)
         PD1 = DETPOL(VECT,KEYF,JJJ,K0,I1,I2,I3,X1,X2,X3,D1)
*----
*  INTERPOLATION IN X AT PLANE Y=J2,Z=K1
*----
         JJJ  = NX*(J2-1)
         PD2 = DETPOL(VECT,KEYF,JJJ,K0,I1,I2,I3,X1,X2,X3,D1)
*----
*  INTERPOLATION IN X AT PLANE Y=J3,Z=K1
*----
         JJJ  = NX*(J3-1)
         PD3 = DETPOL(VECT,KEYF,JJJ,K0,I1,I2,I3,X1,X2,X3,D1)
*----
*  INTERPOLATION IN Y AT PLANE Z=K1
*----
         CALL DETPAR(Y1,Y2,Y3,PD1,PD2,PD3,AH,BE,CE)
         PPD1 = AH*D2*D2 + BE*D2 + CE
*----
*  INTERPOLATION IN X AT PLANE Z=K2
*----
         K0   = (K2-1)*NIJK
*----
*  INTERPOLATION IN X AT PLANE Y=J1,Z=K2
*----
         JJJ  = NX*(J1-1)
         PD1 = DETPOL(VECT,KEYF,JJJ,K0,I1,I2,I3,X1,X2,X3,D1)
*----
*  INTERPOLATION IN X AT PLANE Y=J2,Z=K2
*----
         JJJ  = NX*(J2-1)
         PD2 = DETPOL(VECT,KEYF,JJJ,K0,I1,I2,I3,X1,X2,X3,D1)
*----
*  INTERPOLATION IN X AT PLANE Y=J3,Z=K2
*----
         JJJ  = NX*(J3-1)
         PD3 = DETPOL(VECT,KEYF,JJJ,K0,I1,I2,I3,X1,X2,X3,D1)
*----
*  INTERPOLATION IN Y AT PLANE Z=K2
*----
         CALL DETPAR2(Y1,Y2,Y3,PD1,PD2,PD3,AH,BE,CE)
         PPD2 = AH*D2*D2 + BE*D2 + CE
*----
*  INTERPOLATION IN X AT PLANE Z=K3
*----
         K0   = (K3-1)*NIJK
*----
*  INTERPOLATION IN X AT PLANE Y=J1,Z=K3
*----
         JJJ  = NX*(J1-1)
         PD1 = DETPOL(VECT,KEYF,JJJ,K0,I1,I2,I3,X1,X2,X3,D1)
*----
*  INTERPOLATION IN X AT PLANE Y=J2,Z=K3
*----
         JJJ  = NX*(J2-1)
         PD2 = DETPOL(VECT,KEYF,JJJ,K0,I1,I2,I3,X1,X2,X3,D1)
*----
*  INTERPOLATION IN X AT PLANE Y=J3,Z=K3
*----
         JJJ  = NX*(J3-1)
         PD3 = DETPOL(VECT,KEYF,JJJ,K0,I1,I2,I3,X1,X2,X3,D1)
*----
*  INTERPOLATION IN Y AT PLANE Z=K3
*----
         CALL DETPAR2(Y1,Y2,Y3,PD1,PD2,PD3,AH,BE,CE)
         PPD3 = AH*D2*D2 + BE*D2 + CE
*----
*  INTERPOLATION IN Z
*----
         CALL DETPAR2(Z1,Z2,Z3,PPD1,PPD2,PPD3,AH,BE,CE)
         RESP(III) = AH*D3*D3 + BE*D3 + CE

  10  CONTINUE

      RETURN

 1000 FORMAT(//,57X,'BRACKETING PROCESS',
     >        /,57X,'******************',
     >       //,5X,'DET',4X,'X ',8X,'X1',8X,'X2',8X,'X3',4X,
     >                   4X,'Y ',8X,'Y1',8X,'Y2',8X,'Y3',4X,
     >                   4X,'Z ',8X,'Z1',8X,'Z2',8X,'Z3',4X,/)
 2000 FORMAT(5X,I3.3,12F10.4,/,5X,3X,3(10X,3(2X,I6.6,2X)))

      END
