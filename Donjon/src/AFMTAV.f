*DECK AFMTAV
      SUBROUTINE AFMTAV (NBURN,ITM,XBMAX,XBMIN,YS,NBMIN,NBMAX,XB,SIGAV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Time average calculation using different approximation.
*
*Copyright:
* Copyright (C) 1996 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* M.T. Sissaoui
*
*Parameters: input
* NBURN  total number of steps.
* ITM    type of the approximation (1-Lagrange;  2-spline; 3-Hermite)
* XBMAX  highest value.
* XBMIN  lower value.
* YS     parameter to be integrated
* NBMIN   
* NBMAX   
* XB     steps
*
*Parameters: output
* SIGAV  average value of YS
*
*-----------------------------------------------------------------------
*
      REAL YS(NBURN),XB(NBURN),SIGAV
      REAL UU(2)
      DOUBLE PRECISION DD
      REAL, ALLOCATABLE, DIMENSION(:) :: Y,U
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(Y(NBURN),U(NBURN))
*
      IF(NBMAX.GT.0) THEN
        INMAX=NBMAX-2
        INMIN=NBMIN
      ELSE
        INMAX=ABS(NBMAX)-1
        INMIN=NBMIN
        IF(ITM.EQ.1)
     1  CALL XABORT('AFMTAV: MORE BURNUP STEPS ARE REQUIRED TO USE '
     1 //' LAGRANGE METHOD, CHOOSE HERMIT OR SPLINE METHOD')
      ENDIF
      IF(ABS(NBMAX).GT.NBMIN) THEN
        SIGAV=0.0
*
        IF(ITM.EQ.1) THEN
* TIME AVERAGE CALCULATION USING LAGRANGE APPROXIMATION.
*
          DO 113 IR=INMIN,INMAX
            I1=IR
            I2=IR+1
            I3=IR+2
            XBI=MAX(XBMIN,XB(IR))
            XBF=MIN(XBMAX,XB(IR+1))
            TX=XBF-XBI
            TX2=XBF**2-XBI**2
            TX3=XBF**3-XBI**3
            X12=XB(I1)-XB(I2)
            X13=XB(I1)-XB(I3)
            X23=XB(I2)-XB(I3)
            XA12=XB(I1)+XB(I2)
            XA13=XB(I1)+XB(I3)
            XA23=XB(I2)+XB(I3)
            XM12=XB(I1)*XB(I2)
            XM13=XB(I1)*XB(I3)
            XM23=XB(I2)*XB(I3)
            Y1=YS(I1)/(X12*X13)
            Y2=-YS(I2)/(X12*X23)
            Y3=YS(I3)/(X13*X23)
*
            SIGAV=SIGAV +
     1               Y1*(TX3/3.0-XA23*TX2/2.0+XM23*TX)+
     1               Y2*(TX3/3.0-XA13*TX2/2.0+XM13*TX)+
     1               Y3*(TX3/3.0-XA12*TX2/2.0+XM12*TX)
 113      CONTINUE
*
        ELSE IF(ITM.EQ.2) THEN
* TIME AVERAGE CALCULATION USING SPLINE APPROXIMATION.
*     THE LOWER BOUNDARY CONDITION IS SET TO BE NATURAL
          Y(1)=0.0
          U(1)=0.0
*     THE UPPER BOUNDARY CONDITION IS SET EITHER TO BE NATURAL
          QN=0.0
          UN=0.0
*
          DO 103 IR=2,NBURN-1
            SIG=(XB(IR)-XB(IR-1))/(XB(IR+1)-XB(IR-1))
            P=SIG*Y(IR-1)+2.0
            Y(IR)=(SIG-1.0)/P
            U(IR)=(6.*((YS(IR+1)-YS(IR))/(XB(IR+1)-XB(IR))-
     1            (YS(IR)-YS(IR-1))/(XB(IR)-XB(IR-1)))/(XB(IR+1)-
     1            XB(IR-1))-SIG*U(IR-1))/P
 103      CONTINUE
*
          Y(NBURN)=(UN-QN*U(NBURN-1))/(QN*Y(NBURN-1)+1.0)
*
          DO 104 K=NBURN-1,1,-1
            Y(K)=Y(K)*Y(K+1)+U(K)
 104      CONTINUE
*
*  COMPUTE THE INTEGRAL OF THE X-SECTION
          INMAX=NBMAX-2
          INMIN=NBMIN
          DO 300 IR=INMIN,INMAX
            H=XB(IR+1)-XB(IR)
            XBI=MAX(XBMIN,XB(IR))
            XBF=MIN(XBMAX,XB(IR+1))
*
            DB=XBF-XBI
            HF=XB(IR+1)-XBF
            HI=XB(IR+1)-XBI
*
            AI=-0.5*(HF**2-HI**2)/H
            BI=DB-AI
            CI=-(AI/6)*H**2-(HF**4-HI**4)/(24*H)
            DI=-(BI/6)*H**2-(HF**4-HI**4)/(24*H)
*
            SIGAV=SIGAV+AI*YS(IR)+BI*YS(IR+1)+
     1               CI*Y(IR)+DI*Y(IR+1)
  300     CONTINUE
        ELSE IF(ITM.EQ.3) THEN
* TIME AVERAGE CALCULATION USING HERMIT APPROXIMATION.
          DO 101 I=1,NBURN
            Y(I)=YS(I)
  101     CONTINUE
* TAKE THE DERIVATIVE WITH RESPECT TO BURNUP OR NEUTRON EXPOSURE AT
* TABULATION POINTS.
          CALL ALDERV(NBURN,XB,Y)
*
*  COMPUTE THE INTEGRAL OF THE X-SECTION
          DD=0.0D0
          DO 200 IR=1,NBURN-1
          IF((XBMIN.LT.XB(IR+1)).AND.(XBMAX.GT.XB(IR))) THEN
            DX=XB(IR+1)-XB(IR)
            XBI=MAX(XBMIN,XB(IR))
            XBF=MIN(XBMAX,XB(IR+1))
            CC=0.5*(XBF-XBI)
            U1=(XBI-0.5*(XB(IR)+XB(IR+1)))/DX
            U2=(XBF-0.5*(XB(IR)+XB(IR+1)))/DX
            UU(1)=0.5*(-(U2-U1)*0.577350269189626+U1+U2)
            UU(2)=0.5*((U2-U1)*0.577350269189626+U1+U2)
            DO 190 J=1,2
            H1=3.0*(0.5-UU(J))**2-2.0*(0.5-UU(J))**3
            H2=(0.5-UU(J))**2-(0.5-UU(J))**3
            H3=3.0*(0.5+UU(J))**2-2.0*(0.5+UU(J))**3
            H4=-(0.5+UU(J))**2+(0.5+UU(J))**3
            DD=DD+(H1*YS(IR)+H2*Y(IR)*DX+H3*YS(IR+1)+
     1      H4*Y(IR+1)*DX)*CC
  190       CONTINUE
          ENDIF
  200     CONTINUE
          SIGAV=REAL(DD)
        ENDIF
        SIGAV=SIGAV/(XBMAX-XBMIN)
      ELSE
        SIGAV=YS(NBMIN)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(U,Y)
      RETURN
      END
