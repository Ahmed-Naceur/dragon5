*DECK DETSPL2
      SUBROUTINE DETSPL2(XCNTR ,YCNTR ,NXMAX ,NYMAX ,FXY,
     >                  FP1   ,FP2   ,F2X   ,F2Y   ,FDUMMY,
     >                  XINT  ,YINT  ,FYINT ,FXYINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform spline interpolation.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin
*
*Parameters: 
* XCNTR    
* YCNTR    
* NXMAX    
* NYMAX    
* FXY      
* FP1      
* FP2      
* F2X      
* F2Y      
* FDUMMY   
* XINT     
* YINT     
* FYINT    
* FXYINT   
*
*-----------------------------------------------------------------------
*
      REAL*4     XCNTR(NXMAX),YCNTR(NYMAX),FXY(NXMAX,NYMAX),
     >           XINT        ,YINT        ,FXYINT,
     >           F2X(NXMAX)  ,F2Y(NYMAX)  ,FYINT(NYMAX),
     >           FDUMMY(NXMAX)
*----
*  CALCULATE THE SECOND DERIVATIVES ALONG XCNTR FOR EACH Y
*----
      DO 10 I=1,NYMAX

         DO 20 J=1,NXMAX
            FDUMMY(J) = FXY(J,I)
  20     CONTINUE

         CALL DETSPLI(XCNTR,FDUMMY,NXMAX,FP1,FP2,F2X)
*----
*  INTERPOLATE ALONG THE X COORDINATE FOR EACH Y
*----
         CALL DETSPLI2(XCNTR,FDUMMY,F2X,NXMAX,XINT,FYINT(I))

  10  CONTINUE
*----
*  CALCULATE SECOND DERIVATIVE ALONG Y FOR XINT
*----
      CALL DETSPLI(YCNTR,FYINT,NYMAX,FP1,FP2,F2Y)
*----
*  INTERPOLATE ALONG Y FOR XINT
*----
      CALL DETSPLI2(YCNTR,FYINT,F2Y,NYMAX,YINT,FXYINT)

      RETURN
      END
