*DECK DETSPL3
      SUBROUTINE DETSPL3(XCNTR ,YCNTR ,ZCNTR ,
     >                  NXMAX ,NYMAX ,NZMAX ,
     >                  FXYZ  ,FXY   ,FDUMMY,
     >                  F2X   ,F2Y   ,F2Z   ,
     >                  XINT  ,YINT  ,ZINT  ,
     >                  FP1   ,FP2   ,
     >                  FYINT ,FZINT ,FINTRP,
     >                  N1    ,N2    ,N3    ,ITYPE)
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
*
*Parameters: 
* XCNTR     
* YCNTR     
* ZCNTR     
* NXMAX     
* NYMAX     
* NZMAX     
* FXYZ      
* FXY       
* FDUMMY    
* F2X       
* F2Y       
* F2Z       
* XINT      
* YINT      
* ZINT      
* FP1       
* FP2       
* FYINT     
* FZINT     
* FINTRP    
* N1        
* N2        
* N3        
* ITYPE     
*
*-----------------------------------------------------------------------
*
      REAL*4     XCNTR(NXMAX)  ,YCNTR(NYMAX)     ,ZCNTR(NZMAX),
     >           FXYZ(N1,N2,N3),FXY(NXMAX,NYMAX) ,
     >           FYINT(NYMAX)  ,FZINT(NZMAX)     ,
     >           F2X(NXMAX)    ,F2Y(NYMAX)       ,F2Z(NZMAX),
     >           FDUMMY(NXMAX)
*----
*  INTERPOLATE IN TWO DIMENSIONS AT XINT,YINT FOR EACH Z PLANE
*----
      DO 10 K=1,NZMAX

         DO 20 J=1,NXMAX
            DO 30 I=1,NYMAX

               IF (ITYPE.EQ.1) THEN
                  FXY(J,I) = FXYZ(J,I,K)
               ELSE IF (ITYPE.EQ.2) THEN
                  FXY(J,I) = FXYZ(K,J,I)
               ELSE IF (ITYPE.EQ.3) THEN
                  FXY(J,I) = FXYZ(I,K,J)
               ELSE
                  CALL XABORT('DETSPL3: ERROR IN SPLIN3')
               ENDIF

   30       CONTINUE
  20     CONTINUE

         CALL DETSPL2(XCNTR,YCNTR,NXMAX ,NYMAX ,FXY,
     >               FP1  ,FP2  ,F2X   ,F2Y   ,FDUMMY,
     >               XINT ,YINT ,FYINT ,FXYINT)

         FZINT(K) = FXYINT

   10 CONTINUE
*----
*  CALCULATE SECOND DERIVATIVE ALONG Z AT XINT,YINT
*----
      CALL DETSPLI(ZCNTR,FZINT,NZMAX,FP1,FP2,F2Z)
*----
*  INTERPOLATE ALONG Z FOR XINT,YINT
*----
      CALL DETSPLI2(ZCNTR,FZINT,F2Z,NZMAX,ZINT,FINTRP)

      RETURN
      END
