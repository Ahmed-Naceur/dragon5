*DECK DETSPL
      SUBROUTINE DETSPL(NXMAX,NYMAX,NZMAX,IM,FLUX,FLUXIN,NINT,XCNTR,
     > YCNTR,ZCNTR,COORD,IXX,IPRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the spline interpolation.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin
*
*Parameters: 
* NXMAX     
* NYMAX     
* NZMAX     
* IM        
* FLUX      
* FLUXIN    
* NINT      
* XCNTR     
* YCNTR     
* ZCNTR     
* COORD     
* IXX       
* IPRT      
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NXMAX,NYMAX,NZMAX,IM,NINT,IXX(*),IPRT
      REAL FLUX(*),FLUXIN(NINT),XCNTR(NXMAX),YCNTR(NYMAX),ZCNTR(NZMAX),
     > COORD(*)
*----
*  LOCAL VARIABLES
*----
      LOGICAL     L1DSET
      REAL, ALLOCATABLE, DIMENSION(:) :: FDUMMY,FXINT,FYINT,FZINT,F2X,
     > F2Y,F2Z
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FXY,FYZ,FZX
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: FXYZ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(FXYZ(NXMAX,NYMAX,NZMAX),FDUMMY(IM),FXINT(NXMAX),
     > FYINT(NYMAX),FZINT(NZMAX),FXY(NXMAX,NYMAX),FYZ(NYMAX,NZMAX),
     > FZX(NZMAX,NXMAX),F2X(NXMAX),F2Y(NYMAX),F2Z(NZMAX))
*----
*  SECOND DERIVATIVE IS CALCULATED BASED ON X(I), Y(I) (DEFAULT)
*----
      FP1 = 0.0
      FP2 = 0.0
*----
*  ASSEMBLE THE ARRAY FXYZ OVER THE FULL MESH
*----
      NXNY = NXMAX*NYMAX
*
      DO 10 J=1,NXMAX
         IX = J
         DO 20 I=1,NYMAX
            IY = NXMAX*(I - 1)
            DO 30 K=1,NZMAX
               IZ = NXNY*(K - 1)
*
               IDX = IX + IY + IZ
               IF (IXX(IDX).EQ.0) THEN
                  FXYZ(J,I,K) = 0.0
               ELSE
                  FXYZ(J,I,K) = FLUX(IXX(IDX))
               ENDIF
*
   30       CONTINUE
 20      CONTINUE
10    CONTINUE
*----
*  CALCULATE THE COORDINATES TO INTERPOLATE
*----
      IF(IPRT.GT.4) WRITE(6,1000)
      IF(IPRT.GT.4) WRITE(6,2000)

      N1 = NXMAX
      N2 = NYMAX
      N3 = NZMAX

      DO 40 N=1,NINT
         ININT = 3*(N-1)

         XINT = COORD(ININT + 1)
         YINT = COORD(ININT + 2)
         ZINT = COORD(ININT + 3)
*----
*  INTERPOLATE IN TWO DIMENSIONS AT XINT,YINT FOR EACH Z PLANE
*----
         ITYPE = 1
         CALL DETSPL3(XCNTR ,YCNTR ,ZCNTR ,
     >               NXMAX ,NYMAX ,NZMAX ,
     >               FXYZ  ,FXY   ,FDUMMY,
     >               F2X   ,F2Y   ,F2Z   ,
     >               XINT  ,YINT  ,ZINT  ,
     >               FP1   ,FP2   ,
     >               FYINT ,FZINT ,FINTR1,
     >               N1    ,N2    ,N3    ,ITYPE)

         L1DSET = .TRUE.
         IF (L1DSET) GOTO 1
*----
*  INTERPOLATE IN TWO DIMENSIONS AT YINT,ZINT FOR EACH X PLANE
*----
         ITYPE = 2
         CALL DETSPL3(YCNTR ,ZCNTR ,XCNTR ,
     >               NYMAX ,NZMAX ,NXMAX ,
     >               FXYZ  ,FYZ   ,FDUMMY,
     >               F2Y   ,F2Z   ,F2X   ,
     >               YINT  ,ZINT  ,XINT  ,
     >               FP1   ,FP2   ,
     >               FZINT ,FXINT ,FINTR2,
     >               N1    ,N2    ,N3    ,ITYPE)
*
         IF(IPRT.GT.4) WRITE(6,3000) XINT,YINT,ZINT,FINTR2
*----
*  INTERPOLATE IN TWO DIMENSIONS AT ZINT,XINT FOR EACH Y PLANE
*----
         ITYPE = 3
         CALL DETSPL3(ZCNTR ,XCNTR ,YCNTR ,
     >               NZMAX ,NXMAX ,NYMAX ,
     >               FXYZ  ,FZX   ,FDUMMY,
     >               F2Z   ,F2X   ,F2Y   ,
     >               ZINT  ,XINT  ,YINT ,
     >               FP1   ,FP2   ,
     >               FXINT ,FYINT ,FINTR3,
     >               N1    ,N2    ,N3    ,ITYPE)

         IF(IPRT.GT.4) WRITE(6,3000) XINT,YINT,ZINT,FINTR3
*----
*  GET AVERAGE VALUE
*----
   1     FI = FINTR1

         IF(IPRT.GT.4) WRITE(6,4000) N,XINT,YINT,ZINT,FI

         FLUXIN(N) = FI

  40  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(F2Z,F2Y,F2X,FZX,FYZ,FXY,FZINT,FYINT,FXINT,FDUMMY,FXYZ)
      RETURN
*
 1000 FORMAT(1H1,//,5X,'*** INTERPOLATION PROCESS',/)
 2000 FORMAT(//,1X,'DET NO' ,5X,4X,'XP',4X, 4X,'YP',4X, 4X,'ZP',4X,
     >                       7X,'FI',6X,//)
 3000 FORMAT(   1X,6X       ,5X,F8.3   ,2X, F8.3   ,2X, F8.3,2X,
     >                       3X,1PE12.5)
 4000 FORMAT(   4X,I3.3     ,5X,F8.3   ,2X, F8.3   ,2X, F8.3,2X,
     >                       3X,1PE12.5)

      END
