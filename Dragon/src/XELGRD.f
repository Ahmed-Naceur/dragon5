*DECK XELGRD
      SUBROUTINE XELGRD( IPGEOM,   IPRT,   NDIM, NEXTGE,  ITURN, DMESHO,
     >                     MAXC, RMESHO,  MINDO,  MAXDO, ICORDO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read the geometric input for a specific type of cell.
*
*Copyright:
* Copyright (C) 1987 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPGEOM  pointer to the geometry (L_GEOM).             
* IPRT    intermediate printing level for output.      
* NDIM    number of dimensions.                  
* NEXTGE  rectangular(0)/circular(1) boundary.          
* ITURN   turn index for the geometry (from 1 to 16).  
* DMESHO  dimension of array RMESHO.                    
*
*Parameters: output
* MAXC    number of real meshes to stock in RMESHO.       
* RMESHO  real mesh values (rect/cyl).                 
* MINDO   min index values for all axes (rect/cyl).    
* MAXDO   max index values for all axes (rect/cyl).    
* ICORDO  principal axes direction (X/Y/Z) for meshes. 
*
*-----------------------------------------------------------------------
*
      USE                GANLIB
      IMPLICIT           NONE
*
*     DECLARE DUMMY ARGUMENTS
      TYPE(C_PTR)        IPGEOM 
      INTEGER            IPRT, NDIM, NEXTGE, ITURN, DMESHO, MAXC
      REAL               RMESHO(DMESHO)
      INTEGER            MAXDO(*),MINDO(*),ICORDO(*)
*
*     DECLARE LOCAL VARIABLES
      INTEGER            NSTATE, IOUT, MAXTUR
      PARAMETER        ( NSTATE=40, IOUT=6, MAXTUR=12 )
      REAL               RGAR,CENTER(3),DCENT,RMAX
      INTEGER            ISTATE(NSTATE),ITMIX(2*MAXTUR,3),
     >                   ITXYZ(2*MAXTUR,3),IDCEN,IROT
      DOUBLE PRECISION   GAR,DEL
      CHARACTER          TEDATA*12, TEMESH(4)*1, HSMG*131
      INTEGER            ILEN, ILE2, ITYLCM, IBEGIN, INEW, IOLD,
     >                   I, J, K, ISS, ITYPE, ICI, IDIMEN, ICTYPE
      DOUBLE PRECISION   PI,PIO2,FACT
      PARAMETER        ( PI = 3.14159265358979323846D0, PIO2= 0.5D0*PI)
*
*     ALLOCATABLE ARRAYS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISPLT
*
*     DATA STATEMENTS
      SAVE               TEMESH,ITMIX,ITXYZ
      DATA      TEMESH / 'X','Y','Z','C' /
*                        'X'-AXIS AND ITS SIGN
      DATA      ITMIX  /  1,-2,-1, 2,-1, 2, 1,-2, 0, 0, 0, 0,
     >                    1,-2,-1, 2,-1, 2, 1,-2, 0, 0, 0, 0,
*                        'Y'-AXIS AND ITS SIGN
     >                    2, 1,-2,-1, 2, 1,-2,-1, 0, 0, 0, 0,
     >                    2, 1,-2,-1, 2, 1,-2,-1, 0, 0, 0, 0,
*                        'Z'-AXIS AND ITS SIGN
     >                    3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0,
     >                   -3,-3,-3,-3,-3,-3,-3,-3, 0, 0, 0, 0 /
*
      DATA      ITXYZ  /  1, 2,-1,-2,-1, 2, 1,-2, 0, 0, 0, 0,
     >                    1, 2,-1,-2,-1, 2, 1,-2, 0, 0, 0, 0,
*                        'Y'-AXIS AND ITS SIGN
     >                    2,-1,-2, 1, 2, 1,-2,-1, 0, 0, 0, 0,
     >                    2,-1,-2, 1, 2, 1,-2,-1, 0, 0, 0, 0,
*                        'Z'-AXIS AND ITS SIGN
     >                    3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0,
     >                   -3,-3,-3,-3,-3,-3,-3,-3, 0, 0, 0, 0 /
*
      IF( IPRT.GT.1 )THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(/24H           CELL MESHING     )')
      ENDIF
      IF((ITURN.LE.0).OR.(ITURN.GT.24)) THEN
         WRITE(HSMG,'(24H XELGRD: INVALID ITURN (,I6,2H).)') ITURN
         CALL XABORT(HSMG)
      ENDIF
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      ITYPE= ISTATE(1)
      RMAX= 0.0
      IF( ITYPE.EQ.3 .OR. ITYPE.EQ.6 )THEN
         IF( NEXTGE.NE.1 ) CALL XABORT( 'XELGRD: TYPE IS NOT '//
     >                     'COMPATIBLE WITH CIRCULAR B.C.' )
*        GET MAXIMUM RADIUS &
*        INSCRIBE CIRCULAR REGION IN A SQUARE MESH WITH VOLUME=SURFACE
         TEDATA= 'RADIUS'
         CALL LCMLEN(IPGEOM,TEDATA, ILEN, ITYLCM)
         IF(ILEN.GT.0) THEN
            CALL LCMGET(IPGEOM,TEDATA,RMESHO)
            RMAX= RMESHO(ILEN)
         ENDIF
      ENDIF
      IBEGIN = 1
      DO 10 I= 1, 3
         ICORDO(I)=ITMIX(ITURN,I)
         ICI= ABS(ITXYZ(ITURN,I))
         MINDO(I)= IBEGIN
         IF( NEXTGE.EQ.1.AND.I.LE.2 )THEN
            ILEN= 2
            FACT= 0.5D0 * DSQRT(PI+PIO2*PIO2)
            RMESHO(IBEGIN)=  -REAL(FACT)*RMAX
            RMESHO(IBEGIN+1)= REAL(FACT)*RMAX
            IDIMEN= IBEGIN + 1
         ELSE
            IF( I.LE.NDIM )THEN
               TEDATA= 'MESH'//TEMESH(ICI)
               CALL LCMLEN(IPGEOM,TEDATA, ILEN, ITYLCM)
               CALL LCMGET(IPGEOM,TEDATA,RMESHO(IBEGIN))
               IDIMEN= IBEGIN + ILEN - 1
               TEDATA= 'SPLIT'//TEMESH(ICI)
               CALL LCMLEN(IPGEOM,TEDATA, ILE2, ITYLCM)
               IF( ILE2.NE.0 )THEN
                  IF( ILE2.NE.ILEN-1 )THEN
                     CALL XABORT( 'XELGRD: '//TEDATA//' IS INVALID')
                  ELSE
                     ALLOCATE(ISPLT(ILE2))
                     CALL LCMGET(IPGEOM,TEDATA,ISPLT)
                     INEW=0
                     DO 5 IOLD= 1,ILE2
                        INEW= INEW+ ISPLT(IOLD)
    5                CONTINUE
                     K= INEW+1
                     GAR= DBLE(RMESHO(IDIMEN))
                     DO 7 IOLD= ILE2,1,-1
                        DEL= (GAR-DBLE(RMESHO(IBEGIN+IOLD-1)))
     >                      /DBLE(ISPLT(IOLD))
                        GAR= DBLE(RMESHO(IBEGIN+IOLD-1))
                        DO 6 ISS= ISPLT(IOLD),1,-1
                           RMESHO(IBEGIN+K-1)= REAL(GAR+DEL*REAL(ISS))
                           K=K-1
    6                   CONTINUE
    7                CONTINUE
                     IDIMEN= IBEGIN + INEW
                     ILEN=   INEW + 1
                     DEALLOCATE(ISPLT)
                  ENDIF
               ENDIF
               IF( ITXYZ(ITURN,I).LT.0 )THEN
                  DO 3 IOLD= 1,ILEN/2
                     RGAR=                    -RMESHO(IBEGIN+IOLD-1)
                     RMESHO(IBEGIN+IOLD-1)=   -RMESHO(IBEGIN+ILEN-IOLD)
                     RMESHO(IBEGIN+ILEN-IOLD)= RGAR
    3             CONTINUE
                  IF( 2*(ILEN/2).NE.ILEN )THEN
                     RMESHO(IBEGIN+ILEN/2)= -RMESHO(IBEGIN+ILEN/2)
                  ENDIF
               ENDIF
               IF( IPRT.GT.2 )THEN
                  WRITE(IOUT,'(1H )')
            WRITE(IOUT,'(8X,A1,13H-COORDINATES:/(9X,5(1x,F13.6)))')
     >                  TEMESH(I), (RMESHO(J),J=IBEGIN,IDIMEN)
               ENDIF
            ELSE
               ILEN= 2
               RMESHO(IBEGIN)=   0.0
               RMESHO(IBEGIN+1)= 1.0
               IDIMEN= IBEGIN + 1
            ENDIF
         ENDIF
         MAXDO(I)= IDIMEN
         IBEGIN= IDIMEN + 1
   10 CONTINUE
*
*     DETERMINATE COORDINATES OF CENTER
      IF( ITYPE.GE.20.OR.NEXTGE.EQ.1 )THEN
         IF( NEXTGE.EQ.1 )THEN
            ICTYPE= 3
         ELSE
            ICTYPE= ITYPE-20
            IF( ICTYPE.EQ.0 ) ICTYPE= 3
         ENDIF
*
*        GET OFFCENTER VARIATION
         TEDATA= 'OFFCENTER'
         CALL LCMLEN(IPGEOM,TEDATA, ILEN, ITYLCM)
         IF( ILEN .EQ. 0 )THEN
            CENTER(1)=0.0
            CENTER(2)=0.0
            CENTER(3)=0.0
         ELSE
            CALL LCMGET(IPGEOM,TEDATA,CENTER)
         ENDIF
*
*        GET RADIUS, THE FIRST ONE MUST BE 0.0
         TEDATA= 'RADIUS'
         CALL LCMLEN(IPGEOM,TEDATA, ILEN, ITYLCM)
         IF(ILEN.EQ.0) THEN
            ILEN=1
            RMESHO(IBEGIN+1)=0.0
         ELSE
            CALL LCMGET(IPGEOM,TEDATA,RMESHO(IBEGIN+1))
         ENDIF
         IDIMEN= IBEGIN + ILEN
         TEDATA= 'SPLITR'
         CALL LCMLEN(IPGEOM,TEDATA, ILE2, ITYLCM)
         IF( ILE2.NE.0 )THEN
            IF( ILE2.NE.ILEN-1 )THEN
               CALL XABORT( 'XELGRD: '//TEDATA//' IS INVALID')
            ELSE
               ALLOCATE(ISPLT(ILE2))
               CALL LCMGET(IPGEOM,TEDATA,ISPLT)
               INEW=0
               DO 15 IOLD= 1,ILE2
                  INEW= INEW+ ABS(ISPLT(IOLD))
   15          CONTINUE
               K= INEW+1
               GAR= DBLE(RMESHO(IDIMEN))
               DO 17 IOLD= ILE2,1,-1
                  DEL= (GAR-DBLE(RMESHO(IBEGIN+IOLD)))
     >                /DBLE(ABS(ISPLT(IOLD)))
                  IF(ISPLT(IOLD).LT.0)THEN
                     DEL= DEL*(GAR+RMESHO(IBEGIN+IOLD))
                  ENDIF
                  GAR= DBLE(RMESHO(IBEGIN+IOLD))
                  DO 16 ISS= ABS(ISPLT(IOLD)),1,-1
                   IF( ISPLT(IOLD).GT.0 )THEN
                     RMESHO(IBEGIN+K)= REAL(GAR+DEL*REAL(ISS))
                   ELSE
                     RMESHO(IBEGIN+K)= SQRT(REAL(GAR*GAR+DEL*REAL(ISS)))
                   ENDIF
                   K=K-1
   16             CONTINUE
   17          CONTINUE
               IDIMEN= IBEGIN + INEW + 1
               DEALLOCATE(ISPLT)
            ENDIF
         ENDIF
         IF( RMESHO(IBEGIN+1).NE.0.0 )THEN
            WRITE(IOUT,'(11H           ,17HRADII OF ANNULI: /
     >               (11X,5(1X,F13.6)))')
     >               (RMESHO(J),J=IBEGIN+1,IDIMEN)
            CALL XABORT( 'XELGRD: FIRST RADIUS MUST BE 0.0')
         ENDIF
         ICORDO(4)= ICTYPE
         DO 20 I= 0, 1
            IDCEN=MOD(ICTYPE+I,3)+1
            IROT=ITXYZ(ITURN,IDCEN)
            DCENT=CENTER(ABS(IROT))
            IF( IROT .LT. 0 )THEN
              DCENT=-DCENT
            ENDIF
            RMESHO(IBEGIN+I)= 0.5*(RMESHO(MINDO(IDCEN))
     >                            +RMESHO(MAXDO(IDCEN)))+DCENT
   20    CONTINUE
         MINDO(4)= IBEGIN+2
         MAXDO(4)= IDIMEN
         IF( IPRT.GT.2 )THEN
            WRITE(IOUT,'(1H )')
            IF( NEXTGE.EQ.0 )THEN
               WRITE(IOUT,'(9H        (,A1,1H,,A1,10H)- CENTRE: ,
     >                  2H (,2(1X,F13.6),1H) )')
     >                  TEMESH(MOD(ICTYPE  ,3)+1),
     >                  TEMESH(MOD(ICTYPE+1,3)+1),
     >                  RMESHO(IBEGIN),RMESHO(IBEGIN+1)
            ENDIF
            IF( NDIM.EQ.3 )THEN
               WRITE(IOUT,'(14H              ,A1,7H-RADII:/
     >                  (15X,5(1X,F13.6)))')
     >                  TEMESH(ICTYPE),
     >                  (RMESHO(J),J=IBEGIN+2,IDIMEN)
            ELSE
               WRITE(IOUT,'(24H                 RADII:  /
     >                  (15X,5(1X,F13.6)))')
     >               (RMESHO(J),J=IBEGIN+2,IDIMEN)
            ENDIF
         ENDIF
         DO 40 J  = IBEGIN+2, IDIMEN
            RMESHO(J)= RMESHO(J) * RMESHO(J)
   40    CONTINUE
      ENDIF
      MAXC= IDIMEN
*
      RETURN
      END
