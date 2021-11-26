*DECK PIJWPR
      SUBROUTINE PIJWPR(LOPT,NREG,NSOUT,SIGTAL,PROB,SIGVOL,MSYM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print-out for probability matrices.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy, G. Marleau
*
*Parameters: input
* LOPT    print-out form:
*         <= 0 print all PSS/PVS/PVV;
*         >  0 print only PVV.
* NREG    total number of regions.
* NSOUT   number of outer surface.
* SIGTAL  albedo-cross section vector.
* PROB    collision probabilities.
* SIGVOL  region volumes.
* MSYM    matrix format:
*         = 1  symmetric matrix;
*         = 0  non-symmetric full matrix.
*
*Comments:
*  Format of compress probability matrix
*    NPLEN=(NREG+NSOUT+2)*(NREG+NSOUT+1)/2
*    IND(I,J)=MAX(I+NSOUT+1,J+NSOUT+1)
*            *(MAX(I+NSOUT+1,J+NSOUT+1)-1)/2
*            +MIN(I+NSOUT+1,J+NSOUT+1)
*    IS=-NSOUT,-1; JS=-NSOUT,IS; I=IND(IS,JS)
*      PROB(I)=VOLSUR(IS)*PSS(IS,JS)
*    IV=1,NREG; JS=-NSOUT,-1;    I=IND(IV,JS)
*      SIGT(IV).GT.0.0
*        PROB(I)=SIGT(IV)*VOLSUR(IV)*PVS(IV,JS)
*      SIGT(IV).EQ.0.0
*        PROB(I)=VOLSUR(IV)*PVS(IV,JS)
*     IV=1,NREG; JV=1,IV;       I=IND(IV,JV)
*       SIGT(IV).GT.0.0 AND SIGT(JV).GT.0.0
*         PROB(I)=SIGT(IV)*SIGT(JV)*VOLSUR(IV)*PVV(IV,JV)
*       SIGT(IV).GT.0.0 AND SIGT(JV).EQ.0.0
*         PROB(I)=SIGT(IV)*VOLSUR(IV)*PVV(IV,JV)
*       SIGT(IV).EQ.0.0 AND SIGT(JV).GT.0.0
*         PROB(I)=SIGT(JV)*VOLSUR(IV)*PVV(IV,JV)
*       SIGT(IV).EQ.0.0 AND SIGT(JV).EQ.0.0
*         PROB(I)=VOLSUR(IV)*PVV(IV,JV)
*
*-----------------------------------------------------------------------
*
      IMPLICIT            NONE
      INTEGER           IUNOUT,   LOPT,   NREG,  NSOUT,   MSYM,
     >                     IND,      I,      J,   NSUR,   NVOL,
     >                   NSURC,  NSURM,  NVOLC,  NVOLM,     IP,     IR,
     >                      JR,    III
      PARAMETER        (IUNOUT=6)
      REAL              SIGTAL(-NSOUT:NREG), BILANP(10), XSJR, WPR,
     >                  VPR(10), SIGVOL(NREG), COF
      DOUBLE PRECISION  PROB(*)
*
      IND(I,J) = MAX(I+NSOUT+1,J+NSOUT+1)*(MAX(I+NSOUT+1,J+NSOUT+1)-1)/2
     >         + MIN(I+NSOUT+1,J+NSOUT+1)
*
      WPR(I,J)= REAL(PROB( IND(I,J) ) / PROB( IND(I,0) ))
*
*NOTE:
*     IF( SIGT(I).NE.0.0 )THEN
*        PROB(IND(I,0)= SIGT(I) * VOLSUR(I)
*     ELSE
*        PROB(IND(I,0)= VOLSUR(I)
*     ENDIF
*
      NSUR= -NSOUT
      NVOL=  NREG
*
      WRITE(IUNOUT,'(24H REGIONAL CROSS SECTIONS)')
      WRITE(IUNOUT,'(5(1X,6HREGION,5X,16HCROSS SECTION   ))')
      WRITE(IUNOUT,'(5(1X,I6,3X,E15.7))') (JR,SIGTAL(JR),JR=1,NVOL)
      IF(MSYM .EQ. 0) THEN
        WRITE(IUNOUT,'(5(1X,6HREGION,5X,16HVOLUMES         ))')
        WRITE(IUNOUT,'(5(1X,I6,3X,E15.7))') 
     >  (JR,SIGVOL(JR),JR=1,NVOL)
      ELSE
        WRITE(IUNOUT,'(5(1X,6HREGION,5X,16HSURFACE/VOLUMES ))')
        WRITE(IUNOUT,'(5(1X,I6,3X,E15.7))') 
     >  (JR,PROB(IND(JR,0)),JR=-NSOUT,NVOL)
      ENDIF
      IF( LOPT.LE.0 )THEN
         NSURC = -1
         DO 40 IP  = 1, (9 - NSUR) / 10
            NSURM= MAX( NSUR, NSURC-9 )
            WRITE(IUNOUT,'(30H0  SURFACE CONSERVATION LAWS: ,
     >                     31H( P.S<-S + P.V<-S = 1 + ERR.S ) ,
     >                     31H FOR XS.TOTAL=0, REDUCED P.V<-S ,
     >                     11H IS PRINTED )')
            WRITE(IUNOUT,'(1X,8H(P.S<-S),1X,10( A5,    I6,:)/)')
     >               (' SUR ',-IR,IR= NSURC, NSURM, -1)
            DO 10 IR  =NSURC, NSURM, -1
               BILANP(IR-NSURM+1)= 0.0
   10       CONTINUE
            DO 25 JR = -1,  NSUR, -1
               WRITE(IUNOUT,'(5H SUR ,I4,1X,10F11.8)')
     >                          -JR, (WPR(IR,JR),IR=NSURC,NSURM,-1)
               DO 20 IR  = NSURC, NSURM, -1
                  BILANP(IR-NSURM+1)=  BILANP(IR-NSURM+1)
     >                              +  WPR(IR,JR)
   20          CONTINUE
   25       CONTINUE
            WRITE(IUNOUT,'(1X,8H(P.V<-S) )')
            DO 35 JR  =  1,  NVOL,  1
               IF( SIGTAL(JR).EQ.0.0 ) THEN
                  WRITE(IUNOUT,'(5H VOL ,I4,1X,3H 0*,10(F8.5,:,3H 0*))')
     >                      JR,(WPR(IR,JR),IR=NSURC,NSURM,-1)
                  XSJR= 0.0
               ELSE
                  WRITE(IUNOUT,'(5H VOL ,I4,1X,10F11.8)')
     >                      JR,(WPR(IR,JR),IR=NSURC,NSURM,-1)
                  XSJR= 1.0
               ENDIF
               DO 30 IR  = NSURC, NSURM, -1
                  BILANP(IR-NSURM+1)= BILANP(IR-NSURM+1)
     >                              + XSJR * WPR(IR,JR)
   30          CONTINUE
   35       CONTINUE
            WRITE(IUNOUT,'(1H )')
            WRITE(IUNOUT,'(5H SUM ,5X,10F11.8)')
     >                          (BILANP(IR-NSURM+1),IR=NSURC,NSURM,-1)
            NSURC = NSURC - 10
   40    CONTINUE
      ENDIF
      NVOLC =  1
      DO 90 IP  = 1, (9 + NVOL) / 10
         NVOLM= MIN( NVOL, NVOLC+9 )
         IF( LOPT.LE.0 )THEN
            WRITE(IUNOUT,'(30H0  VOLUME  CONSERVATION LAWS: ,
     >                     31H( P.S<-V + P.V<-V = 1 + ERR.V ) ,
     >                     31H FOR XS.TOTAL=0, REDUCED P.V<-V ,
     >                     11H IS PRINTED )')
         ELSE
            WRITE(IUNOUT,'(30H0  VOLUME  CONSERVATION LAWS: ,
     >                     32H( SUM OF P.V<-V = 1 + ESCAPE.V ) ,
     >                     31H FOR XS.TOTAL=0, REDUCED P.V<-V ,
     >                     11H IS PRINTED )')
         ENDIF
         DO 50 IR  = NVOLC, NVOLM,   1
            BILANP(IR-NVOLC+1)= 0.0
   50    CONTINUE
         IF( LOPT.LE.0 )THEN
            WRITE(IUNOUT,'(1X,8H(P.S<-V),1X,10( A5 ,  I6,:)/)')
     >                    (' VOL ',IR,IR=NVOLC,NVOLM, 1)
            DO 65 JR = -1,  NSUR, -1
               WRITE(IUNOUT,'(5H SUR ,I4,1X,10F11.8)')
     >         -JR, (WPR(IR,JR),IR=NVOLC,NVOLM, 1)
            DO 60 IR  = NVOLC, NVOLM,   1
               BILANP(IR-NVOLC+1)= BILANP(IR-NVOLC+1)
     >                           + WPR(IR,JR)
   60       CONTINUE
   65       CONTINUE
            WRITE(IUNOUT,'(1X,8H(P.V<-V) )')
         ELSE
            WRITE(IUNOUT,'(1X,8H(P.V<-V),1X,10( A5 ,  I6,:)/)')
     >                    (' VOL ',IR,IR=NVOLC,NVOLM, 1)
         ENDIF
*
         IF(LOPT.GT.0.AND.MSYM.EQ.0)THEN
*
*        PRINTING OF PIJK" FULL MATRIX
*
         COF=1.5
         DO 70 JR =  1,  NVOL,  1
            IF( SIGTAL(JR).EQ.0.0 )THEN
              DO 75 IR=NVOLC, NVOLM, 1
               III=JR+NREG*(IR-1)
               VPR(IR-NVOLC+1)=COF*REAL(PROB(III))/SIGVOL(IR)
   75         CONTINUE
               WRITE(IUNOUT,'(5H VOL ,I4,1X,3H 0*,10(F8.5,:,3H 0*))')
     >         JR,(VPR(IR-NVOLC+1),IR=NVOLC,NVOLM,1)
            ELSE
              DO 76 IR=NVOLC, NVOLM, 1
               III=JR+NREG*(IR-1)
               VPR(IR-NVOLC+1)=COF*REAL(PROB(III))/SIGVOL(IR)
               BILANP(IR-NVOLC+1)=BILANP(IR-NVOLC+1)+VPR(IR-NVOLC+1)
   76         CONTINUE
               WRITE(IUNOUT,'(5H VOL ,I4,1X,10F11.8)')
     >         JR,(VPR(IR-NVOLC+1),IR=NVOLC,NVOLM,1)
            ENDIF
   70    CONTINUE
         WRITE(IUNOUT,'(1H )')
         WRITE(IUNOUT,'(5H SUM ,5X,10F11.8)')
     >                         (BILANP(IR-NVOLC+1),IR=NVOLC,NVOLM, 1)
*
         ELSE
*
         DO 85 JR =  1,  NVOL,  1
            IF( SIGTAL(JR).EQ.0.0 )THEN
               WRITE(IUNOUT,'(5H VOL ,I4,1X,3H 0*,10(F8.5,:,3H 0*))')
     >         JR,(WPR(IR,JR),IR=NVOLC,NVOLM,1)
               XSJR= 0.0
            ELSE
               WRITE(IUNOUT,'(5H VOL ,I4,1X,10F11.8)')
     >         JR,(WPR(IR,JR),IR=NVOLC,NVOLM,1)
               XSJR= 1.0
            ENDIF
         DO 80 IR  = NVOLC, NVOLM,   1
            BILANP(IR-NVOLC+1)= BILANP(IR-NVOLC+1)
     >                        + XSJR*WPR(IR,JR)
   80    CONTINUE
   85    CONTINUE
         WRITE(IUNOUT,'(1H )')
         WRITE(IUNOUT,'(5H SUM ,5X,10F11.8)')
     >                         (BILANP(IR-NVOLC+1),IR=NVOLC,NVOLM, 1)
         ENDIF
         NVOLC = NVOLC + 10
   90 CONTINUE
*
      RETURN
      END
