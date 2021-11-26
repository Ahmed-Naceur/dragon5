*DECK PIJSMD
      SUBROUTINE PIJSMD(IMPX,NBMIX,NREGIO,MATCOD,VOLUME,XSSIGW,XSSIGT,
     >                  ILK,PIJSYM,PIJSCT,IOP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluate scattering modified cp matrix.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IMPX    print/check flag (equal to 0 for no print).
* NBMIX   number of mixtures considered.
* NREGIO  number of regions considered.
* MATCOD  material code in region.
* VOLUME  volume of region.
* XSSIGW  within group scattering 0 or 1 harmonic.
* XSSIGT  total macroscopic cross sections.
* ILK     leakage flag (ILK=.true. if leakage exists).
* PIJSYM  group condensed reduced/symmetric pij or pijk matrix.
* IOP     pij (=1) or pijk (=4) collision probability flag.
*
*Parameters: output
* PIJSCT  XSSIGW-modified cp matrix (pij or pijk).
*
*-----------------------------------------------------------------------
*
      INTEGER     IMPX,NBMIX,NREGIO,MATCOD(NREGIO),IOP
      REAL        VOLUME(NREGIO),XSSIGW(0:NBMIX),XSSIGT(0:NBMIX),
     >            PIJSYM(NREGIO*(NREGIO+1)/2)
      DOUBLE PRECISION PIJSCT(NREGIO,2*NREGIO)
      LOGICAL     ILK
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,EPS1=1.0E-4)
      DOUBLE PRECISION WRK,F1
*----
*  INTRINSIC FUNCTION FOR POSITION IN CONDENSE PIJ MATRIX
*----
      INDPOS(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
*----
*  PRINT REDUCED PIJ MATRIX BEFORE SCATTERING REDUCTION
*----
      IF(IMPX.GE.8) THEN
        WRITE(IUNOUT,'(/22H MACROSCOPIC TOTAL XS:/(9X,1P,11E11.3))')
     >  (XSSIGT(MATCOD(I)),I=1,NREGIO)
        WRITE(IUNOUT,'(/40H MACROSCOPIC WITHIN-GROUP SCATTERING XS:/
     >  (9X,1P,11E11.3))') (XSSIGW(MATCOD(I)),I=1,NREGIO)
      ENDIF
      IF(IMPX.GE.10) THEN
        IF(IOP.EQ.1) THEN
          WRITE(IUNOUT,200)
        ELSE
          WRITE(IUNOUT,210)
        ENDIF
        WRITE(IUNOUT,240) (J,J=1,NREGIO)
        DO 10 I=1,NREGIO
          WRITE(IUNOUT,250) I,(PIJSYM(INDPOS(I,J))/VOLUME(I),J=1,NREGIO)
   10   CONTINUE
        WRITE(IUNOUT,'(//)')
      ENDIF
*----
*  COMPUTE SCATTERING MODIFIED PIJ
*----
      DO 30 I=1,NREGIO
        DO 20 J=1,NREGIO
          INDPIJ=INDPOS(I,J)
          PIJSCT(I,J)=-XSSIGW(MATCOD(J))*PIJSYM(INDPIJ)
          PIJSCT(I,NREGIO+J)=PIJSYM(INDPIJ)
   20   CONTINUE
        PIJSCT(I,I)=VOLUME(I)+PIJSCT(I,I)
   30 CONTINUE
      CALL ALSBD(NREGIO,NREGIO,PIJSCT,IERROR,NREGIO)
      IF(IERROR.NE.0) CALL XABORT('PIJSMD: SINGULAR MATRIX.')
      DO 50 I=1,NREGIO
        DO 40 J=1,NREGIO
          PIJSCT(I,J)=PIJSCT(I,NREGIO+J)
   40   CONTINUE
   50 CONTINUE
      IF(IMPX.GE.8) THEN
        IF(IOP.EQ.1) THEN
          WRITE(IUNOUT,220)
        ELSE
          WRITE(IUNOUT,230)
        ENDIF
        WRITE(IUNOUT,240) (J,J=1,NREGIO)
        DO 60 I=1,NREGIO
          WRITE(IUNOUT,250) I,(PIJSCT(I,J),J=1,NREGIO)
   60   CONTINUE
        WRITE(IUNOUT,'(//)')
      ENDIF
      IF((IMPX.GE.10).OR.(IMPX.LT.0).AND.(IOP.EQ.1)) THEN
*----
*  CHECK THE RECIPROCITY CONDITIONS
*----
        VOLTOT=0.0
        DO 70 I=1,NREGIO
          VOLTOT=VOLTOT+VOLUME(I)
   70   CONTINUE
        VOLTOT=VOLTOT/REAL(NREGIO)
        WRK=0.0D0
        DO 90 I=1,NREGIO
          DO 80 J=1,NREGIO
            WRK=MAX(WRK,ABS(PIJSCT(I,J)*VOLUME(I)
     >         -PIJSCT(J,I)*VOLUME(J))/VOLTOT)
   80     CONTINUE
   90   CONTINUE
        IF(WRK.GE.EPS1) WRITE(IUNOUT,260) WRK
*----
*  CHECK THE CONSERVATION CONDITIONS
*----
        IF(.NOT.ILK) THEN
          WRK=0.0D0
          DO 110 I=1,NREGIO
            F1=1.0D0
            DO 100 J=1,NREGIO
              IBM=MATCOD(J)
              F1=F1-PIJSCT(I,J)*(XSSIGT(IBM)-XSSIGW(IBM))
  100       CONTINUE
            WRK=MAX(WRK,ABS(F1))
  110     CONTINUE
          IF(WRK.GE.EPS1) WRITE(IUNOUT,270) WRK
        ENDIF
      ENDIF
      RETURN
*
  200 FORMAT (//51H PIJSMD: REDUCED COLLISION PROBABILITY MATRIX (I --,
     1 6H> J) :/)
  210 FORMAT (//51H PIJSMD: REDUCED DIRECTIONAL COLLISION PROBABILITY ,
     1 18HMATRIX (I --> J) :/)
  220 FORMAT (//51H PIJSMD: SCATTERING-REDUCED COLLISION PROBABILITY M,
     1 17HATRIX (I --> J) :/)
  230 FORMAT (//51H PIJSMD: SCATTERING-REDUCED DIRECTIONAL COLLISION P,
     1 29HROBABILITY MATRIX (I --> J) :/)
  240 FORMAT (11X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,
     1 2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,
     2 5X,2HJ=,I4,:,5X,2HJ=,I4)
  250 FORMAT (3H I=,I4,2H: ,1P,11E11.3/(9X,11E11.3))
  260 FORMAT (/50H PIJSMD: THE SCATTERING MODIFIED CP MATRIX DO NOT ,
     1 40HMEET THE RECIPROCITY CONDITIONS. RECIP =,1P,E10.3/)
  270 FORMAT (/50H PIJSMD: THE SCATTERING MODIFIED CP MATRIX DO NOT ,
     1 40HMEET THE CONSERVATION CONDITIONS. LEAK =,1P,E10.3/)
      END
