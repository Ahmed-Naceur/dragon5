*DECK LIBECT
      SUBROUTINE LIBECT(MAXTRA,LLL,PRI,UUU,DEL,DELTA,NEXT,III,MML,STIS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build the elastic scattering law for neutrons with secondary energy
* in group LLL.
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
* MAXTRA  dimension of array PRI.
* LLL     secondary energy group index.
* PRI     info to rebuild the scat matrix.
* UUU     lethargy limits of the fine groups.
* DEL     elementary lethargy width.
* DELTA   lethargy width of each energy group.
* NEXT    length of x-s structure for the current isotope.
* III     offset in PRI array for the current isotope.
*
*Parameters: output
* MML     number of down-scattering groups (including group LLL).
* STIS    values of the transfer macroscopic cross section:
*         STIS(1):     from group LLL;
*         STIS(2):     from group LLL-1;
*         STIS(LLL):   from group 1.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LLL,NEXT,III,MML
      REAL PRI(MAXTRA),UUU(LLL),DEL,DELTA(LLL),STIS(LLL)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION DAUX
*
      DO 10 I=1,LLL
      STIS(I)=0.0
   10 CONTINUE
      MML=1
      MM=0
      LDELH=INT(UUU(LLL)/DEL+0.1)
      LARGRL=INT(DELTA(LLL)/DEL+0.1)
      LDELB=LDELH-LARGRL+1
      IHM=III+NEXT-1
      LTES=LDELB-NEXT
      INDICE=1
      INTER2=0
      J=0
      DO 70 MM1=1,LLL
      MM=LLL-MM1+1
      MDELH=INT(UUU(MM)/DEL+0.1)
      IF(MDELH.LE.LTES) THEN
         MM=MM+1
         GO TO 80
      ENDIF
      LARGRM=INT(DELTA(MM)/DEL+0.1)
      MDELB=MAX(MDELH-LARGRM+1,LTES+1)
      DAUX=0.0D0
      LARG=MIN(LARGRM,LARGRL)
      IF(LARG.LE.4) THEN
         DO 25 MDEL=MDELB,MDELH
         IBAS=MAX(LDELB-MDEL+III,III)
         IHAUT=MIN(LDELH-MDEL+III,IHM)
         DO 20 I=IBAS,IHAUT
         DAUX=DAUX+PRI(I)
   20    CONTINUE
   25    CONTINUE
         GO TO 60
      ENDIF
      IHAUT=MIN(LDELH-MDELB+III,IHM)
      IF(INDICE.EQ.1) THEN
         INDICE=2
         INTER2=III-1
         J=LARG+1
      ELSE IF(INDICE.EQ.2) THEN
         J=0
         IBAS=MAX(LDELB-MDELH+III,III)
         LARGLI=ABS(LARGRM-LARGRL)
         INTER1=MIN(IBAS+LARG-2,IHAUT)
         DO 30 I=IBAS,INTER1
         J=J+1
         DAUX=DAUX+PRI(I)*REAL(J)
   30    CONTINUE
         INTER1=INTER1+1
         INTER2=MIN(IHAUT,INTER1+LARGLI)
         IF(INTER1.GT.INTER2) GO TO 60
         J=LARG
         DO 40 I=INTER1,INTER2
         DAUX=DAUX+PRI(I)*REAL(LARG)
   40    CONTINUE
      ENDIF
      INTER2=INTER2+1
      IF(INTER2.GT.IHAUT) GO TO 60
      DO 50 I=INTER2,IHAUT
      J=J-1
      DAUX=DAUX+PRI(I)*REAL(J)
   50 CONTINUE
*
   60 STIS(MM1)=REAL(DAUX)
      STIS(MM1)=STIS(MM1)*DEL/DELTA(LLL)
   70 CONTINUE
   80 MML=LLL-MM+1
      RETURN
      END
