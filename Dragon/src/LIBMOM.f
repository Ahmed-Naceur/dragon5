*DECK LIBMOM
      SUBROUTINE LIBMOM(NFS,NDIL,NPAR,DELTA,SIGTF,SIGSF,SIGFF,NOR,
     1 SIGERD,MOMT,MOMP,SEFFER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute a set of total (SIGT**N) and partial (SIGA*(SIGT**N)) moments
* and a set of reference self-shielded flux and cross sections.
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
* NFS     number of fine energy groups.
* NDIL    number of dilutions.
* NPAR    number of partial cross sections (NPAR=0, 1 or 2).
* DELTA   lethargy widths of the fine groups.
* SIGTF   microscopic principal x-sections in the fine groups.
* SIGSF   microscopic partial x-sections 1 in the fine groups.
* SIGFF   microscopic partial x-sections 2 in the fine groups.
* NOR     related to the number of moments to preserve.
*         For the total moments: -NOR+1 <= N <= NOR.
*         For the partial moments: -NOR/2 <= N <= (NOR-1)/2.
* SIGERD  dilutions used to compute SEFFER.
*
*Parameters: output
* MOMT    total moments.
* MOMP    partial moments in absorption.
* SEFFER  Bondarenko self-shielded flux and cross sections at each
*         dilution.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NFS,NDIL,NPAR,NOR
      REAL DELTA(NFS),SIGTF(NFS),SIGSF(NFS),SIGFF(NFS),SIGERD(NDIL),
     1 SEFFER(NPAR+2,NDIL)
      DOUBLE PRECISION MOMT(2*NOR),MOMP(NOR,NPAR)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION T,DEL
*
      DEL=0.0D0
      DO 10 I=1,2*NOR
      MOMT(I)=0.0D0
   10 CONTINUE
      DO 21 I=1,NOR
      DO 20 J=1,NPAR
      MOMP(I,J)=0.0D0
   20 CONTINUE
   21 CONTINUE
      DO 26 IDIL=1,NDIL
      DO 25 IPAR=1,NPAR+2
      SEFFER(IPAR,IDIL)=0.0
   25 CONTINUE
   26 CONTINUE
*
      DO 70 IGF=1,NFS
        DELF=DELTA(IGF)
        SIGT=MAX(0.001,SIGTF(IGF))
        DEL=DEL+DELF
        T=DELF
        DO 30 INOR=NOR,2*NOR
        MOMT(INOR)=MOMT(INOR)+T
        T=T*SIGT
   30   CONTINUE
        T=DELF/SIGT
        DO 40 INOR=NOR-1,1,-1
        MOMT(INOR)=MOMT(INOR)+T
        T=T/SIGT
   40   CONTINUE
        DO 45 IDIL=1,NDIL
        T=SIGERD(IDIL)*DELF
        SEFFER(1,IDIL)=SEFFER(1,IDIL)+REAL(T)/(SIGERD(IDIL)+SIGT)
        SEFFER(2,IDIL)=SEFFER(2,IDIL)+REAL(T)*SIGT/(SIGERD(IDIL)+SIGT)
   45   CONTINUE
        IF(NPAR.GT.0) THEN
          SIGS=MAX(1.E-9,SIGSF(IGF))
          T=DELF*SIGS
          DO 50 INOR=NOR/2+1,NOR
          MOMP(INOR,1)=MOMP(INOR,1)+T
          T=T*SIGT
   50     CONTINUE
          T=DELF*SIGS/SIGT
          DO 60 INOR=NOR/2,1,-1
          MOMP(INOR,1)=MOMP(INOR,1)+T
          T=T/SIGT
   60     CONTINUE
          DO 65 IDIL=1,NDIL
          T=SIGERD(IDIL)*DELF
          SEFFER(3,IDIL)=SEFFER(3,IDIL)+REAL(T)*SIGS/(SIGERD(IDIL)+SIGT)
   65     CONTINUE
          IF(NPAR.EQ.2) THEN
            SIGF=MAX(1.E-9,SIGFF(IGF))
            T=DELF*SIGF
            DO 500 INOR=NOR/2+1,NOR
            MOMP(INOR,2)=MOMP(INOR,2)+T
            T=T*SIGT
  500       CONTINUE
            T=DELF*SIGF/SIGT
            DO 600 INOR=NOR/2,1,-1
            MOMP(INOR,2)=MOMP(INOR,2)+T
            T=T/SIGT
  600       CONTINUE
            DO 650 IDIL=1,NDIL
            T=SIGERD(IDIL)*DELF
            SEFFER(4,IDIL)=SEFFER(4,IDIL)+REAL(T)*SIGF/(SIGERD(IDIL)+
     1      SIGT)
  650       CONTINUE
          ENDIF
        ENDIF
   70 CONTINUE
*
      IF(DEL.EQ.0.0) CALL XABORT('LIBMOM: ALGORITHM FAILURE.')
      DO 80 INOR=1,2*NOR
      MOMT(INOR)=MOMT(INOR)/DEL
   80 CONTINUE
      IF(NPAR.GT.0) THEN
         DO 95 INOR=1,NOR
         DO 90 IPAR=1,NPAR
         MOMP(INOR,IPAR)=MOMP(INOR,IPAR)/REAL(DEL)
   90    CONTINUE
   95    CONTINUE
      ENDIF
      DO 110 IDIL=1,NDIL
      DO 100 IPAR=2,NPAR+2
      SEFFER(IPAR,IDIL)=SEFFER(IPAR,IDIL)/SEFFER(1,IDIL)
  100 CONTINUE
      SEFFER(1,IDIL)=SEFFER(1,IDIL)/REAL(DEL)
  110 CONTINUE
      RETURN
      END
