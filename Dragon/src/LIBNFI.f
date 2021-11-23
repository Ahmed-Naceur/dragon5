*DECK LIBNFI
      SUBROUTINE LIBNFI(IPLIB,NGRO,NBISO,NBMIX,NDEL,NESP,IPISO,MIX,
     1 MAXNFI,NFISSI,LSAME)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the maximum number of fissionable isotopes in a mixture.
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
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NBMIX   number of mixtures present in the calculation domain.
* NDEL    number of delayed precursor groups.
* NESP    number of energy-dependent fission spectra.
* IPISO   pointer array towards microlib isotopes.
* MIX     mixture number of each isotope (can be zero for void).
* MAXNFI  second dimension of array INDFIS.
*
*Parameters: output
* NFISSI  maximum number of fissionable isotopes in a mixture.
* LSAME   fission spectrum mask (=.true. if all the isotopes have the
*         same fission spectrum and the same precursor group decay
*         constants.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NGRO,NBISO,NBMIX,NDEL,NESP,MIX(NBISO),MAXNFI,NFISSI
      LOGICAL LSAME
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIB
      INTEGER MAXGRO,NSTATE
      PARAMETER (MAXGRO=50,NSTATE=40)
      CHARACTER HSMG*131,TEXT12*12
      REAL CHI2(MAXGRO),LAM1(MAXGRO),LAM2(MAXGRO)
      INTEGER IDATA(NSTATE),ISOT,IBM,IFIS,IGR,ILONG,ITYLCM,IWFIS,JBM,
     1 KFIS,LENGT1,LENGT2,LENGTZ
      LOGICAL LFISS
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWRK
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDFIS
      REAL, ALLOCATABLE, DIMENSION(:) :: CHI1
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDFIS(NBMIX,MAXNFI),CHI1(NGRO))
*
      NFISSI=0
      CALL LCMLEN(IPLIB,'MACROLIB',ILONG,ITYLCM)
      LSAME=(NGRO.LE.MAXGRO).AND.(NDEL.LE.MAXGRO)
      IF(ILONG.EQ.-1) THEN
         CALL LCMSIX(IPLIB,'MACROLIB',1)
         CALL LCMGTC(IPLIB,'SIGNATURE',12,1,TEXT12)
         IF(TEXT12.NE.'L_MACROLIB') THEN
            CALL XABORT('LIBNFI: INVALID SIGNATURE ON THE MACROLIB.')
         ENDIF
         CALL LCMGET(IPLIB,'STATE-VECTOR',IDATA)
         IF(IDATA(1).NE.NGRO) THEN
            WRITE(HSMG,'(38HLIBNFI: EXISTING MACROLIB HAVE NGROUP=,I4,
     1      26H NEW MACROLIB HAVE NGROUP=,I4,1H.)') IDATA(1),NGRO
            CALL XABORT(HSMG)
         ELSE IF(IDATA(2).GT.NBMIX) THEN
            WRITE(HSMG,'(37HLIBNFI: EXISTING MACROLIB HAVE NBMIX=,I4,
     1      25H NEW MACROLIB HAVE NBMIX=,I4,1H.)') IDATA(2),NBMIX
            CALL XABORT(HSMG)
         ELSE IF(IDATA(4).GT.NBISO*NESP) THEN
            WRITE(HSMG,'(38HLIBNFI: EXISTING MACROLIB HAVE NFISSI=,I4,
     1      13H GREATER THAN,I5,1H.)') IDATA(4)/NESP,NBISO
            CALL XABORT(HSMG)
         ENDIF
         NFISSI=IDATA(4)/NESP
         LSAME=LSAME.AND.(NFISSI.LE.1)
         IF(NFISSI.GT.0) THEN
            CALL LCMLEN(IPLIB,'FISSIONINDEX',ILONG,ITYLCM)
            IF(ILONG.EQ.0) THEN
*              THE NAMES ARE NOT DEFINED.
               DO 15 IFIS=1,NFISSI
               DO 10 IBM=1,NBMIX
               INDFIS(IBM,IFIS)=0
   10          CONTINUE
   15          CONTINUE
            ELSE IF(ILONG.EQ.NFISSI*NBMIX) THEN
               CALL LCMGET(IPLIB,'FISSIONINDEX',INDFIS)
            ELSE IF(ILONG.LT.NFISSI*NBMIX) THEN
*              REORDER THE 'FISSIONINDEX' MATRIX.
               ALLOCATE(IWRK(ILONG))
               CALL LCMGET(IPLIB,'FISSIONINDEX',IWRK)
               DO 31 IFIS=1,NFISSI
               DO 20 IBM=1,IDATA(2)
               INDFIS(IBM,IFIS)=IWRK((IFIS-1)*IDATA(2)+IBM)
   20          CONTINUE
               DO 30 IBM=IDATA(2)+1,NBMIX
               INDFIS(IBM,IFIS)=0
   30          CONTINUE
   31          CONTINUE
               DEALLOCATE(IWRK)
            ELSE
               CALL XABORT('LIBNFI: INVALID NUMBER OF MIXTURES.')
            ENDIF
         ENDIF
         CALL LCMSIX(IPLIB,' ',2)
      ENDIF
      DO 100 ISOT=1,NBISO
      IBM=MIX(ISOT)
      IF(IBM.GT.0) THEN
         JPLIB=IPISO(ISOT)
         IF(C_ASSOCIATED(JPLIB)) THEN
            CALL LCMLEN(JPLIB,'NUSIGF',ILONG,ITYLCM)
            IF(NESP.EQ.1) THEN
               CALL LCMLEN(JPLIB,'CHI',LENGTZ,ITYLCM)
            ELSE
               CALL LCMLEN(JPLIB,'CHI--01',LENGTZ,ITYLCM)
            ENDIF
            IF((ILONG.GT.0).AND.(LENGTZ.GT.0)) THEN
               IF(NESP.EQ.1) THEN
                  CALL LCMGET(JPLIB,'CHI',CHI1)
               ELSE
                  CALL LCMGET(JPLIB,'CHI--01',CHI1)
               ENDIF
               LFISS=.FALSE.
               DO 35 IGR=1,NGRO
               LFISS=LFISS.OR.(CHI1(IGR).GT.0.0)
   35          CONTINUE
               IF(.NOT.LFISS) GO TO 100
               IF(LSAME) THEN
                  CALL LCMLEN(JPLIB,'LAMBDA-D',LENGT1,ITYLCM)
                  IF((LENGT1.EQ.NDEL).AND.(NDEL.GT.0)) THEN
                     CALL LCMGET(JPLIB,'LAMBDA-D',LAM1)
                  ENDIF
               ENDIF
               DO 40 IFIS=1,NFISSI
               IWFIS=INDFIS(IBM,IFIS)
               IF((IWFIS.EQ.ISOT).OR.(IWFIS.EQ.0)) THEN
                  KFIS=IFIS
                  GO TO 90
               ENDIF
   40          CONTINUE
               IF(LSAME) THEN
                  DO 70 IFIS=1,NFISSI
                  IWFIS=INDFIS(IBM,IFIS)
                  JPLIB=IPISO(IWFIS)
                  CALL LCMGET(JPLIB,'CHI',CHI2)
                  DO 50 IGR=1,NGRO
                  LSAME=LSAME.AND.(ABS(CHI1(IGR)-CHI2(IGR)).LE.1.0E-3)
   50             CONTINUE
                  CALL LCMLEN(JPLIB,'LAMBDA-D',LENGT2,ITYLCM)
                  IF((LENGT1.EQ.NDEL).AND.(LENGT2.EQ.NDEL)
     1                               .AND.(NDEL.GT.0)) THEN
                     CALL LCMGET(JPLIB,'LAMBDA-D',LAM2)
                     DO 60 IGR=1,NDEL
                     LSAME=LSAME.AND.(LAM1(IGR).EQ.LAM2(IGR))
   60                CONTINUE
                  ENDIF
   70             CONTINUE
               ENDIF
               NFISSI=NFISSI+1
               IF(NFISSI.GT.MAXNFI) CALL XABORT('LIBNFI: INDFIS OVERFL'
     1         //'OW.')
               KFIS=NFISSI
               DO 80 JBM=1,NBMIX
               INDFIS(JBM,KFIS)=0
   80          CONTINUE
   90          INDFIS(IBM,KFIS)=ISOT
            ENDIF
         ENDIF
      ENDIF
  100 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(CHI1,INDFIS)
      RETURN
      END
