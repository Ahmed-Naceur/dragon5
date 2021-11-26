*DECK COMRES
      SUBROUTINE COMRES(IMPX,IPISO,IPEDIT,NISOTS,NISO2,ISW,FNORM,ITRES,
     1 PYRES)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute a non-depleting macroscopic residual isotope.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IMPX    print parameter.
* IPISO   pointer to the isotope directory containing the residual
*         isotope.
* IPEDIT  pointer to the edition object (L_EDIT signature).
* NISOTS  number of isotopes in the microlib pointed by IPEDIT.
* NISO2   number of particularized isotopes.
* ISW     selection flag for isotopes in IPEDIT (>0: particularized;
*         <0: included in the macroscopic residual).
* FNORM   flux normalization factor.
*
*Parameters: output
* ITRES   type of the residual isotope (=1: not fissile; =2: fissile).
* PYRES   fission yields of the residual isotope (if ITRES=2).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPISO,IPEDIT
      INTEGER IMPX,NISOTS,NISO2,ISW(NISOTS),ITRES
      REAL FNORM,PYRES(NISO2)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPEDIT,KPEDIT
      INTEGER ISTATE(NSTATE)
      CHARACTER CM*2,TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HVECT
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: HMAKE,HNISO
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIFI,ITYPRO
      REAL, ALLOCATABLE, DIMENSION(:) :: SDEN,PNFTOT,PYIELD,PYDEN,WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAS,SIGS
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WSCAT,PNFIRA,WORK2
*
      CALL LCMGET(IPEDIT,'STATE-VECTOR',ISTATE)
      IF(ISTATE(2).NE.NISOTS) CALL XABORT('COMRES: INVALID EDITION OBJ'
     1 //'ECT.')
      NG=ISTATE(3)
      NL=ISTATE(4)
      NED=ISTATE(13)
      NBESP=ISTATE(16)
      NDEL=ISTATE(19)
      NW=ISTATE(25)
      MAXH=9+3*NW+NL+NED+2*NDEL+NBESP
      ALLOCATE(GAS(NG,MAXH),HMAKE(MAXH+NL),HNISO(NISOTS),SDEN(NISOTS),
     1 WSCAT(NG,NG,NL),PNFIRA(NG,0:NDEL,2),PNFTOT(NISOTS),ITYPRO(NL),
     2 HVECT(NED),WORK(NG+1),WORK2(NG,NG,NL))
      CALL LCMGTC(IPEDIT,'ISOTOPESUSED',12,NISOTS,HNISO)
      CALL LCMGET(IPEDIT,'ISOTOPESDENS',SDEN)
      CALL LCMGTC(IPEDIT,'ADDXSNAME-P0',8,NED,HVECT)
      HMAKE(:MAXH+NL)=' '
      CALL XDRSET(GAS,MAXH*NG,0.0)
      CALL XDRSET(WSCAT,NL*NG*NG,0.0)
      CALL XDRSET(PNFIRA(1,0,2),NG*(NDEL+1),0.0)
      CALL XDRSET(PNFTOT,NISOTS,0.0)
      DENTOT=0.0
      DAWR=0.0
      DECISO=0.0
      ITRES=1
*----
*  SUMMATION OVER NON-PARTICULARIZED ISOTOPES
*----
      IF(NISOTS.GT.0) JPEDIT=LCMGID(IPEDIT,'ISOTOPESLIST')
      DO 170 ISO=1,NISOTS
      IF(ISW(ISO).LT.0) THEN
         ISO3=-ISW(ISO)
         IF(ISO3.GT.NISOTS) CALL XABORT('COMRES: NISOTS OVERFLOW(1).')
         DDEN=SDEN(ISO)
         KPEDIT=LCMGIL(JPEDIT,ISO) ! set ISO-th isotope
         CALL LCMLEN(KPEDIT,'AWR',LENGTH,ITYLCM)
         IF(LENGTH.EQ.1) THEN
*           Compute only heavy-element AWR content
            CALL LCMGET(KPEDIT,'AWR',FLOTT)
            IF(FLOTT.GT.210.0) DAWR=DAWR+DDEN*FLOTT
         ENDIF
         DENTOT=DENTOT+DDEN
         CALL LCMLEN(KPEDIT,'DECAY',LENGTH,ITYLCM)
         IF(LENGTH.EQ.1) THEN
            CALL LCMGET(KPEDIT,'DECAY',FLOTT)
            DECISO=DECISO+FLOTT*DDEN
         ENDIF
*
*        COMPUTE FISSION RATES.
         CALL LCMLEN(KPEDIT,'NFTOT',LENGTH,ITYLCM)
         IF(LENGTH.EQ.NG) THEN
            CALL LCMGET(KPEDIT,'NWT0',GAS(1,1))
            CALL LCMGET(KPEDIT,'NFTOT',WORK)
            DO 30 IGR=1,NG
            DEL=WORK(IGR)*GAS(IGR,1)*DDEN
            PNFTOT(ISO3)=PNFTOT(ISO3)+DEL
   30       CONTINUE
         ENDIF
*
*        SET ARRAY HMAKE.
         DO 40 IW=1,MIN(NW+1,10)
            WRITE(TEXT12,'(3HNWT,I1)') IW-1
            CALL LCMLEN(KPEDIT,TEXT12,LENGTH,ITYLCM)
            IF(LENGTH.EQ.NG) HMAKE(IW)=TEXT12
            WRITE(TEXT12,'(4HNWAT,I1)') IW-1
            CALL LCMLEN(KPEDIT,TEXT12,LENGTH,ITYLCM)
            IF(LENGTH.EQ.NG) HMAKE(1+NW+IW)=TEXT12
            WRITE(TEXT12,'(4HNTOT,I1)') IW-1
            CALL LCMLEN(KPEDIT,TEXT12,LENGTH,ITYLCM)
            IF(LENGTH.EQ.NG) HMAKE(2+2*NW+IW)=TEXT12
   40    CONTINUE
         IOF=3+3*NW
         DO 50 IL=0,NL-1
         IOF=IOF+1
         WRITE (CM,'(I2.2)') IL
         CALL LCMLEN(KPEDIT,'SIGS'//CM,ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF)='SIGS'//CM
   50    CONTINUE
         IOF=IOF+1
         CALL LCMLEN(KPEDIT,'NUSIGF',ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF)='NUSIGF'
         DO 60 IED=1,NED
         IOF=IOF+1
         CALL LCMLEN(KPEDIT,HVECT(IED),ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF)=HVECT(IED)
   60    CONTINUE
         CALL LCMLEN(KPEDIT,'H-FACTOR',ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF+1)='H-FACTOR'
         CALL LCMLEN(KPEDIT,'OVERV',ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF+2)='OVERV'
         CALL LCMLEN(KPEDIT,'TRANC',ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF+3)='TRANC'
         IOF=IOF+3
         DO 70 IDEL=1,NDEL
         IOF=IOF+1
         WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
         CALL LCMLEN(KPEDIT,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF)=TEXT12
   70    CONTINUE
         IOF=IOF+1
         CALL LCMLEN(KPEDIT,'CHI',ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF)='CHI'
         DO 80 IDEL=1,NDEL
         IOF=IOF+1
         WRITE(TEXT12,'(3HCHI,I2.2)') IDEL
         CALL LCMLEN(KPEDIT,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF)=TEXT12
   80    CONTINUE
         DO 85 ISP=1,NBESP
         IOF=IOF+1
         WRITE(TEXT12,'(5HCHI--,I2.2)') ISP
         CALL LCMLEN(KPEDIT,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF)=TEXT12
   85    ENDDO
         IOF=IOF+1
         CALL LCMLEN(KPEDIT,'STRD',ILONG,ITYLCM)
         IF(ILONG.EQ.NG) HMAKE(IOF)='STRD'
         IF(IOF.NE.MAXH) CALL XABORT('COMRES: WRONG OFFSET.')
*
         DO 150 J=1,MAXH
         IF(HMAKE(J).NE.' ') THEN
           CALL LCMLEN(KPEDIT,HMAKE(J),ILONG,ITYLCM)
           IF(ILONG.GT.0) THEN
             CALL LCMGET(KPEDIT,HMAKE(J),WORK)
             IF(HMAKE(J).EQ.'NUSIGF') THEN
               DO 90 IGR=1,NG
               DEL=WORK(IGR)*GAS(IGR,1)*DDEN
               PNFIRA(IGR,0,1)=DEL
               PNFIRA(IGR,0,2)=PNFIRA(IGR,0,2)+DEL
               GAS(IGR,J)=GAS(IGR,J)+WORK(IGR)*DDEN
   90          CONTINUE
             ELSE IF(HMAKE(J)(:3).EQ.'NUS') THEN
               IDEL=J-(7+3*NW+NL+NED)
               DO 100 IGR=1,NG
               DEL=WORK(IGR)*GAS(IGR,1)*DDEN
               PNFIRA(IGR,IDEL,1)=DEL
               PNFIRA(IGR,IDEL,2)=PNFIRA(IGR,IDEL,2)+DEL
               GAS(IGR,J)=GAS(IGR,J)+WORK(IGR)*DDEN
  100          CONTINUE
             ELSE IF(HMAKE(J)(:3).EQ.'NWT') THEN
               DO 110 IGR=1,NG
               GAS(IGR,J)=WORK(IGR)*FNORM
  110          CONTINUE
             ELSE IF((HMAKE(J).EQ.'CHI').OR.
     1               (HMAKE(J)(:5).EQ.'CHI--')) THEN
               DO 120 IGR=1,NG
               GAS(IGR,J)=GAS(IGR,J)+WORK(IGR)*PNFIRA(IGR,0,1)
  120          CONTINUE
             ELSE IF(HMAKE(J)(:3).EQ.'CHI') THEN
               IDEL=J-(8+3*NW+NL+NED+NDEL)
               DO 130 IGR=1,NG
               GAS(IGR,J)=GAS(IGR,J)+WORK(IGR)*PNFIRA(IGR,IDEL,1)
  130          CONTINUE
             ELSE
               DO 140 IGR=1,NG
                 GAS(IGR,J)=GAS(IGR,J)+WORK(IGR)*DDEN
  140          CONTINUE
             ENDIF
           ENDIF
         ENDIF
  150    CONTINUE
         CALL LCMLEN(KPEDIT,'SCAT-SAVED',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
           ALLOCATE(SIGS(NG,NL))
           CALL XDRLGS(KPEDIT,-1,IMPX,0,NL-1,1,NG,SIGS,WORK2,ITYPRO)
           DEALLOCATE(SIGS)
           DO 162 IL=1,NL
           WRITE (CM,'(I2.2)') IL-1
           IF(ITYPRO(IL).NE.0) HMAKE(MAXH+IL)='SCAT'//CM
           DO 161 JGR=1,NG
           DO 160 IGR=1,NG
           WSCAT(IGR,JGR,IL)=WSCAT(IGR,JGR,IL)+WORK2(IGR,JGR,IL)*DDEN
  160          CONTINUE
  161          CONTINUE
  162          CONTINUE
         ENDIF
      ENDIF
  170 CONTINUE
*----
*  SUMMATION OF FISSION YIELDS OVER PARTICULARIZED ISOTOPES
*----
      ALLOCATE(PYDEN(NISO2))
      CALL XDRSET(PYRES,NISO2,0.0)
      CALL XDRSET(PYDEN,NISO2,0.0)
      DO 190 ISO=1,NISOTS
      IF(ISW(ISO).GT.0) THEN
        ISO2=ISW(ISO)
        IF(ISO2.GT.NISO2) CALL XABORT('COMRES: NISO2 OVERFLOW.')
        KPEDIT=LCMGIL(JPEDIT,ISO) ! set ISO-th isotope
        CALL LCMLEN(KPEDIT,'PIFI',NDFI,ITYLCM)
        IF(NDFI.GT.0) THEN
          ALLOCATE(IPIFI(NDFI),PYIELD(NDFI))
          CALL LCMGET(KPEDIT,'PIFI',IPIFI)
          CALL LCMGET(KPEDIT,'PYIELD',PYIELD)
          DO 180 I=1,NDFI
            IFI=IPIFI(I)
            IF(IFI.GT.NISOTS) CALL XABORT('COMRES: NISOTS OVERFLOW(2).')
            IF(ISW(IFI).LT.0) THEN
*             fissile isotope included in the macroscopic residual
              PYRES(ISO2)=PYRES(ISO2)+PYIELD(I)*PNFTOT(-ISW(IFI))
              PYDEN(ISO2)=PYDEN(ISO2)+PNFTOT(-ISW(IFI))
            ENDIF
  180     CONTINUE
          DEALLOCATE(PYIELD,IPIFI)
        ENDIF
      ENDIF
  190 CONTINUE
      DO 200 ISO2=1,NISO2
      IF(PYDEN(ISO2).NE.0.0) PYRES(ISO2)=PYRES(ISO2)/PYDEN(ISO2)
  200 CONTINUE
      DEALLOCATE(PYDEN)
*----
*  SAVE THE RESIDUAL ISOTOPE IN IPISO ISOTOPE OBJECT
*----
      IF(DENTOT.GT.0.0) THEN
         IF(IMPX.GT.0) WRITE (6,600)
         TEXT12='*MAC*RES0001'
         CALL LCMPTC(IPISO,'ALIAS',12,1,TEXT12)
         CALL LCMPUT(IPISO,'AWR',1,2,DAWR)
         DECISO=DECISO/DENTOT
         IF(DECISO.GT.0.0) CALL LCMPUT(IPISO,'DECAY',1,2,DECISO)
         DO 240 J=1,MAXH
         IF(HMAKE(J).EQ.'NUSIGF') ITRES=2
         IF(HMAKE(J).EQ.'OVERV') THEN
            DO 210 IGR=1,NG
            GAS(IGR,J)=GAS(IGR,J)/DENTOT
  210       CONTINUE
         ELSE IF((HMAKE(J).EQ.'CHI').OR.(HMAKE(J)(:5).EQ.'CHI--')) THEN
            DO 220 IGR=1,NG
            IF(GAS(IGR,J).NE.0.0) THEN
               GAS(IGR,J)=GAS(IGR,J)/PNFIRA(IGR,0,2)
            ENDIF
  220       CONTINUE
         ELSE IF(HMAKE(J)(:3).EQ.'CHI') THEN
            IDEL=J-(8+3*NW+NL+NED+NDEL)
            DO 230 IGR=1,NG
            IF(GAS(IGR,J).NE.0.0) THEN
               GAS(IGR,J)=GAS(IGR,J)/PNFIRA(IGR,IDEL,2)
            ENDIF
  230       CONTINUE
         ENDIF
         IF((HMAKE(J).NE.' ').AND.(HMAKE(J)(:4).NE.'SIGS')) THEN
            CALL LCMPUT(IPISO,HMAKE(J),NG,2,GAS(1,J))
         ENDIF
  240    CONTINUE
         IF(ITYPRO(1).NE.0) THEN
            DO 250 IL=1,NL
            ITYPRO(IL)=0
            IF(HMAKE(MAXH+IL).NE.' ') ITYPRO(IL)=1
  250       CONTINUE
            CALL XDRLGS(IPISO,1,IMPX,0,NL-1,1,NG,GAS(1,4+3*NW),WSCAT,
     1      ITYPRO)
         ENDIF
*
         IF(IMPX.GT.3) THEN
           WRITE(6,'(/17H NUMBER DENSITY =,1P,E12.4)') 1.0
           WRITE(6,'(23H WEIGHTED ATOMIC MASS =,1P,E13.5)') DAWR
           DO 260 J=1,MAXH
           IF(HMAKE(J).NE.' ') THEN
             WRITE (6,610) HMAKE(J),(GAS(I,J),I=1,NG)
           ENDIF
  260      CONTINUE
           WRITE (6,610) 'SIGA    ',(GAS(I,3+2*NW)-GAS(I,4+3*NW),I=1,NG)
           WRITE (6,610) 'SIGW00  ',(WSCAT(I,I,1),I=1,NG)
           IF(NL.GT.1) THEN
             WRITE (6,610) 'SIGW01  ',(WSCAT(I,I,2),I=1,NG)
           ENDIF
         ENDIF
      ENDIF
      DEALLOCATE(WORK2,WORK,HVECT,ITYPRO,PNFTOT,PNFIRA,WSCAT,SDEN,HNISO,
     1 HMAKE,GAS)
      RETURN
*
  600 FORMAT (//49H COMRES: CROSS SECTION OF MACRO RESIDUAL ISOTOPE ,
     1 8H*MAC*RES)
  610 FORMAT (/11H REACTION ',A12,2H':/(1X,1P,10E12.4))
      END
