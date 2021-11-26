*DECK AFMCPT
      SUBROUTINE AFMCPT (KENTRY,NBURN,NGRP,NISO,NL,IMPX,
     1 SMACB,XBORB,XPURB,XXENB,XT1FB,XT2FB,XT1CB,
     1 XT2CB,XT1MB,XT2MB,XD1CB,XD2CB,XD1MB,XD2MB,
     1 XSMB,XNP9B,XMFDB,XMMDB,XPF1B,XPF2B,XPF1LB,XPF2LB,
     1 DENSITB,CPW1B,CPW2B,FLUXB,OVERVB,CHIB,
     1 IJJ,NJJ,HISO,CTITRE,NMIX,SIGMA,NTYP,TF,TC,
     1 TM,DC,DM,BOR,XEN,SM,RNP9,XI,TFR,TCR,TMR,XIR,
     1 OVERV,FLUX,CHI,SCAT,MMIX,NPS,PW,XBRH,XBURN,
     1 LTAV,IRAV,IDF,JTAB,IXYZ,ILIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the cross sections
*
*Copyright:
* Copyright (C) 1996 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* M.T. Sissaoui
*
*Parameters: input
* KENTRY  address of the linked list or handle to the xsm file.
* NGRP    number of energy groups.
* NISO    number of extracted isotopes.
* NL      number of Legendre orders (=1 for isotropic scattering).
* NBURN   number of burnup steps.
* NTYP    total number of cross sections type.
* TF      fuel temperature.
* TC      coolant temperature.
* TM      moderator temperature.
* DC      coolant density.
* DM      moderator density.
* BOR     Boron concentration.
* XEN     Xenon concentration.
* SM      Samarium concentration.
* RNP9    Neptunium concentration.
* XI      moderator purity.
* NMIX    mixture number.
* NPS     NPS-2 power shift.
* IXYZ    type of diffusion coefficient (=0: isotropic; =1: directional)
*
*Parameters: 
* IMPX     
* SMACB    
* XBORB    
* XPURB    
* XXENB    
* XT1FB    
* XT2FB    
* XT1CB    
* XT2CB    
* XT1MB    
* XT2MB    
* XD1CB    
* XD2CB    
* XD1MB    
* XD2MB    
* XSMB     
* XNP9B    
* XMFDB    
* XMMDB    
* XPF1B    
* XPF2B    
* XPF1LB   
* XPF2LB   
* DENSITB  
* CPW1B    
* CPW2B    
* FLUXB    
* OVERVB   
* CHIB     
* IJJ      
* NJJ      
* HISO     
* CTITRE   
* SIGMA    
* TFR      
* TCR      
* TMR      
* XIR      
* OVERV      
* FLUX       
* CHI        
* SCAT       
* MMIX       
* PW           
* XBRH     
* XBURN    
* LTAV     
* IRAV     
* IDF      
* JTAB     
* ILIN     
*
*-----------------------------------------------------------------------
*
      CHARACTER HMICRO*12,CTITRE*72
      LOGICAL LTAV
      DOUBLE PRECISION XCOF(3)
      REAL CPF1(3)
      DIMENSION KENTRY(*),
     1 SMACB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XBORB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XXENB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XT1FB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XT2FB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XT1CB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XT2CB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XT1MB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XT2MB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XD1CB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XD2CB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XD1MB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XD2MB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XSMB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XNP9B(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XMFDB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XMMDB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XPF1B(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XPF2B(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XPF1LB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XPF2LB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 XPURB(NGRP*NGRP,NTYP,NISO,NBURN,*),
     1 DENSITB(NISO,NBURN,*),CPW1B(2,NBURN,*),
     1 CPW2B(2,NBURN,*),FLUXB(NGRP,NBURN,*),
     1 CHIB(NGRP,NBURN,*),OVERVB(NGRP,NBURN,*),
     2 HISO(*),JTAB(*),
     6 ELMT(3),PW(*),
     7 SIGMA(MMIX,NGRP,NTYP),XBRH(*),
     8 OVERV(MMIX,*),FLUX(MMIX,*),
     7 CHI(MMIX,*),SCAT(MMIX,NL,NGRP,*),
     1 XBURN(NBURN,*),NJJ(*),IJJ(*)
      REAL, ALLOCATABLE, DIMENSION(:) :: DEL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: CPW1,CPW2
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SMAC
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SMAC(NGRP*NGRP,NTYP,NISO),DEL(NISO),CPW1(3,NPS),
     1 CPW2(3,NPS))
*
*     Remove warning for uninitialized variable
      ISOB=0
      D2MB=0.0
      D2CB=0.0
      D1M=0.0
      D2C=0.0
      CX=0.0
      CSM=0.0
      D1C=0.0
      D2M=0.0
      CNP9=0.0
      CB=0.0
      PUR=0.0
*---------------------------------------------------------------*
      IF(NBURN.EQ.0) CALL XABORT('AFMCPT: ZERO NUMBER OF MIXTURES.')
      IF(NGRP.EQ.0) CALL XABORT('AFMCPT: ZERO NUMBER OF GROUPS.')
*---------------------------------------------------------------*
      NG2=NGRP*NGRP
      NTM=4+2*IXYZ
      IF(IMPX.GT.5) THEN
         WRITE(6,'(9H AFMCPT: ,A)') CTITRE
         WRITE(6,*) ' NTYP ',NTYP
      ENDIF
      XECON=0.0
      TFAV=0.0
      DCAV=0.0
      CPF2=0.0
      CF=0.0
      DO 81 ISO=1,NISO
       DO 80 ITY=1,NTYP
        DO 79 IGR=1,NG2
          SMAC(IGR,ITY,ISO)=0.0
 79     CONTINUE
 80    CONTINUE
 81   CONTINUE
      DO 82 IGR=1,NGRP
         FLUX(NMIX,IGR)=0.0
         OVERV(NMIX,IGR)=0.0
         CHI(NMIX,IGR)=0.0
 82   CONTINUE
      DO 161 ISO=1,NISO
         DEL(ISO)=0.0
 161  CONTINUE
      DO 73 II=1,3
         XCOF(II)=0.0D0
         ELMT(II)=0.0
         CPF1(II)=0.0
 73   CONTINUE
      DO 74 K=1,NPS
         DO II=1,3
            CPW1(II,K)=0.0
            CPW2(II,K)=0.0
         ENDDO
 74   CONTINUE
*
      PREF=PW(1)
      IPFBM=KENTRY(2)
*---------------------------------------------------------------*
*     COMPUTE THE HISTORY COEFFICIENTS
      NPSX=NPS+NPS-1
      DO 101 K=2,NPSX
* XSECTION FOR SNAP-SHOT OR (TIME AVERAGE-HOMOG.)
        IF(LTAV) THEN
          IRMAX=IRAV
          IRMIN=IRAV
          XCOF(1)=1.0D0
        ELSE
*
          IF(K.GT.NPS) THEN
            XIRAD=XBRH(K-NPS+1)
          ELSE
            XIRAD=XBRH(NPS)-XBRH(K-1)
          ENDIF
           CALL AFMLOC(NBURN,NTP,XIRAD,XIRAD,XBURN(1,IDF),
     1                IRMAX,IRMIN,XCOF,ILIN)
        ENDIF
*
* INTERPOLATE THE HISTORY COEFFICIENTS
        IH=0
        DO 910 I = IRMIN,IRMAX
          IH=IH+1
          IF(K.LE.NPS) THEN
            IF(PW(K).GT.PREF) THEN
              CPW1(IH,K)=CPW1B(1,I,IDF) + CPW1(IH,K)
              CPW2(IH,K)=CPW2B(1,I,IDF) + CPW2(IH,K)
            ELSE
              CPW1(IH,K)=CPW1B(2,I,IDF) + CPW1(IH,K)
              CPW2(IH,K)=CPW2B(2,I,IDF) + CPW2(IH,K)
            ENDIF
*
            IF(K.GT.2) THEN
               IF(PW(K-1).GT.PREF) THEN
                 CPW1(IH,K-1)=-CPW1B(1,I,IDF) + CPW1(IH,K-1)
                 CPW2(IH,K-1)=-CPW2B(1,I,IDF) + CPW2(IH,K-1)
               ELSE
                 CPW1(IH,K-1)=-CPW1B(2,I,IDF) + CPW1(IH,K-1)
                 CPW2(IH,K-1)=-CPW2B(2,I,IDF) + CPW2(IH,K-1)
               ENDIF
            ENDIF
          ENDIF
 910    CONTINUE
 101  CONTINUE
*
      YF=0.0
      DO 111 K=2,NPS
        IF(K.EQ.2) THEN
* CORRECTE THE STURATING PSEUDO-FISSILE ISOTOPE
          IF(PW(K).GT.PREF) THEN
            XPW=ALOG(PW(K)/PW(1))
            XPWM=1.0/PW(K)-1.0/PW(1)
            YF=0.0
          ELSE
            XPW=PW(K)-PW(1)
            XPWM=(PW(K)-PW(1))**2
            YF=1.0
          ENDIF
          DO IH=1,3
          CPF1(IH)=CPF1(IH) +CPW1(IH,K)*XPW + CPW2(IH,K)*XPWM
          ENDDO
C
        ELSE
          IF(PW(K).GT.PREF) THEN
            XPW=ALOG(PW(K)/PW(1))
            XPWM=1.0/PW(K)-1.0/PW(1)
          ELSE
            XPW=PW(K)-PW(1)
            XPWM=(PW(K)-PW(1))**2
C
          ENDIF
          DO IH=1,3
          CPF1(IH)=CPF1(IH) +CPW1(IH,K)*XPW + CPW2(IH,K)*XPWM
          ENDDO
        ENDIF
 111  CONTINUE
*---------------------------------------------------------------*
*         COMPUTE THE DEVIATION OF THE PSEUDO-ISOTOPE
*         CPF1 AND CPF2
*      CPF2=CPF1*CPF1
*---------------------------------------------------------------*
*  APPLY  THE FEEDBACK MODEL
      T1F=SQRT(TF)-SQRT(TFR)
      T2F=TF-TFR
      T1C=ALOG(TC/TCR)
      T2C=1.0/TC - 1.0/TCR
      T1M=ALOG(TM/TMR)
      T2M=1.0/TM - 1.0/TMR
*
* RECOVER LOCAL PARAMETER COEFFICIENT AND X-SECTIONS
      II=0
      DMOD=0.0
      DO 900 I = IRMIN,IRMAX
        II=II+1
        RXCOEF=REAL(XCOF(II))
        DO 249 IGR=1,NGRP
          FLUX(NMIX,IGR)=FLUX(NMIX,IGR)+RXCOEF*FLUXB(IGR,I,IDF)
          OVERV(NMIX,IGR)=OVERV(NMIX,IGR)+RXCOEF*OVERVB(IGR,I,IDF)
 249    CONTINUE
* COMPUTE DELTA-CONCENTRATION
        DO 49 ISO=1,NISO
          IF(DENSITB(ISO,I,IDF).EQ.0.0) GO TO 49
          WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
          IF(HMICRO.EQ.'XE135') THEN
            DEL(ISO)=XEN/(IRMAX-IRMIN+1)-DENSITB(ISO,I,IDF)*RXCOEF
     1               +DEL(ISO)
            CX=DEL(ISO)
          ELSE IF(HMICRO.EQ.'BMOD') THEN
            ISOB=ISO
*
*         ERROR IN CALCULATING BORON CONCENTRATION,
*         BUT ADEQUATELY APPROXIMATED IF THE BORON REF EQ 0.0ppm
*             
*            DEL(ISO)=BOR/(IRMAX-IRMIN+1)-
*     1               DENSITB(ISO,I,IDF)*RXCOEF+DEL(ISO)
            DEL(ISO)=-DENSITB(ISO,I,IDF)*RXCOEF+DEL(ISO)
            CB=DEL(ISO)
          ELSE IF(HMICRO.EQ.'CWAT') THEN
            DEL(ISO)=(DC-1.0)*DENSITB(ISO,I,IDF)
            D1C=DEL(ISO)
            D2C=D1C*D1C
            D2CB=D1C*D1C
          ELSE IF(HMICRO.EQ.'MWAT') THEN
            DMOD=DENSITB(ISO,I,IDF)/(IRMAX-IRMIN+1)+DMOD
            DEL(ISO)=(DM-1.0)*DENSITB(ISO,I,IDF)
            D1M=ALOG(DM)
            D2M=1.0/(DM*DENSITB(ISO,I,IDF)) - 1.0/DENSITB(ISO,I,IDF)
            D2MB=DEL(ISO)
*         PURITY
            PUR=(XI-XIR)*DM*DENSITB(ISO,I,IDF)
          ELSE IF(HMICRO.EQ.'SM149') THEN
            DEL(ISO)=SM/(IRMAX-IRMIN+1)-DENSITB(ISO,I,IDF)*RXCOEF
     1               +DEL(ISO)
            CSM=DEL(ISO)
          ELSE IF(HMICRO.EQ.'NP239') THEN
            DEL(ISO)=RNP9/(IRMAX-IRMIN+1)-DENSITB(ISO,I,IDF)*RXCOEF
     1               +DEL(ISO)
            CNP9=DEL(ISO)
          ELSE IF(HMICRO.EQ.'FPC') THEN
*
            CF= CF+DENSITB(ISO,I,IDF)*RXCOEF
            DEL(ISO)=CPF1(II)*RXCOEF+DEL(ISO)
          ELSE IF(HMICRO.EQ.'MACR') THEN
            DEL(ISO)=DENSITB(ISO,I,IDF)
            DO 271 IGR=1,NGRP
              CHI(NMIX,IGR)=CHIB(IGR,I,IDF)*RXCOEF+CHI(NMIX,IGR)
  271       CONTINUE
          ENDIF
   49   CONTINUE
  900 CONTINUE
* R.C. 24/05/2011
* boron unit correction
*    Bnat = 10.811  g/mol
*    O16  = 15.9949 g/mol
*    H1   =  1.0078 g/mol
*    D2   =  2.0141 g/mol
      CB=DMOD*BOR/10.811*((1.0078*(1-XI)+2.0141*XI)*2+15.9949)/3+CB
      CB=CB*DM
      DEL(ISOB)=CB
* R.C.
*  CORRECT THE FUEL TEMPERATURE
      CQ=REAL((CPF1(1)*XCOF(1)+CPF1(2)*XCOF(2)+CPF1(3)*XCOF(3)+CF)/CF)
      RCQ=1.0-CQ
      CQ2=CQ*CQ
      RCQ2=1.0-CQ2
      IF(JTAB(1).EQ.0) THEN
        CQ=0.0
        RCQ=0.0
        CQ2=0.0
        RCQ2=0.0
        CX=0.0
        CSM=0.0
        CNP9=0.0
        D1C=0.0
        D2C=0.0
        D2CB=0.0
        T1F=0.0
        T2F=0.0
        T1C=0.0
        T2C=0.0
      ENDIF
*
* RECOVER MACROSCOPIC X-SECTIONS
      II=0
      DO 901 I = IRMIN,IRMAX
        II=II+1
        RXCOEF=REAL(XCOF(II))
        CPF2=CPF1(II)*CPF1(II)
        DO 98 ISO=1,NISO
         DO 99  ITY=1,NTYP
           IF(ISO.EQ.1) THEN
             ZZ=1.0
             IF(ITY.EQ.1) THEN
               VD2M=D2MB
               VD2C=D2CB
             ELSE IF(ITY.GE.5) THEN
               VD2M=D2MB
               VD2C=D2C
             ELSE
               VD2M=D2M
               VD2C=D2C
             ENDIF
           ELSE
             ZZ=0.0
             VD2M=D2M
             VD2C=D2C
           ENDIF
           DO 100 IGR=1,NG2
             SMAC(IGR,ITY,1)=SMAC(IGR,ITY,1)+
     1        SMACB(IGR,ITY,1,I,IDF)*RXCOEF*ZZ+
     1       (XBORB(IGR,ITY,ISO,I,IDF)*RXCOEF*CB +
     1       XPURB(IGR,ITY,ISO,I,IDF)*RXCOEF*PUR*ZZ +
     1       XXENB(IGR,ITY,ISO,I,IDF)*RXCOEF*CX +
     1       XT1FB(IGR,ITY,ISO,I,IDF)*RXCOEF*T1F*CQ +
     1       XT1FB(IGR,ITY,ISO,1,IDF)*RXCOEF*T1F*RCQ +
     1       XT2FB(IGR,ITY,ISO,I,IDF)*RXCOEF*T2F*CQ2 +
     1       XT2FB(IGR,ITY,ISO,1,IDF)*RXCOEF*T2F*RCQ2 +
     1       XT1CB(IGR,ITY,ISO,I,IDF)*RXCOEF*T1C +
     1       XT2CB(IGR,ITY,ISO,I,IDF)*RXCOEF*T2C +
     1       XT1MB(IGR,ITY,ISO,I,IDF)*RXCOEF*T1M +
     1       XT2MB(IGR,ITY,ISO,I,IDF)*RXCOEF*T2M +
     1       XD1CB(IGR,ITY,ISO,I,IDF)*RXCOEF*D1C +
     1       XD2CB(IGR,ITY,ISO,I,IDF)*RXCOEF*VD2C +
     1       XD1MB(IGR,ITY,ISO,I,IDF)*RXCOEF*D1M +
     1       XD2MB(IGR,ITY,ISO,I,IDF)*RXCOEF*VD2M +
     1       XSMB(IGR,ITY,ISO,I,IDF)*RXCOEF*CSM +
     1       XNP9B(IGR,ITY,ISO,I,IDF)*RXCOEF*CNP9 +
     1       XMFDB(IGR,ITY,ISO,I,IDF)*RXCOEF*T1F*D1C +
     1       XMMDB(IGR,ITY,ISO,I,IDF)*RXCOEF*T1C*D1C +
     1       XPF1B(IGR,ITY,ISO,I,IDF)*RXCOEF*CPF1(II)*(1.-YF) +
     1       XPF2B(IGR,ITY,ISO,I,IDF)*RXCOEF*CPF2*(1.-YF)+
     1       XPF1LB(IGR,ITY,ISO,I,IDF)*RXCOEF*CPF1(II)*YF +
     1       XPF2LB(IGR,ITY,ISO,I,IDF)*RXCOEF*CPF2*YF)*DEL(ISO)
*
 100       CONTINUE
  99     CONTINUE
  98    CONTINUE
 901   CONTINUE
* STORE SCATTERING
       IL= 1
       ITY=5+2*IXYZ+IL
       IGAR=0
       DO 130 JGR=1,NGRP
          DO 120 IGR=IJJ(JGR),IJJ(JGR)-NJJ(JGR)+1,-1
            IGAR=IGAR+1
            SCAT(NMIX,IL,IGR,JGR)=SMAC(IGAR,ITY,1)
* TOTAL OR ABS
            SMAC(IGR,2,1)=SMAC(IGR,2,1)+SCAT(NMIX,IL,IGR,JGR)
 120     CONTINUE
 130  CONTINUE
* STORE X-SECTIONS
      DO 261 ITY=1,NTYP
        DO 260 IGR=1,NGRP
          SIGMA(NMIX,IGR,ITY)=SMAC(IGR,ITY,1)
 260    CONTINUE
 261  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(CPW2,CPW1,DEL,SMAC)
      RETURN
      END
