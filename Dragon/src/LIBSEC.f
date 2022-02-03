*DECK LIBSEC
      SUBROUTINE LIBSEC(MAXTRA,LLL,IANNN,NGRO,IX,UUU,DELTA,SIGS,SIG1,
     1 PRI,NLET,STR,DEL,NRSTR,IANIS,ITY,NEXT,NEXU,NEXV,NEXW,III)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the values of the transfer macroscopic cross section for
* secondary neutrons in group LLL. Component of the APOLIB-1 reader.
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
* MAXTRA  available storage for apollo compacted transfer
*          cross sections.
* LLL     group number for secondary neutrons.
* IANNN   type of transport correction (=-1: transport corrected P0;
*         =0: P0; =1: P1).
* NGRO    number of groups.
* IX      number of groups with up-scattering.
* UUU     groups limits in lethargy units.
* DELTA   groups width in lethargy units.
* SIGS    diffusion P0 microscopic cross sections.
* SIG1    diffusion P1 microscopic cross sections.
* PRI     transfer microscopic cross sections.
* DEL     elementary mesh element in lethargy.
* NRSTR   number of cross section structures own by the isotope.
* IANIS   Legendre order corresponding to each cross section structure.
* ITY     type of each cross section structure.
* NEXT    length of each cross section structure.
* NEXU    information related to each cross section structure.
* NEXV    information related to each cross section structure.
* NEXW    information related to each cross section structure.
* III     offset in vector PRI of each cross section structure.
*
*Parameters: output
* NLET    number of down-scattering groups (including group LLL).
* STR     values of the transfer macroscopic cross section:
*          STR(1)     from group LLL;
*          STR(2)     from group LLL-1;
*          STR(LLL)   from group 1;
*          STR(LLL+1) from group NGRO;
*          STR(LLL+2) from group NGRO-1;
*          STR(NGRO)  from group LLL+1.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXTRA,LLL,IANNN,NGRO,IX,NLET,NRSTR,IANIS(80),ITY(80),
     1 NEXT(80),NEXU(80),NEXV(80),NEXW(80),III(80)
      REAL UUU(NGRO),DELTA(NGRO),SIGS(NGRO),SIG1(NGRO),PRI(MAXTRA),
     1 STR(NGRO),DEL
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION DAUX
      LOGICAL SELF
      EQUIVALENCE(AUX1,K12)
*
      NGROIN=NGRO-IX
      IAN=IANNN
      IF(IANNN.EQ.-1)IAN=0
      DO 71 MM=1,NGRO
      STR(MM)=0.0
   71 CONTINUE
      NLET=1
      SELF=.FALSE.
      IF (NRSTR.GT.0) THEN
         DO 250 K=1,NRSTR
         IF (IANIS(K).NE.IAN) GOTO 250
         MML=0
         IF (ITY(K).EQ.1) THEN
*           ELASTIC SLOWING-DOWN MATRIX.
            IF (LLL.GT.NGROIN) GO TO 250
            SELF=.TRUE.
            LDELH=INT(UUU(LLL)/DEL+0.1)
            LARGRL=INT(DELTA(LLL)/DEL+0.1)
            LDELB=LDELH-LARGRL+1
            NP2=NEXT(K)
            IHM=III(K)+NP2-1
            LTES=LDELB-NP2
            ISOTOP=K
            INDICE=1
            DO 210 MM1=1,LLL
            MM=LLL-MM1+1
            MDELH=INT(UUU(MM)/DEL+0.1)
            IF(MDELH.LE.LTES)GOTO 211
            LARGRM=INT(DELTA(MM)/DEL+0.1)
            MDELB=MDELH-LARGRM+1
            MDELB=MAX0(MDELB,LTES+1)
            DAUX=0.0
            LARG=MIN0(LARGRM,LARGRL)
            IF(LARG.GT.4) THEN
               IHAUT=LDELH-MDELB+III(K)
               IHAUT=MIN0(IHAUT,IHM)
               J=0
               INTER2=0
               IF (INDICE.EQ.2) THEN
                  IBAS=LDELB-MDELH+III(K)
                  IBAS=MAX0(IBAS,III(K))
                  LARGLI=IABS(LARGRM-LARGRL)
                  INTER1=IBAS+LARG-2
                  INTER1=MIN0(INTER1,IHAUT)
                  DO 182 I=IBAS,INTER1
                  J=J+1
                  DAUX=DAUX+PRI(I)*FLOAT(J)
  182             CONTINUE
                  INTER1=INTER1+1
                  INTER2=INTER1+LARGLI
                  INTER2=MIN0(IHAUT,INTER2)
                  IF(INTER1.GT.INTER2) GO TO 1004
                  J=LARG
                  DO 183 I=INTER1,INTER2
                  DAUX=DAUX+PRI(I)*FLOAT(LARG)
  183             CONTINUE
               ELSE IF (INDICE.EQ.1) THEN
                  INDICE=2
                  INTER2=III(K)-1
                  J=LARG+1
               ENDIF
               INTER2=INTER2+1
               DO 184 I=INTER2,IHAUT
               J=J-1
               DAUX=DAUX+PRI(I)*FLOAT(J)
  184          CONTINUE
            ELSE
               DO 83 MDEL=MDELB,MDELH
               IBAS=LDELB-MDEL+III(K)
               IHAUT=LDELH-MDEL+III(K)
               IBAS=MAX0(IBAS,III(K))
               IHAUT=MIN0(IHAUT,IHM)
               DO 82 I=IBAS,IHAUT
               DAUX=DAUX+PRI(I)
   82          CONTINUE
   83          CONTINUE
            ENDIF
 1004       STR(MM1)=STR(MM1)+REAL(DAUX*SIGS(MM)*DEL/DELTA(LLL))
  210       CONTINUE
            MM=MM-1
  211       MML=LLL-MM
         ELSE IF (ITY(K).EQ.4) THEN
*           STANDARD GALOCHE.
            IF (LLL.GT.NGROIN) GO TO 250
            SELF=.TRUE.
            NEX1=NEXU(K)
            NEX2=NEXV(K)
            NEX3=NEXW(K)
            IF(LLL.GT.(NEX2+NEX3)) GO TO 801
            IPR=III(K)-1+(LLL*(LLL-1))/2
            DO 802 I=1,LLL
            IPR=IPR+1
            STR(I)=STR(I)+PRI(IPR)
  802       CONTINUE
            MML=LLL
            GO TO 240
  801       IF(LLL.GT.NEX1) GO TO 803
            IPR=III(K)-1+LLL*(NEX2+NEX3)-((NEX2+NEX3)*(NEX2+NEX3+1))/2
            DO 804 I=1,NEX2
            IPR=IPR+1
            STR(I)=STR(I)+PRI(IPR)
  804       CONTINUE
            LN3=LLL-NEX3+1
            DO 807 I=LN3,LLL
            IPR=IPR+1
            STR(I)=STR(I)+PRI(IPR)
  807       CONTINUE
            MML=LLL
            GO TO 240
  803       IF(NEX2.EQ.0) GO TO 250
            IPR=III(K)-1+NEX1*NEX3-((NEX2+NEX3)*(NEX2+NEX3-1))/2+
     1      (LLL-1)*NEX2
            DO 813 I=1,NEX2
            IPR=IPR+1
            STR(I)=STR(I)+PRI(IPR)
  813       CONTINUE
            MML=NEX2
         ELSE IF (ITY(K).EQ.7) THEN
*           THERMAL TRANSFER MATRIX.
            IF (LLL.LE.NGROIN) GO TO 250
            SELF=.TRUE.
            IPR=III(K)-1+(NGRO-LLL)*IX
            DO 5003 MM=1,LLL-NGROIN
            STR(MM)=STR(MM)+PRI(IPR+MM+NGRO-LLL)
 5003       CONTINUE
            DO 5004 MM=1,NGRO-LLL
            STR(MM+LLL)=STR(MM+LLL)+PRI(IPR+MM)
 5004       CONTINUE
            MML=LLL-NGROIN
         ELSE IF (ITY(K).EQ.8) THEN
*           RECTANGLE SLOWING-DOWN MATRIX.
            IF(LLL.LT.NEXU(K))GO TO 250
            IF(LLL.GT.NEXV(K))GO TO 250
            IPR=III(K)-1+(LLL-NEXU(K))*NEXW(K)
            LN1=LLL-NEXU(K)+2
            DO 355 I=LN1,LN1+NEXW(K)-1
            IPR=IPR+1
            STR(I)=STR(I)+PRI(IPR)
  355       CONTINUE
            MML=NEXW(K)+LLL-NEXU(K)+1
         ELSE IF (ITY(K).EQ.9) THEN
*           GREULING-GOERTZEL SLOWING DOWN MATRIX.
            NEX1=NEXU(K)
            NEX3=NEXW(K)
            IF ((LLL.LT.NEX1).OR.(LLL.GT.NEX3)) GO TO 250
            SELF=.TRUE.
            NEX2=NEXV(K)
            NEX4=NEXT(K)/3
            MML=LLL-NEX1+1
            IPR=NEX3-LLL+III(K)-1
            J=IPR+1+NEX4
            DAUX=PRI(J)
            J=IPR+1+2*NEX4
            IF(LLL.GT.NEX2) THEN
               LN1=LLL-NEX2+1
            ELSE
               LN1=1
               STR(1)=STR(1)+PRI(J)
            ENDIF
            J=IPR+LN1-1
            DO 360 MM=LN1,MML
            J=J+1
            STR(MM)=STR(MM)+REAL(PRI(J)*DAUX)
  360       CONTINUE
         ENDIF
  240    NLET=MAX0(NLET,MML)
  250    CONTINUE
      ENDIF
      IF ((.NOT.SELF).AND.(IAN.EQ.0)) THEN
         STR(1)=STR(1)+SIGS(LLL)
      ELSE IF ((.NOT.SELF).AND.(IAN.EQ.1)) THEN
         STR(1)=STR(1)+3.0*SIG1(LLL)
      ENDIF
      IF (IANNN.EQ.-1) THEN
         STR(1)=STR(1)-SIG1(LLL)
      ENDIF
      RETURN
      END
