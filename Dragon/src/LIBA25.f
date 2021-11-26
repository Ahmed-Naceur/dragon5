*DECK LIBA25
      SUBROUTINE LIBA25(IPLIB,LABS,LDIF,LFIS,L104,NGRO,FGHOMO,NGHOMO,
     1 NSEQHO,NL,SEQHOM,SN,SB,DELTA,ISONAM,TAUX,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute self-shielded flux and cross sections from effective rates.
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
* LABS    absorption reaction flag (=.true. if the absorption reaction
*         is self-shielded).
* LDIF    scattering reaction flag (=.true. if the scattering reaction
*         is self-shielded).
* LFIS    fission reaction flag (=.true. if the fission reaction is
*         self-shielded).
* L104    resonance flux flag (=.true. if the apolib contains dilution
*         /temperature-dependent flux information). If this information
*         is not provided, it will be reconstructed from a balance
*         relation.
* NGRO    number of energy groups.
* FGHOMO  first self-shielded energy group.
* NGHOMO  number of self-shielded energy groups.
* NSEQHO  number of tabulated dilutions.
* NL      number of legendre orders required in the calculation
*         NL=1 or higher.
* SEQHOM  tabulated dilutions.
* SN      dilution of isotope.
* SB      dilution of isotope used in Livolant-Jeanpierre normalization.
* DELTA   lethargy widths.
* ISONAM  alias name of current isotope.
* TAUX    interpolated effective rates:
*         TAUX(I,1) absorption effective rates;
*         TAUX(I,2) diffusion effective rates;
*         TAUX(I,3) nu*fission effective rates;
*         TAUX(I,4) pseudo-absorption effective rates used to
*                   reconstruct the self-shielded flux.
* IMPX    print flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER NGRO,FGHOMO,NGHOMO,NSEQHO,NL,ISONAM(3),IMPX
      REAL SEQHOM(NSEQHO),SN(NGRO),SB(NGRO),DELTA(NGRO),TAUX(NGHOMO,6)
      LOGICAL LABS,LDIF,LFIS,L104
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      CHARACTER HSMG*131,STAR*1
      LOGICAL LFIS2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPRO
      DOUBLE PRECISION TMP,TMP1,AUX,FAC,ZNPHI
      REAL, ALLOCATABLE, DIMENSION(:) :: SECT,VECT,SIGF,SS1
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGS
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ITYPRO(NGRO))
      ALLOCATE(SECT(NGRO),SIGS(NGRO,NL),SCAT(NGRO,NGRO,NL),VECT(NGRO),
     1 SIGF(NGRO),SS1(NGRO))
*----
*  RECOVER INFINITE-DILUTION X-S INFORMATION.
*----
      CALL LCMGET(IPLIB,'NTOT0',SECT)
      CALL XDRLGS(IPLIB,-1,IMPX,0,NL-1,1,NGRO,SIGS,SCAT,ITYPRO)
      CALL LCMLEN(IPLIB,'NUSIGF',LENGT,ITYLCM)
      IF(LENGT.EQ.NGRO) THEN
         CALL LCMGET(IPLIB,'NUSIGF',SIGF)
         LFIS2=LFIS
      ELSE
         CALL XDRSET(SIGF,NGRO,0.0)
         LFIS2=.FALSE.
      ENDIF
*----
*  COMPUTE P0 TRANSFER PROBABILITIES.
*----
      DO 15 IG2=1,NGRO
      VECT(IG2)=SIGS(IG2,1)
      SECT(IG2)=SECT(IG2)-SIGS(IG2,1)
      DO 10 IG1=1,NGRO
      SCAT(IG2,IG1,1)=SCAT(IG2,IG1,1)/SIGS(IG1,1)
   10 CONTINUE
   15 CONTINUE
*----
*  RECOVER THE EFFECTIVE FLUX.
*----
      IF(IMPX.GT.4) WRITE(IOUT,200)
      DO 20 L=1,NGRO
      SS1(L)=1.0
   20 CONTINUE
      DO 50 L=FGHOMO,FGHOMO+NGHOMO-1
      SEIM=MAX(0.0,SN(L))
      IF(SEIM.EQ.0.) CALL XABORT('LIBA25: SELF SHIELDING FAILURE.')
      IF(LABS.AND.LDIF) THEN
*        COMPUTE THE EFFECTIVE FLUX.
         TMP1=0.0D0
         DO 30 IG2=1,FGHOMO-1
         TMP1=TMP1+SCAT(L,IG2,1)*SIGS(IG2,1)*DELTA(IG2)/DELTA(L)
   30    CONTINUE
         IF(TMP1.GT.5.0E-3*TAUX(L-FGHOMO+1,2)/DELTA(L)) THEN
*           USE A SIMPLIFIED MODEL.
            AUX=TAUX(L-FGHOMO+1,1)/DELTA(L)
         ELSE
*           USE A SLOWING-DOWN BALANCE EQUATION.
            TMP=TMP1
            DO 40 IG2=FGHOMO,FGHOMO+NGHOMO-1
            TMP=TMP+SCAT(L,IG2,1)*TAUX(IG2-FGHOMO+1,2)/DELTA(L)
   40       CONTINUE
            AUX=(TAUX(L-FGHOMO+1,1)+TAUX(L-FGHOMO+1,2))/DELTA(L)-TMP
         ENDIF
      ELSE IF(LABS) THEN
*        COMPUTE THE EFFECTIVE FLUX USING A SIMPLIFIED MODEL.
         AUX=TAUX(L-FGHOMO+1,1)/DELTA(L)
      ELSE
         AUX=0.0
      ENDIF
*
      IF((NSEQHO.EQ.1).OR.(SB(L).GE.1.0E10)) THEN
*        USE AN INFINITE DILUTION VALUE.
         ZNPHI=0.0
      ELSE IF(L104.AND.(SEIM.GT.SEQHOM(NSEQHO-1))) THEN
*        USE AN INTERPOLATED VALUE NEAR INFINITE DILUTION.
         FAC=(SEQHOM(NSEQHO-1)/SEIM)**2
         ZNPHI=FAC*TAUX(L-FGHOMO+1,4)+(1.0D0-FAC)*AUX
      ELSE IF(L104) THEN
*        USE AN INTERPOLATED VALUE.
         ZNPHI=TAUX(L-FGHOMO+1,4)
      ELSE
*        USE A CALCULATED VALUE.
         ZNPHI=AUX
      ENDIF
      PHI0=REAL(1.0D0-ZNPHI/SB(L))
      STAR=' '
      IF((PHI0.LE.0.0).OR.(PHI0.GT.1.2)) THEN
         STAR='*'
         IF(IMPX.GT.4) THEN
            WRITE(HSMG,220) PHI0,L,ZNPHI,SEIM,(ISONAM(I0),I0=1,3),STAR
            WRITE(IOUT,'(/1X,A131)') HSMG
         ENDIF
      ENDIF
      IF(PHI0.LE.0.0) PHI0=TAUX(L-FGHOMO+1,2)/DELTA(L)/VECT(L)
      SS1(L)=PHI0
      SECT(L)=TAUX(L-FGHOMO+1,1)/DELTA(L)/PHI0
      IF(LFIS2) SIGF(L)=TAUX(L-FGHOMO+1,3)/DELTA(L)/PHI0
      IF(LDIF) SIGS(L,1)=TAUX(L-FGHOMO+1,2)/DELTA(L)/PHI0
      IF(IMPX.GT.4) WRITE(IOUT,210) L,PHI0,SIGF(L),SECT(L),SIGS(L,1),
     1 SEIM,SB(L),ZNPHI
   50 CONTINUE
      IF(IMPX.GT.4) WRITE(IOUT,'(/)')
*
      CALL LCMPUT(IPLIB,'NWT0',NGRO,2,SS1)
*----
*  SELF-SHIELDING OF THE TRANSFERT CROSS SECTIONS.
*----
      IF(LDIF) THEN
         DO 65 IG1=1,NGRO
         DO 60 IG2=1,NGRO
         SCAT(IG2,IG1,1)=SCAT(IG2,IG1,1)*SIGS(IG1,1)
   60    CONTINUE
   65    CONTINUE
         INGRO=NL-1
         DO 70 IL=NL-1,0,-1
         IF(ITYPRO(IL+1).EQ.0) THEN
            INGRO=INGRO-1
         ELSE
            GO TO 80
         ENDIF
   70    CONTINUE
   80    DO 100 IL=1,NL-1
         IF(ITYPRO(IL+1).GT.0) THEN
            DO 95 IG2=1,NGRO
            SIGS(IG2,IL+1)=SIGS(IG2,IL+1)*SIGS(IG2,1)/VECT(IG2)
            DO 90 IG1=1,NGRO
            SCAT(IG2,IG1,IL+1)=SCAT(IG2,IG1,IL+1)*SIGS(IG1,1)/VECT(IG1)
   90       CONTINUE
   95       CONTINUE
         ENDIF
  100    CONTINUE
*----
*  SAVE SELF-SHIELDED X-S INFORMATION.
*----
         CALL XDRLGS(IPLIB,1,IMPX,0,INGRO,1,NGRO,SIGS,SCAT,ITYPRO)
      ENDIF
*----
*  SELF-SHIELDING OF THE FISSION CROSS SECTIONS.
*----
      IF(LFIS2) THEN
         CALL LCMGET(IPLIB,'NUSIGF',SS1)
         CALL LCMGET(IPLIB,'NFTOT',VECT)
         DO 110 I=1,NGRO
         IF(SS1(I).NE.0.0) VECT(I)=VECT(I)*SIGF(I)/SS1(I)
  110    CONTINUE
         CALL LCMPUT(IPLIB,'NFTOT',NGRO,2,VECT)
         CALL LCMPUT(IPLIB,'NUSIGF',NGRO,2,SIGF)
      ENDIF
*
      DO 120 I=1,NGRO
      SECT(I)=SECT(I)+SIGS(I,1)
  120 CONTINUE
      CALL LCMPUT(IPLIB,'NTOT0',NGRO,2,SECT)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SS1,SIGF,VECT,SCAT,SIGS,SECT)
      DEALLOCATE(ITYPRO)
      RETURN
*
  200 FORMAT(/5X,'GROUP',11X,'PHI0',10X,'SIGF0',10X,'SIGA0',10X,
     1 'SIGS0',10X,'DILUT',13X,'SB',10X,'ZNPHI')
  210 FORMAT(5X,I5,1P,7E15.5)
  220 FORMAT(47HLIBA25: *** WARNING *** INVALID VALUE OF PHI0 (,1P,
     1 E11.3,0P,10H) IN GROUP,I4,8H. ZNPHI=,1P E11.3,2X,5HSEIM=,E11.3,
     2 2X,5HISO=',3A4,2H' ,A1)
      END
