*DECK MCGFCR
      SUBROUTINE MCGFCR(IPRINT,IG,II,NG,NGEFF,KPN,N,NREG,NANI,NFUNL,M,
     1                  LTYPE,KEYFLX,KEYCUR,NZON,NGINDV,FORM,REBAL,FI,
     2                  FIOLD,SC,TAB,NJJ,IJJ,IPOS,XSCAT,AR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute residual of a previous free iterations for ACA method.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IPRINT  print parameter (equal to zero for no print).
* IG      index of group to process in "NG format".
* II      index of group to process in "NGEFF format".
* NG      number of groups.
* NGEFF   number of groups to process.
* KPN     total number of unknowns in vectors SUNKNO and FUNKNO.
* N       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NFUNL   number of moments of the flux (in 2D : NFUNL=NANI*(NANI+1)/2).
* M       number of material mixtures.
* LTYPE   flag to know how the residual vector is organized:
*         .TRUE.  for ACA, with permutation array, only the isotropic
*         moments;
*         .FALSE. for SCR, without permutation array, all the moments.
* KEYFLX  position of flux elements in FI vector.
* KEYCUR  position of current elements in FI vector.
* NZON    index-number of the mixture type assigned to each volume.
* NGINDV  index to pass from "NGEFF format" to "NG format"
* FORM    input flux format flag (.TRUE. same format as output flux ;
*         .FALSE. same format as input source).
* REBAL   type of acceleration (.TRUE. rebalancing ; .FALSE.
*         inner iterations acceleration).
* FI      zonal scalar flux.
* SC      macroscopic "in group" scattering cross section.
* FIOLD   old zonal scalar flux.
* TAB     if LTYPE,  IPERM(N) permutation array,
*         otherwise, KEYANI(NFUNL) 'mode to l' index: l=KEYANI(nu).
* NJJ     scattering information.
* IJJ     scattering information.
* IPOS    scattering information.
* XSCAT   scattering information.
*
*Parameters: output
* AR      residual form previous free iteration.
*
*-----------------------------------------------------------------------
*
*----
* SUBROUTINE ARGUMENTS
*----
      INTEGER IPRINT,IG,II,NG,NGEFF,KPN,N,NREG,NANI,NFUNL,M,
     1 KEYFLX(NREG,NFUNL),KEYCUR(*),TAB(*),NZON(N),NGINDV(NG),
     2 NJJ(0:M),IJJ(0:M),IPOS(0:M)
      REAL FIOLD(KPN,NG),SC(0:M,NANI),XSCAT(0:M*NG)
      DOUBLE PRECISION AR(*),FI(KPN,NGEFF)
      LOGICAL REBAL,FORM,LTYPE
*
      IF(IPRINT.GT.99) WRITE(6,'(23H MCGFCR: PROCESS GROUPS,2I6)') IG,II
      IGG=IG
      IF (FORM) IGG=II
*
      IF (LTYPE) THEN
*---
* ACA RESIDUAL
*---
      DO I=1,N
         J=TAB(I)
         IBM=NZON(J)
         IF(IBM.GE.0) THEN
            SIGC=SC(IBM,1)
            IND=KEYFLX(J,1)
         ELSE
            SIGC=0.5
            IND=KEYCUR(J-NREG)
         ENDIF
         AR(I)=(FI(IND,II)-FIOLD(IND,IGG))*SIGC
         IF ((REBAL).AND.(IBM.GT.0)) THEN
*        rebalancing option on : contribution from other groups.
            JG=IJJ(IBM)
            DO 10 JND=1,NJJ(IBM)
               IF(JG.NE.IG) THEN
               JJ=NGINDV(JG)
               IF (JJ.GT.0) THEN
                  JGG=JG
                  IF (FORM) JGG=JJ
                  AR(I)=AR(I)+XSCAT(IPOS(IBM)+JND-1)*
     1                 (FI(IND,JJ)-FIOLD(IND,JGG))
               ENDIF
               ENDIF
               JG=JG-1
 10         CONTINUE
         ENDIF
      ENDDO
      ELSE
*---
* SCR RESIDUAL
*---
      DO I=1,N
         IBM=NZON(I)
         IF(IBM.GE.0) THEN
            DO INU=1,NFUNL
               IND=KEYFLX(I,INU)
               IL=TAB(INU)
               SIGC=REAL(2*IL+1)*SC(IBM,IL+1)
               AR(IND)=(FI(IND,II)-FIOLD(IND,IGG))*SIGC
            ENDDO
         ELSE
            SIGC=0.5
            IND=KEYCUR(I-NREG)
            AR(IND)=(FI(IND,II)-FIOLD(IND,IGG))*SIGC
         ENDIF
         IF ((REBAL).AND.(IBM.GT.0)) THEN
*        rebalancing option on : contribution from other groups.
            IND=KEYFLX(I,1)
            JG=IJJ(IBM)
            DO 20 JND=1,NJJ(IBM)
               IF(JG.NE.IG) THEN
               JJ=NGINDV(JG)
               IF (JJ.GT.0) THEN
                  JGG=JG
                  IF (FORM) JGG=JJ
                  AR(IND)=AR(IND)+XSCAT(IPOS(IBM)+JND-1)*
     1                 (FI(IND,JJ)-FIOLD(IND,JGG))
               ENDIF
               ENDIF
               JG=JG-1
 20         CONTINUE
         ENDIF
      ENDDO
      ENDIF
      RETURN
      END
