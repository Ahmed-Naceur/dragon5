*DECK MCGSCS
      SUBROUTINE MCGSCS(KPN,K,NREG,M,NANI,NFUNL,NPJJM,KEYFLX,KEYANI,
     1                  PJJIND,NZON,XSW,PJJ,AR,PSI,MATRIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve the SCR anisotropic system.
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
* KPN     total number of unknowns in vectors SUNKNO and FUNKNO.
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes.
* M       number of material mixtures.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NFUNL   number of moments of the flux (in 2D : NFUNL=NANI*(NANI+1)/2).
* NPJJM   second dimension of PJJ.
* KEYFLX  position of flux elements in FI vector.
* KEYANI  'mode to l' index l=KEYANI(nu).
* PJJIND  index for pjj(nu <- nu') modes.
* NZON    index-number of the mixture type assigned to each volume.
* XSW     macroscopic scattering cross section.
* PJJ     used in scr acceleration.
* AR      residuals of the current iteration.
*
*Parameters: output
* PSI     corrective flux.
*
*Parameters: scratch
* MATRIX  undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER KPN,K,NREG,M,NANI,NFUNL,NPJJM,KEYFLX(NREG,NPJJM),
     1 KEYANI(NFUNL),PJJIND(NPJJM,2),NZON(K)
      REAL XSW(0:M,NANI),PJJ(NREG,NPJJM),MATRIX(NFUNL,NFUNL+1,NREG)
      DOUBLE PRECISION AR(KPN),PSI(KPN)
*---
*  LOCAL VARIABLES
*---
      INTEGER IRHS,I,INU,IMOD,INUP,L,LP,L1,LP1,NZI,IND,INDP,IER
      REAL DL,DLP,XSWI,XSWIP,TEMP
*---
* CONSTRUCT LINEAR SYSTEM BY REGION TO SOLVE
*---
      IRHS=NFUNL+1
      CALL XDRSET(MATRIX,(NFUNL*(NFUNL+1)*NREG),0.0)
      DO 10 IMOD=1,NPJJM
         INU=PJJIND(IMOD,1)
         INUP=PJJIND(IMOD,2)
         IF((INU.GT.NFUNL).OR.(INUP.GT.NFUNL)) GOTO 10
         L=KEYANI(INU)
         L1=L+1
         LP=KEYANI(INUP)
         LP1=LP+1
         DL=REAL(2*L+1)
         DLP=REAL(2*LP+1)
         DO I=1,NREG
            NZI=NZON(I)
            IND=KEYFLX(I,INU)
            INDP=KEYFLX(I,INUP)
            XSWI=XSW(NZI,L1)
            XSWIP=XSW(NZI,LP1)
            TEMP=PJJ(I,IMOD)
            MATRIX(INU,IRHS,I)=MATRIX(INU,IRHS,I)+TEMP*REAL(AR(INDP))
            MATRIX(INU,INUP,I)=-DLP*XSWIP*TEMP
            IF(INU.EQ.INUP) THEN
              MATRIX(INU,INUP,I)=MATRIX(INU,INUP,I)+1.0
            ELSE
              MATRIX(INUP,IRHS,I)=MATRIX(INUP,IRHS,I)+TEMP*REAL(AR(IND))
              MATRIX(INUP,INU,I)=-DL*XSWI*TEMP
            ENDIF
         ENDDO
 10   CONTINUE
*---
* SOLVE LINEAR SYSTEM BY REGION
*---
      DO I=1,NREG
         CALL ALSB(NFUNL,1,MATRIX(1,1,I),IER,NFUNL)
         IF(IER.NE.0) CALL XABORT('MCGSCS: PROBLEM')
         DO INU=1,NFUNL
            IND=KEYFLX(I,INU)
            PSI(IND)=MATRIX(INU,IRHS,I)
         ENDDO
      ENDDO
*
      RETURN
      END
