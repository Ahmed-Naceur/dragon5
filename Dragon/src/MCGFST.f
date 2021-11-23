*DECK MCGFST
      SUBROUTINE MCGFST(NGEFF,KPSYS,NCONV,KPN,K,NREG,NANI,NFUNL,NPJJM,
     1           KEYFLX,KEYCUR,PJJIND,NZON,V,S,PHIOUT,IDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Addition of the contribution to the flux of the regional source
* when the 'source term isolation' option is turned on for the 
* method of characteristics integration.
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
* NGEFF   number of groups to process.
* KPSYS   pointer array for each group properties.
* NCONV   logical array of convergence status for each group (.TRUE.
*         not converged).
* KPN     total number of unknowns per group in flux vector.
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NFUNL   number of moments of the flux (in 2D : NFUNL=NANI*(NANI+1)/2).
* NPJJM   number of pjj modes to store for STIS option.
* KEYFLX  position of flux elements in flux vector.
* KEYCUR  position of current elements in flux vector.
* PJJIND  index for pjj(nu <- nu') modes.
* NZON    index-number of the mixture type assigned to each volume.
* V       volumes.
* S       source vector.
* IDIR    direction of fundamental current for TIBERE with MoC
*         (=0,1,2,3).
*
*Parameters: input/output
* PHIOUT  flux vector.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPSYS(NGEFF)
      INTEGER NGEFF,KPN,K,NREG,NANI,NFUNL,NPJJM,KEYFLX(NREG,NFUNL),
     1 KEYCUR(K-NREG),PJJIND(NPJJM,2),NZON(K),IDIR
      REAL V(K)
      DOUBLE PRECISION S(KPN,NGEFF),PHIOUT(KPN,NGEFF)
      LOGICAL NCONV(NGEFF)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPSYS
      INTEGER II,I,IND,IMOD,INU,INUP,INDP
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PJJ,PJJI
      CHARACTER*12 NPJJT,NPJJIT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PJJ(NREG,NPJJM),PJJI(NREG,NPJJM))
*
      IF(IDIR .EQ.0) THEN
        NPJJT='PJJ$MCCG'
        NPJJIT='            '
      ELSEIF(IDIR .EQ. 1) THEN
        NPJJT='PJJX$MCCG'
        NPJJIT='PJJXI$MCCG'
      ELSEIF(IDIR .EQ. 2) THEN
        NPJJT='PJJY$MCCG'
        NPJJIT='PJJYI$MCCG'
      ELSE
        NPJJT='PJJZ$MCCG'
        NPJJIT='PJJZI$MCCG'
      ENDIF
      IF(NANI.LE.0) CALL XABORT('MCGFST: INVALID VALUE OF NANI.')
      DO II=1,NGEFF
      IF (NCONV(II)) THEN
         JPSYS=KPSYS(II)
         CALL LCMGET(JPSYS,NPJJT,PJJ)
         IF(IDIR.GT.0) CALL LCMGET(JPSYS,NPJJIT,PJJI)
         DO I=1,K
         IF(V(I).GT.0.) THEN 
            IF (NZON(I).LT.0) THEN
               IND=KEYCUR(I-NREG)
               PHIOUT(IND,II)=PHIOUT(IND,II)/V(I)
            ELSE
               DO INU=1,NFUNL
                  IND=KEYFLX(I,INU)
                  PHIOUT(IND,II)=PHIOUT(IND,II)/V(I)
               ENDDO
*              DIVIDE THE EXTRA TERMS XI, YI, AND ZI BY THE VOLUME
               IF(IDIR.NE.0) THEN
                 IND=KEYFLX(I,NFUNL)
                 PHIOUT(IND+KPN/2,II)=PHIOUT(IND+KPN/2,II)/V(I)
               ENDIF
               DO IMOD=1,NPJJM
                  INU=PJJIND(IMOD,1)
                  INUP=PJJIND(IMOD,2)
                  IND=KEYFLX(I,INU)
                  INDP=KEYFLX(I,INUP)
                  PHIOUT(IND,II)=PHIOUT(IND,II)+
     1                           PJJ(I,IMOD)*S(INDP,II)
                  IF(IDIR .GT. 0) THEN
                     PHIOUT(IND+KPN/2,II)=PHIOUT(IND+KPN/2,II)+
     1                           PJJI(I,IMOD)*S(INDP,II)
                  ENDIF
                  IF (INU.NE.INUP) THEN
                     PHIOUT(INDP,II)=PHIOUT(INDP,II)+
     1                           PJJ(I,IMOD)*S(IND,II)
                     IF(IDIR .GT. 0) THEN
                       PHIOUT(INDP+KPN/2,II)=PHIOUT(INDP+KPN/2,II)+
     1                              PJJI(I,IMOD)*S(IND,II)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
         ENDDO
      ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PJJI,PJJ)
*
      RETURN
      END
