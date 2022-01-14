*DECK CPOMAR
      SUBROUTINE CPOMAR(IPEDIT,NGROUP,NMERGE,NL    ,NIFISS,NEDMAC,
     >                  HVECT ,IVECT ,NPROC ,ILEAKS,DXSMAC,DSCMAC,
     >                  EMJMAC,DISFC ,IFCDIS,DISFAC )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Get macroscopic cross section from IPEDIT.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input
* IPEDIT  pointer to the edit.
* NGROUP  number of groups condensed.
* NMERGE  number of regions merged.
* NL      number of Legendre orders.
* NIFISS  number of fissile isotopes.
* NEDMAC  number of extra edit vectors.
* HVECT   name of additional xs.
* IVECT   location of additional xs.
* NPROC   number of microscopic xs to process.
* ILEAKS  leak calculation:
*         = 0 no leakage;
*         = 1 homogeneous leakage coefficients;
*         = 2 directional leakage coefficients.
*
*Parameters: output
* DXSMAC  averaged region/group x-s.
* DSCMAC  scattering rates.
* DISFC   disadvantage factor.
* EMJMAC  energy per fission.
* IFCDIS  discontinuity factor present (1) or absent.
* DISFAC  discontinuity factors.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPEDIT
      INTEGER          NGROUP,NMERGE,NL,NIFISS,NEDMAC,IVECT(NEDMAC),
     >                 NPROC,ILEAKS
      CHARACTER        HVECT(NEDMAC)*8
      REAL             DISFC(NGROUP),
     >                 EMJMAC(NMERGE)
      DOUBLE PRECISION DXSMAC(NGROUP,NPROC,NMERGE),
     >                 DSCMAC(NGROUP,NGROUP,NL,NMERGE)
      INTEGER          IFCDIS
      DOUBLE PRECISION DISFAC(2,NGROUP,3)
*----
*  LOCAL PARAMETERS
*----
      TYPE(C_PTR)      JPEDIT,KPEDIT
      INTEGER          NDPROC
      PARAMETER       (NDPROC=20)
      INTEGER          IGR,IED,IXSR,JXSR,KXSR,IMRG,IFIS,ILOCED,
     >                 IL,IPOSIT,JGR1,JGR2,JGR,ILCMLN,ITYLCM
      CHARACTER        CM*2
      INTEGER          IDIR,IPL,IEL
      REAL             TEMP(6)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ
      REAL, ALLOCATABLE, DIMENSION(:) :: SCATC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DNUFI
*----
*  SCRATCH STORAGE ALLOCATION
*   SCATC   compress scattering data.
*   IJJ     position of first diffusion group.
*   NJJ     number of diffusion group.
*   DNUFI   fission source.
*----
      ALLOCATE(IJJ(NMERGE),NJJ(NMERGE))
      ALLOCATE(SCATC(NMERGE*NGROUP))
      ALLOCATE(DNUFI(NMERGE,NIFISS+1))
*----
*  INITIALIZE REACTION RATE VECTOR
*----
      CALL XDDSET(DXSMAC,NGROUP*NPROC*NMERGE,0.0D0)
      CALL XDDSET(DSCMAC,NGROUP*NGROUP*NL*NMERGE,0.0D0)
      CALL XDDSET(DNUFI,NMERGE*(NIFISS+1),0.0D0)
*----
*  READ ALL CROSS SECTION FROM IPEDIT EXCEPT CHI
*----
      CALL LCMLEN(IPEDIT,'FLUXDISAFACT',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.NGROUP) THEN
        CALL LCMGET(IPEDIT,'FLUXDISAFACT',DISFC)
      ELSE
        CALL XDRSET(DISFC,NGROUP,0.0)
      ENDIF
      JPEDIT=LCMGID(IPEDIT,'GROUP')
      DO 100 IGR=1,NGROUP
        KPEDIT=LCMGIL(JPEDIT,IGR)
        IF(NEDMAC.GT.0) THEN
          DO 110 IED=1,NEDMAC
            IXSR=IVECT(IED)
            IF(IXSR.GT.0) THEN
              CALL LCMGET(KPEDIT,HVECT(IED),SCATC)
              DO 111 IMRG=1,NMERGE
                DXSMAC(IGR,IXSR,IMRG)=DBLE(SCATC(IMRG))
 111          CONTINUE
            ENDIF
 110      CONTINUE
        ENDIF
        IXSR=NDPROC+NL+1
        CALL LCMLEN(KPEDIT,'OVERV',ILCMLN,ITYLCM)
        IF(ILCMLN.EQ.NMERGE) THEN
          CALL LCMGET(KPEDIT,'OVERV',SCATC)
          DO 120 IMRG=1,NMERGE
            DXSMAC(IGR,IXSR,IMRG)=DBLE(SCATC(IMRG))
 120      CONTINUE
        ENDIF
        IXSR=1
        CALL LCMGET(KPEDIT,'NTOT0',SCATC)
        DO 130 IMRG=1,NMERGE
          DXSMAC(IGR,IXSR,IMRG)=DBLE(SCATC(IMRG))
 130    CONTINUE
        IXSR=2
        CALL LCMGET(KPEDIT,'TRANC',SCATC)
        DO 170 IMRG=1,NMERGE
          DXSMAC(IGR,IXSR,IMRG)=DBLE(SCATC(IMRG))
 170    CONTINUE
        IF(IFCDIS .EQ. 1) THEN
          CALL LCMLEN(KPEDIT,'ADFGENERAL',ILCMLN,ITYLCM)
          IF(ILCMLN .EQ. 6) THEN
            CALL LCMGET(KPEDIT,'ADFGENERAL',TEMP)
            IEL=0
            DO  IDIR=1,3
              DO  IPL=1,2
                IEL=IEL+1
                DISFAC(IPL,IGR,IDIR)=DBLE(TEMP(IEL))
              ENDDO
            ENDDO
          ELSE
            IFCDIS=0
          ENDIF
        ENDIF
        IXSR=16
        CALL LCMGET(KPEDIT,'FLUX-INTG',SCATC)
        DO 190 IMRG=1,NMERGE
          DXSMAC(IGR,IXSR,IMRG)=DBLE(SCATC(IMRG))
 190    CONTINUE
        IF(NIFISS.GT.0) THEN
          IXSR=3
          JXSR=16
          CALL LCMGET(KPEDIT,'NUSIGF',SCATC)
          ILOCED=1
          DO 150 IFIS=1,NIFISS
            DO 151 IMRG=1,NMERGE
              DXSMAC(IGR,IXSR,IMRG)=DXSMAC(IGR,IXSR,IMRG)
     >          +DBLE(SCATC(ILOCED))
              DNUFI(IMRG,IFIS)=DNUFI(IMRG,IFIS)
     >          +DBLE(SCATC(ILOCED))*DXSMAC(IGR,JXSR,IMRG)
              ILOCED=ILOCED+1
 151        CONTINUE
 150      CONTINUE
          IXSR=4
          CALL LCMGET(KPEDIT,'NFTOT',SCATC)
          DO 153 IMRG=1,NMERGE
            DXSMAC(IGR,IXSR,IMRG)=DXSMAC(IGR,IXSR,IMRG)
     >        +DBLE(SCATC(IMRG))
 153      CONTINUE
        ENDIF
        IXSR=NDPROC
        DO 200 IL=1,NL
          IXSR=IXSR+1
          WRITE (CM,'(I2.2)') IL-1
          CALL LCMGET(KPEDIT,'SCAT'//CM,SCATC)
          CALL LCMGET(KPEDIT,'NJJS'//CM,NJJ)
          CALL LCMGET(KPEDIT,'IJJS'//CM,IJJ)
          IPOSIT=0
          DO 210 IMRG=1,NMERGE
            JGR2=IJJ(IMRG)
            JGR1=JGR2-NJJ(IMRG)+1
            DO 211 JGR=JGR2,JGR1,-1
              IPOSIT=IPOSIT+1
              DSCMAC(IGR,JGR,IL,IMRG)=DBLE(SCATC(IPOSIT))
              DXSMAC(JGR,IXSR,IMRG)=DXSMAC(JGR,IXSR,IMRG)
     >                             +DSCMAC(IGR,JGR,IL,IMRG)
 211        CONTINUE
 210      CONTINUE
 200    CONTINUE
        IF(ILEAKS.EQ.1) THEN
          IXSR=17
          CALL LCMGET(KPEDIT,'DIFF',SCATC)
          DO 180 IMRG=1,NMERGE
            IF(SCATC(IMRG).GT.0.0) THEN
              DXSMAC(IGR,IXSR,IMRG)=1.0D0/(3.0D0*DBLE(SCATC(IMRG)))
            ENDIF
 180      CONTINUE
        ELSE IF(ILEAKS.EQ.2) THEN
          IXSR=17
          CALL LCMGET(KPEDIT,'DIFF',SCATC)
          DO 181 IMRG=1,NMERGE
            IF(SCATC(IMRG).GT.0.0) THEN
              DXSMAC(IGR,IXSR,IMRG)=1.0D0/(3.0D0*DBLE(SCATC(IMRG)))
            ENDIF
 181      CONTINUE
          IXSR=18
          CALL LCMGET(KPEDIT,'DIFFX',SCATC)
          DO 182 IMRG=1,NMERGE
            IF(SCATC(IMRG).GT.0.0) THEN
              DXSMAC(IGR,IXSR,IMRG)=1.0D0/(3.0D0*DBLE(SCATC(IMRG)))
            ENDIF
 182      CONTINUE
          IXSR=19
          CALL LCMGET(KPEDIT,'DIFFY',SCATC)
          DO 183 IMRG=1,NMERGE
            IF(SCATC(IMRG).GT.0.0) THEN
              DXSMAC(IGR,IXSR,IMRG)=1.0D0/(3.0D0*DBLE(SCATC(IMRG)))
            ENDIF
 183      CONTINUE
          IXSR=20
          CALL LCMGET(KPEDIT,'DIFFZ',SCATC)
          DO 184 IMRG=1,NMERGE
            IF(SCATC(IMRG).GT.0.0) THEN
              DXSMAC(IGR,IXSR,IMRG)=1.0D0/(3.0D0*DBLE(SCATC(IMRG)))
            ENDIF
 184      CONTINUE
        ELSE
          IXSR=17
          JXSR=1
          IF(NL.GE.2) THEN
            KXSR=NDPROC+2
            DO 185 IMRG=1,NMERGE
              DXSMAC(IGR,IXSR,IMRG)=DXSMAC(IGR,JXSR,IMRG)
     >                             -DXSMAC(IGR,KXSR,IMRG)
 185        CONTINUE
          ELSE
            DO 186 IMRG=1,NMERGE
              DXSMAC(IGR,IXSR,IMRG)=DXSMAC(IGR,JXSR,IMRG)
 186        CONTINUE
          ENDIF
        ENDIF
 100  CONTINUE
*----
*  PROCESS CHI IF REQUIRED
*----
      IF(NIFISS.GT.0) THEN
        DO 160 IGR=1,NGROUP
          KPEDIT=LCMGIL(JPEDIT,IGR)
          IXSR=5
          CALL LCMGET(KPEDIT,'CHI',SCATC)
          ILOCED=1
          DO 161 IFIS=1,NIFISS
            DO 162 IMRG=1,NMERGE
              DXSMAC(IGR,IXSR,IMRG)=DXSMAC(IGR,IXSR,IMRG)
     >          +DBLE(SCATC(ILOCED))*DNUFI(IMRG,IFIS)
              ILOCED=ILOCED+1
 162        CONTINUE
 161      CONTINUE
 160    CONTINUE
      ENDIF
*----
*  FIND TOTAL ENERGY PRODUCTION
*----
      JXSR=16
      CALL XDRSET(EMJMAC,NMERGE,0.0)
      DO 251 IGR=1,NGROUP
        KPEDIT=LCMGIL(JPEDIT,IGR)
        CALL LCMLEN(KPEDIT,'H-FACTOR',ILCMLN,ITYLCM)
        IF(ILCMLN.EQ.NMERGE) THEN
          CALL LCMGET(KPEDIT,'H-FACTOR',SCATC)
          DO 250 IMRG=1,NMERGE
            EMJMAC(IMRG)=EMJMAC(IMRG)+REAL(DXSMAC(IGR,JXSR,IMRG))*
     >      SCATC(IMRG)*1.0E18
 250      CONTINUE
        ENDIF
 251  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DNUFI)
      DEALLOCATE(SCATC)
      DEALLOCATE(NJJ,IJJ)
      RETURN
      END
