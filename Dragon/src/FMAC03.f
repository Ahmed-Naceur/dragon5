*DECK FMAC03
      SUBROUTINE FMAC03(IPMACR,IG,IPART,NGP,MAXLEN,NANISO,NK,NPART,
     1 HNPRT,NGPRT,NWA,H2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Save a SCAT cross section in the GROUP list of a MACROLIB.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMACR  LCM object address of the MACROLIB.
* IG      secondary energy group.
* IPART   index of the particle type corresponding to the MACROLIB.
* NGP     sum of number of energy groups for all types of particles.
* MAXLEN  second dimension of array H2.
* NANISO  maximum scattering anisotropy.
* NK      number of mixtures.
* NPART   number of particle types.
* HNPRT   character*1 names of particle types.
* NGPRT   number of energy groups per particle type.
* NWA     Legendre order of scattering cross-section information.
* H2      scattering cross-section information.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR
      INTEGER IG,IPART,NGP,MAXLEN,NANISO,NK,NPART,NGPRT(NPART),
     1 NWA(NGP,NK)
      CHARACTER(LEN=1) HNPRT(NPART)
      REAL H2(NGP,MAXLEN)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) KPMACR
      CHARACTER(LEN=2) CM
      CHARACTER(LEN=12) HGROUP
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: JPMACR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT
*----
*  DEFINE GROUP DIRECTORIES PER PARTICLE TYPE
*----
      ALLOCATE(JPMACR(NPART))
      DO JPART=1,NPART
        IF(JPART.EQ.IPART) THEN
          HGROUP='GROUP'
        ELSE
          HGROUP='GROUP-'//HNPRT(JPART)
        ENDIF
        JPMACR(JPART)=LCMLID(IPMACR,HGROUP,NGPRT(IPART))
      ENDDO
      IG1=1
      DO I=1,IPART-1
        IG1=IG1+NGPRT(I)
      ENDDO
      IG2=IG1+NGPRT(IPART)-1
      IGR=IG-IG1+1
*----
*  LOOP OVER PARTICLE TYPES
*----
      DO JPART=1,NPART
        ALLOCATE(SCAT(NK,NGPRT(JPART),NANISO+1))
        SCAT(:NK,:NGPRT(JPART),:NANISO)=0.0
        JG1=1
        DO I=1,JPART-1
          JG1=JG1+NGPRT(I)
        ENDDO
        JG2=JG1+NGPRT(JPART)-1
*----
*  LOOP OVER TRANSITIONS
*----
        DO JG=JG1,JG2
*         Loop over primary energy groups
          DO IBM=1,NK
            IF(NWA(JG,IBM).NE.0) GO TO 10
          ENDDO
          CYCLE
*         Find the primary particle type
   10     JGR=JG-JG1+1
          IOF=0
          DO IBM=1,NK
            IF(NWA(JG,IBM).GT.0) CALL XABORT('FMAC03: POSITIVE NWA NOT'
     1      //' IMPLEMENTED.')
            IF(-NWA(JG,IBM).GT.NANISO+1) CALL XABORT('FMAC03: NWA OVER'
     1      //'FLOW.')
            DO IL=1,-NWA(JG,IBM)
              SCAT(IBM,JGR,IL)=H2(JG,IOF+IL)
            ENDDO
            IOF=IOF-NWA(JG,IBM)
          ENDDO
        ENDDO
*----
*  SAVE SCATTERING INFORMATION ON MACROLIB
*----
        ALLOCATE(NJJ(NK),IJJ(NK),IPOS(NK),GAR(NK*NGPRT(JPART)))
        KPMACR=LCMDIL(JPMACR(JPART),IGR)
        DO IL=1,NANISO
          WRITE (CM,'(I2.2)') IL-1
          IPOSIT=0
          DO IBM=1,NK
            J2=IGR
            J1=IGR
            DO JGR=1,NGPRT(JPART)
              IF(SCAT(IBM,JGR,IL).NE.0.0) THEN
                J2=MAX(J2,JGR)
                J1=MIN(J1,JGR)
              ENDIF
            ENDDO
            NJJ(IBM)=J2-J1+1
            IJJ(IBM)=J2
            IPOS(IBM)=IPOSIT+1
            DO JGR=J2,J1,-1
              IPOSIT=IPOSIT+1
              IF(IPOSIT.GT.NK*NGPRT(JPART)) CALL XABORT('bug')
              GAR(IPOSIT)=SCAT(IBM,JGR,IL)
            ENDDO
          ENDDO
          CALL LCMPUT(KPMACR,'SIGW'//CM,NK,2,SCAT(1,IGR,IL))
          CALL LCMPUT(KPMACR,'SCAT'//CM,IPOSIT,2,GAR)
          CALL LCMPUT(KPMACR,'NJJS'//CM,NK,1,NJJ)
          CALL LCMPUT(KPMACR,'IJJS'//CM,NK,1,IJJ)
          CALL LCMPUT(KPMACR,'IPOS'//CM,NK,1,IPOS)
        ENDDO
        DEALLOCATE(GAR,IPOS,IJJ,NJJ,SCAT)
      ENDDO
      DEALLOCATE(JPMACR)
      RETURN
      END
