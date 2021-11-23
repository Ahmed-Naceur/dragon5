*DECK TRAGRO
      SUBROUTINE TRAGRO(IPMAC1,IPMAC2,NG,NMIL,NL,NF,NDEL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transposition of information in GROUP list.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC1  pointer to the transposed macrolib.
* IPMAC2  pointer to the original macrolib.
* NG      number of energy groups.
* NMIL    number of homogenized mixtures.
* NL      number of Legendre orders.
* NF      number of fissile isotopes.
* NDEL    number of precursor groups.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1,IPMAC2
      INTEGER NG,NMIL,NL,NF,NDEL
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC1,KPMAC1,JPMAC2,KPMAC2
      CHARACTER CM*2,NAMT1*12,NAMT2*12,TEXT12*12
      DOUBLE PRECISION SUMA,SUMB
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR3,XIOF,VOLMER
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: GAR1A,GAR1B,GAR2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(NMIL),NJJ(NMIL),IPOS(NMIL))
      ALLOCATE(GAR1A(NMIL,NF,NG),GAR1B(NMIL,NF,NG),GAR2(NG,NG,NMIL),
     1 GAR3(NMIL*NG))
*----
*  PROCESS TRANSFERT SCATTERING INFORMATION
*----
      JPMAC2=LCMGID(IPMAC2,'GROUP')
      JPMAC1=LCMLID(IPMAC1,'GROUP',NG)
      DO IL=1,NL
        WRITE(CM,'(I2.2)') IL-1
        CALL XDRSET(GAR2,NG*NG*NMIL,0.0)
        DO IG=1,NG
          KPMAC2=LCMGIL(JPMAC2,IG)
          LENGTH=1
          IF(IL.GT.1) CALL LCMLEN(KPMAC2,'SCAT'//CM,LENGTH,ITYLCM)
          IF(LENGTH.GT.0) THEN
            CALL LCMGET(KPMAC2,'SCAT'//CM,GAR3)
            CALL LCMGET(KPMAC2,'NJJS'//CM,NJJ)
            CALL LCMGET(KPMAC2,'IJJS'//CM,IJJ)
            CALL LCMGET(KPMAC2,'IPOS'//CM,IPOS)
            DO IMIL=1,NMIL
              IPOSDE=IPOS(IMIL)
              DO JG=IJJ(IMIL),IJJ(IMIL)-NJJ(IMIL)+1,-1
                GAR2(IG,JG,IMIL)=GAR3(IPOSDE) ! IG <-- JG
                IPOSDE=IPOSDE+1
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        DO IG=1,NG
          KPMAC1=LCMDIL(JPMAC1,IG)
          IPOSDE=0
          DO IMIL=1,NMIL
            IPOS(IMIL)=IPOSDE+1
            IGMIN=IG
            IGMAX=IG
            DO JG=1,NG
              IF(GAR2(NG-JG+1,NG-IG+1,IMIL).NE.0.0) THEN
                IGMIN=MIN(IGMIN,JG)
                IGMAX=MAX(IGMAX,JG)
              ENDIF
            ENDDO
            IJJ(IMIL)=IGMAX
            NJJ(IMIL)=IGMAX-IGMIN+1
            DO JG=IGMAX,IGMIN,-1
              IPOSDE=IPOSDE+1
              GAR3(IPOSDE)=GAR2(NG-JG+1,NG-IG+1,IMIL)
            ENDDO
          ENDDO
          IF(IPOSDE.GT.0) THEN
            CALL LCMPUT(KPMAC1,'SCAT'//CM,IPOSDE,2,GAR3)
            CALL LCMPUT(KPMAC1,'NJJS'//CM,NMIL,1,NJJ)
            CALL LCMPUT(KPMAC1,'IJJS'//CM,NMIL,1,IJJ)
            CALL LCMPUT(KPMAC1,'IPOS'//CM,NMIL,1,IPOS)
          ENDIF
        ENDDO
      ENDDO
*----
*  PROCESS VECTOR CROSS SECTIONS
*----
      ALLOCATE(VOLMER(NMIL))
      CALL LCMLEN(IPMAC2,'VOLUME',ILONG,ITYLCM)
      IF(ILONG.EQ.NMIL) CALL LCMGET(IPMAC2,'VOLUME',VOLMER)
      DO IG=1,NG
        KPMAC2=LCMGIL(JPMAC2,IG)
        KPMAC1=LCMDIL(JPMAC1,NG-IG+1)
        NAMT2=' '
        CALL LCMNXT(KPMAC2,NAMT2)
        TEXT12=NAMT2
   10   IF(NAMT2(:4).EQ.'SCAT') GO TO 20
        IF(NAMT2(:4).EQ.'NJJS') GO TO 20
        IF(NAMT2(:4).EQ.'IJJS') GO TO 20
        IF(NAMT2(:4).EQ.'IPOS') GO TO 20
        CALL LCMLEN(KPMAC2,NAMT2,ILONG,ITYLCM)
        IF(ILONG.GT.0) THEN
          ALLOCATE(XIOF(ILONG))
          CALL LCMGET(KPMAC2,NAMT2,XIOF)
          NAMT1=NAMT2
          IF(NAMT2(:3).EQ.'CHI') THEN
            NAMT1='NUSIGF'//NAMT2(4:)
          ELSE IF(NAMT2(:6).EQ.'NUSIGF') THEN
            NAMT1='CHI'//NAMT2(7:)
          ELSE IF(NAMT2(:9).EQ.'FLUX-INTG') THEN
            NAMT1='NWAT0'
            DO IMIL=1,NMIL
              XIOF(IMIL)=XIOF(IMIL)/VOLMER(IMIL)
            ENDDO
          ELSE IF(NAMT2.EQ.'FLUX-INTG-P1') THEN
            NAMT1='NWAT1'
            DO IMIL=1,NMIL
              XIOF(IMIL)=XIOF(IMIL)/VOLMER(IMIL)
            ENDDO
          ELSE IF(NAMT2(:5).EQ.'NWAT0') THEN
            NAMT1='FLUX-INTG'
            DO IMIL=1,NMIL
              XIOF(IMIL)=XIOF(IMIL)*VOLMER(IMIL)
            ENDDO
          ELSE IF(NAMT2(:5).EQ.'NWAT1') THEN
            NAMT1='FLUX-INTG-P1'
            DO IMIL=1,NMIL
              XIOF(IMIL)=XIOF(IMIL)*VOLMER(IMIL)
            ENDDO
          ENDIF
          CALL LCMPUT(KPMAC1,NAMT1,ILONG,2,XIOF)
          DEALLOCATE(XIOF)
        ENDIF
   20   CALL LCMNXT(KPMAC2,NAMT2)
        IF(NAMT2.NE.TEXT12) GO TO 10
      ENDDO
      DEALLOCATE(VOLMER)
*----
*  FISSION SPECTRUM NORMALIZATION
*----
      DO IDEL=0,NDEL
        IF(IDEL.EQ.0) THEN
          CM=' '
        ELSE
          WRITE(CM,'(I2.2)') IDEL
        ENDIF
        CALL XDRSET(GAR1A,NMIL*NF*NG,0.0)
        CALL XDRSET(GAR1B,NMIL*NF*NG,0.0)
        DO IG=1,NG
          KPMAC1=LCMGIL(JPMAC1,IG)
          CALL LCMGET(KPMAC1,'CHI'//CM,GAR1A(1,1,IG))
          CALL LCMGET(KPMAC1,'NUSIGF'//CM,GAR1B(1,1,IG))
        ENDDO
        DO IFIS=1,NF
          DO IMIL=1,NMIL
            SUMA=0.0D0
            SUMB=0.0D0
            DO IG=1,NG
              SUMA=SUMA+GAR1A(IMIL,IFIS,IG)
              SUMB=SUMB+GAR1B(IMIL,IFIS,IG)
            ENDDO
            IF(SUMA.GT.0.0) THEN
              DO IG=1,NG
                GAR1A(IMIL,IFIS,IG)=GAR1A(IMIL,IFIS,IG)*REAL(SUMB/SUMA)
                GAR1B(IMIL,IFIS,IG)=GAR1B(IMIL,IFIS,IG)*REAL(SUMA/SUMB)
              ENDDO
            ELSE
              DO IG=1,NG
                GAR1A(IMIL,IFIS,IG)=0.0
                GAR1B(IMIL,IFIS,IG)=0.0
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        DO IG=1,NG
          KPMAC1=LCMGIL(JPMAC1,IG)
          CALL LCMPUT(KPMAC1,'CHI'//CM,NMIL*NF,2,GAR1A(1,1,IG))
          CALL LCMPUT(KPMAC1,'NUSIGF'//CM,NMIL*NF,2,GAR1B(1,1,IG))
        ENDDO
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR3,GAR2,GAR1B,GAR1A)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
