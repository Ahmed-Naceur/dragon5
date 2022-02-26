*DECK MCTOUT
      SUBROUTINE MCTOUT(IPOUT,NL,NFM,NDEL,NED,NAMEAD,NBSCO,NMERGE,
     1 NGCOND,ASCORE1,ASCORE2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Reconstruct the macrolib using a macrolib-related tally.
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
* IPOUT   pointer to the MC data structure.
* NL      number of Legendre orders required in the estimations
*         (NL=1 or higher).
* NFM     number of fissile isotopes.
* NDEL    number of delayed precursor groups.
* NED     number of extra edit vectors.
* NAMEAD  names of these extra edits.
* NBSCO   number of macrolib-related scores.
* NMERGE  number of homogenized regions.
* NGCOND  number of condensed energy groups.
* ASCORE1 score for total flux and effective multiplication factor.
* ASCORE2 macrolib score matrix.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPOUT
      INTEGER NL,NFM,NDEL,NED,NAMEAD(2,NED),NBSCO,NMERGE,NGCOND
      DOUBLE PRECISION ASCORE1(3),ASCORE2(NBSCO,NMERGE,NGCOND)
*----
*  LOCAL VARIABLES
*----
      INTEGER NSTATE
      PARAMETER(NSTATE=40)
      TYPE(C_PTR) JPOUT,KPOUT
      INTEGER ISTATE(NSTATE),IGR,JGR,IMIX,IPOSDE,IGMIN,IGMAX,IFM,IOF,
     >        IOF2,IED,I0,IL,IDEL,HSIGN(3)
      CHARACTER TEXT12*12,CM*2
      DOUBLE PRECISION SUM,SUM2
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR1,GAR2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAR3
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(NMERGE),NJJ(NMERGE),IPOS(NMERGE))
      ALLOCATE(GAR1(NMERGE),GAR2(NMERGE*NGCOND),GAR3(NMERGE,NFM))
*----
*  SCATTERING MATRIX NORMALIZATION
*----
      DO IMIX=1,NMERGE
        DO IGR=1,NGCOND
          SUM=0.0D0
          DO JGR=1,NGCOND
            SUM=SUM+ASCORE2(5+JGR,IMIX,IGR) ! JGR <-- IGR
          ENDDO
          DO JGR=1,NGCOND
            ASCORE2(5+JGR,IMIX,IGR)=ASCORE2(5+JGR,IMIX,IGR)*
     1      ASCORE2(3,IMIX,IGR)/SUM
          ENDDO
        ENDDO
      ENDDO
*----
*  FISSION SPECTRUM NORMALIZATION
*----
      IOF=5+NGCOND*NL
      IOF2=IOF+NFM
      DO IMIX=1,NMERGE
        DO IFM=1,NFM
          SUM=0.0D0
          DO IGR=1,NGCOND
            SUM=SUM+ASCORE2(IOF+IFM,IMIX,IGR)
          ENDDO
          SUM2=0.0D0
          DO IGR=1,NGCOND
            SUM2=SUM2+ASCORE2(IOF2+IFM,IMIX,IGR)
          ENDDO
          IF(SUM2.NE.0.0) THEN
            DO IGR=1,NGCOND
              ASCORE2(IOF2+IFM,IMIX,IGR)=ASCORE2(IOF2+IFM,IMIX,IGR)*
     1        SUM/SUM2
            ENDDO
          ENDIF
        ENDDO
      ENDDO
*----
*  REFORMAT CROSS-SECTION INFORMATION INTO MACROLIB FORMAT
*----
      CALL LCMSIX(IPOUT,'MACROLIB',1)
      JPOUT=LCMLID(IPOUT,'GROUP',NGCOND)
      DO IGR=1,NGCOND
        KPOUT=LCMDIL(JPOUT,IGR)
        DO IMIX=1,NMERGE
          GAR1(IMIX)=REAL(ASCORE2(1,IMIX,IGR)/ASCORE1(1))
        ENDDO
        CALL LCMPUT(KPOUT,'NWT0',NMERGE,2,GAR1)
        DO IMIX=1,NMERGE
          GAR1(IMIX)=REAL(ASCORE2(2,IMIX,IGR)/ASCORE2(1,IMIX,IGR))
        ENDDO
        CALL LCMPUT(KPOUT,'NTOT0',NMERGE,2,GAR1)
        DO IMIX=1,NMERGE
          GAR1(IMIX)=REAL((ASCORE2(3,IMIX,IGR)+2.0*ASCORE2(4,IMIX,IGR)+
     1    3.0*ASCORE2(5,IMIX,IGR))/ASCORE2(1,IMIX,IGR))
        ENDDO
        CALL LCMPUT(KPOUT,'SIGS00',NMERGE,2,GAR1)
*----
*  REFORMAT SCATTERING INFORMATION
*----
        IOF=5
        DO IL=1,NL
          WRITE(CM,'(I2.2)') IL-1
          IPOSDE=0
          DO IMIX=1,NMERGE
            IPOS(IMIX)=IPOSDE+1
            IGMIN=IGR
            IGMAX=IGR
            DO JGR=NGCOND,1,-1
              IF(ASCORE2(IOF+IGR,IMIX,JGR).NE.0.0) THEN
                IGMIN=MIN(IGMIN,JGR)
                IGMAX=MAX(IGMAX,JGR)
              ENDIF
            ENDDO
            IJJ(IMIX)=IGMAX
            NJJ(IMIX)=IGMAX-IGMIN+1
            DO JGR=IGMAX,IGMIN,-1
              IPOSDE=IPOSDE+1
              GAR2(IPOSDE)=REAL(ASCORE2(IOF+IGR,IMIX,JGR)/
     1        ASCORE2(1,IMIX,JGR))
              IF(JGR.EQ.IGR) THEN
                GAR2(IPOSDE)=GAR2(IPOSDE)+REAL((2.0*ASCORE2(4,IMIX,IGR)+
     1          3.0*ASCORE2(5,IMIX,IGR))/ASCORE2(1,IMIX,IGR))
              ENDIF
            ENDDO
            GAR1(IMIX)=REAL((ASCORE2(IOF+IGR,IMIX,IGR)
     1      +2.0*ASCORE2(4,IMIX,IGR)+3.0*ASCORE2(5,IMIX,IGR))/
     2      ASCORE2(1,IMIX,IGR))
          ENDDO
          CALL LCMPUT(KPOUT,'SCAT'//CM,IPOSDE,2,GAR2)
          CALL LCMPUT(KPOUT,'NJJS'//CM,NMERGE,1,NJJ)
          CALL LCMPUT(KPOUT,'IJJS'//CM,NMERGE,1,IJJ)
          CALL LCMPUT(KPOUT,'IPOS'//CM,NMERGE,1,IPOS)
          CALL LCMPUT(KPOUT,'SIGW'//CM,NMERGE,2,GAR1)
          IF(IL.GT.1) THEN
            CALL XDRSET(GAR1,NMERGE,0.0)
            DO IMIX=1,NMERGE
              DO JGR=1,NGCOND
                GAR1(IMIX)=GAR1(IMIX)+REAL(ASCORE2(IOF+JGR,IMIX,IGR)/
     1          ASCORE2(1,IMIX,IGR))
              ENDDO
              GAR1(IMIX)=GAR1(IMIX)+REAL((2.0*ASCORE2(4,IMIX,IGR)+3.0*
     1        ASCORE2(5,IMIX,IGR))/ASCORE2(1,IMIX,IGR))
            ENDDO
            CALL LCMPUT(KPOUT,'SIGS'//CM,NMERGE,2,GAR1)
          ENDIF
          IOF=IOF+NGCOND
        ENDDO
*----
*  REFORMAT FISSION INFORMATION
*----
        DO IDEL=1,1+NDEL
          IF(IDEL.EQ.1) THEN
            TEXT12='NUSIGF'
          ELSE
            WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL-1
          ENDIF
          DO IMIX=1,NMERGE
            DO IFM=1,NFM
              GAR3(IMIX,IFM)=REAL(ASCORE2(IOF+IFM,IMIX,IGR)/
     1        ASCORE2(1,IMIX,IGR))
            ENDDO
          ENDDO
          CALL LCMPUT(KPOUT,TEXT12,NMERGE*NFM,2,GAR3)
          IF(IDEL.EQ.1) THEN
            TEXT12='CHI'
          ELSE
            WRITE(TEXT12,'(3HCHI,I2.2)') IDEL-1
          ENDIF
          IOF2=IOF+NFM
          DO IMIX=1,NMERGE
            DO IFM=1,NFM
              SUM=0.0D0
              DO JGR=1,NGCOND
                SUM=SUM+ASCORE2(IOF+IFM,IMIX,JGR)
              ENDDO
              IF(SUM.NE.0.0) THEN
                GAR3(IMIX,IFM)=REAL(ASCORE2(IOF2+IFM,IMIX,IGR)/SUM)
              ELSE
                GAR3(IMIX,IFM)=0.0
              ENDIF
            ENDDO
          ENDDO
          CALL LCMPUT(KPOUT,TEXT12,NMERGE*NFM,2,GAR3)
          IOF=IOF2+NFM
        ENDDO
*----
*  REFORMAT ADDITIONAL EDIT INFORMATION
*----
        DO IED=1,NED
          WRITE(TEXT12,'(2A4)') NAMEAD(1,IED),NAMEAD(2,IED)
          DO IMIX=1,NMERGE
            GAR1(IMIX)=REAL(ASCORE2(IOF+IED,IMIX,IGR)/
     1      ASCORE2(1,IMIX,IGR))
          ENDDO
          CALL LCMPUT(KPOUT,TEXT12,NMERGE,2,GAR1)
        ENDDO
      ENDDO
*----
*  CREATE THE STATE VECTOR AND THE SIGNATURE
*----
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NGCOND
      ISTATE(2)=NMERGE
      ISTATE(3)=NL
      ISTATE(4)=NFM
      ISTATE(5)=NED
      ISTATE(7)=NDEL
      CALL LCMPUT(IPOUT,'STATE-VECTOR',NSTATE,1,ISTATE)
      TEXT12='L_MACROLIB'
      READ(TEXT12,'(3A4)') (HSIGN(I0),I0=1,3)
      CALL LCMPUT(IPOUT,'SIGNATURE',3,3,HSIGN)
      CALL LCMSIX(IPOUT,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR3,GAR2,GAR1)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
