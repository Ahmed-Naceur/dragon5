*DECK FMAC01
      SUBROUTINE FMAC01(IPMACR,IMPX,HPART,LIN,IVERS,NGP,NPART,NGXI,
     1 NEDIT,NUCL,NK,NUFIS,MASM3,N,NGPRT,HNPRT,NPMIN,NPMAX,NANIS,MUFIS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover data and 1D cross sections from the FMAC-M ascii file.
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
* IMPX    print flag.
* HPART   character*1 name of the MACROLIB particle.
* LIN     unit number of the FMAC-M ascii file.
* IVERS   file version number.
* NGP     sum of number of energy groups for all types of particles.
* NPART   number of particle types.
* NGXI    number of groups with non-zero fission spectrum.
* NEDIT   number of additional edit cross sections.
* NUCL    number of nuclides.
* NK      number of mixtures.
* NUFIS   number of fission materials.
* MASM3   length of integer control array.
* N       integer control array.
*
*Parameters: output
* NGPRT   number of energy groups per particle type.
* HNPRT   character*1 names of particle types.
* NPMIN   minimum transition group number.
* NPMAX   maximum transition group number.
* NANIS   number of Legendre orders per energy group.
* MUFIS   fission material number per mixture.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR
      INTEGER IMPX,LIN,IVERS,NGP,NPART,NGXI,NEDIT,NUCL,NK,NUFIS,
     1 MASM3,N(MASM3),NGPRT(NPART),NPMIN(NGP),NPMAX(NGP),NANIS(NGP),
     2 MUFIS(NUFIS)
      CHARACTER(LEN=1) HPART,HNPRT(NPART)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      TYPE(C_PTR) JPMACR
      CHARACTER TEX(50)*6,HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: H1
      REAL, ALLOCATABLE, DIMENSION(:,:) :: H2
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: H3
      CHARACTER(LEN=6), ALLOCATABLE, DIMENSION(:) :: HNGAR
*
      NGROUP=0
      IPART=0
      DO K=1,MASM3
        LTOT=N(K)
        IF(LTOT.EQ.0) CYCLE
        IG1=0
        IG2=0
        SELECT CASE(K)
        CASE(1)
*         title
          READ(LIN,'(12A6)') (TEX(I),I=1,LTOT)
          IF(IMPX.GT.0) WRITE(IOUT,1000) (TEX(I),I=1,LTOT)
        CASE(2)
*         number of energy groups by particle type
          IF(LTOT.NE.NPART) CALL XABORT('FMAC01: BAD RECORD 2.')
          READ(LIN,'(6I12)') (NGPRT(I),I=1,LTOT)
          CALL LCMPUT(IPMACR,'PARTICLE-NGR',NPART,1,NGPRT)
        CASE(3)
*         particle names by particle type
          IF(LTOT.NE.NPART) CALL XABORT('FMAC01: BAD RECORD 3.')
          ALLOCATE(HNGAR(LTOT))
          READ(LIN,'(12A6)') (HNGAR(I),I=1,LTOT)
          DO I=1,LTOT
            IF(HNGAR(I).EQ.'NEUT') THEN
              HNPRT(I)='N'
            ELSE IF(HNGAR(I).EQ.'GAMA') THEN
              HNPRT(I)='G'
            ELSE IF(HNGAR(I).EQ.'BETA') THEN
              HNPRT(I)='B'
            ELSE IF(HNGAR(I).EQ.'POSITR') THEN
              HNPRT(I)='C'
            ELSE IF(HNGAR(I).EQ.'PROT') THEN
              HNPRT(I)='P'
            ELSE
              WRITE(HSMG,'(8HFMAC01: ,A6,26H IS AN INVALID PARTICLE NA,
     1        3HME.)') HNGAR(I)
              CALL XABORT(HSMG)
            ENDIF
          ENDDO
          DEALLOCATE(HNGAR)
          CALL LCMPTC(IPMACR,'PARTICLE-NAM',1,NPART,HNPRT)
          DO I=1,NPART
            IF(HNPRT(I).EQ.HPART) THEN
              IPART=I
              GO TO 30
            ENDIF
          ENDDO
          CALL XABORT('FMAC01: PARTICLE '//HPART//' NOT AVAILABLE IN'
     1    //' FMAC-M FILE.')
   30     CONTINUE
        CASE(4)
*         rest energies by particle type
          IF(LTOT.NE.NPART) CALL XABORT('FMAC01: BAD RECORD 4.')
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          DO I=1,LTOT
            H1(I)=H1(I)*1.0E6
          ENDDO
          CALL LCMPUT(IPMACR,'PARTICLE-MC2',NPART,2,H1)
          DEALLOCATE(H1)
        CASE(5)
*         energy mesh boundaries for all particles
          IF(LTOT.NE.NGP+NPART) CALL XABORT('FMAC01: BAD RECORD 5.')
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          IF(IPART.EQ.0) CALL XABORT('FMAC01: PARTICLE TYPE UNDEFINED.')
          DO I=1,LTOT
            H1(I)=H1(I)*1.0E6
          ENDDO
          IG1=1
          DO I=1,IPART-1
            IG1=IG1+NGPRT(I)+1
          ENDDO
          IG2=IG1+NGPRT(IPART)
          NGROUP=NGPRT(IPART)
          JPMACR=LCMLID(IPMACR,'GROUP',NGROUP)
          CALL FMAC04(NGPRT,NGP,NPART,1,H1)
          CALL LCMPUT(IPMACR,'ENERGY',NGROUP+1,2,H1)
          DEALLOCATE(H1)
        CASE(6)
*         group velocities
          IF(LTOT.NE.NGP) CALL XABORT('FMAC01: BAD RECORD 6.')
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          DEALLOCATE(H1)
        CASE(7)
*         nuclide names
          IF(LTOT.NE.NUCL) CALL XABORT('FMAC01: BAD RECORD 7.')
          READ(LIN,'(12A6)') (TEX(I),I=1,LTOT)
        CASE(8)
*         nuclide nuclear densities and averaged temperatures by
*         materials
          IF(LTOT.NE.(NUCL+1)*NK) CALL XABORT('FMAC01: BAD RECORD 8.')
          ALLOCATE(H2((NUCL+1),NK))
          READ(LIN,'(6E12.0)') ((H2(I,J),I=1,NUCL+1),J=1,NK)
          DEALLOCATE(H2)
        CASE(9)
*         nuclide temperatures by materials
          IF(LTOT.NE.NUCL*NK) CALL XABORT('FMAC01: BAD RECORD 9.')
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          DEALLOCATE(H1)
        CASE(10)
*         additional edit cross section names
          IF(LTOT.NE.NEDIT) CALL XABORT('FMAC01: BAD RECORD 10.')
          IF(IVERS.GE.5) THEN
            READ(LIN,'(12A6)') (TEX(I),I=1,LTOT)
          ELSE
            ALLOCATE(H1(LTOT))
            READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
            DEALLOCATE(H1)
          ENDIF
        CASE(11)
*         PMIN array
          IF(LTOT.NE.NGP) CALL XABORT('FMAC01: BAD RECORD 11.')
          READ(LIN,'(6I12)') (NPMIN(I),I=1,LTOT)
        CASE(12)
*         PMAX array
          IF(LTOT.NE.NGP) CALL XABORT('FMAC01: BAD RECORD 12.')
          READ(LIN,'(6I12)') (NPMAX(I),I=1,LTOT)
        CASE(13)
          CALL XABORT('FMAC01: This record (13) is UNDEFINED in the ve'
     1    //'rsion.ge.3 of format FMAC-M.')
        CASE(14)
*         number of scattering cross-section moments
          IF(LTOT.NE.NGP) CALL XABORT('FMAC01: BAD RECORD 14.')
          READ(LIN,'(6I12)') (NANIS(I),I=1,LTOT)
          IF(IPART.EQ.0) CALL XABORT('FMAC01: PARTICLE TYPE UNDEFINED.')
        CASE(15)
*         total cross sections
          IF(LTOT.NE.NK*NGP) CALL XABORT('FMAC01: BAD RECORD 15.')
          IF(IPART.EQ.0) CALL XABORT('FMAC01: PARTICLE TYPE UNDEFINED.')
          IG1=1
          DO I=1,IPART-1
            IG1=IG1+NGPRT(I)
          ENDDO
          IG2=IG1+NGPRT(IPART)-1
          ALLOCATE(H2(NK,NGP))
          READ(LIN,'(6E12.0)') ((H2(I,J),I=1,NK),J=1,NGP)
          CALL FMAC02(IPMACR,NK,IG2-IG1+1,H2(1,IG1),'NTOT0')
          DEALLOCATE(H2)
        CASE(16)
*         absorption cross sections
          IF(LTOT.NE.NK*NGP) CALL XABORT('FMAC01: BAD RECORD 16.')
          IF(IPART.EQ.0) CALL XABORT('FMAC01: PARTICLE TYPE UNDEFINED.')
          IG1=1
          DO I=1,IPART-1
            IG1=IG1+NGPRT(I)
          ENDDO
          IG2=IG1+NGPRT(IPART)-1
          ALLOCATE(H2(NK,NGP))
          READ(LIN,'(6E12.0)') ((H2(I,J),I=1,NK),J=1,NGP)
          CALL FMAC02(IPMACR,NK,IG2-IG1+1,H2(1,IG1),'ABS')
          DEALLOCATE(H2)
        CASE(17)
*         mixture corresponding to each fissile mixture
          IF(LTOT.NE.NUFIS) CALL XABORT('FMAC01: BAD RECORD 17.')
          READ(LIN,'(6I12)') (MUFIS(I),I=1,LTOT)
        CASE(18)
*         fission cross sections
          IF(LTOT.NE.NUFIS*NGP) CALL XABORT('FMAC01: BAD RECORD 18.')
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          DEALLOCATE(H1)
        CASE(19)
*         neutron production cross sections by fission material
          IF(LTOT.NE.NUFIS*NGP) CALL XABORT('FMAC01: BAD RECORD 19.')
          IF(NUFIS.EQ.0) CALL XABORT('FMAC01: NO FISSILE MIXTURES.')
          IF(NGROUP.EQ.0) CALL XABORT('FMAC01: NGROUP UNDEFINED.')
          ALLOCATE(H2(NUFIS,NGP),H3(NK,NUFIS,NGROUP))
          H3(:NK,:NUFIS,:NGROUP)=0.0
          READ(LIN,'(6E12.0)') ((H2(I,J),I=1,NUFIS),J=1,NGP)
          IF(IPART.EQ.0) CALL XABORT('FMAC01: PARTICLE TYPE UNDEFINED.')
          IG1=1
          DO I=1,IPART-1
            IG1=IG1+NGPRT(I)
          ENDDO
          IG2=IG1+NGPRT(IPART)-1
          DO I=1,NUFIS
            IK=MUFIS(I)
            IF((IK.LE.0).OR.(IK.GT.NK)) CALL XABORT('FMAC01: WRONG MUF'
     1      //'IS VALUE.')
            H3(IK,I,:NGROUP)=H2(I,IG1:IG2)
          ENDDO
          CALL FMAC02(IPMACR,NK*NUFIS,NGROUP,H3,'NUSIGF')
          DEALLOCATE(H3,H2)
        CASE(20)
*         fission spectra by fission material
          IF(LTOT.NE.NUFIS*NGXI) CALL XABORT('FMAC01: BAD RECORD 20.')
          IF(NUFIS.EQ.0) CALL XABORT('FMAC01: NO FISSILE MIXTURES.')
          IF(NGROUP.EQ.0) CALL XABORT('FMAC01: NGROUP UNDEFINED.')
          ALLOCATE(H2(NUFIS,NGXI),H3(NK,NUFIS,NGROUP))
          H3(:NK,:NUFIS,:NGROUP)=0.0
          READ(LIN,'(6E12.0)') ((H2(I,J),I=1,NUFIS),J=1,NGXI)
          IF(IPART.EQ.0) CALL XABORT('FMAC01: PARTICLE TYPE UNDEFINED.')
          IG1=1
          DO I=1,IPART-1
            IG1=IG1+NGPRT(I)
          ENDDO
          IG2=IG1+NGPRT(IPART)-1
          DO I=1,NUFIS
            IK=MUFIS(I)
            IF((IK.LE.0).OR.(IK.GT.NK)) CALL XABORT('FMAC01: WRONG MUF'
     1      //'IS VALUE.')
            H3(IK,I,:NGXI)=H2(I,IG1:IG2+NGXI-NGROUP)
          ENDDO
          CALL FMAC02(IPMACR,NK*NUFIS,NGROUP,H3,'CHI')
          DEALLOCATE(H3,H2)
        CASE(21)
*         fission nuclide numbers
          READ(LIN,'(6I12)') (N(I),I=1,LTOT)
        CASE(22)
*         summary parts of delayed fission neutrons by fission nuclides
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          DEALLOCATE(H1)
        CASE(23)
*         relative delayed group parts of delayed fission neutrons
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          DEALLOCATE(H1)
        CASE(24)
*         decay constants of delayed neutrons
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          DEALLOCATE(H1)
        CASE(25)
*         fission spectra of instantaneous fission neutrons
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          DEALLOCATE(H1)
        CASE(26)
*         fission spectra of delayed fission neutrons
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          DEALLOCATE(H1)
        CASE(27)
*         blocked microscopic neutron production cross sections
          ALLOCATE(H1(LTOT))
          READ(LIN,'(6E12.0)') (H1(I),I=1,LTOT)
          DEALLOCATE(H1)
        CASE(28:34)
          CALL XABORT('FMAC01: This record (28:34) is UNDEFINED in the'
     1    //' version.ge.3 of format FMAC-M.')
        CASE(35)
*         restricted stopping power
          IF(LTOT.NE.NK*(NGP+NPART)) CALL XABORT('FMAC01: INVALID RECO'
     1    //'RD 35.')
          IF(IPART.EQ.0) CALL XABORT('FMAC01: PARTICLE TYPE UNDEFINED.')
          IG1=1
          DO I=1,IPART-1
            IG1=IG1+NGPRT(I)+1
          ENDDO
          IG2=IG1+NGPRT(IPART)
          ALLOCATE(H2(NK,NGP+NPART))
          READ(LIN,'(6E12.0)') ((H2(I,J),I=1,NK),J=1,NGP+NPART)
          IF((HPART.EQ.'N').OR.(HPART.EQ.'G')) THEN
            DEALLOCATE(H2)
            CYCLE
          ENDIF
          CALL FMAC04(NGPRT,NGP,NPART,NK,H2)
          CALL FMAC02(IPMACR,NK,IG2-IG1+1,H2(1,IG1),'ESTOPW')
          DEALLOCATE(H2)
        CASE(36)
*         restricted momentum transfer cross section
          IF(LTOT.NE.NK*NGP) CALL XABORT('FMAC01: BAD RECORD 36.')
          IF(IPART.EQ.0) CALL XABORT('FMAC01: PARTICLE TYPE UNDEFINED.')
          IG1=1
          DO I=1,IPART-1
            IG1=IG1+NGPRT(I)
          ENDDO
          IG2=IG1+NGPRT(IPART)-1
          ALLOCATE(H2(NK,NGP))
          READ(LIN,'(6E12.0)') ((H2(I,J),I=1,NK),J=1,NGP)
          IF((HPART.EQ.'N').OR.(HPART.EQ.'G')) THEN
            DEALLOCATE(H2)
            CYCLE
          ENDIF
          CALL FMAC02(IPMACR,NK,IG2-IG1+1,H2(1,IG1),'EMOMTR')
          DEALLOCATE(H2)
        CASE(37)
*         energy deposition cross section
          IF(LTOT.NE.NK*NGP) CALL XABORT('FMAC01: BAD RECORD 37.')
          IF(IPART.EQ.0) CALL XABORT('FMAC01: PARTICLE TYPE UNDEFINED.')
          IG1=1
          DO I=1,IPART-1
            IG1=IG1+NGPRT(I)
          ENDDO
          IG2=IG1+NGPRT(IPART)-1
          ALLOCATE(H2(NK,NGP))
          READ(LIN,'(6E12.0)') ((H2(I,J),I=1,NK),J=1,NGP)
          CALL FMAC02(IPMACR,NK,IG2-IG1+1,H2(1,IG1),'H-FACTOR')
          DEALLOCATE(H2)
        CASE(38)
*         charge deposition cross section
          IF(LTOT.NE.NK*NGP) CALL XABORT('FMAC01: BAD RECORD 38.')
          IF(IPART.EQ.0) CALL XABORT('FMAC01: PARTICLE TYPE UNDEFINED.')
          IG1=1
          DO I=1,IPART-1
            IG1=IG1+NGPRT(I)
          ENDDO
          IG2=IG1+NGPRT(IPART)-1
          ALLOCATE(H2(NK,NGP))
          READ(LIN,'(6E12.0)') ((H2(I,J),I=1,NK),J=1,NGP)
          CALL FMAC02(IPMACR,NK,IG2-IG1+1,H2(1,IG1),'C-FACTOR')
          DEALLOCATE(H2)
        CASE DEFAULT
          CALL XABORT('FMAC01: This record (>38) is UNDEFINED in the'
     1    //' version.ge.3 of format FMAC-M.')
        END SELECT
      ENDDO
      RETURN
 1000 FORMAT(/9H FMAC01: ,12A6)
      END
