*DECK FMAC
      SUBROUTINE FMAC(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Convert a macroscopic cross section file in ascii FMAC-M format
* towards Version5 macrolib format.
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
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1) create type(L_MACROLIB);
*         HENTRY(2) read-only ascii file containing FMAC-M data.
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPMACR
      PARAMETER(NSTATE=40,IOUT=6,MAS=38)
      CHARACTER TEXT12*12,TEXT18*18,HSIGN*12,TEX(50)*6,HPART*1
      INTEGER ISTATE(NSTATE),N(MAS),IZA(40)
      DOUBLE PRECISION DFLOTT
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NGPRT,NPMIN,NPMAX,NANIS,
     1 MUFIS
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NWA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: H2
      CHARACTER(LEN=1), ALLOCATABLE, DIMENSION(:) :: HNPRT
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.NE.2) CALL XABORT('FMAC: TWO PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('FMAC: LCM'
     1 //' OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0) CALL XABORT('FMAC: ENTRY IN CREATE MODE EXPEC'
     1 //'TED.')
      IPMACR=KENTRY(1)
      HSIGN='L_MACROLIB'
      CALL LCMPTC(IPMACR,'SIGNATURE',12,1,HSIGN)
*----
*  RECOVER FMAC-M FILE
*----
      TEXT12=HENTRY(2)
      IF(IENTRY(2).NE.4) CALL XABORT('FMAC: ASCII FILE NAMED '//TEXT12
     1 //' EXPECTED AT LHS.')
      IF(JENTRY(2).NE.2) CALL XABORT('FMAC: ASCII FILE IN READ-ONLY MO'
     1 //'DE EXPECTED.')
      LIN=FILUNIT(KENTRY(2))
*----
*  READ THE INPUT DATA
*----
      IMPX=1
      HPART=' '
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('FMAC: CHARACTER DATA EXPECTED.')
      IF(TEXT12.EQ.'EDIT') THEN
*       READ THE PRINT INDEX.
        CALL REDGET(INDIC,IMPX,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('FMAC: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT12.EQ.'PARTICLE') THEN
*       READ THE PARTICLE TYPE ('N', 'G', 'C', 'P')
        CALL REDGET(INDIC,NITMA,FLOTT,HPART,DFLOTT)
        IF(INDIC.NE.3) CALL XABORT('FMAC: CHARACTER DATA EXPECTED.')
        CALL LCMPTC(IPMACR,'PARTICLE',1,1,HPART)
      ELSE IF(TEXT12.EQ.';') THEN
        GO TO 20
      ELSE
        CALL XABORT('FMAC: '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
      GO TO 10
*----
*  MACROLIB INITIALIZATION
*----
   20 IPART=0
      NANISO=0
      ITRANC=0
      NALBP=0
      NSTEP=0
      IDF=0
*----
*  PROCESS THE FMAC FILE FOR 1D COUPLED CROSS-SECTION DATA
*----
      READ(LIN,'(A18,I12)') TEXT18,IVERS
      IF(IMPX.GT.0) WRITE(IOUT,1002) TEXT18,IVERS
      READ(LIN,'(A72)')
      READ(LIN,'(A30)')
      READ(LIN,'(6I12)') NTYPE,(N(I),I=1,9)
      READ(LIN,'(12A6)')(TEX(I),I=1,NTYPE)
      IF(TEX(1).NE.'NGCAP ') THEN
        CALL XABORT('FMAC: INVALID DATA TYPE='//TEX(1)//'.')
      ENDIF
      READ(LIN,'(6I12)') (IZA(I),I=1,40)
      WRITE(IOUT,1101) (IZA(I),I=1,40)
      NGP=IZA(1)
      NPART=IZA(2)
      NGXI=IZA(3)
      NGXIF=IZA(4)
      NDELG=IZA(5)
      NEDIT=IZA(8)
      KLEIN=IZA(10)
      IADJ=IZA(12)
      NUCL=IZA(13)
      NK=IZA(15)
      NUFIS=IZA(16)
      MASM3=IZA(17)
      NWO=IZA(20)
      MPOINT=IZA(21)
      ALLOCATE(NGPRT(NPART),HNPRT(NPART),NPMIN(NGP),NPMAX(NGP),
     1 NANIS(NGP),MUFIS(NUFIS))
      IF(MASM3.GT.38) CALL XABORT('FMAC: MASM3 OVERFLOW.')
      READ(LIN,'(6I12)')(N(I),I=1,MASM3)
      CALL FMAC01(IPMACR,IMPX,HPART,LIN,IVERS,NGP,NPART,NGXI,NEDIT,
     1 NUCL,NK,NUFIS,MASM3,N,NGPRT,HNPRT,NPMIN,NPMAX,NANIS,MUFIS)
*----
*  PROCESS SCATTERING INFORMATION
*----
      DO I=1,NPART
        IF(HNPRT(I).EQ.HPART) THEN
          IPART=I
          GO TO 30
        ENDIF
      ENDDO
      CALL XABORT('FMAC: PARTICLE '//HPART//' NOT AVAILABLE IN FMAC-M '
     1 //'FILE.')
   30 IF(IPART.EQ.0) CALL XABORT('FMAC: PARTICLE TYPE NOT DEFINED.')
      IG1=1
      DO I=1,IPART-1
        IG1=IG1+NGPRT(I)
      ENDDO
      IG2=IG1+NGPRT(IPART)-1
      DO IG=IG1,IG2
        NANISO=MAX(NANISO,NANIS(IG))
      ENDDO
      DO IG=1,NGP
        MAXLEN=2*(NANISO+1)*NK
        ALLOCATE(H2(NGP,MAXLEN),NWA(NGP,NK))
        H2(:NGP,:MAXLEN)=0.0
        NWA(:NGP,:NK)=0
        DO NP=NPMIN(IG),NPMAX(IG)
          READ(LIN,'(6I12)') NPP,NQQ,(NWA(NPP,I),I=1,NK),LENGTH
          IF(LENGTH.GT.MAXLEN) CALL XABORT('FMAC: MAXLEN OVERFLOW.')
          READ(LIN,'(6E12.0)') (H2(NPP,J),J=1,LENGTH)
        ENDDO
        IF((IG.GE.IG1).AND.(IG.LE.IG2)) THEN
          CALL FMAC03(IPMACR,IG,IPART,NGP,MAXLEN,NANISO,NK,NPART,HNPRT,
     1    NGPRT,NWA,H2)
        ENDIF
        DEALLOCATE(NWA,H2)
      ENDDO
*----
*  WRITE THE STATE VECTOR AND EXIT
*----
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NGPRT(IPART)
      ISTATE(2)=NK
      ISTATE(3)=NANISO
      ISTATE(4)=NUFIS
      ISTATE(5)=0
      ISTATE(6)=ITRANC
      ISTATE(7)=NDELG
      ISTATE(8)=NALBP
      ISTATE(11)=NSTEP
      ISTATE(12)=IDF
      ISTATE(13)=IADJ
      ISTATE(17)=NPART-1
      CALL LCMPUT(IPMACR,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(IMPX.GT.1) CALL LCMLIB(IPMACR)
      IF(IMPX.GT.0) THEN
        WRITE(IOUT,1010) IMPX,(ISTATE(I),I=1,7),ISTATE(9),ISTATE(11),
     1  ISTATE(12),ISTATE(13),ISTATE(17),HPART
        WRITE(IOUT,1020) (HNPRT(I),I=1,NPART)
        WRITE(IOUT,1030) (NGPRT(I),I=1,NPART)
      ENDIF
      DEALLOCATE(MUFIS,NANIS,NPMAX,NPMIN,HNPRT,NGPRT)
      RETURN
*
 1101 FORMAT(1X,'NG    =',I3,'   number of energy groups;',/
     +,1X,'NPART =',I3,'   number of particle types;',/
     +,1X,'NGXI  =',I3,'   number of groups with non-zero fission spectr
     +um;',/
     +,1X,'NGXIF =',I3,'   number of the first group with non-zero fissi
     +on spectrum;',/
     +,1X,'NGRET =',I3,'   number of delayed neutron groups;',/
     +,1X,'NGXIR =',I3,'   number of groups with non-zero fission spectr
     +um for delayed neutrons;',/
     +,1X,'NGXIFR=',I3,'   number of the first group with non-zero fissi
     +on spectrum for delayed neutrons;',/
     +,1X,'NEDIT =',I3,'   number of additional edit cross-sections;'
     +,/
     +,1X,'MAXGS =',I3,'   not used;',/
     +,1X,'KLEIN =',I3,'   not used;',/
     +,1X,'NGHIGH=',I3,'   total number of groups of cascade region;',/
     +,1X,'IADJ  =',I3,'   0/1 - regular/adjoint cross-section file;',/
     +,1X,'NUCL  =',I3,'   number of nuclides;',/
     +,1X,'NUCLF =',I3,'   number of fission nuclides;',/
     +,1X,'MIX   =',I3,'   number of materials (compositions);',/
     +,1X,'NUFIS =',I3,'   number of fission materials;',/
     +,1X,'MAS   =',I3,'   length of integer control array LL(MAS);',/
     +,1X,'KIN=',I6,'   total number of scattering transitions (not used
     +);',/
     +,1X,'MX    =',I3,'   maximal length of transition array (not used)
     +;',/
     +,1X,'MNW   =',I3,'   order of PL approximation used +1;',/
     +,1X,'MPOINT=',I3,'   order of discrete approximation of scattering
     + indicatrix (number of angular scattering cosine points used);',/
     +,1X,'NUMD(I)=',19I3,' not used')
 1002 FORMAT(1X,A18/' FMAC:  Format FMAC-M version =',I2)
 1010 FORMAT(/8H OPTIONS/8H -------/
     1 7H IMPX  ,I6,30H   (0=NO PRINT/1=SHORT/2=MORE)/
     2 7H NGROUP,I6,28H   (NUMBER OF ENERGY GROUPS)/
     3 7H NBMIX ,I6,39H   (NUMBER OF MIXTURES IN THE MACROLIB)/
     4 7H NANISO,I6,34H   (MAXIMUM SCATTERING ANISOTROPY)/
     5 7H NIFISS,I6,45H   (MAXIMUM NUMBER OF FISSILE ISOTOPES IN A M,
     6 7HIXTURE)/
     7 7H NEDMAC,I6,34H   (NUMBER OF CROSS SECTION EDITS)/
     8 7H ITRANC,I6,45H   (0=NO TRANSPORT CORRECTION/1=APOLLO TYPE/2,
     9 43H=RECOVER FROM LIBRARY/4=LEAKAGE CORRECTION)/
     1 7H NLG   ,I6,39H   (NUMBER OF DELAYED PRECURSOR GROUPS)/
     2 7H ILEAK ,I6,40H   (1=DIFF AVAILABLE; 2=DIFFX AVAILABLE)/
     3 7H NSTEP ,I6,39H   (NUMBER OF PERTURBATION DIRECTORIES)/
     4 7H IDF   ,I6,48H   (=0/2 BOUNDARY FLUXES FOR ADF ABSENT/PRESENT)/
     5 7H IADJ  ,I6,33H   (=0/1 DIRECT/ADJOINT MACROLIB)/
     6 7H NPART0,I6,34H   (NUMBER OF COMPANION PARTICLES)/
     7 6H HPART ,A7,22H   (MACROLIB PARTICLE))
 1020 FORMAT(/22H       PARTICLE NAMES:,10A8)
 1030 FORMAT(22H NB. OF ENERGY GROUPS:,10I8)
      END
