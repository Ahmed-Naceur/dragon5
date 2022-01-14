*DECK LIBLIC
      SUBROUTINE LIBLIC (IPLIB,NBISO,MASKI,IMPX,NGRO,NL,ITRANC,ITIME,
     1 NLIB,NED,HVECT,ISONAM,ISONRF,IPISO,ISHINA,TMPISO,IHLIB,ILLIB,
     2 INAME,NTFG,SN,SB,NIR,GIR,NGF,IGRMAX,NDEL,NBESP,NPART,IPROC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transcription of the useful interpolated microscopic cross section
* data from various format of libraries to lcm. A two dimensional
* interpolation in temperature and dilution is performed (Part B).
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
* NBISO   number of isotopes present in the calculation domain.
* MASKI   isotopic masks. An isotope with index I is processed if
*         MASKI(I)=.true.
* IMPX    print flag.
* NGRO    number of energy groups.
* NL      number of Legendre orders required in the calculation
*         NL=1 (for isotropic scattering) or higher.
* ITRANC  type of transport correction: =0 no transport correction
*         =1 Apollo type transport correction; =2 recover from
*         library; =3 Wims-D type; =4 leakage correction alone.
* ITIME   MATXS type of fission spectrum:
*         =1 steady-state; =2 prompt.
* NLIB    number of independent libraries.
* NED     number of requested vector edits.
* HVECT   names of the requested vector edits.
* ISONAM  alias name of each isotope.
* ISONRF  library reference name of each isotope.
* IPISO   pointer array towards microlib isotopes.
* ISHINA  self-shielding name of each isotope.
* TMPISO  temperature of each isotope.
* IHLIB   isotope options.
* ILLIB   xs library index for each isotope (.le.NLIB).
* INAME   names of the NLIB xs libraries.
* NTFG    number of thermal groups where the thermal inelastic
*         correction is applied.
* SN      dilution cross section in each energy group of each
*         isotope. a value of 1.0E10 is used for infinite dilution.
* SB      dilution cross section as used in Livolant and Jeanpierre
*         normalization.
* NIR     first group index with an imposed IR slowing-down model;
*         =0 for no IR model.
* GIR     value of the imposed Goldstein-Cohen parameter for groups
*         with an IR model.
* NGF     number of fast groups without self-shielding.
* IGRMAX  maximum group index with self-shielding.
* NDEL    number of precursor groups for delayed neutrons.
* NBESP   number of energy-dependent fission spectra.
* NPART   number of particles.
* IPROC   type of library processing.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NBISO,IMPX,NGRO,NL,ITRANC,ITIME,NLIB,NED,NGF,IGRMAX,NDEL,
     > NBESP,NPART,IPROC,ISONAM(3,NBISO),ISONRF(3,NBISO),
     > ISHINA(3,NBISO),IHLIB(2,NBISO,4),ILLIB(NBISO),INAME(16,NLIB),
     > NTFG(NBISO),NIR(NBISO)
      LOGICAL MASKI(NBISO)
      CHARACTER*(*) HVECT(NED)
      REAL TMPISO(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO),GIR(NBISO)
*----
*  INTERNAL PARAMETERS
*----
      TYPE(C_PTR) IPDRL,IPMIC
      INTEGER MAXDIL
      PARAMETER (MAXDIL=65)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,J,IND1,NBIS,NDEL0,NGF0,NGFR0,MAXTRA,ILIB,ILONG,NBESP0,
     1 NPART0
      CHARACTER NAMLBT*8,NAMFIL*64,HSMG*131,NAMLCM*12,NAMMY*12
      LOGICAL LTEST,EMPTY,LCM
*----
*  CHECK FOR DUPLICATE ISOTOPE NAMES.
*----
      DO 40 I=1,NBISO
      IF(MASKI(I).AND.(ILLIB(I).NE.0)) THEN
         DO 30 J=I+1,NBISO
         IF(MASKI(J).AND.(ISONAM(1,I).EQ.ISONAM(1,J)).AND.
     1                   (ISONAM(2,I).EQ.ISONAM(2,J)).AND.
     2                   (ISONAM(3,I).EQ.ISONAM(3,J))) THEN
            WRITE (HSMG,200) ISONAM(1,I),ISONAM(2,I),ISONAM(3,I)
            CALL XABORT(HSMG)
         ENDIF
   30    CONTINUE
      ENDIF
   40 CONTINUE
*
      NPART=1
      NGF0=NGRO+1
      NGFR0=0
      IND1=1
   50 NBIS=1
      LTEST=MASKI(IND1)
      DO 60 I=IND1+1,NBISO
        IF((IHLIB(1,I,1).EQ.IHLIB(1,IND1,1)).AND.
     1     (IHLIB(2,I,1).EQ.IHLIB(2,IND1,1)).AND.
     2     (ILLIB(I).EQ.ILLIB(IND1))) THEN
          NBIS=NBIS+1
          LTEST=LTEST.OR.MASKI(I)
        ELSE
          GO TO 70
        ENDIF
   60 CONTINUE
   70 WRITE(NAMLBT,'(2A4)') IHLIB(1,IND1,1),IHLIB(2,IND1,1)
      ILIB=ILLIB(IND1)
      IF(ILIB.EQ.0) THEN
         NAMFIL=' '
      ELSE
         WRITE(NAMFIL,'(16A4)') (INAME(I,ILIB),I=1,16)
      ENDIF
      NDEL0=0
      IF(LTEST.AND.(NAMLBT.EQ.'DRAGON')) THEN
*        TRANSFER INFORMATION FROM DRAGON LIBRARY TO LCM.
         CALL LCMINF(IPLIB,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
         IF(NAMFIL.EQ.NAMLCM) THEN
           IPDRL=IPLIB
         ELSE
           CALL LCMOP(IPDRL,NAMFIL(:12),2,2,0)
         ENDIF
         CALL LIBDRA (IPLIB,IPDRL,NAMFIL,NGRO,NBIS,NL,ISONAM(1,IND1),
     1   ISONRF(1,IND1),IPISO(IND1),TMPISO(IND1),SN(1,IND1),SB(1,IND1),
     2   MASKI(IND1),NED,HVECT,IMPX,NGF0,NGFR0,NDEL0,NBESP0)
         IF(NAMFIL.NE.NAMLCM) CALL LCMCL(IPDRL,1)
         NBESP=MAX(NBESP,NBESP0)
      ELSE IF(LTEST.AND.(NAMLBT(1:4).EQ.'WIMS')) THEN
*        TRANSFER INFORMATION FROM WIMS LIBRARY FILE TO LCM.
         IF(NAMLBT.EQ.'WIMSD4') THEN
*           WIMS-D4 FORMAT
            CALL LIBWD4(IPLIB,IMPX,NAMFIL,NGRO,NBIS,NL,ISONAM(1,IND1),
     1      ISONRF(1,IND1),IPISO(IND1),ISHINA(1,IND1),TMPISO(IND1),
     2      SN(1,IND1),SB(1,IND1),MASKI(IND1),NGF0,NGFR0)
         ELSE IF(NAMLBT.EQ.'WIMSE') THEN
*           WIMS-E FORMAT
            CALL LIBWE(IPLIB,IMPX,NAMFIL,NGRO,NBIS,NL,ISONAM(1,IND1),
     1      ISONRF(1,IND1),IPISO(IND1),ISHINA(1,IND1),TMPISO(IND1),
     2      SN(1,IND1),SB(1,IND1),MASKI(IND1),NGF0,NGFR0)
         ELSE IF(NAMLBT.EQ.'WIMSAECL') THEN
*           WIMS-AECL FORMAT
            CALL LIBWIM(IPLIB,IMPX,NAMFIL,NGRO,NBIS,NL,ISONAM(1,IND1),
     1      ISONRF(1,IND1),IPISO(IND1),ISHINA(1,IND1),TMPISO(IND1),
     2      SN(1,IND1),SB(1,IND1),MASKI(IND1),NGF0,NGFR0)
         ENDIF
      ELSE IF(LTEST.AND.(NAMLBT.EQ.'MATXS')) THEN
*        TRANSFER INFORMATION FROM MATXS (NJOY-89) TO LCM.
         CALL LIBTR1 (IPLIB,NAMFIL,NGRO,NBIS,NL,ISONAM(1,IND1),
     1   ISONRF(1,IND1),IPISO(IND1),IHLIB(1,IND1,2),IHLIB(1,IND1,3),
     2   NTFG(IND1),TMPISO(IND1),SN(1,IND1),SB(1,IND1),MASKI(IND1),NED,
     3   HVECT,ITIME,IMPX,NGF0,NGFR0)
      ELSE IF(LTEST.AND.(NAMLBT.EQ.'MATXS2')) THEN
*        TRANSFER INFORMATION FROM MATXS (NJOY-91) TO LCM.
         CALL LIBTR2 (IPLIB,NAMFIL,NGRO,NBIS,NL,ISONAM(1,IND1),
     1   ISONRF(1,IND1),IPISO(IND1),IHLIB(1,IND1,2),IHLIB(1,IND1,3),
     2   IHLIB(1,IND1,4),NTFG(IND1),TMPISO(IND1),SN(1,IND1),SB(1,IND1),
     3   MASKI(IND1),NED,HVECT,ITIME,IMPX,NGF0,NGFR0,NPART0)
         NPART=MAX(NPART,NPART0)
      ELSE IF(LTEST.AND.(NAMLBT.EQ.'APLIB1')) THEN
*        TRANSFER INFORMATION FROM APOLIB-1 TO LCM.
         MAXTRA=NL*NGRO**2
         CALL LIBAPL (IPLIB,NAMFIL,MAXTRA,NGRO,NBIS,NL,ISONAM(1,IND1),
     1   ISONRF(1,IND1),IPISO(IND1),ISHINA(1,IND1),MASKI(IND1),
     2   TMPISO(IND1),SN(1,IND1),SB(1,IND1),IMPX,NGF0,NGFR0)
      ELSE IF(LTEST.AND.(NAMLBT.EQ.'APLIB2')) THEN
*        TRANSFER INFORMATION FROM APOLIB-2 TO LCM.
         CALL LIBA20 (IPLIB,NAMFIL,NGRO,NBIS,NL,IPROC,ISONAM(1,IND1),
     1   ISONRF(1,IND1),IPISO(IND1),ISHINA(1,IND1),MASKI(IND1),
     2   TMPISO(IND1),SN(1,IND1),SB(1,IND1),IMPX,NGF0,NGFR0,NDEL0)
      ELSE IF(LTEST.AND.(NAMLBT.EQ.'APXSM')) THEN
*        TRANSFER INFORMATION FROM APOLIB-XSM TO LCM.
         CALL LIBXS4 (IPLIB,NAMFIL,NGRO,NBIS,NL,IPROC,ISONAM(1,IND1),
     1   ISONRF(1,IND1),IPISO(IND1),ISHINA(1,IND1),MASKI(IND1),
     2   TMPISO(IND1),SN(1,IND1),SB(1,IND1),IMPX,NGF0,NGFR0,NDEL0)
      ELSE IF(LTEST.AND.(NAMLBT.EQ.'NDAS')) THEN
         CALL LIBND1 (IPLIB,NAMFIL,NGRO,NBIS,NL,ISONAM(1,IND1),
     1   ISONRF(1,IND1),IPISO(IND1),MASKI(IND1),TMPISO(IND1),SN(1,IND1),
     2   SB(1,IND1),IMPX,NGF0,NGFR0,NDEL0)
      ELSE IF(LTEST.AND.(NAMLBT.EQ.'MICROLIB')) THEN
*        TRANSFER INFORMATION FROM MICROLIB LIBRARY TO LCM.
         CALL LCMINF(IPLIB,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
         IF(NAMFIL.EQ.NAMLCM) THEN
           IPMIC=IPLIB
         ELSE
           CALL LCMOP(IPMIC,NAMFIL(:12),2,2,0)
         ENDIF
         CALL LIBMIC (IPLIB,IPMIC,NAMFIL,NGRO,NBIS,ISONAM(1,IND1),
     1   ISONRF(1,IND1),IPISO(IND1),MASKI(IND1),IMPX,NGF0,NGFR0,NDEL0,
     2   NBESP0)
         IF(NAMFIL.NE.NAMLCM) CALL LCMCL(IPMIC,1)
         NBESP=MAX(NBESP,NBESP0)
      ENDIF
      IF(LTEST) THEN
         NGF=MIN(NGF,NGF0)
         IGRMAX=MAX(IGRMAX,NGFR0)
         IF(NDEL.EQ.0) THEN
            NDEL=NDEL0
         ELSE IF((NDEL0.NE.NDEL).AND.(NDEL0.NE.0)) THEN
            ILIB=ILLIB(IND1)
            IF(ILIB.GT.0) WRITE(6,210) (INAME(I,ILIB),I=1,4),NDEL0,NDEL
            NDEL=MAX(NDEL,NDEL0)
         ENDIF
*
*        COMPUTE THE TRANSPORT XS AND ADD COMPLEMENTARY INFORMATION.
         CALL LIBADD (IPLIB,NBIS,MASKI(IND1),IMPX,NGRO,NL,ITRANC,
     1   ISONAM(1,IND1),IPISO(IND1),NIR(IND1),GIR(IND1))
      ENDIF
*
      IND1=IND1+NBIS
      IF(IND1.LE.NBISO) GO TO 50
      RETURN
*
  200 FORMAT(8HLIBLIC: ,3A4,34H IS A DUPLICATE ISOTOPE/MATERIAL N,
     1 4HAME.)
  210 FORMAT(/51H LIBLIC: INVALID NB OF PRECURSOR GROUPS IN LIBRARY ,
     1 4A4,8H (NDEL0=,I3,6H NDEL=,I3,2H).)
      END
