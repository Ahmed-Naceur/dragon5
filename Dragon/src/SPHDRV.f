*DECK SPHDRV
      SUBROUTINE SPHDRV(IPTRK,IFTRK,IPMACR,IPFLX,IPRINT,IMC,NGCOND,
     1 NMERGE,NALBP,IGRMIN,IGRMAX,SPH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the SPH factors. These factors are used to multiply the
* cross sections and to divide the fluxes. The SPH factors calculation
* is generally application dependent. New SPH algorithms should be
* implemented in this driver.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK   pointer to the macro-tracking LCM object.
* IFTRK   unit of the macro-tracking binary sequential file.
* IPMACR  pointer to the Macrolib (L_MACROLIB signature).
* IPFLX   pointer towards an initialization flux (L_FLUX signature).
* IPRINT  print flag (equal to 0 for no print).
* IMC     type of macro-calculation (=1 diffusion or SPN; 
*         =2 other options;
*         =3 type PIJ with Bell acceleration).
* NGCOND  number of condensed groups.
* NMERGE  number of merged regions.
* NALBP   number of physical albedos.
* IGRMIN  first group to process.
* IGRMAX  last group to process.
*
*Parameters: output
* SPH     SPH homogenization factors.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPMACR,IPFLX
      INTEGER   IFTRK,IPRINT,IMC,NGCOND,NMERGE,NALBP,IGRMIN,IGRMAX
      REAL      SPH(NMERGE+NALBP,NGCOND)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER TEXT12*12,CNDOOR*12,CTITRE*72,SUFF(2)*2,HSMG*131
      INTEGER   ISTATE(NSTATE)
      LOGICAL   ILK
      TYPE(C_PTR) IPSPH
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT2,KEY2,MERG2
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL2
      DATA      SUFF/'00','01'/
*----
*  RECOVER SPH-RELATED INFORMATION
*----
      CALL LCMLEN(IPMACR,'SPH',ILONG,ITYLCM)
      IF(ILONG.EQ.0) CALL XABORT('SPHDRV: MISSING SPH DIRECTORY.')
      IPSPH=LCMDID(IPMACR,'SPH')
      CALL LCMGET(IPSPH,'STATE-VECTOR',ISTATE)
      NSPH=ISTATE(1)
      KSPH=ISTATE(2)
      MAXIT=ISTATE(3)
      MAXNBI=ISTATE(4)
      IF((NSPH.EQ.0).OR.(NSPH.EQ.1)) CALL XABORT('SPHDRV: INVALID VALU'
     > //'E OF NSPH.')
*----
*  RECOVER AND USE AN EXISTING MACRO-TRACKING.
*----
      IF(C_ASSOCIATED(IPTRK)) THEN
         IF(NSPH.GE.2) THEN
            CALL LCMGTC(IPSPH,'SPH$TRK',12,1,CNDOOR)
            CALL LCMGET(IPSPH,'SPH-EPSILON',EPSPH)
         ENDIF
         CALL LCMGTC(IPTRK,'SIGNATURE',12,1,TEXT12)
         IF(TEXT12.NE.'L_TRACK') THEN
            CALL XABORT('SPHDRV: TRACKING DATA STRUCTURE EXPECTED.')
         ENDIF
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
         NREG2=ISTATE(1)
         NUN2=ISTATE(2)
         ILK=ISTATE(3).EQ.0
         CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,CNDOOR)
         IF(CNDOOR.EQ.'MCCG') THEN
            CALL LCMLEN(IPTRK,'KEYFLX',LKFL,ITYLCM)
            NFUNL=LKFL/NREG2
         ELSE
            NFUNL=1
         ENDIF
         IF((CNDOOR.EQ.'SYBIL').AND.(NSPH.EQ.4)) NUN2=NUN2+ISTATE(9)
         IF((CNDOOR.EQ.'EXCELL').OR.(CNDOOR.EQ.'MCCG')) THEN
            ISCAT=ISTATE(6)
         ELSE IF(CNDOOR.EQ.'SN') THEN
            ISCAT=ISTATE(16)
         ELSE IF(CNDOOR.EQ.'BIVAC') THEN
            IF(ISTATE(14).EQ.0) THEN
               ISCAT=1
            ELSE
               ISCAT=ISTATE(16)
            ENDIF
         ELSE IF(CNDOOR.EQ.'TRIVAC') THEN
            IF(ISTATE(30).EQ.0) THEN
               ISCAT=1
            ELSE
               ISCAT=ISTATE(32)
            ENDIF
         ELSE
            ISCAT=1
         ENDIF
         ISCAT=ABS(ISCAT)
         ALLOCATE(VOL2(NREG2),MAT2(NREG2),KEY2(NREG2*NFUNL))
         CALL LCMGET(IPTRK,'VOLUME',VOL2)
         CALL LCMGET(IPTRK,'MATCOD',MAT2)
         CALL LCMGET(IPTRK,'KEYFLX',KEY2)
         CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            CALL LCMGTC(IPTRK,'TITLE',72,1,CTITRE)
         ELSE
            CTITRE='*** NO TITLE PROVIDED ***'
         ENDIF
         NBMIX2=0
         IF(KSPH.EQ.5) THEN
*          HEBERT-BENOIST ALBS TECHNIQUE.
           DO 20 IREG=1,NREG2
           NBMIX2=MAX(NBMIX2,MAT2(IREG))
   20      CONTINUE
           ALLOCATE(MERG2(NBMIX2))
           DO 25 IBM=1,NBMIX2
           MERG2(IBM)=1
   25      CONTINUE
           ILK=.FALSE.
         ELSE
           DO 30 IREG=1,NREG2
           NBMIX2=MAX(NBMIX2,MAT2(IREG))
   30      CONTINUE
           IF(NBMIX2.NE.NMERGE) THEN
              WRITE(HSMG,'(41HSPHDRV: INVALID NUMBER OF MACRO-REGIONS (,
     1        2I6,2H).)') NBMIX2,NMERGE
              CALL XABORT(HSMG)
           ENDIF
           ALLOCATE(MERG2(NBMIX2))
           DO 35 IBM=1,NBMIX2
           MERG2(IBM)=IBM
   35      CONTINUE
         ENDIF
*
*        RECOVER TABULATED FUNCTIONS.
         CALL XDRTA2(IPTRK)
      ELSE
         ISCAT=1
         ILK=.FALSE.
         NBMIX2=0
         NREG2=0
         NUN2=0
      ENDIF
*----
*  GENERAL PROCEDURE FOR COMPUTING THE SPH FACTORS
*----
      CALL SPHEQU(NBMIX2,IPTRK,IFTRK,IPMACR,IPFLX,CNDOOR,NSPH,KSPH,
     1 MAXIT,MAXNBI,EPSPH,IPRINT,IMC,NGCOND,NMERGE,NALBP,ISCAT,NREG2,
     2 NUN2,MAT2,VOL2,KEY2,MERG2,ILK,CTITRE,IGRMIN,IGRMAX,SPH)
      IF(C_ASSOCIATED(IPTRK)) DEALLOCATE(MERG2,KEY2,MAT2,VOL2)
*----
*  PRINT SPH FACTORS
*----
      IF(IPRINT.GT.1) THEN
         WRITE(6,'(/21H SPHDRV: SPH FACTORS:)')
         WRITE(6,200) ((IKK,IGR,SPH(IKK,IGR),IKK=1,NMERGE+NALBP),IGR=1,
     >   NGCOND)
      ENDIF
      RETURN
*
  200 FORMAT(4X,4HSPH(,I5,1H,,I3,2H)=,F9.5,:,4X,4HSPH(,I5,1H,,I3,2H)=,
     > F9.5,:,4X,4HSPH(,I5,1H,,I3,2H)=,F9.5,:,4X,4HSPH(,I5,1H,,I3,2H)=,
     > F9.5,:,4X,4HSPH(,I5,1H,,I3,2H)=,F9.5)
      END
