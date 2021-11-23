*DECK LIBWET
      SUBROUTINE LIBWET(MAXR,NEL,NSTATE,NMDEPL,ITNAM,ISTATE,MATNO,KPAX,
     >                  BPAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Identify the depleting isotopes by type and reorder them (recompute
* the KPAX and BPAX matrices).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input
* MAXR    number of reaction types.
* NEL     number of isotopes on library.
* NSTATE  number of parameters in the state vector.
* NMDEPL  names of reactions:
*           NMDEPL(1)='DECAY'; NMDEPL(2)='NFTOT';
*           NMDEPL(3)='NG'   ; NMDEPL(4)='N2N';
*           etc.
* ITNAM   reactive isotope names in chain.
*
*Parameters: output
* ISTATE  state vector containing the library parameters.
* MATNO   reaction material indice.s
*
*Parameters: input/output
* KPAX    complete reaction type matrix.
* BPAX    complete branching ratio matrix.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      MAXR,NEL,NSTATE,ITNAM(3,NEL),ISTATE(NSTATE),
     >             MATNO(NEL),KPAX(NEL+MAXR,NEL)
      CHARACTER    NMDEPL(MAXR)*8
      REAL         BPAX(NEL+MAXR,NEL)
*----
*  LOCAL VARIABLES
*----
      INTEGER      NDFP,JSO,ISO,ITR,NUNTIL,NHEAVY,NDFI,IUNTIL,NEW,
     >             NLIGHT,NOTHER,NIFP,KPAX0,KPAX1,KPAX2,IEL,JEL,
     >             KREAC,II,NPAR,ICOUNT,NREAC,NSTABL
      REAL         SUMFI
      CHARACTER    TEXT12*12
      LOGICAL      LOGPF
*----
*  INTERNAL PARAMETERS
*----
      INTEGER      KDECAY,KFISSP,KCAPTU,KFISSI,KHEAT
      PARAMETER   (KDECAY=1,KFISSP=2,KCAPTU=3,KFISSI=6,KHEAT=9)
*----
*  COMPUTE NUMBER OF DIRECT FISSION PRODUCT
*  NUMBER DIRECT FISSION PRODUCT
*
*  SET MATNO TO LIGHT (-KFISSP) FOR FISSION PRODUCT
*  SET MATNO TO HEAVY (-KFISSI) FOR FISSILE ISOTOPES
*  SET MATNO TO STABLE (-KHEAT) FOR STABLE ISOTOPES PRODUCING ENERGY
*  SET MATNO TO 0 TO UNUSED ISOTOPES
*----
      NDFI=0
      NDFP=0
      DO 100 JSO=NEL,1,-1
        MATNO(JSO)=0
        KPAX0=KPAX(NEL+KDECAY,JSO)
        KPAX1=KPAX(NEL+KFISSP,JSO)
        KPAX2=KPAX(NEL+KCAPTU,JSO)
        IF((KPAX0.EQ.-9999).OR.(KPAX1.EQ.-9999).OR.(KPAX2.EQ.-9999))
     >  THEN
*         JSO IS A STABLE ISOTOPE
          MATNO(JSO)=-KHEAT
          DO 222 ITR=1,MAXR
            IF(BPAX(NEL+ITR,JSO).GT.0.0) THEN
              KPAX(NEL+ITR,JSO)=2
              KPAX(JSO,:NEL)=0
            ENDIF
 222      CONTINUE
          GO TO 100
        ELSE IF(KPAX1.LT.0) THEN
*         JSO IS A FISSION PRODUCT
          MATNO(JSO)=-KFISSP
          NIFP=0
          DO 101 ISO=NEL,1,-1
            IF((KPAX(JSO,ISO).EQ.KFISSP).AND.
     >         (BPAX(JSO,ISO).GT.0.0)) THEN
              NIFP=NIFP+1
            ENDIF
 101      CONTINUE
          IF(NIFP.GT.0) THEN
            NDFP=NDFP+1
            KPAX(NEL+KFISSP,JSO)=5+100*NDFP
          ELSE
            KPAX(NEL+KFISSP,JSO)=5
          ENDIF
        ELSE IF(KPAX1.GT.0) THEN
*         JSO IS A FISSILE ISOTOPE
          MATNO(JSO)=-KFISSI
          SUMFI=0.0
          DO 221 ISO=1,NEL
            IF(KPAX(ISO,JSO).EQ.KFISSP) SUMFI=SUMFI+BPAX(ISO,JSO)
 221      CONTINUE
          IF(SUMFI.GT.0.0) THEN
            NDFI=NDFI+1
            KPAX(NEL+KFISSP,JSO)=4+100*NDFI
          ELSE
            KPAX(NEL+KFISSP,JSO)=3
          ENDIF
        ENDIF
 100  CONTINUE
*----
*  CHECK IF THE DEPLETION CHAIN IS COHERENT
*----
      DO 200 IEL=1,NEL
        DO 210 JEL=1,NEL
          KREAC=KPAX(JEL,IEL)
          IF(KREAC.EQ.KFISSP) THEN
            IF(MOD(KPAX(NEL+KFISSP,JEL),100).NE.5) THEN
              WRITE(TEXT12,'(3A4)') (ITNAM(II,JEL),II=1,3)
              CALL XABORT('LIBWET: SON '//TEXT12//' IS NOT A FISSION '//
     >                    'PRODUCT')
            ENDIF
            IF((KPAX(NEL+KFISSP,IEL).NE.3).AND.
     >         (MOD(KPAX(NEL+KFISSP,IEL),100).NE.4)) THEN
              WRITE(TEXT12,'(3A4)') (ITNAM(II,IEL),II=1,3)
              CALL XABORT('LIBWET: PARENT '//TEXT12//' IS NOT A FISSI'//
     >                    'LE ISOTOPE')
            ENDIF
          ELSE IF(KREAC.NE.0) THEN
            IF(KPAX(NEL+KREAC,IEL).EQ.0) THEN
              WRITE(TEXT12,'(3A4)') (ITNAM(II,IEL),II=1,3)
              CALL XABORT('LIBWET: PARENT '//TEXT12//' DOES NOT DEPL'//
     >                    'ETE VIA REACTION '//NMDEPL(KREAC))
            ENDIF
          ENDIF
 210    CONTINUE
 200  CONTINUE
*----
*  SET MATNO TO OTHER (-KDECAY) FOR ISOTOPES WITH DAUGHTER OR FATHER
*  AND FOR ISOTOPES WITH DECAY
*----
      DO 112 ISO=1,NEL
        IF(MATNO(ISO).EQ.0) THEN
          IF(KPAX(NEL+KDECAY,ISO).GT.0 .OR.
     >       KPAX(NEL+KCAPTU,ISO).GT.0) THEN
            MATNO(ISO)=-KDECAY
            GO TO 115
          ENDIF
          DO 113 JSO=1,NEL
            IF(KPAX(JSO,ISO).GT.0) THEN
              MATNO(ISO)=-KDECAY
              GO TO 115
            ELSE IF(KPAX(ISO,JSO).GT.0) THEN
              MATNO(ISO)=-KDECAY
              GO TO 115
            ENDIF
 113      CONTINUE
        ENDIF
 115    CONTINUE
 112  CONTINUE
*----
*  SET MATNO TO STABLE (-KHEAT) FOR OTHER ISOTOPES PRODUCING ENERGY
*----
      DO 212 ISO=1,NEL
        IF(MATNO(ISO).EQ.0) THEN
          DO 213 ITR=2,MAXR
            IF(BPAX(NEL+ITR,ISO).NE.0.0) THEN
              MATNO(ISO)=-KHEAT
              GO TO 215
            ENDIF
 213      CONTINUE
        ENDIF
 215    CONTINUE
 212  CONTINUE
*----
*  COMPUTE THE MAXIMUM NUMBER OF PARENTS FOR AN ISOTOPE
*----
      NPAR=0
      DO 116 ISO=1,NEL
        ICOUNT=0
        DO 114 JSO=1,NEL
          LOGPF=((MATNO(ISO).EQ.-KFISSP).AND.(MATNO(JSO).EQ.-KFISSI))
          IF((BPAX(ISO,JSO).NE.0.0).AND.(.NOT.LOGPF)) ICOUNT=ICOUNT+1
 114    CONTINUE
        NPAR=MAX(NPAR,ICOUNT)
 116  CONTINUE
*----
*  COMPUTE THE MAXIMUM NUMBER OF DEPLETION REACTIONS
*----
      NREAC=4
      DO 118 ISO=1,NEL
        DO 117 ITR=1,MAXR
          IF(KPAX(NEL+ITR,ISO).NE.0) NREAC=MAX(NREAC,ITR)
 117    CONTINUE
 118  CONTINUE
*----
*  SET MATNO TO HEAVY (-KFISSI) FOR ISOTOPES RESULTING FROM DECAY
*  OR CAPTURE OF HEAVY ISOTOPE UNTIL ALL HEAVY ISOTOPES IDENTIFIED
*----
      NUNTIL=NEL
      NHEAVY=0
      DO 120 IUNTIL=1,NUNTIL
        NEW=0
        DO 121 ISO=NEL,1,-1
          IF(MATNO(ISO).EQ.-KFISSI) THEN
            NHEAVY=NHEAVY+1
            DO 122 JSO=NEL,1,-1
              IF(MATNO(JSO).NE.-KFISSI) THEN
                IF((KPAX(JSO,ISO).NE.0).AND.(KPAX(JSO,ISO).NE.KFISSP))
     >          THEN
                  NEW=NEW+1
                  MATNO(JSO)=-KFISSI
                ENDIF
              ENDIF
 122        CONTINUE
          ENDIF
 121    CONTINUE
        IF(NEW.EQ.0) GO TO 125
        NHEAVY=0
 120  CONTINUE
 125  CONTINUE
*----
*  SET MATNO TO LIGHT (-KFISSP) FOR ISOTOPES RESULTING FROM DECAY
*  OR CAPTURE OF LIGHT ISOTOPE UNTIL ALL LIGHT ISOTOPES IDENTIFIED
*----
      NUNTIL=NUNTIL-NHEAVY
      NLIGHT=0
      DO 130 IUNTIL=1,NUNTIL
        NEW=0
        DO 131 ISO=NEL,1,-1
          IF(MATNO(ISO).EQ.-KFISSP) THEN
            NLIGHT=NLIGHT+1
            DO 132 JSO=NEL,1,-1
              IF(MATNO(JSO).NE.-KFISSI.AND.
     >           MATNO(JSO).NE.-KFISSP) THEN
                IF(KPAX(JSO,ISO).NE.0) THEN
                  NEW=NEW+1
                  MATNO(JSO)=-KFISSP
                ENDIF
              ENDIF
 132        CONTINUE
          ENDIF
 131    CONTINUE
        IF(NEW.EQ.0) GO TO 135
        NLIGHT=0
 130  CONTINUE
 135  CONTINUE
*----
*  SET MATNO TO OTHER (-KDECAY) FOR ISOTOPES RESULTING FROM DECAY
*  OR CAPTURE OF OTHER ISOTOPE UNTIL ALL OTHER ISOTOPES IDENTIFIED
*----
      NUNTIL=NUNTIL-NLIGHT
      NOTHER=0
      DO 140 IUNTIL=1,NUNTIL
        NEW=0
        DO 141 ISO=NEL,1,-1
          IF(MATNO(ISO).EQ.-KDECAY) THEN
            NOTHER=NOTHER+1
            DO 142 JSO=NEL,1,-1
              IF(MATNO(JSO).NE.-KFISSI.AND.
     >           MATNO(JSO).NE.-KFISSP.AND.
     >           MATNO(JSO).NE.-KDECAY) THEN
                IF(KPAX(JSO,ISO).NE.0) THEN
                  NEW=NEW+1
                  MATNO(JSO)=-KDECAY
                ENDIF
              ENDIF
 142        CONTINUE
          ENDIF
 141    CONTINUE
        IF(NEW.EQ.0) GO TO 145
        NOTHER=0
 140  CONTINUE
 145  CONTINUE
*----
*  COUNT THE NUMBER OF STABLE ISOTOPES AND SET TO ZERO THEIR RADIOACTIVE
*  DECAY CONSTANT.
*----
      NSTABL=0
      DO 150 ISO=1,NEL
        IF(MATNO(ISO).EQ.-KHEAT) THEN
           NSTABL=NSTABL+1
           BPAX(NEL+KDECAY,ISO)=0.0
        ENDIF
 150  CONTINUE
*
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NHEAVY+NLIGHT+NOTHER+NSTABL
      ISTATE(2)=NDFI
      ISTATE(3)=NDFP
      ISTATE(4)=NHEAVY
      ISTATE(5)=NLIGHT
      ISTATE(6)=NOTHER
      ISTATE(7)=NSTABL
      ISTATE(8)=NREAC
      ISTATE(9)=NPAR
*----
*  RETURN
*----
      RETURN
      END
