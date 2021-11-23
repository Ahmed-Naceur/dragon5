*DECK LIBEWR
      SUBROUTINE LIBEWR(CFILNA,MAXR,NEL,ITNAM,KPAX,BPAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read depletion data on a WIMS-AECL formatted library.
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
* CFILNA  file name.
* MAXR    number of reaction types.
* NEL     number of isotopes on library.
*
*Parameters: output
* ITNAM   reactive isotope names in chain.
* KPAX    complete reaction type matrix.
* BPAX    complete branching ratio matrix.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CFILNA*8
      INTEGER MAXR,NEL,ITNAM(3,NEL),KPAX(NEL+MAXR,NEL)
      REAL BPAX(NEL+MAXR,NEL)
*----
*  INTERNAL PARAMETERS
*   CONVE  : ENERGY CONVERSION FACTOR FROM JOULES/(MOLES*10**-24)
*            TO MEV/NUCLIDE = 1.03643526E+13
*   CONVD  : DECAY CONSTANT CONVERSION FACTOR FROM S**(-1) TO
*            10**(-8)*S**(-1) = 1.0+8
*----
      INTEGER      KCAPTU,KDECAY,KFISSP,KN2N,KN3N
      REAL         CONVE,CONVD
      PARAMETER   (KDECAY=1,KFISSP=2,KCAPTU=3,KN2N=4,KN3N=5,
     >             CONVE=1.03643526E+13,CONVD=1.0E+8)
      CHARACTER    TEXT8*8
*----
*  WIMS-AECL LIBRARY PARAMETERS
*   IUTYPE : TYPE OF FILE = 4 (DA)                I
*   LRIND  : LENGHT RECORD ON DA FILE = 256       I
*   IACTO  : OPEN ACTION = 2 (READ ONLY)          I
*   IACTC  : CLOSE ACTION = 2 (KEEP)              I
*   MAXISO : MAX. NB. OF ISO = 246                I
*   MLDEP  : MAXIMUM NUMBER OF REACTION PER       I
*            ISOTOPE IN WIMS-AECL = MAXISO+4
*   LPZ    : LENGTH OF WIMS PARAMETER ARRAY = 9   I
*   LMASTB : LENGTH OF MST TAB = MAXISO+9         I
*   LMASIN : LENGTH OF MST IDX = LMASTB-4         I
*   LGENTB : LENGTH OF GEN TAB = 6                I
*   LGENIN : LENGTH OF GEN IDX = LGENTB           I
*   LSUBTB : LENGTH OF SUB TAB = 6*MAXTEM+21-5+12 I
*   LSUBIN : LENGTH OF SUB IDX = LSUBTB-12        I
*   ICAPTU : WIMS-AECL CAPTURE FLAG = 1           I
*   IDECAY : WIMS-AECL DECAY FLAG = 2             I
*   IFISSP : WIMS-AECL FISSION PRODUCT FLAG = 3   I
*   IFISSI : WIMS-AECL FISSILE ISOTOPE FLAG = 4   I
*   IN2N   : WIMS-AECL N2N FLAG = 5               I
*   IN3N   : WIMS-AECL N3N FLAG = 6               I
*   MASTER : MASTER INDEX ARRAY                   I(LMASTB)
*   GENINX : GENERAL INDEX ARRAY                  I(LGENTB)
*   SUBINX : SUB INDEX ARRAY                      I(LSUBTB)
*   NPZ    : LIST OF MAIN PARAMETERS              I(LPZ)
*   IWISO  : ID OF ISOTOPE                        I(2*MAXISO)
*   IBURN  : INTEGER BURNUP PARAMETERS            I(2,MLDEP)
*   RBURN  : REAL BURNUP PARAMETERS               R(2,MLDEP)
*----
      INTEGER      IUTYPE,LRIND,IACTO,IACTC,MAXISO,MLDEP,LPZ,
     >             MAXTEM,LMASTB,LMASIN,LGENTB,LGENIN,LSUBTB,
     >             LSUBIN,ICAPTU,IDECAY,IFISSP,IFISSI,IN2N,IN3N
      PARAMETER   (IUTYPE=4,LRIND=256,IACTO=2,IACTC=1,MAXISO=246,
     >             MLDEP=MAXISO+4,LPZ=9,MAXTEM=20,LMASTB=MAXISO+9,
     >             LMASIN=LMASTB-4,LGENTB=6,LGENIN=LGENTB,
     >             LSUBTB=6*MAXTEM+28,LSUBIN=LSUBTB-12,ICAPTU=1,
     >             IDECAY=2,IFISSP=3,IFISSI=4,IN2N=5,IN3N=6)
      INTEGER      MASTER(LMASTB),GENINX(LGENTB),SUBINX(LSUBTB),
     >             NPZ(LPZ),IWISO(2*MAXISO)
*----
*  EXTERNAL FUNCTIONS
*----
      INTEGER      KDROPN,LIBWID
*----
*  LOCAL VARIABLES
*----
      INTEGER      IUNIT,IEL2,IEL,ISO,NBURN,NMIN,NFP,JBRN,JSO
      INTEGER      NDECAY,IBURN(2*MLDEP)
      DOUBLE PRECISION TOTLAM
      REAL         RBURN(2*MLDEP)
*----
*  OPEN WIMS-AECL LIBRARY
*  READ INDEX AND GENERAL DIMENSIONING NPZ
*  READ ISOTOPE NAME AND ID NUMBER
*----
      IUNIT=KDROPN(CFILNA,IACTO,IUTYPE,LRIND)
      IF(IUNIT.LE.0) CALL XABORT('LIBEWR: WIMS-AECL LIBRARY '//
     >    CFILNA//' CANNOT BE OPENED FOR DEPLETION')
      CALL OPNIND(IUNIT,MASTER,LMASTB)
      CALL REDIND(IUNIT,MASTER,LMASIN,GENINX,LGENTB,1)
      CALL REDIND(IUNIT,GENINX,LGENIN,NPZ,LPZ,1)
      IF(NPZ(1).NE.NEL) CALL XABORT('LIBEWR: TOO MANY ISOTOPES '//
     >    'ON WIMS-AECL LIBRARY'//CFILNA)
      CALL REDIND(IUNIT,GENINX,LGENIN,IWISO,2*NEL,3)
      IEL2=1
      DO 10 IEL=1,NEL
        CALL UPCKIC(IWISO(IEL2),TEXT8,1)
        READ(TEXT8,'(2A4)') ITNAM(1,IEL),ITNAM(2,IEL)
        IEL2=IEL2+2
 10   CONTINUE
      CALL REDIND(IUNIT,GENINX,LGENIN,IWISO,NEL,2)
*----
*  READ DEPLETION CHAIN FOR EACH ISOTOPES
*----
      DO 100 ISO=1,NEL
*----
*  READ SUB INDEX ASSOCIATED WITH ISOTOPE
*----
        NDECAY=0
        TOTLAM=0.0D0
        CALL REDIND(IUNIT,MASTER,LMASIN,SUBINX,LSUBTB,ISO+4)
        NBURN=SUBINX(LSUBIN+1)
        IF(NBURN.GT.MAXISO) THEN
          CALL XABORT('LIBEWR: NBURN LARGER THAN MAXISO')
        ENDIF
        NMIN=2*MAX0(NBURN,1)
        CALL REDIND(IUNIT,SUBINX,LSUBIN,GENINX,4,1)
        CALL REDIND(IUNIT,GENINX,4,IBURN,NMIN,1)
        CALL REDIND(IUNIT,GENINX,4,RBURN,NMIN,2)
*----
*  STORE REACTION TYPES AND RATES IN KPAX AND BPAX STARTING
*  WITH HEAVIER ISOTOPES
*----
        NFP=0
        DO 101 JBRN=1,NBURN
          IF(IBURN(2*(JBRN-1)+2).EQ.IDECAY.AND.
     >            RBURN(2*(JBRN-1)+1).GE.0.0) THEN
            JSO=LIBWID(NEL,IWISO,IBURN(2*(JBRN-1)+1))
            NDECAY=NDECAY+1
            TOTLAM=TOTLAM+DBLE(RBURN(2*(JBRN-1)+1))
            IF(JSO.GT.0) THEN
              KPAX(JSO,ISO)=KDECAY
              BPAX(JSO,ISO)=RBURN(2*(JBRN-1)+1)
              KPAX(NEL+KCAPTU,JSO)=1
            ENDIF
            KPAX(NEL+KDECAY,ISO)=1
          ELSE IF(IBURN(2*(JBRN-1)+2).EQ.IFISSI.AND.
     >            RBURN(2*(JBRN-1)+1).GE.0.0) THEN
            KPAX(NEL+KFISSP,ISO)=1
            BPAX(NEL+KFISSP,ISO)=RBURN(2*(JBRN-1)+1)*CONVE
          ELSE IF(IBURN(2*(JBRN-1)+2).EQ.ICAPTU) THEN
            JSO=LIBWID(NEL,IWISO,IBURN(2*(JBRN-1)+1))
            IF(JSO.GT.0) THEN
              KPAX(JSO,ISO)=KCAPTU
              BPAX(JSO,ISO)=RBURN(2*(JBRN-1)+1)
              KPAX(NEL+KCAPTU,JSO)=1
            ENDIF
            KPAX(NEL+KCAPTU,ISO)=1
          ELSE IF(IBURN(2*(JBRN-1)+2).EQ.IN2N) THEN
            JSO=LIBWID(NEL,IWISO,IBURN(2*(JBRN-1)+1))
            IF(JSO.GT.0) THEN
              KPAX(JSO,ISO)=KN2N
              BPAX(JSO,ISO)=RBURN(2*(JBRN-1)+1)
              KPAX(NEL+KCAPTU,JSO)=1
            ENDIF
            KPAX(NEL+KN2N,ISO)=1
          ELSE IF(IBURN(2*(JBRN-1)+2).EQ.IN3N) THEN
            JSO=LIBWID(NEL,IWISO,IBURN(2*(JBRN-1)+1))
            IF(JSO.GT.0) THEN
              KPAX(JSO,ISO)=KN3N
              BPAX(JSO,ISO)=RBURN(2*(JBRN-1)+1)
              KPAX(NEL+KCAPTU,JSO)=1
            ENDIF
            KPAX(NEL+KN3N,ISO)=1
          ELSE IF(IBURN(2*(JBRN-1)+2).EQ.IFISSP.AND.
     >            RBURN(2*(JBRN-1)+1).GE.0.0) THEN
            JSO=LIBWID(NEL,IWISO,IBURN(2*(JBRN-1)+1))
            IF(JSO.GT.0) THEN
              NFP=NFP+1
              KPAX(JSO,ISO)=KFISSP
              BPAX(JSO,ISO)=RBURN(2*(JBRN-1)+1)
              KPAX(NEL+KFISSP,JSO)=-1
              KPAX(NEL+KCAPTU,JSO)=1
            ENDIF
          ENDIF
 101    CONTINUE
        IF(NDECAY .EQ. 1) THEN
          BPAX(NEL+KDECAY,ISO)=REAL(TOTLAM)*CONVD
          DO JSO=1,NEL
            IF(KPAX(JSO,ISO).EQ. KDECAY) THEN
              BPAX(JSO,ISO)=1.0
            ENDIF
          ENDDO
        ELSE IF(NDECAY .GT. 1) THEN
          BPAX(NEL+KDECAY,ISO)=REAL(TOTLAM)*CONVD
          DO JSO=1,NEL
            IF(KPAX(JSO,ISO).EQ. KDECAY) THEN
              BPAX(JSO,ISO)=BPAX(JSO,ISO)/REAL(TOTLAM)
            ENDIF
          ENDDO
        ENDIF
 100  CONTINUE
*----
*  CLOSE WIMS-AECL LIBRARY
*----
      CALL CLSIND(IUNIT)
*----
*  RETURN
*----
      RETURN
      END
