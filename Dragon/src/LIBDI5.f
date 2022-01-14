*DECK LIBDI5
      SUBROUTINE LIBDI5 (MAXDIL,NAMFIL,HSHI,NDIL,DILUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the dilutions corresponding to a resonant isotope within a
* library in Apolib-2 format.
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
* MAXDIL  maximum number of dilutions.
* NAMFIL  name of the APOLIB-2 file.
* HSHI    library name of the self-shielding data.
*
*Parameters: output
* NDIL    number of finite dilutions.
* DILUT   dilutions.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXDIL,NDIL
      CHARACTER HSHI*12
      CHARACTER NAMFIL*(*)
      REAL DILUT(MAXDIL)
*----
*  LOCAL VARIABLES
*----
      EXTERNAL LIBA21
      CHARACTER HSMG*131,TEXT8*8,TEXT20*20,NOMOBJ*20,TYPOBJ*8,TYPSEG*8
      LOGICAL LISO,LPTHOM
      INTEGER ISFICH(3)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: VINTE,ITCARO
      TYPE(C_PTR) ICHDIM_PTR,ICHTYP_PTR,ICHDKL_PTR,TSEGM_PTR
      INTEGER, POINTER, DIMENSION(:) :: ICHDIM,ICHTYP,ICHDKL,ITSEGM
      REAL, POINTER, DIMENSION(:) :: RTSEGM
*
      INTEGER TKCARO(31)
      SAVE TKCARO
      DATA TKCARO /
     &   0,   1,   2,   3,  4,   5,  6,  30,   7,  -8,
     &   9, -10,  11, -12, 13, -14, 15,  16, -17,  18,
     & -19,  20, -21,  22, 23, -24, 25, -26,  27, -28,
     &  29   /
*----
*  INDEX THE APOLIB-2 FILE.
*----
      CALL AEXTPA(NAMFIL,ISFICH)
      IADRES=ISFICH(1)
      NBOBJ=ISFICH(2)
      LBLOC=ISFICH(3)
      NIN=KDROPN(NAMFIL,2,4,LBLOC)
      IF(NIN.LE.0) THEN
         WRITE (HSMG,'(35HLIBDI5: UNABLE TO OPEN LIBRARY FILE,1X,A16,
     1   1H.)') NAMFIL
         CALL XABORT(HSMG)
      ENDIF
      ALLOCATE(VINTE(2*NBOBJ))
      CALL AEXDIR(NIN,LBLOC,VINTE,IADRES,2*NBOBJ)
      IDKNO=1-TKCARO(14)
      IDKTY=1-TKCARO(21)
      IDKDS=1-TKCARO(10)
      IDKTS=1-TKCARO(23)
      IDKDA=1-TKCARO(26)
      IDKNS=TKCARO(2)+1
      IDKLS=TKCARO(8)
*
      TEXT20='SSDATA'//HSHI
      LISO=.FALSE.
      DO 70 I=3,NBOBJ
      IDKOBJ=VINTE(2*I-1)
      LGSEG=VINTE(2*I)+1
      ALLOCATE(ITCARO(LGSEG))
      CALL AEXDIR(NIN,LBLOC,ITCARO,IDKOBJ,LGSEG)
      IDK=ITCARO(IDKNO)
      CALL AEXCPC(IDK,20,ITCARO,NOMOBJ)
      IDK=ITCARO(IDKTY)
      CALL AEXCPC(IDK,8,ITCARO,TYPOBJ)
      JDKDS=ITCARO(IDKDS)
      JDKTS=ITCARO(IDKTS)
      NS=ITCARO(IDKNS)
      IDK=ITCARO(IDKDA)
      CALL AEXCPC(IDK,8,ITCARO,TEXT8)
      IF((TYPOBJ.EQ.'APOLIBE').AND.(NOMOBJ.EQ.TEXT20)) THEN
         LISO=.TRUE.
         LPTHOM=.FALSE.
         DO 60 IS=1,NS
         IDK=JDKTS+8*(IS-1)
         CALL AEXCPC(IDK,8,ITCARO,TYPSEG)
         LNGS=ITCARO(IDKLS+IS)
         IF(LNGS.LE.0) GO TO 60
         JDKS=ITCARO(JDKDS+IS)
         CALL AEXTRT(LIBA21,TYPSEG,NBRTYP,ICHDIM_PTR,ICHTYP_PTR,
     1   ICHDKL_PTR)
         CALL C_F_POINTER(ICHDIM_PTR,ICHDIM,(/ NBRTYP /))
         CALL C_F_POINTER(ICHTYP_PTR,ICHTYP,(/ NBRTYP /))
         CALL C_F_POINTER(ICHDKL_PTR,ICHDKL,(/ NBRTYP /))
         TSEGM_PTR=LCMARA(LNGS+1)
         CALL C_F_POINTER(TSEGM_PTR,ITSEGM,(/ LNGS+1 /))
         CALL AEXDIR(NIN,LBLOC,ITSEGM,JDKS,LNGS+1)
         IF(TYPSEG.EQ.'PTHOM1') THEN
            CALL C_F_POINTER(TSEGM_PTR,RTSEGM,(/ LNGS+1 /))
            LPTHOM=.TRUE.
            CALL AEXGNV(16,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NSEQHO)
            NDIL=NSEQHO-1
            IF(NDIL.GT.MAXDIL) CALL XABORT('LIBDI5: INVALID MAXDIL.')
            DILMAX=RTSEGM(IDK+NSEQHO-1)
            IF(DILMAX.LT.1.0E10) THEN
               WRITE(HSMG,'(35HLIBDI5: INVALID INFINITE DILUTION (,1P,
     1         E12.4,14H) FOR ISOTOPE ,A12,1H.)') DILMAX,HSHI
               CALL XABORT(HSMG)
            ENDIF
            DO 50 J=1,NSEQHO
            DILUT(J)=RTSEGM(IDK+J-1)
   50       CONTINUE
         ENDIF
         CALL LCMDRD(TSEGM_PTR)
         CALL LCMDRD(ICHDIM_PTR)
         CALL LCMDRD(ICHTYP_PTR)
         CALL LCMDRD(ICHDKL_PTR)
   60    CONTINUE
         IF(.NOT.LPTHOM) CALL XABORT('LIBDI5: NO PTHOM1 SEGMENT '
     1   //'FOR ISOTOPE '//HSHI//'.')
      ENDIF
      DEALLOCATE(ITCARO)
   70 CONTINUE
      DEALLOCATE(VINTE)
      IF(.NOT.LISO) CALL XABORT('LIBDI5: UNABLE TO FIND ISOTOPE '
     1 //HSHI//'.')
      IER=KDRCLS(NIN,1)
      IF(IER.LT.0) THEN
         WRITE (HSMG,'(36HLIBDI5: UNABLE TO CLOSE LIBRARY FILE,1X,A16,
     1   1H.)') NAMFIL
         CALL XABORT(HSMG)
      ENDIF
      RETURN
      END
