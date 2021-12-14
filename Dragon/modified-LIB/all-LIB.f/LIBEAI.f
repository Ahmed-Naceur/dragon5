*DECK LIBEAI
      SUBROUTINE LIBEAI(CFILNA,NEL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Initialize dimensions for depletion data with APOLIB-2.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input
* CFILNA  APOLIB-2 file name.
*
*Parameters: output
* NEL     number of isotopes on library.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CFILNA*(*)
      INTEGER NEL
*
      EXTERNAL LIBA21
      CHARACTER TYPOBJ*8,TYPSEG*8,TEXT8*8
      INTEGER ISFICH(3)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: VINTE,ITCARO,ITSEGM
      TYPE(C_PTR) ICHDIM_PTR,ICHTYP_PTR,ICHDKL_PTR
      INTEGER, POINTER, DIMENSION(:) :: ICHDIM,ICHTYP,ICHDKL
*
      INTEGER TKCARO(31)
      SAVE TKCARO
      DATA TKCARO /
     &   0,   1,   2,   3,  4,   5,  6,  30,   7,  -8,
     &   9, -10,  11, -12, 13, -14, 15,  16, -17,  18,
     & -19,  20, -21,  22, 23, -24, 25, -26,  27, -28,
     &  29   /
*----
*  PROBE AND OPEN THE APOLIB-2 FILE.
*----
      CALL AEXTPA(CFILNA,ISFICH)
      IADRES=ISFICH(1)
      NBOBJ=ISFICH(2)
      LBLOC=ISFICH(3)
      IUNIT=KDROPN(CFILNA,2,4,LBLOC)
      IF(IUNIT.LE.0) THEN
         TEXT8=CFILNA
         CALL XABORT('LIBEAI: APOLLO-2 LIBRARY '//TEXT8//' CANNOT BE'//
     >   ' OPENED')
      ENDIF
*----
*  INDEX THE APOLIB-2 FILE.
*----
      ALLOCATE(VINTE(2*NBOBJ))
      CALL AEXDIR(IUNIT,LBLOC,VINTE,IADRES,2*NBOBJ)
      IDKNO=1-TKCARO(14)
      IDKTY=1-TKCARO(21)
      IDKDS=1-TKCARO(10)
      IDKTS=1-TKCARO(23)
      IDKNS=TKCARO(2)+1
      IDKLS=TKCARO(8)
*
      DO 70 I=3,NBOBJ
      IDKOBJ=VINTE(2*I-1)
      LGSEG=VINTE(2*I)+1
      ALLOCATE(ITCARO(LGSEG))
      CALL AEXDIR(IUNIT,LBLOC,ITCARO,IDKOBJ,LGSEG)
      IDK=ITCARO(IDKTY)
      CALL AEXCPC(IDK,8,ITCARO,TYPOBJ)
      JDKDS=ITCARO(IDKDS)
      JDKTS=ITCARO(IDKTS)
      NS=ITCARO(IDKNS)
      IF(TYPOBJ.EQ.'APOLIB') THEN
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
         ALLOCATE(ITSEGM(LNGS+1))
         CALL AEXDIR(IUNIT,LBLOC,ITSEGM,JDKS,LNGS+1)
         IF(TYPSEG.EQ.'PHEAD') THEN
            CALL AEXGNV(3,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NV)
            IF(NV.EQ.0) THEN
               TEXT8=CFILNA
               CALL XABORT('LIBEAI: NO ISOTOPES PRESENT ON APOLIB-2 '//
     1         'FILE NAMED: '//TEXT8)
            ENDIF
            NEL=NV/20
         ENDIF
         DEALLOCATE(ITSEGM)
         CALL LCMDRD(ICHDIM_PTR)
         CALL LCMDRD(ICHTYP_PTR)
         CALL LCMDRD(ICHDKL_PTR)
   60    CONTINUE
      ENDIF
      DEALLOCATE(ITCARO)
   70 CONTINUE
      DEALLOCATE(VINTE)
      IERR=KDRCLS(IUNIT,1)
      IF(IERR.LT.0) THEN
         TEXT8=CFILNA
         CALL XABORT('LIBEAI: APOLLO-2 LIBRARY '//TEXT8//' CANNOT BE'//
     >   ' CLOSED')
      ENDIF
      RETURN
      END
