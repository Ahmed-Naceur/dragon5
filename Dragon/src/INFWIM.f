*DECK INFWIM
      SUBROUTINE INFWIM(CFILNA,IPRINT,NBISO,HNAMIS,AWR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover mass for isotopes of WIMS-AECL libraries.
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
* CFILNA  WIMS file name.
* IPRINT  print flag.
* NBISO   number of isotopes.
* HNAMIS  isotope names.
*
*Parameters: output
* AWR     isotope weights
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER    IOUT,IUTYPE,IACTO,LRIND,MAXISO,MAXTEM,NCT,
     >           LPZ,LMASTB,LMASIN,LGENTB,LGENIN,LSUBTB,LSUBIN
      PARAMETER (IOUT=6,IUTYPE=4,IACTO=2,LRIND=256,
     >           MAXISO=246,MAXTEM=20,
     >           NCT=10,LPZ=9,LMASTB=MAXISO+9,
     >           LMASIN=LMASTB-4,LGENTB=6,LGENIN=LGENTB,
     >           LSUBTB=6*MAXTEM+28,LSUBIN=LSUBTB-12)
      INTEGER    MASTER(LMASTB),GENINX(LGENTB),SUBINX(LSUBTB),
     >           IWISO(MAXISO),IPRINT,NBISO,NPZ(LPZ),
     >           ITITLE(2*NCT)
      CHARACTER  CFILNA*64,HNAMIS(NBISO)*8,
     >           CWISO(MAXISO)*8,CTITLE(NCT)*8
      REAL       ZUBINX(LSUBTB),AWR(NBISO)
      INTEGER    IUNIT,IRISO,ISO,JSO,II,KDROPN
      EXTERNAL   KDROPN
      EQUIVALENCE (SUBINX(1),ZUBINX(1))
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ICISO
*
      IF( CFILNA.EQ.' ' )THEN
        CALL XABORT('INFWIM: WIMS LIBRARY HAS NOT BEEN SET')
      ENDIF
*----
*  OPEN WIMSLIB AND READ TITLE
*----
      IRISO=0
      IUNIT=KDROPN(CFILNA,IACTO,IUTYPE,LRIND)
      IF(IUNIT.LE.0) CALL XABORT(
     > 'INFWIM: WIMS-AECL LIBRARY CANNOT BE OPENED FOR MIXS :'//CFILNA)
      CALL OPNIND(IUNIT,MASTER,LMASTB)
      CALL REDIND(IUNIT,MASTER,LMASIN,GENINX,LGENTB,1)
      IF(IPRINT.GT.0) THEN
        CALL REDIND(IUNIT,MASTER,LMASIN,ITITLE,2*NCT,2)
        CALL UPCKIC(ITITLE,CTITLE,NCT)
        WRITE(IOUT,6000) CFILNA
        WRITE(IOUT,'(1X,10A8)') (CTITLE(II),II=1,NCT)
      ENDIF
*----
*  READ GENERAL INDEX, ISOTOPES NAMES AND GROUP STRUCTURE
*----
      CALL REDIND(IUNIT,GENINX,LGENIN,NPZ,LPZ,1)
      CALL REDIND(IUNIT,GENINX,LGENIN,IWISO,NPZ(1),2)
      ALLOCATE(ICISO(2*NPZ(1)))
      CALL REDIND(IUNIT,GENINX,LGENIN,ICISO,2*NPZ(1),3)
      CALL UPCKIC(ICISO,CWISO,NPZ(1))
      IF(IPRINT.GE.100) THEN
        WRITE(IOUT,6200) (CWISO(II),II=1,NPZ(1))
      ENDIF
      DEALLOCATE(ICISO)
*----
* READ THROUGH DRAGON FILE AND ACCUMULATE WEIGHTS.
*----
      DO 120 ISO=1,NBISO
        DO 130 JSO=1,NPZ(1)
          IF(CWISO(JSO).EQ.HNAMIS(ISO)) THEN
            IRISO=JSO
            GO TO 131
          ENDIF
 130    CONTINUE
        WRITE(IOUT,9002) HNAMIS(ISO),CFILNA
        CALL XABORT('INFWIM: ISOTOPE NOT FOUND ON LIBRARY')
 131    CONTINUE
        IF(IPRINT.GT.0) THEN
          WRITE(IOUT,6001) HNAMIS(ISO)
        ENDIF
*----
*  READ SUB INDEX ASSOCIATED WITH ISOTOPE
*----
        CALL REDIND(IUNIT,MASTER,LMASIN,SUBINX,LSUBTB,IRISO+4)
        AWR(ISO)=ZUBINX(LSUBIN+3)
 120  CONTINUE
      CALL CLSIND(IUNIT)
*----
*  RETURN
*----
      RETURN
*----
*  FORMAT
*----
 9002 FORMAT(/' INFWIM: MATERIAL/ISOTOPE ',A8,' IS MISSING ON WIMS',
     >        ' FILE NAME ',A64)
 6000 FORMAT(/' PROCESSING WIMS LIBRARY NAME ',A64)
 6001 FORMAT(/'    PROCESSING ISOTOPE/MATERIAL = ',A12)
 6200 FORMAT(1X,'ISOTOPES ON LIBRARY'/6(4X,A8))
      END
