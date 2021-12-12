*DECK LIBA2G
      SUBROUTINE LIBA2G (NAMFIL,NGRO,IPENER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover energy group information from an APOLIB2 library.
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
* NAMFIL  name of the APOLIB2 file.
*
*Parameters: output
* NGRO    number of energy groups.
* IPENER  pointer of the energy mesh limit array.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  Subroutine arguments
*----
      INTEGER NGRO
      CHARACTER NAMFIL*(*)
      TYPE(C_PTR) IPENER
*----
*  Local variables
*----
      PARAMETER (IACTO=2,IACTC=1,ILIBDA=4)
      EXTERNAL LIBA21
      INTEGER ISFICH(3)
      CHARACTER TEXT80*80,NOMOBJ*20,TYPOBJ*8,TYPSEG*8,HSMG*131
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: VINTE,ITCARO
      INTEGER, POINTER, DIMENSION(:) :: ITSEGM
      REAL, POINTER, DIMENSION(:) :: RTSEGM
      TYPE(C_PTR) ICHDIM_PTR,ICHTYP_PTR,ICHDKL_PTR,TSEGM_PTR
      INTEGER, POINTER, DIMENSION(:) :: ICHDIM,ICHTYP,ICHDKL
      REAL, POINTER, DIMENSION(:) :: ENERG
*
      INTEGER TKCARO(31)
      SAVE TKCARO
      DATA TKCARO /
     &   0,   1,   2,   3,  4,   5,  6,  30,   7,  -8,
     &   9, -10,  11, -12, 13, -14, 15,  16, -17,  18,
     & -19,  20, -21,  22, 23, -24, 25, -26,  27, -28,
     &  29   /
*
      CALL AEXTPA(NAMFIL,ISFICH)
      IADRES=ISFICH(1)
      NBOBJ=ISFICH(2)
      LBLOC=ISFICH(3)
      IUNIT=KDROPN(NAMFIL,IACTO,ILIBDA,LBLOC)
      IF(IUNIT.LE.0) THEN
         WRITE(HSMG,'(26HLIBA2G: APOLLO-2 LIBRARY '',A16,9H'' CANNOT ,
     >   29HBE OPENED BY KDROPN (ERRCODE=,I2,2H).)') NAMFIL,IUNIT
         CALL XABORT(HSMG)
      ENDIF
      ALLOCATE(VINTE(2*NBOBJ))
      CALL AEXDIR(IUNIT,LBLOC,VINTE,IADRES,2*NBOBJ)
      IDKSV=1-TKCARO(12)
      IDKNO=1-TKCARO(14)
      IDKTY=1-TKCARO(21)
      IDKDS=1-TKCARO(10)
      IDKTS=1-TKCARO(23)
      IDKNS=TKCARO(2)+1
      IDKLS=TKCARO(8)
      DO 150 IOBJ=3,NBOBJ
         IDKOBJ=VINTE(2*IOBJ-1)
         LGSEG=VINTE(2*IOBJ)+1
         ALLOCATE(ITCARO(LGSEG))
         CALL AEXDIR(IUNIT,LBLOC,ITCARO,IDKOBJ,LGSEG)
         IDK=ITCARO(IDKSV)
         CALL AEXCPC(IDK,80,ITCARO,TEXT80)
*
         IDK=ITCARO(IDKNO)
         CALL AEXCPC(IDK,20,ITCARO,NOMOBJ)
         IDK=ITCARO(IDKTY)
         CALL AEXCPC(IDK,8,ITCARO,TYPOBJ)
         IF(TYPOBJ.EQ.'APOLIB') THEN
            JDKDS=ITCARO(IDKDS)
            JDKTS=ITCARO(IDKTS)
            NS=ITCARO(IDKNS)
            DO 140 IS=1,NS
              IDK=JDKTS+8*(IS-1)
              CALL AEXCPC(IDK,8,ITCARO,TYPSEG)
              IF(TYPSEG.EQ.'PMAIL') THEN
                LNGS=ITCARO(IDKLS+IS)
                JDKS=ITCARO(JDKDS+IS)
                TSEGM_PTR=LCMARA(LNGS+1)
                CALL C_F_POINTER(TSEGM_PTR,ITSEGM,(/ LNGS+1 /))
                CALL C_F_POINTER(TSEGM_PTR,RTSEGM,(/ LNGS+1 /))
                CALL AEXDIR(IUNIT,LBLOC,ITSEGM,JDKS,LNGS+1)
                CALL AEXTRT(LIBA21,TYPSEG,NBRTYP,ICHDIM_PTR,ICHTYP_PTR,
     1          ICHDKL_PTR)
                CALL C_F_POINTER(ICHDIM_PTR,ICHDIM,(/ NBRTYP /))
                CALL C_F_POINTER(ICHTYP_PTR,ICHTYP,(/ NBRTYP /))
                CALL C_F_POINTER(ICHDKL_PTR,ICHDKL,(/ NBRTYP /))
                CALL AEXGNV(3,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NV)
                NGRO=NV-1
                IPENER=LCMARA(NGRO+1)
                CALL C_F_POINTER(IPENER,ENERG,(/ NGRO+1 /))
                DO 130 IG=1,NV
                ENERG(IG)=RTSEGM(IDK+IG-1)*1.0E6
  130           CONTINUE
                CALL LCMDRD(ICHDIM_PTR)
                CALL LCMDRD(ICHTYP_PTR)
                CALL LCMDRD(ICHDKL_PTR)
                CALL LCMDRD(TSEGM_PTR)
                DEALLOCATE(ITCARO)
                GO TO 160
              ENDIF
  140       CONTINUE
         ENDIF
         DEALLOCATE(ITCARO)
  150 CONTINUE
      CALL XABORT('LIBA2G: NO GROUP STRUCTURE AVAILABLE')
*
  160 IERR=KDRCLS(IUNIT,IACTC)
      IF(IERR.LT.0) THEN
         WRITE(HSMG,'(26HLIBA2G: APOLLO-2 LIBRARY '',A16,9H'' CANNOT ,
     >   29HBE CLOSED BY KDRCLS (ERRCODE=,I2,2H).)') NAMFIL,IERR
         CALL XABORT(HSMG)
      ENDIF
      DEALLOCATE(VINTE)
      RETURN
      END
