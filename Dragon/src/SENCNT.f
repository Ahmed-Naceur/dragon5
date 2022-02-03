*DECK SENCNT
      SUBROUTINE SENCNT(IPLIB,NI,NAMISO,MELISO,NSENS,NSENI,NAMISC,ISOC,
     > NIC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Count the number of sensitivity coefficients.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Laville, G. Marleau
*
*Parameters: input/output
* IPLIB   LCM Library object address.
* NI      number of isotope/mixture.
* NAMISO  name of the isotope/mixture.
* MELISO  mixture associated with the isotope/mixture.
* NSENS   number of sensitivity profiles. 
* NSENI   number of integrated sensitivity profiles. 
* NAMISC  independent isotopes names.
* ISOC    independent isotope number associated with isotope/mixture.
* NIC     number of independent isotopes names.
* 
*-----------------------------------------------------------------------
* 
      USE GANLIB
*----
*  Suboutine arguements
*----
      IMPLICIT      NONE
      TYPE(C_PTR)   IPLIB
      INTEGER       NI
      INTEGER       NAMISO(3,NI),MELISO(NI)
      INTEGER       NSENS,NSENI
      INTEGER       NAMISC(2,NI),ISOC(NI),NIC
*----
*  Parameters
*---- 
      INTEGER       IOUT
      CHARACTER     NAMSBR*6
      PARAMETER    (IOUT=6,NAMSBR='SENCNT')
*----
*  Local variables
*----
      TYPE(C_PTR)   KPISO
      CHARACTER     ISONAM*12,HSMG*131
      INTEGER       ILENG,ITYLCM,II,IJ,ISOMEL
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  Initialize the directory of the isotope/mixture
*  in the library and initialize the SEN information
*  support and the cross section support for the isotope
*----
      NIC=0
      NSENS=0
      NSENI=0
      CALL XDISET(ISOC,NI,0)
*----
*  Find independent isotopes names NAMISC
*  Associate isotope in mixture to 
*  independant isotopes names ISOC
*----
      DO II=1,NI
        IF(ISOC(II) .EQ. 0) THEN
          NIC=NIC+1
          NAMISC(1,NIC)=NAMISO(1,II)
          NAMISC(2,NIC)=NAMISO(2,II)
          ISOC(II)=NIC
          DO IJ=II+1,NI
            IF(NAMISC(1,NIC) .EQ. NAMISO(1,IJ) .AND.
     >         NAMISC(2,NIC) .EQ. NAMISO(2,IJ)) THEN
              ISOC(IJ)=-NIC
            ENDIF
          ENDDO
        ENDIF
      ENDDO
*
      ALLOCATE(IPISO(NI))
      CALL LIBIPS(IPLIB,NI,IPISO)
      DO II=1,NI
        KPISO=IPISO(II) ! set II-th isotope
        IF(.NOT.C_ASSOCIATED(KPISO)) THEN
          WRITE(ISONAM,'(3A4)') NAMISO(1,II),NAMISO(2,II),NAMISO(3,II)
          WRITE(HSMG,'(17HSENCNT: ISOTOPE '',A12,7H'' (ISO=,I8,5H) IS ,
     1    30HNOT AVAILABLE IN THE MICROLIB.)') ISONAM,II
          CALL XABORT(HSMG)
        ENDIF
        ISOMEL=MELISO(II)
*----
*  number of (n,g) sensitivity
*----
        CALL LCMLEN(KPISO,'NG',ILENG,ITYLCM)
        IF(ILENG.GT.0) THEN
          NSENS=NSENS+1
          IF(ISOC(II) .GT. 0) NSENI=NSENI+1
        ENDIF      
*----
*  number of (n,p) sensitivity
*----
        CALL LCMLEN(KPISO,'NP',ILENG,ITYLCM)
        IF(ILENG.GT.0) THEN
          NSENS=NSENS+1
          IF(ISOC(II) .GT. 0) NSENI=NSENI+1
        ENDIF      
*----
*  number of (n,d) sensitivity
*----
        CALL LCMLEN(KPISO,'ND',ILENG,ITYLCM)
        IF(ILENG.GT.0) THEN
          NSENS=NSENS+1
          IF(ISOC(II) .GT. 0) NSENI=NSENI+1
        ENDIF
*----
*  number of (n,a) sensitivity
*----
        CALL LCMLEN(KPISO,'NA',ILENG,ITYLCM)
        IF(ILENG.GT.0) THEN
          NSENS=NSENS+1
          IF(ISOC(II) .GT. 0) NSENI=NSENI+1
        ENDIF
*----
*  number of Capture sensitivity
*----
        NSENS=NSENS+1
        IF(ISOC(II) .GT. 0) NSENI=NSENI+1
*----
*  number of Scattering sensitivity
*----
        NSENS=NSENS+1
        IF(ISOC(II) .GT. 0) NSENI=NSENI+1
*----
*  number of Fissile related sensitivity
*----     
        CALL LCMLEN(KPISO,'NUSIGF',ILENG,ITYLCM)
        IF(ILENG.GT.0) THEN
          NSENS=NSENS+3
          IF(ISOC(II) .GT. 0) NSENI=NSENI+3
        ENDIF
      ENDDO
      DEALLOCATE(IPISO)
      RETURN
      END
