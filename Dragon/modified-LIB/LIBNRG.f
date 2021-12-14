*DECK LIBNRG
      SUBROUTINE LIBNRG(IPLIB,NAMLBT,NAMFIL,NGROUP,NGT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Test for energy mesh compatibility.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input/output
* IPLIB   pointer to the internal library.
* NAMLBT  library type.
* NAMFIL  library file name.
* NGROUP  total number of groups.
* NGT     number of groups to test.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE LIBEEDR
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPLIB
      INTEGER          NGROUP,NGT
      CHARACTER        NAMLBT*8,NAMFIL*(*)
*----
*  LOCAL VARIABLES
*----
      PARAMETER       (IOUT=6,LRIND=256,IACTO=2,IACTC=1,ILIBDA=4,
     >                 NMTYP=8)
      CHARACTER        HSMG*131,HMTYP(NMTYP)*1
*----
*  LIBRARY PARAMETERS
*----
      PARAMETER       (MAXISO=246,NCT=10,LPZ=9,LMASTB=MAXISO+9,
     >                 LMASIN=LMASTB-4,LGENTB=6,LGENIN=LGENTB,
     >                 MAXA=10000,MULT=2)
      TYPE(C_PTR)      IPDRL
      CHARACTER        HPRT*6,NAMLCM*12,NAMMY*12
      LOGICAL          EMPTY,LCM
      INTEGER          ILONG,MASTER(LMASTB),GENINX(LGENTB),NPZ(LPZ),
     >                 IA(MAXA)
      REAL             RA(MAXA)
      DOUBLE PRECISION DA(MAXA/2)
      EQUIVALENCE     (RA(1),IA(1),DA(1))
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, POINTER, DIMENSION(:) :: ENERG
      REAL, ALLOCATABLE, DIMENSION(:) :: UUU,TEMPE,EIER
      TYPE(C_PTR) IPENER
*----
*  DATA STATEMENTS
*----
      SAVE HMTYP
      DATA HMTYP /'N','n','G','g','B','b','C','c'/
*
      NGRI=0
      ILIBIN=2
      !AHMED: IPLIB  = 17971392
      !AHMED: NAMLBT = MATXS2
      !AHMED: NAMFIL =  _test_nc
      !AHMED: NGROUP =  0
      !AHMED: NGT    =  0
*----
*  READ LIBRARY GROUP STRUCTURE
*----
      IF((NAMLBT.EQ.'DRAGON').OR.(NAMLBT.EQ.'MICROLIB')) THEN
*---
*  DRAGON LIBRARY
*----
        CALL LCMINF(IPLIB,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
        IF(NAMFIL.EQ.NAMLCM) THEN
           IPDRL=IPLIB
        ELSE
           CALL LCMOP(IPDRL,NAMFIL,2,2,0)
        ENDIF
        CALL LCMLEN(IPDRL,'ENERGY',LENGT,ITYLCM)
        IF(LENGT.EQ.0) THEN
           CALL LCMLIB(IPDRL)
           CALL XABORT('LIBNRG: NO GROUP STRUCTURE AVAILABLE')
        ENDIF
        NGRI=LENGT-1
        ALLOCATE(ENERG(LENGT))
        CALL LCMGET(IPDRL,'ENERGY',ENERG)
        IF(NAMFIL.NE.NAMLCM) CALL LCMCL(IPDRL,1)
      ELSE IF(NAMLBT.EQ.'WIMSAECL') THEN
*---
*  WIMS-AECL LIBRARY
*----
        IUNIT=KDROPN(NAMFIL,IACTO,ILIBDA,LRIND)
        IF(IUNIT.LE.0) THEN
          WRITE(HSMG,'(27HLIBNRG: WIMS-AECL LIBRARY '',A16,8H'' CANNOT,
     >    30H BE OPENED BY KDROPN (ERRCODE=,I2,2H).)') NAMFIL,IUNIT
          CALL XABORT(HSMG)
        ENDIF
        CALL OPNIND(IUNIT,MASTER,LMASTB)
        CALL REDIND(IUNIT,MASTER,LMASIN,GENINX,LGENTB,1)
        CALL REDIND(IUNIT,GENINX,LGENIN,NPZ,LPZ,1)
        NGRI=NPZ(2)
        ALLOCATE(ENERG(NGRI+1))
        CALL REDIND(IUNIT,GENINX,LGENIN,ENERG,NGRI+1,4)
        CALL CLSIND(IUNIT)
      ELSE IF((NAMLBT.EQ.'WIMSD4').OR.(NAMLBT.EQ.'WIMSE')) THEN
*---
*  WIMSD4 OR WIMSE LIBRARY
*----
        IUNIT=KDROPN(NAMFIL,IACTO,ILIBIN,LRIND)
        IF(IUNIT.LE.0) THEN
          WRITE(HSMG,'(22HLIBNRG: WIMS LIBRARY '',A16,9H'' CANNOT ,
     >    29HBE OPENED BY KDROPN (ERRCODE=,I2,2H).)') NAMFIL,IUNIT
          CALL XABORT(HSMG)
        ENDIF
        READ(IUNIT) (NPZ(II),II=1,LPZ-1)
        NGRI=NPZ(2)
        READ(IUNIT) ITEMP
        ALLOCATE(ENERG(NGRI+1))
        READ(IUNIT) (ENERG(J),J=1,NGRI+1)
        IERR=KDRCLS(IUNIT,IACTC)
        IF(IERR.LT.0) THEN
          WRITE(HSMG,'(22HLIBNRG: WIMS LIBRARY '',A16,9H'' CANNOT ,
     >    29HBE CLOSED BY KDRCLS (ERRCODE=,I2,2H).)') NAMFIL,IERR
          CALL XABORT(HSMG)
        ENDIF
      ELSE IF(NAMLBT.EQ.'APLIB1') THEN
*---
*  APOLLO-1 LIBRARY
*----
        IUNIT=KDROPN(NAMFIL,IACTO,ILIBIN,LRIND)
        IF(IUNIT.LE.0) THEN
          WRITE(HSMG,'(26HLIBNRG: APOLLO-1 LIBRARY '',A16,9H'' CANNOT ,
     >    29HBE OPENED BY KDROPN (ERRCODE=,I2,2H).)') NAMFIL,IUNIT
          CALL XABORT(HSMG)
        ENDIF
        REWIND(IUNIT)
 100    CONTINUE
        READ(IUNIT) INDLOR,NR,NIA,(IA(I),I=1,NIA)
        IF(NIA.GT.MAXA)
     >    CALL XABORT('LIBNRG: DIMENSION MAXA =1000 TOO SMALL')
        IF(INDLOR.EQ.9999)
     >    CALL XABORT('LIBNRG: NO GROUP STRUCTURE AVAILABLE')
        NGRI=IA(1)
        IF(IA(3).EQ.0) THEN
          DO 110 K=1,NR
            READ(IUNIT)
 110      CONTINUE
          GO TO 100
        ELSE
          ALLOCATE(ENERG(NGRI+1),UUU(NGRI))
          READ(IUNIT) E0,DEL,(UUU(I),I=1,NGRI)
          E0=1.0E+6*E0
          ENERG=E0
          DO 120 IG=1,NGRI
            ENERG(IG+1)=E0*EXP(-UUU(IG))
 120      CONTINUE
          DEALLOCATE(UUU)
        ENDIF
        IERR=KDRCLS(IUNIT,IACTC)
        IF(IERR.LT.0) THEN
          WRITE(HSMG,'(26HLIBNRG: APOLLO-1 LIBRARY '',A16,9H'' CANNOT ,
     >    29HBE CLOSED BY KDRCLS (ERRCODE=,I2,2H).)') NAMFIL,IERR
          CALL XABORT(HSMG)
        ENDIF
      ELSE IF(NAMLBT.EQ.'APLIB2') THEN
*---
*  APOLLO-2 LIBRARY
*----
        CALL LIBA2G(NAMFIL,NGRI,IPENER)
        CALL C_F_POINTER(IPENER,ENERG,(/ NGRI+1 /))
      ELSE IF(NAMLBT.EQ.'APXSM') THEN
*---
*  APOLLO-XSM LIBRARY
*----
        CALL LIBXS3(NAMFIL,NGRI,IPENER)
        CALL C_F_POINTER(IPENER,ENERG,(/ NGRI+1 /))
      ELSE IF(NAMLBT.EQ.'MATXS') THEN
*---
*  MATXS LIBRARY
*----
        IUNIT=KDROPN(NAMFIL,IACTO,ILIBIN,LRIND)
        IF(IUNIT.LE.0) THEN
          WRITE(HSMG,'(23HLIBNRG: MATXS LIBRARY '',A16,11H'' CANNOT BE,
     >    27H OPENED BY KDROPN (ERRCODE=,I2,2H).)') NAMFIL,IUNIT
          CALL XABORT(HSMG)
        ENDIF
        NWDS=3
        IREC=2
        CALL XDREED(IUNIT,IREC,RA,NWDS)
        NPART=IA(1)
        NTYPE=IA(2)
        IREC=4
        NWDS=(NPART+NTYPE)*MULT+6*NTYPE+NPART
        IF(NWDS.GT.MAXA)
     >    CALL XABORT('LIBNRG: INSUFFICIENT VALUE OF MAXA(1).')
        CALL XDREED(IUNIT,IREC,RA,NWDS)
        NEX1=(NPART+NTYPE)*MULT+6*NTYPE
        DO 180 I=1,NPART
          NGX=IA(NEX1+I)
          WRITE(HPRT,'(A6)') DA(I)
          IREC=IREC+1
          IF(HPRT.EQ.'NEUT'.OR.HPRT.EQ.'neut'.OR.
     >       HPRT.EQ.'N'.OR.HPRT.EQ.'n') THEN
            IF(NGRI.EQ.0) THEN
              NGRI=NGX
              ALLOCATE(ENERG(NGRI+1))
              CALL XDREED(IUNIT,IREC,ENERG,NGRI+1)
            ELSE
              IF(NGX.NE.NGRI)
     >          CALL XABORT('LIBNRG: INVALID GROUP STRUCTURE.')
              ALLOCATE(TEMPE(NGRI+1))
              CALL XDREED(IUNIT,IREC,TEMPE,NGRI+1)
              DO 170 IG=0,NGRI
                IF(TEMPE(IG+1).NE.ENERG(IG+1))
     >            CALL XABORT('LIBNRG: INVALID GROUP STRUCTURE.')
 170          CONTINUE
              DEALLOCATE(TEMPE)
            ENDIF
          ENDIF
 180    CONTINUE
        CALL XDRCLS(IUNIT)
        IERR=KDRCLS(IUNIT,IACTC)
        IF(IERR.LT.0) THEN
          WRITE(HSMG,'(23HLIBNRG: MATXS LIBRARY '',A16,11H'' CANNOT BE,
     >    27H CLOSED BY KDRCLS (ERRCODE=,I2,2H).)') NAMFIL,IERR
          CALL XABORT(HSMG)
        ENDIF
      ELSE IF(NAMLBT.EQ.'MATXS2') THEN
*---
*  MATXS2 LIBRARY
*----
        !AHMED, NAMFIL= _test_nc 
        IF(NAMFIL(:1).EQ.'_') ILIBIN=3
        IUNIT=KDROPN(NAMFIL,IACTO,ILIBIN,LRIND)
        IF(IUNIT.LE.0) THEN
          WRITE(HSMG,'(24HLIBNRG: MATXS2 LIBRARY '',A16,10H'' CANNOT B,
     >    28HE OPENED BY KDROPN (ERRCODE=,I2,2H).)') NAMFIL,IUNIT
          CALL XABORT(HSMG)
        ENDIF
        NWDS=6
        IREC=2
        IF(ILIBIN.EQ.2) THEN
          CALL XDREED(IUNIT,IREC,RA,NWDS)
        ELSE IF(ILIBIN.EQ.3) THEN
          CALL LIBEED(IUNIT,IREC,RA,NWDS)
        ENDIF
        NPART=IA(1)
        NTYPE=IA(2)
        NMAT=IA(4)
        IREC=4
        NWDS=(NPART+NTYPE+NMAT)*MULT+2*NTYPE+NPART+2*NMAT
        
       ! AHMED NPART =1  NTYPE =1 NMAT =2 NWDS  =15

        IF(NWDS.GT.MAXA)
     >    CALL XABORT('LIBNRG: INSUFFICIENT VALUE OF MAXA(2).')
        NEX1=(NPART+NTYPE+NMAT)*MULT
        IF(ILIBIN.EQ.2) THEN
          CALL XDREED(IUNIT,IREC,RA,NWDS)
        ELSE IF(ILIBIN.EQ.3) THEN
          CALL LIBEED(IUNIT,IREC,RA,NWDS)
        ENDIF
        NGX=IA(NEX1+1) ! use the energy mesh of the first particle
        WRITE(HPRT,'(A6)') DA(1) ! name of the first particle
        !HPRT='b'
        
        IREC=IREC+1
        DO 195 IMTYP=1,NMTYP
          IF(HPRT.EQ.HMTYP(IMTYP)) THEN
            IF(NGRI.EQ.0) THEN
              NGRI=NGX
              ALLOCATE(ENERG(NGRI+1))
              IF(ILIBIN.EQ.2) THEN
                CALL XDREED(IUNIT,IREC,ENERG,NGRI+1)
              ELSE
                CALL LIBEED(IUNIT,IREC,ENERG,NGRI+1)
              ENDIF
            ELSE
              IF(NGX.NE.NGRI)
     >          CALL XABORT('LIBNRG: INVALID GROUP STRUCTURE.')
              ALLOCATE(TEMPE(NGRI+1))
              IF(ILIBIN.EQ.2) THEN
                CALL XDREED(IUNIT,IREC,TEMPE,NGRI+1)
              ELSE
                CALL LIBEED(IUNIT,IREC,TEMPE,NGRI+1)
              ENDIF
              DO 190 IG=0,NGRI
                IF(TEMPE(IG+1).NE.ENERG(IG+1))
     >            CALL XABORT('LIBNRG: INVALID GROUP STRUCTURE.')
 190          CONTINUE
              DEALLOCATE(TEMPE)
            ENDIF
          ENDIF
 195    CONTINUE
        IF(ILIBIN.EQ.2) THEN
          CALL XDRCLS(IUNIT)
        ELSE
          CALL LIBCLS()
        ENDIF
        IERR=KDRCLS(IUNIT,IACTC)
        IF(IERR.LT.0) THEN
          WRITE(HSMG,'(24HLIBNRG: MATXS2 LIBRARY '',A16,10H'' CANNOT B,
     >    28HE CLOSED BY KDRCLS (ERRCODE=,I2,2H).)') NAMFIL,IERR
          CALL XABORT(HSMG)
        ENDIF
      ELSE IF(NAMLBT.EQ.'NDAS') THEN

*---
*  WIMS-NDAS LIBRARY
*----
        CALL LIBND0(NAMFIL,NGRI,IPENER)
        CALL C_F_POINTER(IPENER,ENERG,(/ NGRI+1 /))
      ENDIF
      IF(ENERG(NGRI+1).EQ.0.0) ENERG(NGRI+1)=1.0E-5
      IF(NGT.EQ.0) THEN
*----
*  IF NGT=0 SAVE GROUP STRUCTURE AND SET GROUP PARAMETERS
*----
        NGROUP=NGRI
        CALL LCMPUT(IPLIB,'ENERGY',NGRI+1,2,ENERG)
        JG=0
        DO 210 IG=1,NGROUP
          ENERG(JG+1)=LOG(ENERG(JG+1)/ENERG(JG+2))
          JG=JG+1
 210    CONTINUE
        CALL LCMPUT(IPLIB,'DELTAU',NGROUP,2,ENERG)
        NGT=NGROUP
      ELSE IF(NGRI.EQ.NGT) THEN
*----
*  IF NGT>0 VALIDATE GROUP STRUCTURE
*----
        ALLOCATE(EIER(NGT+1))
        CALL LCMGET(IPLIB,'ENERGY',EIER)
        JG=0
        DO 220 IG=1,NGT
          ERROR=ABS(ENERG(JG+1)-EIER(JG+1))
          IF(ERROR.GT.ABS(ENERG(JG+1))*1.0E-4) THEN
            WRITE(IOUT,'(1X,A20)') 'OLD GROUP STRUCTURE='
            WRITE(IOUT,'(1P,5E15.7)')
     >        (EIER(IPR+1),IPR=0,NGT)
            WRITE(IOUT,'(1X,A20)') 'NEW GROUP STRUCTURE='
            WRITE(IOUT,'(1P,5E15.7)')
     >        (ENERG(IPR+1),IPR=0,NGT)
            WRITE(IOUT,'(7H ERROR=,1P,E10.3,9H IN GROUP,I4)')
     >        ERROR,IG
            CALL XABORT('LIBNRG: INCOMPATIBLE GROUP STRUCTURE')
          ENDIF
          JG=JG+1
 220    CONTINUE
        DEALLOCATE(EIER)
      ELSE
        WRITE(IOUT,'(1X,A20,1X,I10)') 'OLD NUMBER OF GROUPS=',NGT
        WRITE(IOUT,'(1X,A20,1X,I10)') 'NEW NUMBER OF GROUPS=',NGRI
        CALL XABORT('LIBNRG: INCOMPATIBLE NUMBER OF GROUPS')
      ENDIF
      IF((NAMLBT.EQ.'NDAS').OR.(NAMLBT.EQ.'APLIB2').OR.
     >   (NAMLBT.EQ.'APXSM')) THEN
        CALL LCMDRD(IPENER)
      ELSE
        DEALLOCATE(ENERG)
      ENDIF
*----
*  RETURN
*----
      RETURN
      END
