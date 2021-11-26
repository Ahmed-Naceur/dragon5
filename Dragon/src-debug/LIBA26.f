*DECK LIBA26
      SUBROUTINE LIBA26(LGSEG,IG,NGBIN,IUNIT,LBLOC,TKCARO,TCAROB,NSIGF,
     1 TT,NTEMPS,TEMPS,DELTF,SIGTF,SIGAF,DELINF,SGTINF,SGAINF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Temperature interpolation of autolib (bin cross sections) information.
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
* LGSEG   dimension of the directory block.
* IG      coarse energy group under consideration.
* NGBIN   number of coarse energy groups.
* IUNIT   APOLIB-2 file unit number.
* LBLOC   number of words in the direct access buffer.
* TKCARO  index array used to parse tcarob.
* TCAROB  directory block.
* NSIGF   number of fine energy groups.
* TT      temperature of isotope.
* NTEMPS  number of tabulated temperatures.
* TEMPS   tabulated temperatures.
*
*Parameters: output
* DELTF   fine group lethargy widths.
* SIGTF   fine group total x-s.
* SIGAF   fine group absorption x-s.
* DELINF  calculated lethargy width for group IG.
* SGTINF  calculated infinite-dilution total x-s for group IG.
* SGAINF  calculated infinite-dilution absorption x-s for group IG.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LGSEG,TKCARO(31),TCAROB(LGSEG),IG,NGBIN,IUNIT,LBLOC,NSIGF,
     1 NTEMPS
      REAL TT,TEMPS(NTEMPS),DELTF(NSIGF),SIGTF(NSIGF),SIGAF(NSIGF),
     1 DELINF,SGTINF,SGAINF
*----
*  LOCAL VARIABLES
*----
      EXTERNAL LIBA21
      CHARACTER HSMG*131,TYPSEG*8
      PARAMETER (NINT=2,DTMIN=1.0)
      DOUBLE PRECISION D1,D2,D3
      LOGICAL LOK
      TYPE(C_PTR) ICHDIM_PTR,ICHTYP_PTR,ICHDKL_PTR,TSEGM_PTR
      INTEGER, POINTER, DIMENSION(:) :: ICHDIM,ICHTYP,ICHDKL,ITSEGM
      REAL, POINTER, DIMENSION(:) :: RTSEGM
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SQRTEM,WEIJHT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WEIJHT(NTEMPS),SQRTEM(NTEMPS))
*----
*  COMPUTE THE WEIGHTS.
*----
      DO 10 I=1,NTEMPS
      SQRTEM(I)=SQRT(TEMPS(I))
   10 CONTINUE
      IF(NTEMPS.EQ.1) THEN
        IPROX=1
        IGTFIX=1
      ELSE
        STT=SQRT(TT)
        CALL LIBA28(STT,SQRTEM,NTEMPS,NINT,WEIJHT,IORD,IPROX,I0)
        IF(ABS(TT-TEMPS(IPROX)).LE.DTMIN) THEN
          IGTFIX=1
        ELSEIF((STT.LT.SQRTEM(1)).OR.(STT.GT.SQRTEM(NTEMPS))) THEN
          WRITE(HSMG,'(A,F8.2,A,F8.2,A,F8.2)')
     1    'LIBA26: A TEMPSERATURE', TT,'K IS NOT INCLUDED BETWEEN ',
     2    TEMPS(1),' AND ',TEMPS(NTEMPS)
          WRITE(6,'(/1X,A)') HSMG
          IGTFIX=2
        ELSE
          IGTFIX=0
        ENDIF
      ENDIF
*----
*  LOOP OVER TABULATED TEMPERATURES.
*----
      D1=0.0D0
      IDKDS=1-TKCARO(10)
      IDKTS=1-TKCARO(23)
      IDKLS=TKCARO(8)
      JDKDS=TCAROB(IDKDS)
      JDKTS=TCAROB(IDKTS)
      CALL XDRSET(SIGTF,NSIGF,0.0)
      CALL XDRSET(SIGAF,NSIGF,0.0)
      DO 50 J=1,IORD
      IT=I0+J
      IS=(IT-1)*NGBIN+IG
      IDK=JDKTS+8*(IS-1)
      CALL AEXCPC(IDK,8,TCAROB,TYPSEG)
      LNGS=TCAROB(IDKLS+IS)
      IF(LNGS.LE.0) CALL XABORT('LIBA26: INVALID PTHOM5(1).')
      JDKS=TCAROB(JDKDS+IS)
      CALL AEXTRT(LIBA21,TYPSEG,NBRTYP,ICHDIM_PTR,ICHTYP_PTR,ICHDKL_PTR)
      CALL C_F_POINTER(ICHDIM_PTR,ICHDIM,(/ NBRTYP /))
      CALL C_F_POINTER(ICHTYP_PTR,ICHTYP,(/ NBRTYP /))
      CALL C_F_POINTER(ICHDKL_PTR,ICHDKL,(/ NBRTYP /))
      TSEGM_PTR=LCMARA(LNGS+1)
      CALL C_F_POINTER(TSEGM_PTR,ITSEGM,(/ LNGS+1 /))
      CALL C_F_POINTER(TSEGM_PTR,RTSEGM,(/ LNGS+1 /))
      CALL AEXDIR(IUNIT,LBLOC,ITSEGM,JDKS,LNGS+1)
      CALL AEXGNV(1,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDD,NV)
      IF(NV.NE.NSIGF) CALL XABORT('LIBA26: INVALID PTHOM5(2).')
      CALL AEXGNV(3,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDT,NV)
      CALL AEXGNV(5,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDA,NV)
      CALL LCMDRD(ICHDIM_PTR)
      CALL LCMDRD(ICHTYP_PTR)
      CALL LCMDRD(ICHDKL_PTR)
      IF(IT.EQ.I0+1) THEN
         D1=0.0D0
         DO 20 I=1,NSIGF
         DELTF(I)=RTSEGM(IDD+I-1)
         D1=D1+DELTF(I)
   20    CONTINUE
      ELSE
         LOK=.TRUE.
         DO 30 I=1,NSIGF
         LOK=LOK.AND.(DELTF(I).EQ.RTSEGM(IDD+I-1))
   30    CONTINUE
         IF(.NOT.LOK) CALL XABORT('LIBA26: INVALID AUTOLIB MESH.')
      ENDIF
      DO 40 I=1,NSIGF
      SIGTF(I)=SIGTF(I)+REAL(WEIJHT(J)*RTSEGM(IDT+I-1))
      SIGAF(I)=SIGAF(I)+REAL(WEIJHT(J)*RTSEGM(IDA+I-1))
   40 CONTINUE
      CALL LCMDRD(TSEGM_PTR)
   50 CONTINUE
      D2=0.0D0
      D3=0.0D0
      DO 60 I=1,NSIGF
      SIGTF(I)=MAX(SIGTF(I),0.0)
      SIGAF(I)=MAX(SIGAF(I),0.0)
      D2=D2+SIGTF(I)*DELTF(I)
      D3=D3+SIGAF(I)*DELTF(I)
   60 CONTINUE
      DELINF=REAL(D1)
      SGTINF=REAL(D2/D1)
      SGAINF=REAL(D3/D1)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SQRTEM,WEIJHT)
      RETURN
      END
