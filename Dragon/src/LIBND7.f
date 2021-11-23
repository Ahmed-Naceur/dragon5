*DECK LIBND7
      SUBROUTINE LIBND7 (MAXDIL,NGRO,NAMFIL,HNISOR,NDIL,DILUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the dilutions corresponding to a resonant isotope within a
* library in NDAS format.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXDIL  maximum number of dilutions.
* NGRO    number of energy groups.
* NAMFIL  NDAS library name.
* HNISOR  library name of the isotope.
*
*Parameters: output
* NDIL    number of finite dilutions.
* DILUT   dilutions.
*
*Reference:
* Copyright (C) from NDAS Atomic Energy of Canada Limited utility (2006)
*
*-----------------------------------------------------------------------
*
      USE FSDF
      IMPLICIT NONE
*----
*  Subroutine arguments
*----
      CHARACTER NAMFIL*(*),HNISOR*12
      INTEGER MAXDIL,NGRO,NDIL
      REAL DILUT(MAXDIL)
*----
*  Local variables
*----
      INTEGER I,IND,IERR,ISOID,HEADER(16),NISOLB,NF,NTEM,IHEAD(200)
      REAL RHEAD(200)
      CHARACTER HSMG*131,TEXT8*8
      EQUIVALENCE(RHEAD(1),IHEAD(1))
*----
*  Read NDAS library parameters
*----
      CALL XSDOPN(NAMFIL,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND7: XSDOPN could not open Library'
     >  //' files')
      CALL XSDBLD(6001,HEADER,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND7: XSDBLD could not read library'
     > //' parameters')
      IF(NGRO.NE.HEADER(2)) CALL XABORT('LIBND7: Invalid number of e'
     > //'nergy groups')
      NISOLB=HEADER(1)
*----
*  Read through NDAS file
*----
      IND=0
      DO I=1,NISOLB
        CALL XSDNAM(I,ISOID,TEXT8,IERR)
        IF(IERR.NE.0) CALL XABORT('LIBND7: XSDNAM index overflow')
        IF(TEXT8.EQ.HNISOR(:8)) THEN
          IND=I
          GO TO 10
        ENDIF
      ENDDO
      WRITE (HSMG,100) HNISOR,NAMFIL
      CALL XABORT(HSMG)
*     Load nuclide header
   10 CALL XSDISO(7000,6001,IND,RHEAD,IERR)
      NF=IHEAD(5)
      IF((NF.GE.1).AND.(NF.LE.3)) THEN
        CALL XSDRES(IND,IHEAD,IERR)
        NTEM=IHEAD(1)
        NDIL=IHEAD(2)
        IF(NDIL.GT.MAXDIL) CALL XABORT('LIBND7: MAXDIL overflow')
        DO I=1,NDIL
          DILUT(I)=RHEAD(2+NTEM+I)
        ENDDO
        NDIL=NDIL-1
      ELSE
        NDIL=0
      ENDIF
      CALL XSDCL()
      RETURN
*
  100 FORMAT(26HLIBND7: Material/isotope ',A12,20H' is missing on NDAS,
     > 12H file named ,A24,1H.)
      END
