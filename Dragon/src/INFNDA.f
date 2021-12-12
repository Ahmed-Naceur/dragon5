*DECK INFNDA
      SUBROUTINE INFNDA(CFILNA,IPRINT,NBISO,HNAMIS,AWRISO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover isotopic masses for isotopes of NDAS-type libraries.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert 
*
*Parameters: input
* CFILNA  name of the NDAS file.
* IPRINT  print flag.
* NBISO   number of isotopes present in the calculation domain.
* HNAMIS  isotope names.
*
*Parameters: output
* AWRISO  isotopic masses.
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
      INTEGER IPRINT,NBISO
      CHARACTER CFILNA*(*),HNAMIS(NBISO)*8
      REAL AWRISO(NBISO)
*----
*  Local variables
*----
      INTEGER IOUT,MAXISO
      PARAMETER(IOUT=6,MAXISO=500)
      CHARACTER TEXT8*8,HSMG*131
      INTEGER I,ISO,IND,IERR,NEL,ISOID,ISONRF(2),HEADER(16),
     > HNAM(2,MAXISO)
      REAL RHEAD(200)
*----
*  Read NDAS library parameters
*----
      IF(CFILNA.EQ.' ' )THEN
        CALL XABORT('INFNDA: NDAS library has not been set')
      ENDIF
      CALL XSDOPN(CFILNA,IERR)
      IF(IERR.NE.0) CALL XABORT('INFNDA: XSDOPN could not open Library'
     >  //' files')
      CALL XSDBLD(6001,HEADER,IERR)
      IF(IERR.NE.0) CALL XABORT('INFNDA: XSDBLD could not read library'
     > //' parameters')
      NEL=HEADER(1)
      IF(NEL.GT.MAXISO) THEN
        WRITE(IOUT,30) MAXISO,NEL
        CALL XABORT('INFNDA: Invalid number of isotopes')
      ENDIF
*----
*  Recover the isotope names and identifiers from the library
*----
      DO I=1,NEL
        CALL XSDNAM(I,ISOID,TEXT8,IERR)
        IF(IERR.NE.0) CALL XABORT('INFNDA: XSDNAM index overflow')
        READ(TEXT8,'(2A4)') HNAM(1,I),HNAM(2,I)
      ENDDO
*----
*  Read through NDAS file and accumulate isotopic mass values
*----
      DO ISO=1,NBISO
        READ(HNAMIS(ISO),'(2A4)') (ISONRF(I),I=1,2)
        IND=0
        DO I=1,NEL
          IF((ISONRF(1).EQ.HNAM(1,I)).AND.
     >       (ISONRF(2).EQ.HNAM(2,I))) THEN
            IND=I
            GO TO 10
          ENDIF
        ENDDO
        WRITE(HSMG,30) HNAMIS(ISO),CFILNA
        CALL XABORT(HSMG)
*       Load nuclide header
   10   CALL XSDISO(7000,6001,IND,RHEAD,IERR)
        AWRISO(ISO)=RHEAD(3)
        IF(IPRINT.GE.100) WRITE(IOUT,40) HNAMIS(ISO),AWRISO(ISO)
      ENDDO
      CALL XSDCL()
      RETURN
*
   30 FORMAT('INFNDA: MATERIAL/ISOTOPE ',A8,
     >       ' IS MISSING ON NDAS LIBRARY FILE ',A8)
   40 FORMAT('INFNDA: DRAGON ISOTOPE =',A8,
     >       ' HAS ATOMIC WEIGHT RATIO = ',F12.5)
      END
