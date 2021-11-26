*DECK XDRLXS
      SUBROUTINE XDRLXS(IPLIB ,IGS   ,IPRINT,NPROC ,NAMDXS,IORD  ,
     >                  NGROUP,XSREC )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Get/save Legendre-independent cross section data from/on IPLIB.
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
* IPLIB   pointer to the internal library.
* IGS     get or save flag:
*         >0 save;
*         <0 get.
* IPRINT  Print level (cross sections printed if IPRINT>99).
* NPROC   number of Legendre-independent terms to process.
* NAMDXS  names of cross sections to process.
* IORD    cross section order:
*         =1 constant;
*         =2 linear;
*         =3 quadratic.
* NGROUP  number of energy groups.
*
*Parameters: input/output
* XSREC   cross section records for IRPROC=1,NPROC.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPLIB
      INTEGER          IGS,IPRINT,NPROC,IORD,NGROUP
      REAL             XSREC(NGROUP,NPROC)
      CHARACTER        NAMDXS(NPROC)*6,NORD*6
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT
      PARAMETER       (IOUT=6)
      INTEGER          IODIV,IXSR,IG,JG,ILENG,ITYLCM
*
      IF(IORD.EQ.1) THEN
        NORD='      '
        IODIV=1
      ELSE IF(IORD.EQ.2) THEN
        NORD='   LIN'
        IODIV=2
      ELSE IF(IORD.EQ.3) THEN
        NORD='   QUA'
        IODIV=4
      ENDIF
      IF(NPROC.LE.0) THEN
        CALL XABORT('XDRLXS: ZERO OR NEGATIVE VALUE OF NPROC')
      ENDIF
*
      IF(IGS.GT.0) THEN
*----
*  SAVE LOCAL DEFAULT XS IF REQUIRED
*----
*
        DO 100 IXSR=1,NPROC
*----
*  FIND IF XS NOT ALL 0.0
*----
          DO 110 IG=1,NGROUP
            IF(XSREC(IG,IXSR).NE.0.0) GO TO 115
 110      CONTINUE
          GO TO 100
*----
*  SAVE IF XS NOT ALL 0.0
*----
 115      CALL LCMPUT(IPLIB,NAMDXS(IXSR)//NORD,NGROUP,2,XSREC(1,IXSR))
 100    CONTINUE
      ELSE
*----
*  GET LOCAL DEFAULT XS IF REQUIRED
*----
        DO 200 IXSR=1,NPROC
          CALL XDRSET(XSREC(1,IXSR),NGROUP,0.0)
          CALL LCMLEN(IPLIB,NAMDXS(IXSR)//NORD,ILENG,ITYLCM)
          IF(ILENG.GT.0) THEN
            CALL LCMGET(IPLIB,NAMDXS(IXSR)//NORD,XSREC(1,IXSR))
          ENDIF
 200    CONTINUE
      ENDIF
      IF(IPRINT .GE. 100) THEN
*----
*  Print XS
*----
        DO IXSR=1,NPROC
          DO IG=1,NGROUP
            IF(XSREC(IG,IXSR).NE.0.0) THEN
              WRITE(IOUT,6000) NAMDXS(IXSR)//NORD
              WRITE(IOUT,6010) (XSREC(JG,IXSR),JG=1,NGROUP)
              GO TO 210
            ENDIF
          ENDDO
 210      CONTINUE
        ENDDO
      ENDIF
      RETURN
*----
*  Formats
*----
 6000 FORMAT(/' CROSS SECTION TYPE    = ',A12)
 6010 FORMAT(1P,5E16.7)
      END
