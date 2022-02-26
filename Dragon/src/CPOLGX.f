*DECK CPOLGX
      SUBROUTINE CPOLGX(IPLIB ,IGS   ,IPRINT,IORD  ,NGROUP,INDPRO,
     >                  XSREC ,ITYPRO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Get/save standard vectorial cross section data from/on IPLIB.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
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
* IORD    cross section order:
*         =1 constant;
*         =2 linear;
*         =3 quadratic.
* NGROUP  number of energy groups.
* INDPRO  vector for cross section to process:
*         =0 do not process;
*         >0 process.
*
*Parameters: input/output
* XSREC   cross section table.
*
*Parameters: output
* ITYPRO  vector for cross section processed indices:
*         =0 absent  (not processed);
*         >0 present (processed).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          NDPROC
      PARAMETER       (NDPROC=20)
      TYPE(C_PTR)      IPLIB
      INTEGER          IGS,IPRINT,IORD,NGROUP,INDPRO(NDPROC),
     >                 ITYPRO(NDPROC)
      REAL             XSREC(NGROUP,NDPROC)
*----
*  LOCAL PARAMETERS
*  NDPROC = NUMBER OF DEFAULT CROSS SECTIONS = 20
*  NAMDXS = NAME OF NDPROC DEFAULT XS
*----
      INTEGER          IOUT
      PARAMETER       (IOUT=6)
      CHARACTER        NAMDXS(NDPROC)*6,NORD*6,TEXT6*6,TEXT12*12,NAMT*12
      INTEGER          IODIV,LONG,ITYP,IXSR,IXSTN,IG,JG
      SAVE             NAMDXS
      DATA    NAMDXS  /'NTOT0 ','TRANC ','NUSIGF','NFTOT ','CHI   ',
     >                 'NU    ','NG    ','NHEAT ','N2N   ','N3N   ',
     >                 'N4N   ','NP    ','NA    ','GOLD  ','ABS   ',
     >                 'NWT0  ','STRD  ','STRD X','STRD Y','STRD Z'/
      IODIV=0
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
*----
*  READ/INITIALIZE STATE VECTOR
*----
      CALL LCMLEN(IPLIB,'XS-SAVED',LONG,ITYP)
      IF(LONG.EQ.NDPROC) THEN
        CALL LCMGET(IPLIB,'XS-SAVED',ITYPRO)
      ELSE IF(LONG.EQ.0) THEN
        CALL XDISET(ITYPRO,NDPROC,0)
        NAMT=' '
        CALL LCMNXT(IPLIB,NAMT)
        TEXT12=NAMT
   80   CALL LCMLEN(IPLIB,NAMT,LONG,ITYP)
        IF(ITYP.EQ.2) THEN
           DO 90 IXSR=1,NDPROC
           IF(NAMT(:6).EQ.NAMDXS(IXSR)) ITYPRO(IXSR)=1
   90      CONTINUE
        ENDIF
        CALL LCMNXT(IPLIB,NAMT)
        IF(NAMT.NE.TEXT12) GO TO 80
      ELSE
        WRITE(IOUT,9000) NDPROC,LONG
        CALL XABORT('CPOLGX: INVALID VALUE FOR NDPROC')
      ENDIF
      IF(IGS.GT.0) THEN
*----
*  SAVE LOCAL DEFAULT XS IF REQUIRED
*----
        IF(IGS.EQ.1) THEN
          DO 100 IXSR=1,NDPROC
            TEXT6=NAMDXS(IXSR)
            IF(IXSR.EQ.1) TEXT6='TOTAL'
            IF(INDPRO(IXSR).EQ.1) THEN
              IXSTN=MOD(ITYPRO(IXSR)/IODIV,2)
*----
*  FIND IF XS NOT ALL 0.0
*----
              DO 110 IG=1,NGROUP
                IF(XSREC(IG,IXSR).NE.0.0) THEN
                  IF(IXSTN.EQ.0) THEN
                    ITYPRO(IXSR)=ITYPRO(IXSR)+IODIV
                    IXSTN=1
                  ENDIF
                  GO TO 115
                ENDIF
 110          CONTINUE
 115          CONTINUE
              IF((IXSTN.NE.0).OR.(IXSR.EQ.2)) THEN
                CALL LCMPUT(IPLIB,TEXT6//NORD,NGROUP,2,XSREC(1,IXSR))
              ENDIF
            ENDIF
 100      CONTINUE
        ENDIF
        CALL LCMPUT(IPLIB,'XS-SAVED',NDPROC,1,ITYPRO)
      ELSE
*----
*  GET LOCAL DEFAULT XS IF REQUIRED
*----
        IF(IGS.EQ.-1) THEN
          DO 200 IXSR=1,NDPROC
            TEXT6=NAMDXS(IXSR)
            IF(IXSR.EQ.1) TEXT6='NTOT0'
            IF(INDPRO(IXSR).EQ.1) THEN
              IXSTN=MOD(ITYPRO(IXSR)/IODIV,2)
*----
*  READ IF IXSTN = 1
*  INITIALIZE TO 0.0 IF IXSTN = 0
*----
              IF(IXSTN.EQ.1) THEN
                CALL LCMLEN(IPLIB,TEXT6//NORD,LONG,ITYP)
                IF(LONG .EQ. 0) THEN 
                  CALL XDRSET(XSREC(1,IXSR),NGROUP,0.0)
                ELSE
                  CALL LCMGET(IPLIB,TEXT6//NORD,XSREC(1,IXSR))
                ENDIF
              ELSE
                CALL XDRSET(XSREC(1,IXSR),NGROUP,0.0)
              ENDIF
            ENDIF
 200      CONTINUE
        ENDIF
      ENDIF
      IF(IPRINT .GE. 100) THEN
*----
*  Print XS
*----
        DO IXSR=1,NDPROC
          IF(INDPRO(IXSR).EQ.1) THEN
            IXSTN=MOD(ITYPRO(IXSR)/IODIV,2)
            IF(IXSTN.NE.0) THEN
              DO IG=1,NGROUP
                IF(XSREC(IG,IXSR).NE.0.0) THEN
                  WRITE(IOUT,6000) NAMDXS(IXSR)//NORD
                  WRITE(IOUT,6010) (XSREC(JG,IXSR),JG=1,NGROUP)
                  GO TO 210
                ENDIF
              ENDDO
            ENDIF
 210        CONTINUE
          ENDIF
        ENDDO
      ENDIF
      RETURN
*----
*  ABORT FORMAT
*----
 9000 FORMAT(' CPOLGX: ****** ABORT ******'/
     >       ' INVALID LENGTH OF RECORD XS-SAVED '/
     >       ' STORAGE SPACE NDPROC   = ',I10/
     >       ' LENGTH OF RECORD LONG  = ',I10/
     >       ' ***************************')
 6000 FORMAT(/' CROSS SECTION TYPE    = ',A12)
 6010 FORMAT(1P,5E16.7)
      END
