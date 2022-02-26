*DECK CPOISO
      SUBROUTINE CPOISO(IPRINT,IEXTRC,NMERGE,MAXISO,MAXISM,NBMICR,
     >                  NISCPO,NISEXT,ISOCPO,ISOEXT,ISOORD,ISOTMP,
     >                  IMXTMP,IDIMIX,NBIMRG,ICOMIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Identify isotopes to be extracted from macroscopic xs and isotopes
* included in new combined isotopes.
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
* IPRINT  print parameter. Equal to zero for no print.
* IEXTRC  type of extraction: 1 for part 2 for all.
* NMERGE  number of region.
* MAXISO  maximum nunber of isotopes permitted.
* MAXISM  maximum nunber of isotopes per region.
* NBMICR  maximum number of isotopes in EDIT.
* NISCPO  number of Compo isotopes treated.
* NISEXT  number of extracted isotopes treated.
* ISOCPO  Compo name of isotopes.
* ISOEXT  name of extracted isotopes.
* ISOORD  order of extracted isotopes.
* ISOTMP  name of isotopes in EDIT.
* IMXTMP  mixture of isotopes in EDIT.
*
*Parameters: output
* IDIMIX  isotopes identifier in each Compo material.
* NBIMRG  final number of isotope per region.
* ICOMIX  pointer to Compo isotope for region.
*
*-----------------------------------------------------------------------
*
      IMPLICIT      NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER       IPRINT,IEXTRC,NMERGE,MAXISO,MAXISM,NBMICR,NISCPO,
     >              NISEXT,ISOCPO(3,MAXISO),
     >              ISOEXT(3,MAXISO),ISOORD(MAXISO),ISOTMP(3,NBMICR),
     >              IMXTMP(NBMICR),IDIMIX(NMERGE,NBMICR),
     >              NBIMRG(NMERGE),ICOMIX(NMERGE,MAXISM)
*----
*  LOCAL PARAMETERS
*----
      INTEGER       IOUT
      CHARACTER     TEXT4*4
      PARAMETER    (IOUT=6)
      INTEGER       ISOM,ISOE,ISOC,IMRG,ITEXT4,ITC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDETMP
*----
*  SCRATCH STORAGE ALLOCATION
*   IDETMP  extracted isotopes number associated with EDIT isotope.
*----
      ALLOCATE(IDETMP(NBMICR))
      CALL XDISET(IDETMP,NBMICR,0)
*----
*  STORE IN ITEXT4 BLANCK STRING
*----
      TEXT4='    '
      READ(TEXT4,'(A4)') ITEXT4
*----
*  IF EXTRACT ALL USED (IEXTRC=2)
*  GENERATE ISOCPO, ISOEXT AND ISOORD
*  ASSOCIATE WITH ALL ISOTOPES EXTRACTED ISOTOPE NUMBER
*  NAMELY IDETMP(ISOM)=ISOEXT(ISOE)
*----
      IF(IEXTRC.EQ.2) THEN
        NISEXT=0
        DO 100 ISOM=1,NBMICR
          DO 110 ISOE=1,NISEXT
            IF(ISOEXT(1,ISOE).EQ.ISOTMP(1,ISOM).AND.
     >         ISOEXT(2,ISOE).EQ.ISOTMP(2,ISOM)) GO TO 115
 110      CONTINUE
          IF(NISEXT.EQ.MAXISO) THEN
            WRITE(IOUT,7000) MAXISO,ISOTMP(1,ISOM),ISOTMP(2,ISOM)
          ELSE
            NISEXT=NISEXT+1
            ISOEXT(1,NISEXT)=ISOTMP(1,ISOM)
            ISOEXT(2,NISEXT)=ISOTMP(2,ISOM)
            ISOEXT(3,NISEXT)=ITEXT4
            ISOCPO(1,NISEXT)=ISOTMP(1,ISOM)
            ISOCPO(2,NISEXT)=ISOTMP(2,ISOM)
            ISOCPO(3,NISEXT)=ITEXT4
            ISOORD(NISEXT)=NISEXT
            IDETMP(ISOM)=NISEXT
          ENDIF
 115      CONTINUE
 100    CONTINUE
        NISCPO=NISEXT
      ELSE
*----
*  IF SPECIFIC ISOTOPES EXTRACTED (IEXTRC=1)
*  FOR GENERIC EXTRACTED NAME (ISOEXT(3,ISOE)='    ')
*  ASSOCIATE WITH SET OF ISOTOPE EXTRACTED ISOTOPE NUMBER
*  NAMELY IDETMP(ISOM)=ISOEXT(ISOE)
*  FOR EXPLICIT EXTRACTED NAMES
*  ASSOCIATE WITH SPECIFIC ISOTOPE EXTRACTED ISOTOPE NUMBER
*  NAMELY IDETMP(ISOM)=ISOEXT(ISOE)
*----
        DO 120 ISOE=1,NISEXT
          IF(ISOEXT(3,ISOE).EQ.ITEXT4) THEN
            DO 130 ISOM=1,NBMICR
              IF(ISOEXT(1,ISOE).EQ.ISOTMP(1,ISOM).AND.
     >           ISOEXT(2,ISOE).EQ.ISOTMP(2,ISOM)) THEN
                IDETMP(ISOM)=ISOE
              ENDIF
 130        CONTINUE
          ELSE
            DO 140 ISOM=1,NBMICR
              IF(ISOEXT(1,ISOE).EQ.ISOTMP(1,ISOM).AND.
     >           ISOEXT(2,ISOE).EQ.ISOTMP(2,ISOM).AND.
     >           ISOEXT(3,ISOE).EQ.ISOTMP(3,ISOM)) THEN
                IDETMP(ISOM)=ISOE
              ENDIF
 140        CONTINUE
          ENDIF
 120    CONTINUE
      ENDIF
*----
*  IDENTIFY EXTRACTED ISOTOPES
*----
      DO 150 ISOM=1,NBMICR
        IMRG=IMXTMP(ISOM)
        ISOE=IDETMP(ISOM)
        IF(IMRG.NE.0.AND.ISOE.NE.0) THEN
          IDIMIX(IMRG,ISOM)=ISOORD(ISOE)
        ENDIF
 150  CONTINUE
*----
*  COMPUTED NUMBER OF ISOTOPES PER MIXTURE
*----
      DO 160 IMRG=1,NMERGE
        NBIMRG(IMRG)=0
        DO 170 ISOM=1,NBMICR
          ISOC=IDIMIX(IMRG,ISOM)
          IF(ISOC.NE.0) THEN
            DO 180 ISOE=1,NBIMRG(IMRG)
              IF(ISOC.EQ.ICOMIX(IMRG,ISOE)) GO TO 185
 180        CONTINUE
            NBIMRG(IMRG)=NBIMRG(IMRG)+1
            ICOMIX(IMRG,NBIMRG(IMRG))=ISOC
 185        CONTINUE
          ENDIF
 170    CONTINUE
 160  CONTINUE
      IF(IPRINT.GE.1) THEN
        WRITE(IOUT,6000)
        DO 190 IMRG=1,NMERGE
          IF(NBIMRG(IMRG).GT.0) THEN
            DO 191 ISOM=1,NBMICR
              ISOC=IDIMIX(IMRG,ISOM)
              IF(ISOC.NE.0) THEN
                WRITE(IOUT,6001) IMRG,(ISOCPO(ITC,ISOC),ITC=1,3),
     >            (ISOTMP(ITC,ISOM),ITC=1,3)
              ENDIF
 191        CONTINUE
          ENDIF
 190    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IDETMP)
      RETURN
*----
*  PRINT FORMAT
*----
 6000 FORMAT(' CPO: LIST OF EXTRACTED ISOTOPES:'/
     > 10X,'REGION',10X,'CPO NAME    ',10X,'EDIT NAME   ')
 6001 FORMAT(10X,I6,10X,3A4,' CONTAINS ',3A4)
*----
*  WARNING FORMAT
*----
 7000 FORMAT(' CPOISO: ****** WARNING ******'/
     >       ' MAXIMUM NUMBER OF ISOTOPE REACHED = ',I8/
     >       ' SKIP GENERIC ISOTOPE NAME         = ',2A4/
     >       ' *****************************')
      END
