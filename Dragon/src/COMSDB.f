*DECK COMSDB
      SUBROUTINE COMSDB(IMPX,IPCPO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Show the structure of a multicompo file.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Chambon
*
*Parameters: input
* IMPX    print parameter.
* IPCPO   pointer to the multicompo.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX
      TYPE(C_PTR) IPCPO
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXPAR=50,MAXVAL=1000)
      TYPE(C_PTR) JPCPO,KPCPO
      INTEGER ISTATE(NSTATE),NVPO(2),NVALUE(2*MAXPAR)
      CHARACTER RECNAM*12,TEXT12*12,PARFMT(MAXPAR)*8,
     1 VCHAR(MAXVAL)*12,PARKEY(MAXPAR)*12,PARCPO(MAXPAR)*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MUPCPO,JDEBAR,JARBVA,VINTE
      REAL, ALLOCATABLE, DIMENSION(:) :: VREAL
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MUPCPO(2*MAXPAR))
*
      CALL LCMGET(IPCPO,'STATE-VECTOR',ISTATE)
      NMIL=ISTATE(1)
      NCAL=ISTATE(3)
      MAXCAL=ISTATE(4)
      NPAR=ISTATE(5)
      NLOC=ISTATE(6)
      IF(NPAR+NLOC.GT.2*MAXPAR) CALL XABORT('COMSDB: MAXPAR OVERFLOW.')
      IF(NCAL.EQ.0) WRITE(6,*) 'The multi-compo DB is empty.'
*----
*  MAIN LOOP OVER THE HOMOGENEOUS MIXTURES *********************
*----
      JPCPO=LCMGID(IPCPO,'MIXTURES')
      DO 190 IBM=1,NMIL
      KPCPO=LCMDIL(JPCPO,IBM)
*----
*  MAIN LOOP OVER THE NCALR ELEMENTARY CALCULATIONS OF THE COMPO
*----
      WRITE(6,*) 'LIST OF "MUPLET" included in the COMPO'
      DO 170 ICAL=1,NCAL
*----
*  COMPUTE THE MUPLET VECTOR FROM THE COMPO
*----
      CALL LCMSIX(KPCPO,'TREE',1)
      CALL LCMGET(KPCPO,'NVP',NVPO)
      MAXNVP=NVPO(2)
      ALLOCATE(JDEBAR(MAXNVP+1),JARBVA(MAXNVP))
      CALL LCMGET(KPCPO,'NCALS',NCALS)
      CALL LCMGET(KPCPO,'DEBARB',JDEBAR)
      CALL LCMGET(KPCPO,'ARBVAL',JARBVA)
      I0=0
      DO 30 I=NVPO(1)-NCALS+1,NVPO(1)
      IF(JDEBAR(I+1).EQ.ICAL) THEN
         I0=I
         GO TO 40
      ENDIF
   30 CONTINUE
      CALL XABORT('COMSDB: MUPLET ALGORITHM FAILURE 1.')
   40 MUPCPO(NPAR+NLOC)=JARBVA(I0)
      DO 65 IPAR=NPAR+NLOC-1,1,-1
      I0=0
      DO 50 I=1,NVPO(1)-NCALS
      IF(JDEBAR(I+1).GT.I0) THEN
         I0=I
         GO TO 60
      ENDIF
   50 CONTINUE
      CALL XABORT('COMSDB: MUPLET ALGORITHM FAILURE 2.')
   60 MUPCPO(IPAR)=JARBVA(I0)
   65 CONTINUE
      DEALLOCATE(JARBVA,JDEBAR)
      CALL LCMSIX(KPCPO,' ',2)
      WRITE(6,*)'ICAL #',ICAL,': ',(MUPCPO(JM),JM=1,NPAR+NLOC)
      IF(IMPX.LE.1) GOTO 170
*----
*  RECOVER THE GLOBAL PARAMETERS
*----
         CALL LCMSIX(IPCPO,'GLOBAL',1)
         CALL LCMGTC(IPCPO,'PARKEY',12,NPAR,PARCPO)
         CALL LCMGTC(IPCPO,'PARFMT',8,NPAR,PARFMT)
         CALL LCMGET(IPCPO,'NVALUE',NVALUE)
         DO 100 IPAR=1,NPAR
            WRITE(RECNAM,'(''pval'',I8.8)') IPAR
            IVAL=MUPCPO(IPAR)
            IF(PARFMT(IPAR).EQ.'REAL') THEN
               ALLOCATE(VREAL(NVALUE(IPAR)))
               CALL LCMGET(IPCPO,RECNAM,VREAL)
               FLOTT=VREAL(IVAL)
               DEALLOCATE(VREAL)
               write(6,*)'IPAR ',IPAR,'->',PARCPO(IPAR),FLOTT
            ELSE IF(PARFMT(IPAR).EQ.'INTEGER') THEN
               ALLOCATE(VINTE(NVALUE(IPAR)))
               CALL LCMGET(IPCPO,RECNAM,VINTE)
               NITMA=VINTE(IVAL)
               DEALLOCATE(VINTE)
               write(6,*)'IPAR ',IPAR,'->',PARCPO(IPAR),NITMA
            ELSE IF(PARFMT(IPAR).EQ.'STRING') THEN
               IF(NVALUE(IPAR).GT.MAXVAL) CALL XABORT('COMSDB: MAXVAL '
     1         //'OVERFLOW.')
               CALL LCMGTC(IPCPO,RECNAM,12,NVALUE(IPAR),VCHAR)
               TEXT12=VCHAR(IVAL)
               write(6,*)'IPAR ',IPAR,'->',PARCPO(IPAR),TEXT12
            ENDIF
  100    CONTINUE
         CALL LCMSIX(IPCPO,' ',2)
*----
*  RECOVER THE LOCAL PARAMETERS
*----
      CALL LCMSIX(KPCPO,'TREE',1)
      DO 130 ILOC=1,NLOC
         WRITE(RECNAM,'(''pval'',I8.8)') ILOC
         IVAL=MUPCPO(NPAR+ILOC)
         CALL LCMLEN(KPCPO,RECNAM,ILONG,ITYLCM)
         ALLOCATE(VREAL(ILONG))
         CALL LCMGET(KPCPO,RECNAM,VREAL)
         FLOTT=VREAL(IVAL)
         DEALLOCATE(VREAL)
         WRITE(6,*)'ILOC ',ILOC,'->',PARKEY(IPAR),FLOTT
  130 CONTINUE
      CALL LCMSIX(KPCPO,' ',2)
  170 CONTINUE
* END OF LOOP ON CALCULATIONS. *******************************
      IF (IMPX.EQ.0) GOTO 190
      WRITE(6,*) 'Summary of the parameter included in the COMPO'
         CALL LCMSIX(IPCPO,'GLOBAL',1)
         CALL LCMGTC(IPCPO,'PARKEY',12,NPAR,PARCPO)
         CALL LCMGTC(IPCPO,'PARFMT',8,NPAR,PARFMT)
         CALL LCMGET(IPCPO,'NVALUE',NVALUE)
         DO 180 IPAR=1,NPAR
            WRITE(RECNAM,'(''pval'',I8.8)') IPAR
            IF(PARFMT(IPAR).EQ.'REAL') THEN
               ALLOCATE(VREAL(NVALUE(IPAR)))
               CALL LCMGET(IPCPO,RECNAM,VREAL)
               WRITE(6,*)'IPAR ',IPAR,'->',PARKEY(IPAR),
     1                   (VREAL(JM),JM=1,NVALUE(IPAR))
               DEALLOCATE(VREAL)
            ELSE IF(PARFMT(IPAR).EQ.'INTEGER') THEN
               ALLOCATE(VINTE(NVALUE(IPAR)))
               CALL LCMGET(IPCPO,RECNAM,VINTE)
               WRITE(6,*)'IPAR ',IPAR,'->',PARKEY(IPAR),
     1                   (VINTE(JM),JM=1,NVALUE(IPAR))
               DEALLOCATE(VINTE)
            ELSE IF(PARFMT(IPAR).EQ.'STRING') THEN
               IF(NVALUE(IPAR).GT.MAXVAL) CALL XABORT('COMSDB: MAXVAL '
     1         //'OVERFLOW.')
               CALL LCMGTC(IPCPO,RECNAM,12,NVALUE(IPAR),VCHAR)
               WRITE(6,*)'IPAR ',IPAR,'->',PARKEY(IPAR),
     1                   (VCHAR(JM),JM=1,NVALUE(IPAR))
            ENDIF
  180    CONTINUE
         CALL LCMSIX(IPCPO,' ',2)
      CALL LCMSIX(KPCPO,'TREE',1)
      DO 185 ILOC=1,NLOC
         CALL LCMGTC(IPCPO,'PARKEL',12,NLOC,PARCPO)
         CALL LCMGET(IPCPO,'NVALUE',NVALUE)
         WRITE(RECNAM,'(''pval'',I8.8)') ILOC
         CALL LCMLEN(KPCPO,RECNAM,ILONG,ITYLCM)
         ALLOCATE(VREAL(ILONG))
         CALL LCMGET(KPCPO,RECNAM,VREAL)
         WRITE(6,*)'ILOC ',ILOC,'->',PARCPO(ILOC),
     1             (VREAL(JM),JM=1,NVALUE(ILOC))
         DEALLOCATE(VREAL)
  185 CONTINUE
      CALL LCMSIX(KPCPO,' ',2)
      
  190 CONTINUE
* END OF LOOP ON MIXTURES. ***********************************
*
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MUPCPO)
      RETURN
      END
