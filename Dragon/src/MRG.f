*DECK MRG
      SUBROUTINE MRG(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*----------
*
*Purpose:
* Merge EXCELT or NXT geometry.
*
*Copyright:
* Copyright (C) 1997 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1) creation type(L_MC);
*         HENTRY(2) read-only or modification type(L_TRACK);
*         HENTRY(3) read-only type(L_LIBRARY) or type(L_MACROLIB).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*Comments: 
* LINKED LIST / XSM FILE:
*  HENTRY(1) : CREATION OR UPDATE MODE LINKED LIST TYPE(L_TRACK)
*              (CREATION only for EXCELT: type tracking).
*  HENTRY(2) : CREATION MODE SEQUENTIAL BINARY TRACKING FILE
*              (optionnal for NXT: type tracking).
*  HENTRY(3) : READ-ONLY LINKED LIST TYPE(L_TRACK)
*              (optionnal for NXT: type tracking).
*  HENTRY(4) : READ-ONLY SEQUENTIAL BINARY TRACKING FILE
*              (optionnal for NXT: type tracking).
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
      INTEGER          IOUT,NSTATE,NTC,NALB
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NSTATE=40,NTC=18,NALB=6,
     >                 NAMSBR='MRG   ')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
*      INTEGER          KENTRY(NENTRY)
      TYPE(C_PTR)      KENTRY(NENTRY)
      CHARACTER        HENTRY(NENTRY)*12
*----
*  Function for inport/export DDS
*----
      INTEGER          KDROPN,KDRCLS,IRC,IFILE
*----
*  LOCAL VARIABLES
*----
      INTEGER          IDTRKO,IDTRKN,IDTRKE,IDSTRO,IDSTRN,IEN,
     >                 IFTRKO,IFTRKN,IFTRKE
*      INTEGER         IPTRKN,IPTRKO
      TYPE(C_PTR)      IPTRKO,IPTRKN
      INTEGER          ISTATE(NSTATE),ISTATG(NSTATE)
      INTEGER          ITC
      CHARACTER        HSIGN*12
      INTEGER          NREGO,NUNO,NUNN,NSURO,NUNF,ITROP,
     >                 NSOUTO,NVOUTO,NSOUTN,NVOUTN,ILCMLN,ILCMTY,
     >                 NETSUR,NETVOL,NETNUO,NETNUN,NELT,MAXMN
      INTEGER          IPRINT,IUPD(4),NDIM,INDBC
      REAL             ALBEDN(NALB)
*----
*  Tracking file variables
*----
      INTEGER          IFMT,NCOMNT,NBTRK,NSCRP(8),NALBG,NANGL,MXSEG,
     >                 NUNKNO,NUNKNN,IOPTT,MAXMIX
      CHARACTER        CTRK*4,COMENT*80
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMERGE,MIXN,MATO,MATRTO,
     >                 MATN,KEYN,MATRTN,NEXMAT,NEXKEY,NEXMAN,NEXKEN
      REAL, ALLOCATABLE, DIMENSION(:)  :: VOLO,VOLN
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY .LT. 1 .OR. NENTRY .GT. 4) CALL XABORT(NAMSBR//
     >  ': From 1 to 4 data structures required')
*----
*  Find and validate structure types
*----
      IDTRKO=0
      IDTRKN=0
      IDTRKE=0
      IDSTRO=0
      IDSTRN=0
      IOPTT=0
      DO IEN=1,NENTRY
        IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2) THEN
          IF(JENTRY(IEN) .EQ. 0) THEN
            IDSTRN=IEN
          ELSE IF(JENTRY(IEN) .EQ. 1) THEN
            IDSTRN=IEN
            IDSTRO=IEN
          ELSE
            IDSTRO=IEN
          ENDIF
        ELSE IF(IENTRY(IEN) .EQ. 3) THEN
          IF(JENTRY(IEN) .EQ. 0) THEN
            IF(IDTRKN .EQ. 0) THEN
              IDTRKN=IEN
            ELSE
              IDTRKE=IEN
            ENDIF
          ELSE IF(JENTRY(IEN) .EQ. 2) THEN
            IDTRKO=IEN
          ENDIF
        ELSE
          CALL XABORT(NAMSBR//': One data structure has invalid type')
        ENDIF
      ENDDO
      IF(IDSTRO .EQ. 0) CALL XABORT(NAMSBR//
     >': Reference tracking data structure missing')
      IF(IDSTRN .NE. IDSTRO) THEN
*----
*  Make a full copy of the old data structure to the new data structure
*  we will update later
*----
        IPTRKO=KENTRY(IDSTRO)
        IPTRKN=IPTRKO
        IF(IDSTRN .NE. 0) THEN
          IPTRKN=KENTRY(IDSTRN)
          IFILE=KDROPN('DUMMYSQ',0,2,0,0)
          IF(IFILE.LE.0) CALL XABORT(NAMSBR//': KDROPN FAILURE.')
          CALL LCMEXP(IPTRKO,0,IFILE,1,1)
          REWIND(IFILE)
          CALL LCMEXP(IPTRKN,0,IFILE,1,2)
          IRC=KDRCLS(IFILE,2)
        ENDIF
      ELSE
        IPTRKN=KENTRY(IDSTRN)
      ENDIF
*----
*  Test contents of the data structure to update
*  Either old data structure or new data structure that now contains
*  a copy of the old data structure
*----
      CALL LCMGTC(IPTRKN,'SIGNATURE   ',12,1,HSIGN)
      IF(HSIGN .NE. 'L_TRACK') CALL XABORT(NAMSBR//
     >  ': SIGNATURE OF '//HENTRY(3)//' IS '//HSIGN//
     >  '. L_TRACK EXPECTED.')
      CALL LCMGTC(IPTRKN,'TRACK-TYPE  ',12,1,HSIGN)
      IF(HSIGN .NE. 'EXCELL') CALL XABORT(NAMSBR//
     >  ': ILLEGAL TRACKING FORMAT')
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPTRKN,'STATE-VECTOR',ISTATE)
      NREGO=ISTATE(1)
      NUNO=ISTATE(2)
      NUNF=NUNO/NREGO
      MAXMIX=ISTATE(4)
      NSURO=ISTATE(5)
      ITROP=ISTATE(7)
*----
*  Test if NXT: or EXCELT: tracking
*  For EXCELT: -> tracking files are required
*  For NXT:    -> tracking files are optional
*----
      IF(ITROP .GE. 1 .AND. ITROP .LE. 3) THEN
        IF(IDTRKO .EQ. 0 .OR. IDTRKN .EQ. 0) CALL XABORT(NAMSBR//
     >  ': Tracking files required for EXCELT: tracking')
        IOPTT=1
      ELSE IF(ITROP .EQ. 4) THEN
        IF(IDTRKO .EQ. 0 .AND. IDTRKN .EQ. 0) THEN
          IOPTT=2
        ELSE IF(IDTRKO .GT. 0 .AND. IDTRKN .GT. 0) THEN
          IOPTT=3
          IF(IDTRKE .GT. 0) IOPTT=4
        ELSE
          CALL XABORT(NAMSBR//
     >  ': Either 0, 2 or 3 tracking files required for NXT: tracking')
        ENDIF
      ELSE
        CALL XABORT(NAMSBR//
     >  ': Invalid tracking options on tracking data structure')
      ENDIF
      NVOUTO=NREGO
      NSOUTO=NSURO
*----
*  READ MERGE INFORMATION
*----
      ALLOCATE(IMERGE(-NSOUTO:NVOUTO),MIXN(NVOUTO))
      CALL MRGGET(IPRINT,NSOUTO,NVOUTO,NSOUTN,NVOUTN,
     >            IUPD,IMERGE,MIXN,ALBEDN)
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      IF(IOPTT .NE. 4) THEN
        NUNN=NUNF*NVOUTN
        ISTATE(1)=NVOUTN
        ISTATE(2)=NUNN
        ISTATE(5)=NSOUTN
*----
*  Read global records to merge
*----
        ALLOCATE(MATN(NVOUTN),VOLN(NVOUTN),MATRTN(NSOUTN),KEYN(NREGO))
        ALLOCATE(MATO(NVOUTO),VOLO(NVOUTO),MATRTO(NSOUTO))
        CALL LCMGET(IPTRKN,'MATCOD      ',MATO)
        CALL LCMGET(IPTRKN,'VOLUME      ',VOLO)
        CALL LCMLEN(IPTRKN,'BC-REFL+TRAN',ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. NSOUTO) THEN
          INDBC=1
          CALL LCMGET(IPTRKN,'BC-REFL+TRAN',MATRTO)
        ELSE
          INDBC=0
          CALL XDISET(MATRTO,NSOUTO,0)
        ENDIF
        NETNUN=NSOUTN+NVOUTN+1
        IF(ITROP .EQ. 4) THEN
*----
*  Process NXT: Records
*----
          CALL LCMSIX(IPTRKN,'NXTRecords  ',1)
          CALL XDISET(ISTATG,NSTATE,0)
          CALL LCMGET(IPTRKN,'G00000001DIM',ISTATG)
          NETSUR=ISTATG(22)
          NETVOL=ISTATG(23)
          NELT=NETSUR+NETVOL+1
          ALLOCATE(NEXMAT(NELT),NEXKEY(NELT))
          CALL LCMGET(IPTRKN,'MATALB      ',NEXMAT)
          CALL LCMGET(IPTRKN,'KEYMRG      ',NEXKEY)
          CALL MRGVON(IUPD  ,NSOUTO,NVOUTO,NSOUTN,NVOUTN,
     >                NETSUR,NETVOL,NUNN  ,MAXMN ,
     >                IMERGE,MATO  ,VOLO  ,MATRTO,
     >                MATN  ,VOLN  ,KEYN  ,MATRTN,
     >                NEXMAT,NEXKEY)
*----
*  Save NXT: specific records
*----
          CALL LCMPUT(IPTRKN,'MATALB      ',NELT,1,NEXMAT)
          CALL LCMPUT(IPTRKN,'KEYMRG      ',NELT,1,NEXKEY)
          CALL LCMSIX(IPTRKN,'NXTRecords  ',2)
          DEALLOCATE(NEXKEY,NEXMAT)
        ELSE
*----
*  Process EXCELT: Records
*----
          CALL LCMSIX(IPTRKN,'EXCELL      ',1)
          CALL XDISET(ISTATG,NSTATE,0)
          CALL LCMGET(IPTRKN,'STATE-VECTOR',ISTATG)
          NETSUR=ISTATG(2)
          NETVOL=ISTATG(3)
          NETNUO=ISTATG(6)
          ALLOCATE(NEXMAT(-NETSUR:NETVOL),NEXKEY(-NETSUR:NETVOL))
          ALLOCATE(NEXMAN(-NSOUTN:NVOUTN),NEXKEN(-NSOUTN:NVOUTN))
          CALL LCMGET(IPTRKO,'MATALB      ',NEXMAT)
          CALL LCMGET(IPTRKO,'KEYMRG      ',NEXKEY)
          CALL MRGVOL(IUPD  ,NSOUTO,NVOUTO,NSOUTN,NVOUTN,NREGO ,
     >                IMERGE,MIXN  ,MATO  ,VOLO  ,MATN  ,VOLN ,
     >                KEYN  ,MATRTO,MATRTN,MAXMN ,NETVOL,NETSUR,
     >                NEXMAT,NEXKEY,NEXMAN,NEXKEN)
*----
*  Save EXCELT: specific records
*----
          CALL LCMPUT(IPTRKN,'MATALB      ',NETNUN,1,NEXMAN)
          CALL LCMPUT(IPTRKN,'KEYMRG      ',NETNUN,1,NEXKEN)
          CALL LCMPUT(IPTRKN,'STATE-VECTOR',NSTATE,1,ISTATG)
          CALL LCMSIX(IPTRKN,'EXCELL      ',2)
          DEALLOCATE(NEXKEN,NEXMAN)
          DEALLOCATE(NEXKEY,NEXMAT)
        ENDIF
        ISTATE(1)=NVOUTN
        ISTATE(2)=NUNN
        ISTATE(4)=MAXMN
        ISTATE(5)=NSOUTN
*----
*  Save global tracking records
*----
        CALL LCMPUT(IPTRKN,'STATE-VECTOR',NSTATE,1,ISTATE)
        IF(IUPD(4). GT. 0 ) THEN
          CALL LCMPUT(IPTRKN,'ALBEDO      ',NALB,2,ALBEDN)
        ENDIF
        CALL LCMPUT(IPTRKN,'MATCOD      ',NVOUTN,1,MATN)
        CALL LCMPUT(IPTRKN,'VOLUME      ',NVOUTN,2,VOLN)
        IF(IUPD(1) .GT. 0) THEN
          CALL LCMPUT(IPTRKN,'KEYFLX      ',NVOUTN,1,KEYN)
        ENDIF
        CALL LCMPUT(IPTRKN,'BC-REFL+TRAN',NSOUTN,1,MATRTN)
        DEALLOCATE(MATRTO,VOLO,MATO)
        DEALLOCATE(MATN,VOLN,KEYN,MATRTN)
*----
* Processing of tracking data structure finished
* Now process tracking file if required
*----
        IF(IDTRKN .GT. 0) THEN
          IFMT=0
          IFTRKN=FILUNIT(KENTRY(IDTRKN))
          IFTRKO=FILUNIT(KENTRY(IDTRKO))
          READ (IFTRKO) CTRK,NCOMNT,NBTRK
          WRITE(IFTRKN) CTRK,NCOMNT,NBTRK
          DO ITC= 1, NCOMNT
            READ (IFTRKO) COMENT
            IF(COMENT .EQ. 'OPTION      : Extended    ') IFMT=1
            WRITE(IFTRKN) COMENT
          ENDDO
          READ (IFTRKO) (NSCRP(ITC),ITC=1,8)
          NDIM  =NSCRP(1)
          IF(NVOUTO .NE. NSCRP(3)) CALL XABORT(NAMSBR//
     >    ': Number of regions on tracking file inconsistent with '//
     >    'that on tracking data structure')
          IF(NSOUTO .NE. NSCRP(4)) CALL XABORT(NAMSBR//
     >    ': Number of surfaces on tracking file inconsistent with '//
     >    'that on tracking data structure')
          NALBG= NSCRP(5)
          NANGL= NSCRP(7)
          MXSEG= NSCRP(8)
          NSCRP(3)=NVOUTN
          NSCRP(4)=NSOUTN
          WRITE(IFTRKN) (NSCRP(ITC),ITC=1,8)
          NUNKNO=NSOUTO+NVOUTO
          CALL MRGVST(IFTRKO,IFTRKN,IPRINT,IUPD  ,NDIM  ,NALBG,NANGL,
     >                NSOUTO,NVOUTO,NSOUTN,NVOUTN,IMERGE,MIXN  )
          NUNKNN=NSOUTN+NVOUTN
*----
*  TRACKING LINE
*----
          CALL MRGLIN(IPRINT,IFTRKO,NSOUTO,NVOUTO,IFTRKN,
     >                IMERGE,MXSEG)
        ENDIF
      ELSE
*----
*  SPLIT TRACKING FILE
*---
        IF(IUPD(1) .GE. 0) CALL XABORT(NAMSBR//
     >  ': No region specified for EXTR')
        IFTRKN=FILUNIT(KENTRY(IDTRKN))
        IFTRKO=FILUNIT(KENTRY(IDTRKO))
        IFTRKE=FILUNIT(KENTRY(IDTRKE))
        READ (IFTRKO) CTRK,NCOMNT,NBTRK
        WRITE(IFTRKN) CTRK,NCOMNT,NBTRK
        WRITE(IFTRKE) CTRK,NCOMNT,NBTRK
        DO ITC= 1, NCOMNT
          READ (IFTRKO) COMENT
          IF(COMENT .EQ. 'OPTION      : Extended    ') IFMT=1
          WRITE(IFTRKN) COMENT
          WRITE(IFTRKE) COMENT
        ENDDO
        READ (IFTRKO) (NSCRP(ITC),ITC=1,8)
        NDIM  =NSCRP(1)
        IF(NVOUTO .NE. NSCRP(3)) CALL XABORT(NAMSBR//
     >  ': Number of regions on tracking file inconsistent with '//
     >  'that on tracking data structure')
        IF(NSOUTO .NE. NSCRP(4)) CALL XABORT(NAMSBR//
     >  ': Number of surfaces on tracking file inconsistent with '//
     >  'that on tracking data structure')
        NALBG= NSCRP(5)
        NANGL= NSCRP(7)
        MXSEG= NSCRP(8)
        WRITE(IFTRKN) (NSCRP(ITC),ITC=1,8)
        WRITE(IFTRKE) (NSCRP(ITC),ITC=1,8)
        NUNKNO=NSOUTO+NVOUTO
        CALL MRGXTC(IFTRKO,IFTRKN,IFTRKE,IPRINT,IUPD ,NDIM,
     >              NALBG ,NANGL ,NSOUTO,NVOUTO,MXSEG,IMERGE)
      ENDIF
      DEALLOCATE(MIXN,IMERGE)
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
