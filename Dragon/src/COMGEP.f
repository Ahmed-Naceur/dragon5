*DECK COMGEP
      SUBROUTINE COMGEP(IPCPO,IPDEPL,IPLB1,IPLB2,IPEDIT,IMPX,ITIM,NORIG,
     1 NPAR,NLOC,NMIL,MUPLET,LGNEW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover remaining global and local parameters. Update the parameter
* tree (in each mixture) for a new elementary calculation.
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
* IPCPO   pointer to the multicompo.
* IPDEPL  pointer to the burnup object.
* IPLB1   pointer to the first microlib object.
* IPLB2   pointer to the second (optional) microlib object.
* IPEDIT  pointer to the edition object.
* IMPX    print parameter.
* ITIM    index of the current burnup step.
* NORIG   index of the elementary calculation associated to the
*         father node in the parameter tree.
* NPAR    number of global parameters.
* NLOC    number of local parameters.
* NMIL    number of homogenized mixtures.
* MUPLET  tuple of indices associated to each global parameter of the
*         elementary calculation.
* LGNEW   parameter modification flag (=.true. if the I-th global
*         parameter has changed in the new elementary calculation).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPCPO,IPDEPL,IPLB1,IPLB2,IPEDIT
      INTEGER IMPX,ITIM,NORIG,NPAR,NLOC,NMIL,MUPLET(NPAR)
      LOGICAL LGNEW(NPAR)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXPAR=50)
      TYPE(C_PTR) JPCPO,KPCPO,IPLB3
      INTEGER ISTATE(NSTATE),PARMIL(MAXPAR),PARCAD(MAXPAR+1),
     1 PARPAD(MAXPAR+1),NVPO(2)
      CHARACTER PARKEY(MAXPAR)*12,PARCHR(MAXPAR)*8,PARTYP(MAXPAR)*4,
     1 PARBIB(MAXPAR)*12,TEXT4*4,TEXT8*8,TEXT12*12,NAMLCM*12,NAMMY*12
      LOGICAL LGERR,EMPTY,LCM,COMTRE,LAST
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDEBAR,IARBVA,JDEBAR,JARBVA,
     1 IORIGI
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MUPL2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RVALO
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LGNE2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MUPL2(NPAR+NLOC,NMIL),LGNE2(NPAR+NLOC,NMIL))
*----
*  RECOVER INFORMATION FROM THE 'STATE-VECTOR' PARAMETER LIST
*----
      CALL LCMGET(IPCPO,'STATE-VECTOR',ISTATE)
      IF(NPAR.NE.ISTATE(5)) CALL XABORT('COMGEP: WRONG VALUE OF NPAR.')
      IF(NLOC.NE.ISTATE(6)) CALL XABORT('COMGEP: WRONG VALUE OF NLOC.')
      NMIL=ISTATE(1)
      NG=ISTATE(2)
      NCAL=ISTATE(3)
      MAXCAL=ISTATE(4)
      NPCHR=ISTATE(7)
      NPPNT=ISTATE(8)
      NPCHRL=ISTATE(9)
*----
*  RECOVER INFORMATION FROM THE 'GLOBAL' DIRECTORY
*----
      CALL LCMSIX(IPCPO,'GLOBAL',1)
      IF(NPAR.GT.0) CALL LCMGTC(IPCPO,'PARKEY',12,NPAR,PARKEY)
      IF(NPAR.GT.0) CALL LCMGTC(IPCPO,'PARTYP',4,NPAR,PARTYP)
      CALL LCMGET(IPCPO,'PARCAD',PARCAD)
      CALL LCMGET(IPCPO,'PARPAD',PARPAD)
      IF(NPCHR.GT.0) CALL LCMGTC(IPCPO,'PARCHR',8,NPCHR,PARCHR)
      IF(NPPNT.GT.0) CALL LCMGET(IPCPO,'PARMIL',PARMIL)
      IF(NPPNT.GT.0) CALL LCMGTC(IPCPO,'PARBIB',12,NPPNT,PARBIB)
      CALL LCMSIX(IPCPO,' ',2)
*----
*  RECOVER REMAINING GLOBAL PARAMETERS
*----
      DO 10 IPAR=1,NPAR
      IF(PARTYP(IPAR).EQ.'VALU') THEN
         GO TO 10
      ELSE IF(LGNEW(IPAR)) THEN
        CALL XABORT('COMGEP: PARAMETER '//PARTYP(IPAR)//' IS ALREADY D'
     1  //'EFINED (REMOVE THE '//PARKEY(IPAR)//' KEYWORD).')
      ELSE IF((PARTYP(IPAR).EQ.'IRRA').OR.(PARTYP(IPAR).EQ.'TIME').OR.
     1        (PARTYP(IPAR).EQ.'POWR').OR.(PARTYP(IPAR).EQ.'FLUB').OR.
     2        (PARTYP(IPAR).EQ.'FLUX').OR.(PARTYP(IPAR).EQ.'MASL')) THEN
*
*        RECOVER GLOBAL PARAMETER VALUES FROM THE BURNUP OBJECT.
         IF(.NOT.C_ASSOCIATED(IPDEPL)) THEN
            CALL XABORT('COMGEP: NO DEPLETION OBJECT AVAILABLE AMONG T'
     1      //'HE RHS LCM OBJECTS.')
         ENDIF
         CALL LCMGET(IPDEPL,'STATE-VECTOR',ISTATE)
         NBURN=ISTATE(3)
         NBISO=ISTATE(4)
         NREAC=ISTATE(6)
         NVAR=ISTATE(7)
         NBMIX=ISTATE(8)
         CALL COMGEM(IPDEPL,ITIM,PARTYP(IPAR),0,NBURN,NBMIX,NBISO,
     1   NREAC,NVAR,VALPAR)
      ELSE IF((PARTYP(IPAR).EQ.'TEMP').OR.(PARTYP(IPAR).EQ.'CONC'))
     1   THEN
*
*        RECOVER GLOBAL PARAMETER VALUES FROM A MICROLIB OBJECT.
         IF(.NOT.C_ASSOCIATED(IPLB1)) THEN
            CALL XABORT('COMGEP: MICROLIB EXPECTED AT RHS.')
         ENDIF
         IPCAD=PARCAD(IPAR+1)-PARCAD(IPAR)
         IPPAD=PARPAD(IPAR+1)-PARPAD(IPAR)
         IF(IPCAD.EQ.1) IPCAD=PARCAD(IPAR+1)-PARCAD(1)
         IF(IPPAD.EQ.1) IPPAD=PARPAD(IPAR+1)-PARPAD(1)
         TEXT8=' '
         TEXT12=' '
         IMILI=0
         IF(IPCAD.GT.0) TEXT8=PARCHR(IPCAD)
         IF(IPPAD.GT.0) TEXT12=PARBIB(IPPAD)
         IF(IPPAD.GT.0) IMILI=PARMIL(IPPAD)
         CALL LCMGET(IPLB1,'STATE-VECTOR',ISTATE)
         MAXNBI=ISTATE(2)
         IF(C_ASSOCIATED(IPLB2)) THEN
            CALL LCMGET(IPLB2,'STATE-VECTOR',ISTATE)
            MAXNBI=MAX(MAXNBI,ISTATE(2))
         ENDIF
         CALL COMBIB(IPLB1,IPLB2,PARTYP(IPAR),IMILI,TEXT12,TEXT8,MAXNBI,
     1   VALPAR)
      ELSE
         CALL XABORT('COMGEP: '//PARTYP(IPAR)//' IS AN UNKNOWN PARAM'//
     1   'ETER TYPE.')
      ENDIF
      IF(IMPX.GT.0) WRITE(6,100) PARKEY(IPAR),VALPAR
*
      CALL LCMSIX(IPCPO,'GLOBAL',1)
      TEXT8='REAL'
      CALL COMPAV(IPCPO,IPAR,NPAR,TEXT8,VALPAR,NITMA,TEXT12,
     1 MUPLET(IPAR),LGNEW(IPAR))
      CALL LCMSIX(IPCPO,' ',2)
   10 CONTINUE
      DO 25 IBM=1,NMIL
      DO 20 IPAR=1,NPAR
      MUPL2(IPAR,IBM)=MUPLET(IPAR)
      LGNE2(IPAR,IBM)=LGNEW(IPAR)
   20 CONTINUE
   25 CONTINUE
      IF(NLOC.EQ.0) GO TO 50
*----
*  RECOVER INFORMATION FROM THE 'LOCAL' DIRECTORY
*----
      CALL LCMSIX(IPCPO,'LOCAL',1)
      CALL LCMGTC(IPCPO,'PARKEY',12,NLOC,PARKEY)
      CALL LCMGTC(IPCPO,'PARTYP',4,NLOC,PARTYP)
      CALL LCMGET(IPCPO,'PARCAD',PARCAD)
      IF(NPCHRL.GT.0) CALL LCMGTC(IPCPO,'PARCHR',8,NPCHRL,PARCHR)
      CALL LCMSIX(IPCPO,' ',2)
*----
*  RECOVER LOCAL PARAMETERS
*----
      CALL LCMGTC(IPEDIT,'LAST-EDIT',12,1,TEXT12)
      ALLOCATE(RVALO(NLOC,NMIL))
      DO 45 IPAR=1,NLOC
      IF((PARTYP(IPAR).EQ.'IRRA').OR.(PARTYP(IPAR).EQ.'TIME').OR.
     1   (PARTYP(IPAR).EQ.'POWR').OR.(PARTYP(IPAR).EQ.'FLUG').OR.
     2   (PARTYP(IPAR).EQ.'FLUB').OR.(PARTYP(IPAR).EQ.'FLUX').OR.
     3   (PARTYP(IPAR).EQ.'MASL')) THEN
*
*        RECOVER LOCAL PARAMETERS FROM THE BURNUP OBJECT.
         CALL LCMGET(IPDEPL,'STATE-VECTOR',ISTATE)
         NBURN=ISTATE(3)
         NBISO=ISTATE(4)
         NREAC=ISTATE(6)
         NVAR=ISTATE(7)
         NBMIX=ISTATE(8)
         IF(.NOT.C_ASSOCIATED(IPDEPL)) THEN
            CALL XABORT('COMGEP: NO DEPLETION OBJECT AVAILABLE AMONG T'
     1      //'HE RHS LCM OBJECTS.')
         ENDIF
         CALL LCMGET(IPEDIT,'STATE-VECTOR',ISTATE)
         NREG=ISTATE(17)
         CALL COMGEN(IPDEPL,IPEDIT,NREG,NMIL,ITIM,PARTYP(IPAR),NBURN,
     1   NBMIX,NBISO,NREAC,NVAR,IPAR,NLOC,RVALO)
      ELSE IF((PARTYP(IPAR).EQ.'TEMP').OR.(PARTYP(IPAR).EQ.'CONC'))
     1   THEN
*
*        RECOVER LOCAL PARAMETERS FROM THE MICROLIB IN EDIT OBJECT.
         IPCAD=PARCAD(IPAR+1)-PARCAD(IPAR)
         IF(IPCAD.EQ.1) IPCAD=PARCAD(IPAR+1)-PARCAD(1)
         TEXT8=' '
         IF(IPCAD.GT.0) TEXT8=PARCHR(IPCAD)
         CALL LCMSIX(IPEDIT,TEXT12,1)
         CALL LCMGET(IPEDIT,'STATE-VECTOR',ISTATE)
         MAXNBI=ISTATE(2)
         CALL LCMINF(IPEDIT,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
         IPLB3=C_NULL_PTR
         DO 30 IBM=1,NMIL
         CALL COMBIB(IPEDIT,IPLB3,PARTYP(IPAR),IBM,NAMLCM,TEXT8,MAXNBI,
     1   VALPAR)
         RVALO(IPAR,IBM)=VALPAR
   30    CONTINUE
         CALL LCMSIX(IPEDIT,' ',2)
      ELSE
         CALL XABORT('COMGEP: '//PARTYP(IPAR)//' IS AN UNKNOWN LOCAL'//
     1   ' PARAMETER TYPE.')
      ENDIF
      IF(IMPX.GT.1) THEN
         WRITE(6,120) PARKEY(IPAR),(RVALO(IPAR,IBM),IBM=1,NMIL)
      ENDIF
      JPCPO=LCMLID(IPCPO,'MIXTURES',NMIL)
      DO 40 IBM=1,NMIL
      KPCPO=LCMDIL(JPCPO,IBM)
      CALL LCMSIX(KPCPO,'TREE',1)
      FLOTT=RVALO(IPAR,IBM)
      TEXT8='REAL'
      CALL COMPAV(KPCPO,IPAR,NLOC,TEXT8,FLOTT,NITMA,TEXT12,
     1 MUPL2(NPAR+IPAR,IBM),LGNE2(NPAR+IPAR,IBM))
      CALL LCMSIX(KPCPO,' ',2)
   40 CONTINUE
   45 CONTINUE
      DEALLOCATE(RVALO)
*----
*  UPDATE THE PARAMETER TREE IN EACH MIXTURE COMPONENT
*----
   50 JPCPO=LCMLID(IPCPO,'MIXTURES',NMIL)
      DO 90 IBM=1,NMIL
      IF(IMPX.GT.4) THEN
         WRITE(6,'(/17H COMGEP: MIXTURE=,I6)') IBM
         WRITE(6,110) (MUPL2(I,IBM),I=1,NPAR+NLOC)
         WRITE(6,'(/)')
      ENDIF
      DO 55 I=1,NPAR+NLOC
      IF(MUPL2(I,IBM).EQ.0) THEN
         WRITE(6,'(/17H COMGEP: MIXTURE=,I6,23H, UNDEFINED MUPLET ELEM,
     1   4HENT=,I6)') IBM,I
         CALL XABORT('COMGEP: A MUPLET ELEMENT IS NOT ASSIGNED.')
      ENDIF
   55 CONTINUE
      KPCPO=LCMDIL(JPCPO,IBM)
**
** Parameter tree: this tree has a number of stages equal to the
** number of parameters. For each value of the i-th parameter, we
** find the position in the tree corresponding to the value of the
** (i+1)-th parameter.
** NCALAR  Number of elementary calculations stored in the tree.
** NVPO(1) Number of nodes in the parameter tree, including the root.
**         The value corresponding to the root is not used.
** DEBARB  - If the node does not correspond to the last parameter:
**           index in DEBARB of the first daughter of the node.
**         - If the node correspond to the last parameter: index in
**           DEBARB where we recover the index of an elementary
**           calculation.
** ARBVAL  Index of the corresponding parameter in the 'pval'//n
**         record.
*
**     EXEMPLE:   dn = value in DEBARB,  (m) = value in ARBVAL
**
**     Root                          *(0)
**                                     !
**     Param. Nb 1                  d2(1)
**                            -------------------
**                           !                   !
**     Param. Nb 2        d3(1)                4(2)
**                       ---------           ---------
**                      !         !         !    !    !
**     Param. Nb 3   d5(1)      6(3)     d7(1) 8(2) 9(3)   d10
**
**     Calculation Nb:  4         5         1    2    3
**
**     DEBARB:      2  3  5  7 10  4  5  1  2  3
**     ARBVAL:      0  1  1  2  1  3  1  2  3
*
      CALL LCMSIX(KPCPO,'TREE',1)
      CALL LCMLEN(KPCPO,'NVP',ILONG,ITYLCM)
      IF(ILONG.EQ.0) THEN
         MAXNVP=20*(NPAR+NLOC+1)
         ALLOCATE(IDEBAR(MAXNVP+1),IARBVA(MAXNVP))
         CALL XDISET(IDEBAR,MAXNVP+1,0)
         CALL XDISET(IARBVA,MAXNVP,0)
         IARBVA=0
         DO 60 I=1,NPAR+NLOC
         IDEBAR(I)=I+1
         IARBVA(I+1)=1
   60    CONTINUE
         IDEBAR(NPAR+NLOC+1)=NPAR+NLOC+2
         IDEBAR(NPAR+NLOC+2)=1
         NCALAR=1
         NVPNEW=NPAR+NLOC+1
      ELSE
*
*        Find position of the new point and create new PARBRE.
*
*        "II" is the order number of first parameter which recives a
*        "brand new" value.
*        COMTRE returns .TRUE. if the sweep throught the tree reaches
*        its bottom, otherwise it returns "KK" value: level of the
*        first new node to be introduced.
*
         CALL LCMGET(KPCPO,'NVP',NVPO)
         MAXNVP=NVPO(2)
         ALLOCATE(JDEBAR(NVPO(1)+1),JARBVA(NVPO(1)))
         CALL LCMGET(KPCPO,'NCALS',NCALAR)
         CALL LCMGET(KPCPO,'DEBARB',JDEBAR)
         CALL LCMGET(KPCPO,'ARBVAL',JARBVA)
         DO 70 IPAR=1,NPAR+NLOC
         IF(LGNE2(IPAR,IBM)) THEN
            II=IPAR
            GO TO 80
         ENDIF
   70    CONTINUE
         II=NPAR+NLOC+1
   80    LGERR=COMTRE(NPAR+NLOC,NVPO(1),JARBVA,JDEBAR,
     1   MUPL2(1,IBM),KK,I0,IORD,JJ,LAST)
         IF((II.GT.NPAR+NLOC).AND.LGERR) THEN
            WRITE(TEXT4,'(I4)') IORD
            CALL XABORT('COMGEP: ELEMENTARY CALCULATION HAS THE SAME'//
     1      ' PARAMETERS AS ELEMENTARY CALCULATION NB '//TEXT4)
         ENDIF
*
*        Size of the new tree.
*
         NVPNEW=NVPO(1)+NPAR+NLOC+1-MIN(II,KK)
         IF(NVPNEW.GT.MAXNVP) MAXNVP=NVPNEW+MAXNVP
         ALLOCATE(IDEBAR(MAXNVP+1),IARBVA(MAXNVP))
         CALL XDISET(IDEBAR(NVPNEW+2),MAXNVP-NVPNEW,0)
         CALL XDISET(IARBVA(NVPNEW+1),MAXNVP-NVPNEW,0)
*
*        Update values and suppress old PARBRE.
*
         CALL COMARB(NPAR+NLOC,NVPO(1),NVPNEW,JDEBAR,JARBVA,
     1   LGNE2(1,IBM),MUPL2(1,IBM),NCALAR,IDEBAR,IARBVA)
         DEALLOCATE(JARBVA,JDEBAR)
      ENDIF
      IF(NCALAR.NE.NCAL) CALL XABORT('COMGEP: INVALID NCALAR.')
      NVPO(1)=NVPNEW
      NVPO(2)=MAXNVP
      CALL LCMPUT(KPCPO,'NVP',2,1,NVPO)
      CALL LCMPUT(KPCPO,'NCALS',1,1,NCALAR)
      CALL LCMPUT(KPCPO,'DEBARB',NVPNEW+1,1,IDEBAR)
      CALL LCMPUT(KPCPO,'ARBVAL',NVPNEW,1,IARBVA)
      DEALLOCATE(IARBVA,IDEBAR)
      IF(NCALAR.EQ.1) THEN
         ALLOCATE(IORIGI(MAXCAL))
         CALL XDISET(IORIGI,MAXCAL,0)
      ELSE
         CALL LCMLEN(KPCPO,'ORIGIN',MAXOLD,ITYLCM)
         IF(MAXOLD.GT.MAXCAL) CALL XABORT('COMGEP: ORIGIN OVERFLOW(1).')
         ALLOCATE(IORIGI(MAXCAL))
         CALL XDISET(IORIGI,MAXCAL,0)
         CALL LCMGET(KPCPO,'ORIGIN',IORIGI)
      ENDIF
      IF(NCALAR.GT.MAXCAL) CALL XABORT('COMGEP: ORIGIN OVERFLOW(2).')
      IORIGI(NCALAR)=NORIG
      CALL LCMPUT(KPCPO,'ORIGIN',MAXCAL,1,IORIGI)
      DEALLOCATE(IORIGI)
      CALL LCMSIX(KPCPO,' ',2)
   90 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LGNE2,MUPL2)
      RETURN
*
  100 FORMAT(31H COMGEP: SET GLOBAL PARAMETER ',A,3H' =,1P,E12.4)
  110 FORMAT(17H COMGEP: MUPLET =,10I6:/(17X,10I6))
  120 FORMAT(30H COMGEP: SET LOCAL PARAMETER ',A,3H' =,1P,5E12.4/(37X,
     1 5E12.4))
      END
