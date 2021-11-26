*DECK SAPGEP
      SUBROUTINE SAPGEP(IPSAP,IPDEPL,IPLB1,IPLB2,IPEDIT,IMPX,ITIM,NORIG,
     1 NPAR,MUPLET,LGNEW,NVPNEW,NCALAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover remaining global parameters and local values. Update the
* parameter tree for a new elementary calculation.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSAP   pointer to the Saphyb.
* IPDEPL  pointer to the burnup object.
* IPLB1   pointer to the first microlib object.
* IPLB2   pointer to the second (optional) microlib object.
* IPEDIT  pointer to the edition object.
* IMPX    print parameter.
* ITIM    index of the current burnup step.
* NORIG   index of the elementary calculation associated to the
*         father node in the parameter tree.
* NPAR    number of global parameters.
* MUPLET  tuple of indices associated to each global parameter of the
*         elementary calculation.
* LGNEW   parameter modification flag (.TRUE. only if the I-th global
*         parameter has changed in the new elementary calculation).
*
*Parameters: output
* NVPNEW  number of nodes in the global parameter tree.
* NCALAR  index of the new elementary calculation.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP,IPDEPL,IPLB1,IPLB2,IPEDIT
      INTEGER IMPX,ITIM,NORIG,NPAR,MUPLET(NPAR),NVPNEW,NCALAR
      LOGICAL LGNEW(NPAR)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPLB3
      PARAMETER (NDIMSA=50,MAXPAR=50)
      INTEGER IDATA(NDIMSA),PARMIL(MAXPAR),
     1 PARCAD(MAXPAR+1),PARPAD(MAXPAR+1),LOCADR(MAXPAR+1)
      CHARACTER PARKEY(MAXPAR)*4,PARCHR(MAXPAR)*8,PARTYP(MAXPAR)*4,
     1 PARFMT(MAXPAR)*8,PARBIB(MAXPAR)*12,PARNAM(MAXPAR)*80,TEXT4*4,
     2 TEXT8*8,TEXT12*12,NAMLCM*12,NAMMY*12,HSMG*131
      LOGICAL LGERR,EMPTY,LCM,COMTRE,LAST
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDEBAR,IARBVA,JDEBAR,JARBVA,
     1 IORIGI
      REAL, ALLOCATABLE, DIMENSION(:) :: RVALO
*----
*  RECOVER INFORMATION FROM THE 'DIMSAP' PARAMETER LIST.
*----
      NVPNEW=0
      CALL LCMGET(IPSAP,'DIMSAP',IDATA)
      IF(NPAR.NE.IDATA(8)) CALL XABORT('SAPGEP: WRONG VALUE OF NPAR.')
      NMIL=IDATA(7)
      NPCHR=IDATA(9)
      NPPNT=IDATA(10)
      NLOC=IDATA(11)
      NPCHRL=IDATA(12)
      NPPNTL=IDATA(13)
      NVPO=IDATA(17)
      NCALAR=IDATA(19)
      NG=IDATA(20)
*----
*  RECOVER INFORMATION FROM THE 'paramdescrip' DIRECTORY.
*----
      IF(NPAR.EQ.0) GO TO 45
      CALL LCMSIX(IPSAP,'paramdescrip',1)
      CALL LCMGTC(IPSAP,'PARKEY',4,NPAR,PARKEY)
      CALL LCMGTC(IPSAP,'PARTYP',4,NPAR,PARTYP)
      CALL LCMGTC(IPSAP,'PARFMT',8,NPAR,PARFMT)
      CALL LCMGET(IPSAP,'PARCAD',PARCAD)
      CALL LCMGET(IPSAP,'PARPAD',PARPAD)
      IF(NPCHR.GT.0) CALL LCMGTC(IPSAP,'PARCHR',8,NPCHR,PARCHR)
      IF(NPPNT.GT.0) CALL LCMGET(IPSAP,'PARMIL',PARMIL)
      IF(NPPNT.GT.0) CALL LCMGTC(IPSAP,'PARBIB',12,NPPNT,PARBIB)
      CALL LCMSIX(IPSAP,' ',2)
*----
*  RECOVER REMAINING GLOBAL PARAMETERS.
*----
      DO 10 IPAR=1,NPAR
      IF(PARTYP(IPAR).EQ.'VALE') THEN
         GO TO 10
      ELSE IF((PARTYP(IPAR).EQ.'IRRA').OR.(PARTYP(IPAR).EQ.'TIME').OR.
     1        (PARTYP(IPAR).EQ.'PUIS').OR.(PARTYP(IPAR).EQ.'FLUB').OR.
     2        (PARTYP(IPAR).EQ.'FLUX').OR.(PARTYP(IPAR).EQ.'MASL')) THEN
*
*        RECOVER GLOBAL PARAMETER VALUES FROM THE DEPLETION OBJECT.
         IF(.NOT.C_ASSOCIATED(IPDEPL)) CALL XABORT('SAPGEP: NO DEPLETI'
     1   //'ON OBJECT AVAILABLE AMONG THE RHS LCM OBJECTS.')
         CALL LCMGET(IPDEPL,'STATE-VECTOR',IDATA)
         NBURN=IDATA(3)
         NBISO=IDATA(4)
         NREAC=IDATA(6)
         NVAR=IDATA(7)
         NBMIX=IDATA(8)
         CALL COMGEM(IPDEPL,ITIM,PARTYP(IPAR),0,NBURN,NBMIX,NBISO,
     1   NREAC,NVAR,VALPAR)
      ELSE IF((PARTYP(IPAR).EQ.'TEMP').OR.(PARTYP(IPAR).EQ.'CONC'))
     1   THEN
*
*        RECOVER GLOBAL PARAMETER VALUES FROM A MICROLIB OBJECT.
         IF(.NOT.C_ASSOCIATED(IPLB1)) CALL XABORT('SAPGEP: MICROLIB EX'
     1   //'PECTED AT RHS.')
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
         CALL LCMGET(IPLB1,'STATE-VECTOR',IDATA)
         MAXNBI=IDATA(2)
         IF(C_ASSOCIATED(IPLB2)) THEN
            CALL LCMGET(IPLB2,'STATE-VECTOR',IDATA)
            MAXNBI=MAX(MAXNBI,IDATA(2))
         ENDIF
         CALL COMBIB(IPLB1,IPLB2,PARTYP(IPAR),IMILI,TEXT12,TEXT8,MAXNBI,
     1   VALPAR)
         IF(PARTYP(IPAR).EQ.'TEMP') VALPAR=VALPAR-273.16
      ELSE
         CALL XABORT('SAPGEP: '//PARTYP(IPAR)//' IS AN UNKNOWN PARAM'//
     1   'ETER TYPE.')
      ENDIF
      IF(IMPX.GT.0) WRITE(6,100) PARKEY(IPAR),VALPAR
*
      CALL SAPPAV(IPSAP,IPAR,NPAR,'FLOTTANT',VALPAR,NITMA,TEXT12,
     1 MUPLET(IPAR),LGNEW(IPAR))
   10 CONTINUE
      IF(IMPX.GT.2) THEN
         WRITE(6,110) (MUPLET(I),I=1,NPAR)
         WRITE(6,'(/)')
      ENDIF
      DO 15 I=1,NPAR
      IF(MUPLET(I).EQ.0) THEN
         WRITE(HSMG,'(33HSAPGEP: UNDEFINED MUPLET ELEMENT=,I6)') I
         CALL XABORT(HSMG)
      ENDIF
   15 CONTINUE
*----
*  INTRODUCE VALUES INTO GLOBAL PARAMETER TREE.
*----
**
** Parameter tree: this tree has a number of stages equal to the
** number of parameters. For each value of the i-th parameter, we
** find the position in the tree corresponding to the value of the
** (i+1)-th parameter.
** NCALAR  Number of elementary calculations stored in the tree.
** NVP     Number of nodes in the parameter tree, including the root.
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
      CALL LCMSIX(IPSAP,'paramarbre',1)
      CALL LCMLEN(IPSAP,'ARBVAL',MAXNVP,ITYLCM)
      IF(MAXNVP.EQ.0) THEN
         MAXNVP=100*(NPAR+1)
         ALLOCATE(IDEBAR(MAXNVP+1),IARBVA(MAXNVP))
         CALL XDISET(IDEBAR,MAXNVP+1,0)
         CALL XDISET(IARBVA,MAXNVP,0)
         IARBVA=0
         DO 20 I=1,NPAR
         IDEBAR(I)=I+1
         IARBVA(I+1)=1
   20    CONTINUE
         IDEBAR(NPAR+1)=NPAR+2
         IDEBAR(NPAR+2)=1
         NCALAR=1
         NVPNEW=NPAR+1
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
         ALLOCATE(JDEBAR(MAXNVP+1),JARBVA(MAXNVP))
         CALL LCMGET(IPSAP,'DEBARB',JDEBAR)
         CALL LCMGET(IPSAP,'ARBVAL',JARBVA)
         DO 30 IPAR=1,NPAR
         IF(LGNEW(IPAR)) THEN
            II=IPAR
            GO TO 40
         ENDIF
   30    CONTINUE
         II=NPAR+1
   40    LGERR=COMTRE(NPAR,NVPO,JARBVA,JDEBAR,MUPLET,KK,
     1   I0,IORD,JJ,LAST)
         IF((II.GT.NPAR).AND.LGERR) THEN
            WRITE(TEXT4,'(I4)') IORD
            CALL XABORT('SAPGEP: ELEMENTARY CALCULATION HAS THE SAME'//
     1      ' GLOBAL PARAMETERS AS ELEMENTARY CALCULATION NB '//TEXT4)
         ENDIF
*
*        Size of the new tree.
*
         NVPNEW=NVPO+NPAR+1-MIN(II,KK)
         IF(NVPNEW.GT.MAXNVP) MAXNVP=NVPNEW+MAXNVP
         ALLOCATE(IDEBAR(MAXNVP+1),IARBVA(MAXNVP))
         CALL XDISET(IDEBAR(NVPNEW+2),MAXNVP-NVPNEW,0)
         CALL XDISET(IARBVA(NVPNEW+1),MAXNVP-NVPNEW,0)
*
*        Update values and suppress old PARBRE.
*
         CALL COMARB(NPAR,NVPO,NVPNEW,JDEBAR,JARBVA,LGNEW,MUPLET,NCALAR,
     1   IDEBAR,IARBVA)
         DEALLOCATE(JARBVA,JDEBAR)
      ENDIF
      CALL LCMPUT(IPSAP,'NCALS',1,1,NCALAR)
      CALL LCMPUT(IPSAP,'DEBARB',NVPNEW+1,1,IDEBAR)
      CALL LCMPUT(IPSAP,'ARBVAL',NVPNEW,1,IARBVA)
      DEALLOCATE(IARBVA,IDEBAR)
      IF(NCALAR.EQ.1) THEN
         MAXNCA=1000
         ALLOCATE(IORIGI(MAXNCA))
         CALL XDISET(IORIGI,MAXNCA,0)
      ELSE
         CALL LCMLEN(IPSAP,'ORIGIN',MAXNCA,ITYLCM)
         IF(NCALAR.GT.MAXNCA) MAXNCA=NCALAR+MAXNCA
         ALLOCATE(IORIGI(MAXNCA))
         CALL XDISET(IORIGI,MAXNCA,0)
         CALL LCMGET(IPSAP,'ORIGIN',IORIGI)
      ENDIF
      IORIGI(NCALAR)=NORIG
      CALL LCMPUT(IPSAP,'ORIGIN',MAXNCA,1,IORIGI)
      DEALLOCATE(IORIGI)
      CALL LCMSIX(IPSAP,' ',2)
*----
*  RECOVER INFORMATION FROM THE 'varlocdescri' DIRECTORY.
*----
   45 IF(NLOC.EQ.0) RETURN
      CALL LCMSIX(IPSAP,'varlocdescri',1)
      CALL LCMGTC(IPSAP,'PARNAM',80,NPAR,PARNAM)
      CALL LCMGTC(IPSAP,'PARKEY',4,NLOC,PARKEY)
      CALL LCMGTC(IPSAP,'PARTYP',4,NLOC,PARTYP)
      CALL LCMGTC(IPSAP,'PARFMT',8,NLOC,PARFMT)
      CALL LCMGET(IPSAP,'PARCAD',PARCAD)
      IF(NPCHRL.GT.0) CALL LCMGTC(IPSAP,'PARCHR',8,NPCHRL,PARCHR)
      CALL LCMSIX(IPSAP,' ',2)
*
      CALL LCMGTC(IPEDIT,'LAST-EDIT',12,1,TEXT12)
*----
*  INITIALIZE LOCADR AND ALLOCATE RVALOC.
*----
      IADR=1
      LOCADR(1)=1
      DO 50 IPAR=1,NLOC
      IF((PARTYP(IPAR).EQ.'EQUI').OR.(PARTYP(IPAR).EQ.'VITE')) THEN
         IADR=IADR+NG
      ELSE IF(PARTYP(IPAR).EQ.'COUR') THEN
         IADR=IADR+2*NG
      ELSE
         IADR=IADR+1
      ENDIF
      LOCADR(IPAR+1)=IADR
   50 CONTINUE
      NVLC=LOCADR(NLOC+1)-1
      ALLOCATE(RVALO(NVLC*NMIL))
*----
*  RECOVER LOCAL VARIABLES.
*----
      DO 70 IPAR=1,NLOC
      IF((PARTYP(IPAR).EQ.'IRRA').OR.(PARTYP(IPAR).EQ.'TIME').OR.
     1   (PARTYP(IPAR).EQ.'PUIS').OR.(PARTYP(IPAR).EQ.'FLUG').OR.
     2   (PARTYP(IPAR).EQ.'FLUB').OR.(PARTYP(IPAR).EQ.'FLUX').OR.
     3   (PARTYP(IPAR).EQ.'MASL')) THEN
*
*        RECOVER LOCAL VARIABLES FROM THE DEPLETION OBJECT.
         IF(.NOT.C_ASSOCIATED(IPDEPL)) CALL XABORT('SAPGEP: NO DEPLET'
     1   //'ION OBJECT AVAILABLE AMONG THE RHS LCM OBJECTS.')
         CALL LCMGET(IPDEPL,'STATE-VECTOR',IDATA)
         NBURN=IDATA(3)
         NBISO=IDATA(4)
         NREAC=IDATA(6)
         NVAR=IDATA(7)
         NBMIX=IDATA(8)
         CALL LCMGET(IPEDIT,'STATE-VECTOR',IDATA)
         NREG=IDATA(17)
         CALL COMGEN(IPDEPL,IPEDIT,NREG,NMIL,ITIM,PARTYP(IPAR),NBURN,
     1   NBMIX,NBISO,NREAC,NVAR,LOCADR(IPAR),NVLC,RVALO)
      ELSE IF((PARTYP(IPAR).EQ.'TEMP').OR.(PARTYP(IPAR).EQ.'CONC'))
     1   THEN
*
*        RECOVER LOCAL VARIABLES FROM THE MICROLIB IN EDIT OBJECT.
         IPCAD=PARCAD(IPAR+1)-PARCAD(IPAR)
         IF(IPCAD.EQ.1) IPCAD=PARCAD(IPAR+1)-PARCAD(1)
         TEXT8=' '
         IF(IPCAD.GT.0) TEXT8=PARCHR(IPCAD)
         CALL LCMSIX(IPEDIT,TEXT12,1)
         CALL LCMGET(IPEDIT,'STATE-VECTOR',IDATA)
         MAXNBI=IDATA(2)
         CALL LCMINF(IPEDIT,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
         IPLB3=C_NULL_PTR
         DO 60 IBM=1,NMIL
         CALL COMBIB(IPEDIT,IPLB3,PARTYP(IPAR),IBM,NAMLCM,TEXT8,MAXNBI,
     1   VALPAR)
         IF(PARTYP(IPAR).EQ.'TEMP') VALPAR=VALPAR-273.16
         RVALO((IBM-1)*NVLC+LOCADR(IPAR))=VALPAR
   60    CONTINUE
         CALL LCMSIX(IPEDIT,' ',2)
      ELSE IF(PARTYP(IPAR).EQ.'EQUI') THEN
*        RECOVER A SET OF SPH EQUIVALENCE FACTORS.
         CALL SAPSPH(IPEDIT,NG,NMIL,LOCADR(IPAR),NVLC,RVALO)
      ELSE
         CALL XABORT('SAPGEP: '//PARTYP(IPAR)//' IS AN UNKNOWN LOCAL'//
     1   ' VARIABLE TYPE.')
      ENDIF
      IF(IMPX.GT.1) WRITE(6,120) PARKEY(IPAR),
     1            (RVALO((IBM-1)*NVLC+LOCADR(IPAR)),IBM=1,NMIL)
   70 CONTINUE
      WRITE(TEXT12,'(''calc'',I8)') NCALAR
      CALL LCMSIX(IPSAP,TEXT12,1)
      CALL LCMSIX(IPSAP,'info',1)
      CALL LCMPUT(IPSAP,'NLOC',1,1,NLOC)
      CALL LCMPTC(IPSAP,'LOCNAM',80,NLOC,PARNAM)
      CALL LCMPTC(IPSAP,'LOCKEY',4,NLOC,PARKEY)
      CALL LCMPTC(IPSAP,'LOCTYP',4,NLOC,PARTYP)
      CALL LCMPUT(IPSAP,'LOCADR',NLOC+1,1,LOCADR)
      CALL LCMSIX(IPSAP,' ',2)
      DO 80 IBM=1,NMIL
      WRITE(TEXT12,'(''mili'',I8)') IBM
      CALL LCMSIX(IPSAP,TEXT12,1)
      CALL LCMPUT(IPSAP,'RVALOC',NVLC,2,RVALO((IBM-1)*NVLC+1))
      CALL LCMSIX(IPSAP,' ',2)
   80 CONTINUE
      CALL LCMSIX(IPSAP,' ',2)
      DEALLOCATE(RVALO)
      RETURN
*
  100 FORMAT(31H SAPGEP: SET GLOBAL PARAMETER ',A,3H' =,1P,E12.4)
  110 FORMAT(/16H SAPGEP: MUPLET=,10I6:/(16X,10I6))
  120 FORMAT(29H SAPGEP: SET LOCAL VARIABLE ',A,3H' =,1P,5E12.4/(36X,
     1 5E12.4))
      END
