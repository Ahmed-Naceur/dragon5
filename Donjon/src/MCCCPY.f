*DECK MCCCPY
      SUBROUTINE MCCCPY(IMPX,IPMAP,IPMAP2,NPARAM,NCH,NB,REC1,REC2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Modification of the data stored in the PARAM folder of a fuel map
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal.
*
*Author(s): 
* M. Cordiez
*
*Parameters: input
* IMPX    printing index (=0 for no print).
* IPMAP   pointer of the fuel map 1
* IPMAP2  pointer of the fuel map 2
* NPARAM  number of parameters in the PARAM folder
* NCH     number of fuel channels
* NB      number of fuel bundles per channel
* REC1    record to be updated
* REC2    record to be copied
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,NPARAM,NCH,NB
      CHARACTER REC1*40,REC2*40
      TYPE(C_PTR) IPMAP,IPMAP2
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE2(NSTATE)
      INTEGER VALSIZE,PTYPE
      INTEGER PTYPE2,NPARAM2,NCH2,NB2,VALSIZE2
      CHARACTER PNAME*12,PNAME2*12
      REAL, ALLOCATABLE, DIMENSION(:) :: VALMOD,VALCOP
      TYPE(C_PTR) JPMAP,KPMAP,JPMAP2,KPMAP2
*----
*  RECOVER L_MAP 2 STATE-VECTOR
*----
      CALL LCMGET(IPMAP2,'STATE-VECTOR',ISTATE2)
      NB2=ISTATE2(1)
      NCH2=ISTATE2(2)
      NPARAM2=ISTATE2(8)
      IF(NB.NE.NB2) CALL XABORT('@MCCCCPY: THE NUMBER OF FUEL'
     >            //' BUNDLES PER CHANNEL IS DIFFERENT BETWEEN'
     >            //' FLMAP1 AND FLMAP2.')
      IF(NCH.NE.NCH2) CALL XABORT('@MCCCCPY: THE NUMBER OF FUEL'
     >            //' CHANNELS IS DIFFERENT BETWEEN'
     >            //' FLMAP1 AND FLMAP2.')
*----
*  RECOVERY OF L_MAP PARAMETERS
*----
* L_MAP1 (to be updated)
      JPMAP=LCMGID(IPMAP,'PARAM')
      DO IPAR=1,NPARAM,1
        KPMAP=LCMGIL(JPMAP,IPAR)
        CALL LCMGTC(KPMAP,'P-NAME',12,1,PNAME)
        CALL LCMGET(KPMAP,'P-TYPE',PTYPE)
        IF(PNAME.EQ.REC1) THEN
          IF(IMPX.GE.3) CALL LCMLIB(KPMAP)
          EXIT
        ENDIF
      ENDDO
* L_MAP2 (to be copied)
      JPMAP2=LCMGID(IPMAP2,'PARAM')
      DO IPAR=1,NPARAM2,1
        KPMAP2=LCMGIL(JPMAP2,IPAR)
        CALL LCMGTC(KPMAP2,'P-NAME',12,1,PNAME2)
        CALL LCMGET(KPMAP2,'P-TYPE',PTYPE2)
        IF(PNAME2.EQ.REC2) THEN
          IF(IMPX.GE.3) CALL LCMLIB(KPMAP2)
          EXIT
        ENDIF
      ENDDO
*
* Checking of the type (local or global) of REC1
      IF(PTYPE.EQ.1) THEN
        IF(IMPX.GE.1) WRITE(6,210) 'PARAMETER ',PNAME,' IS GLOBAL'
        VALSIZE=1
        ALLOCATE(VALMOD(VALSIZE))
      ELSE IF(PTYPE.EQ.2) THEN
        IF(IMPX.GE.1) WRITE(6,210) 'PARAMETER ',PNAME,' IS LOCAL'
        VALSIZE=NCH*NB
        ALLOCATE(VALMOD(VALSIZE))
      ELSE 
        CALL XABORT('@MCCCPY: '//PNAME//'IS NEITHER LOCAL NOR GLOBAL'
     >         //'AND THAT IS IMPOSSIBLE.') 
      ENDIF
* Checking of the type (local or global) of REC2
      IF((PTYPE2.NE.1).AND.(PTYPE2.NE.2)) THEN
        CALL XABORT('@MCCCPY: '//PNAME2//'IS NEITHER LOCAL NOR GLOBAL'
     >         //'AND THAT IS IMPOSSIBLE.') 
      ENDIF
      IF(PTYPE2.EQ.1) THEN
        IF(IMPX.GE.1) WRITE(6,210) 'PARAMETER ',PNAME2,' IS GLOBAL'
        VALSIZE2=1
        ALLOCATE(VALCOP(VALSIZE2))
        CALL LCMGET(KPMAP2,'P-VALUE',VALCOP)
      ELSE
        IF(IMPX.GE.1) WRITE(6,210) 'PARAMETER ',PNAME2,' IS LOCAL'
        VALSIZE2=NCH2*NB2
        ALLOCATE(VALCOP(VALSIZE2))
        CALL LCMGET(KPMAP2,'P-VALUE',VALCOP)
      ENDIF
      IF(PTYPE.EQ.1.AND.PTYPE2.EQ.2) CALL XABORT('@MCCCPY: '//PNAME
     >         //'IS GLOBAL ON THE CORE AND '//PNAME2//'IS LOCAL.')
      IF(PTYPE.EQ.2.AND.PTYPE2.EQ.1) CALL XABORT('@MCCCPY: '//PNAME
     >         //'IS LOCAL ON THE CORE AND '//PNAME2//'IS GLOBAL.')
*
* Modification of REC1

      VALMOD(:)=VALCOP(:)

      CALL LCMPUT(KPMAP,'P-VALUE',VALSIZE,2,VALMOD)

      IF(IMPX.GT.0.AND.C_ASSOCIATED(IPMAP,IPMAP2)) THEN
        WRITE(6,220) 'EVERY VALUE OF THE RECORD ',REC1,' HAS BEEN '
     >             //'UPDATED WITH THE ONES FROM ',REC2,' (SAME '
     >             //'FUEL MAP).'
      ELSEIF(IMPX.GT.0.AND.C_ASSOCIATED(IPMAP,IPMAP2)) THEN
        WRITE(6,220) 'EVERY VALUE OF THE RECORD ',REC1,' HAS BEEN '
     >             //'UPDATED WITH THE ONES FROM ',REC2,' (FUEL '
     >             //'MAP FLMAP2).'
      ENDIF

* Array deallocation
      DEALLOCATE(VALMOD)
      DEALLOCATE(VALCOP)

      RETURN

  210 FORMAT(1X,A,A6,A/)
  220 FORMAT(1X,A,A6,A,A6,A/)
      END
