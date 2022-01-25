*DECK NCRTRP
      SUBROUTINE NCRTRP(IPCPO,LCUB2,IMPX,IBMOLD,NPAR,NLOC,NCAL,MUPLET,
     1 MUTYPE,VALR,VARVAL,TERP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the TERP interpolation/derivation/integration factors using
* table-of-content information of the multicompo for mixture IBMOLD.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert, R. Chambon
*
*Parameters: input
* IPCPO   address of the multidimensional multicompo object.
* LCUB2   interpolation type for each parameter (=.TRUE.: cubic Ceschino
*         interpolation; =.FALSE: linear Lagrange interpolation).
* IMPX    print parameter (equal to zero for no print).
* IBMOLD  material mixture index in the multicompo.
* NPAR    number of global parameters.
* NLOC    number of local parameters.
* NCAL    number of elementary calculations in the multicompo.
* MUPLET  tuple used to identify an elementary calculation.
* MUTYPE  type of interpolation (=1: interpolation; =2: delta-sigma).
* VALR    real values of the interpolated point.
* VARVAL  exit burnup used if MUTYPE(IPAR(ID))=3.
*
*Parameters: output
* TERP    interpolation factors.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER, PARAMETER::MAXPAR=50
      TYPE(C_PTR) IPCPO
      INTEGER IMPX,IBMOLD,NPAR,NLOC,NCAL,MUPLET(NPAR+NLOC),
     1 MUTYPE(NPAR+NLOC)
      REAL VALR(2*MAXPAR,2),VARVAL,TERP(NCAL)
      LOGICAL LCUB2(NPAR+NLOC)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXDIM=10
      INTEGER, PARAMETER::MAXVAL=200
      INTEGER IPAR(MAXDIM),NVALUE(MAXPAR),NVAL(MAXDIM),IDDIV(MAXDIM),
     1 NVPO(2)
      REAL VREAL(MAXVAL),T1D(MAXVAL,MAXDIM),WORK(MAXVAL)
      REAL BURN0, BURN1, DENOM, TERTMP
      INTEGER ICAL, IDTMP, IDTOT, ID, ILONG, ITYLCM, I, JD, MAXNVP,
     & NDELTA, NDIM, NID, NTOT, NCRCAL
      CHARACTER HSMG*131,RECNAM*12,PARKEY(MAXPAR)*12
      LOGICAL LCUBIC,LSINGL
      TYPE(C_PTR) JPCPO,KPCPO,LPCPO
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JDEBAR,JARBVA
      REAL, ALLOCATABLE, DIMENSION(:) :: TERPA
*----
*  RECOVER TREE INFORMATION
*----
      JPCPO=LCMGID(IPCPO,'GLOBAL')
      CALL LCMGTC(JPCPO,'PARKEY',12,NPAR,PARKEY)
      JPCPO=LCMGID(IPCPO,'MIXTURES')
      KPCPO=LCMGIL(JPCPO,IBMOLD)
      LPCPO=LCMGID(KPCPO,'TREE')
      CALL LCMGET(LPCPO,'NVP',NVPO)
      CALL LCMLEN(LPCPO,'ARBVAL',MAXNVP,ITYLCM)
      IF(NVPO(1).GT.MAXNVP) CALL XABORT('NCRTRP: NVP OVERFLOW.')
      ALLOCATE(JDEBAR(MAXNVP+1),JARBVA(MAXNVP))
      CALL LCMGET(LPCPO,'DEBARB',JDEBAR)
      CALL LCMGET(LPCPO,'ARBVAL',JARBVA)
*----
*  COMPUTE TERP FACTORS
*----
      TERP(:NCAL)=0.0
      IPAR(:MAXDIM)=0
      NDIM=0
      NDELTA=0
      DO 10 I=1,NPAR+NLOC
        IF(MUPLET(I).EQ.-1) THEN
          NDIM=NDIM+1
          IF(MUTYPE(I).NE.1) NDELTA=NDELTA+1
          IF(NDIM.GT.MAXDIM) THEN
            WRITE(HSMG,'(7HNCRTRP:,I4,29H-DIMENSIONAL INTERPOLATION NO,
     1      14HT IMPLEMENTED.)') NDIM
            CALL XABORT(HSMG)
          ENDIF
          IPAR(NDIM)=I
        ENDIF
   10 CONTINUE
      IF(IMPX.GT.2) THEN
        WRITE(IOUT,'(16H NCRTRP: MUPLET=,10I4/(16X,10I4))')
     1  (MUPLET(I),I=1,NPAR)
        WRITE(IOUT,'(8H NCRTRP:,I4,31H-DIMENSIONAL INTERPOLATION IN C,
     1  12HOMPO MIXTURE,I5,1H.)') NDIM,IBMOLD
      ENDIF
      IF(NDIM.EQ.0) THEN
        ICAL=NCRCAL(1,NVPO(1),NPAR+NLOC,JDEBAR,JARBVA,MUPLET)
        IF(ICAL.GT.NCAL) CALL XABORT('NCRTRP: TERP OVERFLOW(1).')
        IF(ICAL.EQ.0) GO TO 200
        IF(ICAL.EQ.-1) GO TO 210
        TERP(ICAL)=1.0
      ELSE
        NTOT=1
        IDDIV(:MAXDIM)=1
        DO 70 ID=1,NDIM
        IF(IPAR(ID).LE.NPAR) THEN
          LPCPO=LCMGID(IPCPO,'GLOBAL')
          WRITE(RECNAM,'(''pval'',I8.8)') IPAR(ID)
          CALL LCMGET(LPCPO,'NVALUE',NVALUE)
          NID=NVALUE(IPAR(ID))
        ELSE
          JPCPO=LCMGID(IPCPO,'MIXTURES')
          KPCPO=LCMGIL(JPCPO,IBMOLD)
          LPCPO=LCMGID(KPCPO,'TREE')
          WRITE(RECNAM,'(''pval'',I8.8)') IPAR(ID)-NPAR
          CALL LCMGET(LPCPO,'NVALUE',NVALUE)
          NID=NVALUE(IPAR(ID)-NPAR)
        ENDIF
        NTOT=NTOT*NID
        DO 15 IDTMP=1,NDIM-ID
        IDDIV(IDTMP)=IDDIV(IDTMP)*NID
   15   CONTINUE
        CALL LCMLEN(LPCPO,RECNAM,ILONG,ITYLCM)
        IF(ILONG.GT.MAXVAL) CALL XABORT('NCRTRP: MAXVAL OVERFLOW.')
        CALL LCMGET(LPCPO,RECNAM,VREAL)
        BURN0=VALR(IPAR(ID),1)
        BURN1=VALR(IPAR(ID),2)
        LSINGL=(BURN0.EQ.BURN1)
        LCUBIC=LCUB2(IPAR(ID))
        IF((MUTYPE(IPAR(ID)).EQ.1).AND.LSINGL) THEN
          CALL ALTERP(LCUBIC,NID,VREAL,BURN0,.FALSE.,T1D(1,ID))
        ELSE IF(MUTYPE(IPAR(ID)).EQ.1) THEN
          IF(BURN0.GE.BURN1) CALL XABORT('@NCRTRP: INVALID BURNUP'
     1     //' LIMITS(1).')
          CALL ALTERI(LCUBIC,NID,VREAL,BURN0,BURN1,T1D(1,ID))
          DO 20 I=1,NID
          T1D(I,ID)=T1D(I,ID)/(BURN1-BURN0)
   20     CONTINUE
        ELSE IF((MUTYPE(IPAR(ID)).EQ.2).AND.(.NOT.LSINGL)) THEN
          CALL ALTERP(LCUBIC,NID,VREAL,BURN0,.FALSE.,WORK)
          CALL ALTERP(LCUBIC,NID,VREAL,BURN1,.FALSE.,T1D(1,ID))
          DO 30 I=1,NID
          T1D(I,ID)=T1D(I,ID)-WORK(I)
   30     CONTINUE
        ELSE IF((MUTYPE(IPAR(ID)).EQ.2).AND.(LSINGL)) THEN
          T1D(:NID,ID)=0.0
        ELSE IF(MUTYPE(IPAR(ID)).EQ.3) THEN
*          DERIVATIVE WITH RESPECT TO A SINGLE EXIT BURNUP. USE
*          EQ.(3.3) OF RICHARD CHAMBON'S THESIS.
          IF(BURN0.GE.BURN1) CALL XABORT('@NCRTRP: INVALID BURNUP'
     1    //' LIMITS(2).')
          IF(PARKEY(IPAR(ID)).NE.'BURN') THEN
            CALL XABORT('@NCRTRP: BURN EXPECTED.')
          ENDIF
          ALLOCATE(TERPA(NID))
          CALL ALTERI(LCUBIC,NID,VREAL,BURN0,BURN1,TERPA)
          DO 40 I=1,NID
          T1D(I,ID)=-TERPA(I)
   40     CONTINUE
          CALL ALTERP(LCUBIC,NID,VREAL,BURN0,.FALSE.,TERPA)
          DO 50 I=1,NID
          T1D(I,ID)=T1D(I,ID)-TERPA(I)*BURN0
   50     CONTINUE
          CALL ALTERP(LCUBIC,NID,VREAL,BURN1,.FALSE.,TERPA)
          DENOM=VARVAL*(BURN1-BURN0)
          DO 60 I=1,NID
          T1D(I,ID)=(T1D(I,ID)+TERPA(I)*BURN1)/DENOM
   60     CONTINUE
          DEALLOCATE(TERPA)
        ELSE
          CALL XABORT('NCRTRP: INVALID OPTION.')
        ENDIF
        NVAL(ID)=NID
   70   CONTINUE

* Example: NDIM=3, NVALUE=(3,2,2)
* IDTOT 1  2  3  4  5  6  7  8  9 10 11 12
* ID(1) 1  2  3  1  2  3  1  2  3  1  2  3
* ID(2) 1  1  1  2  2  2  1  1  1  2  2  2
* ID(3) 1  1  1  1  1  1  2  2  2  2  2  2
* (NTOT=12, IDDIV=(6,3,1))
        DO 100 IDTOT=1,NTOT               ! Ex.: IDTOT       = 9
          TERTMP=1.0
          IDTMP=IDTOT
          DO 80 JD=1,NDIM                 ! Ex.: JD          = 1,2,3
            ID=(IDTMP-1)/IDDIV(JD)+1      ! Ex.: ID(NDIM...1)= 2,1,3
            IDTMP=IDTMP-(ID-1)*IDDIV(JD)  ! Ex.: IDTMP       = 3,3,1
            MUPLET(IPAR(NDIM-JD+1))=ID
            TERTMP=TERTMP*T1D(ID,NDIM-JD+1)
   80     CONTINUE
          ICAL=NCRCAL(1,NVPO(1),NPAR+NLOC,JDEBAR,JARBVA,MUPLET)
          IF(ICAL.GT.NCAL) CALL XABORT('NCRTRP: TERP OVERFLOW(2).')
          IF(ICAL.EQ.0) GO TO 200
          IF(ICAL.EQ.-1) GO TO 210
          TERP(ICAL)=TERP(ICAL)+TERTMP
  100   CONTINUE
      ENDIF
      IF(IMPX.GT.3) THEN
        WRITE(IOUT,'(35H NCRTRP: TERP PARAMETERS IN MIXTURE,I4,1H:/(1X,
     1  1P,10E12.4))') IBMOLD,(TERP(I),I=1,NCAL)
      ENDIF
      DEALLOCATE(JARBVA,JDEBAR)
      RETURN
*----
*  MISSING ELEMENTARY CALCULATION EXCEPTION.
*----
  200 WRITE(IOUT,'(16H NCRTRP: MUPLET=,10I4/(16X,10I4))')
     1 (MUPLET(I),I=1,NPAR+NLOC)
      CALL XABORT('NCRTRP: MISSING ELEMENTARY CALCULATION.')
  210 WRITE(IOUT,'(16H NCRTRP: MUPLET=,10I4/(16X,10I4))')
     1 (MUPLET(I),I=1,NPAR+NLOC)
      CALL XABORT('NCRTRP: DEGENERATE ELEMENTARY CALCULATION.')
      END
