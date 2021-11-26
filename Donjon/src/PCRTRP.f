*DECK PCRTRP
      SUBROUTINE PCRTRP(LCUB2,IMPX,NPAR,NCAL,NVALUE,MUPLET,MUTYPE,
     1 VALR,VARVAL,MUBASE,VREAL,TERP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the TERP interpolation/derivation/integration factors using
* table-of-content information of the PMAXS file.
*
*Copyright:
* Copyright (C) 2018 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* LCUB2   interpolation type for each parameter (=.TRUE.: cubic Ceschino
*         interpolation; =.FALSE: linear Lagrange interpolation).
* IMPX    print parameter (equal to zero for no print).
* NPAR    number of parameters.
* NCAL    number of elementary calculations in the PMAXS file.
* NVALUE  number of tabulation values for each parameter.
* MUPLET  tuple used to identify an elementary calculation.
* MUTYPE  type of interpolation (=1: interpolation; =2: delta-sigma).
* VALR    real values of the interpolated point.
* VARVAL  exit burnup used if MUTYPE(IPAR(ID))=3.
* MUBASE  muplet database.
* VREAL   local parameter values at tabulation points.
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
      INTEGER, PARAMETER::MAXVAL=200
      INTEGER, PARAMETER::MAXPAR=50
      INTEGER IMPX,NPAR,NCAL,NVALUE(NPAR),MUPLET(NPAR),MUTYPE(NPAR),
     1 MUBASE(NPAR,NCAL)
      REAL VALR(MAXPAR,2),VARVAL,VREAL(MAXVAL,MAXPAR),TERP(NCAL)
      LOGICAL LCUB2(NPAR)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXDIM=10
      INTEGER IPAR(MAXDIM),NVAL(MAXDIM),IDDIV(MAXDIM)
      REAL T1D(MAXVAL,MAXDIM),WORK(MAXVAL)
      REAL BURN0, BURN1, DENOM, TERTMP
      INTEGER ICAL, IDTMP, IDTOT, ID, I, JD, NDELTA, NDIM, NID, NTOT,
     1 IIPAR, PCRCAL
      CHARACTER HSMG*131,RECNAM*12
      LOGICAL LCUBIC,LSINGL
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: TERPA
*----
*  COMPUTE TERP FACTORS
*----
      TERP(:NCAL)=0.0
      IPAR(:MAXDIM)=0
      NDIM=0
      NDELTA=0
      DO 10 I=1,NPAR
        IF(MUPLET(I).EQ.-1) THEN
          NDIM=NDIM+1
          IF(MUTYPE(I).NE.1) NDELTA=NDELTA+1
          IF(NDIM.GT.MAXDIM) THEN
            WRITE(HSMG,'(7HPCRTRP:,I4,29H-DIMENSIONAL INTERPOLATION NO,
     1      14HT IMPLEMENTED.)') NDIM
            CALL XABORT(HSMG)
          ENDIF
          IPAR(NDIM)=I
        ENDIF
   10 CONTINUE
      IF(IMPX.GT.2) THEN
        WRITE(IOUT,'(16H PCRTRP: MUPLET=,10I4/(16X,10I4))')
     1  (MUPLET(I),I=1,NPAR)
        WRITE(IOUT,'(8H PCRTRP:,I4,27H-DIMENSIONAL INTERPOLATION.)')
     1  NDIM
      ENDIF
      IF(NDIM.EQ.0) THEN
        ICAL=PCRCAL(NPAR,NCAL,MUPLET,MUBASE)
        IF(ICAL.GT.NCAL) CALL XABORT('PCRTRP: TERP OVERFLOW(1).')
        IF(ICAL.EQ.0) GO TO 200
        IF(ICAL.EQ.-1) GO TO 210
        TERP(ICAL)=1.0
      ELSE
        NTOT=1
        IDDIV(:MAXDIM)=1
        DO 70 ID=1,NDIM
        IIPAR=IPAR(ID)
        NID=NVALUE(IIPAR)
        NTOT=NTOT*NID
        DO 15 IDTMP=1,NDIM-ID
        IDDIV(IDTMP)=IDDIV(IDTMP)*NID
   15   CONTINUE
        BURN0=VALR(IIPAR,1)
        BURN1=VALR(IIPAR,2)
        LSINGL=(BURN0.EQ.BURN1)
        LCUBIC=LCUB2(IIPAR)
        IF((MUTYPE(IIPAR).EQ.1).AND.LSINGL) THEN
          CALL ALTERP(LCUBIC,NID,VREAL(1,IIPAR),BURN0,.FALSE.,
     1    T1D(1,ID))
        ELSE IF(MUTYPE(IIPAR).EQ.1) THEN
          IF(BURN0.GE.BURN1) CALL XABORT('@PCRTRP: INVALID BURNUP'
     1     //' LIMITS(1).')
          CALL ALTERI(LCUBIC,NID,VREAL(1,IIPAR),BURN0,BURN1,T1D(1,ID))
          DO 20 I=1,NID
          T1D(I,ID)=T1D(I,ID)/(BURN1-BURN0)
   20     CONTINUE
        ELSE IF((MUTYPE(IIPAR).EQ.2).AND.(.NOT.LSINGL)) THEN
          CALL ALTERP(LCUBIC,NID,VREAL(1,IIPAR),BURN0,.FALSE.,WORK)
          CALL ALTERP(LCUBIC,NID,VREAL(1,IIPAR),BURN1,.FALSE.,T1D(1,ID))
          DO 30 I=1,NID
          T1D(I,ID)=T1D(I,ID)-WORK(I)
   30     CONTINUE
        ELSE IF((MUTYPE(IIPAR).EQ.2).AND.(LSINGL)) THEN
          T1D(:NID,ID)=0.0
        ELSE IF(MUTYPE(IIPAR).EQ.3) THEN
*          DERIVATIVE WITH RESPECT TO A SINGLE EXIT BURNUP. USE
*          EQ.(3.3) OF RICHARD CHAMBON'S THESIS.
          IF(BURN0.GE.BURN1) CALL XABORT('@PCRTRP: INVALID BURNUP'
     1    //' LIMITS(2).')
          IF(RECNAM.NE.'BURN') CALL XABORT('@PCRTRP: BURN EXPECTED.')
          ALLOCATE(TERPA(NID))
          CALL ALTERI(LCUBIC,NID,VREAL(1,IIPAR),BURN0,BURN1,TERPA)
          DO 40 I=1,NID
          T1D(I,ID)=-TERPA(I)
   40     CONTINUE
          CALL ALTERP(LCUBIC,NID,VREAL(1,IIPAR),BURN0,.FALSE.,TERPA)
          DO 50 I=1,NID
          T1D(I,ID)=T1D(I,ID)-TERPA(I)*BURN0
   50     CONTINUE
          CALL ALTERP(LCUBIC,NID,VREAL(1,IIPAR),BURN1,.FALSE.,TERPA)
          DENOM=VARVAL*(BURN1-BURN0)
          DO 60 I=1,NID
          T1D(I,ID)=(T1D(I,ID)+TERPA(I)*BURN1)/DENOM
   60     CONTINUE
          DEALLOCATE(TERPA)
        ELSE
          CALL XABORT('PCRTRP: INVALID OPTION.')
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
          ICAL=PCRCAL(NPAR,NCAL,MUPLET,MUBASE)
          IF(ICAL.GT.NCAL) CALL XABORT('PCRTRP: TERP OVERFLOW(2).')
          IF(ICAL.EQ.0) GO TO 200
          IF(ICAL.EQ.-1) GO TO 210
          TERP(ICAL)=TERP(ICAL)+TERTMP
  100   CONTINUE
      ENDIF
      IF(IMPX.GT.3) THEN
        WRITE(IOUT,'(25H PCRTRP: TERP PARAMETERS:/(1X,1P,10E12.4))')
     1  (TERP(I),I=1,NCAL)
      ENDIF
      RETURN
*----
*  MISSING ELEMENTARY CALCULATION EXCEPTION.
*----
  200 WRITE(IOUT,'(16H PCRTRP: MUPLET=,10I4/(16X,10I4))')
     1 (MUPLET(I),I=1,NPAR)
      CALL XABORT('PCRTRP: MISSING ELEMENTARY CALCULATION.')
  210 WRITE(IOUT,'(16H PCRTRP: MUPLET=,10I4/(16X,10I4))')
     1 (MUPLET(I),I=1,NPAR)
      CALL XABORT('PCRTRP: DEGENERATE ELEMENTARY CALCULATION.')
      END
