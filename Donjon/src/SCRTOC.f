*DECK SCRTOC
      SUBROUTINE SCRTOC(IPSAP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print the table of content of a Saphyb.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPSAP   address of the multidimensional Saphyb object.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXLAM=20
      INTEGER, PARAMETER::MAXPAR=50
      INTEGER, PARAMETER::MAXVAL=200
      INTEGER I, ILENG,ILONG, IPAR, ITYLCM, 
     & NADRX, NCALS, NGROUP, NISO, NISOTS, NLAM, NMAC, NMIL, NPAR, 
     & NPARL, NPRC, NREA, NSURFD
      INTEGER DIMSAP(50),NVALUE(MAXPAR),VINTE(MAXVAL)
      REAL VREAL(MAXVAL)
      CHARACTER PARKEY(MAXPAR)*4,PARTYP(MAXPAR)*4,PARFMT(MAXPAR)*8,
     1 VCHAR(MAXVAL)*12,RECNAM*12,NAMLAM(MAXLAM)*8
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: TEXT8
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: TEXT12
*----
*  DIMSAP INFORMATION
*----
      CALL LCMLEN(IPSAP,'DIMSAP',ILENG,ITYLCM)
      IF(ILENG.EQ.0) CALL XABORT('SCRTOC: INVALID SAPHYB.')
      CALL LCMGET(IPSAP,'DIMSAP',DIMSAP)
      NLAM=DIMSAP(3)   ! number of radioactive decay reactions
      NREA=DIMSAP(4)   ! number of neutron-induced reactions
      NISO=DIMSAP(5)   ! number of particularized isotopes
      NMAC=DIMSAP(6)   ! number of macroscopic sets
      NMIL=DIMSAP(7)   ! number of mixtures
      NPAR=DIMSAP(8)   ! number of global parameters
      NPARL=DIMSAP(11) ! number of local variables
      NADRX=DIMSAP(18) ! number of address sets
      NCALS=DIMSAP(19) ! number of elementary calculations in the Saphyb
      NGROUP=DIMSAP(20) ! number of energy groups
      NPRC=DIMSAP(31)   ! number of delayed neutron precursor groups
      NISOTS=DIMSAP(32) ! maximum number of isotopes in output tables
      WRITE(IOUT,'(/38H SCRTOC: table of content information:)')
      WRITE(IOUT,'(42H   number of radioactive decay reactions =,I3)')
     1 NLAM
      WRITE(IOUT,'(40H   number of neutron-induced reactions =,I3)')
     1 NREA
      WRITE(IOUT,'(38H   number of particularized isotopes =,I4)') NISO
      WRITE(IOUT,'(31H   number of macroscopic sets =,I2)') NMAC
      WRITE(IOUT,'(23H   number of mixtures =,I5)') NMIL
      WRITE(IOUT,'(32H   number of global parameters =,I4)') NPAR
      WRITE(IOUT,'(30H   number of local variables =,I4)') NPARL
      WRITE(IOUT,'(27H   number of address sets =,I4)') NADRX
      WRITE(IOUT,'(27H   number of calculations =,I7)') NCALS
      WRITE(IOUT,'(28H   number of energy groups =,I4)') NGROUP
      WRITE(IOUT,'(31H   number of precursor groups =,I4)') NPRC
      WRITE(IOUT,'(48H   maximum number of isotopes in output tables =,
     1 I4/)') NISOTS
      IF(NLAM.GT.0) THEN
        CALL LCMSIX(IPSAP,'constphysiq',1)
          IF(NLAM.GT.MAXLAM) CALL XABORT('SCRTOC: MAXLAM OVERFLOW')
          CALL LCMGTC(IPSAP,'NOMLAM',8,NLAM,NAMLAM)
          WRITE(IOUT,'(40H   names of radioactive decay reactions:/
     1    (5X,5A10))') (NAMLAM(I),I=1,NLAM)
        CALL LCMSIX(IPSAP,' ',2)
      ENDIF
      CALL LCMSIX(IPSAP,'contenu',1)
        IF(NREA.GT.0) THEN
          ALLOCATE(TEXT12(NREA))
          CALL LCMGTC(IPSAP,'NOMREA',12,NREA,TEXT12)
          WRITE(IOUT,'(38H   names of neutron-induced reactions:/
     1    (5X,A12,2X,A12,2X,A12,2X,A12,2X,A12))') (TEXT12(I),I=1,NREA)
          DEALLOCATE(TEXT12)
        ENDIF
        IF(NISO.GT.0) THEN
          ALLOCATE(TEXT8(NISO))
          CALL LCMGTC(IPSAP,'NOMISO',8,NISO,TEXT8)
          WRITE(IOUT,'(36H   names of particularized isotopes:/
     1    (5X,A8,2X,A8,2X,A8,2X,A8,2X,A8))') (TEXT8(I),I=1,NISO)
          DEALLOCATE(TEXT8)
        ENDIF
        IF(NMAC.GT.0) THEN
          ALLOCATE(TEXT8(NMAC))
          CALL LCMGTC(IPSAP,'NOMMAC',8,NMAC,TEXT8)
          WRITE(IOUT,'(29H   names of macroscopic sets:/
     1    (5X,A8,2X,A8,2X,A8,2X,A8,2X,A8))') (TEXT8(I),I=1,NMAC)
          DEALLOCATE(TEXT8)
        ENDIF
      CALL LCMSIX(IPSAP,' ',2)
      CALL LCMSIX(IPSAP,'geom',1)
      CALL LCMLEN(IPSAP,'outgeom',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
        CALL LCMSIX(IPSAP,'outgeom',1)
        CALL LCMLEN(IPSAP,'SURF',NSURFD,ITYLCM)
        WRITE(IOUT,'(36H   number of discontinuity factors =,I4/)')
     1  NSURFD
        CALL LCMSIX(IPSAP,' ',2)
      ENDIF
      CALL LCMSIX(IPSAP,' ',2)
*----
*  GLOBAL PARAMETERS INFORMATION
*----
      IF(NPAR.GT.MAXPAR) CALL XABORT('SCRTOC: MAXPAR OVERFLOW')
      CALL LCMSIX(IPSAP,'paramdescrip',1)
        CALL LCMGET(IPSAP,'NVALUE',NVALUE)
        CALL LCMGTC(IPSAP,'PARKEY',4,NPAR,PARKEY)
        CALL LCMGTC(IPSAP,'PARTYP',4,NPAR,PARTYP)
        CALL LCMGTC(IPSAP,'PARFMT',8,NPAR,PARFMT)
      CALL LCMSIX(IPSAP,' ',2)
      CALL LCMSIX(IPSAP,'paramvaleurs',1)
        DO IPAR=1,NPAR
          WRITE(IOUT,'(25H SCRTOC: global parameter,A5,8H of type,A5,
     1    1H:)') PARKEY(IPAR),PARTYP(IPAR)
          IF(NVALUE(IPAR).GT.MAXVAL) CALL XABORT('SCRTOC: MAXVAL OVERF'
     1    //'LOW')
          WRITE(RECNAM,'(''pval'',I8)') IPAR
          IF(PARFMT(IPAR).EQ.'ENTIER') THEN
            CALL LCMGET(IPSAP,RECNAM,VINTE)
            WRITE(IOUT,'(20H   TABULATED POINTS=,1P,6I12/(20X,6I12))')
     1      (VINTE(I),I=1,NVALUE(IPAR))
          ELSE IF(PARFMT(IPAR).EQ.'FLOTTANT') THEN
            CALL LCMGET(IPSAP,RECNAM,VREAL)
            WRITE(IOUT,'(20H   TABULATED POINTS=,1P,6E12.4/(20X,
     1      6E12.4))') (VREAL(I),I=1,NVALUE(IPAR))
          ELSE IF(PARFMT(IPAR).EQ.'CHAINE') THEN
            CALL LCMGTC(IPSAP,RECNAM,12,NVALUE(IPAR),VCHAR)
            WRITE(IOUT,'(20H   TABULATED POINTS=,2X,6A12/(22X,6A12))')
     1      (VCHAR(I),I=1,NVALUE(IPAR))
          ENDIF
        ENDDO
      CALL LCMSIX(IPSAP,' ',2)
*----
*  LOCAL VARIABLES INFORMATION
*----
      IF(NPARL.GT.0) THEN
        IF(NPARL.GT.MAXPAR) CALL XABORT('SCRTOC: MAXPAR OVERFLOW')
        CALL LCMSIX(IPSAP,'varlocdescri',1)
          CALL LCMGTC(IPSAP,'PARKEY',4,NPARL,PARKEY)
          CALL LCMGTC(IPSAP,'PARTYP',4,NPARL,PARTYP)
          CALL LCMGTC(IPSAP,'PARFMT',8,NPARL,PARFMT)
          DO IPAR=1,NPARL
            WRITE(IOUT,'(23H SCRTOC: local variable,A5,8H of type,A5,
     1      11H and format,A9,1H:)') PARKEY(IPAR),PARTYP(IPAR),
     2      PARFMT(IPAR)
          ENDDO
        CALL LCMSIX(IPSAP,' ',2)
      ENDIF
      WRITE(IOUT,'(/)')
      RETURN
      END
