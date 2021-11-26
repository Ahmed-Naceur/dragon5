!
!---------------------------------------------------------------------
!
!Purpose:
! To read surfacic file.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
!
!Author(s):
! X. Warin
!
!---------------------------------------------------------------------
!
MODULE SALGET_FUNS_MOD
  
  USE CONSTUTILES, ONLY : FORMATR,FORMATI
  !
  !       GENERIC INTERFACES
  !
  INTERFACE SALGET
     MODULE PROCEDURE &
          SALRIN, SALRIN_0, &
          SALRRE, SALRRE_0, &
          SALRDB, SALRDB_0, &
          SALRCH, SALRCH_0
  END INTERFACE
  !
CONTAINS
  !
  SUBROUTINE SALRIN(DATAIN,N,FIN,FOUT,TEXT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! reading text comment cards followed by N integers
    ! first card before integer data must have a '*' in first column
    !
    !Parameters: input
    ! N         number of integer to be read
    ! FOUT      index of output file (if =0 no printing)
    ! TEXT      description of sought input (for debug)
    !
    !Parameters: input/output
    ! FIN       index of input file (FIN < 0 => do not call SALTIT)
    !
    !Parameters: output
    ! DATAIN    integer array of dimension >= N
    !
    !---------------------------------------------------------------------
    !
    INTEGER,           INTENT(IN)                :: N,FOUT
    INTEGER,           INTENT(INOUT)             :: FIN
    INTEGER,           INTENT(OUT), DIMENSION(N) :: DATAIN
    CHARACTER (LEN=*), INTENT(IN)                :: TEXT
    !**
    INTEGER                           :: I
    LOGICAL                           :: LGFIN
    !**
    LGFIN=FIN.LT.0
    IF(LGFIN)THEN
       FIN=-FIN
    ELSE
       CALL SALTIT('*',FIN,FOUT,TEXT)
    ENDIF
    READ(FIN,*)(DATAIN(I),I=1,N)
    IF(FOUT.NE.0)WRITE(FOUT,'(1X,10I8)')(DATAIN(I),I=1,N)
    IF(LGFIN)FIN=-FIN
    !
  END SUBROUTINE SALRIN
  !
  SUBROUTINE SALRIN_0(DATAIN,FIN,FOUT,TEXT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! reading text comment cards followed by a single integer
    ! first card before integer data must have a '*' in first column
    !
    !Parameters: input
    ! FOUT      index of output file (if =0 no printing)
    ! TEXT      description of sought input (for debug)
    !
    !Parameters: input/output
    ! FIN       index of input file (FIN < 0 => do not call SALTIT)
    !
    !Parameters: output
    ! DATAIN    integer value
    !
    !---------------------------------------------------------------------
    !
    INTEGER,           INTENT(IN)     :: FOUT
    INTEGER,           INTENT(INOUT)  :: FIN
    INTEGER,           INTENT(INOUT)  :: DATAIN
    CHARACTER (LEN=*), INTENT(IN)     :: TEXT
    !****
    LOGICAL                           :: LGFIN
    !****
    LGFIN=FIN.LT.0
    IF(LGFIN)THEN
       FIN=-FIN
    ELSE
       CALL SALTIT('*',FIN,FOUT,TEXT)
    ENDIF
    READ(FIN,*) DATAIN
    IF(FOUT.NE.0)WRITE(FOUT,'(1X,10I8)') DATAIN
    IF(LGFIN)FIN=-FIN
    !
  END SUBROUTINE SALRIN_0
  !
  SUBROUTINE SALRRE(DATAIN,N,FIN,FOUT,TEXT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! reading text comment cards followed by N reals
    ! first card before integer data must have a '*' in first column
    !
    !Parameters: input
    ! N         number of integer to be read (format 1P5E15.6)
    ! FOUT      index of output file (if =0 no printing)
    ! TEXT      description of sought input (for debug)
    !
    !Parameters: input/output
    ! FIN       index of input file (FIN < 0 => do not call SALTIT)
    !
    !Parameters: output
    ! DATAIN    real array of dimension >= N
    !
    !---------------------------------------------------------------------
    !
    INTEGER,           INTENT(IN)                :: N,FIN,FOUT
    REAL,              INTENT(OUT), DIMENSION(N) :: DATAIN
    CHARACTER (LEN=*), INTENT(IN)                :: TEXT
    !****
    INTEGER                                      :: I
    !****
    CALL SALTIT('*',FIN,FOUT,TEXT)
    READ(FIN,*)(DATAIN(I),I=1,N)
    IF(FOUT.NE.0)WRITE(FOUT,'(1X,1P,5'//FORMATR//')')(DATAIN(I),I=1,N)
    !
  END SUBROUTINE SALRRE
  !
  SUBROUTINE SALRRE_0(DATAIN,FIN,FOUT,TEXT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! reading text comment cards followed by a single real
    ! first card before integer data must have a '*' in first column
    !
    !Parameters: input
    ! FOUT      index of output file (if =0 no printing)
    ! TEXT      description of sought input (for debug)
    !
    !Parameters: input/output
    ! FIN       index of input file (FIN < 0 => do not call SALTIT)
    !
    !Parameters: output
    ! DATAIN    real value
    !
    !---------------------------------------------------------------------
    !
    INTEGER,           INTENT(IN)  :: FIN,FOUT
    REAL,              INTENT(OUT) :: DATAIN
    CHARACTER (LEN=*), INTENT(IN)  :: TEXT
    !****
    CALL SALTIT('*',FIN,FOUT,TEXT)
    READ(FIN,*)DATAIN
    IF(FOUT.NE.0)WRITE(FOUT,'(1X,1P,5'//FORMATR//')')DATAIN
    !
  END SUBROUTINE SALRRE_0
  !
  SUBROUTINE SALRDB(DATAIN,N,FIN,FOUT,PREC,TEXT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! reading text comment cards followed by N reals*8
    ! first card before integer data must have a '*' in first column
    !
    !Parameters: input
    ! N         number of real*8 to be read
    ! FOUT      index of output file (if =0 no printing)
    ! PREC      = 0 read with 4E20.12, otherwise use 5E12.6
    ! TEXT      description of sought input (for debug)
    !
    !Parameters: input/output
    ! FIN       index of input file (FIN < 0 => do not call SALTIT)
    !
    !Parameters: output
    ! DATAIN    real*8 array of dimension >= N
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS, ONLY : PDB
    !****
    INTEGER,           INTENT(IN)                :: N,FIN,FOUT,PREC
    REAL(PDB),         INTENT(OUT), DIMENSION(N) :: DATAIN
    CHARACTER (LEN=*), INTENT(IN)                :: TEXT
    !****
    INTEGER                                      :: I
    !****
    CALL SALTIT('*',FIN,FOUT,TEXT)
    IF(PREC.EQ.0)THEN
       READ(FIN,'(4E20.12)')(DATAIN(I),I=1,N)
       IF(FOUT.NE.0)WRITE(FOUT,'(1X,1P,4E20.12)')(DATAIN(I),I=1,N)
    ELSE
       READ(FIN,*)(DATAIN(I),I=1,N)
       IF(FOUT.NE.0)WRITE(FOUT,'(1X,1P,5'//FORMATR//')')(DATAIN(I),I=1,N)
    ENDIF
    !
  END SUBROUTINE SALRDB
  !
  SUBROUTINE SALRDB_0(DATAIN,FIN,FOUT,PREC,TEXT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! reading text comment cards followed by a single real*8
    ! first card before integer data must have a '*' in first column
    !
    !Parameters: input
    ! FOUT      index of output file (if =0 no printing)
    ! PREC      = 0 read with 4E20.12, otherwise use 5E12.6
    ! TEXT      description of sought input (for debug)
    !
    !Parameters: input/output
    ! FIN       index of input file (FIN < 0 => do not call SALTIT)
    !
    !Parameters: output
    ! DATAIN    real*8 value
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS, ONLY : PDB
    !****
    INTEGER,           INTENT(IN)  :: FIN,FOUT,PREC
    REAL(PDB),         INTENT(OUT) :: DATAIN
    CHARACTER (LEN=*), INTENT(IN)  :: TEXT
    !****
    !****
    CALL SALTIT('*',FIN,FOUT,TEXT)
    IF(PREC.EQ.0)THEN
       READ(FIN,'(4E20.12)')DATAIN
       IF(FOUT.NE.0)WRITE(FOUT,'(1X,1P,4E20.12)')DATAIN
    ELSE
       READ(FIN,*)DATAIN
       IF(FOUT.NE.0)WRITE(FOUT,'(1X,1P,5'//FORMATR//')')DATAIN
    ENDIF
    !
  END SUBROUTINE SALRDB_0
  !
  SUBROUTINE SALRCH(DATAIN,N,FIN,FOUT,TEXT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! reading text comment cards followed by N chains
    ! first card before integer datain must have a '*' in first column
    !
    !Parameters: input
    ! N         number of strings to be read (3X,A12)
    ! FOUT      index of output file (if =0 no printing)
    ! TEXT      description of sought input (for debug)
    !
    !Parameters: input/output
    ! FIN       index of input file (FIN < 0 => do not call SALTIT)
    !
    !Parameters: output
    ! DATAIN    character array of dimension >= N
    !
    !---------------------------------------------------------------------
    !
    INTEGER,           INTENT(IN)                :: N,FIN,FOUT
    CHARACTER (LEN=*), INTENT(OUT), DIMENSION(N) :: DATAIN
    CHARACTER (LEN=*), INTENT(IN)                :: TEXT
    !****
    INTEGER                                      :: I
    !****
    CALL SALTIT('*',FIN,FOUT,TEXT)
    READ(FIN,'(4(3X,A12))')(DATAIN(I),I=1,N)
    IF(FOUT.NE.0)WRITE(FOUT,'(1X,4(3X,A12))')(DATAIN(I),I=1,N)
    !
  END SUBROUTINE SALRCH
  !
  SUBROUTINE SALRCH_0(DATAIN,FIN,FOUT,TEXT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! reading a single chain of characters.
    ! first card before integer datain must have a '*' in first column
    !
    !Parameters: input
    ! FOUT      index of output file (if =0 no printing)
    ! TEXT      description of sought input (for debug)
    !
    !Parameters: input/output
    ! FIN       index of input file (FIN < 0 => do not call SALTIT)
    !
    !Parameters: output
    ! DATAIN    character chain
    !
    !---------------------------------------------------------------------
    !
    INTEGER,           INTENT(IN)  :: FIN,FOUT
    CHARACTER (LEN=*), INTENT(OUT) :: DATAIN
    CHARACTER (LEN=*), INTENT(IN)  :: TEXT
    !****
    !****
    CALL SALTIT('*',FIN,FOUT,TEXT)
    READ(FIN,'(4(3X,A12))')DATAIN
    IF(FOUT.NE.0)WRITE(FOUT,'(1X,4(3X,A12))')DATAIN
    !
  END SUBROUTINE SALRCH_0
  !
  SUBROUTINE SALTIT(WORD,FIN,FOUT,SEEK)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! reads and prints lines of length 80 char until the line that
    ! begins with word or with 'end'
    !
    !Parameters: input
    ! WORD   character string (mask)
    ! FIN    logical number of input file
    ! FOUT   logical number of output file (if =0 no printing)
    ! SEEK   description of sought input (for debug)
    !
    !---------------------------------------------------------------------
    !
    !**    reads and prints lines of length 80 char until the line that
    !      begins with word or with 'end'
    !>     WORD     = character string (mask)
    !>     FIN      = logical number of input file
    !>     FOUT     = logical number of output file (if =0 no printing)
    !>     SEEK     = description of sought input (for debug)
    !**
    INTEGER,           INTENT(IN)     :: FOUT,FIN
    CHARACTER (LEN=*), INTENT(IN)     :: WORD,SEEK
    !**
    CHARACTER (LEN=80)    :: TEXT
    INTEGER               :: LL
    LOGICAL               :: LGOUT
    !**
    LL=LEN(WORD)
    DO
       READ(FIN,'(A80)')TEXT
       IF(TEXT(1:3).EQ.'END')THEN
          IF(FOUT.NE.0)WRITE(FOUT,'(1X,''SEEKS => '',A)')SEEK
          WRITE(FOUT,'(5X,''READ END AND STOP IN TITLE'')')
          CALL XABORT('SALTIT: FAILURE')
       ENDIF
       IF(TEXT(1:5).EQ.'%SKIP')THEN
          DO
             READ(FIN,'(A80)')TEXT
             IF(FOUT.NE.0)WRITE(FOUT,'(1X,A80)')TEXT
             IF(TEXT(1:5).EQ.'%SKIP')EXIT
          ENDDO
       ELSE
          LGOUT=TEXT(1:LL).EQ.WORD
          IF(FOUT.NE.0)THEN
             IF(LGOUT)WRITE(FOUT,'(1X,''SEEKS => '',A)')SEEK
             WRITE(FOUT,'(1X,A80)')TEXT
          ENDIF
          IF(LGOUT)EXIT
       ENDIF
    ENDDO
    !
  END SUBROUTINE SALTIT
END MODULE SALGET_FUNS_MOD
