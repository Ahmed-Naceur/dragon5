*DECK LCMNAN
      SUBROUTINE LCMNAN(IPLIST)
*
*----------------------------------------------------------------------
*
*Purpose:
* Scan a LCM object for NaN.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIST  address of the table or handle to the XSM file.
*
*----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXLEV=50)
      CHARACTER NAMT*12,MYNAME*12,PATH(MAXLEV)*12,FIRST(MAXLEV)*12,
     1 NAMLCM*12,HSMG*131
      INTEGER IVEC(MAXLEV),KJLON(MAXLEV),IGO(MAXLEV)
      TYPE(C_PTR) KDATA(MAXLEV)
      LOGICAL EMPTY,LCM
*----
*  POINTER VARIABLES
*----
      TYPE(C_PTR) PT_DATA
      REAL, POINTER :: RRR(:)
      DOUBLE PRECISION, POINTER :: DDD(:)
      COMPLEX, POINTER :: CCC(:)
*
      ILEV=1
      KDATA(1)=IPLIST
      KJLON(1)=-1
      IVEC(1)=1
      IGO(1)=5
      CALL LCMVAL(IPLIST,' ')
      CALL LCMINF(IPLIST,NAMLCM,MYNAME,EMPTY,ILONG,LCM)
      IF(EMPTY) GO TO 65
*
*     ASSOCIATIVE TABLE.
   10 NAMT=' '
      CALL LCMNXT(IPLIST,NAMT)
      LENNAM=12
      IF(NAMT.EQ.' ') LENNAM=0
*
      FIRST(ILEV)=NAMT
   20 CALL LCMLEN(IPLIST,NAMT,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 60
      IF(ITYLCM.EQ.0) THEN
*        ASSOCIATIVE TABLE DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMNAN: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',NAMLCM,'''(1).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=-1
         KDATA(ILEV)=LCMGID(IPLIST,NAMT)
         PATH(ILEV)=NAMT
         IPLIST=KDATA(ILEV)
         IVEC(ILEV)=1
         IGO(ILEV)=1
         GO TO 10
      ELSE IF(ITYLCM.EQ.10) THEN
*        LIST DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMNAN: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',NAMLCM,'''(2).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=ILONG
         KDATA(ILEV)=LCMGID(IPLIST,NAMT)
         PATH(ILEV)=NAMT
         IPLIST=KDATA(ILEV)
         IVEC(ILEV)=0
         IGO(ILEV)=2
         GO TO 70
      ELSE IF(ITYLCM.LE.6) THEN
         CALL LCMGPD(IPLIST,NAMT,PT_DATA)
         IF(ITYLCM.EQ.2) THEN
*           SINGLE PRECISION DATA.
            CALL C_F_POINTER(PT_DATA, RRR, (/ ILONG /))
            DO I=1,ILONG
              IF(RRR(I).NE.RRR(I)) THEN
                WRITE(HSMG,'(36HLCMNAN: NAN DETECTED IN REAL ARRAY: ,
     1          A12)') NAMT
                CALL XABORT(HSMG)
              ENDIF
            ENDDO
         ELSE IF(ITYLCM.EQ.4) THEN
*           DOUBLE PRECISION DATA.
            CALL C_F_POINTER(PT_DATA, DDD, (/ ILONG /))
            DO I=1,ILONG
              IF(DDD(I).NE.DDD(I)) THEN
                WRITE(HSMG,'(38HLCMNAN: NAN DETECTED IN DOUBLE ARRAY: ,
     1          A12)') NAMT
                CALL XABORT(HSMG)
              ENDIF
            ENDDO
         ELSE IF(ITYLCM.EQ.6) THEN
*           COMPLEX DATA.
            CALL C_F_POINTER(PT_DATA, CCC, (/ ILONG /))
            DO I=1,ILONG
              IF(CCC(I).NE.CCC(I)) THEN
                WRITE(HSMG,'(39HLCMNAN: NAN DETECTED IN COMPLEX ARRAY: ,
     1          A12)') NAMT
                CALL XABORT(HSMG)
              ENDIF
            ENDDO
         ENDIF
      ELSE
         WRITE(HSMG,'(34HLCMNAN: UNKNOWN TYPE RECORD NAMED ,A12,
     1   5H (1).)') NAMLCM
         CALL XABORT(HSMG)
      ENDIF
      GO TO 60
*
   55 NAMT=PATH(ILEV)
      ILEV=ILEV-1
      IPLIST=KDATA(ILEV)
*
   60 CALL LCMNXT(IPLIST,NAMT)
      IF(NAMT.NE.FIRST(ILEV)) GO TO 20
   65 GO TO (55,55,95,95,100),IGO(ILEV)
*
*     LIST.
   70 IVEC(ILEV)=IVEC(ILEV)+1
      IF(IVEC(ILEV).GT.KJLON(ILEV)) THEN
         GO TO (55,55,95,95,100),IGO(ILEV)
      ENDIF
      CALL LCMLEL(KDATA(ILEV),IVEC(ILEV),ILONG,ITYLCM)
      IF((ILONG.NE.0).AND.(ITYLCM.EQ.0)) THEN
*        ASSOCIATIVE TABLE DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMNAN: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',NAMLCM,'''(3).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=-1
         KDATA(ILEV)=LCMGIL(IPLIST,IVEC(ILEV-1))
         IPLIST=KDATA(ILEV)
         IVEC(ILEV)=1
         IGO(ILEV)=3
         GO TO 10
      ELSE IF((ILONG.NE.0).AND.(ITYLCM.EQ.10)) THEN
*        LIST DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMNAN: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',NAMLCM,'''(4).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=ILONG
         KDATA(ILEV)=LCMGIL(IPLIST,IVEC(ILEV-1))
         IPLIST=KDATA(ILEV)
         IVEC(ILEV)=0
         IGO(ILEV)=4
         GO TO 70
      ELSE IF((ILONG.NE.0).AND.(ITYLCM.LE.6)) THEN
         CALL LCMGPL(IPLIST,IVEC(ILEV),PT_DATA)
         IF(ITYLCM.EQ.2) THEN
*           SINGLE PRECISION DATA.
            CALL C_F_POINTER(PT_DATA, RRR, (/ ILONG /))
            DO I=1,ILONG
              IF(RRR(I).NE.RRR(I)) THEN
                WRITE(HSMG,'(36HLCMNAN: NAN DETECTED IN REAL ARRAY: ,
     1          I12)') IVEC(ILEV)
                CALL XABORT(HSMG)
              ENDIF
            ENDDO
         ELSE IF(ITYLCM.EQ.4) THEN
*           DOUBLE PRECISION DATA.
            CALL C_F_POINTER(PT_DATA, DDD, (/ ILONG /))
            DO I=1,ILONG
              IF(DDD(I).NE.DDD(I)) THEN
                WRITE(HSMG,'(38HLCMNAN: NAN DETECTED IN DOUBLE ARRAY: ,
     1          I12)') IVEC(ILEV)
                CALL XABORT(HSMG)
              ENDIF
            ENDDO
         ELSE IF(ITYLCM.EQ.6) THEN
*           COMPLEX DATA.
            CALL C_F_POINTER(PT_DATA, CCC, (/ ILONG /))
            DO I=1,ILONG
              IF(CCC(I).NE.CCC(I)) THEN
                WRITE(HSMG,'(39HLCMNAN: NAN DETECTED IN COMPLEX ARRAY: ,
     1          I12)') IVEC(ILEV)
                CALL XABORT(HSMG)
              ENDIF
            ENDDO
         ENDIF
      ELSE IF(ILONG.NE.0) THEN
         WRITE(HSMG,'(34HLCMNAN: UNKNOWN TYPE RECORD NAMED ,A12,
     1   5H (2).)') NAMLCM
         CALL XABORT(HSMG)
      ENDIF
      GO TO 70
*
   95 ILEV=ILEV-1
      IPLIST=KDATA(ILEV)
      GO TO 70
*
  100 WRITE(6,'(25H LCMNAN: NO NaN DETECTED.)')
      RETURN
      END

