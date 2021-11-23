*DECK LIBND6
      SUBROUTINE LIBND6(CFILNA,MAXR,NEL,ITNAM,KPAX,BPAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read depletion data on a NDAS formatted library.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: input
* CFILNA  NDAS file name.
* MAXR    number of reaction types.
* NEL     number of isotopes on library.
*
*Parameters: output
* ITNAM   reactive isotope names in chain.
* KPAX    complete reaction type matrix.
* BPAX    complete branching ratio matrix.
*
*Reference:
* Copyright (C) from NDAS Atomic Energy of Canada Limited utility (2006)
*
*-----------------------------------------------------------------------
*
      USE FSDF
      IMPLICIT NONE
*----
*  Subroutine arguments
*----
      INTEGER MAXR,NEL,ITNAM(3,NEL),KPAX(NEL+MAXR,NEL)
      CHARACTER CFILNA*(*)
      REAL BPAX(NEL+MAXR,NEL)
*----
*  Local variables
*----
      CHARACTER TEXT8*8,TEXT12*12
      INTEGER IND,J,IERR,HEADER(16),IHEAD(200),ISO,JSO,ISOID,NBCHIL,
     > LIBWID
*----
*  INTERNAL PARAMETERS
*   CONVE  : ENERGY CONVERSION FACTOR FROM JOULES/(MOLES*10**-24)
*            TO MEV/NUCLIDE = 1.03643526E+13
*   CONVD  : DECAY CONSTANT CONVERSION FACTOR FROM S**(-1) TO
*            10**(-8)*S**(-1) = 1.0+8
*----
      INTEGER KCAPTU,KDECAY,KFISSP,KN2N
      REAL CONVE,CONVD
      PARAMETER(KCAPTU=3,KDECAY=1,KFISSP=2,KN2N=4,CONVE=1.03643526E+13,
     > CONVD=1.0E+8)
      INTEGER NDECAY
      DOUBLE PRECISION TOTLAM
      EXTERNAL LIBWID
      INTEGER, ALLOCATABLE, DIMENSION(:) :: CHILDR,IWISO
      REAL, ALLOCATABLE, DIMENSION(:) :: BURNDA
*----
*  Scratch storage allocation
*----
      ALLOCATE(CHILDR(2*NEL),IWISO(NEL))
      ALLOCATE(BURNDA(2*NEL))
*----
*  Open and probe the NDAS file
*----
      CALL XSDOPN(CFILNA,IERR)
      IF(IERR.NE.0) THEN
         TEXT12=CFILNA
         CALL XABORT('LIBND6: NDAS library '//TEXT12//' cannot be'//
     >   ' opened')
      ENDIF
      CALL XSDBLD(6001,HEADER,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND6: XSDBLD could not read library'
     > //' parameters')
      ISO=0
      DO IND=1,HEADER(1)
*       Load nuclide header
        CALL XSDISO(7000,6001,IND,IHEAD,IERR)
        NBCHIL=IHEAD(1)
        IF(NBCHIL.GT.NEL) CALL XABORT('LIBND6: Children overflow')
        IF(NBCHIL.NE.0) THEN
          ISO=ISO+1
          IF(ISO.GT.NEL) CALL XABORT('LIBND6: NEL overflow')
          CALL XSDNAM(IND,IWISO(ISO),TEXT8,IERR)
          IF(IERR.NE.0) CALL XABORT('LIBND6: XSDNAM index overflow')
        ENDIF
      ENDDO
      ISO=0
      DO IND=1,HEADER(1)
*       Load nuclide header
        CALL XSDISO(7000,6001,IND,IHEAD,IERR)
        NBCHIL=IHEAD(1)
        IF(NBCHIL.NE.0) THEN
          ISO=ISO+1
          NDECAY=0
          TOTLAM=0.0D0
          CALL XSDNAM(IND,ISOID,TEXT8,IERR)
          READ(TEXT8,'(2A4)') ITNAM(1,ISO),ITNAM(2,ISO)
*         Load burnup children data
          CALL XSDISO(7000,5002,IND,CHILDR,IERR)
*         Load burnup coefficients
          CALL XSDISO(7000,5003,IND,BURNDA,IERR)
          DO J=1,2*NBCHIL,2
            JSO=LIBWID(NEL,IWISO,CHILDR(J))
            IF(CHILDR(J+1).EQ.1) THEN
              IF(JSO.GT.0) THEN
                KPAX(JSO,ISO)=KCAPTU
                BPAX(JSO,ISO)=BURNDA(J)
                KPAX(NEL+KCAPTU,JSO)=1
              ENDIF
              KPAX(NEL+KCAPTU,ISO)=1
            ELSE IF(CHILDR(J+1).EQ.2) THEN
              NDECAY=NDECAY+1
              TOTLAM=TOTLAM+DBLE(BURNDA(J))
              IF(JSO.GT.0) THEN
                KPAX(JSO,ISO)=KDECAY
                BPAX(JSO,ISO)=BURNDA(J)
                KPAX(NEL+KCAPTU,JSO)=1
              ENDIF
              KPAX(NEL+KDECAY,ISO)=1
            ELSE IF(CHILDR(J+1).EQ.3) THEN
              IF(JSO.GT.0) THEN
                KPAX(JSO,ISO)=KFISSP
                BPAX(JSO,ISO)=BURNDA(J)
                KPAX(NEL+KFISSP,JSO)=-1
                KPAX(NEL+KCAPTU,JSO)=1
              ENDIF
            ELSE IF(CHILDR(J+1).EQ.4) THEN
              KPAX(NEL+KFISSP,ISO)=1
              BPAX(NEL+KFISSP,ISO)=BURNDA(J)*CONVE
            ELSE IF(CHILDR(J+1).EQ.5) THEN
              IF(JSO.GT.0) THEN
                KPAX(JSO,ISO)=KN2N
                BPAX(JSO,ISO)=BURNDA(J)
                KPAX(NEL+KCAPTU,JSO)=1
              ENDIF
              KPAX(NEL+KN2N,ISO)=1
            ENDIF
          ENDDO
          IF(NDECAY .EQ. 1) THEN
            BPAX(NEL+KDECAY,ISO)=REAL(TOTLAM)*CONVD
            DO JSO=1,NEL
              IF(KPAX(JSO,ISO).EQ. KDECAY) THEN
                BPAX(JSO,ISO)=1.0
              ENDIF
            ENDDO
          ELSE IF(NDECAY .GT. 1) THEN
            BPAX(NEL+KDECAY,ISO)=REAL(TOTLAM)*CONVD
            DO JSO=1,NEL
              IF(KPAX(JSO,ISO).EQ. KDECAY) THEN
                BPAX(JSO,ISO)=BPAX(JSO,ISO)/REAL(TOTLAM)
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      CALL XSDCL()
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(BURNDA)
      DEALLOCATE(IWISO,CHILDR)
      RETURN
      END
