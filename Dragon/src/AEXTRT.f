*DECK AEXTRT
      SUBROUTINE AEXTRT(AEXFAP,TYPSEG,NBRTYP,IPCHDIM,IPCHTYP,
     1 IPCHDKL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Extraction of component information for a segment of type typseg.
* Component of a FORTRAN-77 emulator of the SAPHYR archive system.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* AEXFAP  external subroutine providing segment type information.
* TYPSEG  character type of segment.
*
*Parameters: output
* NBRTYP  number of components in segment.
* IPCHDIM pointer of vector tchdim containing the
*         dimensions of segment components.
* IPCHTYP pointer of vector tchtyp containing the
*         character types of segment components.
* IPCHDKL pointer of vector tchdkl containing the
*         positions of segment components.
*
*Comments:
* The setara pointers ichdim, ichtyp and ichdkl should be
* deallocated at completion of work.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT INTEGER(A-Z)
      EXTERNAL AEXFAP
      CHARACTER TYPSEG*(*),TEXT4*4,CHAIN*80,HSMG*131
      TYPE(C_PTR) IPCHDIM,IPCHTYP,IPCHDKL
      INTEGER, POINTER, DIMENSION(:) :: ICHDIM,ICHTYP,ICHDKL
      
*
      LONMOT=4
      CHAIN=' '
      IRETOU=1
      IF(TYPSEG(1:1) .EQ. '.') THEN
        CHAIN = TYPSEG(2:)
      ELSE
        CALL AEXFAP(TYPSEG,CHAIN)
      ENDIF
*
* DECODING THE COMPONENT FIELDS.
      MMOT = INDEX(CHAIN, ' ') - 1
      IF(MMOT .GT. 0) THEN
        IRETOU = 0
        NBRTYP = 0
        IP = 0
        DO 10 I = 1, MMOT
          II = INDEX('0123456789', CHAIN(I:I))
          IF(II .NE. 0) THEN
            IP = IP * 10 + II
            IP = IP - 1
          ELSE IF(CHAIN(I:I) .EQ. 'C') THEN
            NBRTYP = NBRTYP + IP + 2
            IP = 0
          ELSE
            NBRTYP = NBRTYP + IP + 1
            IP = 0
          ENDIF
   10   CONTINUE
        NBRTYP = NBRTYP + IP
        MMOTC = MMOT
        IPCHDIM=LCMARA(NBRTYP)
        IPCHTYP=LCMARA(NBRTYP)
        IPCHDKL=LCMARA(NBRTYP)
        CALL C_F_POINTER(IPCHDIM,ICHDIM,(/ NBRTYP /))
        CALL C_F_POINTER(IPCHTYP,ICHTYP,(/ NBRTYP /))
        CALL C_F_POINTER(IPCHDKL,ICHDKL,(/ NBRTYP /))
        IRETOU = 0
        NM = 0
        NN = 0
        NU = 0

        IC = 0
        IDIM = 0
        IP = 0
        JC = 0
        KC = -1

        DO 20 I = 1, MMOTC + 1
          IF(I .LE. MMOTC) THEN
            II = INDEX('1234567890RICLD', CHAIN(I:I))
          ELSE
            II = 1024
          ENDIF
          IF(II .EQ. 0) THEN
            IRETOU = I
            GO TO 30
          ELSE IF(II .LE. 10) THEN
            IP = IP * 10 + II
          ELSE
            IF(JC .NE. 0) THEN
              IF(CHAIN(JC:JC) .EQ. 'C') THEN
                IP = IP + 1
              ENDIF
              IC = IC + 1
              KC = KC + 1
              TEXT4 = CHAIN(JC:JC)
              CALL LCMCAR(TEXT4,.TRUE.,ICHTYP(IC))
              ICHDIM(IC) = IP
              IF(IP .EQ. 0) THEN
                IF(CHAIN(JC:JC) .EQ. 'D') THEN
                 IF(LONMOT .EQ. 8) THEN
                  ICHDKL(IC) = KC
                 ELSE
                  IF(MOD (KC, 2) .EQ. 0) THEN
                    ICHDKL(IC) = KC / 2
                    KC = KC + 1
                  ELSE
                    ICHDKL(IC) = (KC + 1) / 2
                    KC = KC + 2
                  ENDIF
                 ENDIF
                ELSE
                  ICHDKL(IC) = KC
                ENDIF
              ELSE
                IF(NM .EQ. 0) THEN
                  KC = KC - 1
                ENDIF
                ICHDKL(IC) = - KC
                DO JC = 1, IP
                  IC = IC + 1
                  KC = KC + 1
                  TEXT4 = 'I'
                  CALL LCMCAR(TEXT4,.TRUE.,ICHTYP(IC))
                  ICHDIM(IC) = 0
                  ICHDKL(IC) = KC
                ENDDO
              ENDIF

              NN = NN + IP + 1
              IF(IP .GT. 0) THEN
                NM = NM + 1
              ELSE IF(NM .GT. 0) THEN
                IRETOU = I
                GO TO 30
              ELSE
                NU = NU + 1
              ENDIF
              IP = 0
            ELSE IF(IP .NE. 0) THEN
              IRETOU = I
              GO TO 30
            ENDIF

            JC = I
          ENDIF
   20   CONTINUE
        IF(NM .GT. 0) THEN
          CALL LCMCAR(TEXT4,.FALSE.,ICHTYP(NU+1))
          IF(TEXT4(1:1) .EQ. 'C') THEN
            ICHDKL(NU+1) = (NN - 1) * LONMOT
          ELSE IF(TEXT4(1:1) .EQ. 'D') THEN
           IF(LONMOT .EQ. 8) THEN
             ICHDKL(NU+1) = NN - 1
           ELSE
            IF(MOD (NN, 2) .EQ. 0) THEN
              ICHDKL(NU+1) = NN / 2
            ELSE
              ICHDKL(NU+1) = (NN - 1) / 2
            ENDIF
           ENDIF
          ELSE
            ICHDKL(NU+1) = NN - 1
          ENDIF
        ENDIF

        IF(IC .NE. NBRTYP) THEN
          IRETOU = NBRTYP + IC
          GO TO 30
        ENDIF
      ENDIF
*
   30 IF(IRETOU.NE.0) THEN
        WRITE(HSMG,'(22HAEXTRT: FAILURE NUMBER,I5,18H FOR SEGMENT TYPE ,
     1  A)') IRETOU,TYPSEG
        CALL XABORT(HSMG)
      ENDIF
      RETURN
      END
