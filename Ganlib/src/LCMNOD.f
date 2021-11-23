*DECK LCMNOD
      SUBROUTINE LCMNOD(NUNIT,IMODE,IDIR,JLONG,ITYLCM,PT_DATA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Export/import a single node. Called by LCMEXP.
*
*Copyright:
* Copyright (C) 1993 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NUNIT   file unit number where the export/import is performed.
* IMODE   type of export/import file (=1 sequential unformatted;
*         =2 sequential formatted -- ascii).
* IDIR    export-import flag(=1 to export ; =2 to import).
* JLONG   node length.
* ITYLCM  node type.
* PT_DATA c_ptr address of data.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER :: NUNIT,IMODE,IDIR,JLONG,ITYLCM
      TYPE(C_PTR) :: PT_DATA
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NBLK=24,LENDAT=4)
      INTEGER, POINTER :: III(:)
      REAL, POINTER :: RRR(:)
      LOGICAL, POINTER :: LLL(:)
      DOUBLE PRECISION, POINTER :: DDD(:)
      COMPLEX, POINTER :: CCC(:)
*
      IF(IDIR.EQ.1) THEN
*        EXPORT A NODE.
         IF(ITYLCM.EQ.1) THEN
*           INTEGER DATA.
            CALL C_F_POINTER(PT_DATA, III, (/ JLONG /))
            DO 40 I=1,1+(JLONG-1)/NBLK
            JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               WRITE(NUNIT) (III((I-1)*NBLK+J),J=1,JMIN)
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               WRITE(NUNIT,'(8I10)') (III((I-1)*NBLK+J),J=1,JMIN)
            ENDIF
   40       CONTINUE
         ELSE IF(ITYLCM.EQ.2) THEN
*           SINGLE PRECISION DATA.
            CALL C_F_POINTER(PT_DATA, RRR, (/ JLONG /))
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               DO 50 I=1,1+(JLONG-1)/NBLK
               JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
               WRITE(NUNIT) (RRR((I-1)*NBLK+J),J=1,JMIN)
   50          CONTINUE
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               WRITE(NUNIT,'(1P,5E16.8)') (RRR(I),I=1,JLONG)
            ENDIF
         ELSE IF(ITYLCM.EQ.3) THEN
*           CHARACTER*4 DATA.
            CALL C_F_POINTER(PT_DATA, III, (/ JLONG /))
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               DO 60 I=1,1+(JLONG-1)/NBLK
               JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
               WRITE(NUNIT) (LENDAT,K=1,JMIN)
   60          CONTINUE
               DO 70 I=1,1+(JLONG-1)/NBLK
               JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
               WRITE(NUNIT) (III((I-1)*NBLK+J),J=1,JMIN)
   70          CONTINUE
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               WRITE(NUNIT,'(8I10)') (LENDAT,I=1,JLONG)
               WRITE(NUNIT,'(20A4)') (III(I),I=1,JLONG)
            ENDIF
         ELSE IF(ITYLCM.EQ.4) THEN
*           DOUBLE PRECISION DATA.
            CALL C_F_POINTER(PT_DATA, DDD, (/ JLONG /))
            DO 90 I=1,1+(JLONG-1)/NBLK
            JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               WRITE(NUNIT) (DDD((I-1)*NBLK+J),J=1,JMIN)
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               WRITE(NUNIT,'(1P,4D20.12)') (DDD((I-1)*NBLK+J),J=1,JMIN)
            ENDIF
   90       CONTINUE
         ELSE IF(ITYLCM.EQ.5) THEN
*           LOGICAL DATA.
            CALL C_F_POINTER(PT_DATA, LLL, (/ JLONG /))
            DO 110 I=1,1+(JLONG-1)/NBLK
            JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               WRITE(NUNIT) (LLL((I-1)*NBLK+J),J=1,JMIN)
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               WRITE(NUNIT,'(8L10)') (LLL((I-1)*NBLK+J),J=1,JMIN)
            ENDIF
  110       CONTINUE
         ELSE IF(ITYLCM.EQ.6) THEN
*           COMPLEX DATA.
            CALL C_F_POINTER(PT_DATA, CCC, (/ JLONG /))
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               DO 120 I=1,1+(JLONG-1)/NBLK
               JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
               WRITE(NUNIT) (CCC((I-1)*NBLK+J),J=1,JMIN)
  120          CONTINUE
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               WRITE(NUNIT,'(1P,5E16.8)') (CCC(I),I=1,JLONG)
            ENDIF
         ENDIF
      ELSE IF(IDIR.EQ.2) THEN
*        IMPORT A NODE.
         IF(ITYLCM.EQ.1) THEN
*           INTEGER DATA.
            CALL C_F_POINTER(PT_DATA, III, (/ JLONG /))
            DO 190 I=1,1+(JLONG-1)/NBLK
            JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               READ(NUNIT) (III((I-1)*NBLK+J),J=1,JMIN)
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               READ(NUNIT,'(8I10)') (III((I-1)*NBLK+J),J=1,JMIN)
            ENDIF
  190       CONTINUE
         ELSE IF(ITYLCM.EQ.2) THEN
*           SINGLE PRECISION DATA.
            CALL C_F_POINTER(PT_DATA, RRR, (/ JLONG /))
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               DO 200 I=1,1+(JLONG-1)/NBLK
               JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
               READ(NUNIT) (RRR((I-1)*NBLK+J),J=1,JMIN)
  200          CONTINUE
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               READ(NUNIT,'(5E16.0)') (RRR(I),I=1,JLONG)
            ENDIF
         ELSE IF(ITYLCM.EQ.3) THEN
*           CHARACTER*4 DATA.
            CALL C_F_POINTER(PT_DATA, III, (/ JLONG /))
            KLONG=0
            DO 215 I=1,1+(JLONG-1)/NBLK
            JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               READ(NUNIT) (III((I-1)*NBLK+J),J=1,JMIN)
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               READ(NUNIT,'(8I10)') (III((I-1)*NBLK+J),J=1,JMIN)
            ENDIF
            DO 210 J=1,JMIN
            KLONG=KLONG+III(J)
  210       CONTINUE
  215       CONTINUE
            IF((KLONG+3)/4.GT.JLONG) CALL XABORT('LCMNOD: OVERFLOW.')
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               DO 220 I=1,1+(JLONG-1)/NBLK
               JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
               READ(NUNIT) (III((I-1)*NBLK+J),J=1,JMIN)
  220          CONTINUE
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               READ(NUNIT,'(20A4)') (III(I),I=1,JLONG)
            ENDIF
         ELSE IF(ITYLCM.EQ.4) THEN
*           DOUBLE PRECISION DATA.
            CALL C_F_POINTER(PT_DATA, DDD, (/ JLONG /))
            DO 230 I=1,1+(JLONG-1)/NBLK
            JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               READ(NUNIT) (DDD((I-1)*NBLK+J),J=1,JMIN)
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               READ(NUNIT,'(4D20.0)') (DDD((I-1)*NBLK+J),J=1,JMIN)
            ENDIF
  230       CONTINUE
         ELSE IF(ITYLCM.EQ.5) THEN
*           LOGICAL DATA.
            CALL C_F_POINTER(PT_DATA, LLL, (/ JLONG /))
            DO 240 I=1,1+(JLONG-1)/NBLK
            JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               READ(NUNIT) (LLL((I-1)*NBLK+J),J=1,JMIN)
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               READ(NUNIT,'(8L10)') (LLL((I-1)*NBLK+J),J=1,JMIN)
            ENDIF
  240       CONTINUE
         ELSE IF(ITYLCM.EQ.6) THEN
*           COMPLEX DATA.
            CALL C_F_POINTER(PT_DATA, CCC, (/ JLONG /))
            IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
               DO 250 I=1,1+(JLONG-1)/NBLK
               JMIN=MIN(NBLK,JLONG-(I-1)*NBLK)
               READ(NUNIT) (CCC((I-1)*NBLK+J),J=1,JMIN)
  250          CONTINUE
            ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
               READ(NUNIT,'(5E16.0)') (CCC(I),I=1,JLONG)
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
