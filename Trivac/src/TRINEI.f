*DECK TRINEI
      SUBROUTINE TRINEI(IOPT,IDIR,ICAS,ISPLH,ICR,I,KK1,KK2,KK3,KEL,
     >                  IQF,NUM1,NTPH,NTPL,NVT1,NVT2,NVT3,IVAL,KN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the three neighbours of triangle I.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Benaboud
*
*Parameters: input
* IDIR    axis index: W: 1 ; X: 2 ; Y: 3 ; Z: 1.
* ISPLH   used to compute the numbrt of triangles per hexagon using
*         (6*(ISPLH-1)**2).
* ICAS    type of calculation: = 1 (with KK3); = 2 (without KK3).
* I       number of triangles.
* KN      element-ordered unknown list.
*
*Parameters: output
* KK1     first neighbours of triangle I.
* KK2     second neighbours of triangle I.
* KK3     third neighbours of triangle I.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IOPT,IDIR,ICAS,ISPLH,ICR,I,KK1,KK2,KK3,KEL,IQF,NUM1,NTPH,
     > NTPL,NVT1,NVT2,NVT3,IVAL,KN(*)
*----
*  LOCAL VARIABLES
*----
      LOGICAL LPAIR
      INTEGER IPER(180,3),ICF(6,3)
      DATA IPER  /1,2,3,4,5,6, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     >   16,17,18,19,20,21,22,23,24, 1,2,3,4,5,6,7,8,9,10,11,12,13,
     >   14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,
     >   34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,
     >   54, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     >   22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,
     >   42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,
     >   62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,
     >   82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,
     >   2,3,6,1,4,5, 4,5,11,12,19,2,3,9,10,17,18,24,1,7,8,
     >  15,16,22,23,6,13,14,20,21, 6,7,15,16,26,27,38,4,5,13,14,24,
     >  25,36,37,47,2,3,11,12,22,23,34,35,45,46,54,1,9,10,20,21,32,
     >  33,43,44,52,53,8,18,19,30,31,41,42,50,51,17,28,29,39,40,48,
     >  49, 8,9,19,20,32,33,47,48,63,6,7,17,18,30,31,45,46,61,62,76,
     >  4,5,15,16,28,29,43,44,59,60,74,75,87,2,3,13,14,26,27,41,42,
     >  57,58,72,73,85,86,96,1,11,12,24,25,39,40,55,56,70,71,83,84,
     >  94,95,10,22,23,37,38,53,54,68,69,81,82,92,93,21,35,36,51,52,
     >  66,67,79,80,90,91,34,49,50,64,65,77,78,88,89,
     >  3,6,5,2,1,4, 12,19,18,24,23,5,11,10,17,16,22,21,4,3,9,8,15,
     >  14,20,2,1,7,6,13, 27,38,37,47,46,54,53,16,26,25,36,35,45,44,
     >  52,51,7,15,14,24,23,34,33,43,42,50,49,6,5,13,12,22,21,32,31,
     >  41,40,48,4,3,11,10,20,19,30,29,39,2,1,9,8,18,17,28,
     >  48,63,62,76,75,87,86,96,95,33,47,46,61,60,74,73,85,84,94,93,
     >  20,32,31,45,44,59,58,72,71,83,82,92,91,9,19,18,30,29,43,42,
     >  57,56,70,69,81,80,90,89,8,7,17,16,28,27,41,40,55,54, 68,67,79,
     >  78,88,6,5,15,14,26,25,39,38,53,52,66,65,77,4,3,13,12,24,23,
     >  37,36,51,50,64,2,1,11,10,22,21,35,34,49/
      DATA ICF  /6,2,1,5,3,4,1,3,2,6,4,5,2,4,3,1,5,6/
*
      PATRI = REAL(I)/2.
      LPAIR = (AINT(PATRI).EQ.PATRI)
*
      IF(I.LE.NTPL) THEN
         IF(I.EQ.1) THEN
            IQF = ICF(1,IDIR)
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IOPT*NTPH+IQF)
            IF(KK2.GT.0) KK2 = KN((KK2-1)*
     >                          IVAL+IPER(ICR+NVT1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.2) KK3 = KN(NUM1+IPER(ICR+I+NTPL,IDIR))
               IF(ISPLH.GT.2) KK3 = KN(NUM1+IPER(ICR+I+NTPL+1,IDIR))
            ENDIF
         ELSE IF(I.EQ.NTPL) THEN
            IQF = ICF(2,IDIR)
            KK1 = KN(NUM1+IOPT*NTPH+IQF)
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(KK1.GT.0) KK1 = KN((KK1-1)*
     >                          IVAL+IPER(ICR+1+NTPH/2,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.2) KK3 = KN(NUM1+IPER(ICR+I+NTPL,IDIR))
               IF(ISPLH.GT.2) KK3 = KN(NUM1+IPER(ICR+I+NTPL+1,IDIR))
            ENDIF
         ELSE
            IQF = ICF(3,IDIR)
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.2) THEN
                  KK3 = KN(NUM1+IOPT*NTPH+IQF)
                  IF(KK3.GT.0) KK3 = KN((KK3-1)*
     >                               IVAL+IPER(ICR+I+NTPL,IDIR))
               ELSE
                  IF(.NOT.LPAIR) KK3 =KN(NUM1+IPER(ICR+I+NTPL+1,IDIR))
                  IF(LPAIR) THEN
                     KK3 = KN(NUM1+IOPT*NTPH+IQF)
                     IF(KK3.GT.0) KK3 = KN((KK3-1)*
     >                           IVAL+IPER(ICR+I+NTPH-NTPL,IDIR))
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ELSE IF(((I.GT.NTPL).AND.(I.LE.(2*NTPL+2)))
     >                .AND.ISPLH.GE.3) THEN
         IF(I.EQ.(NTPL+1)) THEN
            IQF = ICF(1,IDIR)
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IOPT*NTPH+IQF)
            IF(KK2.GT.0) KK2 = KN((KK2-1)*
     >                              IVAL+IPER(ICR+NVT2,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.3) KK3 = KN(NUM1+IPER(ICR+I+NTPL+2,IDIR))
               IF(ISPLH.GT.3) KK3 = KN(NUM1+IPER(ICR+I+NTPL+3,IDIR))
            ENDIF
         ELSE IF(I.EQ.(2*NTPL+2)) THEN
            IQF = ICF(2,IDIR)
            KK1 = KN(NUM1+IOPT*NTPH+IQF)
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(KK1.GT.0) KK1 = KN((KK1-1)*
     >                            IVAL+IPER(ICR+NVT1+1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.3) KK3 = KN(NUM1+IPER(ICR+I+NTPL+2,IDIR))
               IF(ISPLH.GT.3) KK3 = KN(NUM1+IPER(ICR+I+NTPL+3,IDIR))
            ENDIF
         ELSE
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(.NOT.LPAIR)            KK3 = KN(NUM1+
     >                                        IPER(ICR+I-NTPL-1,IDIR))
               IF(LPAIR.AND.ISPLH.EQ.3) KK3 = KN(NUM1+
     >                                        IPER(ICR+I+NTPL+2,IDIR))
               IF(LPAIR.AND.ISPLH.GT.3) KK3 = KN(NUM1+
     >                                        IPER(ICR+I+NTPL+3,IDIR))
            ENDIF
         ENDIF
      ELSE IF(((I.GT.(2*NTPL+2)).AND.(I.LE.(3*NTPL+6)))
     >                .AND.ISPLH.GE.4) THEN
         IF(I.EQ.(2*NTPL+3)) THEN
            IQF = ICF(1,IDIR)
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IOPT*NTPH+IQF)
            IF(KK2.GT.0) KK2 = KN((KK2-1)*
     >                           IVAL+IPER(ICR+NVT3,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.4) KK3 = KN(NUM1+IPER(ICR+I+NTPL+4,IDIR))
               IF(ISPLH.EQ.5) KK3 = KN(NUM1+IPER(ICR+I+NTPL+5,IDIR))
            ENDIF
         ELSE IF(I.EQ.(3*NTPL+6)) THEN
            IQF = ICF(2,IDIR)
            KK1 = KN(NUM1+IOPT*NTPH+IQF)
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(KK1.GT.0) KK1 = KN((KK1-1)*
     >                         IVAL+IPER(ICR+NVT2+1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.4) KK3 = KN(NUM1+IPER(ICR+I+NTPL+4,IDIR))
               IF(ISPLH.EQ.5) KK3 = KN(NUM1+IPER(ICR+I+NTPL+5,IDIR))
            ENDIF
         ELSE
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(LPAIR)                      KK3 = KN(NUM1+
     >                                        IPER(ICR+I-NTPL-3,IDIR))
               IF(.NOT.LPAIR.AND.ISPLH.EQ.4) KK3 = KN(NUM1+
     >                                        IPER(ICR+I+NTPL+4,IDIR))
               IF(.NOT.LPAIR.AND.ISPLH.EQ.5) KK3 = KN(NUM1+
     >                                        IPER(ICR+I+NTPL+5,IDIR))
            ENDIF
         ENDIF
      ELSE IF(((I.GT.(3*NTPL+6)).AND.(I.LE.(4*NTPL+12)))
     >                .AND.ISPLH.EQ.5) THEN
         IF(I.EQ.(3*NTPL+7)) THEN
            IQF = ICF(1,IDIR)
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IOPT*NTPH+IQF)
            IF(ICAS.EQ.1) KK3 = KN(NUM1+IPER(ICR+I+NTPL+6,IDIR))
            IF(KK2.GT.0) KK2 = KN((KK2-1)*
     >                          IVAL+IPER(ICR+NTPH,IDIR))
         ELSE IF(I.EQ.(4*NTPL+12)) THEN
            IQF = ICF(2,IDIR)
            KK1 = KN(NUM1+IOPT*NTPH+IQF)
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ICAS.EQ.1) KK3 = KN(NUM1+IPER(ICR+I+NTPL+6,IDIR))
            IF(KK1.GT.0) KK1 = KN((KK1-1)*
     >                          IVAL+IPER(ICR+NTPH-NTPL+1,IDIR))
         ELSE
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(LPAIR)       KK3 = KN(NUM1+IPER(ICR+I+NTPL+6,IDIR))
               IF(.NOT.LPAIR) KK3 = KN(NUM1+IPER(ICR+I-NTPL-5,IDIR))
            ENDIF
         ENDIF
      ELSE IF(((I.GT.(NTPH-4*NTPL-12)).AND.(I.LE.(NTPH-3*NTPL-6)))
     >                .AND.ISPLH.EQ.5) THEN
         IF(I.EQ.(NTPH-4*NTPL-11)) THEN
            IQF = ICF(4,IDIR)
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IOPT*NTPH+IQF)
            IF(ICAS.EQ.1) KK3 = KN(NUM1+IPER(ICR+I-NTPL-6,IDIR))
            IF(KK2.GT.0) KK2 = KN((KK2-1)*
     >                          IVAL+IPER(ICR+NTPL,IDIR))
         ELSE IF(I.EQ.(NTPH-3*NTPL-6)) THEN
            IQF = ICF(5,IDIR)
            KK1 = KN(NUM1+IOPT*NTPH+IQF)
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ICAS.EQ.1) KK3 = KN(NUM1+IPER(ICR+I-NTPL-6,IDIR))
            IF(KK1.GT.0) KK1 = KN((KK1-1)*
     >                          IVAL+IPER(ICR+1,IDIR))
         ELSE
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(LPAIR)       KK3 = KN(NUM1+IPER(ICR+I+NTPL+5,IDIR))
               IF(.NOT.LPAIR) KK3 = KN(NUM1+IPER(ICR+I-NTPL-6,IDIR))
            ENDIF
         ENDIF
      ELSE IF(((I.GT.(NTPH-3*NTPL-6)).AND.(I.LE.(NTPH-2*NTPL-2)))
     >                .AND.ISPLH.GE.4) THEN
         IF(I.EQ.(NTPH-3*NTPL-5)) THEN
            IQF = ICF(4,IDIR)
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IOPT*NTPH+IQF)
            IF(KK2.GT.0) KK2 = KN((KK2-1)*
     >                          IVAL+IPER(ICR+NTPH-NVT2,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.4) KK3 = KN(NUM1+IPER(ICR+I-NTPL-4,IDIR))
               IF(ISPLH.EQ.5) KK3 = KN(NUM1+IPER(ICR+I-NTPL-5,IDIR))
            ENDIF
         ELSE IF(I.EQ.(NTPH-2*NTPL-2)) THEN
            IQF = ICF(5,IDIR)
            KK1 = KN(NUM1+IOPT*NTPH+IQF)
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ISPLH.EQ.4) THEN
               IF(KK1.GT.0) KK1 = KN((KK1-1)*
     >                             IVAL+IPER(ICR+1,IDIR))
               IF(ICAS.EQ.1) KK3 = KN(NUM1+IPER(ICR+I-NTPL-4,IDIR))
            ELSE
               IF(KK1.GT.0) KK1 = KN((KK1-1)*
     >                             IVAL+IPER(ICR+NTPL+1,IDIR))
               IF(ICAS.EQ.1) KK3 = KN(NUM1+IPER(ICR+I-NTPL-5,IDIR))
            ENDIF
         ELSE
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(.NOT.LPAIR)            KK3 = KN(NUM1+
     >                                   IPER(ICR+I+NTPL+3,IDIR))
               IF(LPAIR.AND.ISPLH.EQ.4) KK3 = KN(NUM1+
     >                                   IPER(ICR+I-NTPL-4,IDIR))
               IF(LPAIR.AND.ISPLH.EQ.5) KK3 = KN(NUM1+
     >                                   IPER(ICR+I-NTPL-5,IDIR))
            ENDIF
         ENDIF
      ELSE IF(((I.GT.(NTPH-2*NTPL-2)).AND.(I.LE.(NTPH-NTPL)))
     >                .AND.ISPLH.GE.3) THEN
         IF(I.EQ.(NTPH-2*NTPL-1)) THEN
            IQF = ICF(4,IDIR)
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IOPT*NTPH+IQF)
            IF(KK2.GT.0) KK2 = KN((KK2-1)*
     >                          IVAL+IPER(ICR+NTPH-NVT1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.3) KK3 = KN(NUM1+IPER(ICR+I-NTPL-2,IDIR))
               IF(ISPLH.GT.3) KK3 = KN(NUM1+IPER(ICR+I-NTPL-3,IDIR))
            ENDIF
         ELSE IF(I.EQ.(NTPH-NTPL)) THEN
            IQF = ICF(5,IDIR)
            KK1 = KN(NUM1+IOPT*NTPH+IQF)
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ISPLH.EQ.3) THEN
               IF(KK1.GT.0) KK1 = KN((KK1-1)*
     >                             IVAL+IPER(ICR+1,IDIR))
               IF(ICAS.EQ.1) KK3 = KN(NUM1+IPER(ICR+I-NTPL-2,IDIR))
            ELSE
               IF(KK1.GT.0) KK1 = KN((KK1-1)*
     >                       IVAL+IPER(ICR+NTPH-NVT2+1,IDIR))
               IF(ICAS.EQ.1) KK3 = KN(NUM1+IPER(ICR+I-NTPL-3,IDIR))
            ENDIF
         ELSE
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(LPAIR)                      KK3 = KN(NUM1+
     >                                        IPER(ICR+I+NTPL+1,IDIR))
               IF(.NOT.LPAIR.AND.ISPLH.EQ.3) KK3 = KN(NUM1+
     >                                        IPER(ICR+I-NTPL-2,IDIR))
               IF(.NOT.LPAIR.AND.ISPLH.GT.3) KK3 = KN(NUM1+
     >                                        IPER(ICR+I-NTPL-3,IDIR))
            ENDIF
         ENDIF
      ELSE IF(((I.GT.(NTPH-NTPL)).AND.(I.LE.NTPH))) THEN
         IF(I.EQ.(NTPH-NTPL+1)) THEN
            IQF = ICF(4,IDIR)
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IOPT*NTPH+IQF)
            IF(KK2.GT.0) KK2 = KN((KK2-1)*
     >                          IVAL+IPER(ICR+NTPH/2,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.2) KK3 = KN(NUM1+IPER(ICR+I-NTPL,IDIR))
               IF(ISPLH.GT.2) KK3 = KN(NUM1+IPER(ICR+I-NTPL-1,IDIR))
            ENDIF
         ELSE IF(I.EQ.NTPH) THEN
            IQF = ICF(5,IDIR)
            KK1 = KN(NUM1+IOPT*NTPH+IQF)
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(KK1.GT.0) KK1 = KN((KK1-1)*
     >                          IVAL+IPER(ICR+NTPH-NVT1+1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.2) KK3 = KN(NUM1+IPER(ICR+I-NTPL,IDIR))
               IF(ISPLH.GT.2) KK3 = KN(NUM1+IPER(ICR+I-NTPL-1,IDIR))
            ENDIF
         ELSE
            IQF = ICF(6,IDIR)
            KK1 = KN(NUM1+IPER(ICR+I+1,IDIR))
            KK2 = KN(NUM1+IPER(ICR+I-1,IDIR))
            IF(ICAS.EQ.1) THEN
               IF(ISPLH.EQ.2) THEN
                  KK3 = KN(NUM1+IOPT*NTPH+IQF)
                  IF(KK3.GT.0) KK3 = KN((KK3-1)*
     >                               IVAL+IPER(ICR+I-NTPL,IDIR))
               ELSE
                  IF(LPAIR) KK3 = KN(NUM1+IPER(ICR+I-NTPL-1,IDIR))
                  IF(.NOT.LPAIR) THEN
                     KK3 = KN(NUM1+IOPT*NTPH+IQF)
                     IF(KK3.GT.0) KK3 = KN((KK3-1)*
     >                           IVAL+IPER(ICR+I+NTPL-NTPH,IDIR))
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      KEL = KN(NUM1+IPER(ICR+I,IDIR))
      RETURN
      END
