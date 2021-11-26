*DECK SYB7TR
      SUBROUTINE SYB7TR (MNA,NRD,NZIS,NZRS,IFAC,ISYM,NUMREG,ZZIS,ZZRS,
     1 NZIR,NZRR,ZZIR,ZZRR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Unfold the tracking information related to an hexagonal sectorized
* heterogeneous cell.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MNA     number of angles in (0,$\\pi$/6).
* NRD     one plus the number of tubes in the cell.
* NUMREG  tubes indices.
* NZIS    undefined.
* NZRS    undefined.
* IFAC    undefined.
* ISYM    undefined.
* ZZIS    undefined.
* ZZRS    undefined.
*
*Parameters: input/output
* NZRR    length if the original/unfolded real tracking information.
* ZZRR    original/unfolded real tracking information.
* NZIR    length if the original/unfolded integer tracking information.
* ZZIR    original/unfolded integer tracking information.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER        MNA,NRD,NZIS,NZRS,IFAC,ISYM,NUMREG(0:5,NRD),
     &               ZZIS(NZIS),NZIR,NZRR,ZZIR(*)
      REAL           ZZRS(NZRS),ZZRR(*)
*
      IZRR = 0
      IZRS = 0
      IZIR = 0
      IZIS = 0
      DO IA = 1, MNA
        ISS = IFAC
        DO IST = 1, 3
          ISC = 2 * MOD(ISS, 6)
          ZZRR(IZRR+ISC+1) = ZZRS(IZRS+1)
          ZZRR(IZRR+ISC+2) = ISYM * ZZRS(IZRS+2)
          ISC = 2 * MOD(ISS+3, 6)
          ZZRR(IZRR+ISC+1) = ZZRS(IZRS+1)
          ZZRR(IZRR+ISC+2) = ISYM * ZZRS(IZRS+2)
          ISS = ISS + ISYM
          IZRS = IZRS + 2
        ENDDO
      IZRR = IZRR + 12
*
      IZRS = IZRS + 1
      IZRR = IZRR + 1
      W = ZZRS(IZRS)
      ZZRR(IZRR) = W
*
      ISSDEB = IFAC
      IZIS = IZIS + 1
      MNT  = ZZIS(IZIS)
      IZIR = IZIR + 1
      ZZIR(IZIR) = MNT - 1
      IZIR = IZIR + 1
      ZZIR(IZIR) = IFAC * ISYM
*
      DO ITT = 1, MNT
        IZIS = IZIS + 1
        NH  = ZZIS(IZIS)
        IZIS = IZIS + 1
        NX  = ZZIS(IZIS)
        IF (NX .LT. 0) CALL XABORT('SYB7TR: NEGATIVE TRACKS.')
        IF (NH .EQ. 0) THEN
          ISSDEB = ISSDEB + ISYM
        ELSE
          ZZIR(IZIR+2) = NX
          IZIR = IZIR + 2
*
        IZIR = IZIR + 1
        ZZIR(IZIR) = MOD(ISSDEB, 6)
*
        ISS = ISSDEB
        ISR = 0
        DO IHS = 1, NH
          ISP = ISR
          ISR = ZZIS(IZIS+IHS)
          IF (ISR .EQ. ISP) THEN
            ISS = ISS + ISYM
          ENDIF
          ISC = MOD(ISS, 6)
          IRC = NUMREG(ISC, ISR) + 5
          ZZIR(IZIR+IHS) = IRC
        ENDDO
          IZIS = IZIS + NH
*
        DO ITX = 1, NX
          IZRS = IZRS + 1
          IZRR = IZRR + 1
          W = ZZRS(IZRS)
          ZZRR(IZRR) = W
          IRC = 0
        DO IHS = 1, NH
          IRP = IRC
          W   = ZZRS(IZRS+IHS)
          IRC = ZZIR(IZIR+IHS)
          IF (IRC .EQ. IRP ) THEN
            ZZRR(IZRR) = W + ZZRR(IZRR)
          ELSE
            IZRR = IZRR + 1
            ZZRR(IZRR) = W
          ENDIF
        ENDDO
          IZRS = IZRS + NH
        ENDDO
*
          IRC = 0
          JZIR = IZIR
        DO IHS = 1, NH
          IRP = IRC
          IRC = ZZIR(IZIR+IHS)
          IF (IRC .NE. IRP) THEN
            JZIR = JZIR + 1
            ZZIR(JZIR) = IRC
          ENDIF
        ENDDO
          NHR = JZIR - IZIR
          ZZIR(IZIR-2) = NHR
          IZIR = JZIR + 1
          ISC = MOD(ISS, 6)
          ZZIR(IZIR) = ISC
        ENDIF
*
      ENDDO
      ENDDO
*
      NZIR = IZIR
      NZRR = IZRR
*
      RETURN
      END
