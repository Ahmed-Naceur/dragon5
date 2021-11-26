*DECK MOCIK3
      SUBROUTINE MOCIK3(NANI,NFUNL,NMOD,ISGNR,KEYANI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Generate all signs ISGNR(L,M,K) for spherical harmonics R(L,M) for
* $0 \\le L \\le$ NANI (and for $-L \\le M \\le L$) on the 8 
* octant angular modes for $1 \\le K \\le 8$.
* All these ISGNR values are compressed to be used according to the
* rectangular dimension.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NANI    scattering anisotropy (=0 for isotropic scattering).
* NFUNL   number of moments of the flux.
* NMOD    first dimension of ISGNR.
*
*Parameters: output
* ISGNR   array of the spherical harmonics signs for the different
*         reflections.
* KEYANI  mode to l index: l=KEYANI(nu).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER  NANI,NFUNL,NMOD,ISGNR(NMOD,NFUNL),KEYANI(NFUNL)
*----
*  LOCAL VARIABLES
*----
      INTEGER  NEWMOD(8,4),K,L,M,IND3,KNEW,NSELEC
      LOGICAL  LROK
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISIWRK
      DATA NEWMOD    / 1,0,0,0, 0,0,2,0,
     >                 1,2,3,4, 0,0,0,0,
     >                 0,0,0,0, 0,0,0,0,
     >                 1,2,3,4, 5,6,7,8 /
*
*     INDEX FOR SIGN ISIWRK
      IND3(L,M,K)=    L*(L+1) + M + 1 + (K-1)*((NANI+1)*(NANI+1))
*
* Definition of signs:
* ISIWRK(L,M,1)= +1
* ISIWRK(L,M,2)= SIGN(M)*(-1)**M
* ISIWRK(L,M,3)= SIGN(M)
* ISIWRK(L,M,4)= (-1)**M
* ISIWRK(L,M,5)= (-1)**(L+M)
* ISIWRK(L,M,6)= SIGN(M)*(-1)**L
* ISIWRK(L,M,7)= SIGN(M)*(-1)**(L+M)
* ISIWRK(L,M,8)= (-1)**L
* where SIGN(M)= +1 for 0 <= M
*                -1 for M <  0
*
      ALLOCATE(ISIWRK(8*(NANI+1)*(NANI+1)))
      DO 20 L= 0, NANI
         DO 10 M= -L, L
            ISIWRK(IND3(L,M,1))= 1
            ISIWRK(IND3(L,M,2))= ISIGN(1,M)*(-1)**M
            ISIWRK(IND3(L,M,3))= ISIGN(1,M)
            ISIWRK(IND3(L,M,4))= (-1)**M
            ISIWRK(IND3(L,M,5))= (-1)**(L+M)
            ISIWRK(IND3(L,M,6))= ISIGN(1,M)*(-1)**L
            ISIWRK(IND3(L,M,7))= ISIGN(1,M)*(-1)**(L+M)
            ISIWRK(IND3(L,M,8))= (-1)**L
   10    CONTINUE
   20 CONTINUE
*
***** SELECTS THE GOOD SIGN ISIWRK(L,M) FUNCTIONS
*             FOR NMOD=2(SLAB),4(TWO-D RECT),8(THREE-D).
*     COMPRESSES ISIWRK INTO ISGNR.
*
      DO 50 K= 1, 8
         NSELEC= 0
         KNEW= NEWMOD(K,NMOD/2)
         IF(KNEW.GT.NMOD) CALL XABORT('MOCIK3: NMOD OVERFLOW')
         IF( KNEW.NE.0 )THEN
            DO 40 L= 0, NANI
               DO 30 M= -L, L
                  LROK=.FALSE.
                  IF( NMOD.EQ.2 )THEN
                     LROK= M.EQ.0
                  ELSEIF( NMOD.EQ.4 )THEN
                     LROK= MOD(L+M,2).EQ.0
                  ELSEIF( NMOD.EQ.8 )THEN
                     LROK= .TRUE.
                  ENDIF
                  IF( LROK )THEN
                     NSELEC= NSELEC+1
                     ISGNR(KNEW,NSELEC)= ISIWRK(IND3(L,M,K))
                     KEYANI(NSELEC) = L
                  ENDIF
   30          CONTINUE
   40       CONTINUE
            IF(NSELEC.NE.NFUNL) CALL XABORT('MOCIK3: INVALID NSELEC')
         ENDIF
   50 CONTINUE
      DEALLOCATE(ISIWRK)
*
      RETURN
      END
