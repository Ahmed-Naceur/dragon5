*DECK FLDORT
      SUBROUTINE FLDORT(IPSYS,IPFLUX,NUN,NGRP,LMOD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Test the biorthogonality of the direct-CADjoint eigenvectors.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSYS   L_SYSTEM pointer to system matrices.
* IPFLUX  L_FLUX pointer to the solution.
* NUN     number of unknowns in each energy group.
* NGRP    number of energy groups.
* LMOD    number of orthogonal harmonics to compute.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPFLUX
      INTEGER NUN,NGRP,LMOD
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER TEXT12*12,HSMG*131
      TYPE(C_PTR) JPFLUX,KPFLUX,MPFLUX
      REAL, DIMENSION(:), POINTER :: AGARM
      TYPE(C_PTR) AGARM_PTR
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, DIMENSION(:), ALLOCATABLE :: GAR
      COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: CEV,CAD
      COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: DWORK,ORTHO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DWORK(NUN,NGRP),CEV(NUN,NGRP,LMOD),CAD(NUN,NGRP,LMOD),
     1 ORTHO(LMOD,LMOD),GAR(NUN))
*----
*  FLUX RECOVERY
*----
      CALL LCMLEN(IPFLUX,'MODE',ILONG,ITYLCM)
      IF((ILONG.EQ.0).AND.(LMOD.EQ.1)) THEN
        MPFLUX=LCMGID(IPFLUX,'AFLUX')
        DO IGR=1,NGRP
          CALL LCMGDL(MPFLUX,IGR,GAR)
          CAD(:NUN,IGR,1)=GAR(:NUN)
        ENDDO
        MPFLUX=LCMGID(IPFLUX,'FLUX')
        DO IGR=1,NGRP
          CALL LCMGDL(MPFLUX,IGR,GAR)
          CEV(:NUN,IGR,1)=GAR(:NUN)
        ENDDO
      ELSE IF(ILONG.GT.0) THEN
        DO IMOD=1,LMOD
          JPFLUX=LCMGID(IPFLUX,'MODE')
          CALL LCMLEL(JPFLUX,IMOD,ILONG,ITYLCM)
          IF(ILONG.EQ.0) THEN
            WRITE(6,'(20HFLDORT: MISSING MODE,I4,1H.)') IMOD
            CALL XABORT(HSMG)
          ENDIF
          KPFLUX=LCMGIL(JPFLUX,IMOD)
          MPFLUX=LCMGID(KPFLUX,'AFLUX')
          DO IGR=1,NGRP
            CALL LCMLEL(MPFLUX,IGR,ILONG,ITYLCM)
            IF(ITYLCM.EQ.2) THEN
              CALL LCMGDL(MPFLUX,IGR,GAR)
              CAD(:NUN,IGR,IMOD)=GAR(:NUN)
            ELSE IF(ITYLCM.EQ.6) THEN
              CALL LCMGDL(MPFLUX,IGR,CAD(1,IGR,IMOD))
            ENDIF
          ENDDO
          MPFLUX=LCMGID(KPFLUX,'FLUX')
          DO IGR=1,NGRP
            CALL LCMLEL(MPFLUX,IGR,ILONG,ITYLCM)
            IF(ITYLCM.EQ.2) THEN
              CALL LCMGDL(MPFLUX,IGR,GAR)
              CEV(:NUN,IGR,IMOD)=GAR(:NUN)
            ELSE IF(ITYLCM.EQ.6) THEN
              CALL LCMGDL(MPFLUX,IGR,CEV(1,IGR,IMOD))
            ENDIF
          ENDDO
        ENDDO
      ELSE
        CALL XABORT('FLDORT: MODE INFORMATION MISSING.')
      ENDIF
*----
*  MULTIPLY FLUX WITH B MATRIX
*----
      CALL LCMGET(IPSYS,'STATE-VECTOR',ISTATE)
      LL4=ISTATE(2)
      DO JMOD=1,LMOD
        DWORK(:NUN,:NGRP)=0.0D0
        DO IGR=1,NGRP
          DO JGR=1,NGRP
            WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
            CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
            IF(ILONG.EQ.0) CYCLE
            CALL LCMGPD(IPSYS,TEXT12,AGARM_PTR)
            CALL C_F_POINTER(AGARM_PTR,AGARM,(/ ILONG /))
            DO I=1,ILONG
              DWORK(I,IGR)=DWORK(I,IGR)+CMPLX(AGARM(I)*CEV(I,JGR,JMOD))
            ENDDO
          ENDDO
        ENDDO
*----
*  COMPUTE ORTHONORMAL MATRIX
*----
        DO IMOD=1,LMOD
          ORTHO(IMOD,JMOD)=0.0D0
          DO I=1,LL4
            DO IGR=1,NGRP
              ORTHO(IMOD,JMOD)=ORTHO(IMOD,JMOD)+CAD(I,IGR,IMOD)*
     1        DWORK(I,IGR)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
*----
*  PRINT ORTHONORMAL MATRIX
*----
      WRITE(6,'(/28H FLDORT: ORTHONORMAL MATRIX:)')
      DO IMOD=1,LMOD
        WRITE(6,'(3X,1P,15E12.4)') REAL(ORTHO(IMOD,:LMOD))
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR,ORTHO,CAD,CEV,DWORK)
      RETURN
      END
