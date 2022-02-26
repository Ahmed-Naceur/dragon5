*DECK B1TODO
      SUBROUTINE B1TODO(OPTION,TYPE,IMPX,NGRO,IJJ1,NJJ1,IDEL,FLXIN,ST,
     1 SFNU,SCAT1,OLDBIL,B2,D,ALAM1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the leakage coefficients using the Todorova approximation.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* OPTION  type of leakage coefficients; can be 'P0' (P-0), 'P1' (P-1)
*         or 'P0TR' (P-0 with transport correction).
* TYPE    type of iteration; can be 'DIFF' (compute D and exit);
*         'K' (keff search) or 'B' (buckling search).
* IMPX    print flag.
* NGRO    number of groups.
* IJJ1    most thermal group in band for P1 scattering.
* NJJ1    number of groups in band for P1 scattering.
* IDEL    dimension of matrices SCAT0 and SCAT1.
* FLXIN   integrated fluxes.
* ST      total macroscopic cross sections.
* SFNU    nu * macroscopic fission cross-sections.
* SCAT1   packed diffusion P1 macroscopic cross sections.
* OLDBIL  previous norm of the flux.
* B2      buckling.
*
*Parameters: output
* D       diffusion coefficients.
* ALAM1   effective multiplication factor.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER*4 OPTION,TYPE
      INTEGER IMPX,NGRO,IJJ1(NGRO),NJJ1(NGRO),IDEL(2)
      REAL ST(NGRO),SFNU(NGRO),SCAT1(IDEL(2)),D(NGRO),B2
      CHARACTER HSMG*131
      DOUBLE PRECISION FLXIN(NGRO),OLDBIL,ALAM1
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION PROD
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: STOD
*----
*  COMPUTE THE LEAKAGE COEFFICIENTS USING THE TODOROVA APPROXIMATION.
*----
      IF((IMPX.GT.0).AND.(TYPE.EQ.'DIFF')) THEN
         WRITE (6,100)
      ELSE IF(IMPX.GT.0) THEN
         WRITE (6,110) OPTION,TYPE,B2
      ENDIF
      IF((OPTION.EQ.'LKRD').OR.(OPTION.EQ.'RHS')) GO TO 90
*----
*  PROCESS PARTICULAR TYPE='DIFF'
*----
      IF(TYPE.EQ.'DIFF') THEN
         DO IGR=1,NGRO
           BETA=1.0/(3.0*ST(IGR)*ST(IGR)+B2)
           D(IGR)=BETA*ST(IGR)/(1.0-B2*BETA)
         ENDDO
         RETURN
      ENDIF
*
      IF((OPTION.EQ.'P0').OR.(OPTION.EQ.'P0TR')) THEN
        DO IGR=1,NGRO
          D(IGR)=1.0/(3.0*ST(IGR))
        ENDDO
      ELSE IF(OPTION.EQ.'P1') THEN
*       Inscatter approximation
        ALLOCATE(STOD(NGRO,NGRO+1))
        STOD(:NGRO,:NGRO)=0.0D0
        IPOSDE=0
        DO IGR=1,NGRO
          STOD(IGR,IGR)=ST(IGR)
          DO JGR=IJJ1(IGR),IJJ1(IGR)-NJJ1(IGR)+1,-1 ! IGR <- JGR
            IPOSDE=IPOSDE+1
            STOD(IGR,JGR)=STOD(IGR,JGR)-SCAT1(IPOSDE)
          ENDDO
          STOD(IGR,NGRO+1)=FLXIN(IGR)/3.0D0
        ENDDO
        CALL ALSBD(NGRO,1,STOD,IER,NGRO)
        IF(IER.NE.0) CALL XABORT('B1TODO: SINGULAR MATRIX.')
        DO IGR=1,NGRO
          D(IGR)=REAL(STOD(IGR,NGRO+1)/FLXIN(IGR))
        ENDDO
        DEALLOCATE(STOD)
      ELSE
        WRITE(HSMG,'(15HB1TODO: OPTION ,A,25H IS INVALID WITH TODOROVA,
     1  15H APPROXIMATION.)') OPTION
        CALL XABORT(HSMG)
      ENDIF
*----
*  COMPUTE THE EFFECTIVE MULTIPLICATION FACTOR.
*----
   90 IF(TYPE.NE.'K') CALL XABORT('B1TODO: TYPE K EXPECTED.')
      PROD=0.0D0
      IPOSDE=0
      DO IGR=1,NGRO
        PROD=PROD+SFNU(IGR)*FLXIN(IGR)
      ENDDO
      ALAM1=ALAM1*PROD/OLDBIL
      OLDBIL=PROD
      IF(IMPX.GT.1) WRITE (6,120) B2,ALAM1
      RETURN
  100 FORMAT(/43H B1TODO: DIFFUSION COEFFICIENT CALCULATION.)
  110 FORMAT(/21H B1TODO: SOLUTION OF ,A4,21H EQUATIONS WITH TYPE ,A4/
     1 9X,17HINITIAL BUCKLING=,1P,E13.5)
  120 FORMAT(/18H B1TODO: BUCKLING=,1P,E13.5,15H K-EFFECTIVE  =,E13.5)
      END
