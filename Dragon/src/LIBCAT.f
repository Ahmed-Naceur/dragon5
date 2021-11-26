*DECK LIBCAT
      SUBROUTINE LIBCAT (MAXNOR,NPAR,NDIL,DEMT,DEMP,IPRECI,LNORAJ,
     1 SIGERD,SEFFER,NOR,PROSIG,ERRBST)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the weights and base points for a principal cross-section
* type and the partial base points for NPAR partial reactions.
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
* MAXNOR  maximum order of a probability table.
* NPAR    number of partial cross sections types.
* NDIL    number of dilutions used to test the accuracy of the table.
* DEMT    moments of the principal cross section.
* DEMP    moments for each partial cross section.
* IPRECI  accuracy criteria for the table (=1/=2/=3).
* LNORAJ  algorithm flag (=.true.: find an order NOR.le.MAXNOR
*         corresponding to accuracy IPRECI; =.false.: compute the
*         table at order NOR. if this is impossible, try an order
*         smaller than NOR).
* SIGERD  list of dilutions used to test the accuracy of the table.
* SEFFER  list of reference self-shielded cross sections corresponding
*         to each cross-section type and each dilution.
*
*Parameters: input/output
* NOR     input order of the table if LNORAJ=.false.and 
*         output order of the table.
*
*Parameters: output
* PROSIG  probability table.
*         PROSIG(inor,1): weights;
*         PROSIG(inor,2): base points for the principal x-s;
*         PROSIG(inor,3): base points for a partial x-s;
*         etc.
* ERRBST  probability table error.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXNOR,NPAR,NDIL,IPRECI,NOR
      REAL SIGERD(NDIL),SEFFER(NPAR+2,NDIL),PROSIG(MAXNOR,2+NPAR),
     1 ERRBST
      DOUBLE PRECISION DEMT(1-MAXNOR:MAXNOR),
     1                 DEMP(-MAXNOR/2:(MAXNOR-1)/2,NPAR)
      LOGICAL LNORAJ
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SEFF
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: PROSIC
*
      EPS=0.2**(1+IPRECI)
      IF(.NOT.LNORAJ) THEN
*        COMPUTE A TABLE FOR AN IMPOSED VALUE OF NOR.
         IF(NOR.GT.MAXNOR) CALL XABORT('LIBCAT: INVALID INPUT ORDER.')
   10    CALL ALPRTB(NOR,1-NOR,DEMT(1-NOR),IER,PROSIG(1,1),PROSIG(1,2))
         IF(IER.NE.0) THEN
            NOR=NOR-1
            GO TO 10
         ENDIF
         JINI=-NOR/2
         DO 20 IPAR=1,NPAR
         CALL LIBMPA(NOR,JINI,PROSIG(1,1),PROSIG(1,2),DEMP(JINI,IPAR),
     1   PROSIG(1,IPAR+2))
   20    CONTINUE
         RETURN
      ENDIF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PROSIC(MAXNOR,2+NPAR,MAXNOR),SEFF(NPAR+2,NDIL))
*----
* COMPUTE THE TABLE FOR EACH AVAILABLE ORDER. THE MAXIMUM ORDER IS
* LIMITED TO 10 TO AVOID NUMERICAL INSTABILITIES.
*----
      ERRBST=1.0E10
      DO 90 INOR=1,10
      CALL ALPRTB(INOR,1-INOR,DEMT(1-INOR),IER,PROSIC(1,1,INOR),
     1 PROSIC(1,2,INOR))
      IF(IER.NE.0) GO TO 90
      JINI=-INOR/2
      DO 30 IPAR=1,NPAR
      CALL LIBMPA(INOR,JINI,PROSIC(1,1,INOR),PROSIC(1,2,INOR),
     1 DEMP(JINI,IPAR),PROSIC(1,IPAR+2,INOR))
   30 CONTINUE
*----
* COMPUTE THE SELF-SHIELDED CROSS SECTIONS FROM THE TABLE.
*----
      DO 70 IDIL=1,NDIL
      DO 40 IPAR=1,NPAR+2
      SEFF(IPAR,IDIL)=0.0
   40 CONTINUE
      DO 55 IOR=1,INOR
      ASTPD=SIGERD(IDIL)*PROSIC(IOR,1,INOR)/(PROSIC(IOR,2,INOR)+
     1 SIGERD(IDIL))
      SEFF(1,IDIL)=SEFF(1,IDIL)+ASTPD
      DO 50 IPAR=2,NPAR+2
      SEFF(IPAR,IDIL)=SEFF(IPAR,IDIL)+ASTPD*PROSIC(IOR,IPAR,INOR)
   50 CONTINUE
   55 CONTINUE
      DO 60 IPAR=2,NPAR+2
      SEFF(IPAR,IDIL)=SEFF(IPAR,IDIL)/SEFF(1,IDIL)
   60 CONTINUE
   70 CONTINUE
*----
* COMPUTE THE TABLE ACCURACY.
*----
      ERROR=0.0
      DO 85 IDIL=1,NDIL
      DO 80 I=2,NPAR+2
      ERROR=MAX(ERROR,ABS(SEFF(I,IDIL)-SEFFER(I,IDIL))/
     1 ABS(SEFFER(I,NDIL)))
   80 CONTINUE
   85 CONTINUE
      IF(1.2*ERROR.LT.ERRBST) THEN
         NOR=INOR
         IF(ERROR.LT.EPS) THEN
            ERRBST=ERROR
            GO TO 100
         ENDIF
         ERRBST=ERROR
      ENDIF
   90 CONTINUE
*----
* SELECT THE ORDER NOR TABLE.
*----
  100 DO 115 IPAR=1,2+NPAR
      DO 110 I=1,NOR
      PROSIG(I,IPAR)=PROSIC(I,IPAR,NOR)
  110 CONTINUE
  115 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SEFF,PROSIC)
      RETURN
      END
