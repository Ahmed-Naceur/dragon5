*DECK PIJABC
      SUBROUTINE PIJABC(NREG,NSOUT,NPRB,SIGTAL,MATRT,PROB,PSST,PSVT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Reconstruct collision probabilities (CP) for  all zones eliminating 
* surfaces from the system.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau, R. Roy, I. petrovic
*
*Parameters: input
* NREG    number of zones for geometry.
* NSOUT   number of surfaces for geometry. 
* NPRB    number of probabilities.
* SIGTAL  albedo-sigt vector.
* MATRT   reflection/transmission vector.
*
*Parameters: output
* PROB    directional cp matrix for all types.
* PSST    PSST=(A**(-1)-PSS)**(-1).
* PSVT    PSVT=PSST*PSV.
*
*-----------------------------------------------------------------------
*

      IMPLICIT         NONE
*----
* VARIABLES
*-----
      INTEGER          NREG,NSOUT,NPRB,NSP1,IPSS,ISUR,ISX,IS,JSX,JS,IER,
     1                 IV,JV,IVSI,IVVI,ISV,IVS,IVV
      INTEGER          MATRT(NSOUT)
*
      REAL             SIGTAL(-NSOUT:NREG)
      DOUBLE PRECISION PROB(NPRB),PSST(NSOUT,NSOUT),PSVT(NSOUT,NREG)
*----
*  EVALUATE MATRIX (A**(-1)-PSS)
*----
      NSP1=NSOUT+1
      IPSS=0
      ISUR=(NSOUT*NSP1)/2
      ISX=0
      DO 100 IS=-NSOUT,-1,1
        ISX=ISX+1
        JSX=0
        ISUR=ISUR+1
        DO 101 JS=-NSOUT,IS,1
          JSX=JSX+1
          IPSS=IPSS+1
          IF((SIGTAL(IS).EQ.0.0).OR.(SIGTAL(JS).EQ.0.0)) THEN
            PSST(ISX,JSX)= 0.0D0
          ELSE
            PSST(ISX,JSX)=-PROB(IPSS)
          ENDIF
          IF(JS.NE.IS) THEN
            PSST(JSX,ISX)=PSST(ISX,JSX)
          ENDIF
 101    CONTINUE
        IF(SIGTAL(IS) .EQ. 0.0)THEN
          PSST(ISX,ISX)=PROB(ISUR)
        ELSE
          JS=-MATRT(-IS)
          IF(JS .EQ. IS) THEN
            PSST(ISX,ISX)=PSST(ISX,ISX)+PROB(ISUR)/SIGTAL(IS)
          ELSE IF(JS .LT. IS) THEN
            JSX=NSOUT+JS+1
            PSST(ISX,JSX)=PSST(ISX,JSX)+PROB(ISUR)/SIGTAL(IS)
            PSST(JSX,ISX)=PSST(ISX,JSX)
          ENDIF
        ENDIF
 100  CONTINUE
*----
*  INVERSE MATRIX PSST=(A**(-1)-PSS)
*----
      CALL ALINVD(NSOUT,PSST,NSOUT,IER)
*----
*  CHECK IF INVERSE IS VALID
*----
      IF(IER .NE. 0 ) CALL XABORT
     >  ('PIJABC: IMPOSSIBLE TO INVERT PSS COUPLING MATRIX')
      IVSI=(NSP1*(NSP1+1))/2
      IVVI=IVSI+NSP1
      DO 110 IV=1,NREG
*----
*    PSVT(IS,IV)=SUM(JSS) PSST(ISS,JSS)*PSV(JSS,IV)
*----
        DO 111 IS=1,NSOUT
          PSVT(IS,IV)=0.0D0
 111    CONTINUE
        DO 120 IS=1,NSOUT
          DO 121 JS=1,NSOUT
            ISV=IVSI+JS
            PSVT(IS,IV)=PSVT(IS,IV)+PSST(IS,JS)*PROB(ISV)
 121      CONTINUE
 120    CONTINUE
        IVV=IVVI
        DO 130 JV=1,IV
          IVV=IVV+1
          ISV=0
          IVS=IVSI
          DO 131 IS=-NSOUT,-1,1
            ISV=ISV+1
            IVS=IVS+1
            IF(SIGTAL(IS).NE.0.0) THEN
              PROB(IVV)=PROB(IVV)+PROB(IVS)*PSVT(ISV,JV)
            ENDIF
 131      CONTINUE
 130    CONTINUE
        IVSI=IVSI+NSP1+IV
        IVVI=IVVI+NSP1+IV
 110  CONTINUE
      RETURN
      END
