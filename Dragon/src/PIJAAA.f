*DECK PIJAAA
      SUBROUTINE PIJAAA(NREG,NSOUT,SIGTAL,PROB,PSVT,PROBS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculates directional collision probabilities for all zones 
* eliminating surfaces from the system: 
* PIJK=PIJK+PISK*((1-PSS)**(-1))*PSJ.
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): I. petrovic
*
*Parameters: input
* NREG    number of zones for geometry.
* NSOUT   number of surfaces for geometry. 
* SIGTAL  albedo-sigt vector.
* PROB    directional cp matrix for all types.
* PSVT    PSST matrix:
*         PSVT=(A**(-1)-PSS)**(-1)*PSV.
*
*Parameters: output
* PROBS   directional CP matrix
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* INTERFACE VARIABLES
*----
      INTEGER          NREG,NSOUT
      REAL             SIGTAL(-NSOUT:NREG),PROBS(*)
      DOUBLE PRECISION PROB(*),PSVT(NSOUT,NREG)
*----
* LOCAL VARIABLES
*----
      INTEGER NSP1,IVSI,IDPSV,IV,IPRL,IPRU,JV,ISV,IPSV,IVS,ISU
*
      NSP1=NSOUT+1
      IVSI=(NSP1*(NSP1+1))/2
      IDPSV=IVSI
      DO 100 IV=1,NREG
        IPRL=NREG*(IV-1)+1
        IPRU=IV
        IPSV=IDPSV
        DO 110 JV=1,IV
          ISV=0
          IVS=IVSI
          DO 120 ISU=-NSOUT,-1,1
            ISV=ISV+1
            IVS=IVS+1
            IPSV=IPSV+1
            IF(SIGTAL(ISU).NE.0.0) THEN
              PROBS(IPRL)=PROBS(IPRL)+REAL(PROB(IVS)*PSVT(ISV,JV))
              IF(IPRL.NE.IPRU) THEN
                PROBS(IPRU)=PROBS(IPRU)+REAL(PROB(IPSV)*PSVT(ISV,IV))
              ENDIF
            ENDIF
 120      CONTINUE
          IPSV=IPSV+JV+1
          IPRL=IPRL+1
          IPRU=NREG*JV+IV
 110    CONTINUE
        IVSI=IVSI+NSP1+IV
 100  CONTINUE
*
      RETURN
      END
