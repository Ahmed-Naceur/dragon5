*DECK RODMOD
      SUBROUTINE RODMOD(IPMAP,NGRP,MAXMIX,NCH,NB,RNAME,INS,INSS,
     > HMIX,RMIX,NREB,RODSIZE,INTROD,INI,NUMMIX,NCALL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Modification of the rod data stored in the PARAM folder of a fuel map
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Tixier
*
*Parameters: input
* IPMAP   pointer to the fuel map
* NGRP    number of rod groups
* MAXMIX  maximum number of rod mix
* NCH     number of fuel channels
* NB      number of fuel bundles per channel
* RNAME   name of rod group
* INS     rod insertion for each rod group
* INSS    rod insertion step (in cm)
* HMIX    height of rod mix (if more than 2 rod mix are defined)
* RMIX    number associated to rod mix
* NREB    number of bottom-reflective meshes
* RODSIZE total number of meshes for the fuel zone (=NCH*NB)
* INTROD  final value of data stored in the PARAM folder after 
*         calculation
* INI     initial value of data stored in the PARAM folder (no rod
*         inserted)
* NUMMIX  number of rod mix for each rod group
* NCALL   number to distinguish first or other calls to the ROD: 
*         module
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRP,MAXMIX,INS(NGRP),NREB,NCH,NB,RODSIZE,NUMMIX(NGRP)
      INTEGER INSM,NCALL
      REAL HMIX(NGRP*MAXMIX),RMIX(NGRP*MAXMIX),INTROD(RODSIZE)
      REAL INI,INSS
      CHARACTER(LEN=3) RNAME(NGRP),RASS(NCH)
      TYPE(C_PTR) IPMAP
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT*3
      INTEGER I,J,K,L,M
      REAL FLOT,F1,F2,ICM,HMB,HMT,ZMI1,ZMI2
      DOUBLE PRECISION DFLOT
      TYPE(C_PTR) JPMAP,MPMAP
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: GEOZZ

* Recover axial meshing
      JPMAP=LCMGID(IPMAP,'GEOMAP')
      CALL LCMLEN(JPMAP,'MESHZ',ILONG,ITYLCM)
      ALLOCATE(GEOZZ(ILONG))
      CALL LCMGET(JPMAP,'MESHZ',GEOZZ)

* Recover rod parameters
      MPMAP=LCMGID(IPMAP,'ROD-INFO')
      CALL LCMGET(MPMAP,'ROD-INIT',INI)
      CALL LCMGET(MPMAP,'INS-MAX',INSM)
      CALL LCMGET(MPMAP,'STEP-CM',INSS)
      CALL LCMGET(MPMAP,'REFL-BOTTOM',NREB)      
      CALL LCMGTC(MPMAP,'ROD-NAME',3,NGRP,RNAME)
      CALL LCMGET(MPMAP,'ROD-INSERT',INS)
      CALL LCMGET(MPMAP,'ROD-RIN',RMIX)
      CALL LCMGET(MPMAP,'ROD-NBZONE',NUMMIX)
      CALL LCMGET(MPMAP,'ROD-HEIGHT',HMIX)
      IF(NB+NREB+1.GT.ILONG) CALL XABORT('RODMOD: GEOZZ OVERFLOW.')

      IF(NCALL.EQ.1) THEN
* If it is the first call to the ROD: module, the rod map is stored in
* the fuel map
      N=1
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@RODMOD: CHARACTER DATA EXPECTED.')
      DO WHILE(N.LE.NCH)
        RASS(N)=TEXT
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        N=N+1
      END DO
      CALL LCMPTC(MPMAP,'ROD-MAP',3,NCH,RASS)
      ENDIF
* Recover rod map
      CALL LCMGTC(MPMAP,'ROD-MAP',3,NCH,RASS)
      I=1
      DO WHILE(I.LE.NCH)
        K=NB
        DO WHILE(K.GE.1)
        IF((RASS(I).EQ.'|').OR.(RASS(I).EQ.'-').OR.(RASS(I).EQ.'-|-'
     >  )) THEN
* If no control rod is defined
          INTROD((K-1)*NCH+I)=INI
          GO TO 10
        ELSE
* Recover control rod information
          J=1
          DO WHILE(J.LE.NGRP)
            IF(RASS(I).EQ.RNAME(J)) THEN
              EXIT
            ELSE
              J=J+1
            ENDIF 
          END DO
          M=NUMMIX(J)
          DO WHILE(M.GE.1)
            IF(M.EQ.1) THEN
              ICM=INS(J)*INSS
            ELSE
              ICM=INS(J)*INSS-HMIX(J+(M-2)*NGRP)
            ENDIF
* Mesh size calculations
            HMB=GEOZZ(NB+NREB+1)-GEOZZ(K+NREB)
            HMT=GEOZZ(NB+NREB+1)-GEOZZ(K+NREB+1)
            IF(ICM.LT.0.0) THEN
              M=M-1
            ELSE
              ZMI1=HMB-ICM
              IF(ZMI1.LE.0.0) THEN
                INTROD((K-1)*NCH+I)=RMIX(J+(M-1)*NGRP)
                EXIT
              ELSE
                ZMI2=HMT-ICM
                IF(ZMI2.GE.0.0) THEN
                  IF(M.EQ.1) THEN
                    INTROD((K-1)*NCH+I)=INI
                    GO TO 10
                  ELSE
                    M=M-1
                  ENDIF
                ELSE
* Calculation of the proportion of control rod for the mesh considered
                    F1=HMB-ICM
                    F2=GEOZZ(K+NREB+1)-GEOZZ(K+NREB)-F1
                  IF(M.EQ.1) THEN
                    INTROD((K-1)*NCH+I)=(F2*RMIX(J+(M-1)*NGRP)
     >              +F1*INI)/(F1+F2)
                    GO TO 10
                  ELSE
                    IF(M.GE.3) THEN
* It is not possible to have two interfaces in one mesh
                    IF(((INS(J)*INSS-HMIX(J+(M-3)*NGRP)).GE.HMT).AND.
     >              ((INS(J)*INSS-HMIX(J+(M-3)*NGRP)).LE.HMB))
     >              CALL XABORT('@RODMOD: IT IS NOT POSSIBLE TO HAVE'
     1   //' A MESH WITH MORE THAN TWO INTERFACES. HMIX HAS TO BE'
     1   //' BIGGER')
                    ENDIF
* If two mixtures fill one mesh, we have to compute a fraction of
* insertion for each rod mixture
                    INTROD((K-1)*NCH+I)=(F2*RMIX(J+(M-1)*NGRP)+
     >              F1*RMIX(J+(M-2)*NGRP))/(F1+F2)
                    EXIT
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          END DO
        ENDIF
        K=K-1
        END DO
   10   L=K-1
* If a control rod does not fill one mesh entirely, the meshes located
* at the bottom of this mesh do not contain control rod
        DO WHILE(L.GE.1)
          INTROD((L-1)*NCH+I)=INI
          L=L-1
        END DO
      I=I+1
      END DO
      DEALLOCATE(GEOZZ)
      RETURN
      END
