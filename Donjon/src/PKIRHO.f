*DECK PKIRHO
      SUBROUTINE PKIRHO(IPMAP,NALPHA,T,H,PARAMI,PARAMB,RHO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the reactivity during a Runge-Kutta time step taking into
* account feedback effects.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
* This library is free software; you can redIribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMAP   pointer to the point kinetic directory
* NALPHA  number of feedback parameters
* T       time at beggining of step
* H       time-step duration
* PARAMI  initial values of the global parameters corresponding to
*         RHO=0
* PARAMB  values of global parameters at beginning of stage
*
*Parameters: ouput
* PARAMB  values of global parameters at end of Runge-Kutta time step
* RHO     reactivity during Runge-Kutta time step
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  Subroutine arguments
*----
      TYPE(C_PTR) IPMAP
      INTEGER NALPHA
      REAL PARAMI(NALPHA),PARAMB(NALPHA)
      DOUBLE PRECISION T,H,RHO(3)
*----
*  Local variables
*----
      TYPE(C_PTR) JPPAR,KPPAR
      TYPE(C_PTR) X_PTR,Y_PTR
      DOUBLE PRECISION TS(3),DSUM
      LOGICAL LCUBIC
*----
*  Allocatable arrays
*----
      REAL, POINTER, DIMENSION(:) :: X,Y
      REAL, ALLOCATABLE, DIMENSION(:) :: TERP,GAR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PARAM
*----
*  Compute the values of the global parameters during the time step
*----
      TS(1)=T
      TS(2)=T+H/2.0D0
      TS(3)=T+H
      ALLOCATE(PARAM(3,NALPHA))
      DO IAL=1,NALPHA
         CALL XDRSET(PARAM(1,IAL),3,PARAMB(IAL))
      ENDDO
      JPPAR=LCMGID(IPMAP,'ALPHA')
      DO IAL=1,NALPHA
        KPPAR=LCMGIL(JPPAR,IAL)
        CALL LCMLEN(KPPAR,'TIME-LAW-T',NXY,ITYLCM)
        IF(NXY.NE.0) THEN
          ALLOCATE(TERP(NXY))
          CALL LCMGPD(KPPAR,'TIME-LAW-T',X_PTR)
          CALL C_F_POINTER(X_PTR,X,(/ NXY /))
          CALL LCMGPD(KPPAR,'TIME-LAW-P',Y_PTR)
          CALL C_F_POINTER(Y_PTR,Y,(/ NXY /))
          CALL LCMGET(KPPAR,'TIME-LAW-I',LCUBIC)
          DO I=1,3
            CALL ALTERP(LCUBIC,NXY,X,REAL(TS(I)),.FALSE.,TERP)
            DSUM=0.0D0
            DO J=1,NXY
              DSUM=DSUM+TERP(J)*Y(J)
            ENDDO
            PARAM(I,IAL)=REAL(DSUM)
          ENDDO
          PARAMB(IAL)=PARAM(3,IAL)
          DEALLOCATE(TERP)
        ENDIF
      ENDDO
*----
*  Compute the reactivity
*----
      CALL XDDSET(RHO,3,0.0D0)
      JPPAR=LCMGID(IPMAP,'ALPHA')
      DO IAL=1,NALPHA
        KPPAR=LCMGIL(JPPAR,IAL)
        CALL LCMLEN(KPPAR,'ALPHA-LAW-P',NXY,ITYLCM)
        IF(NXY.NE.0) THEN
          ALLOCATE(TERP(NXY),GAR(NXY))
          CALL LCMGPD(KPPAR,'ALPHA-LAW-P',X_PTR)
          CALL C_F_POINTER(X_PTR,X,(/ NXY /))
          CALL LCMGPD(KPPAR,'ALPHA-LAW-R',Y_PTR)
          CALL C_F_POINTER(Y_PTR,Y,(/ NXY /))
          CALL LCMGET(KPPAR,'ALPHA-LAW-T',ITYPE)
          CALL LCMGET(KPPAR,'ALPHA-LAW-I',LCUBIC)
          DO I=1,3
            IF(ITYPE.EQ.1) THEN
              CALL ALTERP(LCUBIC,NXY,X,PARAM(I,IAL),.FALSE.,TERP)
              CALL ALTERP(LCUBIC,NXY,X,PARAMI(IAL),.FALSE.,GAR)
              DSUM=0.0D0
              DO J=1,NXY
                DSUM=DSUM+(TERP(J)-GAR(J))*Y(J)
              ENDDO
            ELSE IF((ITYPE.EQ.2).AND.(PARAMI(IAL).LT.PARAM(I,IAL))) THEN
              CALL ALTERI(LCUBIC,NXY,X,PARAMI(IAL),PARAM(I,IAL),TERP)
              DSUM=0.0D0
              DO J=1,NXY
                DSUM=DSUM+TERP(J)*Y(J)
              ENDDO
            ELSE IF((ITYPE.EQ.2).AND.(PARAMI(IAL).GT.PARAM(I,IAL))) THEN
              CALL ALTERI(LCUBIC,NXY,X,PARAM(I,IAL),PARAMI(IAL),TERP)
              DSUM=0.0D0
              DO J=1,NXY
                DSUM=DSUM-TERP(J)*Y(J)
              ENDDO
            ELSE IF(ITYPE.EQ.3) THEN
              DO J=1,NXY
                GAR(J)=SQRT(X(J))
              ENDDO
              GAR1=SQRT(PARAMI(IAL))
              GAR2=SQRT(PARAM(I,IAL))
              IF(GAR1.LT.GAR2) THEN
                CALL ALTERI(LCUBIC,NXY,GAR,GAR1,GAR2,TERP)
                DSUM=0.0D0
                DO J=1,NXY
                  DSUM=DSUM+TERP(J)*Y(J)
                ENDDO
              ELSE IF(GAR2.LT.GAR1) THEN
                CALL ALTERI(LCUBIC,NXY,GAR,GAR2,GAR1,TERP)
                DSUM=0.0D0
                DO J=1,NXY
                  DSUM=DSUM-TERP(J)*Y(J)
                ENDDO
              ELSE
                CYCLE
              ENDIF
            ELSE
              CYCLE
            ENDIF
            RHO(I)=RHO(I)+DSUM
           ENDDO
           DEALLOCATE(GAR,TERP)
         ENDIF
      ENDDO
      DEALLOCATE(PARAM)
      RETURN
      END
