*DECK MCGACA
      SUBROUTINE MCGACA(LFORW,PACA,N,NG,NFIRST,NGEFF,M,LC,NGIND,NGINDV,
     1                  NCONV,KPSYS,JPMACR,NZON,IPERM,IM,MCU,JU,XIN,
     2                  XOUT,TEMP,LC0,IM0,MCU0)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the product of the left-hand side ACA matrix in its multigroup
* form with a vector and apply group per group left preconditioner.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* LFORW   flag set to .false. to transpose the coefficient matrix.
* PACA    type of preconditioner to solve the ACA corrective system.
* N       number of unknowns per group.
* NG      total number of groups.
* NFIRST  first group to proceed.
* NGEFF   number of  unconverged groups.
* M       number of material mixtures.
* LC      dimension of profiled matrices MCU and CQ.
* NGIND   index of the groups to process.
* NGINDV  index to pass from "NGEFF format" to "NG format".
* NCONV   logical array of convergence status for each group (.TRUE.
*         not converged).
* KPSYS   pointer array for each group properties.
* JPMACR  pointer to the macrolib LCM object ('GROUP' directory).
* NZON    index-number of the mixture type assigned to each volume.
* IPERM   permutation array for ACA.
* IM      connection matrix.
* MCU     connection matrix.
* JU      used for ilu0 preconditioner.
* XIN     undefined.
* LC0     used in ILU0-ACA acceleration.
* IM0     used in ILU0-ACA acceleration.
* MCU0    used in ILU0-ACA acceleration.
*
*Parameters: output
* XOUT    undefined.
*
*Parameters: scratch
* TEMP    preconditioning coefficients.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPSYS(NGEFF),JPMACR
      INTEGER PACA,N,NFIRST,NGEFF,NG,M,LC,NGIND(NGEFF),NGINDV(NG),
     1 NZON(N),IPERM(N),IM(N+1),MCU(LC),JU(N),LC0,IM0(*),MCU0(*)
      DOUBLE PRECISION XIN(N,NGEFF),XOUT(N,NGEFF),TEMP(N)
      LOGICAL LFORW,NCONV(NGEFF)
*----
* LOCAL VARIABLES
*----
      TYPE(C_PTR) JPSYS,KPMACR
      INTEGER I,J,II,IG,JG,JJ,JND,IBM
      REAL, TARGET, SAVE, DIMENSION(1) :: DUMMY
*----
* ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,IJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT
*
      TYPE(C_PTR) DIAGF_PTR,CF_PTR,LUDF_PTR,LUCF_PTR,DIAGQ_PTR,CQ_PTR
      REAL, POINTER, DIMENSION(:) :: DIAGF,CF,LUDF,LUCF,DIAGQ,CQ
*----
*  INITIALIZE POINTERS
*----
      LUDF=>DUMMY
      LUCF=>DUMMY
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NJJ(0:M),IJJ(0:M),IPOS(0:M),XSCAT(0:M*NG))
*
      DO II=NFIRST,NGEFF
         IF(NCONV(II)) THEN
*           compute temp=sum_{g' ne g} Sigma_s^{g<-g'} XIN^{g'}
            IG=NGIND(II)
            JPSYS=KPSYS(II)
            CALL LCMGPD(JPSYS,'DIAGF$MCCG',DIAGF_PTR)
            CALL LCMGPD(JPSYS,'CF$MCCG',CF_PTR)
            CALL C_F_POINTER(DIAGF_PTR,DIAGF,(/ N /))
            CALL C_F_POINTER(CF_PTR,CF,(/ LC /))
            IF(PACA.GE.2) THEN
               CALL LCMGPD(JPSYS,'ILUDF$MCCG',LUDF_PTR)
               CALL C_F_POINTER(LUDF_PTR,LUDF,(/ N /))
               IF(PACA.LT.4) THEN
                  CALL LCMGPD(JPSYS,'ILUCF$MCCG',LUCF_PTR)
                  CALL C_F_POINTER(LUCF_PTR,LUCF,(/ LC /))
               ENDIF
            ENDIF
            KPMACR=LCMGIL(JPMACR,IG)
            CALL LCMGPD(JPSYS,'DIAGQ$MCCG',DIAGQ_PTR)
            CALL LCMGPD(JPSYS,'CQ$MCCG',CQ_PTR)
            CALL C_F_POINTER(DIAGQ_PTR,DIAGQ,(/ N /))
            CALL C_F_POINTER(CQ_PTR,CQ,(/ LC /))
            CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
            CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
            CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
            CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
            DO I=1,N
               TEMP(I)=0.0
               J=IPERM(I)
               IBM=NZON(J)
               IF(IBM.GT.0) THEN
                  JG=IJJ(IBM)
                  DO 10 JND=1,NJJ(IBM)
                  IF(JG.NE.IG) THEN
                     JJ=NGINDV(JG)
                     IF(JJ.GE.NFIRST) THEN
                        TEMP(I)=TEMP(I)+XSCAT(IPOS(IBM)+JND-1)*XIN(I,JJ)
                     ENDIF
                  ENDIF
                  JG=JG-1
 10               CONTINUE
               ENDIF
            ENDDO
*           compute E^{g}*temp
            CALL MCGPRA(LFORW,1,PACA,.FALSE.,N,LC,IM,MCU,JU,DIAGQ,CQ,
     1           LUDF,LUCF,DIAGF,TEMP,XOUT(1,II),LC0,IM0,MCU0,CF)
*           compute D^{g}*XIN^{g}
            CALL MCGPRA(LFORW,1,PACA,.FALSE.,N,LC,IM,MCU,JU,DIAGF,CF,
     1           LUDF,LUCF,DIAGF,XIN(1,II),TEMP,LC0,IM0,MCU0,CF)
*           temp=D^{g}*XIN^{g}-E^{g}*sum_{g' ne g} Sigma_s^{g<-g'} XIN^{g'}
            DO I=1,N
               TEMP(I)=TEMP(I)-XOUT(I,II)
            ENDDO
*           apply single-group preconditioner XOUT^{g}=P^{g}*temp
            CALL MCGPRA(LFORW,2,PACA,.TRUE.,N,LC,IM,MCU,JU,DIAGF,CF,
     1           LUDF,LUCF,DIAGF,XOUT(1,II),TEMP,LC0,IM0,MCU0,CF)
         ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XSCAT,IPOS,IJJ,NJJ)
      RETURN
      END
