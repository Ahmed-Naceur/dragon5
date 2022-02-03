*DECK MCGCDD
      SUBROUTINE MCGCDD(IPRINT,IPMACR,IG,NG,NGEFF,NGIND,NCONV,M,NLONG,
     1                  NUN,NREG,LC,LFORW,PACA,NZON,KEYFLX,KEYCUR,IPERM,
     2                  IM,MCU,JU,EPSINT,MAXINT,FLUX,Q,DIAGQ,CQ,DIAGF,
     3                  CF,ILUDF,ILUCF,LC0,IM0,MCU0)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of the CDD equations (ACA method) for a synthetic diffusion
* flux calculation.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): I. Suslov and R. Le Tellier
*
*Parameters: input
* IPRINT  print parameter.
* IPMACR  pointer to the macrolib LCM object ('GROUP' directory)
* IG      group being processed.
* NG      number of groups.
* NGEFF   number of groups to process.
* NGIND   index of the groups to process.
* NCONV   logical array of convergence status for each group (.TRUE.
*         not converged).
* M       number of material mixtures.
* NLONG   size of the corrective system.
* NUN     number of unknowns per group.
* NREG    number of volume regions.
* LC      dimension of profiled matrices MCU and CQ.
* LFORW   flag set to .false. to transpose the coefficient matrix.
* PACA    type of preconditioner to solve the ACA corrective system.
* NZON    index-number of the mixture type assigned to each volume.
* KEYFLX  position of flux elements in FLUX vector.
* KEYCUR  position of current elements in FLUX vector.
* IPERM   permutation array.
* IM      used in cdd acceleration.
* MCU     used in cdd acceleration.
* JU      used in ACA  acceleration for ilu0.
* EPSINT  stopping criterion for BICGSTAB in ACA resolution.
* MAXINT  maximum number of iterations allowed for BICGSTAB in ACA
*         resolution.
* Q       source vector.
* DIAGQ   used in cdd acceleration.
* CQ      used in cdd acceleration.
* DIAGF   used in cdd acceleration.
* CF      used in cdd acceleration. 
* ILUDF   used in cdd acceleration.
* ILUCF   used in cdd acceleration.
* LC0     used in ILU0-ACA acceleration.
* IM0     used in ILU0-ACA acceleration.
* MCU0    used in ILU0-ACA acceleration.
*
*Parameters: output
* FLUX    zonal scalar flux.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
* SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR
      INTEGER IPRINT,IG,NG,NGEFF,NGIND(NGEFF),M,NLONG,NUN,NREG,LC,PACA,
     1 NZON(NLONG),KEYFLX(NREG),KEYCUR(*),IPERM(NLONG),IM(NLONG+1),
     2 MCU(LC),JU(NLONG),MAXINT,LC0,IM0(*),MCU0(*)
      REAL EPSINT,FLUX(NUN),Q(NUN,NG),DIAGQ(NLONG),CQ(LC),DIAGF(NLONG),
     1 CF(LC),ILUDF(NLONG),ILUCF(LC)
      LOGICAL LFORW,NCONV(NGEFF)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMACR,KPMACR
      DOUBLE PRECISION FF
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,IJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: AR,PHI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      IF(C_ASSOCIATED(IPMACR)) THEN
         JPMACR=LCMGID(IPMACR,'GROUP')
         ALLOCATE(NJJ(0:M),IJJ(0:M),IPOS(0:M),XSCAT(0:M*NG))
      ELSE
         JPMACR=C_NULL_PTR
      ENDIF
      ALLOCATE(PHI(NLONG),AR(NLONG))
*----
* CONSTRUCT RHS OF THE CDD SYSTEM
*----
      DO I=1,NLONG
         J=IPERM(I)
         IF(NZON(J).GE.0) THEN
            FF=DIAGQ(I)*Q(KEYFLX(J),IG)
         ELSE
            FF=0.0
         ENDIF
         DO IL=IM(I)+1,IM(I+1)
         IF(MCU(IL).GT.0) THEN
            L=IPERM(MCU(IL))
            IF(NZON(L).GE.0) FF=FF+CQ(IL)*Q(KEYFLX(L),IG)
         ENDIF
         ENDDO
         PHI(I)=FF
      ENDDO
*----
* INVERSE THE SYSTEM BY THE ITERATIVE METHOD BICGSTAB
*----
*     apply preconditioner to RHS
      CALL MCGPRA(LFORW,2,PACA,.TRUE.,NLONG,LC,IM,MCU,JU,DIAGF,CF,
     1     ILUDF,ILUCF,DIAGF,AR,PHI,LC0,IM0,MCU0,CF)
*
      CALL MCGABG(IPRINT,LFORW,PACA,NLONG,LC,EPSINT,MAXINT,IM,MCU,JU,
     1     DIAGF,CF,ILUDF,ILUCF,AR,PHI,1.0,LC0,IM0,MCU0)
*
      IF(C_ASSOCIATED(JPMACR)) THEN
*---
* MODIFY THE CONTRIBUTION FORM THIS GROUP TO OTHER GROUP ISOTROPIC SOURCES
* (JACOBI -> GAUSS-SEIDEL) 
*---
         DO JJ=1,NGEFF
         IF(NCONV(JJ)) THEN
            JG=NGIND(JJ)
            IF(JG.GT.IG) THEN
            KPMACR=LCMGIL(JPMACR,JG)
            CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
            CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
            CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
            CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
            DO 10 I=1,NLONG
               J=IPERM(I)
               IBM=NZON(J)
               IF(IBM.GT.0) THEN
                  IND=KEYFLX(J)
                  IGG=IJJ(IBM)
                  DO 20 JND=1,NJJ(IBM)
                  IF(IG.EQ.IGG) THEN
                     Q(IND,JG)=Q(IND,JG)+
     1               XSCAT(IPOS(IBM)+JND-1)*(REAL(PHI(I))-FLUX(IND))
                     GOTO 10
                  ENDIF
                  IGG=IGG-1
 20               CONTINUE
               ENDIF
 10         CONTINUE
            ENDIF
         ENDIF
         ENDDO
      ENDIF
*---
* REORDER FLUX VECTOR
*---
      DO I=1,NLONG
         J=IPERM(I)
         IF(NZON(J).GE.0) THEN
            FLUX(KEYFLX(J))=REAL(PHI(I))
         ELSE
            FLUX(KEYCUR(J-NREG))=REAL(PHI(I))
         ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(AR,PHI)
      IF(C_ASSOCIATED(IPMACR)) DEALLOCATE(XSCAT,IPOS,IJJ,NJJ)
      RETURN
      END
