*DECK MCGPJJ
      SUBROUTINE MCGPJJ(IPTRK,IPRINT,NDIM,NANI,MAXNU,NPJJM,KEYANI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the effective number of pjj intermodes to be stored in 2D or
* 3D and the corresponding index when an expansion up to order L 
* of the scattering cross-section is considered in order to
* construct the source term of the scalar flux moments for
* a method of characteristics iteration.
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
* IPTRK  pointer to the tracking (L_TRACK signature).
* IPRINT print flag (> 1 for print).
* NDIM   number of dimensions for the geometry.
* NANI   scattering anisotropy (=1 for isotropic scattering).
* MAXNU  number of angular modes nu=(l,m).
* KEYANI 'mode to l' index: l=KEYANI(nu):
*        Cartesian 2D KEYANI(NU)=INT(0.5*(SQRT(REAL(1+8*(NU-1)))-1.0));
*        Cartesian 3D KEYANI(NU)=INT(SQRT(REAL(NU-1))).
*
*Parameters: output
* NPJJM  number of non-vanishing pjj modes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER IPRINT,NDIM,NANI,MAXNU,NPJJM,KEYANI(MAXNU)
*----
*  LOCAL VARIABLES
*----
      INTEGER MAXPJJM,IL,INU,ILP,INUP,L,LT,DT,NMODE,NMODO,NPJJM0
      LOGICAL EVEN
      CHARACTER CDIM*2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPJJM
*--- 
*  Compute the number of effective modes
*---
      MAXPJJM=MAXNU*(MAXNU+1)/2
      L=NANI-1
      LT=L/2
      DT=L-2*LT
      IF (NDIM.EQ.2) THEN
         CDIM='2D'
         NMODE=(LT+1)**2        ! number of 'even l' modes
         NMODO=(LT+1)*(LT+2*DT) ! number of 'odd l' modes
      ELSE ! NDIM.EQ.3
         CDIM='3D'
         NMODE=(L+1-DT)*(LT+1)
         NMODO=(L+1+DT)*(LT+DT)
      ENDIF
*     'even l' and 'odd l' modes are unconnected
*            (i.e. pjj('even l' <- 'odd l') = 0)
*     and pjj(nu <- nu') = pjj(nu' <- nu)
*     so the effective number of pjj is
      NPJJM=NMODE*(NMODE+1)/2+NMODO*(NMODO+1)/2
*
      IF (IPRINT.GT.1) THEN
         WRITE(*,*) '--------------------'
         WRITE(*,*) 'ANISOTROPY INDEX FOR L=',NANI-1,' IN ',CDIM
         CALL PRINIM('NU->L ',KEYANI,MAXNU)
         WRITE(*,*) '--------------------'
         WRITE(*,*) NPJJM,
     1      ' PJJ(NU <- NU'') MODES OUT OF',MAXPJJM,' TO BE STORED'
      ENDIF
*---
*  Compute and store the indexes for these modes
*---
      ALLOCATE(IPJJM(2*NPJJM))
      NPJJM0=-1
      DO INU=1,MAXNU
         IL=KEYANI(INU)
         DO INUP=1,INU
            ILP=KEYANI(INUP)
            EVEN=(MOD((IL+ILP),2).EQ.0)
            IF (EVEN) THEN
               NPJJM0=NPJJM0+1
               IPJJM(NPJJM0+1)=INU
               IPJJM(NPJJM+NPJJM0+1)=INUP
            ENDIF
         ENDDO
      ENDDO
      NPJJM0=NPJJM0+1
      IF(NPJJM0.NE.NPJJM) CALL XABORT('MCGPIJ: bug.')
      CALL LCMPUT(IPTRK,'PJJIND$MCCG',2*NPJJM,1,IPJJM)
      IF (IPRINT.GT.1) THEN
         WRITE(*,*) 'INDEXES FOR THE CORRESPONDING',NPJJM,' MODES:'
         CALL PRINIM('-> NU ',IPJJM,NPJJM)
         CALL PRINIM('NU''-> ',IPJJM(NPJJM+1),NPJJM)
      ENDIF
      DEALLOCATE(IPJJM)
*
      RETURN
      END
