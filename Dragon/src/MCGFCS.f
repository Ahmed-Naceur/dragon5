*DECK MCGFCS
      SUBROUTINE MCGFCS(N,NDIM,NZON,QN,FI,M,NANI,NLIN,NFUNL,SC,S,KPN,
     1                  NREG,IPRINT,KEYFLX,KEYCUR,IBC,SIGAL,STIS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of source for collision at iteration iter.
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
* N       number of spatial unknowns.
* NDIM    number of dimensions for the geometry.
* NZON    index-number of the mixture type assigned to each volume.
* QN      input source (fission-other groups) vector.
* FI      unknown vector.
* M       number of material mixtures.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NLIN    linear discontinuous flag (=1 SC/DD0; =3 LDC/DD1).
* NFUNL   number of spherical harmonics components.
* SC      macroscopic scattering cross section.
* KPN     total number of unknowns in vectors QN and FI.
* NREG    number of volumes.
* IPRINT  print parameter (equal to zero for no print).
* KEYFLX  position of flux elements in FI vector.
* KEYCUR  position of current elements in FI.
* IBC     index for boundary condition to connect a surface to another.
* STIS    integration strategy flag.
* SIGAL   total cross-section and albedo array.
*
*Parameters: output
* S       source elements vector.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,NDIM,NZON(N),M,NANI,NLIN,NFUNL,KPN,NREG,IPRINT,
     1 KEYFLX(NREG,NLIN,NFUNL),KEYCUR(N-NREG),IBC(N-NREG),STIS
      REAL QN(KPN),FI(KPN),SC(0:M,NANI),SIGAL(-6:M)
      DOUBLE PRECISION S(KPN)
*
      IF(NDIM.EQ.2) THEN
*     2D geometry
      DO IR=1,N
         IBM=NZON(IR)
         IF(IBM.LT.0) THEN
*        Boundary condition
            ISUR=IR-NREG
            ISUR2=IBC(ISUR)
            IND=KEYCUR(ISUR)
            IND2=KEYCUR(ISUR2)
            IF(IND.GT.0) S(IND)=SIGAL(IBM)*FI(IND2)
         ELSEIF(IBM.GE.0) THEN
*        Volume cell
            DO IL=0,NANI-1
               XSC=REAL(2*IL+1)*SC(IBM,IL+1)
               DO IM=0,IL
                  DO IE=1,NLIN
                     IND=KEYFLX(IR,IE,1+IL*(IL+1)/2+IM)
                     IF(IND.GT.0) THEN
                       S(IND)=QN(IND)+XSC*FI(IND)
                       IF(STIS.EQ.-1) S(IND)=S(IND)/SIGAL(IBM)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      ELSE ! NDIM.EQ.3   
*     3D geometry     
      DO IR=1,N
         IBM=NZON(IR)
         IF(IBM.LT.0) THEN
*        Boundary condition
            ISUR=IR-NREG
            ISUR2=IBC(ISUR)
            IND=KEYCUR(ISUR)
            IND2=KEYCUR(ISUR2)
            IF(IND.GT.0) S(IND)=SIGAL(IBM)*FI(IND2)
         ELSEIF(IBM.GE.0) THEN
*        Volume cell
            INDA=0
            DO IL=0,NANI-1
               XSC=REAL(2*IL+1)*SC(IBM,IL+1)
               DO IM=-IL,IL
                  DO IE=1,NLIN
                     INDA=INDA+1
                     IND=KEYFLX(IR,IE,INDA)
                     IF(IND.GT.0) THEN
                       S(IND)=QN(IND)+XSC*FI(IND)
                       IF(STIS.EQ.-1) S(IND)=S(IND)/SIGAL(IBM)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      ENDIF
*
      IF(IPRINT.GT.6) CALL PRINDM ('S     ',S,KPN)
      RETURN
      END
