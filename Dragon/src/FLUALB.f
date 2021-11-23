*DECK FLUALB
      SUBROUTINE FLUALB(IPSYS,NREGIO,NUNKNO,IR,MATCOD,VOLUME,KEYFLX,
     > FUNKNO,SUNKNO,SIGS0,SIGT0,F1,F2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Computes information related to an albedo search.
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
* IPSYS   pointer to the pij object (L_PIJ signature).
* NREGIO  total number of volumes in the domain.
* NUNKNO  number of unknown in the system.
* IR      number of mixtures.
* MATCOD  mixture index in each volume.
* VOLUME  volumes.
* KEYFLX  index pointing to the average fluxes in vector FUNKNO.
* FUNKNO  unknowns.
* SUNKNO  sources.
* SIGS0   within-group scattering macroscopic cross sections of each
*         mixture.
* SIGT0   total macroscopic cross sections of each mixture.
*
*Parameters: output
* F1      first part of the neutron flux.
* F2      second part of the neutron flux.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS
      INTEGER   NREGIO,NUNKNO,IR,MATCOD(NREGIO),KEYFLX(NREGIO)
      REAL      VOLUME(NREGIO),FUNKNO(NUNKNO),SUNKNO(NUNKNO),
     >          SIGS0(0:IR),SIGT0(0:IR),F1(NREGIO),F2(NREGIO)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(NREGIO**2))
*----
*  READ PIS MATRIX
*----
      CALL LCMLEN(IPSYS,'DRAGON-WIS',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.NREGIO) THEN
        CALL LCMGET(IPSYS,'DRAGON-WIS',F2)
      ELSE
        CALL LCMLIB(IPSYS)
        CALL XABORT('FLUALB: THE ALBS OPTION OF THE ASM: MODULE HAVE N'
     >  //'OT BEEN ACTIVATED.')
      ENDIF
*
      ZNUM=0.0
      ZDEN=0.0
      DO 10 I=1,NREGIO
      ZNUM=ZNUM+VOLUME(I)*FUNKNO(KEYFLX(I))
      ZDEN=ZDEN+VOLUME(I)*(SIGT0(MATCOD(I))-SIGS0(MATCOD(I)))*F2(I)
   10 CONTINUE
      ZNUM=-ZNUM/ZDEN
*----
*  READ SCATTERING MODIFIED CP MATRIX.
*----
      CALL LCMLEN(IPSYS,'DRAGON-PCSCT',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.NREGIO**2) THEN
        CALL LCMGET(IPSYS,'DRAGON-PCSCT',WORK)
      ELSE
        CALL LCMLIB(IPSYS)
        CALL XABORT('FLUALB: THE SCATTERING MODIFIED PIJ ARE ABSENT FR'
     >  //'OM LCM.')
      ENDIF
*
      DO 20 I=1,NREGIO
      F2(I)=F2(I)*ZNUM
      F1(I)=0.0
   20 CONTINUE
      DO 40 J=1,NREGIO
      SSS=SUNKNO(KEYFLX(J))
      IOF=(J-1)*NREGIO
      DO 30 I=1,NREGIO
      F1(I)=F1(I)+WORK(IOF+I)*SSS
   30 CONTINUE
   40 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK)
      RETURN
      END
