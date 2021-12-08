*DECK DOORFVR
      SUBROUTINE DOORFVR(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,NMAT,
     1 IDIR,NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUNKNO,FUNKNO,
     2 IPMACR,REBFLG,DCUTOFF,IPSOUR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the flux. Vectorial version. Multigroup rebalancing
* option.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* CDOOR   name of the geometry/solution operator.
* IPSYS   pointer to the assembly LCM object (L_PIJ signature). IPSYS is
*         a list of directories.
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  unit of the sequential binary tracking file.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* NMAT    number of mixtures in the internal library.
* IDIR    directional collision probability flag:
*         =0 for pij or wij;
*         =k for pijk or wijk k=1,2,3.
*         direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3).  
* NREG    total number of merged blocks for which specific values
*         of the neutron flux and reactions rates are required.
* NUN     total number of unknowns in vectors SUNKNO and FUNKNO.
* IPHASE  type of flux solution (=1: use a native flux solution door;
*         =2: use collision probabilities).
* LEXAC   type of exponential function calculation (=.false. to compute
*         exponential functions using tables).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  index of L-th order flux components in unknown vector.
* TITR    title.
* SUNKNO  input source vector.
* FUNKNO  unknown vector.
* DCUTOFF energy deposition under the energy cutoff.
* IPMACR  pointer to the macrolib LCM object.
* REBFLG  ACA or SCR rebalancing flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK,IPMACR,IPSOUR
      CHARACTER CDOOR*12,TITR*72
      INTEGER NPSYS(NGRP),IFTRAK,IMPX,NGRP,NMAT,IDIR,NREG,NUN,IPHASE,
     > MAT(NREG),KEYFLX(NREG)
      REAL VOL(NREG),SUNKNO(NUN,NGRP),FUNKNO(NUN,NGRP),DCUTOFF(NREG)
      LOGICAL LEXAC,REBFLG
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER IPAR(NSTATE)
*
      IF(IPHASE.EQ.1) THEN
*     USE A NATIVE DOOR
         IF(CDOOR.EQ.'SYBIL') THEN
            CALL SYBILF(IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,IDIR,NREG,
     >                  NMAT,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR)
         ELSE IF(CDOOR.EQ.'BIVAC') THEN
            CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
            NLF=IPAR(14)
            IF(NLF.EQ.0) THEN
               CALL BIVAF(IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,IDIR,NREG,
     >                    NMAT,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR)
            ELSE
               CALL PNF(IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,IDIR,NREG,
     >                  NMAT,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR)
            ENDIF
         ELSE IF(CDOOR.EQ.'TRIVAC') THEN
            CALL TRIVAF(IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,IDIR,NREG,
     >                  NMAT,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR)
         ELSE IF(CDOOR.EQ.'SN') THEN
            CALL SNF(IPSYS,NPSYS,IPTRK,IPSOUR,IFTRAK,IMPX,NGRP,IDIR,
     >               NREG,NMAT,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR,
     >               DCUTOFF)
         ELSE IF(CDOOR.EQ.'MCCG') THEN
            CALL MCCGF(IPSYS,NPSYS,IPTRK,IFTRAK,IPMACR,IMPX,NGRP,IDIR,
     >                 NREG,NMAT,NUN,LEXAC,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,
     >                 TITR,REBFLG)
         ENDIF
      ELSE IF(IPHASE.EQ.2) THEN
         CALL TRFICF(IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,IDIR,NREG,NMAT,
     >               NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR)
      ENDIF
*
      RETURN
      END
