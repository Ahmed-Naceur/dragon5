*DECK DMA
      SUBROUTINE DMA(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set the source of an adjoint fixed source eigenvalue problem. The
* source is the gradient of a macrolib.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): creation type(L_SOURCE);
*         HENTRY(2): read-only type(L_FLUX);
*         HENTRY(3): read-only type(L_MACROLIB);
*         HENTRY(4): read-only type(L_TRACKING).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) IPDMA,IPFLX,IPMAC,IPTRK
      CHARACTER HSIGN*12,TEXT12*12
      INTEGER ISTATE(NSTATE)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NAMEAD,NMIX,NDGRP,IMERGE,
     1 IGCR,MAT,KEY
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.NE.4) CALL XABORT('DMA: FOUR PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('DMA: LI'
     1 //'NKED LIST OR XSM FILE EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0) CALL XABORT('DMA: ENTRY IN CREATE MODE EXPE'
     1 //'CTED.')
      DO I=2,4
         IF((JENTRY(I).NE.2).OR.((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2)))
     1   CALL XABORT('DMA: LINKED LIST OR XSM FILE IN READ-ONLY MODE E'
     2   //'XPECTED AT RHS.')
      ENDDO
      IPDMA=KENTRY(1)
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_FLUX') THEN
         TEXT12=HENTRY(2)
         CALL XABORT('DMA: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_FLUX EXPECTED.')
      ENDIF
      IPFLX=KENTRY(2)
      CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.EQ.'L_MACROLIB') THEN
         IPMAC=KENTRY(3)
      ELSE IF(HSIGN.EQ.'L_LIBRARY') THEN
         IPMAC=LCMGID(KENTRY(3),'MACROLIB')
      ELSE
         TEXT12=HENTRY(3)
         CALL XABORT('DMA: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_MACROLIB OR L_LIBRARY EXPECTED.')
      ENDIF
      CALL LCMGTC(KENTRY(4),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_TRACK') THEN
         TEXT12=HENTRY(4)
         CALL XABORT('DMA: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_TRACK EXPECTED.')
      ENDIF
      IPTRK=KENTRY(4)
*----
*  RECOVER STATE VECTOR INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NREG=ISTATE(1)
      NBMIX=ISTATE(4)
      CALL LCMGET(IPFLX,'STATE-VECTOR',ISTATE)
      NG=ISTATE(1)
      NUN=ISTATE(2)
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.NG) CALL XABORT('DMA: INVALID NUMBER OF GROUPS.')
      NMIL=ISTATE(2)
      NL=ISTATE(3)
      NFM=ISTATE(4)
      NED=ISTATE(5)
      NDEL=ISTATE(7)
      ALLOCATE(NAMEAD(2*NED))
      IF(NED.GT.0) CALL LCMGET(IPMAC,'ADDXSNAME-P0',NAMEAD)
*----
*  READ INPUT PARAMETERS
*----
      ALLOCATE(NMIX(NREG))
      CALL LCMGET(IPTRK,'MATCOD',NMIX)
      CALL DMAGET(IPDMA,NG,NREG,NBMIX,NMIX,IPRINT,NMERGE,NGCOND)
      DEALLOCATE(NMIX)
      NCST=(5+NGCOND*NL+2*NFM*(1+NDEL)+NED)*NMERGE*NGCOND
*----
*  COMPUTE THE GPT SOURCE
*----
      ALLOCATE(NDGRP(NG))
      CALL XDISET(NDGRP,NG,0)
      ALLOCATE(IMERGE(NREG),IGCR(NG),MAT(NREG),KEY(NREG),VOL(NREG))
      CALL LCMGET(IPDMA,'REF:IMERGE',IMERGE)
      CALL LCMGET(IPDMA,'REF:IGCOND',IGCR)
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'KEYFLX',KEY)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      IOF=1
      JJ=IGCR(1)
      DO IND=1,NG
        IF(IND.GT.JJ) THEN
          IOF=IOF+1
          IF(IOF.GT.NGCOND) CALL XABORT('DMA: NGCOND OVERFLOW.')
          JJ=IGCR(IOF)
        ENDIF
        NDGRP(IND)=IOF
      ENDDO
      CALL DMASOU(IPRINT,IPDMA,IPMAC,IPFLX,NG,NREG,NMIL,NL,NDEL,
     1 NED,NAMEAD,NUN,NMERGE,NGCOND,NCST,IMERGE,NDGRP,MAT,KEY,VOL)
      DEALLOCATE(VOL,KEY,MAT,IGCR,IMERGE)
      DEALLOCATE(NDGRP,NAMEAD)
*----
*  SAVE THE SIGNATURE AND STATE VECTOR
*----
      HSIGN='L_SOURCE'
      CALL LCMPTC(IPDMA,'SIGNATURE',12,1,HSIGN)
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NG
      ISTATE(2)=NUN
      ISTATE(3)=0
      ISTATE(4)=NCST
      ISTATE(5)=NMERGE
      ISTATE(6)=NGCOND
      IF(IPRINT.GT.0) WRITE(6,100) (ISTATE(I),I=1,6)
      CALL LCMPUT(IPDMA,'STATE-VECTOR',NSTATE,1,ISTATE)
      RETURN
*
  100 FORMAT(/8H OPTIONS/8H -------/
     1 7H NG    ,I8,28H   (NUMBER OF ENERGY GROUPS)/
     2 7H NUN   ,I8,40H   (NUMBER OF UNKNOWNS PER ENERGY GROUP)/
     3 7H NDIR  ,I8,35H   (NUMBER OF DIRECT FIXED SOURCES)/
     4 7H NCST  ,I8,36H   (NUMBER OF ADJOINT FIXED SOURCES)/
     5 7H NMERGE,I8,34H   (NUMBER OF HOMOGENIZED REGIONS)/
     6 7H NGCOND,I8,38H   (NUMBER OF CONDENSED ENERGY GROUPS))
      END
