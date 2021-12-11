*DECK MCGFLS
      SUBROUTINE MCGFLS(IMPX,IPTRK,IPMACR,NUN,K,NREG,NLONG,M,NG,NGEFF,
     1                  LC,LFORW,PACA,NZON,KEYFLX,KEYCUR,NGIND,KPSYS,
     2                  NCONV,EPSI,MAXI,FIMEM,QFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Synthetic diffusion (ACA) flux calculation.
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
* IMPX   print flag (equal to zero for no print).
* IPTRK   pointer to the tracking (L_TRACK signature).
* IPMACR  pointer to the macrolib LCM object.
* NUN     number of unknowns per group.
* K       number of volumes and outer surfaces.
* NREG    number of volume regions.
* NLONG   size of the corrective system.
* M       number of material mixtures.
* NG      number of groups.
* NGEFF   number of groups to process.
* LC      dimension of profiled matrices MCU and CQ.
* LFORW   flag set to .false. to transpose the coefficient matrix.
* PACA    type of preconditioner to solve the ACA corrective system.
* NZON    index-number of the mixture type assigned to each volume.
* KEYFLX  position of flux elements in FIMEM vector.
* KEYCUR  position of current elements in FIMEM vector.
* NGIND   index of the groups to process.
* KPSYS   pointer to system groups.
* NCONV   array of convergence flag for each group.
* EPSI    stopping criterion for BICGSTAB in ACA resolution.
* MAXI    maximum number of iterations allowed for BICGSTAB in ACA
*         resolution.
* QFR     input source vector.
*
*Parameters: input/output
* FIMEM   unknown vector.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
* SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPSYS(NGEFF),IPTRK,IPMACR
      INTEGER IMPX,NUN,K,NREG,M,NG,NGEFF,LC,PACA,NZON(NLONG),
     1 KEYFLX(NREG),KEYCUR(NLONG-NREG),NGIND(NGEFF),MAXI
      REAL EPSI,FIMEM(NUN,NGEFF),QFR(NUN,NG)
      LOGICAL LFORW,NCONV(NGEFF)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPSYS
      INTEGER, TARGET, SAVE, DIMENSION(1) :: IDUMMY
      REAL, TARGET, SAVE, DIMENSION(1) :: DUMMY
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) IM_PTR,MCU_PTR,JU_PTR,IM0_PTR,MCU0_PTR,IPERM_PTR,
     1 DIAGF_PTR,CF_PTR,CQ_PTR,LUDF_PTR,LUCF_PTR,DIAGQ_PTR
      INTEGER, POINTER, DIMENSION(:) :: IM,MCU,IPERM,JU,IM0,MCU0
      REAL, POINTER, DIMENSION(:) :: DIAGQ,CQ,LUDF,LUCF,CF,DIAGF
*----
*  INITIALIZE POINTERS
*----
      JU=>IDUMMY
      IM0=>IDUMMY
      MCU0=>IDUMMY
      LUDF=>DUMMY
      LUCF=>DUMMY
      CF=>DUMMY
      DIAGF=>DUMMY
*
      IF(K.LE.0) CALL XABORT('MCGFLS: INVALID VALUE OF K.')
      M1=M+1
*     recover connection matrices
      CALL LCMGPD(IPTRK,'IM$MCCG',IM_PTR)
      CALL LCMGPD(IPTRK,'MCU$MCCG',MCU_PTR)
      CALL C_F_POINTER(IM_PTR,IM,(/ NLONG+1 /))
      CALL C_F_POINTER(MCU_PTR,MCU,(/ LC /))
*     recover permutation array
      CALL LCMGPD(IPTRK,'PI$MCCG',IPERM_PTR)
      CALL C_F_POINTER(IPERM_PTR,IPERM,(/ NLONG /))
      IF(PACA.GE.2) THEN
         CALL LCMGPD(IPTRK,'JU$MCCG',JU_PTR)
         CALL C_F_POINTER(JU_PTR,JU,(/ NLONG /))
      ENDIF
      IF(PACA.EQ.3) THEN
         CALL LCMLEN(IPTRK,'IM0$MCCG',LIM0,ITYLCM)
         CALL LCMLEN(IPTRK,'MCU0$MCCG',LMCU0,ITYLCM)
         CALL LCMGPD(IPTRK,'IM0$MCCG',IM0_PTR)
         CALL LCMGPD(IPTRK,'MCU0$MCCG',MCU0_PTR)
         CALL C_F_POINTER(IM0_PTR,IM0,(/ LIM0 /))
         CALL C_F_POINTER(MCU0_PTR,MCU0,(/ LMCU0 /))
      ELSE
         LMCU0=0
      ENDIF
      DO II=1,NGEFF
      IF(NCONV(II)) THEN
         IG=NGIND(II)
         JPSYS=KPSYS(II)
         CALL LCMGPD(JPSYS,'DIAGF$MCCG',DIAGF_PTR)
         CALL LCMGPD(JPSYS,'CF$MCCG',CF_PTR)
         CALL C_F_POINTER(DIAGF_PTR,DIAGF,(/ NLONG /))
         CALL C_F_POINTER(CF_PTR,CF,(/ LC /))
         IF(PACA.GE.2) THEN
            CALL LCMGPD(JPSYS,'ILUDF$MCCG',LUDF_PTR)
            CALL C_F_POINTER(LUDF_PTR,LUDF,(/ NLONG /))
            IF(PACA.LT.4) THEN
               CALL LCMGPD(JPSYS,'ILUCF$MCCG',LUCF_PTR)
               CALL C_F_POINTER(LUCF_PTR,LUCF,(/ LC /))
            ENDIF
         ENDIF
         CALL LCMGPD(JPSYS,'DIAGQ$MCCG',DIAGQ_PTR)
         CALL LCMGPD(JPSYS,'CQ$MCCG',CQ_PTR)
         CALL C_F_POINTER(DIAGQ_PTR,DIAGQ,(/ NLONG /))
         CALL C_F_POINTER(CQ_PTR,CQ,(/ LC /))
         CALL MCGCDD(IMPX,IPMACR,IG,NG,NGEFF,NGIND,NCONV,M,NLONG,NUN,
     1        NREG,LC,LFORW,PACA,NZON,KEYFLX,KEYCUR,IPERM,IM,MCU,JU,
     2        EPSI,MAXI,FIMEM(1,II),QFR,DIAGQ,CQ,DIAGF,CF,LUDF,LUCF,
     3        LMCU0,IM0,MCU0)
      ENDIF
      ENDDO
*     
      RETURN
      END
