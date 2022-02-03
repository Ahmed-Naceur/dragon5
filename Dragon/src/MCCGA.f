*DECK MCCGA
      SUBROUTINE MCCGA(IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,NBMIX,NANI,
     1 NALBP,ISTRM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of PJJ for flux integration when isotropic scattering is
* considered and calculation of preconditioning matrices for 
* Algebraic Collapsing Acceleration or Self-Collision Probability 
* acceleration of inner iterations (vectorial version).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and R. Le Tellier
*
*Parameters: input/output
* IPSYS   pointer to the PIJ LCM object (L_PIJ signature). IPSYS is a
*         list of NGRP directories.
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  tracking file unit number.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* NBMIX   number of mixtures.
* NANI    number of Legendre orders.
* NALBP   number of physical albedos.
* ISTRM   type of streaming effect:
*         =1 no streaming effect;
*         =2 isotropic streaming effect;
*         =3 anisotropic streaming effect.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK
      INTEGER IFTRAK,IMPX,NGRP,NBMIX,NANI,NALBP,ISTRM,NPSYS(NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER JPAR(NSTATE),TRTY,PACA,STIS,IGB(8)
      CHARACTER*4 TEXT4
      REAL ZREAL(4),DELU,FACSYM
      LOGICAL LEXA,LEXF,CYCLIC,LTMT,LACA,LPJJ,LPJJAN,LVOID,LPRISM,
     1 LBIHET
      TYPE(C_PTR) JPSYS
      EXTERNAL MCGDSCA,MCGDSCE,MCGDDDF,MCGDS2,MOCDS2,MCGDSP,MOCDSP,
     1 MCGDS2A,MCGDS2E
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NGIND
      REAL, ALLOCATABLE, DIMENSION(:) :: CPO,SIGAL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CAZ0,CAZ1,CAZ2
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: KPSYS
*----
*  RECOVER MCCG3D SPECIFIC PARAMETERS
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      IF(JPAR(4).GT.NBMIX) CALL XABORT('MCCGA: INVALID NBMIX.')
      IF(IFTRAK.LE.0) CALL XABORT('MCCGA: INVALID TRACKING FILE.')
*     recover state-vector information
      LBIHET=JPAR(40).NE.0
      IF(LBIHET) THEN
         CALL LCMSIX(IPTRK,'BIHET',1)
         CALL LCMGET(IPTRK,'PARAM',IGB)
         NREG=IGB(3)
         CALL LCMSIX(IPTRK,' ',2)
      ELSE
         NREG=JPAR(1)
      ENDIF
      NFI=NREG+JPAR(5)
      IF(JPAR(6).NE.NANI) CALL XABORT('MCCGA: INVALID NANI.')
      TRTY=JPAR(9)
      IF(TRTY.EQ.1) THEN
         IF(JPAR(5).EQ.0) NFI=NREG+1
         CYCLIC=.TRUE.
         NLONG=NREG
      ELSE
         CYCLIC=.FALSE.
         NLONG=NFI
      ENDIF
      NZP=JPAR(39)
      LPRISM=(NZP.NE.0)
      CALL LCMGET(IPTRK,'MCCG-STATE',JPAR)
      NMU=JPAR(2)
      NMAX=JPAR(5)
      IAAC=JPAR(7)
      STIS=JPAR(15)
      ISCR=JPAR(8)
      LC=JPAR(6)
      LPS=JPAR(9)
      PACA=JPAR(10)
      LC0=JPAR(17)
      LTMT=(JPAR(14).EQ.1)
      LEXA=(JPAR(11).EQ.1)
      LEXF=(JPAR(12).EQ.1)
      NPJJM=JPAR(16)
*     recover real parameters
      CALL LCMGET(IPTRK,'REAL-PARAM',ZREAL)
      HDD=ZREAL(2)
      DELU=ZREAL(3)
      FACSYM=ZREAL(4)
*     recover tracking file information
      REWIND IFTRAK
      READ(IFTRAK) TEXT4,NCOMNT,NBTR,IFMT
      DO ICOM=1,NCOMNT
         READ(IFTRAK)
      ENDDO
      READ(IFTRAK) NDIM,ISPEC,N2REG,N2SOU,NALBG,NCOR,NANGL,MXSUB,MXSEG
      IF(NCOR.NE.1) 
     1 CALL XABORT('MCCGA: INVALID TRACKING FILE: NCOR.NE.1')
      READ(IFTRAK)
      READ(IFTRAK)
      READ(IFTRAK)
      READ(IFTRAK)
      ALLOCATE(CAZ0(NANGL),CAZ1(NANGL),CAZ2(NANGL),CPO(NMU))
      IF(NDIM.EQ.2) THEN
         CALL LCMGET(IPTRK,'XMU$MCCG',CPO)
         READ(IFTRAK) (CAZ1(JJ),CAZ2(JJ),JJ=1,NANGL)
      ELSE ! NDIM.EQ.3
**        correction Sylvie Musongela, december 2019
         READ(IFTRAK) (CAZ1(JJ),CAZ2(JJ),CAZ0(JJ),JJ=1,NANGL)
         DO JJ=1,NANGL
            CAZ1(JJ)=CAZ1(JJ)/SQRT(1.0D0-CAZ0(JJ)*CAZ0(JJ))
            CAZ2(JJ)=CAZ2(JJ)/SQRT(1.0D0-CAZ0(JJ)*CAZ0(JJ))
         ENDDO
      ENDIF
*---
* DETERMINE THE NUMBER OF GROUPS TO BE PROCESSED
* RECOVER POINTERS TO EACH GROUP PROPERTIES
* CREATE AN INDEX FOR THE GROUPS TO BE PROCESSED
*---
      NGEFF=0
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) NGEFF=NGEFF+1
      ENDDO
      ALLOCATE(NGIND(NGEFF),KPSYS(NGEFF))
      II=1
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) THEN
            NGIND(II)=IG
            IF(LBIHET) THEN
               JPSYS=LCMGIL(IPSYS,IOFSET)
               KPSYS(II)=LCMGID(JPSYS,'BIHET')
            ELSE
               KPSYS(II)=LCMGIL(IPSYS,IOFSET)
            ENDIF
            II=II+1
         ENDIF
      ENDDO
*----
*  CONSTRUCT TOTAL CROSS SECTIONS ARRAY AND CHECK FOR ZERO CROSS SECTION
*----
      ALLOCATE(SIGAL((NBMIX+7)*NGEFF))
      CALL MCGSIG(IPTRK,NBMIX,NGEFF,NALBP,KPSYS,SIGAL,LVOID)
      IF((LVOID).AND.(STIS.EQ.-1)) THEN
         IF(IMPX.GT.0) 
     1       WRITE(6,*) 'VOID EXISTS -> STIS SET TO 1 INSTEAD OF -1'
         STIS=1
      ENDIF
*---
* IS THERE SOMETHING TO DO ?
*---
      LACA=(IAAC.GT.0)
      LPJJ=((STIS.EQ.1).OR.(ISCR.GT.0))
      IF(.NOT.(LACA.OR.LPJJ)) GOTO 10
      LPJJAN=(LPJJ.AND.(NANI.GT.1))
      IF(HDD.GT.0.0) THEN
         ISCH=0
      ELSEIF(LEXF) THEN
         ISCH=-1
      ELSE
         ISCH=1
      ENDIF
*----
*  PRECONDITIONING MATRICES CALCULATION
*----
      IF(ISCH.EQ.1) THEN
*     PJJ/SCR: Step-Characteristics Scheme with Tabulated Exponentials
         IF(CYCLIC) THEN
*        ACA: cyclic tracking
            IF(LEXA) THEN
*           ACA: Exact Exponentials
               CALL MCGASM(MCGDSCA,MOCDS2,MOCDSP,MCGDS2E,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ELSE
*           ACA: Tabulated Exponentials
               CALL MCGASM(MCGDSCA,MOCDS2,MOCDSP,MCGDS2A,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ENDIF
         ELSE
*        ACA: non-cyclic tracking
            IF(LEXA) THEN
*           ACA: Exact Exponentials
               CALL MCGASM(MCGDSCA,MCGDS2,MCGDSP,MCGDS2E,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ELSE
*           ACA: Tabulated Exponentials
               CALL MCGASM(MCGDSCA,MCGDS2,MCGDSP,MCGDS2A,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ENDIF
         ENDIF        
      ELSEIF(ISCH.EQ.0) THEN
*     PJJ/SCR: Diamond-Differencing Scheme
         IF(CYCLIC) THEN
*        ACA: cyclic tracking
            IF(LEXA) THEN
*           ACA: Exact Exponentials
               CALL MCGASM(MCGDDDF,MOCDS2,MOCDSP,MCGDS2E,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ELSE
*           ACA: Tabulated Exponentials
               CALL MCGASM(MCGDDDF,MOCDS2,MOCDSP,MCGDS2A,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ENDIF
         ELSE
*        ACA: non-cyclic tracking
            IF(LEXA) THEN
*           ACA: Exact Exponentials
               CALL MCGASM(MCGDDDF,MCGDS2,MCGDSP,MCGDS2E,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ELSE
*           ACA: Tabulated Exponentials
               CALL MCGASM(MCGDDDF,MCGDS2,MCGDSP,MCGDS2A,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ENDIF
         ENDIF
      ELSEIF(ISCH.EQ.-1) THEN
*     PJJ/SCR: Step-Characteristics Scheme with Exact Exponentials
         IF(CYCLIC) THEN
*        ACA: cyclic tracking
            IF(LEXA) THEN
*           ACA: Exact Exponentials
               CALL MCGASM(MCGDSCE,MOCDS2,MOCDSP,MCGDS2E,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ELSE
*           ACA: Tabulated Exponentials
               CALL MCGASM(MCGDSCE,MOCDS2,MOCDSP,MCGDS2A,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ENDIF
         ELSE
*        ACA: non-cyclic tracking
            IF(LEXA) THEN
*           ACA: Exact Exponentials
               CALL MCGASM(MCGDSCE,MCGDS2,MCGDSP,MCGDS2E,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ELSE
*           ACA: Tabulated Exponentials
               CALL MCGASM(MCGDSCE,MCGDS2,MCGDSP,MCGDS2A,IPTRK,
     1              KPSYS,IMPX,IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,
     2              NBMIX,NMU,NANGL,NMAX,LC,NDIM,NGIND,CYCLIC,ISCR,
     3              CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,
     4              LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,
     5              FACSYM,ISTRM)
            ENDIF
         ENDIF
      ENDIF
*
   10 DEALLOCATE(SIGAL,KPSYS,NGIND,CPO,CAZ2,CAZ1,CAZ0)
      RETURN
      END
