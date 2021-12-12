*DECK PIJXL3
      SUBROUTINE PIJXL3(IPTRK,IPRT,NGRP,NANI,NBMIX,NPSYS,NRENOR,LEAKSW,
     > XSSIGT,XSSIGW,NELPIJ,PIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the collision probabilities in EXCELL without producing
* a tracking file. Based on subroutine XL3TRK in DRAGON 3.4.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and A. Hebert
*
*Parameters: input
* IPTRK   pointer to the tracking (L_TRACK signature).
* IPRT    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* NANI    number of Legendre orders (usually equal to one).
* NBMIX   number of mixtures.
* NPSYS   index set to zero if a group is not to be processed. Usually,
*         NPSYS(I)=I.
* NRENOR  normalization scheme for PIJ matrices.
* LEAKSW  leakage flag (=.true. if neutron leakage through external
*         boundary is present).
* XSSIGT  total macroscopic cross sections ordered by mixture.
* XSSIGW  P0 within-group scattering macroscopic cross sections
*         ordered by mixture.
* NELPIJ  number of elements in symmetrized pij matrix.
*
*Parameters: output
* PIJ     reduced and symmetrized collision probabilities.
*
*-----------------------------------------------------------------------
*
      USE             GANLIB
      IMPLICIT        NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL         LEAKSW
      TYPE(C_PTR)     IPTRK
      INTEGER         IPRT,NGRP,NANI,NBMIX,NPSYS(NGRP),NRENOR,NELPIJ
      REAL            XSSIGT(0:NBMIX,NGRP),XSSIGW(0:NBMIX,NANI,NGRP),
     >                PIJ(NELPIJ,NGRP)
*----
*  LOCAL VARIABLES
*----
      INTEGER         IOUT,NALB,NSTATE,ICPALL,ICPEND
      PARAMETER      (IOUT=6,NALB=6,NSTATE=40,ICPALL=4,ICPEND=3)
      INTEGER         NDIM  ,ISPEC ,NANGLE,NANGL ,ISYMM, NORE
      INTEGER         NALBG ,NC    ,NTR   ,NTZ   ,ICL   ,NSOUT ,
     >                ITGEO ,NRMV  ,NTY   ,LTRK  ,LINMAX,NUNK  ,
     >                MAXR  ,NEXTGE,NTX   ,NCOR  ,NSUR  ,NTOTCL,
     >                NVOL  ,NV    ,NS    ,IGRP  ,ISOUT ,ILONG ,
     >                ITYPE ,INDPIJ,IIN   ,IBM   ,I     ,J     ,
     >                NPIJ  ,NREG  ,NUNKMR
      INTEGER         ISTATE(NSTATE),LCLSYM(3)
      INTEGER         MXANGL,ICODE(NALB)
      LOGICAL         SWVOID,SWNZBC,LSKIP
      REAL            ALBOLD(NALB),EXTKOP(NSTATE),DENUSR, RCUTOF,
     >                CUTOFX, FACT
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MINDIM,MAXDIM,ICORD,INDEL,
     > KEYMRG,MATALB,MATMRG,ICUR,INCR,NUMERO,MATRT
      REAL, ALLOCATABLE, DIMENSION(:) :: REMESH,VOLSUR,VOLMRG,CONV,
     > TRKBEG,TRKDIR,FFACT,ANGLES,DENSTY
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGVOL,SIGTAL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LENGHT,VOLTRK,
     > PSST,PSVT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DBLPIJ,PCSCT
*----
*  INTRINSIC FUNCTION FOR POSITION IN CONDENSE PIJ MATRIX
*----
      INTEGER INDPOS
      INDPOS(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
*----
*  READ THE GEOMETRY INFORMATION STORED ON IPTRK
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NORE=ISTATE(8)
      LTRK=ISTATE(9)+1
      NANGLE=ISTATE(11)
      ISYMM=ISTATE(12)
      CALL LCMGET(IPTRK,'EXCELTRACKOP',EXTKOP)
      CUTOFX=EXTKOP(1)
      DENUSR=EXTKOP(2)
      RCUTOF=EXTKOP(3)
      CALL LCMSIX(IPTRK,'EXCELL',1)
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NDIM     =ISTATE(1)
      NSUR    =-ISTATE(2)
      NVOL     =ISTATE(3)
      NTOTCL   =ISTATE(4)
      MAXR     =ISTATE(5)
      NUNK     =ISTATE(6)
      NEXTGE   =ISTATE(7)
*----
*  Intrinsic symmetries used in geometry
*  Use these to simplify tracking unless 
*  NOSYMM tracking option activated
*----
      LCLSYM(1) =ISTATE(8)
      LCLSYM(2) =ISTATE(9)
      LCLSYM(3) =ISTATE(10)
      IF(ISYMM .NE. 0) THEN
        ISYMM=0
        IF(NDIM .EQ. 2) THEN
          IF(LCLSYM(1) .NE. 0) THEN
*----
*  X SYMMETRY
*----
            ISYMM=2
          ENDIF
          IF(LCLSYM(2) .NE. 0) THEN
            IF(ISYMM .EQ. 0) THEN
*----
*  Y SYMMETRY
*----
              ISYMM=4
            ELSE
*----
*  X AND Y SYMMETRY
*----
              ISYMM=8 
            ENDIF
          ENDIF
*----
*  X-Y DIAGONAL SYMMETRY
*---- 
        ELSE
          IF(LCLSYM(1) .NE. 0) THEN
*----
*  X SYMMETRY
*----
            ISYMM=2
          ENDIF
          IF(LCLSYM(2) .NE. 0) THEN
            IF(ISYMM .EQ. 0) THEN
*----
*  Y SYMMETRY
*----
              ISYMM=4
            ELSE
*----
*  X AND Y SYMMETRY
*----
              ISYMM=8 
            ENDIF
          ENDIF 
          IF(LCLSYM(3) .NE. 0) THEN
*----
*  Z SYMMETRY
*----
            ISYMM=ISYMM+16
          ENDIF
        ENDIF
        IF(ISYMM .EQ. 0) ISYMM=1
      ENDIF
      ALLOCATE(MINDIM(NTOTCL),MAXDIM(NTOTCL),ICORD(NTOTCL),
     > INDEL(4*NUNK),KEYMRG(NUNK),MATALB(NUNK))
      ALLOCATE(REMESH(MAXR),VOLSUR(NUNK))
      CALL LCMGET(IPTRK,'MINDIM      ',MINDIM)
      CALL LCMGET(IPTRK,'MAXDIM      ',MAXDIM)
      CALL LCMGET(IPTRK,'ICORD       ',ICORD )
      CALL LCMGET(IPTRK,'INDEX       ',INDEL )
      CALL LCMGET(IPTRK,'REMESH      ',REMESH)
      CALL LCMGET(IPTRK,'KEYMRG      ',KEYMRG)
      CALL LCMGET(IPTRK,'MATALB      ',MATALB)
      CALL LCMGET(IPTRK,'VOLSUR      ',VOLSUR)
      CALL LCMSIX(IPTRK,'EXCELL      ',2)
      CALL LCMGET(IPTRK,'ICODE       ',ICODE )
      CALL LCMGET(IPTRK,'ALBEDO      ',ALBOLD)
*----
*  VERIFY SYMMETRY AND STUDY TRACKING PARAMETERS. ARE THEY BASICALLY
*  POSSIBLE ?
*----
      MXANGL=0
      IF(LTRK .EQ. 1)THEN
        NCOR= 1
        IF(NDIM .EQ. 2) THEN
          MXANGL=NANGLE
          IF(ISYMM .GE. 2) THEN
            NANGL = (NANGLE+1)/2
          ELSE 
            NANGL = NANGLE
          ENDIF
          IF( RCUTOF.GT.0.0 ) NCOR= 2
        ELSE IF(NDIM .EQ. 3) THEN
          IF(MOD(NANGLE,2) .EQ. 1)THEN
            NANGLE=NANGLE+1
            WRITE(IOUT,'(/31H MESS = ONLY EVEN # EQN ANGLES )')
          ENDIF
          IF(NANGLE .GT. 16)THEN
            NANGLE=16
            WRITE(IOUT,'(/31H MESS = 16 IS MAX # EQN ANGLES )')
          ENDIF
          MXANGL=(NANGLE * (NANGLE+2)) / 2 
          IF(NEXTGE .EQ. 1) THEN
            NANGL = (NANGLE * (NANGLE+2)) / 8
          ELSE
            IF(ISYMM .EQ. 8 .OR. ISYMM .EQ. 24) THEN
              NANGL = (NANGLE * (NANGLE+2)) / 8 
            ELSE IF(ISYMM .EQ. 2  .OR. ISYMM .EQ. 4 .OR. 
     >              ISYMM .EQ. 18 .OR. ISYMM .EQ. 20 ) THEN
              NANGL = (NANGLE * (NANGLE+2)) / 4 
            ELSE
              NANGL = (NANGLE * (NANGLE+2)) / 2 
            ENDIF
          ENDIF
          IF(RCUTOF .GT. 0.0) NCOR= 4
        ENDIF
      ELSEIF( LTRK.EQ.2 )THEN
        NCOR  = 1
        MXANGL=NANGLE
        IF( NDIM.EQ.2 )THEN
           NANGL = NANGLE
        ELSEIF( NDIM.EQ.3 )THEN
           CALL XABORT('PIJXL3: *TSPC* NOT AVAILABLE FOR 3-D GEOMETRY')
        ENDIF
        CUTOFX= RCUTOF
      ENDIF
      IF( IPRT.GE.1 ) THEN
         WRITE(IOUT,6002) NANGL,ISYMM,CUTOFX,DENUSR,RCUTOF
      ENDIF
      IF( IPRT.GT.1 .AND. NEXTGE.EQ.0 )THEN
*
*        IF PRINT REQUIRED AND OVERALL CARTESIAN GEOMETRY
*        PRINT CARTESIAN REGION MAP
         NTX= MAXDIM(1)-MINDIM(1)
         NTY= MAXDIM(2)-MINDIM(2)
         NTZ= MAXDIM(3)-MINDIM(3)
         NTR=0
         DO 103 ICL=4,NTOTCL
            NTR= MAX(NTR,MAXDIM(ICL)-MINDIM(ICL)+1)
  103    CONTINUE
         CALL XELGPR(NDIM,NTX,NTY,NTZ,NTR,ISYMM,
     >               NSUR,NVOL,NTOTCL,MINDIM,MAXDIM,
     >               KEYMRG,INDEL,MATALB)
      ENDIF
      ALLOCATE(VOLTRK((NANGL+1)*NUNK))
*
      NV=  NVOL
      NS= -NSUR
      ALLOCATE(VOLMRG(NUNK),MATMRG(NUNK))
      ITGEO=3
      CALL XELCMP(   NS,    NV, VOLSUR, MATALB, KEYMRG,
     >            NSOUT,  NREG, VOLMRG, MATMRG, ITGEO,ICODE)
      NUNKMR= NREG+NSOUT+1
      NPIJ= (NUNKMR*(NUNKMR+1))/2
      IF( IPRT .GT. 1 ) WRITE(IOUT,6000) (NGRP*NPIJ/128)
      ALLOCATE(DBLPIJ(NPIJ,NGRP))
      IF( IPRT .GT. 1 ) WRITE(IOUT,6001)
*
*     ALLOCATE AND CHARGE TOTAL XS PER REGION
      ALLOCATE(SIGTAL(NUNKMR,NGRP),SIGVOL(NREG,NGRP))
*
*  3) DO THE TRACKING OF THE EXACT GEOMETRY FOR *NEWT* OPTION.
      IF( LTRK.NE.0 )THEN
         NC= NTOTCL - 3
         IF( IPRT.GE.1 )THEN
            WRITE(IOUT,'(1H )')
            IF( NC.EQ.0 )THEN
               WRITE(IOUT,'(/38H NOW, TRACKING GEOMETRY WITH NO CYLIND,
     >         2HER/)')
            ELSEIF( NC.EQ.1 )THEN
               WRITE(IOUT,'(/38H NOW, TRACKING GEOMETRY WITH ONE CYLIN,
     >         3HDER/)')
            ELSE
               WRITE(IOUT,'(/28H NOW, TRACKING GEOMETRY WITH,I4,
     >         10H CYLINDERS/)') NC
            ENDIF
         ENDIF
         ALLOCATE(ICUR(NTOTCL),INCR(NTOTCL))
         ALLOCATE(CONV(NTOTCL),TRKBEG(NTOTCL),TRKDIR(NTOTCL))
*
*  3.0)  WRITE FIRST RECORDS OF THE UNNORMALIZED TRACKING FILE
         IF( LTRK.EQ.1 )THEN
            LINMAX= 2*NVOL + 10
         ELSE
            LINMAX= 8*NANGL*(2*NVOL + 8)
         ENDIF
         ISPEC = LTRK-1
         NALBG = 6
         ALLOCATE(NUMERO(LINMAX))
         ALLOCATE(LENGHT(LINMAX),ANGLES(3*MXANGL),DENSTY(MXANGL))
*
         NRMV=1
         CALL XL3TI3( IPRT,   NANGLE, DENUSR, ISYMM,  ANGLES, DENSTY,
     >                NTOTCL, NEXTGE, MAXR,   REMESH, LINMAX, RCUTOF,
     >                NSUR,   NVOL,   INDEL,  MINDIM, MAXDIM, ICORD,
     >                INCR,   ICUR,   TRKBEG, CONV,   TRKDIR, LENGHT, 
     >                NUMERO, NPIJ,   NGRP,   SIGTAL, SWVOID, NORE,
     >                NRMV,   VOLTRK, KEYMRG,-NSOUT,  NREG,   NPSYS,
     >                DBLPIJ )
*
         CALL XL3NTR( IPRT, NDIM, ISPEC, NS, NV, NORE, VOLSUR, KEYMRG,
     >                MATALB, NANGL, VOLTRK, DENSTY )
*
         CALL XL3SIG( NGRP, NBMIX,  XSSIGT, ALBOLD, NPSYS, NGRP, -NSOUT,
     >                NREG, MATMRG, VOLMRG(NSOUT+2), SIGTAL, SIGVOL,
     >                SWVOID, SWNZBC)
*
         NRMV=0
         CALL XL3TI3( IPRT,   NANGLE, DENUSR, ISYMM,  ANGLES, DENSTY,
     >                NTOTCL, NEXTGE, MAXR,   REMESH, LINMAX, RCUTOF,
     >                NSUR,   NVOL,   INDEL,  MINDIM, MAXDIM, ICORD ,
     >                INCR,   ICUR,   TRKBEG, CONV,   TRKDIR, LENGHT, 
     >                NUMERO, NPIJ,   NGRP,   SIGTAL, SWVOID, NORE,
     >                NRMV,   VOLTRK, KEYMRG, -NSOUT, NREG  , NPSYS,
     >                DBLPIJ )
*
         CALL QIJCMP(NREG,-NSOUT,NPIJ,NGRP,NCOR,VOLMRG,SIGTAL,DBLPIJ,
     >               NPSYS)
*----
*  RENORMALIZE ALL ISOTROPIC PROBS WITH VARIOUS OPTIONS
*----
         DO 2060 IGRP=1,NGRP
         IF(NPSYS(IGRP).EQ.0) GO TO 2060
         IF( NRENOR.EQ.1 )THEN
*
*           NORMALIZATION USING GELBARD SCHEME
            CALL PIJRGL(IPRT,NREG,NSOUT,SIGTAL(1,IGRP),DBLPIJ(1,IGRP))
         ELSEIF( NRENOR.EQ.2 )THEN
*
*           NORMALIZATION WORKING ON DIAGONAL COEFFICIENTS
            CALL PIJRDG(NREG,NSOUT,SIGTAL(1,IGRP),DBLPIJ(1,IGRP) )
         ELSEIF( NRENOR.EQ.3 )THEN
*
*           NORMALIZATION WORKING ON WEIGHT FACTORS TO KEEP DIAG = 0.0
            CALL PIJRNL(IPRT,NREG,NSOUT,SIGTAL(1,IGRP),DBLPIJ(1,IGRP))
         ELSEIF( NRENOR .EQ. 4 )THEN  ! ATTENTION
*
*           NORMALIZATION WORKING ON WEIGHT FACTORS ADDITIVE (HELIOS)
            CALL PIJRHL(IPRT,NREG,NSOUT,SIGTAL(1,IGRP),DBLPIJ(1,IGRP))
         ENDIF
         IF( IPRT.GE.ICPALL )THEN
            WRITE(IOUT,'(1H )')
            WRITE(IOUT,'(35H   COLLISION PROBABILITIES OUTPUT: ,
     >                   35H *BEFORE* ALBEDO REDUCTION          )')
            CALL PIJWPR(0,NREG,NSOUT,SIGTAL(1,IGRP),DBLPIJ(1,IGRP),
     >      SIGVOL(1,IGRP),1)
         ENDIF
 2060    CONTINUE
*----
*  ELIMINATION OF SURFACES FOR PIJ
*----
         IF( SWNZBC )THEN
            ALLOCATE(PSST(NSOUT*NSOUT),PSVT(NSOUT*NREG))
            ALLOCATE(MATRT(NSOUT))
            CALL LCMLEN(IPTRK,'BC-REFL+TRAN',ILONG,ITYPE)
            IF(ILONG.EQ.NSOUT) THEN
              CALL LCMGET(IPTRK,'BC-REFL+TRAN',MATRT)
            ELSE
               DO 130 ISOUT=1,NSOUT
                 MATRT(ISOUT)=ISOUT
 130           CONTINUE
            ENDIF
            DO 2080 IGRP=1,NGRP
              IF(NPSYS(IGRP).EQ.0) GO TO 2080
              CALL PIJABC(NREG,NSOUT,NPIJ,SIGTAL(1,IGRP),MATRT,
     >                    DBLPIJ(1,IGRP),PSST,PSVT)
 2080       CONTINUE
*
            DEALLOCATE(MATRT)
            DEALLOCATE(PSVT,PSST)
         ENDIF
*
         ALLOCATE(FFACT(NREG))
         DO 2090 IGRP=1,NGRP
         IF(NPSYS(IGRP).EQ.0) GO TO 2090
         IF( IPRT.GE.ICPEND )THEN
            WRITE(IOUT,'(1H )')
            WRITE(IOUT,'(35H   COLLISION PROBABILITIES OUTPUT: ,
     >                   35H *AFTER* ALBEDO REDUCTION          )')
            CALL PIJWPR(1,NREG,NSOUT,SIGTAL(1,IGRP),DBLPIJ(1,IGRP),
     >      SIGVOL(1,IGRP),1)
         ENDIF
*----
*  CHARGE PIJ MATRIX IN THE DRAGON SYMMETRIZED FORMAT
*----
         DO 160 IIN=1,NREG
            IF(SIGTAL(NSOUT+IIN+1,IGRP).EQ.0.0) THEN
               FFACT(IIN)=1.0
            ELSE
               FFACT(IIN)=1.0/SIGTAL(NSOUT+IIN+1,IGRP)
            ENDIF
  160    CONTINUE
         CALL PIJD2R(NREG,NSOUT,DBLPIJ(1,IGRP),FFACT,.FALSE.,NELPIJ,
     >               NPIJ,PIJ(1,IGRP))
 2090    CONTINUE
         DEALLOCATE(FFACT)
*
         DEALLOCATE(DENSTY,ANGLES,LENGHT,NUMERO,TRKDIR,TRKBEG,CONV,
     >   INCR,ICUR)
      ENDIF
      DEALLOCATE(INDEL,ICORD,MAXDIM,MINDIM,REMESH,DBLPIJ,SIGTAL,SIGVOL,
     > VOLSUR,VOLTRK,KEYMRG,MATALB)
*----
*  CHECK IF SCATTERING REDUCTION IS REQUIRED
*----
      ALLOCATE(PCSCT(NREG,2*NREG))
      DO 3000 IGRP=1,NGRP
      IF(NPSYS(IGRP).EQ.0) GO TO 3000
      LSKIP=.TRUE.
      DO 200 IBM=1,NBMIX
        LSKIP=LSKIP.AND.(XSSIGW(IBM,1,IGRP).EQ.0.0)
  200 CONTINUE
*----
*  COMPUTE THE SCATTERING-REDUCED CP MATRICES
*----
      IF(.NOT.LSKIP) THEN
        CALL PIJSMD(IPRT,NBMIX,NREG,MATMRG(NSOUT+2),VOLMRG(NSOUT+2),
     >              XSSIGW(0,1,IGRP),XSSIGT(0,IGRP),LEAKSW,PIJ(1,IGRP),
     >              PCSCT,1)
        DO 220 I=1,NREG
          FACT=VOLMRG(NSOUT+I+1)
          DO 210 J=1,NREG
            INDPIJ=INDPOS(I,J)
            PIJ(INDPIJ,IGRP)=REAL(PCSCT(I,J))*FACT
  210     CONTINUE
  220   CONTINUE
      ENDIF
 3000 CONTINUE
      DEALLOCATE(PCSCT,VOLMRG,MATMRG)
      RETURN 
*----
*  FORMAT
*----
 6000 FORMAT(' *** SPACE REQUIRED FOR CP MATRICES = ',I10,' K ***')
 6001 FORMAT(' *** CP MATRICES ALLOCATED            ',10X,'   ***')
 6002 FORMAT(
     > ' -----------------'/' RECOMPUTED PARAMETERS '/
     > ' NANGL  =',I10  ,' (NUMBER OF TRACKING ANGLES)'/
     > ' ISYMM  =',I10  ,' (TRACKING SYMMETRY FACTOR)'/
     > ' CUTOFX =',F10.5,' (CUTOFF FOR TRACK LENGTH)'/
     > ' DENS   =',F10.5,' (TRACK DENSITY)'/
     > ' PCORN  =',F10.5,' (CORNER DUPLICATION DISTANCE)'/
     > ' -----------------'/)
      END
