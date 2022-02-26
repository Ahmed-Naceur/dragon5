*DECK SHIDRV
      SUBROUTINE SHIDRV (IPLIB,IPTRK,IFTRAK,LEVEL,NGRO,NBISO,NBMIX,
     1 NREG,NUN,CDOOR,NRES,IMPX,ISONRF,ISONAM,MIX,DEN,SN,SB,LSHI,
     2 IPHASE,IPROB,MAT,VOL,KEYFLX,LEAKSW,TITR,IGRMIN,IGRMAX,MAXX0,
     3 IBIEFF,IGC,ITRANZ,EPS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform a multidimensional self-shielding calculation in order to
* compute the dilution cross section of each resonant isotope present
* in the domain.
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
* IPLIB   pointer to the internal microscopic cross section library
*         (L_LIBRARY signature).
* IPTRK   pointer to the tracking. (L_TRACK signature).
* IFTRAK  unit number of the sequential binary tracking file.
* LEVEL   type of self-shielding model (=0 original Stamm'ler model;
*         =1 original Stamm'ler model with Nordheim approximation;
*         =2 Stamm'ler model with Nordheim approximation and Riemann
*         integration method).
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NBMIX   number of mixtures in the macrolib.
* NREG    number of regions.
* NUN     number of unknowns in the flux or source vector in one
*         energy group.
* CDOOR   name of the geometry/solution module.
* NRES    number of resonant mixtures.
* IMPX    print flag.
* ISONRF  reference name of isotopes.
* ISONAM  alias name of isotopes.
* MIX     mix number of each isotope (can be zero).
* DEN     density of each isotope.
* LSHI    resonant region number associated with each isotope.
*         Infinite dilution will be assumed if LSHI(i)=0.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* IPROB   adjoint macrolib flag (=0 direct; =1 adjoint).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* LEAKSW  leakage flag (=.TRUE. only if leakage is present on the outer
*         surface).
* TITR    title.
* IGRMIN  first group where the self-shielding is applied.
* IGRMAX  most thermal group where the self-shielding is applied.
* MAXX0   maximum number of self-shielding iterations.
* IBIEFF  Livolant-Jeanpierre normalization flag (=1 to activate).
* IGC     Goldstein-Cohen approximation flag (=1 to activate).
*         The Goldstein-Cohen approximation is activated only in cases
*         where they are available in the internal library.
* ITRANZ  type of transport correction used in the self-shielding
*         calculations.
* EPS     convergence criterion for the self-shielding iterations.
*
*Parameters: input/output
* SN      on input, estimate of the dilution cross section in each 
*         energy group of each isotope. A value of 1.0e10 is used 
*         for infinite dilution and;
*         on output, computed dilution cross section in each energy 
*         group of each isotope.
*
*Parameters: output
* SB      dilution cross section as used in Livolant-Jeanpierre
*         normalization.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPTRK
      INTEGER IFTRAK,LEVEL,NGRO,NBISO,NBMIX,NREG,NUN,NRES,IMPX,
     1 ISONRF(3,NBISO),ISONAM(3,NBISO),MIX(NBISO),LSHI(NBISO),
     2 IPHASE,IPROB,MAT(NREG),KEYFLX(NREG),IGRMIN,IGRMAX,MAXX0,IBIEFF,
     3 IGC,ITRANZ
      REAL DEN(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO),VOL(NREG),EPS
      LOGICAL LEAKSW
      CHARACTER CDOOR*12,TITR*72
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NALPHA=9,NRAT=(NALPHA+1)/2,NSTATE=40)
      TYPE(C_PTR) JPLIB,KPLIB
      INTEGER IPAR(NSTATE)
      REAL TMPDAY(3)
      CHARACTER HSMG*130
      LOGICAL START,BIEFF,LGC,LOGDO
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LSHI2
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGE
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGT1,SIGT2,SIGT3
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK,MASKL
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: NOCONV
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SIGT1(NBMIX,NGRO),SIGT2(NBMIX,NGRO),SIGT3(NBMIX,NGRO))
      ALLOCATE(MASK(NBMIX),MASKL(NGRO),NOCONV(NBMIX,NGRO))
*
      IF(IMPX.GE.5) THEN
         WRITE (6,'(//23H SHIDRV: VALUES OF MAT:)')
         I1=1
         KI=(NREG-1)/11+1
         DO 10 I=1,KI
         I2=I1+10
         IF(I2.GT.NREG) I2=NREG
         WRITE (6,340) (J,J=I1,I2)
         WRITE (6,350) (MAT(J),J=I1,I2)
         I1=I1+11
   10    CONTINUE
         WRITE (6,'(//)')
      ENDIF
*----
*  RECOVER SELF SHIELDING DATA
*----
      IF(LEAKSW) CALL XABORT('SHIDRV: NEUTRON LEAKAGE IS FORBIDDEN.')
      IF(CDOOR.EQ.' ') CALL XABORT('SHIDRV: THE GEOMETRY IS NOT YET '
     1 //'DEFINED.')
      BIEFF=(IBIEFF.EQ.1)
      LGC=(IGC.EQ.1)
*
      IF(IMPX.GT.0) THEN
         WRITE (6,400) TITR,CDOOR
         WRITE (6,'(25H STAMM''LER APPROXIMATION./)')
         WRITE (6,405) IGRMIN,IGRMAX,MAXX0,IBIEFF,IGC,ITRANZ,LEVEL,
     1   IPHASE,EPS
      ENDIF
      IF(NRES.EQ.0) THEN
         WRITE (6,410)
         RETURN
      ENDIF
      DO 30 I=1,NREG
      IF(MAT(I).GT.NBMIX) THEN
         WRITE (HSMG,380) NBMIX
         CALL XABORT(HSMG)
      ENDIF
   30 CONTINUE
      IGRMAX=MIN(IGRMAX,NGRO)
      DO 55 LLL=1,NGRO
      DO 40 I=1,NBISO
      SB(LLL,I)=SN(LLL,I)
   40 CONTINUE
      DO 50 IBM=1,NBMIX
      NOCONV(IBM,LLL)=.FALSE.
   50 CONTINUE
   55 CONTINUE
      CALL LCMSIX(IPLIB,'MACROLIB',1)
      JPLIB=LCMGID(IPLIB,'GROUP')
      DO 70 LLL=IGRMIN,IGRMAX
      IF(IPROB.EQ.0) LL=LLL
      IF(IPROB.EQ.1) LL=NGRO-LLL+1
      KPLIB=LCMGIL(JPLIB,LL)
      CALL LCMGET(KPLIB,'NTOT0',SIGT2(1,LLL))
*----
*  TRANSPORT CORRECTION
*----
      IF(ITRANZ.NE.0) THEN
         CALL LCMGET(KPLIB,'TRANC',SIGT3(1,LLL))
      ELSE
         CALL XDRSET(SIGT3(1,LLL),NBMIX,0.0)
      ENDIF
*
      CALL XDLSET(NOCONV(1,LLL),NBMIX,.TRUE.)
   70 CONTINUE
      CALL LCMSIX(IPLIB,' ',2)
      IF(IMPX.GE.5) THEN
         WRITE (6,'(/20H SHIBA INPUT VALUES:/)')
         DO 80 L=IGRMIN,IGRMAX
         WRITE(6,420) L
         WRITE(6,460) (SN(L,J),J=1,NBISO)
         WRITE(6,480) (SIGT2(IBM,L),IBM=1,NBMIX)
   80    CONTINUE
         WRITE(6,490)
      ENDIF
*----
*  ELIMINATE ISOTOPE ABSENT FROM GEOMETRY
*----
      DO IBM=1,NBMIX
        DO IREG=1,NREG
          IF(MAT(IREG).EQ.IBM) GO TO 85
        ENDDO
        DO ISO=1,NBISO
          IF(MIX(ISO).EQ.IBM) LSHI(ISO)=0
        ENDDO
   85   CONTINUE
      ENDDO
*----
*  RECOMPUTE THE VECTOR LSHI
*----
      IF(LEVEL.EQ.0) THEN
         NRES2=NRES
      ELSE
         ALLOCATE(LSHI2(NBISO))
         NRES1=0
         NRES2=0
         DO 90 ISO=1,NBISO
         LSHI2(ISO)=0
   90    CONTINUE
         DO 140 INRS=1,NRES
  100    DENMAX=0.0
         KSOT=0
         DO 120 ISO=1,NBISO
         IF(LSHI2(ISO).EQ.0) THEN
            VOLISO=0.0
            DO 110 I=1,NREG
            IF(MAT(I).EQ.MIX(ISO)) VOLISO=VOLISO+VOL(I)
  110       CONTINUE
            IF((LSHI(ISO).EQ.INRS).AND.(DEN(ISO)*VOLISO.GT.DENMAX)) THEN
               KSOT=ISO
               DENMAX=DEN(ISO)*VOLISO
            ENDIF
         ENDIF
  120    CONTINUE
         IF(KSOT.GT.0) THEN
           NRES2=NRES2+1
           DO 130 ISO=1,NBISO
           IF((ISONRF(1,ISO).EQ.ISONRF(1,KSOT)).AND.
     1        (ISONRF(2,ISO).EQ.ISONRF(2,KSOT)).AND.
     2        (ISONRF(3,ISO).EQ.ISONRF(3,KSOT)).AND.
     3        (LSHI(ISO).EQ.INRS)) LSHI2(ISO)=NRES2
           IF((ISONAM(1,ISO).EQ.ISONAM(1,KSOT)).AND.
     1        (ISONAM(2,ISO).EQ.ISONAM(2,KSOT)).AND.
     2        (LSHI(ISO).EQ.INRS)) LSHI2(ISO)=NRES2
  130      CONTINUE
           GO TO 100
         ENDIF
         IF(NRES2.EQ.NRES1) THEN
           WRITE(HSMG,'(43HSHIDRV: NO RESONANT ISOTOPES IN RESONANT RE,
     1     11HGION NUMBER,I4,5H (1).)') INRS
           CALL XABORT(HSMG)
         ENDIF
         NRES1=NRES2
  140    CONTINUE
      ENDIF
*----
*  DETERMINE THE AMOUNT OF SCRATCH STORAGE REQUIRED
*----
      NBMIX2=0
      DO 150 ISO=1,NBISO
      IF(LSHI(ISO).GT.0) NBMIX2=NBMIX2+1
  150 CONTINUE
*----
*  ITERATION LOOP
*----
      IF(LEVEL.EQ.0) ALLOCATE(SIGE(NRES2*NGRO))
      NITER=0
  160 NITER=NITER+1
      START=(NITER.EQ.1)
      IF(IMPX.GT.5) WRITE (6,430) NITER
      DO 175 L=IGRMIN,IGRMAX
      DO 170 IBM=1,NBMIX
      SIGT1(IBM,L)=SIGT2(IBM,L)
  170 CONTINUE
  175 CONTINUE
      IF(LEVEL.EQ.0) THEN
         CALL SHISN2 (IPLIB,IPTRK,IFTRAK,NGRO,NBISO,NBMIX,NREG,NUN,
     1   CDOOR,NRES,NBMIX2,IMPX,ISONAM,MIX,DEN,SN,SB,LSHI,IPHASE,
     2   MAT,VOL,KEYFLX,LEAKSW,TITR,START,SIGT2,SIGT3,NOCONV,BIEFF,
     3   LGC,SIGE)
      ELSE
         DO 210 INRS=1,NRES2
         NBNRS=0
         DO 200 IBM=1,NBMIX
         LOGDO=.FALSE.
         DO 180 I=1,NREG
         LOGDO=LOGDO.OR.(MAT(I).EQ.IBM)
  180    CONTINUE
         IF(.NOT.LOGDO) GO TO 200
         DO 190 ISO=1,NBISO
         IF((MIX(ISO).EQ.IBM).AND.(LSHI2(ISO).EQ.INRS)) THEN
            NBNRS=NBNRS+1
            GO TO 200
         ENDIF
  190    CONTINUE
  200    CONTINUE
         IF(NBNRS.EQ.0) THEN
            IF(START.AND.(IMPX.GE.1)) WRITE(6,385) 'SHIDRV',INRS
            GO TO 210
         ELSE IF(START.AND.(NBNRS.GT.1).AND.(IMPX.GE.1)) THEN
            WRITE (6,370) NBNRS,INRS
         ENDIF
         CALL SHISN3 (IPLIB,IPTRK,IFTRAK,LEVEL,NGRO,NBISO,NBMIX,NREG,
     1   NUN,CDOOR,INRS,NBNRS,IMPX,ISONAM,MIX,DEN,SN,SB,LSHI2,IPHASE,
     2   MAT,VOL,KEYFLX,LEAKSW,TITR,START,SIGT2,SIGT3,NOCONV,BIEFF,LGC)
  210    CONTINUE
      ENDIF
      ZZMAX=0.0
      LNGRO=0
      ICOUNT=0
      DO 240 L=IGRMIN,IGRMAX
      ZNORM=0.0
      DO 220 IBM=1,NBMIX
      ZNORM=MAX(ZNORM,ABS(SIGT2(IBM,L)))
  220 CONTINUE
      ZMAX=0.0
      MASKL(L)=.FALSE.
      DO 230 IBM=1,NBMIX
      YMAX=ABS(SIGT1(IBM,L)-SIGT2(IBM,L))/ZNORM
      ZMAX=MAX(ZMAX,YMAX)
      NOCONV(IBM,L)=(NOCONV(IBM,L).AND.(YMAX.GT.EPS))
      MASKL(L)=MASKL(L).OR.NOCONV(IBM,L)
  230 CONTINUE
      IF(MASKL(L)) ICOUNT=ICOUNT+1
      IF(ZMAX.GT.ZZMAX) THEN
         ZZMAX=ZMAX
         LNGRO=L
      ENDIF
  240 CONTINUE
      IF(IMPX.GE.2) WRITE (6,440) NITER,ICOUNT,ZZMAX,LNGRO
      IF(IMPX.GE.10) THEN
         WRITE (6,450) (L,MASKL(L),L=IGRMIN,IGRMAX)
         WRITE (6,'(/31H INPUT MACROSCOPIC X-S IN GROUP,I4,1H:)') LNGRO
         WRITE (6,'(1X,1P,10E12.4)') (SIGT1(IBM,LNGRO),IBM=1,NBMIX)
         WRITE (6,'(/32H OUTPUT MACROSCOPIC X-S IN GROUP,I4,1H:)') LNGRO
         WRITE (6,'(1X,1P,10E12.4)') (SIGT2(IBM,LNGRO),IBM=1,NBMIX)
      ENDIF
      IF(IMPX.GT.3) THEN
         WRITE (6,'(/29H OUTPUT DILUTION X-S IN GROUP,I4,1H:)') LNGRO
         WRITE (6,'(1X,1P,10E12.4)') (SN(LNGRO,J),J=1,NBISO)
      ENDIF
      IF((NITER.GE.MAXX0).AND.(ICOUNT.GT.0)) THEN
         WRITE (6,390)
         GO TO 250
      ELSE IF(ICOUNT.GT.0) THEN
         GO TO 160
      ENDIF
*----
*  CONVERGENCE IS OBTAINED
*----
  250 IF(LEVEL.GT.0) DEALLOCATE(LSHI2)
      IF(LEVEL.EQ.0) DEALLOCATE(SIGE)
      IF(IMPX.GE.3) THEN
         WRITE (6,'(/21H SHIBA OUTPUT VALUES:/)')
         DO 260 L=IGRMIN,IGRMAX
         WRITE(6,420) L
         WRITE(6,460) (SN(L,J),J=1,NBISO)
         IF(BIEFF) WRITE(6,470) (SB(L,J),J=1,NBISO)
         IF(IMPX.GE.5) WRITE(6,480) (SIGT2(IBM,L),IBM=1,NBMIX)
  260    CONTINUE
         WRITE(6,490)
      ENDIF
*----
*  COMPUTE THE NEW SELF-SHIELDED MACROSCOPIC CROSS SECTIONS
*----
      CALL XDLSET(MASKL,NGRO,.FALSE.)
      DO 280 LLL=IGRMIN,IGRMAX
      MASKL(LLL)=.TRUE.
  280 CONTINUE
      DO 300 IBM=1,NBMIX
      DO 290 ISO=1,NBISO
      MASK(IBM)=(MIX(ISO).EQ.IBM).AND.(LSHI(ISO).GT.0)
      IF(MASK(IBM)) GO TO 300
  290 CONTINUE
  300 CONTINUE
      ITSTMP=0
      TMPDAY(1)=0.0
      TMPDAY(2)=0.0
      TMPDAY(3)=0.0
      CALL LIBMIX(IPLIB,NBMIX,NGRO,NBISO,ISONAM,MIX,DEN,MASK,MASKL,
     > ITSTMP,TMPDAY)
      IF(IMPX.GT.0) WRITE (6,500) NITER,ZZMAX
*----
*  STORE THE GENERAL SHIBA PARAMETERS
*----
      CALL XDISET(IPAR,NSTATE,0)
      IPAR(1)=IGRMIN
      IPAR(2)=IGRMAX
      IPAR(3)=MAXX0
      IPAR(4)=IBIEFF
      IPAR(5)=IGC
      IPAR(6)=ITRANZ
      IPAR(7)=LEVEL
      IPAR(8)=IPHASE
      CALL LCMSIX(IPLIB,'SHIBA',1)
      CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,IPAR)
      CALL LCMPUT(IPLIB,'EPS-SHIBA',1,2,EPS)
      CALL LCMSIX(IPLIB,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(NOCONV,MASKL,MASK)
      DEALLOCATE(SIGT3,SIGT2,SIGT1)
      RETURN
*
  340 FORMAT(//26H VOLUME NB.              :,11(I5,3X,1HI))
  350 FORMAT(  26H MIXTURE (MAT)           :,11(I5,3X,1HI))
  370 FORMAT(/42H SHIDRV: USE THE NORDHEIM MODEL TO PROCESS,I3,5H RESO,
     1 39HNANT MIXTURES IN RESONANT REGION NUMBER,I3,1H.)
  380 FORMAT(32HSHIDRV: INVALID VALUE OF NBMIX (,I5,2H).)
  385 FORMAT(A6,1X,': RESONANT REGION =',I10,1X,'NOT USED.')
  390 FORMAT(/1X,61(1H*)/42H SHIDRV: MAXIMUM NUMBER OF SELF-SHIELDING ,
     1 20HITERATIONS EXCEEDED./1X,61(1H*)/)
  400 FORMAT(
     > 1X,'SHIBA MULTIDIMENSIONAL SELF-SHIELDING CALCULATION',
     > 1X,'-> A. HEBERT'/
     > 1X,A72/
     > 1X,'COLLISION PROBABILITY MODULE: ',A12/)
  405 FORMAT(/8H OPTIONS/8H -------/
     1  7H IGRMIN,I8,27H   (FIRST GROUP TO PROCESS)/
     2  7H IGRMAX,I8,34H   (MOST THERMAL GROUP TO PROCESS)/
     3  7H MAXX0 ,I8,33H   (MAXIMUM NUMBER OF ITERATIONS)/
     4  7H IBIEFF,I8,46H   (=1: USE LIVOLANT-JEANPIERRE NORMALIZATION)/
     5  7H IGC   ,I8,42H   (=1: USE GOLDSTEIN-COHEN APPROXIMATION)/
     6  7H ITRANZ,I8,45H   (0=NO TRANSPORT CORRECTION/1=APOLLO TYPE/2,
     7  57H=RECOVER FROM LIBRARY/3=WIMS-D TYPE/4=LEAKAGE CORRECTION)/
     8  7H LEVEL ,I8,46H   (=0: STAMM'LER; =1: STAMM'LER/NORDHEIM; =2:,
     9  18H RIEMANN/NORDHEIM)/
     1  7H IPHASE,I8,37H   (=1: NATIVE ASSEMBLY; =2: USE PIJ)/
     2  7H EPS   ,1P,E8.1,22H   (STOPING CRITERION)/)
  410 FORMAT(/52H SHIDRV: THERE IS NO REQUEST TO PROCESS ANY RESONANT,
     1 9H ISOTOPE./)
  420 FORMAT(/1X,131(1H-)//8H GROUP =,I4/)
  430 FORMAT(/40H PERFORMING SELF-SHIELDING ITERATION NB.,I5)
  440 FORMAT(/27H SELF-SHIELDING ITERATION =,I4,5X,14HNUMBER OF NON ,
     1 18HCONVERGED GROUPS =,I4,5X,7HERROR =,1P,E13.4,0P,9H IN GROUP,
     2 I4/)
  450 FORMAT(7H MASKL(,I3,2H)=,L1,:,8H  MASKL(,I3,2H)=,L1,:,8H  MASKL(,
     1 I3,2H)=,L1,:,8H  MASKL(,I3,2H)=,L1,:,8H  MASKL(,I3,2H)=,L1,:,
     2 8H  MASKL(,I3,2H)=,L1,:,8H  MASKL(,I3,2H)=,L1,:,8H  MASKL(,I3,
     3 2H)=,L1,:,8H  MASKL(,I3,2H)=,L1)
  460 FORMAT(/37H MICROSCOPIC DILUTION CROSS SECTIONS:/(9X,1P,11E11.3))
  470 FORMAT(/53H LIVOLANT AND JEANPIERRE MICROSCOPIC DILUTION CROSS S,
     1 8HECTIONS:/(9X,1P,11E11.3))
  480 FORMAT(/34H MACROSCOPIC TOTAL CROSS SECTIONS:/(9X,1P,11E11.3))
  490 FORMAT(/1X,131(1H-)/)
  500 FORMAT(/41H CONVERGENCE REACHED AT SHIBA ITERATION =,I4,7H  ERROR,
     1 2H =,1P,E11.3/)
      END
