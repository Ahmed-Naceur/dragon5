*DECK CPODRV
      SUBROUTINE CPODRV(IPCPO ,IPEDIT,IPDEPL,IPRINT,CURNAM,CTITRE,
     >                  NAMCPO,NGROUP,NMERGE,NBMICR,NIFISS,MXBURN,
     >                  NL    ,NISCPO,NPROC ,ILEAKS,NXXXZ ,NEDMAC,
     >                  HVECT ,NSBS  ,ILOCAL,ISOCPO,ISOTMP,IDIMIX,
     >                  NBIMRG,ICOMIX,VOLMER,ENERGY,TIME  ,BURN  ,
     >                  WIRRAD,IBSTEP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover cross section information located on directory CURNAM or on
* directory family with prefix CURNAM.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input
* IPCPO   pointer to the compo (L_COMPO signature).
* IPEDIT  pointer to edit information (L_EDIT signature).
* IPDEPL  pointer to depletion information (L_BURNUP signature).
* IPRINT  print parameter. Equal to zero for no print.
* CURNAM  name of the output directory (or prefix of output
*         directory in burnup cases).
* CTITRE  character*72 title.
* NAMCPO  character*8 name of the material mixture sub-directory.
* NGROUP  number of energy groups in output data.
* NMERGE  number of output regions.
* NBMICR  maximum number of isotopes.
* NIFISS  number of fissile isotopes.
* MXBURN  maximum number of output burnup sets.
* NL      number of Legendre orders (=1 for isotropic scattering).
* NISCPO  number of Compo isotopes treated.
* NPROC   number of microscopic xs to process.
* ILEAKS  leak option: 0 no leakage ; 1 homogeneous leakage ;
*         2 heterogeneous leakage.
* NXXXZ   maximum dimension of ISO dependent vector = max(nbmicr,1).
* NEDMAC  number of edit xs.
* HVECT   name of edit xs.
* NSBS    number of sub-burnup step considered.
* ILOCAL  local parameter flag (0: global; 1:local).
* ISOCPO  Compo name of isotopes.
* ISOTMP  name of isotopes in EDIT.
* IDIMIX  isotopes identifier in each Compo material.
* NBIMRG  final number of isotope per region.
* ICOMIX  pointer to Compo isotope for region.
* VOLMER  merge volume.
* ENERGY  energy.
* TIME    time steps.
* BURN    burnup.
* WIRRAD  irradiation.
* IBSTEP  sub-burnup step considered.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPCPO,IPEDIT,IPDEPL
      INTEGER          IPRINT,NGROUP,NMERGE,CTITRE(18),NBMICR,NIFISS,
     >                 MXBURN,NL,NISCPO,NPROC ,ILEAKS,NXXXZ,NEDMAC,
     >                 NSBS  ,ILOCAL,ISOCPO(3,NXXXZ),ISOTMP(3,NXXXZ),
     >                 IDIMIX(NMERGE,NXXXZ),NBIMRG(NMERGE),
     >                 ICOMIX(NMERGE,NXXXZ),IBSTEP(MXBURN)
      CHARACTER        CURNAM*12,NAMCPO*8,HVECT(NEDMAC)*8
      REAL             VOLMER(NMERGE),ENERGY(NGROUP+1),
     >                 TIME(MXBURN),BURN(MXBURN),WIRRAD(MXBURN)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INDPRO,ITYPRO,NAMI
      REAL, ALLOCATABLE, DIMENSION(:) :: DENTMP,EMJMAC,VECT,XSREC,XSCAT,
     1 DISFC,DENSI,EMJI
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: RVALOC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DENCPO,XSREM,SCREM
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DXSMIC,DMJCPO
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DSCMIC,DXSMAC,
     1 DISFAC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DSCMAC
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT,NSTATE,NDPROC,NPARAM
      REAL             CUTOFF
      PARAMETER       (IOUT=6,NSTATE=40,NDPROC=20,NPARAM=4,
     >                 CUTOFF=1.0E5)
      INTEGER          NEFBRN,IBR,IBR2,IBURN,NDUM1,NDUM2,MAXDM,
     >                 IMRG,ISOC,ISOR,ITC,IDFLU,
     >                 ISTATE(NSTATE),IPARAM(NPARAM),MXISOS,
     >                 NBISO,NREAC,NVAR,NBMIX,NREG,NLOC
      CHARACTER        NAMMIX*12,NAMBRN*12,NAMISO*12,NAMMAC*12
      REAL             DELTA(2),TMPDAY(3),DELERR(3)
      DOUBLE PRECISION DMJMAC
      INTEGER          IFCDIS
*----
*  SCRATCH STORAGE ALLOCATION
*   DENTMP  density of EDI isotopes.
*   EMJMAC  fission energy for macroscopic data.
*   DENCPO  density of Compo isotopes.
*   DMJCPO  fission energy for macroscopic data.
*   INDPRO  identifier for xs processing.
*   ITYPRO  identifier for xs processed.
*   DXSMIC  micro vector xs.
*   DSCMIC  micro scattering matrix xs.
*   DXSMAC  macro vector xs.
*   DSCMAC  macro scattering matrix xs.
*   DISFAC  discontinuity factors.
*   RVALOC  local burnup and irradiation values.
*----
      ALLOCATE(INDPRO(NPROC),ITYPRO(NPROC))
      ALLOCATE(DENTMP(NXXXZ),EMJMAC(NMERGE),RVALOC(2,NMERGE,MXBURN))
      ALLOCATE(DENCPO(NXXXZ),DXSMIC(NGROUP,NPROC),
     > DSCMIC(NGROUP,NGROUP,NL),DXSMAC(NGROUP,NPROC,NMERGE),
     > DSCMAC(NGROUP,NGROUP,NL,NMERGE),DMJCPO(2,NXXXZ),
     > DISFAC(2,NGROUP,3))
*----
*  GET GLOBAL BURNUP AND IRRADIATION
*----
      IFCDIS=1
      NAMMAC='MACR       '
      IF(NSBS.EQ.0) THEN
        BURN(1)=0.0
        WIRRAD(1)=0.0
        NEFBRN=1
        IBSTEP(NEFBRN)=0
      ELSE
        DO 100 IBR=1,NSBS
          IBURN=IBSTEP(IBR)
          IF(IBURN.GT.0.AND. IBURN.LE.MXBURN) THEN
            WRITE(NAMBRN,'(8HDEPL-DAT,I4.4)') IBURN
            CALL LCMSIX(IPDEPL,NAMBRN,1)
            CALL LCMGET(IPDEPL,'BURNUP-IRRAD',DELTA)
            BURN(IBR)=DELTA(1)
            WIRRAD(IBR)=DELTA(2)
            TIME(IBR)=TIME(IBURN)
            CALL LCMSIX(IPDEPL,NAMBRN,2)
          ENDIF
 100    CONTINUE
        NEFBRN=NSBS
      ENDIF
*----
*  GET LOCAL BURNUP AND IRRADIATION
*----
      IF((NSBS.EQ.0).OR.(ILOCAL.EQ.0).OR.(.NOT.C_ASSOCIATED(IPDEPL)))
     1 THEN
         CALL XDRSET(RVALOC,2*NMERGE*NEFBRN,0.0)
      ELSE
         NLOC=2
         CALL LCMGET(IPDEPL,'STATE-VECTOR',ISTATE)
         IF(ISTATE(3).NE.MXBURN) CALL XABORT('CPODRV: INVALID STATE-VE'
     1   //'CTOR.')
         NBISO=ISTATE(4)
         NREAC=ISTATE(6)
         NVAR=ISTATE(7)
         NBMIX=ISTATE(8)
         CALL LCMGET(IPEDIT,'STATE-VECTOR',ISTATE)
         NREG=ISTATE(17)
         DO 105 IBR=1,NSBS
           IBURN=IBSTEP(IBR)
           CALL COMGEN(IPDEPL,IPEDIT,NREG,NMERGE,IBURN,'FLUB',MXBURN,
     1     NBMIX,NBISO,NREAC,NVAR,1,NLOC,RVALOC(1,1,IBR))
           CALL COMGEN(IPDEPL,IPEDIT,NREG,NMERGE,IBURN,'IRRA',MXBURN,
     1     NBMIX,NBISO,NREAC,NVAR,2,NLOC,RVALOC(1,1,IBR))
  105    CONTINUE
      ENDIF
*----
*  INITIALIZE INDPRO FOR MICROSCOPIC XS TO PROCESS
*----
      ALLOCATE(VECT(NEDMAC))
      CALL XDRSET(VECT,NEDMAC,0.0)
      CALL CPONED(NPROC ,0,NL-1,MAX(1,ILEAKS),NEDMAC,HVECT,VECT,INDPRO)
      INDPRO(6)=0
      INDPRO(16)=0
      ALLOCATE(XSREC(NGROUP*NPROC),XSCAT(NGROUP*NGROUP*NL))
*----
*  LOOP OVER BURNUP STEPS
*----
      NDUM1=NMERGE*MAX(NIFISS,NGROUP)
      NDUM2=NGROUP*NGROUP
      MAXDM=MAX(NDUM1,NDUM2)
      ALLOCATE(DISFC(NGROUP))
      ALLOCATE(XSREM(NGROUP*NPROC),SCREM(NGROUP*NGROUP*NL))
      IDFLU=16
      MXISOS=0
      DO 110 IBR=1,NEFBRN
        WRITE(NAMBRN,'(A8,I4)') 'BURN    ',IBR
        IBURN=IBSTEP(IBR)
        IF(IBURN.GT.0) WRITE(CURNAM(9:12),'(I4.4)') IBURN
        CALL LCMSIX(IPEDIT,CURNAM,1)
        IF(NISCPO.GT.0) CALL LCMGET(IPEDIT,'ISOTOPESDENS',DENTMP)
        IF(IPRINT.GE.10) WRITE(IOUT,6000) CURNAM
        CALL LCMSIX(IPEDIT,'MACROLIB',1)
        CALL LCMGET(IPEDIT,'TIMESTAMP',TMPDAY)
        DELERR(1)=CUTOFF*ABS(TMPDAY(1)-TIME(IBR)/8.64E-4)
        DELERR(2)=CUTOFF*ABS(TMPDAY(2)-BURN(IBR))
        DELERR(3)=CUTOFF*ABS(TMPDAY(3)-WIRRAD(IBR))
        IF( (DELERR(1).GT.TMPDAY(1)) .OR.
     >      (DELERR(2).GT.TMPDAY(2)) .OR.
     >      (DELERR(3).GT.TMPDAY(3)) ) THEN
          IF(TIME(IBR)   .EQ. 0.0 .AND.
     >       BURN(IBR)   .EQ. 0.0 .AND.
     >       WIRRAD(IBR) .EQ. 0.0) THEN
            WRITE(IOUT,7001)
          ELSE
            WRITE(IOUT,7000)
     >        TMPDAY(1),TIME(IBR)/8.64E-4,DELERR(1)/CUTOFF,
     >        TMPDAY(2),BURN(IBR),DELERR(2)/CUTOFF,
     >        TMPDAY(3),WIRRAD(IBR)      ,DELERR(3)/CUTOFF
          ENDIF
          TIME(IBR)=TMPDAY(1)*8.64E-4
          BURN(IBR)=TMPDAY(2)
          WIRRAD(IBR)=TMPDAY(3)
        ENDIF
*----
*  READ MACROSCOPIC XS FOR ALL GROUP AND ALL REGIONS
*----
        CALL CPOMAR(IPEDIT,NGROUP,NMERGE,NL    ,NIFISS,NEDMAC,
     >              HVECT ,VECT ,NPROC ,ILEAKS,DXSMAC,
     >              DSCMAC,EMJMAC,DISFC,IFCDIS,DISFAC)
        CALL LCMSIX(IPEDIT,'MACROLIB',2)
        CALL XDDSET(DENCPO,NISCPO,0.0D0)
        DO 120 IMRG=1,NMERGE
          DMJMAC=DBLE(EMJMAC(IMRG))
          WRITE(NAMMIX,'(A8,I4)') NAMCPO,IMRG
          CALL LCMSIX(IPCPO,NAMMIX,1)
          CALL LCMSIX(IPCPO,NAMBRN,1)
          IF(IPRINT.GE.10) WRITE(IOUT,6001) NAMMIX
          CALL XDDSET(DMJCPO,2*NBMICR,0.0D0)
          CALL XDDSET(XSREM,NGROUP*NPROC,0.0D0)
          CALL XDDSET(SCREM,NGROUP*NGROUP*NL,0.0D0)
          DO 130 ISOC=1,NBIMRG(IMRG)
            ISOR=ICOMIX(IMRG,ISOC)
            WRITE(NAMISO,'(3A4)') (ISOCPO(ITC,ISOR),ITC=1,3)
            IF(IPRINT.GE.10) WRITE(IOUT,6002) NAMISO
*----
*   CREATE AND SAVE XS FOR A CPO ISOTOPE IN CURRENT REGION
*----
            CALL LCMSIX(IPCPO,NAMISO,1)
            CALL XDDSET(DXSMIC,NGROUP*NPROC,0.0D0)
            CALL XDDSET(DSCMIC,NGROUP*NGROUP*NL,0.0D0)
            CALL CPOMIC(IPCPO ,IPEDIT,IPRINT,NGROUP,NMERGE,NBMICR,
     >                  NL    ,IMRG  ,ISOR  ,NPROC ,ISOTMP,IDIMIX,
     >                  INDPRO,ITYPRO,DENCPO,DENTMP,DXSMIC,DSCMIC,
     >                  DMJCPO,DXSMAC(1,IDFLU,IMRG))
            CALL LCMSIX(IPCPO,NAMISO,2)
*----
*  REMOVE CONTRIBUTION OF CPO ISOTOPE FROM MACROSCOPIC.
*----
            IF(DENCPO(ISOR).GT.0.0D0) THEN
              CALL CPOREM(NGROUP,NL    ,NPROC ,INDPRO,DENCPO(ISOR),
     >                    DXSMIC,DSCMIC,XSREM ,SCREM )
            ENDIF
            DMJMAC=DMJMAC-DMJCPO(1,ISOR)
 130      CONTINUE
*----
*  WRITE MACROSCOPIC XS FOR ALL GROUP IN THIS REGION REGIONS
*----
          CALL CPOMAW(IPCPO ,IPRINT,NGROUP,NL    ,NPROC ,INDPRO,
     >                ITYPRO,DXSMAC(1,1,IMRG),DSCMAC(1,1,1,IMRG),
     >                XSREM,SCREM,DISFC,DMJMAC,IFCDIS,DISFAC)
          ALLOCATE(DENSI(NBIMRG(IMRG)+1),EMJI(NBIMRG(IMRG)+1))
          DENSI(1)=1.0
          EMJI=REAL(DMJMAC)*1.0E-18
          DO 140 ISOC=1,NBIMRG(IMRG)
            ISOR=ICOMIX(IMRG,ISOC)
            DENSI(ISOC+1)=REAL(DENCPO(ISOR))
            IF(DMJCPO(2,ISOR).GT.0.0D0) THEN
              EMJI(ISOC+1)=1.0E-18*REAL(DMJCPO(1,ISOR)/DMJCPO(2,ISOR))
            ELSE
              EMJI(ISOC+1)=0.0
            ENDIF
 140      CONTINUE
          CALL LCMPUT(IPCPO,'ISOTOPESDENS',(NBIMRG(IMRG)+1),2,DENSI)
          CALL LCMPUT(IPCPO,'ISOTOPES-EFJ',(NBIMRG(IMRG)+1),2,EMJI)
          DEALLOCATE(EMJI,DENSI)
          CALL LCMSIX(IPCPO,NAMBRN,2)
*----
*  PUT REMAINING INFORMATION ON CPO FOR THIS MIXTURE
*----
          CALL LCMPUT(IPCPO,'TITLE',18,3,CTITRE)
          CALL LCMPUT(IPCPO,'VOLUME',1,2,VOLMER(IMRG))
          CALL LCMPUT(IPCPO,'ENERGY',NGROUP+1,2,ENERGY)
          IF(IBR.EQ.NEFBRN) THEN
            IF(ILOCAL.EQ.1) THEN
              DO 145 IBR2=1,NEFBRN
                WIRRAD(IBR2)=RVALOC(1,IMRG,IBR2)
                BURN(IBR2)=RVALOC(2,IMRG,IBR2)
 145          CONTINUE
            ENDIF
            IF(IPRINT.GT.1) THEN
               WRITE(IOUT,7002) IMRG,'IRRA',(WIRRAD(IBR2),IBR2=1,NEFBRN)
               WRITE(IOUT,7002) IMRG,'BURN',(BURN(IBR2),IBR2=1,NEFBRN)
            ENDIF
            CALL LCMPUT(IPCPO,'N/KB  ',NEFBRN,2,WIRRAD)
            CALL LCMPUT(IPCPO,'BURNUP',NEFBRN,2,BURN)
            ALLOCATE(NAMI(3*(NBIMRG(IMRG)+1)))
            READ(NAMMAC,'(3A4)') (NAMI(ITC+1),ITC=0,2)
            ITC=3
            DO 150 ISOC=1,NBIMRG(IMRG)
              ISOR=ICOMIX(IMRG,ISOC)
              NAMI(ITC+1)=ISOCPO(1,ISOR)
              NAMI(ITC+2)=ISOCPO(2,ISOR)
              NAMI(ITC+3)=ISOCPO(3,ISOR)
              ITC=ITC+3
 150        CONTINUE
            CALL LCMPUT(IPCPO,'ISOTOPESNAME',3*(NBIMRG(IMRG)+1),3,NAMI)
            CALL XDISET(NAMI,(NBIMRG(IMRG)+1),0)
            CALL LCMPUT(IPCPO,'JTAB',(NBIMRG(IMRG)+1),1,NAMI)
            DEALLOCATE(NAMI)
            CALL XDISET(IPARAM,NPARAM,0)
            IPARAM(1)=NGROUP
            IPARAM(2)=NBIMRG(IMRG)+1
            IPARAM(3)=NL
            IPARAM(4)=NEFBRN
            MXISOS=MAX(MXISOS,NBIMRG(IMRG)+1)
            CALL LCMPUT(IPCPO,'PARAM',NPARAM,1,IPARAM)
          ENDIF
          CALL LCMSIX(IPCPO,NAMMIX,2)
 120    CONTINUE
        CALL LCMSIX(IPEDIT,CURNAM,2)
 110  CONTINUE
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NMERGE
      ISTATE(2)=NGROUP
      ISTATE(3)=MXISOS
      ISTATE(4)=NL
      ISTATE(5)=NEFBRN
      ISTATE(6)=NPARAM
      ISTATE(7)=IFCDIS
      CALL LCMPUT(IPCPO,'STATE-VECTOR',NSTATE,1,ISTATE)
      DEALLOCATE(SCREM,XSREM,DISFC)
      IF(NISCPO.GT.0) DEALLOCATE(XSCAT,XSREC,VECT)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DISFAC,DMJCPO,DSCMAC,DXSMAC,DSCMIC,DXSMIC,DENCPO)
      DEALLOCATE(RVALOC,EMJMAC,DENTMP)
      DEALLOCATE(ITYPRO,INDPRO)
      RETURN
*----
*  PRINT FORMAT
*----
 6000 FORMAT(' CPODRV: STEPPING UP ON DIRECTORY = ',A12)
 6001 FORMAT(' CPODRV: CREATING MIXTURE         = ',A12)
 6002 FORMAT(' CPODRV: CREATING ISOTOPE         = ',A12)
*----
*  WARNING FORMAT
*----
 7000 FORMAT(
     > ' CPODRV: WARNING -> BURNUP AND EDIT DATA DIFFER',1P/
     > ' TIME: EDIT=',E15.7,5X,' BURNUP=',E15.7,' DIFF=',E15.7/
     > ' BURN: EDIT=',E15.7,5X,' BURNUP=',E15.7,' DIFF=',E15.7/
     > ' WIRR: EDIT=',E15.7,5X,' BURNUP=',E15.7,' DIFF=',E15.7/
     > ' USE EDIT DATA   ')
 7001 FORMAT(
     > ' CPODRV: WARNING -> 0 BURNUP STEP, USE EDIT DATA')
 7002 FORMAT(/13H CPODRV: MIX=,I4,3X,A,1H=,1P,6E12.4/(25X,6E12.4))
      END
