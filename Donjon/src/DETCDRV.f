*DECK DETCDRV
      SUBROUTINE DETCDRV(IPDET,NGRP,NEL,NUN,NX,NY,NZ,MESHX,MESHY,MESHZ,
     1           KEYF,FLUX,IPRT,KC,DT,LHEX,LSIMEX,LNORM,VNORM,LPARAB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the module DETECT:
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin, M. Guyot
*
*Parameters:
* IPDET  pointer to the library object
* NGRP   number of energy groups
* NEL    number of finite elements
* NUN    number of unknowns
* NX     number of x mesh-splitted elements 
* NY     number of y mesh-splitted elements 
* NZ     number of z mesh-splitted elements
* MESHX  
* MESHY  
* MESHZ  
* KEYF   keyflux recover from L_TRACk object
* FLUX   flux for each mesh-splitted elements
* IPRT   printing index
* KC     calculation type reference
* DT     time step
* LHEX   =.TRUE. if hexagonal detectors are present
* LSIMEX =.TRUE. if keyword SIMEX is present
* LNORM  =.TRUE. if keyword NORM is present
* VNORM  real used for normalization
* LPARAB =.TRUE. if parabolic interpolation is performed
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDET
      INTEGER IPRT,KC,NGRP,NEL,NX,NY,NZ,KEYF(NEL),NUN
      LOGICAL LHEX,LSIMEX,LNORM,LPARAB
      REAL    FLUX(NUN,NGRP),DT,MESHX(NX+1),MESHY(NY+1),MESHZ(NZ+1),
     1        VNORM
*----
*  LOCAL VARIABLES
*----
      INTEGER   NSTATE,IOUT
      PARAMETER (NSTATE=40,IOUT=6) 
      INTEGER   ILONG,ITYLCM,INFO(2),NREP,ITHEX,NHEX,J
      REAL      DEVPOS(6),PLNL,REF,RESP,VLAMDA,XMULT,TLG
      CHARACTER NXTYP*12,FIRST*12,NXDET*12,FIRST1*12
      LOGICAL   LREGUL
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IHEX
      REAL, ALLOCATABLE, DIMENSION(:) :: SPEC,REP,FRACT,NVCST,PDD,APD,
     1          BPD
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SPEC(NGRP))
*
      IF(IPRT.GE.1) THEN
        IF(KC.EQ.0) THEN
          WRITE(IOUT,*) 'DETECT: CALCULATION TYPE REFERENCE, KC = ',KC
        ELSE
          WRITE(IOUT,*) 'DETECT: CALCULATION TYPE NORMAL, KC = ',KC
        ENDIF
        WRITE(IOUT,*) 'DETECT: TIME STEP USED, DT = ',DT
      ENDIF
      CALL LCMSIX(IPDET,' ',0)
      NXTYP = ' '
      CALL LCMNXT(IPDET,NXTYP)
      FIRST = NXTYP
 10   CALL LCMNXT(IPDET,NXTYP)
      CALL LCMLEN(IPDET,NXTYP,ILONG,ITYLCM)
      IF (ITYLCM.EQ.0) THEN
        CALL LCMSIX(IPDET,NXTYP,1)
        CALL LCMGET(IPDET,'INFORMATION',INFO)
        CALL LCMGET(IPDET,'SPECTRAL',SPEC)
        NREP = INFO(2)
        ALLOCATE(REP(NREP))
        IF( NXTYP(1:5).EQ.'PLATN' )THEN
          ALLOCATE(FRACT(NREP-1),NVCST(NREP-2),PDD(NREP-2),APD(NREP-2),
     1    BPD(NREP-2))
          CALL LCMGET(IPDET,'FRACTION',FRACT)
          CALL LCMGET(IPDET,'INV-CONST',NVCST)
        ENDIF
        NXDET = ' '
        CALL LCMNXT(IPDET,NXDET)
        FIRST1 = NXDET
 20     CALL LCMNXT(IPDET,NXDET)
        CALL LCMLEN(IPDET,NXDET,ILONG,ITYLCM)
        IF (ITYLCM.EQ.0) THEN
           IF(IPRT.GT.3) WRITE(IOUT,*) 'NAME DETECTOR ',NXDET
           CALL LCMSIX(IPDET,NXDET,1)
           IF(LHEX) THEN
             CALL LCMLEN(IPDET,'NHEX',NHEX,ITHEX)
             IF(NHEX.EQ.0) CALL XABORT('@DETCDRV: HEXAGON NUMBERS'
     +             //' NOT PRESENT IN DETECT')
             ALLOCATE(IHEX(NHEX))
             CALL LCMGET(IPDET,'NHEX',IHEX)
           ENDIF
           CALL LCMGET(IPDET,'POSITION',DEVPOS)
           CALL LCMGET(IPDET,'RESPON',REP)
            IF(LSIMEX.AND.NXTYP.EQ.'VANAD_REGUL') THEN
              CALL DETINT(NX,NY,NZ,NEL,NUN,LPARAB,MESHX,MESHY,MESHZ,
     +        KEYF,FLUX,NGRP,DEVPOS,RESP,IPRT)
            ELSE
              CALL DETFLU(LHEX,NX,NY,NZ,NEL,NUN,MESHX,MESHY,MESHZ,KEYF,
     +        FLUX,NGRP,SPEC,DEVPOS,NHEX,IHEX,RESP,IPRT)
            ENDIF
*----
*  DETECTOR RESPONSE CALCULATION
*----
            IF(.NOT.LNORM)THEN
              PLNL = REP(1)
              REF  = REP(2)
              IF(NXTYP.EQ.'VANAD_REGUL')THEN
*----
*  VANADIUM RESPONSE CALCULATION
*----
                IF(LSIMEX) THEN
                  REF = RESP
                ELSE
                  IF(KC.EQ.1) THEN
                    VLAMDA = 1./225.
                    XMULT = 1.0 + VLAMDA*DT
                    XMULT = 1.0/XMULT
                    RESP = XMULT*(PLNL+DT*VLAMDA*RESP)
                    REF = PLNL
                  ELSE
                    REF = RESP
                  ENDIF
                ENDIF
              ELSEIF(NXTYP(1:5).EQ.'PLATN')THEN
*----
*  PLATINIUM RESPONSE CALCULATION
*----
                LREGUL = .FALSE.
                DO 30 J=1,NREP-2
                  PDD(J) = REP(J)
  30            CONTINUE
                IF(NXTYP.EQ.'PLATN_REGUL')THEN
                  LREGUL = .TRUE.
                ENDIF
                CALL DETPLAT(DT,RESP,REF,KC,PDD,LREGUL,FRACT,NVCST,
     +          NREP-2,APD,BPD)
                DO 40 J=1,NREP-2
                  REP(J) = PDD(J)
  40            CONTINUE
              ELSEIF (NXTYP(1:5).EQ.'CHION') THEN
*----
*  LECTURE DE CHAMBRES D'ION
*----
                IF(NREP.NE.3)CALL XABORT('@DETCDRV: ION CHAMBERS MUST '
     +            //'HAVE THREE STORED VALUES FOR RESPONSES')
                IF (KC.EQ.1) THEN
                  REF = REP(3)
                  RESP = LOG10(RESP/REF)
                  TLG = (RESP-PLNL)/DT
                  REF = TLG
                ELSE
                  REP(3) = RESP
                  RESP =LOG10(RESP/REP(3))
                ENDIF
              ENDIF
            ELSE
              REF = VNORM/RESP
              RESP = VNORM
            ENDIF
*----
*  DETECTOR RESPONSE STORAGE
*----
           REP(1) = RESP
           REP(2) = REF
           IF(IPRT.GT.4) WRITE(6,*) 'RESP, REF ',RESP, REF
           CALL LCMPUT(IPDET,'RESPON',NREP,2,REP)
           CALL LCMSIX(IPDET,' ',2)
        ENDIF
        IF(LHEX) DEALLOCATE(IHEX)
        IF(NXDET.EQ.FIRST1) GOTO 45
        GOTO 20
   45   CALL LCMSIX(IPDET,' ',2)
        DEALLOCATE(REP)
        IF(NXTYP(1:5).EQ.'PLATN')THEN
          DEALLOCATE(FRACT,NVCST,PDD,APD,BPD)
        ENDIF
      ENDIF
      IF (NXTYP.EQ.FIRST) GOTO 50
        GOTO 10
  50  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SPEC)
      RETURN
      END
