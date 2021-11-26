*DECK RESPAR
      SUBROUTINE RESPAR(IPMAP,NCH,NB,NFUEL,NCOMB,NPARM,NX,NY,NZ,NSTATE,
     1 ISTATE,IMPX,NASB,LMAP2,IPMP2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read and store the data related to global and local parameters.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki, R. Chambon, M. Guyot, V. Descotes, B. Toueg
*
*Parameters: input/output
* IPMAP   pointer to fuel-map information.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* NFUEL   number of fuel types.
* NCOMB   number of combustion zones.
* NX      number of elements along x-axis in fuel map.
* NY      number of elements along y-axis in fuel map.
* NZ      number of elements along z-axis in fuel map.
* NSTATE  maximum number of state-vector records.
* IMPX    printing index (=0 for no print).
* NASB    total number of assembly
* LMAP2   flag to set if second fuel-map information is used to
*         recover burnup information
* IPMP2   pointer to the second fuel-map information.
*
*Parameters: output
* ISTATE  updated state-vector.
* NPARM   total number of recorded parameters.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP,IPMP2
      INTEGER NCH,NB,NFUEL,NCOMB,NPARM,NX,NY,NZ,NSTATE,ISTATE(NSTATE),
     1 IMPX
      LOGICAL LMAP2
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER INAME(3),IZONE(NCH),IVECT(NCOMB,NB),NSCH(NCH),
     1 IBSH(NCOMB),MIX(NX*NY*NZ),FMIX(NCH,NB),IAZ(NCH),
     2 ISTAT2(NSTATE)
      REAL VALUE(NCH,NB),POWER(NCH,NB),FPOWER(NB)
      CHARACTER CVALUE(NCH,NB)*12
      DOUBLE PRECISION DFLOT
      CHARACTER TEXT*12,TEXT12*12,PNAME*12,KEYN*12,PNAME2*12
      LOGICAL LRSCH,LBURN
      TYPE(C_PTR) JPMAP,KPMAP,ZPMAP,JPMP2,KPMP2
*----
*  ALLOCATABLE STATEMENTS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ALCH,NUM,IND
      REAL, ALLOCATABLE, DIMENSION(:) :: DENSMOD,BRN,BASS,VAL2,ZZ,VB
*----
*  READ INPUT DATA
*----
      LRSCH=.FALSE.
      LBURN=.FALSE.
      PTOT=0.0
      CALL LCMGET(IPMAP,'FLMIX',FMIX)
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT12,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@RESPAR: CHARACT'
     1 //'ER DATA EXPECTED ('//TEXT12//').')
      IF(TEXT12.EQ.';')THEN
        GOTO 500
*     PRINTING INDEX
      ELSEIF(TEXT12.EQ.'EDIT')THEN
        CALL REDGET(ITYP,NITMA,FLOT,TEXT12,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@RESPAR: INTEGER DA'
     1   //'TA FOR EDIT EXPECTED.')
        IMPX=MAX(0,NITMA)
*----
*  ADD NEW PARAMETER
*----
      ELSEIF(TEXT12.EQ.'ADD-PARAM')THEN
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(TEXT.NE.'PNAME')CALL XABORT('@RESPAR: KEY'
     1   //'WORD PNAME EXPECTED.')
*       READ PARAMETER NAME
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESPAR: CHARACTER'
     1   //' DATA FOR PARAMETER NAME EXPECTED.')
        IF(IMPX.GT.0)WRITE(IOUT,1000)TEXT
        IF(NPARM.GT.0)THEN
          JPMAP=LCMGID(IPMAP,'PARAM')
          DO IPAR=1,NPARM
            KPMAP=LCMGIL(JPMAP,IPAR)
            CALL LCMGET(KPMAP,'P-NAME',INAME)
            WRITE(PNAME,'(3A4)') (INAME(I),I=1,3)
            IF(PNAME.EQ.TEXT)CALL XABORT('@RESPAR: THE '
     1       //'PARAMETER '//TEXT//' ALREADY EXISTS.')
          ENDDO
        ENDIF
        NPARM=NPARM+1
        JPMAP=LCMLID(IPMAP,'PARAM',NPARM)
        KPMAP=LCMDIL(JPMAP,NPARM)
        READ(TEXT,'(3A4)') (INAME(I),I=1,3)
        CALL LCMPUT(KPMAP,'P-NAME',3,3,INAME)
*       READ PARKEY NAME
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(TEXT.NE.'PARKEY')CALL XABORT('@RESPAR: KEY'
     1   //'WORD PARKEY EXPECTED.')
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESPAR: CHARACTER'
     1   //' DATA FOR PARKEY NAME EXPECTED.')
        READ(TEXT,'(3A4)') (INAME(I),I=1,3)
        CALL LCMPUT(KPMAP,'PARKEY',3,3,INAME)
*       READ PARAMETER TYPE
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(TEXT.EQ.'GLOBAL')THEN
          IPTYP=1
        ELSEIF(TEXT.EQ.'LOCAL')THEN
          IPTYP=2
        ELSE
          CALL XABORT('@RESPAR: INVALID KEYWORD '//TEXT)
        ENDIF
        CALL LCMPUT(KPMAP,'P-TYPE',1,1,IPTYP)
        ISTATE(8)=NPARM
*----
*  SET PARAMETER VALUES
*----
      ELSEIF(TEXT12.EQ.'SET-PARAM')THEN
        IF(NPARM.EQ.0)CALL XABORT('@RESPAR: PARAM'
     1   //'ETER NOT YET DEFINED NPARM=0')
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESPAR: CHARACT'
     1   //'ER DATA FOR PARAMETER NAME EXPECTED.')
*       RECOVER PARAMETER
        JPMAP=LCMGID(IPMAP,'PARAM')
        DO IPAR=1,NPARM
          KPMAP=LCMGIL(JPMAP,IPAR)
          CALL LCMGET(KPMAP,'P-NAME',INAME)
          WRITE(PNAME,'(3A4)') (INAME(I),I=1,3)
          IF(PNAME.EQ.TEXT)THEN
            CALL LCMGET(KPMAP,'P-TYPE',IPTYP)
            GOTO 30
          ENDIF
        ENDDO
        CALL XABORT('@RESPAR: UNABLE TO FIND PARAME'
     1   //'TER WITH PNAME '//TEXT)
   20   IF(IMPX.GT.0)WRITE(IOUT,1001)TEXT
   30   IF(IPTYP.EQ.1)THEN
*         GLOBAL PARAMETER
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
          IF((ITYP.EQ.3).AND.(TEXT.EQ.'OLDMAP')) THEN
            IPTYP=11
            GOTO 20
          ENDIF
          IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL'
     1     //' DATA or OLDMAP keyword FOR VALUE EXPECTED.')
          CALL LCMPUT(KPMAP,'P-VALUE',1,2,FLOT)
        ELSE
*         LOCAL PARAMETER
          CALL XDRSET(VALUE,NCH*NB,0.)
          IF(IPTYP.NE.11) CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
          IF(TEXT.EQ.'SAME')THEN
            CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
            DO ICH=1,NCH
              DO IB=1,NB
                IF(FMIX(ICH,IB).NE.0) THEN
                  IF(ITYP.EQ.2)VALUE(ICH,IB)=FLOT
                  IF(ITYP.EQ.3)CVALUE(ICH,IB)=TEXT
                ENDIF
              ENDDO
            ENDDO
*
          ELSEIF(TEXT.EQ.'CHAN')THEN
            DO ICH=1,NCH
              CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
              IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL'
     1         //' DATA FOR VALUE EXPECTED.')
              DO 40 IB=1,NB
              IF(FMIX(ICH,IB).NE.0) VALUE(ICH,IB)=FLOT
   40         CONTINUE
            ENDDO
*
          ELSEIF(TEXT.EQ.'BUND')THEN
            DO 55 IB=1,NB
            DO 50 ICH=1,NCH
              IF(FMIX(ICH,IB).EQ.0) GO TO 50
              CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
              IF(ITYP.EQ.2)VALUE(ICH,IB)=FLOT
              IF(ITYP.EQ.3)CVALUE(ICH,IB)=TEXT
   50       CONTINUE
   55       CONTINUE
          ELSEIF(TEXT.EQ.'TIMES')THEN
! try to find the parameters called DMOD
            CALL REDGET(ITYP,NITMA,FLOT,KEYN,DFLOT)
            IF(ITYP.NE.3)CALL XABORT('@RESPAR: CHARACTER'
     1       //' DATA FOR VALUE EXPECTED.')
            JPMAP=LCMGID(IPMAP,'PARAM')
            DO IPAR=1,NPARM
              ZPMAP=LCMGIL(JPMAP,IPAR)
              CALL LCMGET(ZPMAP,'P-NAME',INAME)
              WRITE(PNAME,'(3A4)') (INAME(I),I=1,3)
              IF(PNAME.EQ.KEYN)THEN
                CALL LCMGET(ZPMAP,'P-TYPE',IPTYP)
                GOTO 60
              ENDIF
            ENDDO
            CALL XABORT('@RESPAR: UNABLE TO FIND PARAME'
     1        //'TER WITH PNAME '//KEYN)
   60       CONTINUE
            ALLOCATE(DENSMOD(NCH*NB))
            CALL LCMGET(ZPMAP,'P-VALUE',DENSMOD)
            CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
            IF(ITYP.NE.3.OR.TEXT.NE.'SAME')CALL XABORT('@RESPAR:'
     1       //' KEYWORD SAME EXPECTED.')
            CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
            IF(ITYP.NE.2) CALL XABORT('@RESPAR: REAL DATA EXPECTED.')
            DO IB=1,NB
              DO ICH=1,NCH
                IF(FMIX(ICH,IB).NE.0) THEN
                  VALUE(ICH,IB)=FLOT*DENSMOD(ICH+(IB-1)*NCH)
                ENDIF
              ENDDO
            ENDDO
            DEALLOCATE(DENSMOD)
*  R. Chambon - begin
          ELSEIF(TEXT.EQ.'OLDMAP')THEN
            IF(.NOT.LMAP2) CALL XABORT('@RESPAR: SECOND'
     1        //' L_MAP EXPECTED.')
            CALL LCMGET(IPMP2,'STATE-VECTOR',ISTAT2)
            NPARM2=ISTAT2(8)
            JPMP2=LCMGID(IPMP2,'PARAM')
            DO IPAR=1,NPARM2
              KPMP2=LCMGIL(JPMP2,IPAR)
              CALL LCMGET(KPMP2,'P-NAME',INAME)
              WRITE(PNAME2,'(3A4)') (INAME(I),I=1,3)
              IF(PNAME.EQ.PNAME2)THEN
                GOTO 70
              ENDIF
            ENDDO
            CALL XABORT('@RESPAR: UNABLE TO FIND PARAME'
     1       //'TER WITH PNAME in second L_MAP '//TEXT)
   70       CALL LCMLEN(KPMP2,'P-VALUE',NITMA,INDIC)
            IF(NITMA.EQ.0) CALL XABORT('@RESPAR: Record BURN-INST in '
     1        //'SECOND L_MAP EXPECTED.')
            ALLOCATE(VAL2(NITMA))
            CALL LCMGET(KPMP2,'P-VALUE',VAL2)
*         global parameter
            IF(NITMA.EQ.1) THEN
              VALUE(1,1)=VAL2(1)
*         recovered from previous calculation with the same geometry
*         but not the same initialization part
*         example: homogeneous calculation followed by a pin power
*         reconstruction
            ELSEIF(NITMA.EQ.NCH*NB) THEN
              DO ICH=1,NCH
                DO IB=1,NB
                  I=ICH+(IB-1)*NCH
                  VALUE(ICH,IB)=VAL2(I)
                ENDDO
              ENDDO
*         recovered from previous calculation with a different geometry
*         the second geometry must correspond to the assembly geometry
*         of the new geometry
*         examples: homogeneous calculation followed by a heterogeneous
*                   calculation
*                   homogeneous calculation followed by a pin power
*                   calculation
            ELSEIF(NITMA.EQ.NASB*NB) THEN
              CALL LCMGET(IPMAP,'A-ZONE',IAZ)
              DO ICH=1,NCH
                DO IB=1,NB
                  VALUE(ICH,IB)=VAL2(IAZ(ICH)+(IB-1)*NCH)
                ENDDO
              ENDDO
            ENDIF
            DEALLOCATE(VAL2)
*  R. Chambon - End
          ELSEIF(TEXT.EQ.'LEVEL')THEN
*           move a control rod over each channel
            ITOP=1
   75       CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
            IF(ITYP.NE.3)CALL XABORT('@RESPAR: CHARACTER DATA H+, H-,'
     1       //'SAME OR CHAN EXPECTED.')
            IF(TEXT.EQ.'H+')THEN
*             PWR-type moving rod
              ITOP=1
              GO TO 75
            ELSEIF(TEXT.EQ.'H-')THEN
*             BWR-type moving rod
              ITOP=-1
              GO TO 75
            ELSEIF(TEXT.EQ.'SAME') THEN
              CALL REDGET(ITYP,NITMA,ZLEVEL,TEXT,DFLOT)
              IF(ITYP.NE.2) CALL XABORT('@RESPAR: REAL DATA EXPECTED.')
            ENDIF
            JPMAP=LCMGID(IPMAP,'GEOMAP')
            CALL LCMGET(JPMAP,'STATE-VECTOR',ISTAT2)
            NX=ISTAT2(3)
            NY=ISTAT2(4)
            NZ=ISTAT2(5)
            NEL=ISTAT2(6)
            IF((ISTAT2(1).EQ.9).AND.(NY.EQ.0)) NY=1
            ALLOCATE(ZZ(NZ+1),NUM(NEL),IND(NZ),VB(NB))
            CALL LCMGET(JPMAP,'MESHZ',ZZ)
            CALL LCMGET(IPMAP,'BMIX',NUM)
            ICH=0
            DO 105 IY=1,NY
            DO 100 IX=1,NX
            IEL=(IY-1)*NX+IX
            DO 80 IZ=1,NZ
            IF(NUM((IZ-1)*NX*NY+IEL).NE.0) GO TO 90
   80       CONTINUE
            GO TO 100
   90       ICH=ICH+1
            IF(TEXT.EQ.'CHAN') THEN
              CALL REDGET(ITYP,NITMA,ZLEVEL,TEXT,DFLOT)
              IF(ITYP.NE.2) CALL XABORT('@RESPAR: REAL DATA EXPECTED.')
            ENDIF
            IF((ZLEVEL.LT.0.0).OR.(ZLEVEL.GT.1.0)) THEN
              CALL XABORT('@RESPAR: 0<=LEVEL<=1 EXPECTED.')
            ENDIF
            IB=0
            DO IZ=1,NZ
              IND(IZ)=0
              IF(NUM((IZ-1)*NX*NY+IEL).EQ.0) CYCLE
              IB=IB+1
              IND(IZ)=IB
            ENDDO
            IF(IB.NE.NB) CALL XABORT('@RESPAR: INVALID NUMBER OF BUNDL'
     1      //'ES.')
            CALL RESROD(NB,NZ,ZZ,IND,ZLEVEL,ITOP,VB)
            DO IB=1,NB
              IF(FMIX(ICH,IB).NE.0) VALUE(ICH,IB)=VB(IB)
            ENDDO
  100       CONTINUE
  105       CONTINUE
            IF(ICH.NE.NCH) CALL XABORT('@RESPAR: INVALID NUMBER OF CHA'
     1      //'NNELS.')
            DEALLOCATE(VB,IND,NUM,ZZ)
          ELSE
            CALL XABORT('@RESPAR: INVALID KEYWORD '//TEXT)
          ENDIF
          IF(ITYP.EQ.2)CALL LCMPUT(KPMAP,'P-VALUE',NCH*NB,2,VALUE)
          IF(ITYP.EQ.3)CALL LCMPTC(KPMAP,'P-VALUE',12,NCH*NB,CVALUE)
          IF(ITYP.EQ.11)CALL LCMPUT(KPMAP,'P-VALUE',1,2,VALUE(1,1))
        ENDIF
*----
*  CHANNEL REFUELLING SCHEMES
*----
      ELSEIF(TEXT12.EQ.'REF-SHIFT')THEN
        IF(IMPX.GT.0)WRITE(IOUT,1002)
*       BUNDLE-SHIFT NUMBERS
        CALL XDISET(IBSH,NCOMB,0)
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.EQ.1)THEN
          IF(NITMA.LE.0.OR.NITMA.GT.NB)CALL XABORT('@RESPAR:'
     1     //' BUNDLE-SHIFT MUST BE POSITIVE AND NON-ZERO NUMBER'
     1     //' AND MAX EQUAL TO NUMBER OF FUEL BUNDLES PER CHANNEL')
          DO 110 ICZ=1,NCOMB
          IBSH(ICZ)=NITMA
  110     CONTINUE
        ELSEIF((ITYP.EQ.3).AND.(TEXT.EQ.'COMB'))THEN
          DO ICZ=1,NCOMB
            CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
            IF(ITYP.NE.1)CALL XABORT('@RESPAR: INTEGER BUNDLE'
     1      //'-SHIFT NUMBER PER COMBUSTION-ZONE EXPECTED.')
            IF(NITMA.LE.0.OR.NITMA.GT.NB)CALL XABORT('@RESPAR:'
     1      //' BUNDLE-SHIFT MUST BE POSITIVE AND NON-ZERO NUMBER.'
     1      //' AND MAX EQUAL TO NUMBER OF FUEL BUNDLES PER CHANNEL')
            IBSH(ICZ)=NITMA
          ENDDO
        ELSE
          CALL XABORT('@RESPAR: INVALID INPUT FOR REF-SHIFT.')
        ENDIF
        CALL LCMPUT(IPMAP,'REF-SHIFT',NCOMB,1,IBSH)
*       REFUELLING VECTOR
        CALL XDISET(IVECT,NCOMB*NB,0)
        DO 120 ICZ=1,NCOMB
        ISHIFT=IBSH(ICZ)
        IF(ISHIFT.EQ.NB)GOTO 120
        NREF=NB-ISHIFT
        DO IREF=1,NREF
          IPOS=ISHIFT+IREF
          IVECT(ICZ,IPOS)=IREF
        ENDDO
  120   CONTINUE
        CALL LCMPUT(IPMAP,'REF-VECTOR',NCOMB*NB,1,IVECT)
*       CHANNEL REFUELLING SCHEMES
        CALL LCMGET(IPMAP,'B-ZONE',IZONE)
        CALL LCMGET(IPMAP,'BMIX',MIX)
        CALL XDISET(NSCH,NCH,0)
        IEL=0
        ICH=0
        DO 135 IY=1,NY
        DO 130 IX=1,NX
        IEL=IEL+1
        IF(MIX(IEL).EQ.0)GOTO 130
        ICH=ICH+1
        ISHIFT=IBSH(IZONE(ICH))
        NSCH(ICH)=((-1)**(IEL+IY-1))*ISHIFT
  130   CONTINUE
  135   CONTINUE
        CALL LCMPUT(IPMAP,'REF-SCHEME',NCH,1,NSCH)
        LRSCH=.TRUE.
*----
*  BURNUP INTERPOLATION TYPE
*----
      ELSEIF(TEXT12.EQ.'BTYPE')THEN
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESPAR: BURN'
     1   //'UP INTERPOLATION OPTION EXPECTED.')
        IBTYP=0
        IF(TEXT.EQ.'TIMAV-BURN')THEN
          IBTYP=1
        ELSEIF(TEXT.EQ.'INST-BURN')THEN
          IBTYP=2
        ELSE
          CALL XABORT('@RESPAR: INVALID INPUT FOR BTYPE.')
        ENDIF
        ISTATE(5)=IBTYP
*----
*  AVERAGE EXIT BURNUPS
*----
      ELSEIF(TEXT12.EQ.'TIMAV-BVAL')THEN
        IF(IMPX.GT.0)WRITE(IOUT,1003)
        ALLOCATE(BRN(NCOMB))
        CALL XDRSET(BRN,NCOMB,0.)
        DO ICZ=1,NCOMB
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL DATA'
     1     //' FOR BURNUP VALUE EXPECTED(1).')
          IF(FLOT.LE.0.)CALL XABORT('@RESPAR: INVALID'
     1     //' DATA FOR AVERAGE BURNUP VALUE =0.')
          BRN(ICZ)=FLOT
        ENDDO
        CALL LCMPUT(IPMAP,'BURN-AVG',NCOMB,2,BRN)
        DEALLOCATE(BRN)
        LBURN=.TRUE.
*----
*  INSTANTANEOUS BURNUPS
*----
      ELSEIF(TEXT12.EQ.'INST-BVAL')THEN
        IF(IMPX.GT.0)WRITE(IOUT,1004)
        CALL XDRSET(VALUE,NCH*NB,0.)
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESPAR: KEYWORD'
     1   //' SAME|CHAN|BUND EXPECTED (1).')
        IF(TEXT.EQ.'BUND')THEN
          DO IB=1,NB
            DO ICH=1,NCH
              IF(FMIX(ICH,IB).NE.0) THEN
                CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
                IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL DATA'
     1            //' FOR BURNUP VALUE EXPECTED(2).')
                IF(FLOT.LT.0.)CALL XABORT('@RESPAR: INVALID DA'
     1            //'TA FOR BURNUP VALUE <0.')
                VALUE(ICH,IB)=FLOT
              ENDIF
            ENDDO
          ENDDO
        ELSEIF(TEXT.EQ.'CHAN')THEN
          DO ICH=1,NCH
            CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
            IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL DATA'
     1        //' FOR BURNUP VALUE EXPECTED(2).')
            IF(FLOT.LT.0.)CALL XABORT('@RESPAR: INVALID DA'
     1        //'TA FOR BURNUP VALUE <0.')
            DO IB=1,NB
              IF(FMIX(ICH,IB).NE.0) VALUE(ICH,IB)=FLOT
            ENDDO
          ENDDO
        ELSEIF(TEXT.EQ.'ASBL')THEN
          IF(NASB.EQ.0)CALL XABORT('@RESPAR: ASSEMBLY'
     1        //' NOT DEFINED.')
          ALLOCATE(BASS(NASB))
          DO IASS=1,NASB
            CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
            IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL DATA'
     1        //' FOR BURNUP VALUE EXPECTED(2).')
            IF(FLOT.LT.0.)CALL XABORT('@RESPAR: INVALID DA'
     1        //'TA FOR BURNUP VALUE <0.')
            BASS(IASS)=FLOT
          ENDDO
          CALL LCMGET(IPMAP,'A-ZONE',IAZ)
          DO ICH=1,NCH
            DO IB=1,NB
              VALUE(ICH,IB)=BASS(IAZ(ICH))
            ENDDO
          ENDDO
          DEALLOCATE(BASS)
        ELSEIF(TEXT.EQ.'SAME')THEN
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL DATA'
     1      //' FOR BURNUP VALUE EXPECTED(2).')
          IF(FLOT.LT.0.)CALL XABORT('@RESPAR: INVALID DA'
     1      //'TA FOR BURNUP VALUE <0.')
          DO ICH=1,NCH
            DO IB=1,NB
              IF(FMIX(ICH,IB).NE.0) VALUE(ICH,IB)=FLOT
            ENDDO
          ENDDO
*  R. Chambon - begin
        ELSEIF(TEXT.EQ.'OLDMAP')THEN
          IF(.NOT.LMAP2) CALL XABORT('@RESPAR: SECOND'
     1      //' L_MAP EXPECTED.')
          CALL LCMLEN(IPMP2,'BURN-INST',NITMA,INDIC)
          IF(NITMA.EQ.0) CALL XABORT('@RESPAR: Record BURN-INST in '
     1      //'SECOND L_MAP EXPECTED.')
          ALLOCATE(VAL2(NITMA))
          CALL LCMGET(IPMP2,'BURN-INST',VAL2)
*         recovered from previous calculation with the same geometry but
*         not the same initialization part
*         example: homogeneous calculation followed by a pin power
*         reconstruction
          IF(NITMA.EQ.NCH*NB) THEN
            DO ICH=1,NCH
              DO IB=1,NB
                I=ICH+(IB-1)*NCH
                VALUE(ICH,IB)=VAL2(I)
              ENDDO
            ENDDO
*         recovered from previous calculation with a different geometry
*         the second geometry must correspond to the assembly geometry
*         of the new geometry
*         examples: homogeneous calculation followed by a heterogeneous
*                   calculation
*                   homogeneous calculation followed by a pin power
*                   calculation
          ELSEIF(NITMA.EQ.NASB*NB) THEN
            CALL LCMGET(IPMAP,'A-ZONE',IAZ)
            DO ICH=1,NCH
              DO IB=1,NB
                VALUE(ICH,IB)=VAL2(IAZ(ICH)+(IB-1)*NCH)
              ENDDO
            ENDDO
          ENDIF
          DEALLOCATE(VAL2)
*  R. Chambon - End
        ELSEIF(TEXT.EQ.'SMOOTH')THEN
* EACH 'BURN-INST' WILL HAVE THE SAME BURNUP AS THEIR FIRST INDEX IN 'FLMIX'
          CALL LCMGET(IPMAP,'BURN-INST',VALUE)
          DO ICH=1,NCH
            DO IB=1,NB
              JBKEEP=0
              DO JCH=1,NCH
                DO JB=1,NB
* FIRST INDEX OF FMIX(ICH,IB) IS AT JCH,JB
                  JBKEEP=JB
                  IF(FMIX(ICH,IB).EQ.FMIX(JCH,JB)) GOTO 140
                ENDDO
              ENDDO
              CALL XABORT('@RESPAR: ASSERTION ERROR (NO FIRST INDEX)')
  140         VALUE(ICH,IB)=VALUE(JCH,JBKEEP)
            ENDDO
          ENDDO
        ELSE
          CALL XABORT('@RESPAR: KEYWORD'
     1      //' SAME|CHAN|BUND|ASBL|OLDMAP|SMOOTH EXPECTED (2).')
        ENDIF
        CALL LCMPUT(IPMAP,'BURN-INST',NCH*NB,2,VALUE)
*----
*  BUNDLE POWERS IN KW
*----
      ELSEIF(TEXT12.EQ.'BUNDLE-POW')THEN
        IF(IMPX.GT.0)WRITE(IOUT,1006)
        CALL XDRSET(POWER,NCH*NB,0.)
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESPAR: KEYWORD'
     1   //' BUND|CHAN|SAME EXPECTED (3).')
        IF(TEXT.EQ.'BUND')THEN
          DO IB=1,NB
            DO ICH=1,NCH
              IF(FMIX(ICH,IB).NE.0) THEN
                CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
                IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL DATA'
     1            //' FOR POWER VALUE EXPECTED(1).')
                IF(FLOT.LT.0.)CALL XABORT('@RESPAR: INVALID DA'
     1            //'TA FOR POWER VALUE <0.')
                POWER(ICH,IB)=FLOT
              ENDIF
            ENDDO
          ENDDO
        ELSEIF(TEXT.EQ.'CHAN')THEN
          DO ICH=1,NCH
            CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
            IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL DATA'
     1        //' FOR POWER VALUE EXPECTED(2).')
            IF(FLOT.LT.0.)CALL XABORT('@RESPAR: INVALID DA'
     1        //'TA FOR POWER VALUE <0.')
            DO IB=1,NB
              IF(FMIX(ICH,IB).NE.0) POWER(ICH,IB)=FLOT
            ENDDO
          ENDDO
        ELSEIF(TEXT.EQ.'SAME')THEN
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL DATA'
     1      //' FOR BURNUP VALUE EXPECTED(2).')
          IF(FLOT.LT.0.)CALL XABORT('@RESPAR: INVALID DA'
     1      //'TA FOR POWER VALUE <0.')
          DO ICH=1,NCH
            DO IB=1,NB
              IF(FMIX(ICH,IB).NE.0) POWER(ICH,IB)=FLOT
            ENDDO
          ENDDO
        ELSE
          CALL XABORT('@RESPAR: KEYWORD SAME|CHAN|BUND EXPECTED (4).')
        ENDIF
        CALL LCMPUT(IPMAP,'BUND-PW',NCH*NB,2,POWER)
        PTOT=0.0
        DO ICH=1,NCH
          DO IB=1,NB
            PTOT=PTOT+POWER(ICH,IB)
          ENDDO
        ENDDO
        PTOT=PTOT/1.0E3
        CALL LCMPUT(IPMAP,'REACTOR-PW',1,2,PTOT)
*----
*  AXIAL POWERS FORM FACTORS
*----
      ELSEIF(TEXT12.EQ.'AXIAL-PFORM')THEN
        IF(IMPX.GT.0)WRITE(IOUT,1007)
        IF(PTOT.EQ.0.0)CALL XABORT('@RESPAR: FULL REACTOR POWER NOT S'
     1  //'ET.')
        CALL XDRSET(FPOWER,NB,0.)
        DO IB=1,NB
          CALL REDGET(ITYP,NITMA,FPOWER(IB),TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL DATA FOR POWERS FOR'
     1      //'M FACTORS VALUE EXPECTED.')
          IF(FPOWER(IB).LT.0.)CALL XABORT('@RESPAR: INVALID DATA FOR '
     1      //'POWERS FORM FACTORS VALUE <0.')
        ENDDO
        CALL LCMPUT(IPMAP,'AXIAL-FPW',NB,2,FPOWER)
        DSUM=0.0
        DO IB=1,NB
          DSUM=DSUM+FPOWER(IB)
        ENDDO
        DO ICH=1,NCH
          DO IB=1,NB
            POWER(ICH,IB)=FPOWER(IB)*PTOT*1.0E3/(DSUM*REAL(NCH))
          ENDDO
        ENDDO
        CALL LCMPUT(IPMAP,'BUND-PW',NCH*NB,2,POWER)
*----
*  FULL REACTOR POWER IN MW
*----
      ELSEIF(TEXT12.EQ.'REACTOR-POW')THEN
        IF(IMPX.GT.0)WRITE(IOUT,1008)
        CALL REDGET(ITYP,NITMA,PTOT,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL DATA'
     1    //' FOR FULL REACTOR POWER VALUE EXPECTED.')
        IF(PTOT.LT.0.)CALL XABORT('@RESPAR: INVALID DA'
     1    //'TA FOR FULL REACTOR POWER VALUE <0.')
        CALL LCMPUT(IPMAP,'REACTOR-PW',1,2,PTOT)
*----
*  FUEL-TYPE DATA
*----
      ELSEIF(TEXT12.EQ.'FUEL')THEN
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESPAR: KEYWO'
     1   //'RD FOR FUEL-TYPE PARAMETER EXPECTED.')
        IF((TEXT.NE.'WEIGHT').AND.(TEXT.NE.'ENRICH').AND.
     1   (TEXT.NE.'POISON'))CALL XABORT('@RESPAR: INVAL'
     2   //'ID INPUT FOR FUEL.')
        IF(IMPX.GT.0)WRITE(IOUT,1005)TEXT
        JPMAP=LCMLID(IPMAP,'FUEL',NFUEL)
        DO IFUEL=1,NFUEL
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@RESPAR: REAL'
     1     //' DATA PER EACH FUEL-TYPE EXPECTED.')
          KPMAP=LCMDIL(JPMAP,IFUEL)
          CALL LCMPUT(KPMAP,TEXT,1,2,FLOT)
          IF(IMPX.GT.0)CALL LCMLIB(KPMAP)
        ENDDO
      ELSEIF(TEXT12.EQ.'CELL')THEN
        ALLOCATE(ALCH(NCH))
        DO 150 I=1,NCH
          CALL REDGET(INDIC,ALCH(I),FLOTT,TEXT12,DFLOT)
          IF(INDIC.NE.1) CALL XABORT('@RESPAR: INTEGER DATA EXPECTED.')
  150   CONTINUE
        CALL RESCEL(IPMAP,NCH,NB,ALCH)
        DEALLOCATE(ALCH)
        ISTATE(5)=2
      ELSE
        CALL XABORT('@RESPAR: WRONG KEYWORD '//TEXT12)
      ENDIF
      GOTO 10
  500 IF(LRSCH.OR.LBURN)CALL RESBRN(IPMAP,NCH,NB,NCOMB,
     1 NX,NY,NZ,LRSCH,IMPX)
      RETURN
*
 1000 FORMAT(/1X,'INPUT OF NEW PARAMETER: ',A12)
 1001 FORMAT(/1X,'READING VALUES FOR PARAMETER: ',A12)
 1002 FORMAT(/1X,'READING INPUT FOR REF-SHIFT')
 1003 FORMAT(/1X,'READING AVERAGE EXIT BURNUPS')
 1004 FORMAT(/1X,'READING INSTANTANEOUS BURNUPS')
 1005 FORMAT(/1X,'READING DATA FOR FUEL-TYPE PARAMETER: ',A12)
 1006 FORMAT(/1X,'READING BUNDLE POWERS IN KW')
 1007 FORMAT(/1X,'READING BUNDLE POWERS FORM FACTORS')
 1008 FORMAT(/1X,'READING FULL REACTOR POWER IN MW')
      END
