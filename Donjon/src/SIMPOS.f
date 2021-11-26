*DECK SIMPOS
      SUBROUTINE SIMPOS(LX,LY,NCH,NB,HCYC,HOLD,HHX,IHY,ZONE,INFMIX,
     > NIS,CYCLE,NAME,BURNUP,FMIX,RFOLLO,ONAME,OBURNU,OFMIX,OFOLLO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set the correspondance between assembly indices during refuelling.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input/output
* LX      number of assemblies along the X axis.
* LY      number of assemblies along the Y axis.
* NCH     number of assemblies or number of quart-of-assemblies.
* NB      number of axial burnup subdivisions in an assembly.
* HCYC    name of cycle.
* HOLD    name of previous cycle.
* HHX     naval battle indices along X axis.
* IHY     naval battle indices along Y axis.
* ZONE    default assembly or quart-of-assembly names as defined in
*         the fuel map.
* INFMIX  assembly types as defined in the fuel map.
* NIS     number of particularized isotopes.
* CYCLE   shuffling matrix for refuelling as provided by the plant
*         operator. The name "|" is reserved for empty locations. 
* NAME    names of each assembly or of each quart-of assembly during
*         a refuelling cycle. All quart-of-assembly belonging to the
*         same assembly have the same name.
* BURNUP  burnups during a refuelling cycle. A value of -999.0 means
*         a non-initialized value.
* FMIX    assembly mixtures after refuelling.
* RFOLLO  number densities of the particularized isotopes after
*         refuelling.
* ONAME   names of each assembly or of each quart-of assembly during
*         a previous refuelling cycle.
* OBURNU  burnups during a previous refuelling cycle.
* OFMIX   assembly types in a previous refuelling cycle.
* OFOLLO  number densities of the particularized isotopes at the end
*         of a previous refuelling cycle.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LX,LY,NCH,NB,IHY(LY),INFMIX(NCH,NB),NIS,FMIX(NCH,NB),
     > OFMIX(NCH,NB)
      CHARACTER HCYC*12,HOLD*12,HHX(LX)*1,ZONE(NCH)*4,CYCLE(LX,LY)*4,
     > NAME(NCH)*12,ONAME(NCH)*12
      REAL BURNUP(NCH,NB),RFOLLO(NCH,NB,NIS),OBURNU(NCH,NB),
     > OFOLLO(NCH,NB,NIS)
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT4*4,TEXT1*1,HSMG*131
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: ZONE2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ZONE2(NCH))
*
      MAXINF=0
      DO 10 ICH=1,NCH
      MAXINF=MAX(MAXINF,MAXVAL(INFMIX(ICH,:NB)))
      ZONE2(ICH)=ZONE(ICH)
   10 CONTINUE
      DO ICH=1,NCH
        TEXT4=ZONE(ICH)
        READ(TEXT4,'(A1,I2)') TEXT1,INTG2
        INDX=0
        DO IX=1,LX
          IF(TEXT1.EQ.HHX(IX)) INDX=IX
        ENDDO
        IF(INDX.EQ.0) CALL XABORT('@SIMPOS: UNABLE TO FIND INDX.')
        INDY=0
        DO IY=1,LY
          IF(INTG2.EQ.IHY(IY)) INDY=IY
        ENDDO
        IF(INDY.EQ.0) CALL XABORT('@SIMPOS: UNABLE TO FIND INDY.')
        TEXT4=CYCLE(INDX,INDY)
        IF((TEXT4.EQ.'|').OR.(TEXT4.EQ.'-').OR.(TEXT4.EQ.'-|-')) THEN
          WRITE(HSMG,'(16H@SIMPOS: CHANNEL,I4,21H REFERS TO LOCATION (,
     >    I4,1H,,I4,37H) WHICH IS OUTSIDE THE CORE AT CYCLE ,A12,1H.)')
     >    ICH,INDX,INDY,HCYC
          CALL XABORT(HSMG)
        ELSE IF(TEXT4.EQ.'SPC') THEN
          DO IB=1,NB
            BURNUP(ICH,IB)=-999.0
            FMIX(ICH,IB)=INFMIX(ICH,IB)
            DO ISO=1,NIS
              RFOLLO(ICH,IB,ISO)=0.0
            ENDDO
          ENDDO
          WRITE(NAME(ICH),'(A3,1H/,A8)') TEXT4(:3),HCYC(:8)
        ELSE IF(TEXT4.EQ.'NEW') THEN
          DO IB=1,NB
            BURNUP(ICH,IB)=0.0
            FMIX(ICH,IB)=INFMIX(ICH,IB)
            DO ISO=1,NIS
              RFOLLO(ICH,IB,ISO)=0.0
            ENDDO
          ENDDO
          WRITE(NAME(ICH),'(A3,1H/,A8)') TEXT4(:3),HCYC(:8)
        ELSE IF(TEXT4(4:).EQ.'@') THEN
          READ(TEXT4,'(I3,1X)') NITMA
          IF(NITMA.GT.MAXINF) CALL XABORT('@SIMPOS: MAXINF OVERFLOW.')
          DO IB=1,NB
            BURNUP(ICH,IB)=0.0
            FMIX(ICH,IB)=INFMIX(ICH,IB)
            IF(INFMIX(ICH,IB).NE.0) FMIX(ICH,IB)=NITMA
            DO ISO=1,NIS
              RFOLLO(ICH,IB,ISO)=0.0
            ENDDO
          ENDDO
          WRITE(NAME(ICH),'(A3,1H/,A8)') 'NEW',HCYC(:8)
        ELSE
          IF(HOLD.EQ.' ') CALL XABORT('@SIMPOS: NO PREVIOUS CYCLE.')
          IOLD=0
          DO ICH2=1,NCH
            IF(ZONE2(ICH2).EQ.TEXT4) THEN
              IOLD=ICH2
              ZONE2(ICH2)=' '
              GO TO 20
            ENDIF
          ENDDO
          WRITE(HSMG,'(33H@SIMPOS: UNABLE TO FIND ASSEMBLY ,A4,
     >    25HIN THE FUEL MAP AT CYCLE ,A12,1H.)') TEXT4,HCYC
          CALL XABORT(HSMG)
   20     DO IB=1,NB
            BURNUP(ICH,IB)=OBURNU(IOLD,IB)
            FMIX(ICH,IB)=OFMIX(IOLD,IB)
            DO ISO=1,NIS
              RFOLLO(ICH,IB,ISO)=OFOLLO(IOLD,IB,ISO)
            ENDDO
          ENDDO
          NAME(ICH)=ONAME(IOLD)
        ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ZONE2)
      RETURN
      END
