*DECK SIMSET
      SUBROUTINE SIMSET(NCH,NB,HCYC,NASMB1,ASMB1,BURN,IFUEL,ZONE,NAME,
     > BURNUP,FMIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set the burnup and fuel type of a group of assemblies at positions
* ASMB1.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input/output
* NCH     number of assemblies or number of quart-of-assemblies.
* NB      number of axial burnup subdivisions in an assembly.
* HCYC    name of cycle.
* NASMB1  number of assemblies to set.
* ASMB1   group of assembly names, as defined in the fuel map, to set
*         at specific burnup or fuel type.
* BURN    burnup in MW-day/tonne. Burnup must be set if .ne.-999.0.
* IFUEL   fuel type. Fuel type must be set if .ne.0.
* ZONE    default assembly or quart-of-assembly names as defined in
*         the fuel map.
* NAME    names of each assembly or of each quart-of assembly during
*         a refuelling cycle. All quart-of-assembly belonging to the
*         same assembly have the same name.
* BURNUP  burnups during a refuelling cycle. A value of -999.0 means
*         a non-initialized value.
* FMIX    assembly types after refuelling.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NCH,NB,NASMB1,FMIX(NCH,NB)
      CHARACTER HCYC*12,ZONE(NCH)*4,ASMB1(NASMB1)*4,NAME(NCH)*12
      REAL BURN,BURNUP(NCH,NB)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131
*
      DO IASMB1=1,NASMB1
        DO ICH=1,NCH
          IF(ASMB1(IASMB1).EQ.ZONE(ICH)) THEN
            DO IB=1,NB
              IF(BURN.NE.-999.0) THEN
                IF(BURNUP(ICH,IB).NE.-999.0) THEN
                  WRITE(HSMG,'(36H@SIMSET: BURNUP ALREADY DEFINED IN C,
     >            6HHANNEL,I4,10HAND BUNDLE,I4,10H AT CYCLE ,A12,1H.)')
     >            ICH,IB,HCYC
                ENDIF
                BURNUP(ICH,IB)=BURN
              ENDIF
              IF(IFUEL.NE.0) THEN
                FMIX(ICH,IB)=IFUEL
              ENDIF
            ENDDO
            NAME(ICH)=ASMB1(IASMB1)(:3)//HCYC(:9)
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END
