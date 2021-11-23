*DECK SIMCPY
      SUBROUTINE SIMCPY(NCH,NB,HCYC,NASMB1,ASMB1,ASMB1B,ZONE,NIS,NAME,
     > BURNUP,FMIX,RFOLLO,ONAME,OBURNU,OFMIX,OFOLLO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the burnup of an assembly in another cycle.
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
*         at specific burnup.
* ASMB1B  assembly name, as defined in the fuel map, to which we
*         want to copy burnup.
* ZONE    default assembly or quart-of-assembly names as defined in
*         the fuel map.
* NIS     number of particularized isotopes.
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
* OBURNU  burnups at the end of a previous refuelling cycle.
* OFMIX   assembly types in a previous refuelling cycle.
* OFOLLO  number densities of the particularized isotopes at the end
*         of a previous refuelling cycle.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NCH,NB,NASMB1,NIS,FMIX(NCH,NB),OFMIX(NCH,NB)
      CHARACTER HCYC*12,ASMB1(NASMB1)*4,ASMB1B*4,ZONE(NCH)*4,
     > NAME(NCH)*12,ONAME(NCH)*12
      REAL BURNUP(NCH,NB),RFOLLO(NCH,NB,NIS),OBURNU(NCH,NB),
     > OFOLLO(NCH,NB,NIS)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: ZONE2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ZONE2(NCH))
*
      DO IASMB1=1,NASMB1
        DO 10 ICH=1,NCH
        ZONE2(ICH)=ZONE(ICH)
   10   CONTINUE
        DO ICH=1,NCH
          IF(ZONE(ICH).EQ.ASMB1(IASMB1)) THEN
            IOLD=0
            DO ICH2=1,NCH
              IF(ZONE2(ICH2).EQ.ASMB1B) THEN
                IOLD=ICH2
                ZONE2(ICH2)=' '
                GO TO 20
              ENDIF
            ENDDO
            WRITE(HSMG,'(33H@SIMCPY: UNABLE TO FIND ASSEMBLY ,A4,
     >      25HIN THE FUEL MAP AT CYCLE ,A12,1H.)') ASMB1B,HCYC
            CALL XABORT(HSMG)
   20       DO IB=1,NB
              IF(BURNUP(ICH,IB).NE.-999.0) THEN
                WRITE(HSMG,'(38H@SIMCPY: BURNUP ALREADY DEFINED IN CHA,
     >          4HNNEL,I4,10HAND BUNDLE,I4,10H AT CYCLE ,A12,1H.)')
     >          ICH,IB,HCYC
              ENDIF
              BURNUP(ICH,IB)=OBURNU(IOLD,IB)
              FMIX(ICH,IB)=OFMIX(ICH,IB)
              DO ISO=1,NIS
                RFOLLO(ICH,IB,ISO)=OFOLLO(IOLD,IB,ISO)
              ENDDO
            ENDDO
            NAME(ICH)=ONAME(IOLD)
            CYCLE
          ENDIF
        ENDDO
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ZONE2)
      RETURN
      END
