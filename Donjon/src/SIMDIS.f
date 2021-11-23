*DECK SIMDIS
      SUBROUTINE SIMDIS(LSET,NCH,NB,HCYC,NASMB1,ASMB1,FORM,ASMB1B,ZONE,
     > BURNUP,OBURNU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Axial normalization of the burnup distribution using information from
* another cycle.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input/output
* LSET    type of normalization (=.true.: use FORM info; =.false: use
*         an existing assembly).
* NCH     number of assemblies or number of quart-of-assemblies.
* NB      number of axial burnup subdivisions in an assembly.
* HCYC    name of cycle.
* NASMB1  number of assemblies to set.
* ASMB1   group of assembly names, as defined in the fuel map, to set
*         at specific burnup.
* FORM    axial form factor used if LSET=.true.
* ASMB1B  assembly name, as defined in the fuel map, to which we
*         want to use the burnup distribution if LSET=.false.
* ZONE    default assembly or quart-of-assembly names as defined in
*         the fuel map.
* BURNUP  burnups during a refuelling cycle.
* OBURNU  burnups during a previous refuelling cycle.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL LSET
      INTEGER NCH,NB,NASMB1
      CHARACTER HCYC*12,ASMB1(NASMB1)*4,ASMB1B*4,ZONE(NCH)*4
      REAL FORM(NB),BURNUP(NCH,NB),OBURNU(NCH,NB)
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
      ZNUM=0.0
      ZDEN=0.0
      DO IASMB1=1,NASMB1
        DO 10 ICH=1,NCH
        ZONE2(ICH)=ZONE(ICH)
   10   CONTINUE
        DO ICH=1,NCH
          IF(ZONE(ICH).EQ.ASMB1(IASMB1)) THEN
            IF(LSET) THEN
              ZNUM=0.0
              ZDEN=0.0
              DO IB=1,NB
                ZNUM=ZNUM+BURNUP(ICH,IB)
                ZDEN=ZDEN+FORM(IB)
              ENDDO
              DO IB=1,NB
                BURNUP(ICH,IB)=FORM(IB)*ZNUM/ZDEN
              ENDDO
            ELSE
              IOLD=0
              DO ICH2=1,NCH
                IF(ZONE2(ICH2).EQ.ASMB1B) THEN
                  IOLD=ICH2
                  ZONE2(ICH2)=' '
                  GO TO 20
                ENDIF
              ENDDO
              WRITE(HSMG,'(33H@SIMDIS: UNABLE TO FIND ASSEMBLY ,A4,
     >        25HIN THE FUEL MAP AT CYCLE ,A12,1H.)') ASMB1(IASMB1),
     >        HCYC
              CALL XABORT(HSMG)
   20         ZNUM=0.0
              ZDEN=0.0
              DO IB=1,NB
                ZNUM=ZNUM+BURNUP(ICH,IB)
                ZDEN=ZDEN+OBURNU(IOLD,IB)
              ENDDO
              DO IB=1,NB
                BURNUP(ICH,IB)=OBURNU(IOLD,IB)*ZNUM/ZDEN
              ENDDO
            ENDIF
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
