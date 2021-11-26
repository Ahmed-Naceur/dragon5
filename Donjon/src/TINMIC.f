*DECK TINMIC
      SUBROUTINE TINMIC(IPMIC,IPMIC2,IPMIC3,NB,NCH,NW,ICH,NISO,NISO2,
     1 IWORK,BSH,NDENS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Microscopic refueling to update the microlib (micro-depletion)
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
*
*Author(s): 
* M. Guyot
*
*Parameters: input/output
* IPMIC   Adress of the L_LIBRARY in creation mode.
* IPMIC2  Adress of the fuel-map L_LIBRARY in read-only mode.
* IPMIC3  Adress of the L_LIBRARY in read-only mode, containing new
*         fuel properties.
* NB      Number of bundles
* NCH     Number of channels
* NW      Vector containing new index for the refuelling
* ICH     Number of the channel to refuel
* NISO    Number of isotopes in the fuel-map microlib
* NISO2   Number of isotopes in the third microlib
* IWORK   Useful vector for refueling
* BSH     Vector containing new mixtures after shifting
* NDENS   New isotopic densities after refueling
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMIC,IPMIC2,IPMIC3
      INTEGER NB,NCH,NW(NB),ICH,NISO,NISO2,IWORK(NB,2),BSH(NB)
      REAL NDENS(NISO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IOUT=6,MAXISO=100)
      INTEGER IB,ISO,ISO2,IMIX,I,SHT(NB),IND(MAXISO),IND2(MAXISO),I1,I2
      CHARACTER TEXT*12,TEXT2*12
      LOGICAL LMIX
      TYPE(C_PTR) JPMIC2,JPMIC3
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX,MIX2,TODO,TODO2,TYP,TYP2
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NAME,NAME2,USED,USED2
      REAL, ALLOCATABLE, DIMENSION(:) :: DENS,DENS2,TEMP,TEMP2,VOL,VOL2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MIX(NISO),MIX2(NISO2),TODO(NISO),TODO2(NISO2),TYP(NISO),
     1 TYP2(NISO2))
      ALLOCATE(NAME(3,NISO),NAME2(3,NISO2),USED(3,NISO),USED2(3,NISO2))
      ALLOCATE(DENS(NISO),DENS2(NISO2),TEMP(NISO),TEMP2(NISO2),
     1 VOL(NISO),VOL2(NISO2))
*----
*  RECOVER INFORMATION
*----
      CALL LCMGET(IPMIC2,'ISOTOPESMIX',MIX)
      CALL LCMGET(IPMIC3,'ISOTOPESMIX',MIX2)
      CALL LCMGET(IPMIC2,'ISOTOPERNAME',NAME)
      CALL LCMGET(IPMIC3,'ISOTOPERNAME',NAME2)
      CALL LCMGET(IPMIC2,'ISOTOPESUSED',USED)
      CALL LCMGET(IPMIC3,'ISOTOPESUSED',USED2)
      CALL LCMGET(IPMIC2,'ISOTOPESDENS',DENS)
      CALL LCMGET(IPMIC3,'ISOTOPESDENS',DENS2)
      CALL LCMGET(IPMIC2,'ISOTOPESTODO',TODO)
      CALL LCMGET(IPMIC3,'ISOTOPESTODO',TODO2)
      CALL LCMGET(IPMIC2,'ISOTOPESTYPE',TYP)
      CALL LCMGET(IPMIC3,'ISOTOPESTYPE',TYP2)
      CALL LCMGET(IPMIC2,'ISOTOPESTEMP',TEMP)
      CALL LCMGET(IPMIC3,'ISOTOPESTEMP',TEMP2)
      CALL LCMGET(IPMIC2,'ISOTOPESVOL',VOL)
      CALL LCMGET(IPMIC3,'ISOTOPESVOL',VOL2)
*----
*  CHECK IF THE MIXTURES TO SHIFT EXIST IN THE MICROLIB
*----
      DO 10 IB=1,NB
        LMIX=.FALSE.
        DO 15 ISO2=1,NISO2
          IF(MIX2(ISO2).EQ.IWORK(IB,2)) THEN
            LMIX=.TRUE.
          ENDIF
 15     CONTINUE
        IMIX=MIX2(ISO2)
        IF(.NOT.LMIX) THEN
          WRITE(IOUT,*) '@TINMIC: THE MIXTURE ',IMIX,' IS NOT PRESENT '
     +      //'IN THE MICROLIB FOR THE REFUEL. '
          CALL XABORT('@TINMIC: REFUELING ERROR. ')
        ENDIF
 10   CONTINUE
*----
*  COMPUTE THE VECTORS FOR THE REFUELING
*----
*     SHT CONTAINS THE MIXTURES OF THE CHANNEL TO SHIFT
      CALL XDISET(SHT,NB,0)
      DO I=1,NB
        SHT(I)=ICH+(I-1)*NCH
      ENDDO
*     BSH CONTAINS THE NEW MIXTURE AFTER SHIFTING
      DO I=1,NB
        IF(NW(I).EQ.0) THEN
          BSH(I)=0
        ELSE
          BSH(I)=SHT(NW(I))
        ENDIF
      ENDDO

      CALL LCMGET(IPMIC,'ISOTOPESDENS',NDENS)

      DO 20 IB=1,NB
        CALL XDISET(IND,MAXISO,0)
        CALL XDISET(IND2,MAXISO,0)
        I1=0
        I2=0
        DO 25 ISO=1,NISO
          IF(MIX(ISO).EQ.SHT(IB)) THEN
            I1=I1+1
            IF(I1.GE.MAXISO) CALL XABORT('@TINMIC: NUMBER OF ISOTOPES'
     +        //' OVERFLOW(1). ')
            IND(I1)=ISO
          ENDIF
 25     CONTINUE
        IF(BSH(IB).EQ.0) THEN
*       THE PROPERTIES ARE RECOVERED FROM THE THIRD LIBRARY
          DO 30 ISO2=1,NISO2
            IF(MIX2(ISO2).EQ.IWORK(IB,2)) THEN
              I2=I2+1
              IF(I2.GE.MAXISO) CALL XABORT('@TINMIC: NUMBER OF ISOTOPES'
     +          //' OVERFLOW(2). ')
              IND2(I2)=ISO2
            ENDIF
 30       CONTINUE
          IF(I1.NE.I2) CALL XABORT('@TINMIC: WRONG NUMBER OF ISOTOPES '
     +      //'IN THE NEW MIXTURE(1). ')
          DO 35 J=1,I1
            NDENS(IND(J))=DENS2(IND2(J))
            WRITE(TEXT,'(3A4)') (USED(I0,IND(J)),I0=1,3)
            WRITE(TEXT2,'(3A4)') (USED2(I0,IND2(J)),I0=1,3)
            JPMIC3=LCMGID(IPMIC3,TEXT2)
            CALL LCMSIX(IPMIC,TEXT,1)
            CALL LCMEQU(JPMIC3,IPMIC)
            CALL LCMSIX(IPMIC,' ',2)
 35       CONTINUE
*       THE PROPERTIES ARE RECOVERED FROM THE FUEL MAP LIBRARY
        ELSE
          DO 40 ISO=1,NISO
            IF(MIX(ISO).EQ.BSH(IB)) THEN
              I2=I2+1
              IF(I2.GE.MAXISO) CALL XABORT('@TINMIC: NUMBER OF ISOTOPES'
     +          //' OVERFLOW(3). ')
              IND2(I2)=ISO
            ENDIF
 40       CONTINUE
          IF(I1.NE.I2) CALL XABORT('@TINMIC: WRONG NUMBER OF ISOTOPES '
     +      //'IN THE NEW MIXTURE(2). ')
          DO 45 J=1,I1
            NDENS(IND(J))=DENS(IND2(J))
            WRITE(TEXT,'(3A4)') (USED(I0,IND(J)),I0=1,3)
            WRITE(TEXT2,'(3A4)') (USED(I0,IND2(J)),I0=1,3)
            JPMIC2=LCMGID(IPMIC2,TEXT2)
            CALL LCMSIX(IPMIC,TEXT,1)
            CALL LCMEQU(JPMIC2,IPMIC)
            CALL LCMSIX(IPMIC,' ',2)
 45       CONTINUE
        ENDIF
 20   CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VOL2,VOL,TEMP2,TEMP,DENS2,DENS)
      DEALLOCATE(USED2,USED,NAME2,NAME)
      DEALLOCATE(TYP2,TYP,TODO2,TODO,MIX2,MIX)
      RETURN
      END
