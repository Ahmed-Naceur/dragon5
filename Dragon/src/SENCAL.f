*DECK SENCAL
      SUBROUTINE SENCAL(IPSENS,IPLIB,IPRINT,NR,NG,NI,NANIS,NAMISO,
     >                  MELISO,MAT,DENISO,KEFF,P,D,NAMISC,ISOC,NIC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute sensitivity profiles.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Laville, G. Marleau
*
*Parameters: input/output
* IPSENS  LCM Sensitivity object address.
* IPLIB   LCM Library object address.
* IPRINT  print level.
* NR      number of region in Tracking object.
* NG      number of energy group in Library object.
* NI      number of isotope/mixture.
* NANIS   anisotropy order kept.
* NAMISO  name of the isotope/mixture.
* MELISO  mixture associated with the isotope/mixture.
* MAT     mixture of each region.
* DENISO  density of each mixtures.
* KEFF    keff.
* P       matrix calcution for sensitivity analysis.
* D       weighting coefficient for sensitivity analysis.
* NAMISC  independent isotopes names.
* ISOC    independent isotope number associated with isotope/mixture.
* NIC     number of independent isotopes names.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT      NONE
*----
*  Parameters
*---- 
      INTEGER       IOUT,NS
      CHARACTER     NAMSBR*6
      PARAMETER    (IOUT=6,NS=9,NAMSBR='SENCAL')      
*----
*  Suboutine arguements
*----
      TYPE(C_PTR)   IPLIB
      INTEGER       IPSENS,IPRINT,NR,NG,NI,NANIS,NAMISO(3,NI),
     >              MELISO(NI),MAT(NR)
      REAL          DENISO(NI),KEFF,P(NR,NANIS,NG,NG),D
      INTEGER       NAMISC(2,NI),ISOC(NI),NIC
*----
*  Local variables
*----
      TYPE(C_PTR)   KPISO
      REAL          ZERO
      CHARACTER     ISONAM*12,ISONAC*8,CL*2,HSMG*131
      INTEGER       ILENG,ITYLCM,IR,IS,IG,IP,JG,ISF,ISN,ISC,IL,II,
     >              ISOMEL,IIC,NGG
      DOUBLE PRECISION DD,DDD,SENRIG
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IDS
      REAL, ALLOCATABLE, DIMENSION(:) :: NUSIGF,CHI,CAPT,SIGS,SCAT,
     > NUBAR,NFTOT,SIGD,SIGA,SIGG,SIGP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  SENT,SENC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::  SENG,SENTI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) ::  SENRG,SENGI
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  Scratch storage allocation
*   NUSIGF  Vector for nubar multiplied by the fission cross-section
*   CHI     Vector for fission spectra.
*   SENRG   Sensitivities for each region and group.
*   SENG    Sensitivities for each group integrated on each region.
*   SENT    Integrated sensitivities.
*   CAPT    Vector for capture cross-sections.
*   SIGS    Vector for scattering cross-sections.
*   SCAT    Vector for transfter sections.
*   NUBAR   Vector for nubar.
*   NFTOT   Vector for fission cross-sections.
*   SIGD    Vector for (n,d) cross-sections.
*   SIGA    Vector for (n,a) cross-sections.
*   SIGG    Vector for (n,g) cross-sections.
*   SIGP    Vector for (n,p) cross-sections.
*   SENC    Vector used for sensitivities to the fission spectra CHI.
*   IJJ     Highest energy group for which the scattering does no
*           vanish.
*   NJJ     Number of energy group for which the scattering does not
*           vanish.
*   IDS     Used to compute integrated isotope sensitivities.
*   SENGI   Sensitivities for integrated isotope for each group.
*   SENTI   Sensitivities for integrated isotope.
*----
      ALLOCATE(IJJ(NG),NJJ(NG),IDS(NS,NIC))
      ALLOCATE(CAPT(NG),NUSIGF(NG),CHI(NG),SIGS(NG),SCAT(NG*NG),
     < NUBAR(NG),NFTOT(NG),SIGD(NG),SIGA(NG),SIGG(NG),SIGP(NG))
      ALLOCATE(SENRG(NR,NG,NS),SENG(NG,NS),SENT(NS),SENC(NR),
     < SENGI(NG,NS,NIC),SENTI(NS,NIC))
      ALLOCATE(IPISO(NI))
*----
*  Initialize the directory of the isotope/mixture
*  in the library and initialize the SEN information
*  support and the cross section support for the isotope
*----
      ZERO=0.0
      NGG=NG*NG
*----
*  Loop over isotopes
*----
      CALL XDISET(IDS,NS*NIC,0)
      CALL XDDSET(SENGI,NG*NS*NIC,0.0D0)
      CALL XDDSET(SENTI,NS*NIC,0.0D0)
      CALL LIBIPS(IPLIB,NI,IPISO)
      DO II=1,NI
      WRITE(ISONAM,'(3A4)') NAMISO(1,II),NAMISO(2,II),NAMISO(3,II)
      KPISO=IPISO(II) ! set II-th isotope
      IF(.NOT.C_ASSOCIATED(KPISO)) THEN
         WRITE(HSMG,'(17HSENCAL: ISOTOPE '',A12,7H'' (ISO=,I8,5H) IS ,
     1   30HNOT AVAILABLE IN THE MICROLIB.)') ISONAM,II
         CALL XABORT(HSMG)
      ENDIF
      IIC=ABS(ISOC(II))
      ISOMEL=MELISO(II)
      DD=DBLE(DENISO(II)/D)
      DDD=DBLE(DD/KEFF)
*----
*  Process isotope
*----
      CALL XDDSET(SENRG,NR*NG*NS,0.0D0)
      CALL XDDSET(SENG,NG*NS,0.0D0)
      CALL XDDSET(SENT,NS,0.0D0)
      CALL XDDSET(SENC,NR,0.0D0)
      CALL XDRSET(SIGD,NG,0.0)
      CALL XDRSET(NUSIGF,NG,0.0)
      CALL XDRSET(NUBAR,NG,0.0)
      CALL XDRSET(NFTOT,NG,0.0)
      CALL XDRSET(CHI,NG,0.0)
      CALL XDRSET(SIGA,NG,0.0)     
      CALL XDRSET(SIGP,NG,0.0)
      CALL XDRSET(SIGG,NG,0.0)
      CALL XDRSET(CAPT,NG,0.0)
      CALL XDRSET(SIGS,NG,0.0)
      CALL XDRSET(SCAT,NGG,0.0)
      CALL XDISET(IJJ,NG,0)
      CALL XDISET(NJJ,NG,0)   
*----
*  (n,g) sensitivity calculation
*----
      IS=1
      CALL LCMLEN(KPISO,'NG',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(KPISO,'NG',SIGG)
        DO IR=1,NR
          IF(MAT(IR).EQ.ISOMEL) THEN
            DO IG=1,NG
              DO IL=0,NANIS-1
                 SENRG(IR,IG,IS)=SENRG(IR,IG,IS)-DD*
     >           DBLE(((2*IL)+1)*SIGG(IG)*P(IR,IL+1,IG,IG))
              ENDDO
              SENG(IG,IS)=SENG(IG,IS)+SENRG(IR,IG,IS)
              SENT(IS)=SENT(IS)+SENRG(IR,IG,IS)
            ENDDO
          ENDIF
        ENDDO
*----
*  Add contribution to integrated isotope
*----
        DO IG=1,NG
          SENGI(IG,IS,IIC)=SENGI(IG,IS,IIC)+SENG(IG,IS)
        ENDDO
        SENTI(IS,IIC)=SENTI(IS,IIC)+SENT(IS)
        IDS(IS,IIC)=IDS(IS,IIC)+1
*----
*  Send to IPSENS
*----
        WRITE(IPSENS,7000) ISONAM(1:8),'n,g    ',102,-ISOMEL,ZERO
        WRITE(IPSENS,7001) SENT(IS),ZERO,ZERO
        WRITE(IPSENS,7002) (SENG(IG,IS),IG=1,NG)
*----
*  Print information if required for (n,g)
*----
        IF(IPRINT .GE. 5) THEN
          WRITE(IOUT,6000) '  (n,g)',ISONAM,SENT(IS)
          WRITE(IOUT,6001) (SENG(IG,IS),IG=1,NG)
        ENDIF
      ENDIF      
*----
*  (n,p) sensitivity calculation
*----
      IS=2
      CALL LCMLEN(KPISO,'NP',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(KPISO,'NP',SIGP)
        DO IR=1,NR
          IF(MAT(IR).EQ.ISOMEL) THEN
            DO IG=1,NG
              DO IL=0,NANIS-1
                 SENRG(IR,IG,IS)=SENRG(IR,IG,IS)-DD*
     >           DBLE(((2*IL)+1)*SIGP(IG)*P(IR,IL+1,IG,IG))
              ENDDO
              SENG(IG,IS)=SENG(IG,IS)+SENRG(IR,IG,IS)
              SENT(IS)=SENT(IS)+SENRG(IR,IG,IS)
            ENDDO
          ENDIF
        ENDDO      
*----
*  Add contribution to integrated isotope
*----
        DO IG=1,NG
          SENGI(IG,IS,IIC)=SENGI(IG,IS,IIC)+SENG(IG,IS)
        ENDDO
        SENTI(IS,IIC)=SENTI(IS,IIC)+SENT(IS)
        IDS(IS,IIC)=IDS(IS,IIC)+1
*----
*  Send to IPSENS
*----
        WRITE(IPSENS,7000) ISONAM(1:8),'n,p    ',103,-ISOMEL,ZERO
        WRITE(IPSENS,7001) SENT(IS),ZERO,ZERO
        WRITE(IPSENS,7002) (SENG(IG,IS),IG=1,NG)
*----
*  Print information if required for (n,p)
*----
        IF(IPRINT .GE. 5) THEN
          WRITE(IOUT,6000) '  (n,p)',ISONAM,SENT(IS)
          WRITE(IOUT,6001) (SENG(IG,IS),IG=1,NG)
        ENDIF
      ENDIF      
*----
*  (n,d) sensitivity calculation
*----
      IS=3
      CALL LCMLEN(KPISO,'ND',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(KPISO,'ND',SIGD)
        DO IR=1,NR
          IF(MAT(IR).EQ.ISOMEL) THEN
            DO IG=1,NG
              DO IL=0,NANIS-1
                 SENRG(IR,IG,IS)=SENRG(IR,IG,IS)-DD*
     >           DBLE(((2*IL)+1)*SIGD(IG)*P(IR,IL+1,IG,IG))
              ENDDO            
              SENG(IG,IS)=SENG(IG,IS)+SENRG(IR,IG,IS)
              SENT(IS)=SENT(IS)+SENRG(IR,IG,IS)
            ENDDO
          ENDIF
        ENDDO
*----
*  Add contribution to integrated isotope
*----
        DO IG=1,NG
          SENGI(IG,IS,IIC)=SENGI(IG,IS,IIC)+SENG(IG,IS)
        ENDDO
        SENTI(IS,IIC)=SENTI(IS,IIC)+SENT(IS)
        IDS(IS,IIC)=IDS(IS,IIC)+1
*----
*  Send to IPSENS
*----
        WRITE(IPSENS,7000) ISONAM(1:8),'n,d    ',104,-ISOMEL,ZERO
        WRITE(IPSENS,7001) SENT(IS),ZERO,ZERO
        WRITE(IPSENS,7002) (SENG(IG,IS),IG=1,NG)
*----
*  Print information if required for (n,d)
*----
        IF(IPRINT .GE. 5) THEN
          WRITE(IOUT,6000) '  (n,d)',ISONAM,SENT(IS)
          WRITE(IOUT,6001) (SENG(IG,IS),IG=1,NG)
        ENDIF
      ENDIF
*----
*  (n,a) sensitivity calculation
*----
      IS=4
      CALL LCMLEN(KPISO,'NA',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(KPISO,'NA',SIGA)
        DO IR=1,NR
          IF(MAT(IR).EQ.ISOMEL) THEN
            DO IG=1,NG
              DO IL=0,NANIS-1
                 SENRG(IR,IG,IS)=SENRG(IR,IG,IS)-DD*
     >           DBLE(((2*IL)+1)*SIGA(IG)*P(IR,IL+1,IG,IG))
              ENDDO            
              SENG(IG,IS)=SENG(IG,IS)+SENRG(IR,IG,IS)
              SENT(IS)=SENT(IS)+SENRG(IR,IG,IS)
            ENDDO
          ENDIF
        ENDDO
*----
*  Add contribution to integrated isotope
*----
        DO IG=1,NG
          SENGI(IG,IS,IIC)=SENGI(IG,IS,IIC)+SENG(IG,IS)
        ENDDO
        SENTI(IS,IIC)=SENTI(IS,IIC)+SENT(IS)
        IDS(IS,IIC)=IDS(IS,IIC)+1
*----
*  Send to IPSENS
*----
        WRITE(IPSENS,7000) ISONAM(1:8),'n,a    ',107,-ISOMEL,ZERO
        WRITE(IPSENS,7001) SENT(IS),ZERO,ZERO
        WRITE(IPSENS,7002) (SENG(IG,IS),IG=1,NG)
*----
*  Print information if required for (n,a)
*----
        IF(IPRINT .GE. 5) THEN
          WRITE(IOUT,6000) '  (n,a)',ISONAM,SENT(IS)
          WRITE(IOUT,6001) (SENG(IG,IS),IG=1,NG)
        ENDIF
      ENDIF
*----
*  Capture sensitivity calculation
*----
      IS=5
      DO IG=1,NG
        CAPT(IG)=SIGD(IG)+SIGA(IG)+SIGP(IG)+SIGG(IG)
      ENDDO
      DO IR=1,NR
        IF(MAT(IR).EQ.ISOMEL) THEN
          DO IG=1,NG
              DO IL=0,NANIS-1
                 SENRG(IR,IG,IS)=SENRG(IR,IG,IS)-DD*
     >           DBLE(((2*IL)+1)*CAPT(IG)*P(IR,IL+1,IG,IG))
              ENDDO          
            SENG(IG,IS)=SENG(IG,IS)+SENRG(IR,IG,IS)
            SENT(IS)=SENT(IS)+SENRG(IR,IG,IS)
          ENDDO
        ENDIF
      ENDDO
*----
*  Add contribution to integrated isotope
*----
      DO IG=1,NG
        SENGI(IG,IS,IIC)=SENGI(IG,IS,IIC)+SENG(IG,IS)
      ENDDO
      SENTI(IS,IIC)=SENTI(IS,IIC)+SENT(IS)
      IDS(IS,IIC)=IDS(IS,IIC)+1
*----
*  Send to IPSENS
*----
      WRITE(IPSENS,7000) ISONAM(1:8),'capture',101,-ISOMEL,ZERO
      WRITE(IPSENS,7001) SENT(IS),ZERO,ZERO
      WRITE(IPSENS,7002) (SENG(IG,IS),IG=1,NG)
*----
*  Print information if required for capture
*----
      IF(IPRINT .GE. 5) THEN
        WRITE(IOUT,6000) 'capture',ISONAM,SENT(IS)
        WRITE(IOUT,6001) (SENG(IG,IS),IG=1,NG)
      ENDIF
*----
*  Scattering sensitivity calculation
*----
      IS=6
      DO IL=1,NANIS
        WRITE(CL,'(I2.2)') IL-1
        CALL LCMLEN(KPISO,'SIGS'//CL,ILENG,ITYLCM)
        IF(ILENG.GT.0) THEN
          IF( IL .EQ. 1) THEN
            CALL LCMGET(KPISO,'SIGS'//CL,SIGS)
          ENDIF
          CALL LCMGET(KPISO,'SCAT'//CL,SCAT)
          CALL LCMGET(KPISO,'NJJS'//CL,NJJ)
          CALL LCMGET(KPISO,'IJJS'//CL,IJJ)
*----
*  Decompress scattering matrix
*  SCAT(JG,IG) is from IG to JG
*----
          DO IR=1,NR
            IF(MAT(IR).EQ.ISOMEL) THEN
              IP=1
              DO JG=1,NG
                DO IG=IJJ(JG),IJJ(JG)-NJJ(JG)+1,-1
                  IF(IG .EQ. JG) THEN
                    SENRIG=DBLE(SCAT(IP)*P(IR,IL,IG,JG))-
     >           DBLE(((2*IL)-1)*SIGS(IG)*P(IR,IL,IG,IG))
                    SENRG(IR,IG,IS)=SENRG(IR,IG,IS)+SENRIG
                  ELSE
                    SENRG(IR,IG,IS)=SENRG(IR,IG,IS)+
     >              DBLE(SCAT(IP)*P(IR,IL,IG,JG))
                  ENDIF
                  IP=IP+1
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DO IR=1,NR
        IF(MAT(IR).EQ.ISOMEL) THEN
          DO IG=1,NG
            SENRG(IR,IG,IS)=DD*SENRG(IR,IG,IS)
            SENG(IG,IS)=SENG(IG,IS)+SENRG(IR,IG,IS)
            SENT(IS)=SENT(IS)+SENRG(IR,IG,IS)
          ENDDO
        ENDIF
      ENDDO
*----
*  Add contribution to integrated isotope
*----
      DO IG=1,NG
        SENGI(IG,IS,IIC)=SENGI(IG,IS,IIC)+SENG(IG,IS)
      ENDDO
      SENTI(IS,IIC)=SENTI(IS,IIC)+SENT(IS)
      IDS(IS,IIC)=IDS(IS,IIC)+1
*----
*  Send to IPSENS
*----
      WRITE(IPSENS,7000) ISONAM(1:8),'scatter',0,-ISOMEL,ZERO
      WRITE(IPSENS,7001) SENT(IS),ZERO,ZERO
      WRITE(IPSENS,7002) (SENG(IG,IS),IG=1,NG)
*----
*  Print information if required for scattering
*----
      IF(IPRINT .GE. 5) THEN
        WRITE(IOUT,6000) 'scatter',ISONAM,SENT(IS)
        WRITE(IOUT,6001) (SENG(IG,IS),IG=1,NG)
      ENDIF
*----
*  Check if the isotope is fissile and get the
*  cross section informations
*----     
      CALL LCMLEN(KPISO,'NUSIGF',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(KPISO,'NUSIGF',NUSIGF)
        CALL LCMGET(KPISO,'CHI',CHI)
        CALL LCMGET(KPISO,'NFTOT',NFTOT)
        DO IG=1,NG
          IF(NFTOT(IG).NE.0.0) THEN
            NUBAR(IG)=NUSIGF(IG)/NFTOT(IG)
          ENDIF
        ENDDO
*----
*  Fission sensitivity calculation if fissile isotope
*  Fission (ISF), nubar (ISN) and chi (ISC)
*---- 
        ISF=7
        ISN=8
        ISC=9
        DO IR=1,NR
          IF(MAT(IR).EQ.ISOMEL) THEN
            DO IG=1,NG
              DO IL=0,NANIS-1
*  Fission
                SENRG(IR,IG,ISF)=SENRG(IR,IG,ISF)-DD*
     >           DBLE(((2*IL)+1)*NFTOT(IG)*P(IR,IL+1,IG,IG))
              ENDDO            
              DO JG=1,NG
* Fission         
                SENRG(IR,IG,ISF)=SENRG(IR,IG,ISF)+DDD*
     >           DBLE(NUSIGF(IG)*CHI(JG)*P(IR,1,IG,JG))
* Nubar
                SENRG(IR,IG,ISN)=SENRG(IR,IG,ISN)+
     >           DBLE(NUSIGF(IG)*CHI(JG)*P(IR,1,IG,JG))
* Chi
                SENRG(IR,IG,ISC)=SENRG(IR,IG,ISC)+
     >           DBLE(NUSIGF(JG)*CHI(IG)*P(IR,1,JG,IG))
              ENDDO
* Fission         
              SENG(IG,ISF)=SENG(IG,ISF)+SENRG(IR,IG,ISF)
              SENT(ISF)=SENT(ISF)+SENRG(IR,IG,ISF)
* Nubar
              SENRG(IR,IG,ISN)=DDD*SENRG(IR,IG,ISN)
              SENG(IG,ISN)=SENG(IG,ISN)+SENRG(IR,IG,ISN)
              SENT(ISN)=SENT(ISN)+SENRG(IR,IG,ISN)
* Chi
              SENRG(IR,IG,ISC)=DDD*SENRG(IR,IG,ISC)
              SENC(IR)=SENC(IR)+SENRG(IR,IG,ISC)     
            ENDDO
          ENDIF
        ENDDO
*----
* Modification of sensitivty to Chi
*----        
        DO IR=1,NR
          IF(MAT(IR).EQ.ISOMEL) THEN
            DO IG=1,NG
              SENRG(IR,IG,ISC)=SENRG(IR,IG,ISC)-
     >           DBLE(CHI(IG))*SENC(IR)
              SENG(IG,ISC)=SENG(IG,ISC)+SENRG(IR,IG,ISC)
              SENT(ISC)=SENT(ISC)+SENRG(IR,IG,ISC)   
            ENDDO
          ENDIF
        ENDDO
*----
*  Add contribution to integrated isotope
*----
        DO IG=1,NG
          SENGI(IG,ISF,IIC)=SENGI(IG,ISF,IIC)+SENG(IG,ISF)
          SENGI(IG,ISN,IIC)=SENGI(IG,ISN,IIC)+SENG(IG,ISN)
          SENGI(IG,ISC,IIC)=SENGI(IG,ISC,IIC)+SENG(IG,ISC)
        ENDDO
        SENTI(ISF,IIC)=SENTI(ISF,IIC)+SENT(ISF)
        SENTI(ISN,IIC)=SENTI(ISN,IIC)+SENT(ISN)
        SENTI(ISC,IIC)=SENTI(ISC,IIC)+SENT(ISC)
        IDS(ISF,IIC)=IDS(ISF,IIC)+1
        IDS(ISN,IIC)=IDS(ISN,IIC)+1
        IDS(ISC,IIC)=IDS(ISC,IIC)+1
*----
*  Send to IPSENS
*----
        WRITE(IPSENS,7000) ISONAM(1:8),'fission',18,-ISOMEL,ZERO
        WRITE(IPSENS,7001) SENT(ISF),ZERO,ZERO
        WRITE(IPSENS,7002) (SENG(IG,ISF),IG=1,NG)
        WRITE(IPSENS,7000) ISONAM(1:8),'nubar  ',452,-ISOMEL,ZERO
        WRITE(IPSENS,7001) SENT(ISN),ZERO,ZERO
        WRITE(IPSENS,7002) (SENG(IG,ISN),IG=1,NG)
        WRITE(IPSENS,7000) ISONAM(1:8),'chi    ',1018,-ISOMEL,ZERO
        WRITE(IPSENS,7001) SENT(ISC),ZERO,ZERO
        WRITE(IPSENS,7002) (SENG(IG,ISC),IG=1,NG)
*----
*  Print information if required 
*----
        IF(IPRINT .GE. 5) THEN
          WRITE(IOUT,6000) 'fission',ISONAM,SENT(ISF)
          WRITE(IOUT,6001) (SENG(IG,ISF),IG=1,NG)
          WRITE(IOUT,6000) '  nubar',ISONAM,SENT(ISN)
          WRITE(IOUT,6001) (SENG(IG,ISN),IG=1,NG)
          WRITE(IOUT,6000) '    chi',ISONAM,SENT(ISC)
          WRITE(IOUT,6001) (SENG(IG,ISC),IG=1,NG)
        ENDIF
*----
*  End of the loop for fissile isotope
*----        
      ENDIF 
*----
*  End loop for isotopes
*----
      ENDDO
*----
*  Save integrated contributions
*----
      DO IIC=1,NIC
        WRITE(ISONAC,'(2A4)') NAMISC(1,IIC),NAMISC(2,IIC)
*----
*  (n,g) sensitivity calculation
*----
        DO IS=1,NS
          IF(IDS(IS,IIC).GT.0) THEN
            IF(IS .EQ.1) THEN
              WRITE(IPSENS,7000) ISONAC,'n,g    ', 102,0,ZERO
            ELSE IF(IS .EQ.2) THEN
              WRITE(IPSENS,7000) ISONAC,'n,p    ', 103,0,ZERO
            ELSE IF(IS .EQ.3) THEN
              WRITE(IPSENS,7000) ISONAC,'n,d    ', 104,0,ZERO
            ELSE IF(IS .EQ.4) THEN
              WRITE(IPSENS,7000) ISONAC,'n,a    ', 107,0,ZERO
            ELSE IF(IS .EQ.5) THEN
              WRITE(IPSENS,7000) ISONAC,'capture', 101,0,ZERO
            ELSE IF(IS .EQ.6) THEN
              WRITE(IPSENS,7000) ISONAC,'scatter',   0,0,ZERO
            ELSE IF(IS .EQ.7) THEN
              WRITE(IPSENS,7000) ISONAC,'fission',  18,0,ZERO
            ELSE IF(IS .EQ.8) THEN
              WRITE(IPSENS,7000) ISONAC,'nubar  ', 452,0,ZERO
            ELSE IF(IS .EQ.9) THEN
              WRITE(IPSENS,7000) ISONAC,'chi    ',1018,0,ZERO
            ENDIF
            WRITE(IPSENS,7001) SENTI(IS,IIC),ZERO,ZERO
            WRITE(IPSENS,7002) (SENGI(IG,IS,IIC),IG=1,NG)
*----
*  Print information if required for (n,g)
*----
            IF(IPRINT .GE. 5) THEN
              IF(IS .EQ.1) THEN
                WRITE(IOUT,6005) '  (n,g)',ISONAC,SENTI(IS,IIC)
              ELSE IF(IS .EQ.2) THEN
                WRITE(IOUT,6005) '  (n,p)',ISONAC,SENTI(IS,IIC)
              ELSE IF(IS .EQ.3) THEN
                WRITE(IOUT,6005) '  (n,d)',ISONAC,SENTI(IS,IIC)
              ELSE IF(IS .EQ.4) THEN
                WRITE(IOUT,6005) '  (n,a)',ISONAC,SENTI(IS,IIC)
              ELSE IF(IS .EQ.5) THEN
                WRITE(IOUT,6005) 'capture',ISONAC,SENTI(IS,IIC)
              ELSE IF(IS .EQ.6) THEN
                WRITE(IOUT,6005) 'scatter',ISONAC,SENTI(IS,IIC)
              ELSE IF(IS .EQ.7) THEN
                WRITE(IOUT,6005) 'fission',ISONAC,SENTI(ISF,IIC)
              ELSE IF(IS .EQ.8) THEN
                WRITE(IOUT,6005) '  nubar',ISONAC,SENTI(ISN,IIC)
              ELSE IF(IS .EQ.9) THEN
                WRITE(IOUT,6005) '    chi',ISONAC,SENTI(ISC,IIC)
              ENDIF
              WRITE(IOUT,6001) (SENGI(IG,IS,IIC),IG=1,NG)
            ENDIF
          ENDIF
        ENDDO
*----
*  End loop for isotopes
*----
      ENDDO
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(IPISO)
      DEALLOCATE(SENTI,SENGI,SENC,SENT,SENG,SENRG)
      DEALLOCATE(SIGP,SIGG,SIGA,SIGD,NFTOT,NUBAR,SCAT,SIGS,CHI,NUSIGF,
     < CAPT)
      DEALLOCATE(IDS,NJJ,IJJ)
      RETURN
*----
*  Formats
*----
 6000 FORMAT('Name of the isotope/mixture: ',7X,
     > 'Integrated sensitivity to :',
     >A7/A12,44X,1P,E14.6/'Sensitivity profiles')
 6001 FORMAT(1P,5E14.6)
 6005 FORMAT('Name of the isotope/mixture: ',7X,
     > 'Integrated sensitivity to :',
     >A7/A8,48X,1P,E14.6/'Sensitivity profiles')
 7000 FORMAT(A8,4X,1X,A7,'  -0000',I7,I7,1P,E14.6)
 7001 FORMAT(1P,3E14.6)
 7002 FORMAT(1P,5E14.6)
      END
