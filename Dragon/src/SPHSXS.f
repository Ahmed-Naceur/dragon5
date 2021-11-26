*DECK SPHSXS
      SUBROUTINE SPHSXS(NREA,IDIM2,NADRX,NGROUP,NL,NDATAX,NDATAP,INDX,
     1 IAD,ADRX,RDATAX,IDATAP,NOMREA,SIGS,SS2D,XS,LXS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the cross sections of an elementary calculation and single
* mixture in a Saphyb.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NREA    number of reactions in the Saphyb.
* IDIM2   second dimension of ADRX array.
* NADRX   number of address sets.
* NGROUP  number of energy groups.
* NL      maximum Legendre order (NL=1 is for isotropic scattering).
* NDATAX  number of components in RDATAX.
* NDATAP  number of components in IDATAP.
* INDX    position of isotopic set in current mixture.
* IAD     last index in ADRX.
* ADRX    index for RDATAX in the Saphyb.
* RDATAX  main cross section container in the Saphyb.
* IDATAP  index for scattering matrix information in the Saphyb.
* NOMREA  names of reactions in the Saphyb.
*
*Parameters: output
* SIGS    scattering cross sections.
* SS2D    complete scattering matrix.
* XS      cross sections per reaction.
* LXS     existence flag of each reaction.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NREA,IDIM2,NADRX,NGROUP,NL,NDATAX,NDATAP,INDX,IAD,
     1 ADRX(NREA+2,IDIM2,NADRX),IDATAP(NDATAP)
      REAL SS2D(NGROUP,NGROUP,NL),SIGS(NGROUP,NL),XS(NGROUP,NREA),
     1 RDATAX(NDATAX)
      LOGICAL LXS(NREA)
      CHARACTER NOMREA(NREA)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER FAGG,LAGG,FDGG,WGAL,FAG,LAG,FDG(NGROUP),ADR(NGROUP+1)
*----
*  FILL OUTPUT ARRAYS
*----
      IREAPR=0
      IRENTO=0
      DO IREA=1,NREA
        IF(NOMREA(IREA).EQ.'PROFIL') IREAPR=IREA
        IF(NOMREA(IREA).EQ.'TOTALE') IRENTO=IREA
      ENDDO
      SIGS(:NGROUP,:NL)=0.0
      SS2D(:NGROUP,:NGROUP,:NL)=0.0
      XS(:NGROUP,:NREA)=0.0
      NL1=ADRX(NREA+1,INDX,IAD)
      NL2=ADRX(NREA+2,INDX,IAD)
      IF((NL1.GT.NL).OR.(NL2.GT.NL)) THEN
        CALL XABORT('SPHSXS: NL OVERFLOW.')
      ENDIF
      DO IREA=1,NREA
        IOF=ADRX(IREA,INDX,IAD)
        IF(IOF.EQ.0) CYCLE
        IF(NOMREA(IREA).EQ.'DIFFUSION') THEN
          DO IL=1,NL1
            DO IGR=1,NGROUP
              SIGS(IGR,IL)=RDATAX(IOF+(IL-1)*NGROUP+IGR-1)
              LXS(IREA)=LXS(IREA).OR.(SIGS(IGR,IL).NE.0.0)
            ENDDO
          ENDDO
          IF(ADRX(IRENTO,INDX,IAD).EQ.0) THEN
            DO IGR=1,NGROUP
              XS(IGR,IRENTO)=XS(IGR,IRENTO)+RDATAX(IOF+IGR-1)
              LXS(IRENTO)=LXS(IRENTO).OR.(XS(IGR,IRENTO).NE.0.0)
            ENDDO
          ENDIF
        ELSE IF(NOMREA(IREA).EQ.'ABSORPTION') THEN
          DO IGR=1,NGROUP
            XS(IGR,IREA)=RDATAX(IOF+IGR-1)
            LXS(IREA)=LXS(IREA).OR.(XS(IGR,IREA).NE.0.0)
          ENDDO
          IF(ADRX(IRENTO,INDX,IAD).EQ.0) THEN
            DO IGR=1,NGROUP
              XS(IGR,IRENTO)=XS(IGR,IRENTO)+RDATAX(IOF+IGR-1)
              LXS(IRENTO)=LXS(IRENTO).OR.(XS(IGR,IRENTO).NE.0.0)
            ENDDO
          ENDIF
        ELSE IF(NOMREA(IREA).EQ.'PROFIL') THEN
          CYCLE
        ELSE IF(NOMREA(IREA).EQ.'TRANSFERT') THEN
          IF(IREAPR.EQ.0) CALL XABORT('SPHSXS: MISSING PROFILE INFO.')
          IPROF=ADRX(IREAPR,INDX,IAD)
          FAGG=IDATAP(IPROF)
          LAGG=IDATAP(IPROF+1)
          FDGG=IDATAP(IPROF+2)
          WGAL=IDATAP(IPROF+3)
          FAG=IDATAP(IPROF+4)
          LAG=IDATAP(IPROF+5)
          DO IGR=1,NGROUP
            FDG(IGR)=IDATAP(IPROF+5+IGR)
            ADR(IGR)=IDATAP(IPROF+5+NGROUP+IGR)
          ENDDO
          ADR(NGROUP+1)=IDATAP(IPROF+6+2*NGROUP)
          JOFS=0
          DO IL=1,NL2
            ZIL=REAL(2*IL-1)
            IF(WGAL.NE.0) THEN
              DO IGR=FAGG,LAGG
                DO JGR=FDGG,FDGG+WGAL-1
                  SS2D(IGR,JGR,IL)=RDATAX(IOF+JOFS)/ZIL
                  JOFS=JOFS+1
                  LXS(IREA)=LXS(IREA).OR.(SS2D(IGR,JGR,IL).NE.0.0)
                ENDDO
              ENDDO
            ENDIF
            DO IGR=FAG,LAG
              DO JGR=FDG(IGR),FDG(IGR)+(ADR(IGR+1)-ADR(IGR))-1
                SS2D(IGR,JGR,IL)=RDATAX(IOF+JOFS)/ZIL
                JOFS=JOFS+1
                LXS(IREA)=LXS(IREA).OR.(SS2D(IGR,JGR,IL).NE.0.0)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO IGR=1,NGROUP
            XS(IGR,IREA)=RDATAX(IOF+IGR-1)
            LXS(IREA)=LXS(IREA).OR.(XS(IGR,IREA).NE.0.0)
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
