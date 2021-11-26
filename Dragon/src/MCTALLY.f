*DECK MCTALLY
      SUBROUTINE MCTALLY(ITALLY,NFREG,NMIX,NGRP,NL,NFM,NDEL,NED,NBSCO,
     1 NMERGE,NGCOND,IREG,IGR,NU,MATCOD,IMERGE,INDGRP,XSM,XST,XSS,
     2 XSN2N,XSN3N,XSSNN,XSNUSI,XSCHI,XSEDI,SCORE1,SCORE2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Score for effective multiplication factor and macrolib information.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* ITALLY  type of tally (=1 score effective multiplication factor;
*         =2 also score macrolib information).
* NFREG   number of regions.
* NMIX    number of material mixtures.
* NGRP    number of energy groups.
* NL      number of Legendre orders required in the estimations
*         (NL=1 or higher).
* NFM     number of fissile isotopes.
* NDEL    number of delayed precursor groups.
* NED     number of extra edit vectors.
* NBSCO   number of macrolib-related scores.
* NMERGE  number of homogenized regions.
* NGCOND  number of condensed energy groups.
* IREG    index of the region where the particle is located.
* IGR     index of the energy group of the particle.
* NU      particle weight.
* MATCOD  region material.
* IMERGE  homogenized regions indices.
* INDGRP  condensed groups indices.
* XSM     maximum macroscopic total cross section.
* XST     total macroscopic cross sections for each mixture and energy
*         group.
* XSS     total scattering cross sections for each mixture and energy
*         group.
* XSN2N   N2N macroscopic cross sections for each mixture and energy
*         group.
* XSN3N   N3N macroscopic cross sections for each mixture and energy
*         group.
* XSSNN   in-group and out-of-group macroscopic transfert cross sections
*         for each mixture.
* XSNUSI  the values of Nu time the fission cross sections for each
*         isotope per mixture and energy group.
* XSCHI   fission spectrum for isotopes per mixture and energy group.
* XSEDI   extra edit cross sections for each mixture and energy group.
*
*Parameters: input/output
* SCORE1  score for total flux and effective multiplication factor.
* SCORE2  macrolib score matrix.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITALLY,NFREG,NMIX,NGRP,NL,NFM,NDEL,NED,NBSCO,NMERGE,
     1 NGCOND,IREG,IGR,MATCOD(NFREG),IMERGE(NFREG),INDGRP(NGRP)
      REAL NU,XSM,XST(NMIX,NGRP),XSS(NMIX,NGRP,NL),XSN2N(NMIX,NGRP),
     1 XSN3N(NMIX,NGRP),XSSNN(NMIX,NGRP,NGRP,NL),
     2 XSNUSI(NMIX,NFM,NGRP,1+NDEL),XSCHI(NMIX,NFM,NGRP,1+NDEL),
     3 XSEDI(NMIX,NGRP,NED),SCORE1(3),SCORE2(NBSCO,NMERGE,NGCOND)
*----
*  LOCAL VARIABLES
*----
      INTEGER IBM,IFM,IED,JGR,IBMCD,IGRCD,JGRCD,IOF,IL,IDEL
      REAL WW
*
      WW=NU/XSM
      IF(IREG.GT.0) THEN
        SCORE1(1)=SCORE1(1)+WW
        IBM=MATCOD(IREG)
        IF(IBM.GT.0) THEN
          DO IFM=1,NFM
            SCORE1(2)=SCORE1(2)+WW*XSNUSI(IBM,IFM,IGR,1)
          ENDDO
          SCORE1(3)=SCORE1(3)+WW*(XST(IBM,IGR)-XSS(IBM,IGR,1)-2.0*
     1    XSN2N(IBM,IGR)-3.0*XSN3N(IBM,IGR))
        ENDIF
        IF(ITALLY.EQ.2) THEN
          IBMCD=IMERGE(IREG)
          IF(IBMCD.EQ.0) GO TO 10
          IGRCD=INDGRP(IGR)
          IF(IGRCD.EQ.0) GO TO 10
          SCORE2(1,IBMCD,IGRCD)=SCORE2(1,IBMCD,IGRCD)+WW
          IF(IBM.EQ.0) GO TO 10
          SCORE2(2,IBMCD,IGRCD)=SCORE2(2,IBMCD,IGRCD)+WW*XST(IBM,IGR)
          SCORE2(3,IBMCD,IGRCD)=SCORE2(3,IBMCD,IGRCD)+WW*XSS(IBM,IGR,1)
          SCORE2(4,IBMCD,IGRCD)=SCORE2(4,IBMCD,IGRCD)+WW*XSN2N(IBM,IGR)
          SCORE2(5,IBMCD,IGRCD)=SCORE2(5,IBMCD,IGRCD)+WW*XSN3N(IBM,IGR)
          IOF=5
          DO IL=1,NL
            DO JGR=1,NGRP
             JGRCD=INDGRP(JGR)
             SCORE2(IOF+JGRCD,IBMCD,IGRCD)=SCORE2(IOF+JGRCD,IBMCD,IGRCD)
     1       +WW*XSSNN(IBM,JGR,IGR,IL)
            ENDDO
            IOF=IOF+NGCOND
          ENDDO
          DO IDEL=1,1+NDEL
            DO IFM=1,NFM
              SCORE2(IOF+IFM,IBMCD,IGRCD)=SCORE2(IOF+IFM,IBMCD,IGRCD)+
     1        WW*XSNUSI(IBM,IFM,IGR,IDEL)
            ENDDO
            IOF=IOF+NFM
            DO IFM=1,NFM
              DO JGR=1,NGRP
                JGRCD=INDGRP(JGR)
                SCORE2(IOF+IFM,IBMCD,JGRCD)=SCORE2(IOF+IFM,IBMCD,JGRCD)+
     1          WW*XSNUSI(IBM,IFM,IGR,IDEL)*XSCHI(IBM,IFM,JGR,IDEL)
              ENDDO
            ENDDO
            IOF=IOF+NFM
          ENDDO
          DO IED=1,NED
            SCORE2(IOF+IED,IBMCD,IGRCD)=SCORE2(IOF+IED,IBMCD,IGRCD)+WW*
     1      XSEDI(IBM,IGR,IED)
          ENDDO
        ENDIF
      ENDIF
   10 RETURN
      END
