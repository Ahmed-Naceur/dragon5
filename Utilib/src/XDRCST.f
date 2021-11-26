*DECK XDRCST
      FUNCTION XDRCST(CSTNAM,CSTUNT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To return the physical constants in the required units or get constant
* for converting between units.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
*
*Author(s): G. Marleau.
*
*Parameters: input
* CSTNAM  constant name or initial units for conversion where
*         \begin{itemize}
*         \item \moc{CSTNAM}=\moc{Avogadro} is for Avogadro number
*         with units in N/moles;
*         \item \moc{CSTNAM}=\moc{Plank} is for Plank constant
*         with units in J$\times$s, MeV$\times$s or eV$\times$s;
*         \item \moc{CSTNAM}=\moc{Boltzmann} is for Boltzmann constant
*         with units in J/K, MeV/K or eV/K ;
*         \item \moc{CSTNAM}=\moc{Neutron mass} is for neutron mass
*         with units in kg, amu, MeV or eV;
*         \item \moc{CSTNAM}=\moc{Proton mass} is for proton mass
*         with units in kg, amu, MeV or eV;
*         \item \moc{CSTNAM}=\moc{kg} is the factor to transform kg 
*         into amu, MeV, eV or J;
*         \item \moc{CSTNAM}=\moc{amu} is the factor to transform amu 
*         into kg, MeV, eV or J;
*         \item \moc{CSTNAM}=\moc{eV} is the factor to transform eV 
*         into J or K;
*         \item \moc{CSTNAM}=\moc{K} is the factor to transform K 
*         into J or eV;
*         \item \moc{CSTNAM}=\moc{J} is the factor to transform J 
*         into eV or K;
*         \item \moc{CSTNAM}=\moc{Pi} is for $\pi$
*         without units.
*         \end{itemize}
* CSTUNT  units for the constant or final units as described for
*         \moc{CSTNAM}.
*
*Parameters: input
* XDRCST  numerical value of the constant in required unit.
*
*References: 
*  Peter J. Mohr and Barry N. Taylor, CODATA Recommended Values of the 
*  Fundamental Physical Constants: 2002. (to be published)
*  http://physics.nist.gov/constants
*  Last visit: September 04, 2004
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      CHARACTER        CSTNAM*(*)
      CHARACTER        CSTUNT*(*)
      DOUBLE PRECISION XDRCST
*----
*  Parameters
*----
      INTEGER          IOUT
      PARAMETER       (IOUT=6)
      CHARACTER        NAMSBR*6
      PARAMETER       (NAMSBR='XDRCST')
      DOUBLE PRECISION EVJOULE,AMUKG,SPEEDL
      PARAMETER       (EVJOULE=1.60217653D-19,AMUKG=1.66053886D-27,
     >                 SPEEDL=2.99792458D+08)
      DOUBLE PRECISION PI
      PARAMETER       (PI=3.14159265358979323846D+00)
      CHARACTER        CSTNT*12
      CHARACTER        CSTUT*12
*----
*  Various constants
*---- 
      XDRCST=0.0D0
      CSTNT=CSTNAM
      CSTUT=CSTUNT
      IF(INDEX(CSTNT,'Avogadro') .NE. 0) THEN
        XDRCST=6.0221415D+23
        IF(INDEX(CSTUT,'N/moles') .EQ. 0) THEN
          CALL XABORT(NAMSBR//': Invalid units for Avogadro number')
        ENDIF
      ELSE IF(INDEX(CSTNT,'Plank') .NE. 0) THEN
        XDRCST=6.6260693D-34
        IF(INDEX(CSTUT,'MeV s') .NE. 0) THEN
          XDRCST=1.0D-06*XDRCST/EVJOULE
        ELSE IF(INDEX(CSTUT,'eV s') .NE. 0) THEN
          XDRCST=XDRCST/EVJOULE
        ELSE IF(INDEX(CSTUT,'J s') .EQ. 0) THEN
          CALL XABORT(NAMSBR//': Invalid units for Plank constant')
        ENDIF
      ELSE IF(INDEX(CSTNT,'Boltzmann') .NE. 0) THEN
        XDRCST=1.3806505D-23
        IF(INDEX(CSTUT,'MeV/K') .NE. 0) THEN
          XDRCST=1.0D-06*XDRCST/EVJOULE
        ELSE IF(INDEX(CSTUT,'eV/K') .NE. 0) THEN
          XDRCST=XDRCST/EVJOULE
        ELSE IF(INDEX(CSTUT,'J/K') .EQ. 0) THEN
          CALL XABORT(NAMSBR//': Invalid units for Boltzmann constant')
        ENDIF
*----
*  Various mass
*----
      ELSE IF(INDEX(CSTNT,'Neutron mass') .NE. 0) THEN
        XDRCST=1.67492728D-27
        IF(INDEX(CSTUT,'amu') .NE. 0) THEN
          XDRCST=XDRCST/AMUKG
        ELSE IF(INDEX(CSTUT,'MeV') .NE. 0) THEN
          XDRCST=1.0D-06*XDRCST*SPEEDL*SPEEDL/EVJOULE
        ELSE IF(INDEX(CSTUT,'eV') .NE. 0) THEN
          XDRCST=XDRCST*SPEEDL*SPEEDL/EVJOULE
        ELSE IF(INDEX(CSTUT,'kg') .EQ. 0) THEN
          CALL XABORT(NAMSBR//': Invalid units for neutron mass')
        ENDIF
      ELSE IF(INDEX(CSTNT,'Proton mass') .NE. 0) THEN
        XDRCST=1.67262171D-27
        IF(INDEX(CSTUT,'amu') .NE. 0) THEN
          XDRCST=XDRCST/AMUKG
        ELSE IF(INDEX(CSTUT,'MeV') .NE. 0) THEN
          XDRCST=1.0D-06*XDRCST*SPEEDL*SPEEDL/EVJOULE
        ELSE IF(INDEX(CSTUT,'eV') .NE. 0) THEN
          XDRCST=XDRCST*SPEEDL*SPEEDL/EVJOULE
        ELSE IF(INDEX(CSTUT,'kg') .EQ. 0) THEN
          CALL XABORT(NAMSBR//': Invalid units for neutron mass')
        ENDIF
*----
*  Various mass energy conversion units
*----
      ELSE IF(INDEX(CSTNT,'kg') .NE. 0) THEN
        IF(INDEX(CSTUT,'amu') .NE. 0) THEN
          XDRCST=1.0D0/AMUKG
        ELSE IF(INDEX(CSTUT,'MeV') .NE. 0) THEN
          XDRCST=1.0D-06*SPEEDL*SPEEDL/EVJOULE
        ELSE IF(INDEX(CSTUT,'eV') .NE. 0) THEN
          XDRCST=SPEEDL*SPEEDL/EVJOULE
        ELSE IF(INDEX(CSTUT,'J') .NE. 0) THEN
          XDRCST=SPEEDL*SPEEDL
        ELSE
          CALL XABORT(NAMSBR//': No relations between '//CSTNT//
     >                        ' and '//CSTUT)
        ENDIF
      ELSE IF(INDEX(CSTNT,'amu') .NE. 0) THEN
        IF(INDEX(CSTUT,'kg') .NE. 0) THEN
          XDRCST=AMUKG
        ELSE IF(INDEX(CSTUT,'MeV') .NE. 0) THEN
          XDRCST=1.0D-06*AMUKG*SPEEDL*SPEEDL/EVJOULE
        ELSE IF(INDEX(CSTUT,'eV') .NE. 0) THEN
          XDRCST=AMUKG*SPEEDL*SPEEDL/EVJOULE
        ELSE IF(INDEX(CSTUT,'J') .NE. 0) THEN
          XDRCST=AMUKG*SPEEDL*SPEEDL
        ELSE
          CALL XABORT(NAMSBR//': No relations between '//CSTNT//
     >                        ' and '//CSTUT)
        ENDIF
      ELSE IF(INDEX(CSTNT,'eV') .NE. 0) THEN
        IF(INDEX(CSTUT,'J') .NE. 0) THEN
          XDRCST=EVJOULE
        ELSE IF(INDEX(CSTUT,'K') .NE. 0) THEN
          XDRCST=EVJOULE/1.3806505D-23
        ELSE
          CALL XABORT(NAMSBR//': No relations between '//CSTNT//
     >                        ' and '//CSTUT)
        ENDIF
      ELSE IF(INDEX(CSTNT,'K') .NE. 0) THEN
        IF(INDEX(CSTUT,'J') .NE. 0) THEN
          XDRCST=1.3806505D-23
        ELSE IF(INDEX(CSTUT,'eV') .NE. 0) THEN
          XDRCST=1.3806505D-23/EVJOULE
        ELSE
          CALL XABORT(NAMSBR//': No relations between '//CSTNT//
     >                        ' and '//CSTUT)
        ENDIF
      ELSE IF(INDEX(CSTNT,'J') .NE. 0) THEN
        IF(INDEX(CSTUT,'eV') .NE. 0) THEN
          XDRCST=1.0D0/EVJOULE
        ELSE IF(INDEX(CSTUT,'K') .NE. 0) THEN
          XDRCST=1.0D0/1.3806505D-23
        ELSE
          CALL XABORT(NAMSBR//': No relations between '//CSTNT//
     >                        ' and '//CSTUT)
        ENDIF
      ELSE IF(INDEX(CSTNT,'Pi') .NE. 0) THEN
        IF(CSTUT .EQ. ' ') THEN
          XDRCST=PI
        ELSE
          CALL XABORT(NAMSBR//': No units for Pi ')
        ENDIF
      ELSE
        CALL XABORT(NAMSBR//': '//CSTNT//
     >             ' is an invalid constant or unit') 
      ENDIF
      RETURN
      END
