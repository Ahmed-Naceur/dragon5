*DECK INFWAT
      SUBROUTINE INFWAT(TEMPC,PURWGT,DENSTY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute water density as a function of temperature and pressure.
*
*Copyright:
* Copyright (C) 1995 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): R. Roy 
*
*Parameters: input
* TEMPC   temperature (celcius).
* PURWGT  D2O purity (in wgt%).
*
*Parameters: output
* DENSTY  density (G/CM**3).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*
*     TEMPERATURE DEPENDENT DENSITY CALCULATION IS DONE ACCORDING TO:
*
*         DENSITY = D(H2O) * D(D2O) /
*                   ( WGT%(H2O) * D(D2O) + WGT%(D2O) * D(H2O) )
*
      REAL TEMPC,  PURWGT, DENSTY
      REAL DEND2O, DENH2O, WGTD2O, WGTH2O
      REAL AIKINT
      REAL TDND2O(20), TDNH2O(85), TMPD2O(20), TMPH2O(85)
*
*------------------------------------------------------------------
*
* >>> TABLES ORIGINALLY CAME FROM ROUTINE dmats.f (WIMS)
*
*
* * D2O DATA CONSISTENT WITH AECL 7531, TABLE FROM J. PHYS. CHEM. REF.
*      DATA, VOL11, NO.1, 1982, P6 (SAME AUTHORS)
*
       DATA TMPD2O /3.8, 6.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0,
     >      49.99, 100.0, 111.02, 150.02, 200.0, 250.0, 275.,
     >     300., 325., 350.127, 360.057/
*
       DATA TDND2O / 0.90464, 0.90439, 0.90419, 0.90428, 0.90472,
     >     0.90545, 0.90645, 0.90771, 0.90918, 0.91274, 0.94057,
     >     0.94866, 0.98296, 1.04354, 1.13149, 1.19270, 1.2740,
     >     1.3917, 1.6044, 1.7709 /
*
* * H2O DATA FROM CRNL-1533
*
       DATA TMPH2O /3.98,20.0,30.0,40.0,50.0,
     >     60.0,70.0,80.0,90.0,99.6,
     >     120.2,133.5,143.6,151.8,158.8,165.0,
     >     170.4,175.4,179.9,188.0,195.0,201.4,
     >     207.1,212.4,217.2,221.8,226.0,230.0,
     >     233.8,237.4,240.9,244.2,247.3,250.3,
     >     253.2,256.0,258.8,261.4,263.9,266.4,
     >     268.8,271.1,273.3,275.6,277.7,279.8,
     >     281.8,283.8,285.8,287.7,289.6,291.4,
     >     293.2,295.0,296.7,298.4,300.1,301.7,
     >     303.3,307.2,311.0,314.6,318.0,321.4,
     >     324.6,327.8,330.8,333.8,336.6,339.4,
     >     342.1,344.8,347.3,349.8,352.3,354.6,
     >     357.0,359.2,361.4,363.6,365.7,367.8,
     >     369.8,371.8,373.7/
*
       DATA TDNH2O /0.999973,0.998418,0.995848,0.992385,0.988164,
     >     0.983274,0.977785,0.971753,0.965218,0.958479,
     >     0.942737,0.931621,0.922712,0.915136,0.908466,0.902459,
     >     0.896961,0.891869,0.887107,0.878373,0.870456,0.863165,
     >     0.856371,0.849982,0.843931,0.838165,0.832644,0.827336,
     >     0.822215,0.817259,0.812449,0.807771,0.803211,0.798756,
     >     0.794401,0.790133,0.785946,0.781833,0.777787,0.773804,
     >     0.769879,0.766007,0.762183,0.758405,0.754669,0.750972,
     >     0.747310,0.743682,0.740085,0.736516,0.732973,0.729455,
     >     0.725958,0.722482,0.719025,0.715585,0.712160,0.708750,
     >     0.705352,0.696902,0.688503,0.680135,0.671780,0.663420,
     >     0.655039,0.646619,0.638141,0.629586,0.620932,0.612153,
     >     0.603220,0.594098,0.584739,0.575085,0.565162,0.554167,
     >     0.543564,0.531796,0.519245,0.505736,0.490945,0.474241,
     >     0.454256,0.427351,0.374304/
*
*------------------------------------------------------------------
*
*     COMPUTE DENSITIES FOR LIGHT AND HEAVY WATER (PURE)
*      AT SUCH TEMPERATURE:
*
*      1. USE  DIRECT INTERPOLATION FOR LIGHT WATER:
*
          DENH2O =     AIKINT(TEMPC,TMPH2O,TDNH2O,85,1.0E-5)
*
*      2. INVERSE THE INTERPOLATION FOR HEAVY WATER:
*
          DEND2O = 1.0/AIKINT(TEMPC,TMPD2O,TDND2O,20,1.0E-5)
*
      WGTD2O = 0.01 * PURWGT
      WGTH2O = 1.00 - WGTD2O
*
*     COMPUTE GLOBAL DENSITY A MIX OF LIGHT AND HEAVY WATER:
*
      DENSTY = DENH2O * DEND2O /( WGTH2O * DEND2O + WGTD2O * DENH2O )
*
      RETURN
      END
