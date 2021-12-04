*DECK XDRCRE
      SUBROUTINE XDRCRE(NAMMOD,IBEAF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To print DRAGON credits.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* NAMMOD  name of DRAGON module.
* IBEAF   flag for beginning or finishing module where:
*         =1  before module execution;
*         =-1 after module execution.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER        NAMMOD*12
      INTEGER          IBEAF
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='XDRCRE')
*----
*  LOCAL VARIABLES
*----
      CHARACTER        USE*72,AUT*72
      INTEGER          ICOPYR
*----
*  PRINT CREDITS
*----
      ICOPYR=0
      IF(IBEAF .EQ. 1) THEN
        WRITE(IOUT,6000) NAMMOD
        IF     (NAMMOD .EQ.  'ASM:        ') THEN      
          USE='To built system matrices (CP and IC)'
          AUT='A. Hebert, G. Marleau, R. Roy'
        ELSE IF(NAMMOD .EQ.  'COMPO:      ') THEN      
          USE='Create multiparameter reactor composition database'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'EDI:        ') THEN      
          USE='Editing module'
          AUT='A. Hebert, G. Marleau'
        ELSE IF(NAMMOD .EQ.  'EVO:        ') THEN      
          USE='Isotopic depletion and fuel burnup'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'EXCELT:     ') THEN      
          USE='Excell tracking in 2- and 3-D'
          AUT='G. Marleau, M. Ouisloumen, R. Roy'
        ELSE IF(NAMMOD .EQ.  'NXT:        ') THEN      
          USE='New excell tracking in 2- and 3-D'
          AUT='G. Marleau'
        ELSE IF(NAMMOD .EQ.  'MCCGT:      ') THEN      
          USE='Method of characteristics in 2- and 3-D'
          AUT='I. Suslov, R. Roy, R. Le Tellier, A. Hebert'
        ELSE IF(NAMMOD .EQ.  'BIVACT:     ') THEN      
          USE='2-D diffusion or SPN finite element tracking'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'TRIVAT:     ') THEN      
          USE='3-D diffusion or SPN finite element tracking'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'SNT:        ') THEN      
          USE='Discrete ordinates tracking'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'FLU:        ') THEN      
          USE='Solve the flux equations'
          AUT='R. Roy, A. Hebert, G. Marleau'
        ELSE IF(NAMMOD .EQ.  'GEO:        ') THEN      
          USE='Geometry definition'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'INFO:       ') THEN      
          USE='Information on water, UO2 and ThO2'
          AUT='R. Roy'
        ELSE IF(NAMMOD .EQ.  'LIB:        ') THEN      
          USE='Microscopic xs-library processing'
          AUT='A. Hebert, G. Marleau'
        ELSE IF(NAMMOD .EQ.  'MAC:        ') THEN      
          USE='Macroscopic xs processor'
          AUT='G. Marleau'
        ELSE IF(NAMMOD .EQ.  'MRG:        ') THEN      
          USE='Merge excell tracking file'
          AUT='G. Marleau'
        ELSE IF(NAMMOD .EQ.  'PSP:        ') THEN      
          USE='Generates ps graphics for dragon'
          AUT='K.E. Kohler, G. Marleau'
          ICOPYR=1
        ELSE IF(NAMMOD .EQ.  'SHI:        ') THEN      
          USE='Self-shielding by improved Stammler method'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'USS:        ') THEN      
          USE='Self-shielding by subgroup method'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'TONE:       ') THEN      
          USE='Self-shielding by Tone method'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'SYBILT:     ') THEN      
          USE='Sybil 2-D tracking'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'TLM:        ') THEN      
          USE='Generate matlab line tracking file'
          AUT='C. Plamondon, G. Marleau'
        ELSE IF(NAMMOD .EQ.  'M2T:        ') THEN      
          USE='Generate an apotrim interface file'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'FMAC:       ') THEN      
          USE='Recover information from a FMAC-M interface file'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'PSOUR:      ') THEN      
          USE='Compute a fixed source from companion particles'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'SOUR:      ') THEN      
          USE='Definie external fixed sources'
          AUT='C. Bienvenue'
        ELSE IF(NAMMOD .EQ.  'HEAT:       ') THEN      
          USE='Compute the energy and charge deposition values'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'CHAB:       ') THEN      
          USE='Modify and renormalize a microlib'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'CPO:        ') THEN      
          USE='Create Version3 reactor composition database'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'SAP:        ') THEN      
          USE='Create a Saphyb composition database'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'MC:        ') THEN
          USE='Multigroup Monte-Carlo calculation'
          AUT='R. Le Tellier, B. Arsenault'
        ELSE IF(NAMMOD .EQ.  'T:        ') THEN
          USE='Transpose a macrolib'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'DMAC:      ') THEN
          USE='Set the GPT adjoint sources (Macrolib gradient)'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'FMT:        ') THEN      
          USE='Transfer data structure to formatted output files'
          AUT='G. Marleau'
        ELSE IF(NAMMOD .EQ.  'EPC:        ') THEN      
          USE='Error propagation module'
          AUT='G. Marleau'
        ELSE IF(NAMMOD .EQ.  'SPH:        ') THEN      
          USE='Superhomogenization (SPH) calculation'
          AUT='A. Hebert'
        ELSE IF(NAMMOD .EQ.  'CFC:        ') THEN      
          USE='Construction of a feedback database for CANDU reactors'
          AUT='M. T. Sissaoui'
        ELSE IF(NAMMOD .EQ.  'SENS:       ') THEN      
          USE='Sensitivity analysis to cross-section on the reactivity'
          AUT='C. Laville, G. Marleau'          
        ELSE IF(NAMMOD .EQ.  'DUO:        ') THEN      
          USE='Perturbative analysis using the Clio formula'
          AUT='A. Hebert'          
        ELSE IF(NAMMOD .EQ.  'BREF:       ') THEN      
          USE='Discontinuity factors calculation in a 1D reflector'
          AUT='A. Hebert'          
        ELSE IF(NAMMOD .EQ.  'SALT:       ') THEN      
          USE='Track calculations from a SALOME surfacic file'
          AUT='A. Hebert, X. Warin'          
        ELSE IF(NAMMOD .EQ.  'G2S:        ') THEN      
          USE='Generate a surfacic file in SALOME format from a DRAGON'
     >    //' geometry'
          AUT='G. Civario'          
        ELSE IF(NAMMOD .EQ.  'G2MC:       ') THEN
          USE='Generate a surfacic file in Monte Carlo format from a '
     >    //'DRAGON geometry'
          AUT='G. Civario'          
        ELSE IF(NAMMOD .EQ.  'MRG:        ') THEN      
          USE='Merge regions in tracking data structure'
          AUT='G. Marleau'          
        ELSE
          USE='No description available for this module'
          AUT='No author provided for this module'
        ENDIF
        WRITE(IOUT,6100) USE,AUT
        IF(ICOPYR .EQ. 1) THEN
          WRITE(IOUT,6111)
        ELSE
          WRITE(IOUT,6110)
        ENDIF
      ELSE
        WRITE(IOUT,6001) NAMMOD
      ENDIF
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT('->@BEGIN MODULE : ',A12)
 6001 FORMAT('->@END MODULE   : ',A12)
 6100 FORMAT('->@DESCRIPTION  : ',A72/
     >       '->@CREDITS      : ',A72)
 6110 FORMAT('->@COPYRIGHTS   : ECOLE POLYTECHNIQUE DE MONTREAL  '/
     >       '                  GNU LESSER GENERAL PUBLIC LICENSE')
 6111 FORMAT('->@COPYRIGHTS   : ECOLE POLYTECHNIQUE DE MONTREAL  '/
     >       '                  GNU LESSER GENERAL PUBLIC LICENSE'/
     >       '                  K.E. KOHLER FOR PSPLOT           ')
      END  
