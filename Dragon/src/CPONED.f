*DECK CPONED
      SUBROUTINE CPONED(NPROC ,MINLEG,MAXLEG,ILEAKS ,NED   ,HVECT ,
     >                  IVECT ,INDPRO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set up INDPRO for cross section to read on IPLIB.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input
* NPROC   number of terms to process.
* MINLEG  mimimum Legendre order to process for scattering.
* MAXLEG  maximum Legendre order to process for scattering.
* ILEAKS  leakage calculation: = 1 STRD; = 2 STRDX, STRDY and STRDZ.
* NED     number of extra vector edits.
* HVECT   names of the extra vector edits.
*
*Parameters: output
* IVECT   pointer to additional xs possible.
* INDPRO  vector for cross section to process:
*           = 0 do not process;
*           > 0 process.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          NPROC ,MINLEG,MAXLEG,ILEAKS,NED,IVECT(NED),
     >                 INDPRO(NPROC)
      CHARACTER        HVECT(NED)*8
*----
*  LOCAL PARAMETERS
*  NDPROC = NUMBER OF DEFAULT CROSS SECTIONS = 20
*  NAMDXS = NAME OF NDPROC DEFAULT XS
*  SCATTERING CROSS SECTIONS START AT NDPROC+1 WITH
*  NAME NAMSCT='SIGS'//NAMLEG AND NAMSCT='SCAT'//NAMLEG
*  WITH NAMLEG DEFINED BY
*  WRITE(NAMLEG ,'(I2.2)') ILEG
*  FOR ILEG=0 TO NDPROC-NPROC-1
*----
      INTEGER          NDPROC,IOUT,NEDOTH,IED,IXSR
      PARAMETER       (NDPROC=20,IOUT=6)
      CHARACTER        NAMDXS(NDPROC)*6
      SAVE             NAMDXS
      DATA    NAMDXS  /'NTOT0 ','TRANC ','NUSIGF','NFTOT ','CHI   ',
     >                 'NU    ','NG    ','NHEAT ','N2N   ','N3N   ',
     >                 'N4N   ','NP    ','NA    ','GOLD  ','ABS   ',
     >                 'NWT0  ','STRD  ','STRD X','STRD Y','STRD Z'/
*----
*  SCAN FOR ADDITIONAL AND STANDARD CROSS SECTIONS TO BE SAVED
*----
      CALL XDRSET(IVECT,NED,0)
      CALL XDRSET(INDPRO,NPROC,0)
      NEDOTH=NED
      DO 100 IED=1,NED
        IF(HVECT(IED).EQ.'        ') THEN
          NEDOTH=NEDOTH-1
        ELSE
          DO 110 IXSR=1,NDPROC
            IF(HVECT(IED)(:6).EQ.NAMDXS(IXSR)) THEN
              NEDOTH=NEDOTH-1
              INDPRO(IXSR)=1
              IF(HVECT(IED).EQ.'NFTOT') GO TO 115
              IVECT(IED)=IXSR
              GO TO 115
            ENDIF
 110      CONTINUE
 115      CONTINUE
        ENDIF
 100  CONTINUE
      IF(NEDOTH.GE.1) THEN
        WRITE(IOUT,9000)
        DO 120 IED=1,NED
          IF(IVECT(IED).EQ.0.AND.
     >       HVECT(IED).NE.'NFTOT'.AND.HVECT(IED).NE.' ') THEN
            WRITE(IOUT,9001) HVECT(IED)
          ENDIF
 120    CONTINUE
        WRITE(IOUT,9002)
      ENDIF
      DO 130 IXSR=1,7
        INDPRO(IXSR)=1
 130  CONTINUE
      INDPRO(16)=1
      IF(ILEAKS.EQ.1) THEN
        INDPRO(17)=1
      ELSE IF(ILEAKS.EQ.2) THEN
        INDPRO(18)=1
        INDPRO(19)=1
        INDPRO(20)=1
      ENDIF
      DO 140 IXSR=NDPROC+MINLEG+1,NDPROC+MAXLEG+1
        INDPRO(IXSR)=1
 140  CONTINUE
      RETURN
*----
*  FORMAT
*----
 9000 FORMAT(' CPONED: ************ WARNING ************')
 9001 FORMAT(' CROSS-SECTION TYPE NOT RECOVERED : ',A8)
 9002 FORMAT(' *****************************************')
      END
