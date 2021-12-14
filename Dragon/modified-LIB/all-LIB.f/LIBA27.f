*DECK LIBA27
      SUBROUTINE LIBA27(NAMFIL,NBISO,NISOT,NSEGM,NL,ISONRF,ISHINA,
     1 MASKI,NOM,NOMOB,IPR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Probe the APOLIB-2 file and compute the IPR main index vector.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NAMFIL  name of the APOLIB-2 file.
* NBISO   number of isotopes present in the calculation domain.
* NISOT   number of isotopes in the PHEAD segment.
* NSEGM   number of APOLIBE segments in the APOLIB-2 file.
* NL      number of Legendre orders required in the calculation
*         NL=1 or higher.
* ISONRF  library name of isotopes.
* ISHINA  self shielding name.
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* NOM     isotope names in the PHEAD segment.
* NOMOB   APOLIBE segment names.
*
*Parameters: output
* IPR     main index vector:
*         IPR(1,I)  index in PHEAD segment table;
*         IPR(2,I)  segment index of main data (ISOTOP);
*         IPR(3,I)  segment index of production data (PHYSIQ);
*         IPR(4,I)  segment index of delayed neutron data (BETAEF);
*         IPR(5,I)  segment index of main ss data (SSDATA);
*         IPR(6,I)  segment index of 104 flux data (SSPOND);
*         IPR(7,I)  segment index of autolib data (SSSECT);
*         IPR(..,I) segment index of pn diff xs data (DIFF..);
*         IPR(..,I) segment index of pn transfer data (TRAN..).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER NAMFIL*(*)
      INTEGER NBISO,NISOT,NSEGM,NL,ISONRF(3,NBISO),ISHINA(3,NBISO),
     1 NOM(5,NISOT),NOMOB(7,NSEGM),IPR(7+2*(NL-1),NBISO)
      LOGICAL MASKI(NBISO)
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT20*20,HNISOR*12,HNISSS*12,HSMG*131
      INTEGER NITCA(5)
*
      CALL XDISET(IPR,(7+2*(NL-1))*NBISO,0)
      DO 200 IMX=1,NBISO
      IF(MASKI(IMX)) THEN
         WRITE(HNISOR,'(3A4)') (ISONRF(I0,IMX),I0=1,3)
         WRITE(HNISSS,'(3A4)') (ISHINA(I0,IMX),I0=1,3)
         CALL LCMCAR(HNISOR,.TRUE.,NITCA)
         KISO=0
         DO 10 ISO=1,NISOT
         IF(NITCA(1).EQ.NOM(1,ISO)) THEN
            IF(NITCA(2).EQ.NOM(2,ISO)) THEN
               IF(NITCA(3).EQ.NOM(3,ISO)) THEN
                  KISO=ISO
                  GO TO 20
               ENDIF
            ENDIF
         ENDIF
   10    CONTINUE
         WRITE (HSMG,300) HNISOR,NAMFIL
         CALL XABORT(HSMG)
   20    IPR(1,IMX)=KISO
*
         TEXT20='ISOTOP'//HNISOR(:12)
         CALL LCMCAR(TEXT20,.TRUE.,NITCA)
         KISEG=0
         DO 30 ISEG=1,NSEGM
         IF(NITCA(1).EQ.NOMOB(1,ISEG)) THEN
            IF(NITCA(2).EQ.NOMOB(2,ISEG)) THEN
               IF(NITCA(3).EQ.NOMOB(3,ISEG)) THEN
                  IF(NITCA(4).EQ.NOMOB(4,ISEG)) THEN
                     IF(NITCA(5).EQ.NOMOB(5,ISEG)) THEN
                        KISEG=ISEG
                        GO TO 40
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
   30    CONTINUE
         WRITE (HSMG,300) HNISOR,NAMFIL
         CALL XABORT(HSMG)
   40    IPR(2,IMX)=KISEG
*
         TEXT20='PHYSIQ'//HNISOR(:12)
         CALL LCMCAR(TEXT20,.TRUE.,NITCA)
         KISEG=0
         DO 50 ISEG=1,NSEGM
         IF(NITCA(1).EQ.NOMOB(1,ISEG)) THEN
            IF(NITCA(2).EQ.NOMOB(2,ISEG)) THEN
               IF(NITCA(3).EQ.NOMOB(3,ISEG)) THEN
                  IF(NITCA(4).EQ.NOMOB(4,ISEG)) THEN
                     IF(NITCA(5).EQ.NOMOB(5,ISEG)) THEN
                        KISEG=ISEG
                        GO TO 60
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
   50    CONTINUE
   60    IPR(3,IMX)=KISEG
*
         TEXT20='BETAEF'//HNISOR(:12)
         CALL LCMCAR(TEXT20,.TRUE.,NITCA)
         KISEG=0
         DO 70 ISEG=1,NSEGM
         IF(NITCA(1).EQ.NOMOB(1,ISEG)) THEN
            IF(NITCA(2).EQ.NOMOB(2,ISEG)) THEN
               IF(NITCA(3).EQ.NOMOB(3,ISEG)) THEN
                  IF(NITCA(4).EQ.NOMOB(4,ISEG)) THEN
                     IF(NITCA(5).EQ.NOMOB(5,ISEG)) THEN
                        KISEG=ISEG
                        GO TO 80
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
   70    CONTINUE
   80    IPR(4,IMX)=KISEG
*
         IF(HNISSS.NE.' ') THEN
            TEXT20='SSDATA'//HNISSS(:12)
            CALL LCMCAR(TEXT20,.TRUE.,NITCA)
            KISEG=0
            DO 90 ISEG=1,NSEGM
            IF(NITCA(1).EQ.NOMOB(1,ISEG)) THEN
               IF(NITCA(2).EQ.NOMOB(2,ISEG)) THEN
                  IF(NITCA(3).EQ.NOMOB(3,ISEG)) THEN
                     IF(NITCA(4).EQ.NOMOB(4,ISEG)) THEN
                        IF(NITCA(5).EQ.NOMOB(5,ISEG)) THEN
                           KISEG=ISEG
                           GO TO 100
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
   90       CONTINUE
            WRITE (HSMG,310) HNISSS,NAMFIL
            CALL XABORT(HSMG)
  100       IPR(5,IMX)=KISEG
*
            TEXT20='SSPOND'//HNISSS(:12)
            CALL LCMCAR(TEXT20,.TRUE.,NITCA)
            KISEG=0
            DO 110 ISEG=1,NSEGM
            IF(NITCA(1).EQ.NOMOB(1,ISEG)) THEN
               IF(NITCA(2).EQ.NOMOB(2,ISEG)) THEN
                  IF(NITCA(3).EQ.NOMOB(3,ISEG)) THEN
                     IF(NITCA(4).EQ.NOMOB(4,ISEG)) THEN
                        IF(NITCA(5).EQ.NOMOB(5,ISEG)) THEN
                           KISEG=ISEG
                           GO TO 120
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
  110       CONTINUE
  120       IPR(6,IMX)=KISEG
*
            TEXT20='SSSECT'//HNISSS(:12)
            CALL LCMCAR(TEXT20,.TRUE.,NITCA)
            KISEG=0
            DO 130 ISEG=1,NSEGM
            IF(NITCA(1).EQ.NOMOB(1,ISEG)) THEN
               IF(NITCA(2).EQ.NOMOB(2,ISEG)) THEN
                  IF(NITCA(3).EQ.NOMOB(3,ISEG)) THEN
                     IF(NITCA(4).EQ.NOMOB(4,ISEG)) THEN
                        IF(NITCA(5).EQ.NOMOB(5,ISEG)) THEN
                           KISEG=ISEG
                           GO TO 140
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
  130       CONTINUE
  140       IPR(7,IMX)=KISEG
         ENDIF
*
         DO 190 IL=2,NL
         WRITE(TEXT20,'(4HDIFF,I2.2,A12)') IL-1,HNISOR
         CALL LCMCAR(TEXT20,.TRUE.,NITCA)
         KISEG=0
         DO 150 ISEG=1,NSEGM
         IF(NITCA(1).EQ.NOMOB(1,ISEG)) THEN
            IF(NITCA(2).EQ.NOMOB(2,ISEG)) THEN
               IF(NITCA(3).EQ.NOMOB(3,ISEG)) THEN
                  IF(NITCA(4).EQ.NOMOB(4,ISEG)) THEN
                     IF(NITCA(5).EQ.NOMOB(5,ISEG)) THEN
                        KISEG=ISEG
                        GO TO 160
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  150    CONTINUE
  160    IPR(7+(IL-1),IMX)=KISEG
*
         WRITE(TEXT20,'(4HTRAN,I2.2,A12)') IL-1,HNISOR
         CALL LCMCAR(TEXT20,.TRUE.,NITCA)
         KISEG=0
         DO 170 ISEG=1,NSEGM
         IF(NITCA(1).EQ.NOMOB(1,ISEG)) THEN
            IF(NITCA(2).EQ.NOMOB(2,ISEG)) THEN
               IF(NITCA(3).EQ.NOMOB(3,ISEG)) THEN
                  IF(NITCA(4).EQ.NOMOB(4,ISEG)) THEN
                     IF(NITCA(5).EQ.NOMOB(5,ISEG)) THEN
                        KISEG=ISEG
                        GO TO 180
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  170    CONTINUE
  180    IPR(7+(NL-1)+(IL-1),IMX)=KISEG
  190    CONTINUE
      ENDIF
  200 CONTINUE
      RETURN
*
  300 FORMAT(26HLIBA27: MATERIAL/ISOTOPE ',A12,20H' IS MISSING ON APOL,
     1 15HIB-2 FILE NAME ,A12,1H.)
  310 FORMAT(49HLIBA27: SELF-SHIELDING DATA OF MATERIAL/ISOTOPE ',A12,
     1 35H' IS MISSING ON APOLIB-2 FILE NAME ,A12,1H.)
      END
