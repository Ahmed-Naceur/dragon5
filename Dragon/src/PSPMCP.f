*DECK PSPMCP
      SUBROUTINE PSPMCP(ISPSP,OFFC,FACT,N,COORD,REGI,EVENT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Add MC: neutron paths to the graphics of a 2-D NXT geometry.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* ISPSP   pointer to the POSTSCRIPT file.
* OFFC    offset vector. 
* FACT    scaling factor.
* N       number of points.
* COORD   points coordinates.
* REGI    regions indexes.
* EVENT   event indexes.
* 
*-----------------------------------------------------------------------
* 
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ISPSP,N,REGI(N),EVENT(N)
      DOUBLE PRECISION COORD(3,N),OFFC(2),FACT
*----
*  LOCAL VARIABLES
*----
      REAL WLINE,HTEX,HCRO
      PARAMETER(WLINE=0.002,HTEX=0.11,HCRO=0.03)
      INTEGER I,IDIR,NCHAR,IREG
      REAL POS(2),POSO(2),SEG(2,2),CENTER(2),RADANG(2,2)
      CHARACTER TEXT*5,FORM*4
      LOGICAL START
      INTEGER IORDER(2)
      DATA IORDER /-2,-1 /
*
      START=.TRUE.
      DO I=1,N
*        CALCULATE POSITION IN GRAPHICS COORDINATES
         POS(1)=REAL(FACT*(COORD(1,I)-OFFC(1)))
         POS(2)=REAL(FACT*(COORD(2,I)-OFFC(2)))
         IF (START) THEN
*        STARTING POINT: DRAW A CIRCLE
            CALL PSMOVE(ISPSP,POS,-3)
            RADANG(1,1)=HCRO
            RADANG(2,1)=0.0
            RADANG(1,2)=HCRO
            RADANG(2,2)=6.30
            CALL PSDRAI(ISPSP,2,IORDER,POS,RADANG)
            CALL PSSTRK(ISPSP,WLINE,0,0)
            CENTER(1)=-POS(1)
            CENTER(2)=-POS(2)
            CALL PSMOVE(ISPSP,CENTER,-3)
         ELSE
*        DRAW A SEGMENT FROM PREVIOUS POINT TO THIS ONE
            DO IDIR=1,2
               SEG(IDIR,1)=POSO(IDIR)
               SEG(IDIR,2)=POS(IDIR)
            ENDDO
            CALL PSDREG(ISPSP,2,SEG)
            CALL PSSTRK(ISPSP,WLINE,0,0)           
         ENDIF
         IF (REGI(I).GT.0) THEN
            IREG=REGI(I)
            START=.FALSE.
         ELSE
            IREG=-REGI(I)
*           ENDING POINT: DRAW A CROSS
            DO IDIR=1,2
               SEG(IDIR,1)=POS(IDIR)-HCRO
               SEG(IDIR,2)=POS(IDIR)+HCRO
            ENDDO
            CALL PSDREG(ISPSP,2,SEG)
            CALL PSSTRK(ISPSP,WLINE,0,0)
            SEG(1,1)=SEG(1,1)+2.0*HCRO
            SEG(1,2)=SEG(1,2)-2.0*HCRO
            CALL PSDREG(ISPSP,2,SEG)
            CALL PSSTRK(ISPSP,WLINE,0,0)
            START=.TRUE.
         ENDIF
*        SAVE PREVIOUS POSITION
         POSO(1)=POS(1)
         POSO(2)=POS(2)
*        INDICATE REGION/SURFACE INDEX
         NCHAR=1
         IF ((IREG.GE.10.).AND.(IREG.LT.100)) THEN
            NCHAR=2
         ELSEIF ((IREG.GE.100.).AND.(IREG.LT.1000)) THEN
            NCHAR=3
         ELSEIF ((IREG.GE.1000.).AND.(IREG.LT.10000)) THEN
            NCHAR=4
         ELSEIF ((IREG.GE.10000.).AND.(IREG.LT.100000)) THEN
            NCHAR=5
         ENDIF
*        WHICH EVENT TOOK PLACE? 
         IF (EVENT(I).LT.0) THEN
*        ENCOUNTERING A SURFACE: indicated by a minus in front of the
*        region index
            NCHAR=NCHAR+1
            IREG=-IREG
            IF (EVENT(I).EQ.-1) THEN
*           X- surface
            POS(1)=POS(1)-0.5*NCHAR*HTEX
            POS(2)=POS(2)-0.5*HTEX
            ELSEIF (EVENT(I).EQ.-2) THEN
*           X+ surface
            POS(1)=POS(1)+0.5*NCHAR*HTEX
            POS(2)=POS(2)-0.5*HTEX
            ELSEIF (EVENT(I).EQ.-3) THEN
*           Y- surface
            POS(2)=POS(2)-1.2*HTEX
            ELSEIF (EVENT(I).EQ.-4) THEN
*           Y+ surface
            POS(2)=POS(2)+0.2*HTEX
            ENDIF
         ELSE
*        INTERACTION IN REGION
*
*        etc ...
*
            POS(2)=POS(2)+0.2*HTEX
         ENDIF
         WRITE(FORM,'(1H(,A1,I1,1H))') 'I',NCHAR
         WRITE(TEXT,FORM) IREG
         CALL PSTEXT(ISPSP,NCHAR,TEXT(1:NCHAR),POS,HTEX,1,0)
      ENDDO
*
      RETURN
      END
