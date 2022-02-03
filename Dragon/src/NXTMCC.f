*DECK NXTMCC
      SUBROUTINE NXTMCC(IPTRK,NAMCEL,NREGC,NSURC,NREGF,NSURF,INDEX,
     1                  IDSUR,IDREG)
*-----------------------------------------------------------------------
*
*Purpose:
* Calculate and store the compressed index and region/surface ids for
* an elementary geometry.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): Romain Le Tellier
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* NAMCEL  name of the elementary geometry to be treated.
* NREGC   number of regions (uncompressed).
* NSURC   number of surfaces (uncompressed).
*
*Parameters: output
* NREGF   number of regions (compressed).
* NSURF   number of surfaces (compressed).
*
*Parameters: input/output
* INDEX   index vector (uncompressed and compressed).
* IDSUR   surface identificator (uncompressed and compressed).
* IDREG   region identificator (uncompressed and compressed).
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPTRK
*      INTEGER IPTRK
      INTEGER NREGC,NSURC,NREGF,NSURF,INDEX(5,-NSURC:NREGC,2),
     1 IDSUR(NSURC,2),IDREG(NREGC,2)
      CHARACTER NAMCEL*9
*----
*  LOCAL VARIABLES
*----
      INTEGER NSTATE,IOUT
      PARAMETER(NSTATE=40,IOUT=6)
      INTEGER ESTATE(NSTATE)
      INTEGER I,ISUR,INDF,JJ,ITMP,IREG
      CHARACTER NAMREC*12
*----
*  SCAN THE SURFACES AND FILL IN THE SURFACE ID AND CORREPONDING INDEX
*----
      NSURF=0
      INDF=-NSURC-1
      DO I=NSURC,1,-1
         ISUR=IDSUR(I,1)
         IF (IDSUR(I,1).NE.0) THEN
            NSURF=NSURF+1
            IDSUR(NSURF,2)=ABS(ISUR)
            INDF=INDF+1
            DO JJ=1,4
               INDEX(JJ,INDF,2)=INDEX(JJ,-I,1)
            ENDDO
         ENDIF
      ENDDO
*----
*  REVERSE SURFACE ID IN SUCH A WAY THAT
*  IDSUR(I,2) CORRESPONDS TO INDEX(:,-NSURC+NSURF-I,2)
*----
      DO I=1,NSURF/2
         ITMP=IDSUR(NSURF+1-I,2)
         IDSUR(NSURF+1-I,2)=IDSUR(I,2)
         IDSUR(I,2)=ITMP
      ENDDO
      INDF=INDF+1
      DO JJ=1,4
         INDEX(JJ,INDF,2)=0
      ENDDO
*----
*  SCAN THE REGIONS AND FILL IN THE SURFACE ID AND CORREPONDING INDEX
*----
      NREGF=0
      DO I=1,NREGC
         IREG=IDREG(I,1)
         IF (IDREG(I,1).NE.0) THEN
            NREGF=NREGF+1
            IDREG(NREGF,2)=ABS(IREG)
            INDF=INDF+1
            DO JJ=1,4
               INDEX(JJ,INDF,2)=INDEX(JJ,I,1)
            ENDDO
         ENDIF
      ENDDO
*----
*  STORE THE FINAL NUMBER OF REGIONS/SURFACES
*  AND THE COMPRESSED IDS AND INDEX
*----
      NAMREC=NAMCEL//'DIM'
      CALL XDISET(ESTATE,NSTATE,0)
      CALL LCMGET(IPTRK,NAMREC,ESTATE)
      ESTATE(39)=NREGF
      ESTATE(40)=NSURF
      CALL LCMPUT(IPTRK,NAMREC,NSTATE,1,ESTATE)
      IF (NREGF.GT.0) THEN
         NAMREC=NAMCEL//'RIC'
         CALL LCMPUT(IPTRK,NAMREC,NREGF,1,IDREG(1,2))
      ENDIF
      IF (NSURF.GT.0) THEN
         NAMREC=NAMCEL//'SIC'
         CALL LCMPUT(IPTRK,NAMREC,NSURF,1,IDSUR(1,2))
      ENDIF
      INDF=NREGF+NSURF+1
      IF (INDF.GT.0) THEN
         NAMREC=NAMCEL//'VSC'
         CALL LCMPUT(IPTRK,NAMREC,5*INDF,1,INDEX(1,-NSURC,2))
      ENDIF
*
      RETURN
      END
