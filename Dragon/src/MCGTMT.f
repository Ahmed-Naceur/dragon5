*DECK MCGTMT
      SUBROUTINE MCGTMT(NMERG,NTRTMT,NSETMT,NSEG,NSEG0,NOM,NOM0,WEIGHT,
     1                  WEIGHT0,H,H0,LFORC,NTPROC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Track merging.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* NMERG   number of merged tracks for the track under contruction.
* NTRTMT  total number of finalized tracks.
* NSETMT  total number of segments in finalized tracks.
* NSEG    number of segments in the track to be processed.
* NSEG0   number of segments in the track under construction.
* NOM0    integer tracking elements for the under construction.
* WEIGHT0 weight of the under construction.
* H       real tracking elements for the track to be processed.
* LFORC   flag to force a merged track to be finalized.
*
*Parameters: input/output 
* NOM     integer tracking elements for the finalized track.
* WEIGHT  weight of the finalized track.
* H0      real tracking elements for the finalized track.
*
*Parameters: output 
* NTPROC  number of merged tracks for the finalized track.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NMERG,NTRTMT,NSETMT,NSEG,NSEG0,NOM(NSEG),
     1 NOM0(NSEG),NTPROC
      DOUBLE PRECISION WEIGHT,H(NSEG),WEIGHT0,H0(NSEG)
      LOGICAL LFORC
*----
*  LOCAL VARIABLES
*----
      INTEGER II,NTEMP
      DOUBLE PRECISION TEMP
*
      IF (NMERG.GT.0) THEN
         IF ((LFORC).OR.(NSEG.NE.NSEG0)) GOTO 20
         DO II=1,NSEG
            IF (NOM(II).NE.NOM0(II)) GOTO 20
         ENDDO
*        merge this track with the previous one
         DO II=1,NSEG
            H0(II)=H0(II)+WEIGHT*H(II)
         ENDDO
         WEIGHT0=WEIGHT0+WEIGHT
      ELSE
*     start
         DO II=1,NSEG
            H0(II)=WEIGHT*H(II)
            NOM0(II)=NOM(II)
         ENDDO
         WEIGHT0=WEIGHT
      ENDIF
      NMERG=NMERG+1
      NSEG0=NSEG
      NTPROC=0
      RETURN
 20   CONTINUE 
*     finalize this "merged" track and start a new one
      NTRTMT=NTRTMT+1
      NSETMT=NSETMT+NSEG0
      DO II=1,MAX(NSEG0,NSEG)
         TEMP=H(II)
         H(II)=H0(II)/WEIGHT0
         H0(II)=WEIGHT*TEMP
         NTEMP=NOM(II)
         NOM(II)=NOM0(II)
         NOM0(II)=NTEMP
      ENDDO
      NTPROC=NMERG
      NMERG=1
      TEMP=WEIGHT
      WEIGHT=WEIGHT0
      WEIGHT0=TEMP
      NTEMP=NSEG
      NSEG=NSEG0
      NSEG0=NTEMP
      RETURN
      END
