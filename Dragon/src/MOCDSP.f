*DECK MOCDSP
      SUBROUTINE MOCDSP(N,NFI,NLONG,LC,NZON,NOM,KM,MCU,IM,PREV,NEXT,H)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the position of the coefficients relative to a track
* in ACA matrices. Cyclic tracking version.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): I. Suslov and R. Le Tellier
*
*Parameters: input
* NFI     total number of volumes and surfaces.
* NLONG   total number of cells with unknowns quantities.
* LC      dimension of vector MCU.
* NZON    index-number of the mixture type assigned to each volume.
* KM      used in CDD acceleration.
* MCU     used in CDD acceleration.
* IM      used in CDD acceleration.
*
*Parameters: output
* PREV    PREV(I): location of non diagonal element (NOM(I),NOM(I-1))
*                  of preconditioning matrices in vector CF and CQ.
* NEXT    NEXT(I): location of non diagonal element (NOM(I),NOM(I+1))
*                  of preconditioning matrices in vector CF and CQ.
*
*Parameters: input/output
* N       number of elements in the current track.
* NOM     integer tracking elements.
* H       real tracking elements.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,NFI,NLONG,LC,NZON(NFI),NOM(N),KM(NLONG),MCU(LC),
     1 IM(NLONG),PREV(N),NEXT(N)
      DOUBLE PRECISION H(N)
*---- 
*  LOCAL VARIABLES
*----
      INTEGER I,I1,I2,II,FOUN,FOUP,NOMI,NOMIN,NOMIP,NZI
*----
*  UNFOLD CYCLIC TRACKING LINE
*----
      CALL MCGTRK(NFI,NZON,N,NOM,H)
*----
*  CONSTRUCT PREV & NEXT
*----
      DO I=1,N
         FOUN=-1
         FOUP=-1
*        current cell
         NOMI=NOM(I)
         IF((NOMI.LT.1).OR.(NOMI.GT.NFI)) THEN
            WRITE(6,*) 'NFI= ',NFI,' NOMI= ',NOMI
            CALL XABORT('MOCDSP: KM OVERFLOW.')
         ENDIF
         NZI=NZON(NOMI)
         IF (NZI.LT.0) THEN
            NEXT(I)=0
            PREV(I)=0
         ELSE
            I1=IM(NOMI)
*           next cell 
            IF (I.EQ.N) THEN
               NOMIN=NOM(1)
            ELSE
               NOMIN=NOM(I+1)
            ENDIF
*           previous cell
            IF(I.EQ.1) THEN
               NOMIP=NOM(N)
            ELSE
               NOMIP=NOM(I-1)
            ENDIF
            IF ((NOMI.EQ.NOMIN).OR.(NZON(NOMIN).LT.0)) THEN
               FOUN=0
            ENDIF
            IF ((NOMI.EQ.NOMIP).OR.(NZON(NOMIP).LT.0)) THEN
               FOUP=0
            ENDIF
*
            I2=I1+KM(NOMI)
            I1=I1+1
            DO II=I1,I2
               IF ((FOUN.LT.0).AND.(MCU(II).EQ.NOMIN)) THEN
                  FOUN=II
               ENDIF
               IF ((FOUP.LT.0).AND.(MCU(II).EQ.NOMIP)) THEN
                  FOUP=II
               ENDIF
               IF ((FOUN.GE.0).AND.(FOUP.GE.0)) GOTO 10
            ENDDO
*           connectivity between NOMI and NOMIN and/or NOMIP not found
            WRITE(6,100) I,NOMI,NOMIN,NOMIP
            CALL PRINIM('NOM   ',NOM,N)
            CALL PRINIM('MCU   ',MCU(I1),KM(NOMI))
            CALL XABORT('MOCDSP: FAILURE 1.')
 10         IF ((FOUN.LE.LC).AND.(FOUP.LE.LC)) THEN
               PREV(I)=FOUP
               NEXT(I)=FOUN
            ELSE             
               CALL XABORT('MOCDSP: CQ/CF OVERFLOW.')
            ENDIF
         ENDIF
      ENDDO
*
 100  FORMAT(1X,'I=',I3,' NOMI=',I5,' NOMIN=',I5,' NOMIP=',I5)
*
      RETURN
      END
