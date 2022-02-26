*DECK MCGTRK
      SUBROUTINE MCGTRK(NFI,NZON,NSEG,NOM,H)
*
*-----------------------------------------------------------------------
*
*Purpose: 
* Unfold a tracking line (MOCC) for connection matrices.
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
*Parameters: output
* NFI     number of regions (volumes and surfaces).
* NZON    index-number of the mixture type assigned to each volume.
*Parameters: input/output
* NSEG    number of segment of this track folded/unfolded.
* NOM     segment index array of this track folded/unfolded.
* H       segment length array of this track folded/unfolded.
*
*Comments:
*   Preconditioner calculation.
*   ex: -5 2 3 1 -4 -4 2 1 -1 -1 1 -5
*        r        v  v      r  r    r
*   -->    2 3 1 -4    2 1
*   where r stands for reflective boundary condition,
*         v     \\     void             \\          .
* 
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NFI,NZON(NFI),NSEG,NOM(NSEG)
      DOUBLE PRECISION H(NSEG)
*----
*  LOCAL VARIABLES
*----
      INTEGER IBCV
      PARAMETER(IBCV=-7)
      INTEGER I,IP,NOMI,NZI
*----
      I=1
      IP=0
      DO WHILE (I.LE.NSEG)
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         IF ((NZI.GE.0).OR.(NZI.EQ.IBCV)) THEN
*        this cell is a volume or a voided boundary
            IF ((IP.GE.1).AND.(NOM(IP).EQ.NOMI)) THEN
               H(IP)=H(IP)+H(I)
            ELSE
               IP=IP+1
               NOM(IP)=NOM(I)
               H(IP)=H(I)
            ENDIF
         ENDIF
         I=I+1
      ENDDO
      IF ((NOM(IP).EQ.NOM(1)).AND.(IP.GT.1)) THEN
         H(1)=H(1)+H(IP)
         IP=IP-1
      ENDIF
      NSEG=IP
*     
      RETURN
      END
