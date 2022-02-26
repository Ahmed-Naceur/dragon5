*DECK XCWROD
      SUBROUTINE XCWROD(NRIN,NRODS,NRODR,RODR,RODP,RADC,NFSEG,NLSEG,
     >                  SEGLEN,NRSEG,NNSEG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform rod tracking for 2-D cluster geometry.
*
*Copyright:
* Copyright (C) 1992 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G.Marleau
*
*Parameters: input
* NRIN    current region number.
* NRODS   integer description of rod type:
*         NRODS(1) = number of rod;
*         NRODS(2) = number of subrods in rod.
* NRODR   subrod region.
* RODR    subrod radius.
* RODP    rod position:
*         RODP(1,IRD) = X-position;
*         RODP(2,IRD) = Y-position.
* RADC    Y-position of track.
*
*Parameters: output
* NFSEG   initial segment position.
* NLSEG   final segment position.
* SEGLEN  length of track.
* NRSEG   region crossed by track.
* NNSEG   region crossed by track (left).
*
*----------------------------------------------------------------------
*
      INTEGER    NRIN,NRODS(2),NRODR,NFSEG,NLSEG,NRSEG(*),NNSEG(*)
      REAL       RODR(*),RODP(2,*)
      DOUBLE PRECISION SEGLEN(*),RADC,RADR,RADR2
*----
*  FILL IN SEGLEN FROM THE END STARTING WITH ROD FURTHER FROM
*  TRACK STARTING POINT UNTIL CENTER OF TRACK REACHED
*----
      NPROD=(NRODS(1)+3)/2
      NSBR=NRODS(2)
      IF(RADC.GE.0.0D0) THEN
        IPDEB=1
        IPFIN=NPROD
        IPSTP=1
        IMDEB=NPROD
        IMFIN=1
        IMSTP=-1
      ELSE
        RADR=RODP(2,1)-RADC
        IF(ABS(RADR).LT.RODR(NSBR)) THEN
          IPDEB=NRODS(1)+1
          IPFIN=MAX(2,NRODS(1)+1-NPROD)
        ELSE
          IPDEB=NRODS(1)
          IPFIN=MAX(1,NRODS(1)-NPROD)
        ENDIF
        IPSTP=-1
        IMDEB=IPFIN
        IMFIN=IPDEB
        IMSTP=1
      ENDIF
      NXSEG=NLSEG
      DO 100 IRZ=IPDEB,IPFIN,IPSTP
        IF(IRZ.EQ.NRODS(1)+1) THEN
          IRD=1
        ELSE
          IRD=IRZ
        ENDIF
        RADR=RODP(2,IRD)-RADC
        RADR2=RADR*RADR
        NREG=NRIN
        IF( ABS(RADR).LT.RODR(NSBR) ) THEN
*----
*  ROD INTERCEPS
*----
          XTRA=SQRT(RODR(NSBR)*RODR(NSBR)-REAL(RADR2))
          XLST=RODP(1,IRD)+XTRA
          XFST=RODP(1,IRD)-XTRA
          IF(XLST.LT.0.0) THEN
*----
*  CENTER OF TRACK REACHED/EXIT
*----
            GO TO 1000
          ELSE
*----
*  SET POINTERS TO SEGLEN VECTOR W.R.T. LAST POSITION FREE
*----
            NFLSEG=NXSEG-2*NSBR
            NLLSEG=NXSEG
            NXSEG=NFLSEG
          ENDIF
          SEGLEN(NLLSEG)=XLST
          NRSEG(NLLSEG)=NREG
          NNSEG(NFLSEG+1)=-NREG
          NLLSEG=NLLSEG-1
          NREG=NRODR
          NFLSEG=NFLSEG+1
          SEGLEN(NFLSEG)=XFST
          NRSEG(NFLSEG)=NREG
          NNSEG(NLLSEG+1)=-NREG
          DO 110 ISBR=NSBR-1,1,-1
            IF( ABS(RADR).LT.RODR(ISBR) ) THEN
*----
*  SUBROD INTERCEPS
*----
              XTRA=SQRT(RODR(ISBR)*RODR(ISBR)-REAL(RADR2))
              SEGLEN(NLLSEG)=RODP(1,IRD)+XTRA
              NRSEG(NLLSEG)=NREG
              NNSEG(NFLSEG+1)=-NREG
              NLLSEG=NLLSEG-1
              NREG=NREG-1
              NFLSEG=NFLSEG+1
              SEGLEN(NFLSEG)=RODP(1,IRD)-XTRA
              NRSEG(NFLSEG)=NREG
              NNSEG(NLLSEG+1)=-NREG
            ENDIF
 110      CONTINUE
        ENDIF
 100  CONTINUE
 1000 CONTINUE
      NLSEG=NXSEG
*----
*  FILL IN SEGLEN FROM THE BEGINNING STARTING WITH ROD CLOSEST FROM
*  TRACK STARTING POINT UNTIL CENTER OF TRACK REACHED
*----
      NXSEG=NFSEG
      DO 200 IRZ=IMDEB,IMFIN,IMSTP
        IF(IRZ.EQ.NRODS(1)+1) THEN
          IRD=1
        ELSE
          IRD=IRZ
        ENDIF
        RADR=RODP(2,IRD)-RADC
        RADR2=RADR*RADR
        NREG=NRIN
        IF( ABS(RADR).LT.RODR(NSBR) ) THEN
*----
*  ROD INTERCEPS
*----
          XTRA=SQRT(RODR(NSBR)*RODR(NSBR)-REAL(RADR2))
          XLST=RODP(1,IRD)+XTRA
          XFST=RODP(1,IRD)-XTRA
          IF(XLST.LT.0.0) THEN
*----
*  SET POINTERS TO SEGLEN VECTOR W.R.T. FIRST POSITION FREE
*----
            NLLSEG=NXSEG+2*NSBR
            NFLSEG=NXSEG
            NXSEG=NLLSEG
          ELSE
*----
*  CENTER OF TRACK REACHED/EXIT
*----
            GO TO 2000
          ENDIF
          SEGLEN(NLLSEG)=XLST
          NRSEG(NLLSEG)=NREG
          NNSEG(NFLSEG+1)=-NREG
          NLLSEG=NLLSEG-1
          NREG=NRODR
          NFLSEG=NFLSEG+1
          SEGLEN(NFLSEG)=XFST
          NRSEG(NFLSEG)=NREG
          NNSEG(NLLSEG+1)=-NREG
          DO 210 ISBR=NSBR-1,1,-1
            IF( ABS(RADR).LT.RODR(ISBR) ) THEN
*----
*  SUBROD INTERCEPS
*----
              XTRA=SQRT(RODR(ISBR)*RODR(ISBR)-REAL(RADR2))
              SEGLEN(NLLSEG)=RODP(1,IRD)+XTRA
              NRSEG(NLLSEG)=NREG
              NNSEG(NFLSEG+1)=-NREG
              NLLSEG=NLLSEG-1
              NREG=NREG-1
              NFLSEG=NFLSEG+1
              SEGLEN(NFLSEG)=RODP(1,IRD)-XTRA
              NRSEG(NFLSEG)=NREG
              NNSEG(NLLSEG+1)=-NREG
            ENDIF
 210      CONTINUE
        ENDIF
 200  CONTINUE
 2000 CONTINUE
      NFSEG=NXSEG
      RETURN
      END
