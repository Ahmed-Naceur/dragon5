*DECK PSPLEG
      SUBROUTINE PSPLEG(IPRINT,ISPSP ,ITYPE ,ICOLR ,NSUR  ,NVOL  ,
     >                  NAMLEG,NUNKNO,FLUX  ,NREGT ,
     >                  MATALB,KEYMRG,KEYFLX,COLREG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Associate a color to a region and print legend.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPRINT  print level.
* ISPSP   PSP file unit.
* ITYPE   type of graphic:
*         =  0 color per region number;
*         =  1 color per material;
*         =  2 color for flux (one group);
*         =  3 color for flux (multigroup);
*         =  4 color per material for homogenization (HMIX).
* ICOLR   color set used:
*         = -4 fill hsb with no-contour;
*         = -3 fill cmyk with no-contour;
*         = -2 fill rgb with no-contour;
*         = -1 fill bw with no-contour;
*         =  0 no fill contour only;
*         =  1 fill bw and contour;
*         =  2 fill rgb and contour;
*         =  3 fill cmyk and contour;
*         =  4 fill hsb and contour.
* NSUR    number of outer surface.
* NVOL    maximum number of regions.
* NAMLEG  legend name.
* NUNKNO  number of unknowns.
* FLUX    unknown vector.
* NREGT   dimension of KEYFLX vector.
* MATALB  albedo-material of regions.
* KEYMRG  merge index.
* KEYFLX  flux location.
* COLREG  region color.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      REAL             WLINE
      PARAMETER       (IOUT=6,WLINE=0.002,NAMSBR='PSPLEG')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IPRINT,ISPSP,ITYPE,ICOLR,NSUR,NVOL,
     >                 NUNKNO,NREGT
      INTEGER          MATALB(NSUR:NVOL),KEYMRG(NSUR:NVOL),
     >                 KEYFLX(NVOL)
      REAL             FLUX(NUNKNO),COLREG(4,NVOL)
      REAL             COLTMP(4)
      CHARACTER        NAMLEG*24
*----
*  LOCAL PARAMETERS
*----
      CHARACTER        COLNAM*4,LEGTXT*48,FLXTXT*80
      INTEGER          MXMIX,MREG,IVOL,IMX,IRG,ICOLA,
     >                 ILEG,IFRM,MXCOL,ICOLF,IKEY
      INTEGER          KMX,ICT
      REAL             XYPOS(2),POSL,POSB,DELX,DELY,DELXC,DELYC,
     >                 XYPTS(2,4),FLXMIN,FLXMAX,DELFLX,COLFLX(4)
      INTEGER          KFS,KFR,KSS,KSR
*----
*  INITIALIZE LEGEND
*----
      KFS=0
      KFR=0
      KSS=0
      KSR=0
      ICOLA=ABS(ICOLR)
      IF(ICOLA .GT. 0) THEN
        KFS=1
        KSR=1
      ENDIF
      IF(ICOLA .GE. 2) THEN
        LEGTXT='Color by '//NAMLEG
      ELSE
        LEGTXT='Graylevel by '//NAMLEG
      ENDIF
      ILEG=1
      IF(IPRINT .LE. 0) THEN
        ILEG=0
      ENDIF
*----
*  GENERATE RANDOM COLOR
*  FOR RGB USE ALL THREE COLORS
*  FOR BW USE ONLY FIRST COLOR
*  SKIP FOR NONE
*----
      IF(ICOLA .GT. 0) THEN
        POSL=0.0
        POSB=10.0
        XYPOS(1)=POSL
        XYPOS(2)=POSB
        IF(ILEG .EQ. 1) THEN
          CALL PSTEXT(ISPSP,6,'Legend',
     >      XYPOS,0.1,0,0.0)
        ENDIF
        IF(ITYPE .EQ. 0) THEN
*----
*  COMPUTE NUMBER OF REGIONS AFTER MERGE
*----
          MREG=0
          DO 100 IVOL=1,NVOL
            MREG=MAX(MREG,KEYMRG(IVOL))
 100      CONTINUE
*----
*  GENERATE ONE COLOR PER REGION
*----
          POSB=POSB-0.2
          XYPOS(2)=POSB
          IF(ILEG .EQ. 1) THEN
            CALL PSTEXT(ISPSP,48,LEGTXT,XYPOS,0.1,0,0.0)
          ENDIF
          POSB=POSB-0.2
          IF(MREG .GT. 10000) THEN
            ILEG=0
          ENDIF
          DELX=0.2
          DELY=DELX/2.0
          DELXC=DELY
          DELYC=DELXC/4.0
          DO 110 IRG=1,MREG
            IFRM=0
            IF(MOD(IRG-1,30) .EQ. 0 .AND. ILEG .EQ. 1) THEN
              POSB=POSB-DELY
            ENDIF
            DO 111 IVOL=1,NVOL
              IF(KEYMRG(IVOL) .EQ. IRG) THEN
                CALL PSPCOL(ICOLA,MREG,IRG,COLREG(1,IVOL))
                IF(IFRM .EQ. 0 .AND. ILEG .EQ.1) THEN
                  IFRM=IFRM+1
                  POSL=MOD(IRG-1,30)*DELX
                  XYPTS(1,1)=POSL
                  XYPTS(2,1)=POSB
                  XYPTS(1,2)=POSL+DELX
                  XYPTS(2,2)=POSB
                  XYPTS(1,3)=POSL+DELX
                  XYPTS(2,3)=POSB+DELY
                  XYPTS(1,4)=POSL
                  XYPTS(2,4)=POSB+DELY
                  CALL PSDREG(ISPSP,4,XYPTS)
                  IF(ICOLA .GT. 0) THEN
                    CALL PSFILL(ISPSP,ICOLA,COLREG(1,IVOL),KFS,KFR)
                  ENDIF
                  CALL PSSTRK(ISPSP,WLINE,KSS,KSR)
                  WRITE(COLNAM,'(I4)') IRG
                  XYPOS(1)=POSL+DELXC
                  XYPOS(2)=POSB+DELYC
                  CALL PSTEXT(ISPSP,4,COLNAM,XYPOS,0.05,1,0.0)
                ENDIF
              ENDIF
 111        CONTINUE
 110      CONTINUE
        ELSE IF(ITYPE .EQ. 1 .OR. ITYPE .EQ. 4) THEN
*----
*  COMPUTE NUMBER OF MIXTURES
*----
          MXMIX=0
          DO 120 IVOL=1,NVOL
            MXMIX=MAX(MXMIX,MATALB(IVOL))
 120      CONTINUE
          POSB=POSB-0.2
          XYPOS(2)=POSB
          IF(ILEG .EQ. 1) THEN
            CALL PSTEXT(ISPSP,32,LEGTXT,XYPOS,0.1,0,0.0)
          ENDIF
          POSB=POSB-0.2
          IF(MXMIX .GT. 10000) THEN
            ILEG=0
          ENDIF
          KMX=0
          DELX=0.2
          DELY=DELX/2.0
          DELXC=DELY
          DELYC=DELXC/4.0
*----
*  GENERATE ONE COLOR PER MIXTURE
*----
          DO 130 IMX=0,MXMIX
            KMX=KMX+1
            IFRM=0
            IF(MOD(KMX-1,30).EQ.0 .AND. ILEG .EQ. 1) THEN
              POSB=POSB-DELY
            ENDIF
            CALL PSPCOL(ICOLA,MXMIX,IMX,COLTMP)
            IF (ILEG.EQ.1) THEN
              POSL=MOD(KMX-1,30)*DELX
              XYPTS(1,1)=POSL
              XYPTS(2,1)=POSB
              XYPTS(1,2)=POSL+DELX
              XYPTS(2,2)=POSB
              XYPTS(1,3)=POSL+DELX
              XYPTS(2,3)=POSB+DELY
              XYPTS(1,4)=POSL
              XYPTS(2,4)=POSB+DELY
              CALL PSDREG(ISPSP,4,XYPTS)
              IF(ICOLA .GT. 0) THEN
                 CALL PSFILL(ISPSP,ICOLA,COLTMP,KFS,KFR)
              ENDIF
              CALL PSSTRK(ISPSP,WLINE,KSS,KSR)
              WRITE(COLNAM,'(I4)') IMX
              XYPOS(1)=POSL+DELXC
              XYPOS(2)=POSB+DELYC
              CALL PSTEXT(ISPSP,4,COLNAM,XYPOS,0.05,1,0.0)
            ENDIF
*----
*  ASSOCIATE MIXTURE COLOR WITH REGION
*----
            DO 131 IVOL=1,NVOL
              IF(MATALB(IVOL) .EQ. IMX) THEN
                DO 132 ICT=1,4
                  COLREG(ICT,IVOL)=COLTMP(ICT)
 132            CONTINUE
              ENDIF
 131        CONTINUE
 130      CONTINUE
        ELSE IF(ITYPE .EQ. 2 .OR. ITYPE .EQ. 3 .OR.
     >          ITYPE .EQ. 5 .OR. ITYPE .EQ. 6) THEN
*----
*  COMPUTE NUMBER OF REGIONS AFTER MERGE
*----
          POSB=POSB-0.2
          XYPOS(2)=POSB
          IF(ILEG .EQ. 1) THEN
            CALL PSTEXT(ISPSP,32,LEGTXT,XYPOS,0.1,0,0.0)
          ENDIF
          POSB=POSB-0.2
*----
*  FIND MAXIMUM AND MINIMUM FLUX
*----
          FLXMAX=FLUX(KEYFLX(1))
          FLXMIN=FLUX(KEYFLX(1))
          DO 150 IRG=2,NREGT
            IKEY=KEYFLX(IRG)
            FLXMAX=MAX(FLXMAX,FLUX(IKEY))
            FLXMIN=MIN(FLXMIN,FLUX(IKEY))
 150      CONTINUE
          MXCOL=20
          DELFLX=(FLXMAX-FLXMIN)/REAL(MXCOL)
          WRITE(FLXTXT,5000) FLXMIN,DELFLX,FLXMIN,DELFLX
          XYPOS(2)=POSB
          IF(ILEG .EQ. 1) THEN
            CALL PSTEXT(ISPSP,80,FLXTXT,XYPOS,0.1,0,0.0)
          ENDIF
          POSB=POSB-0.2
          DELX=0.2
          DELY=DELX/2.0
          DELXC=DELY
          DELYC=DELXC/4.0
*----
*  GENERATE ONE COLOR PER FLUX LEVEL
*  COLOR I IS GIVEN BY:
*  I=MIN(INT((FLUX-FLXMIN)/DELFLX)+1,MXCOL)
*----
          POSB=POSB-DELY
          DO 160 ICOLF=1,MXCOL
            CALL PSPCOL(ICOLA,MXCOL,ICOLF,COLFLX)
            POSL=MOD(ICOLF-1,30)*DELX
            XYPTS(1,1)=POSL
            XYPTS(2,1)=POSB
            XYPTS(1,2)=POSL+DELX
            XYPTS(2,2)=POSB
            XYPTS(1,3)=POSL+DELX
            XYPTS(2,3)=POSB+DELY
            XYPTS(1,4)=POSL
            XYPTS(2,4)=POSB+DELY
            CALL PSDREG(ISPSP,4,XYPTS)
            IF(ICOLA .GT. 0) THEN
              CALL PSFILL(ISPSP,ICOLA,COLFLX,KFS,KFR)
            ENDIF
            CALL PSSTRK(ISPSP,WLINE,KSS,KSR)
            WRITE(COLNAM,'(I4)') ICOLF
            XYPOS(1)=POSL+DELXC
            XYPOS(2)=POSB+DELYC
            CALL PSTEXT(ISPSP,4,COLNAM,XYPOS,0.05,1,0.0)
 160      CONTINUE
          DO 170 IRG=1,NREGT
            IKEY=KEYFLX(IRG)
            ICOLF=INT((FLUX(IKEY)-FLXMIN)/DELFLX)+1
            ICOLF=MIN(ICOLF,MXCOL)
            DO 171 IVOL=1,NVOL
              IF(KEYMRG(IVOL) .EQ. IRG) THEN
                CALL PSPCOL(ICOLA,MXCOL,ICOLF,COLREG(1,IVOL))
              ENDIF
 171        CONTINUE
 170      CONTINUE
        ENDIF
      ENDIF
      RETURN
*----
*  FORMAT
*----
 5000 FORMAT(1P,E9.2,'+(i-1)*',E9.2,
     >       ' < Flux(i) <= ',E9.2,'+i*',E9.2)
      END
