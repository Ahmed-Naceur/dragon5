*DECK TRKHEX
      SUBROUTINE TRKHEX(IPRT,NCEL,FVOL,REMESH,MESH,PAS1,A,COS1,COS2,
     +        COS3,POP,STAIRS,IPLANZ,FACST,NDIM,NCYL,MAT,IFILE,
     +        IANGL,POIDS,MAT2,NSECT,T0,TSEC,V0,VSEC,PAS2,RAYON,
     +        ZMIN,ZMAX,FACB,NVOL,SECTOR,NSMIN,SURB,RAUX,NSURF,CORN,
     +        ICOR,VOISIN,NCEL2,NSOUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the tracking information related to a
* hexagonal heterogeneous assembly for a given
* angle in a 3D or 2D geometry.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* M. Ouisloumen
*
*Parameters: input
* IPRT    print option.                           
* NCEL    number of cells of the assembly.
* FVOL    first zone in the cell. 
* NVOL    first volume number in the cell. 
* NCYL    number of cylinder in the cell.
* REMESH  contains the coordinates of the hexagon centers (height for
*         Z), the location of the cylinders and their radius.
*         The order is:
*         XH(1),XH(2),...,YH(1),YH(2),...,ZH(1),ZH(2),...,XC(1),YC(1),
*         ZC(1),RC(1),XC(2),YC(2),ZC(2),RC(2)...
* A       lenght of one of the hexagone.
* COS1    X director cosine.
* COS2    Y director cosine.
* COS3    Z director cosine.
* PAS1    line spacing in Y.
* PAS2    line spacing in Z.
* STAIRS  cell maximum number in plane.
* FACST   surface maximum number in plane.
* NDIM    number of dimensions of problem.
* IFILE   file track unit number.
* NSECT   number of sectors in volume. 
* FACB    first bottom or top surface number.
* SECTOR  flag for sector.
* SURB    initial surface numbering. 
* VOISIN  6 neighbors of each cell in the X-Y plane.
* NSOUT   Total number of surfaces in global geometry.
*
*Parameters: input/output
* IANGL   angle number.
*
*Parameters: scratch
* POP     undefined.
* MAT     undefined.
* T0      undefined.
* TSEC    undefined.
* V0      undefined.
* VSEC    undefined.
* RAUX    undefined.
* MESH    undefined.
* IPLANZ  undefined.
* POIDS   undefined.
* MAT2    undefined.
* RAYON   undefined.
* ZMIN    undefined.
* ZMAX    undefined.
* NSMIN   undefined.
* NSURF   undefined.
* CORN    undefined.
* ICOR    undefined.
* NCEL2   undefined.
*
*Comments:
*   Assembly coordinate system
*                                            .Z
*                                            .
*                                            .    . Y
*                    ++++++++                .   .
*                 +           +              .  .
*                +              +            . .
*                  +          +              ..
*                I   ++++++++   I          O ..............   X
*                I   I      I   I
*                I   I      I   I           
*                I   I      I   I           
*                  + I      I +
*                    ++++++++
*
*---------------------------  TRKHEX  ---------------------------------
*
      IMPLICIT NONE
      INTEGER          IOUT
      DOUBLE PRECISION SQRT3,PI
      CHARACTER*6      NAMSBR
      PARAMETER       (IOUT=6,SQRT3=1.732050807568877D0,
     >                 PI=3.141592653589793D0,NAMSBR='TRKHEX')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IPRT,NCEL,MESH,IPLANZ,NDIM,ICOR,IFILE,NSMIN,
     >                 NSURF,NCEL2,IANGL,NSOUT
      REAL             REMESH(MESH)
      DOUBLE PRECISION POP(*),A,COS1,COS2,COS3,PAS1,PAS2,RAYON,POIDS,
     >                 ZMIN,ZMAX
      DOUBLE PRECISION T0(*),TSEC(*),RAUX(*)
      INTEGER          FVOL(NCEL),STAIRS(IPLANZ),FACST(IPLANZ),
     >                 NCYL(NCEL),MAT(*),MAT2(*),NSECT(*),V0(*),
     >                 VSEC(*),FACB(*),NVOL(*),SURB(NSMIN:*),
     >                 CORN(ICOR),VOISIN(6,NCEL2)
      LOGICAL          SECTOR
*----
*  LOCAL VECTORS
*----
      INTEGER          SURC(6),SURF1(8),SURF2(8),FACEM(5),
     >                 FACES(5),FACEM1(5),FACEM2(5),FACEMD(5),
     >                 SURFX(8),FACEMU(5),MAT3(10)
      DOUBLE PRECISION T(9),XDR(8),YDR(8),ZDR(8)
      LOGICAL          START,DIRUP,LGSTOR,LGDIM,LGPG1,LGFAC,
     >                 LGMAT3,LGDIR,LGPER,LGFAC8,LGFAC7,
     >                 LTHROU,LGDEB,LGPASS,PASSU,PASSD,LGPAS0,
     >                 LGPASU,LGPASD,LGOUT1,LGOUTU,LGOUTD,LGG,LGG1
*----
*  ADDITIONNAL VARIABLES
*----
      INTEGER          IFONC,IFCOUR,JSEC2,JSEC,JCC,J,IVS1,IVS0,JC,
     >                 JCC0,KMAX,N,IZONE,ISECTR,L,LMIN,ILINE,
     >                 IZZ,IZ0,LMAX,L0,LL,LSTEP,ISV,JSV,
     >                 JVOUT2,JVOUT1,IFV,IVSEC1,IS,IVS,IVSEC2,JAUX,
     >                 ISEC,IS1,IAUX,LPOP,IMAM,IFOUT2,KSECT,JMAX,JMIN,
     >                 IXY,KFACE,ISURF2,ISURF1,NFACX,ISTAIR,IMAX,IW,
     >                 IFACST,LFOUT2,ISXX,KS,IP,KOUT,KFXX,MSECT,LVOUT,
     >                 IFOUT1,NSURF2,NSURF1,IVOUT1,IMAT2,IVOL,KPOP,
     >                 IFDOWN,MFACD,IMATD,ICDOWN,I3,IVAUX,ISC,ISS,
     >                 ISSX,I2,I1,JS,MFACU,IMAT0,KPC,IVOUT2,NFAC,ICXX,
     >                 IYY,NFXX,ICELLD,IXX,KPLD,IFACEU,IVOUT,IFACE,
     >                 ICELLX,IMAT,LSV,I,IX,NCOUR,LIM2,ICELL,LIM1,
     >                 ITAB,ICELL0,KFACE2,KFACED,KFACE1,KPER,
     >                 IPOP,NCELZ,LFACE,KAM,KMA,KCEL1,IPP,ICELL2,
     >                 IVDOWN,IPOS,IVUP,KUP,KDOWN,KMAT2,LVOIS,KC,
     >                 IFF1,IXF,IFF0,IAC,MM,KPER3,IFG,IC,LAUX,KVOL2,
     >                 KVOL1,IV01,IV02,K,ICYL,ISTEP,IFACE2,IFACE1,
     >                 ICC,KPL,MFAC2,MFAC,MFAC1,IFF,IPAR,ICEL0,IC8,
     >                 IPPP,IPX,KCEL2,JFF,ISURM,IVOLC,IETAG,KCELD,ISURC
      DOUBLE PRECISION ANG0,DIST,YAUX,XAUX,TTHYP0,TTHYP1,THYP,XT,YT,
     >                 YRF,YRH,PENTE,XPT,ATAN0,ANG2,ANG1,YPT,DEL0,
     >                 YPH2,TGA,TT1,TT0,TTHYP,YPF,DEL1,ATAN1,YPH,
     >                 COSR,SINR,R2AZ,ZZ,COAZ,
     >                 CY0MIN,CY0MAX,DIV2,COS1I,DIV1,Y0MIN,Y0MAX,SAZ,
     >                 TERM3,COEF1,TERM2,ZLA,TERM1,Z0,Z00,Z0MAX,COEF2,
     >                 Z0MIN,COS2I,XZTMIN,ZTMAX,XZTMAX,SSS,Y0COS,ZTMIN,
     >                 Y00,SSQ,CZMAX,ZZ0,CZMIN,CZX1,CZX2,XY02,Y0,XY01,
     >                 DELTA,SDELTA,ALP,Z2,PASY,Z1,Y,X,ACOS6,PASZ,EPS1,
     >                 EPS2,EPS3,EPS4,EPS5,EPS6,EPSY,EPSZ,EPSX,EPSS,
     >                 YDEN2,R2,YDEN1,YDDD,YDEN,RSCOS2,POID1,AY,SCOS3,
     >                 AX,TGDIR,SCOS2,AZ,SCOS1,TTT,ZP,ZM,T0MIN,T0MAX,
     >                 TM,R,YYC,YCYL,XCYL,SDEL,TP,DEL,BZ,DEE,YDROIT,
     >                 XDROIT,XMA2,XMA,XPA2,Y6,Y1,Y3,ZDROIT,Y4,XPA,
     >                 YY,Y2,Y5,ZZZ,YUP,ZDOWN,ZUP,FACY,FACX,YDOWN,
     >                 FTX,FTY
      DOUBLE PRECISION AUX1,WEIGHT
      DOUBLE PRECISION YT1,YT2,YT3,YT4,YT5,YT6,XT1,XT2,XT3,XT4,EPST,DZ0
      CHARACTER        HSMG*131
*
*         FONCTION NECESSAIRE POUR LA ROTATION ET LA RECHERCHE
*         DES SURFACES EXTERNES
*
      IFONC(N,L)= 2+(N-1)*(L+3*(N-2))
      IFCOUR(N)=NINT( (4.+SQRT(1.+4.*FLOAT(N-1)/3.)
     +                 +SQRT(1.+4.*FLOAT(N-2)/3.))*.25)
*
*         TRACKING PAR RAPPORT AU PLAN YOZ. LES DIMENSIONS DU PLAN
*         A TRACKER SONT DETERMINEES EN FONCTION DES ANGLES
*
      DO 767 ITAB=1,5
        SURF1(ITAB)=0
        SURF2(ITAB)=0
        FACEM(ITAB)=0
 767  CONTINUE
      T(9)=1.0D30
      EPS1=1.0D-6
      EPST=EPS1*EPS1
      DZ0=0.0D0
      EPS2=2.0D0*EPS1 
      EPS3=5.0D0*EPS1 
      EPS4=1.0D-2*EPS1 
      EPS5=9.0D0*EPS1 
      EPS6=10.0D0*EPS1 
      EPSS=MAX(EPS1*ABS(COS1),EPS1*ABS(COS2))
      EPSS=MAX(EPSS,EPS1*ABS(COS3))
      EPSX=EPSS
      EPSY=EPSS
      EPSZ=EPSS
      PASZ=PAS2
      PASY=PAS1
      Z1=0.0
      Z2=0.0
      NCOUR=1
      IF(STAIRS(1).GT.1)NCOUR=IFCOUR(STAIRS(1))
      LIM2=STAIRS(1)
      LIM1=IFONC(NCOUR,0)
      ACOS6=.5*SQRT3*A
      X=0.
      Y=0.
      ICELL=1
      ICELL0=0
      KPER=0
      IPOP=0
      KFACE1=9
      KFACE2=9
      KFACED=9
      LGOUT1=.FALSE.
      LGOUTD=.FALSE.
      LGOUTU=.FALSE.
      LGDEB=.TRUE.
      LGPAS0=.FALSE.
      LGMAT3=.FALSE.
      START=.TRUE.
      DIRUP=.TRUE.
      NCELZ=2*NCEL
      LGDIM=NDIM.EQ.3
      LFACE=6
      SCOS1=COS1*COS1
      SCOS2=COS2*COS2
      AZ=SCOS1+SCOS2
      LGPG1=.FALSE.
      TGDIR=COS2/COS1
      YDDD=SQRT3*TGDIR
      YDEN=1.+YDDD
      YDEN1=-1.+YDDD
      YDEN2=1.-YDDD
      ICDOWN=0
      ICELLD=0
      ICELLX=0
      IETAG=0
      IFDOWN=0
      IMAT2=0
      IMATD=0
      ISURC=0
      ISURM=0
      IVOLC=0
      KCEL1=0
      KCEL2=0
      KCELD=0
      KMAX=0
      KPER3=0
      KPL=0
      KPOP=0
      KSECT=0
      LPOP=0
      MFACD=0
      NSURF1=0
      IF(LGDIM) THEN
        R2=RAYON*RAYON
        SCOS3=COS3*COS3
        AX=SCOS2+SCOS3
        AY=SCOS1+SCOS3
        RSCOS2=R2*SCOS2
        LGPG1=IPLANZ.GT.1
        LFACE=8
        POID1=POIDS
*
*--L'ASSEMBLAGE A TRAITER EST CONFINE DANS UN CYLINDRE DE RAYON=RAYON
*  FERME PAR DEUX CALOTTES SPHERIQUES A CES EXTREMITES
*
        ALP=SCOS1+SCOS3
        ZLA=SCOS1/R2
*
*--EXTEMUMS DE Z0: LIMITES SELON L'AXE Z DU PLAN A TRACKER
*
        TERM1=SCOS1+SCOS2*SCOS3
        TERM2=-SCOS1+SCOS2*SCOS3
        TERM3=RAYON/(COS1*SQRT((1.-SCOS3)*TERM1))
        COEF1=TERM3*TERM1
        COEF2=TERM3*TERM2
        Z0MIN=ZMIN+COEF1
        Z0MIN=MIN(Z0MIN,ZMIN-COEF1)
        Z0MIN=MIN(Z0MIN,ZMIN+COEF2)
        Z0MIN=MIN(Z0MIN,ZMIN-COEF2)
        Z0MAX=ZMAX+COEF1
        Z0MAX=MAX(Z0MAX,ZMAX-COEF1)
        Z0MAX=MAX(Z0MAX,ZMAX+COEF2)
        Z0MAX=MAX(Z0MAX,ZMAX-COEF2)
        PASZ=(Z0MAX-Z0MIN)*PAS2
        Z0=Z0MIN+PASZ*.5
        Z00=Z0
        SAZ=SQRT(AZ)
        CY0MIN=-RAYON*SAZ/ABS(COS1)
        CY0MAX=-CY0MIN
        COAZ=COS3/AZ
        R2AZ=R2*AZ
        ZZ=RAYON*ABS(COS3*COS2/(COS1*SAZ))
        Y0MIN=999999.0
        Y0MAX=-Y0MIN
      ELSE
        R2=0.0D0
        POID1=0.0D0
        RSCOS2=0.0D0
        ALP=0.0D0
        ZLA=0.0D0
        CY0MIN=0.0D0
        CY0MAX=0.0D0
        COAZ=0.0D0
        R2AZ=0.0D0
        ZZ=0.0D0
        Z0=0.0D0
        Z00=0.0D0
        Y0MIN=-RAYON/ABS(COS1)
        Y0MAX=-Y0MIN
        PASY=(Y0MAX-Y0MIN)*PAS1
        POIDS=POIDS*PASY
        Z0MIN=999999.0
        Z0MAX=-Z0MIN
      ENDIF
      DIV1=1./(COS2+COS1*SQRT3)
      DIV2=1./(COS2-COS1*SQRT3)
      COS1I=1./COS1
      COS2I=1./COS2
  777 CONTINUE
      ILINE=0
      IF(LGDIM) THEN
*
*--DOMAINE DE VARIATION DE Y0 POUR Z0 FIXE
*
        CZX1=Z0-ZZ
        CZX2=Z0+ZZ
        CZMIN=MIN(CZX1,CZX2)
        CZMAX=MAX(CZX1,CZX2)
        IF(CZMAX.LT.ZMIN)THEN
          ZZ0=Z0-ZMIN
          DELTA=-ZLA*ZZ0*ZZ0+ALP
          IF(DELTA.GE.0.) THEN
            SDELTA=SQRT(DELTA)
            XY01=ZZ0*COS2*COS3
            XY02=RAYON*SDELTA
            Y0MIN=(XY01-XY02)/ALP
            Y0MAX=(XY01+XY02)/ALP
          ELSE
            CALL XABORT('TRKHEX: ALGORITHME FAILURE -A')
          ENDIF
        ELSEIF(CZMIN.LT.ZMAX)THEN
          Y0MIN=CY0MIN
          Y0MAX=CY0MAX
        ELSE
          ZZ0=Z0-ZMAX
          DELTA=-ZLA*ZZ0*ZZ0+ALP
          IF(DELTA.GE.0.) THEN
            SDELTA=SQRT(DELTA)
            XY01=ZZ0*COS2*COS3
            XY02=RAYON*SDELTA
            Y0MIN=(XY01-XY02)/ALP
            Y0MAX=(XY01+XY02)/ALP
          ELSE
            CALL XABORT('TRKHEX: ALGORITHME FAILURE - B')
          ENDIF
        ENDIF
        PASY=(Y0MAX-Y0MIN)*PAS1
        POIDS=POID1*PASZ*PASY
        Y0=Y0MIN
        ILINE=ILINE+1
  611   SSQ=SCOS1*(R2-Y0*Y0)+RSCOS2
        IF(SSQ.GT.0.) THEN
        SSS=SQRT(SSQ)
        Y0COS=Y0*COS2
        XZTMAX=Z0+(SSS-Y0COS)*COAZ
        XZTMIN=Z0+(-SSS-Y0COS)*COAZ
        ZTMAX=MAX(XZTMAX,XZTMIN)
        ZTMIN=MIN(XZTMAX,XZTMIN)
        IF(ZTMAX.LT.ZMIN) THEN
          Y0=Y0+PASY
          IF(Y0.GT.Y0MAX) GOTO 53
          GOTO 611
        ELSEIF(ZTMIN.GT.ZMAX) THEN
          Y0=Y0+PASY
          IF(Y0.GT.Y0MAX) GOTO 53
          GOTO 611
        ENDIF
        ENDIF
      ELSE
        Y0=Y0MIN+PASY*.5
        ILINE=ILINE+1
      ENDIF
      Y00=Y0
      LGDEB=.TRUE.
      IX=1
      I=CORN(1)
      LGPER=.FALSE.
      GOTO 42
*
*         INTERSECTION DES PLANS DE L'HEXAGONE AVEC LA DROITE
*
   10 YY=Y-Y0
      TERM1=SQRT3*(X+A)
      TERM2=SQRT3*(X-A)
      T(1)=(TERM1+YY)*DIV1
      T(2)=(YY+ACOS6)*COS2I
      T(3)=(YY-TERM2)*DIV2
      T(4)=(YY+TERM2)*DIV1
      T(5)=(YY-ACOS6)*COS2I
      T(6)=(YY-TERM1)*DIV2
      IF(LGDIM) THEN
        Z2=REMESH(NCELZ+ICELL) 
        IF(LGPG1) THEN
        IF(ICELL.GT.STAIRS(1))THEN
          KPL=INT(AINT(REAL(ICELL)/REAL(STAIRS(1))))
          IF(MOD(ICELL,STAIRS(1)).NE.0)KPL=KPL+1
          ICC=ICELL-STAIRS(KPL-1)
          IF(KPL.GT.2)ICC=ICC+STAIRS(KPL-2)
          Z1=REMESH(NCELZ+ICC)
          GOTO 14
        ENDIF
        ENDIF
   14   T(7)=(Z1-Z0)/COS3
        T(8)=(Z2-Z0)/COS3
        ZDR(1)=COS3*T(1)+Z0
        ZDR(2)=COS3*T(2)+Z0
        ZDR(3)=COS3*T(3)+Z0
        ZDR(4)=COS3*T(4)+Z0
        ZDR(5)=COS3*T(5)+Z0
        ZDR(6)=COS3*T(6)+Z0
        ZDR(7)=Z1
        XDR(7)=COS1*T(7)
        YDR(7)=COS2*T(7)+Y0
        XDR(8)=COS1*T(8)
        YDR(8)=COS2*T(8)+Y0
        ZDR(8)=Z2
      ENDIF
      XDR(1)=COS1*T(1)
      YDR(1)=SQRT3*(X-XDR(1)+A)+Y
      XDR(2)=COS1*T(2)
      YDR(2)=Y+ACOS6
      XDR(3)=COS1*T(3)
      YDR(3)=SQRT3*(-X+A+XDR(3))+Y
      XDR(4)=COS1*T(4)
      YDR(4)=SQRT3*(X-XDR(4)-A)+Y
      XDR(5)=COS1*T(5)
      YDR(5)=Y-ACOS6
      XDR(6)=COS1*T(6)
      YDR(6)=SQRT3*(-X+XDR(6)-A)+Y
*
*         RECHERCHE DES FACES DE SORTIE
*
   15 IFACE2=9
      IFACE1=9
      MFAC=0
      MFAC1=0
      MFAC2=0
      LGFAC=.FALSE.
      LGPASU=.FALSE.
      LGPASD=.FALSE.
      Y2=Y+ACOS6
      Y5=Y-ACOS6
      XPA=X+A
      XMA=X-A
      XPA2=X+.5*A
      XMA2=X-.5*A 
      DO 20 I=1,LFACE
         YDROIT=YDR(I)
         XDROIT=XDR(I)
         IF(LGDIM) THEN
           ZDROIT=ZDR(I)
           IF((ZDROIT.LE.Z2.OR.ABS(ZDROIT-Z2).LE.EPSZ).AND.
     +        (ZDROIT.GE.Z1.OR.ABS(ZDROIT-Z1).LE.EPSZ)) GOTO 17
           GOTO 19
         ENDIF
   17      IF((YDROIT.GE.Y5.OR.ABS(YDROIT-Y5).LE.EPSY).AND.
     +        (YDROIT.LE.Y2.OR.ABS(YDROIT-Y2).LE.EPSY))THEN
            IF((XDROIT.GE.XMA.OR.ABS(XDROIT-XMA).LE.EPSX).AND.
     +         (XDROIT.LE.XPA.OR.ABS(XDROIT-XPA).LE.EPSX))THEN
               Y4=SQRT3*(X-XDROIT-A)+Y
               Y3=SQRT3*(-X+A+XDROIT)+Y
               Y6=SQRT3*(-X+XDROIT-A)+Y
               Y1=SQRT3*(X-XDROIT+A)+Y
               IF(ABS(YDROIT-Y2).LT.EPSY) THEN
        IF(((XDROIT.GE.XMA2).OR.ABS(XDROIT-XMA2).LE.EPSX).AND.
     +  ((XDROIT.LE.XPA2).OR.ABS(XDROIT-XPA2).LE.EPSX))GOTO 18
                 GOTO 19
               ELSEIF(ABS(YDROIT-Y5).LT.EPSY) THEN
        IF(((XDROIT.GE.XMA2).OR.ABS(XDROIT-XMA2).LE.EPSX).AND.
     +  ((XDROIT.LE.XPA2).OR.ABS(XDROIT-XPA2).LE.EPSX))GOTO 18
                 GOTO 19
               ELSEIF(ABS(YDROIT-Y1).LT.EPSY) THEN
                 GOTO 18
               ELSEIF(ABS(YDROIT-Y3).LT.EPSY) THEN
                 GOTO 18
               ELSEIF(ABS(YDROIT-Y4).LT.EPSY) THEN
                 GOTO 18
               ELSEIF(ABS(YDROIT-Y6).LT.EPSY) THEN
                 GOTO 18
               ELSEIF(I.EQ.7) THEN
                 IF(XDROIT.GE.XMA2.AND.XDROIT.LE.XPA2)THEN
                   GOTO 18
                 ELSEIF(XDROIT.GE.XMA.AND.XDROIT.LE.XMA2)THEN
             IF(((YDROIT.GE.Y4).OR.ABS(YDROIT-Y4).LE.EPSY).AND.
     +        ((YDROIT.LE.Y3).OR.ABS(YDROIT-Y3).LE.EPSY))GOTO 18
                   GOTO 19
                 ELSEIF(XDROIT.GE.XPA2.AND.XDROIT.LE.XPA)THEN
             IF(((YDROIT.GE.Y6).OR.ABS(YDROIT-Y6).LE.EPSY).AND.
     +        ((YDROIT.LE.Y1).OR.ABS(YDROIT-Y1).LE.EPSY))GOTO 18
                   GOTO 19
                 ENDIF
               ELSEIF(I.EQ.8) THEN
                 IF(XDROIT.GE.XMA2.AND.XDROIT.LE.XPA2)THEN
                   GOTO 18
                 ELSEIF(XDROIT.GE.XMA.AND.XDROIT.LE.XMA2)THEN
             IF(((YDROIT.GE.Y4).OR.ABS(YDROIT-Y4).LE.EPSY).AND.
     +        ((YDROIT.LE.Y3).OR.ABS(YDROIT-Y3).LE.EPSY))GOTO 18
                   GOTO 19
                 ELSEIF(XDROIT.GE.XPA2.AND.XDROIT.LE.XPA)THEN
             IF(((YDROIT.GE.Y6).OR.ABS(YDROIT-Y6).LE.EPSY).AND.
     +        ((YDROIT.LE.Y1).OR.ABS(YDROIT-Y1).LE.EPSY))GOTO 18
                   GOTO 19
                 ENDIF
               ELSE
                 CALL XABORT('TRKHEX: ALGORITHME FAILURE - C')
               ENDIF
               GOTO 19
   18          IF(IFACE1.NE.9) THEN
                 MFAC=MFAC+1
                 FACEM(MFAC)=IFACE1
               ENDIF
               IFACE1=IFACE2
               IFACE2=I
             ENDIF
           ENDIF
   19      CONTINUE
   20 CONTINUE
*
* CAS OU LA FACE N'A PAS ETE COMPTE
*
      DO 720 IFF=1,LFACE
        IF(IFF.NE.IFACE1.AND.IFACE2.NE.IFF) THEN
          FTX=ABS(T(IFACE1)-T(IFF))
          FTY=ABS(T(IFACE2)-T(IFF))
          IF(FTX.LE.EPS2) THEN
            DO 715 JFF=1,MFAC
              IF(FACEM(JFF).EQ.IFF) GOTO 718
 715        CONTINUE
            IF(IFACE2.EQ.9)THEN
             IFACE2=IFF
            ELSE
             MFAC=MFAC+1
             FACEM(MFAC)=IFF
            ENDIF
          ELSEIF(FTY.LE.EPS2) THEN
            DO 717 JFF=1,MFAC
              IF(FACEM(JFF).EQ.IFF) GOTO 718
 717        CONTINUE
            IF(IFACE1.EQ.9)THEN
             IFACE1=IFF
            ELSE
             MFAC=MFAC+1
             FACEM(MFAC)=IFF
            ENDIF
          ENDIF
        ENDIF
 718  CONTINUE
 720  CONTINUE
      IF(IFACE1.EQ.8) THEN
        IPX=IFACE1
        IFACE1=IFACE2
        IFACE2=IPX
      ELSEIF(IFACE2.EQ.7) THEN
        IPX=IFACE1
        IFACE1=IFACE2
        IFACE2=IPX
      ENDIF
      LGFAC=MFAC.GT.0
      IF(IFACE1.GT.8.OR.IFACE2.GT.8) THEN
*
*        CAS OU LA DROITE EST TANGENTE A UNE CELLULE PERIPHERIQUE
*
        IF(LGDEB) THEN
          I=ICELL+1
          IF(I.GT.NCEL) GOTO 58
          GOTO 42
        ENDIF
        KPER=KPER+1
        MAT2(KPER)=ICELL
        IF(DIRUP) THEN
          IF(LGOUTU) THEN
            ICELL=KCEL2
            IFACE2=KFACE2
            GOTO 739
          ENDIF
        ELSE
          IF(LGOUTD) THEN
            ICELL=KCELD
            IFACE1=KFACED
            GOTO 739
          ENDIF
        ENDIF
        IF(LGPAS0) GOTO 56
        ISURC=1
        IETAG=0
        ISURM=6
        SURC(1)=1
        SURC(2)=2
        SURC(3)=3
        SURC(4)=4
        SURC(5)=5
        SURC(6)=6
        IF(LGDIM) THEN
          IF(ICELL.LE.STAIRS(1)) THEN
            IVOLC=ICELL
            I=VOISIN(1,IVOLC)
          ELSE
            IVOLC=ICELL-STAIRS(KPL-1)
            ISURC=0
            I=IVOLC
            IF(DIRUP) THEN
              IF(KPL.GT.2) IETAG=STAIRS(KPL-2)
            ELSE
              IF(KPL.EQ.IPLANZ) THEN
              IETAG=STAIRS(KPL-1)
              ELSE
              IETAG=STAIRS(KPL)
              ENDIF
            ENDIF
          ENDIF
        ELSE
          IVOLC=ICELL
          I=VOISIN(1,IVOLC)
        ENDIF
        LGPAS0=.TRUE.
        IF(I.GT.NCEL2) THEN
          I=NCEL+10
          GOTO 56
        ENDIF
        I=I+IETAG
        LGPER=.TRUE.
        GOTO 42
      ENDIF
*
*--RECHERCHE DU PLUS GRAND PARCOURS
*
      IF(LGFAC) THEN
      IF(T(IFACE1).GT.T(IFACE2))THEN
          IPPP=IFACE2
          IFACE2=IFACE1
          IFACE1=IPPP
      ENDIF
        DO 730 IPAR=1,MFAC
          IF(T(FACEM(IPAR)).GT.T(IFACE2))THEN
             IPPP=IFACE2
             IFACE2=FACEM(IPAR)
             FACEM(IPAR)=IPPP
          ELSEIF(T(FACEM(IPAR)).LT.T(IFACE1))THEN
             IPPP=IFACE1
             IFACE1=FACEM(IPAR)
             FACEM(IPAR)=IPPP
          ENDIF
 730    CONTINUE
      ENDIF
*
      IF(LGDIM) THEN
        IF(ABS(COS3).GT.EPS1) THEN
        ZUP=COS3*T(IFACE2)+Z0
        ZDOWN=COS3*T(IFACE1)+Z0
*
*    CAS OU IFACE1 ET IFACE2  FORMENT LE MEME COIN
*
        IF(ZUP.LT.ZDOWN) THEN
          IFF=IFACE1
          IFACE1=IFACE2
          IFACE2=IFF
          ZZZ=ZUP
          ZUP=ZDOWN
          ZDOWN=ZUP
        ENDIF
        ELSE
        YUP=COS2*T(IFACE2)+Y0
        YDOWN=COS2*T(IFACE1)+Y0
        IF(YUP.LT.YDOWN) THEN
          IFF=IFACE1
          IFACE1=IFACE2
          IFACE2=IFF
        ENDIF
        ENDIF
*
*  ON PREND LA FACE 7 OU 8 DANS LE CAS D'UN COIN
*
        IF(LGFAC) THEN
        DO 21 IC8=1,MFAC
           IF(FACEM(IC8).EQ.8) THEN
             FACEM(IC8)=IFACE2
             IFACE2=8
           ELSEIF(FACEM(IC8).EQ.7) THEN
             FACEM(IC8)=IFACE1
             IFACE1=7
           ENDIF
   21   CONTINUE
        ENDIF
        IF(IFACE2.EQ.7)THEN
          IFF=IFACE2
          IFACE2=IFACE1
          IFACE1=IFF
        ENDIF
        IF(IFACE1.EQ.8)THEN
          IFF=IFACE2
          IFACE2=IFACE1
          IFACE1=IFF
        ENDIF
      ELSE
        YUP=COS2*T(IFACE2)+Y0
        YDOWN=COS2*T(IFACE1)+Y0
        IF(YUP.LT.YDOWN) THEN
          IFF=IFACE1
          IFACE1=IFACE2
          IFACE2=IFF
        ENDIF
      ENDIF
      MFAC1=0
      MFAC2=0
*
*  STOCKAGE POUR LE TRAITEMENT DES FACES DU COIN
*
      IF(LGFAC)THEN
        DO 301 IFF=1,MFAC
          FACX=ABS(T(FACEM(IFF))-T(IFACE2))
          FACY=ABS(T(FACEM(IFF))-T(IFACE1))
          IF(FACY.LE.EPS1) THEN
            MFAC1=MFAC1+1
            FACEM1(MFAC1)=FACEM(IFF)
            IF(FACX.LE.EPS1) THEN
              MFAC2=MFAC2+1
              FACEM2(MFAC2)=FACEM(IFF)
            ENDIF
          ELSEIF(FACX.LE.EPS1) THEN
            MFAC2=MFAC2+1
            FACEM2(MFAC2)=FACEM(IFF)
          ELSE
            IF(FACX.LT.FACY) THEN
              MFAC2=MFAC2+1
              FACEM2(MFAC2)=FACEM(IFF)
            ELSE
              MFAC1=MFAC1+1
              FACEM1(MFAC1)=FACEM(IFF)
            ENDIF
          ENDIF
  301   CONTINUE
      ENDIF
*
*   DETERMINATION DES FACES D'ENTREE OU DE SORTIE
*
      IF(.NOT.LGDEB) THEN
        ICEL0=ICELL
        ISTEP=0
        IF(ICELL.GT.STAIRS(1)) THEN
          ISTEP=STAIRS(KPL-1)
          ICEL0=ICELL-ISTEP
        ENDIF
        IXF=0
        LGG=.TRUE.
        LGG1=.FALSE.
        IF(DIRUP) THEN
          IFF0=IFACE1
          IFF1=IFACE2
          IFF=IFACE1
  554     CONTINUE
          IF(IFF.LT.7) THEN
            LVOIS=VOISIN(IFF,ICEL0)+ISTEP
          ELSE
            LVOIS=ICELL-STAIRS(1)
          ENDIF
          IF(LVOIS.EQ.ICELL0) THEN
            IFACE1=IFF
            IF(LGG1)FACEM1(IXF)=IFF0
            GOTO 556
          ENDIF
          IXF=IXF+1
          LGG1=.FALSE.
          IF(IXF.LE.MFAC1) THEN
            IFF=FACEM1(IXF)
            LGG1=.TRUE.
            GOTO 554
          ENDIF
          IF(LGG) THEN
            LGG=.FALSE.
            IF(ABS(T(IFF0)-T(IFACE2)).LE.EPS6) THEN
             IFF=IFACE2
             IFACE2=IFACE1
             GOTO 554
            ENDIF
          ENDIF
          IFACE1=IFF0
          IFACE2=IFF1
        ELSE
          IFF0=IFACE2
          IFF1=IFACE1
          IFF=IFACE2
  557     CONTINUE
          IF(IFF.LT.7) THEN
            LVOIS=VOISIN(IFF,ICEL0)+ISTEP
          ELSE
            LVOIS=ICELL+STAIRS(1)
          ENDIF
          IF(LVOIS.EQ.ICELL0) THEN
            IFACE2=IFF
            IF(LGG1)FACEM2(IXF)=IFF0
            GOTO 556
          ENDIF
          IXF=IXF+1
          LGG1=.FALSE.
          IF(IXF.LE.MFAC2) THEN
            LGG1=.TRUE.
            IFF=FACEM2(IXF)
            GOTO 557
          ENDIF
          IF(LGG) THEN
            LGG=.FALSE.
            IF(ABS(T(IFF0)-T(IFACE1)).LE.EPS6) THEN
            IFF=IFACE1
            IFACE1=IFACE2
            GOTO 557
            ENDIF
          ENDIF
          IFACE2=IFF0
          IFACE1=IFF1
        ENDIF
      ENDIF
 556  CONTINUE
      ICELL0=ICELL
*
*  ELEMINATION DES PARCOURS NULS SAUF LE DERNIER.
*
      IF(ABS(T(IFACE1)-T(IFACE2)).GT.EPS1) GO TO 930
*
      KPER=KPER+1
      MAT2(KPER)=ICELL
      IF(LGDEB) THEN
       I=ICELL+1
       IF(I.GT.NCEL) GOTO 58
       GOTO 42
      ENDIF
      ISURC=1
      IETAG=0
      ISURM=6
      SURC(1)=1
      SURC(2)=2
      SURC(3)=3
      SURC(4)=4
      SURC(5)=5
      SURC(6)=6
      IF(LGDIM) THEN
          IF(ICELL.LE.STAIRS(1)) THEN
            IVOLC=ICELL
            IF(.NOT.LGFAC) THEN
              IF(DIRUP) THEN
              IF(IFACE2.LT.8) THEN
                I=VOISIN(IFACE2,ICELL)
                IF(I.GT.NCEL2) THEN
                 KPER=KPER-1
                 ISURC=0
                 GOTO 931
                ELSEIF(IFACE1.LT.7) THEN
                 IF(VOISIN(IFACE1,ICELL).GT.NCEL2)THEN
                   IFF=IFACE1
                   IFACE1=IFACE2
                   IFACE2=IFF
                   KPER=KPER-1
                   ISURC=0
                   GOTO 931
                 ENDIF
                ENDIF
                ISURM=1
                SURC(1)=IFACE2
                IF(IFACE1.LT.7) THEN
                  ISURM=2
                  SURC(2)=IFACE1
                ENDIF
              ELSE
                IF(IPLANZ.EQ.1) THEN
                 KPER=KPER-1
                 ISURC=0
                 GOTO 930
                ENDIF
                IF(IFACE1.LT.7) THEN
                 IF(VOISIN(IFACE1,ICELL).GT.NCEL2)THEN
                   IFF=IFACE1
                   IFACE1=IFACE2
                   IFACE2=IFF
                   KPER=KPER-1
                   ISURC=0
                   GOTO 931
                 ENDIF
                ENDIF
                ISURM=1
                SURC(1)=IFACE1
              ENDIF
              ELSE
              IF(IFACE1.LT.7) THEN
                I=VOISIN(IFACE1,ICELL)
                IF(I.GT.NCEL2) THEN
                 KPER=KPER-1
                 ISURC=0
                 GOTO 931
                ENDIF
                IF(IFACE2.LT.7) THEN
                 IF(VOISIN(IFACE2,ICELL).GT.NCEL2)THEN
                   IFF=IFACE1
                   IFACE1=IFACE2
                   IFACE2=IFF
                   KPER=KPER-1
                   ISURC=0
                   GOTO 931
                 ENDIF
                ENDIF
                ISURM=1
                SURC(1)=IFACE1
                IF(IFACE2.LT.7) THEN
                  ISURM=2
                  SURC(2)=IFACE2
                ENDIF
              ELSE
                 KPER=KPER-1
                 ISURC=0
                 GOTO 930
              ENDIF
              ENDIF
            ELSE
              ISURM=1
              IF(IFACE2.EQ.8) THEN
                ISURM=ISURM-1
              ELSE
                SURC(ISURM)=IFACE2
              ENDIF
              ISURM=ISURM+1
              IF(IFACE1.EQ.7) THEN
                ISURM=ISURM-1
              ELSE
                SURC(ISURM)=IFACE1
              ENDIF
              DO 111 KC=1,MFAC
                IF(FACEM(KC).LT.7) THEN
                 ISURM=ISURM+1
                 SURC(ISURM)=FACEM(KC)
                ENDIF
  111         CONTINUE
              I=VOISIN(SURC(1),IVOLC)
            ENDIF
          ELSE
            IVOLC=ICELL-STAIRS(KPL-1)
            IF(.NOT.LGFAC) THEN
              IF(DIRUP) THEN
              IF(IFACE2.LT.8) THEN
                IF(IFACE2.EQ.7) THEN
                  IFG=IFACE1
                  IFACE1=IFACE2
                  IFACE2=IFG
                ENDIF
                IF(VOISIN(IFACE2,IVOLC).GT.NCEL2) THEN
                 KPER=KPER-1
                 ISURC=0
                 GOTO 931
                ELSEIF(IFACE1.LT.7) THEN
                 IF(VOISIN(IFACE1,IVOLC).GT.NCEL2) THEN
                   IFF=IFACE2
                   IFACE2=IFACE1
                   IFACE1=IFF
                   KPER=KPER-1
                   ISURC=0
                   GOTO 931
                 ENDIF
                ENDIF
                ISURM=1
                SURC(1)=IFACE2
                IF(IFACE1.LT.7) THEN
                 ISURM=2
                 SURC(2)=IFACE1
                ENDIF
              ELSEIF(KPL.EQ.IPLANZ) THEN
                 KPER=KPER-1
                 ISURC=0
                 GOTO 930
              ELSE
              ISURM=1
              SURC(1)=IFACE1
              ENDIF
              ELSE
              IF(IFACE1.LT.7) THEN
                IF(VOISIN(IFACE1,IVOLC).GT.NCEL2) THEN
                 KPER=KPER-1
                 ISURC=0
                 GOTO 931
                ELSEIF(IFACE2.LT.7) THEN
                 IF(VOISIN(IFACE2,IVOLC).GT.NCEL2) THEN
                   IFF=IFACE2
                   IFACE2=IFACE1
                   IFACE1=IFF
                   KPER=KPER-1
                   ISURC=0
                   GOTO 931
                 ENDIF
                ENDIF
                ISURM=1
                SURC(1)=IFACE1
                IF(IFACE2.LT.7) THEN
                  ISURM=2
                  SURC(2)=IFACE2
                ENDIF
              ELSEIF(KPL.EQ.1) THEN
                 KPER=KPER-1
                 ISURC=0
                 GOTO 930
              ELSE
                ISURM=1
                SURC(1)=IFACE2
              ENDIF
              ENDIF
            ELSE
              ISURM=1
              IF(IFACE2.EQ.8) THEN
                ISURM=ISURM-1
              ELSE
                SURC(ISURM)=IFACE2
              ENDIF
              ISURM=ISURM+1
              IF(IFACE1.EQ.7) THEN
                ISURM=ISURM-1
              ELSE
                SURC(ISURM)=IFACE1
              ENDIF
              DO 222 KC=1,MFAC
                IF(FACEM(KC).LT.7) THEN
                 ISURM=ISURM+1
                 SURC(ISURM)=FACEM(KC)
                ENDIF
  222         CONTINUE
            ENDIF
            ISURC=0
            I=IVOLC
            IF(DIRUP) THEN
              IF(KPL.GT.2) THEN
              IETAG=STAIRS(KPL-2)
              ENDIF
            ELSE
              IF(KPL.EQ.IPLANZ) THEN
              IETAG=STAIRS(KPL-1)
              ELSE
              IETAG=STAIRS(KPL)
              ENDIF
            ENDIF
          ENDIF
      ELSE
          IF(.NOT.LGFAC) THEN
            IF(DIRUP) THEN
              IF(VOISIN(IFACE2,ICELL).GT.NCEL2) THEN
                KPER=KPER-1
                ISURC=0
                GOTO 931
              ENDIF
            ELSE
              IF(VOISIN(IFACE1,ICELL).GT.NCEL2) THEN
                KPER=KPER-1
                ISURC=0
                GOTO 931
              ENDIF
            ENDIF
            ISURM=2
            SURC(1)=IFACE2
            SURC(2)=IFACE1
          ELSE
            DO 333 IC=1,MFAC
              ISURM=ISURM+1
              SURC(ISURM)=FACEM(IC)
  333       CONTINUE
          ENDIF
          IVOLC=ICELL
          I=VOISIN(SURC(1),IVOLC)
      ENDIF
      LGPAS0=.TRUE.
      LGPER=.TRUE.
      IF(I.GT.NCEL2) THEN
          I=NCEL+10
          GOTO 56
      ENDIF
      I=I+IETAG
      ICELL0=IVOLC+IETAG
      GOTO 42
 931  CONTINUE
  930 LGDEB=.FALSE.
*
* STOKAGE DES ANCIENS VOLUMES PARCOURUS
*
      IF(ABS(T(IFACE1)-T(IFACE2)).GE.EPS3) THEN
      IF(LGMAT3) THEN
        DO 555 IAC=1,KPER3
          KPER=KPER+1
          MAT2(KPER)=MAT3(IAC)
  555   CONTINUE
      ENDIF
      LGMAT3=.FALSE.
      KPER3=0
      ENDIF
*
*        CALCUL DES PARCOURS OPTIQUE. LA CELLULE PEUT CONTENIR DES
*        CYLINDRES D'AXES DIFFERENTS
*
        KPER=KPER+1
        IF(KPER.GT.NCEL+10) THEN
         GOTO 384
        ENDIF
        MAT2(KPER)=ICELL
*
*  COMPTAGE DES CELLULES DU COIN EN PERIPHERIE DANS MAT2
*
      LGOUT1=.FALSE.
      LGOUTU=.FALSE.
        IF(LGFAC) THEN
          DO 812 MM=1,MFAC1
            IF(FACEM1(MM).LT.7) THEN
              IF(ICELL.LE.STAIRS(1)) THEN
                KMAT2=VOISIN(FACEM1(MM),ICELL)
                IF(KMAT2.LE.NCEL2) THEN
                IF(IFACE1.LT.7) THEN
                IF(VOISIN(IFACE1,ICELL).GT.NCEL2) THEN
                  IPP=IFACE1
                  IFACE1=FACEM1(MM)
                  FACEM1(MM)=IPP
                  IF(IPP.LT.7) THEN
                    KMAT2=VOISIN(IPP,ICELL)
                    IF(KMAT2.LE.LIM2) THEN
                      IF(KMAT2.GE.LIM1)THEN
                        KPER=KPER+1
                        MAT2(KPER)=KMAT2
                      ELSEIF(LGDIM)THEN
                        KPER=KPER+1
                        MAT2(KPER)=KMAT2
                      ENDIF
                    ENDIF
                  ENDIF
                ELSE
                  IF(KMAT2.LE.LIM2) THEN
                    IF(KMAT2.GE.LIM1)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ELSEIF(LGDIM)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ENDIF
                  ELSEIF(LGDIM) THEN
                    LGMAT3=.TRUE.
                    KPER3=KPER3+1
                    MAT3(KPER3)=KMAT2
                  ENDIF
                ENDIF
                ELSE
                  IF(KMAT2.LE.LIM2) THEN
                    IF(KMAT2.GE.LIM1)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ELSEIF(LGDIM)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ENDIF
                  ELSEIF(LGDIM) THEN
                    LGMAT3=.TRUE.
                    KPER3=KPER3+1
                    MAT3(KPER3)=KMAT2
                  ENDIF
                ENDIF
                ENDIF
              ELSE
                ICELL2=ICELL-STAIRS(KPL-1)
                KMAT2=VOISIN(FACEM1(MM),ICELL2)
                IF(KMAT2.LE.NCEL2) THEN
                IF(IFACE1.LT.7) THEN
                IF(VOISIN(IFACE1,ICELL2).GT.NCEL2) THEN
                  IPP=IFACE1
                  IFACE1=FACEM1(MM)
                  FACEM1(MM)=IPP
                  IF(IPP.EQ.7) THEN
                    KPER=KPER+1
                    MAT2(KPER)=ICELL-STAIRS(KPL-1)
                    IF(KPL.GT.2)MAT2(KPER)=MAT2(KPER)+STAIRS(KPL-2)
                  ELSE
                    KMAT2=VOISIN(IPP,ICELL2)
                    IF(KMAT2.LE.LIM2) THEN
                      IF(KMAT2.GE.LIM1)THEN
                        KPER=KPER+1
                        MAT2(KPER)=KMAT2
                      ELSEIF(KPL.EQ.1)THEN
                        KPER=KPER+1
                        MAT2(KPER)=KMAT2
                      ENDIF
                    ENDIF
                  ENDIF
                ELSE
                  IF(KMAT2.LE.LIM2) THEN
                    IF(KMAT2.GE.LIM1)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ELSEIF(KPL.EQ.1)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ENDIF
                  ELSEIF(KPL.EQ.1) THEN
                    LGMAT3=.TRUE.
                    KPER3=KPER3+1
                    MAT3(KPER3)=KMAT2
                  ENDIF
                ENDIF
                ELSE
                  IF(KMAT2.LE.LIM2) THEN
                    IF(KMAT2.GE.LIM1)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ELSEIF(KPL.EQ.1)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ENDIF
                  ELSEIF(KPL.EQ.1) THEN
                    LGMAT3=.TRUE.
                    KPER3=KPER3+1
                    MAT3(KPER3)=KMAT2
                  ENDIF
                ENDIF
                ELSE
                  LGOUT1=.TRUE.
                  KCEL1=ICELL
                  KFACE1=FACEM1(MM)
                ENDIF
              ENDIF
            ENDIF
 812      CONTINUE
          DO 813 MM=1,MFAC2
            IF(FACEM2(MM).LT.7) THEN
              IF(ICELL.LE.STAIRS(1)) THEN
                KMAT2=VOISIN(FACEM2(MM),ICELL)
                IF(KMAT2.LE.NCEL2) THEN
                IF(IFACE2.LT.7) THEN
                IF(VOISIN(IFACE2,ICELL).GT.NCEL2) THEN
                  IPP=IFACE2
                  IFACE2=FACEM2(MM)
                  FACEM2(MM)=IPP
                  IF(IPP.EQ.8) THEN
                    KPER=KPER+1
                    MAT2(KPER)=ICELL+STAIRS(1)
                  ELSE
                    KAM=VOISIN(IPP,ICELL)
                    IF(KAM.LE.LIM2) THEN
                      IF(KAM.GE.LIM1)THEN
                        KPER=KPER+1
                        MAT2(KPER)=KAM
                      ELSEIF(LGDIM) THEN
                        KPER=KPER+1
                        MAT2(KPER)=KAM
                      ENDIF
                    ENDIF
                  ENDIF
                ELSE
                    IF(KMAT2.GE.LIM1)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ELSEIF(LGDIM) THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ENDIF
                ENDIF
                ELSE
                    IF(KMAT2.GE.LIM1)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ELSEIF(LGDIM) THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2
                    ENDIF
                ENDIF
                ELSE
                  KFACE2=FACEM2(MM)
                  KCEL2=ICELL
                  LGOUTU=.TRUE.
                ENDIF
              ELSE
                ICELL2=ICELL-STAIRS(KPL-1)
                KMAT2=VOISIN(FACEM2(MM),ICELL2)
                IF(KMAT2.LE.NCEL2) THEN
                IF(IFACE2.LT.7) THEN
                IF(VOISIN(IFACE2,ICELL2).GT.NCEL2) THEN
                  IPP=IFACE2
                  IFACE2=FACEM2(MM)
                  FACEM2(MM)=IPP
                  IF(IPP.EQ.8) THEN
                    IF(ICELL2.LE.LIM2) THEN
                      IF(ICELL2.GE.LIM1) THEN
                        LGMAT3=.TRUE.
                        KPER3=KPER3+1
                        MAT3(KPER3)=ICELL2+STAIRS(KPL)
                      ELSEIF(KPL.EQ.IPLANZ-1) THEN
                        LGMAT3=.TRUE.
                        KPER3=KPER3+1
                        MAT3(KPER3)=ICELL2+STAIRS(KPL)
                      ENDIF
                    ELSEIF(KPL.EQ.IPLANZ-1)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=ICELL2+STAIRS(KPL)
                    ENDIF
                  ELSE
                    KMA=VOISIN(IPP,ICELL2)
                    IF(KMA.LE.LIM2) THEN
                      IF(KMA.GE.LIM1) THEN
                        LGMAT3=.TRUE.
                        KPER3=KPER3+1
                        MAT3(KPER3)=KMA+STAIRS(KPL-1)
                      ELSEIF(KPL.EQ.IPLANZ) THEN
                        LGMAT3=.TRUE.
                        KPER3=KPER3+1
                        MAT3(KPER3)=KMA+STAIRS(KPL-1)
                      ENDIF
                    ELSEIF(KPL.EQ.IPLANZ)THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMA+STAIRS(KPL-1)
                    ENDIF
                  ENDIF
                ELSE
                  IF(KMAT2.LE.LIM2) THEN
                    IF(KMAT2.GE.LIM1) THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2+STAIRS(KPL-1)
                    ELSEIF(KPL.EQ.IPLANZ) THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2+STAIRS(KPL-1)
                    ENDIF
                  ELSEIF(KPL.EQ.IPLANZ)THEN
                    LGMAT3=.TRUE.
                    KPER3=KPER3+1
                    MAT3(KPER3)=KMAT2+STAIRS(KPL-1)
                  ENDIF
                ENDIF
                ELSE
                  IF(KMAT2.LE.LIM2) THEN
                    IF(KMAT2.GE.LIM1) THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2+STAIRS(KPL-1)
                    ELSEIF(KPL.EQ.IPLANZ) THEN
                      LGMAT3=.TRUE.
                      KPER3=KPER3+1
                      MAT3(KPER3)=KMAT2+STAIRS(KPL-1)
                    ENDIF
                  ELSEIF(KPL.EQ.IPLANZ)THEN
                    LGMAT3=.TRUE.
                    KPER3=KPER3+1
                    MAT3(KPER3)=KMAT2+STAIRS(KPL-1)
                  ENDIF
                ENDIF
                ELSE
                  KFACE2=FACEM2(MM)
                  KCEL2=ICELL
                  LGOUTU=.TRUE.
                ENDIF
              ENDIF
            ENDIF
 813      CONTINUE
          IF(IPOP.LT.1)LGMAT3=.FALSE.
        ENDIF
****
      IF(ICELL.GT.NCEL) THEN
        WRITE(HSMG,'(29HTRKHEX: CELL OVERFLOW (ICELL=,I6,6H NCYL=,I6,
     +  6H NCEL=,I6,1H))') ICELL,NCYL(ICELL),NCEL
        CALL XABORT(HSMG)
      ENDIF
****
      KUP=1
      KDOWN=2
      IVUP=0
      IVDOWN=0
      LTHROU=.FALSE.
      IF(NCYL(ICELL).GT.0) THEN
        IPOS=2*NCEL
        IF(LGDIM) IPOS=3*NCEL
        DO 22 I=1,ICELL-1
           IPOS=IPOS+NCYL(I)
   22   CONTINUE
        KUP=NCYL(ICELL)+1
        KDOWN=KUP+1
        LGFAC8=.TRUE.
        LGFAC7=.TRUE.
        LAUX=0
        DO 25 ICYL=1,NCYL(ICELL)
           IPOS=IPOS+1
           XCYL=REMESH(ICELL)
           YCYL=REMESH(NCEL+ICELL)
           R   =REMESH(IPOS)
           LAUX=LAUX+1
           RAUX(LAUX)=R*R
*
*--TOUS LES CYLINDRES SONT SUPPOSES D'AXES Z
*
             YYC=Y0-YCYL
             BZ=-COS1*XCYL+YYC*COS2
             DEE=YYC*COS1+XCYL*COS2
             DEL=RAUX(LAUX)*AZ-DEE*DEE
             IF(DEL.GT.0.) THEN
               SDEL=SQRT(DEL)
               TP=(-BZ+SDEL)/AZ
               TM=(-BZ-SDEL)/AZ
               IF(LGDIM) THEN
                 ZP=COS3*TP+Z0
                 ZM=COS3*TM+Z0
                  IF(ZP.LT.ZM) THEN
                    ZZZ=ZP
                    ZP=ZM
                    ZM=ZZZ
                    TTT=TP
                    TP=TM
                    TM=TTT
                  ENDIF
                 IF(ZM.GT.Z2.OR.ZP.LT.Z1) GOTO 24
                 LTHROU=.TRUE.
                 IF(ZP.GE.Z2.OR.ABS(ZP-Z2).LE.EPS1) THEN
                   IF(IFACE2.EQ.8)THEN
                     IF(LGFAC8)THEN
                       LGFAC8=.FALSE.
                       IVUP=ICYL
                     ENDIF
                     IF(ABS(TM-T(8)).LT.EPS6) THEN
                        LTHROU=.FALSE.
                        IF(.NOT.LGFAC8)IVUP=ICYL+1
                        GO TO 24
                     ENDIF
                   ENDIF
                   IF(ZM.LE.Z1) THEN
                     IF(IFACE1.EQ.7)THEN
                       IF(LGFAC7)THEN
                         LGFAC7=.FALSE.
                         IVDOWN=ICYL
                       ENDIF
                     ENDIF
                     GOTO 24
                   ENDIF
                   IVDOWN=ICYL
                   IF(ABS(TM-T(7)).LT.EPS4) THEN
                      LGFAC7=.FALSE.
                      GO TO 24
                   ENDIF
                   IVDOWN=ICYL+1
                   KDOWN=KDOWN-1
                   T0(KDOWN)=TM
                   V0(KDOWN)=ICYL+1
                   GO TO 24
                 ELSEIF(ZM.LE.Z1.OR.ABS(ZM-Z1).LE.EPS1) THEN
                   IF(IFACE1.EQ.7)THEN
                      IF(LGFAC7)THEN
                        LGFAC7=.FALSE.
                        IVDOWN=ICYL
                      ENDIF
                      IF(ABS(TP-T(7)).LT.EPS4) THEN
                         LTHROU=.FALSE.
                         IF(.NOT.LGFAC7) IVDOWN=ICYL+1
                         GO TO 24
                      ENDIF
                   ENDIF
                   IVUP=ICYL
                   IF(ABS(TP-T(8)).LT.EPS4) THEN
                      LGFAC8=.FALSE.
                      GO TO 24
                   ENDIF
                   IVUP=ICYL+1
                   KUP=KUP+1
                   T0(KUP)=TP
                   V0(KUP)=ICYL
                   GO TO 24
                 ENDIF
               ENDIF
               IF(LGFAC8)IVUP=ICYL+1
               IF(LGFAC7)IVDOWN=ICYL+1
               KUP=KUP+1
               KDOWN=KDOWN-1
               T0(KUP)=TP
               T0(KDOWN)=TM
               V0(KUP)=ICYL
               V0(KDOWN)=ICYL+1
             ENDIF
   24      CONTINUE
   25   CONTINUE
      ENDIF
*
*        VOLUME A BORDURE HEXAGONALE
*
      KVOL2=NCYL(ICELL)+1
      IF(IFACE2.EQ.8) THEN
        IF(LTHROU) KVOL2=IVUP
      ENDIF
      KVOL1=NCYL(ICELL)+1
      IF(IFACE1.EQ.7) THEN
        IF(LTHROU) KVOL1=IVDOWN
      ENDIF
      KUP=KUP+1
      KDOWN=KDOWN-1
      T0(KUP)=T(IFACE2)
      T0(KDOWN)=T(IFACE1)
      V0(KUP)=KVOL2
      V0(KDOWN)=KVOL1
      IF(COS3.LT.0.) THEN
       T0MIN=T0(KUP)
       T0MAX=T0(KDOWN)
      ELSE
       T0MAX=T0(KUP)
       T0MIN=T0(KDOWN)
      ENDIF
*
*        STOCKAGE POUR SECTEURS
*
      K=0
      IV01=V0(KDOWN)
      IV02=V0(KDOWN)
      LSV=FVOL(ICELL)+NCYL(ICELL)
      IFV=FVOL(ICELL)-1
      JVOUT1=0
      JVOUT2=0
      IVSEC1=0
      IVSEC2=0
      DO 125 I=KDOWN,KUP
         K=K+1
         TSEC(K)=T0(I)
         IV01=MIN(IV01,V0(I))
         IV02=MAX(IV02,V0(I))
         IVS=0
         DO 233 IS=1,V0(I)-1
            IF(NSECT(IFV+IS).GT.1) THEN
              IVS=IVS+6*(NSECT(IFV+IS)-1)
            ELSE
              IVS=IVS+1
            ENDIF
  233    CONTINUE
         IF(NSECT(IFV+V0(I)).GT.1) THEN
           YPT=COS2*T0(I)+Y0
           XPT=COS1*T0(I)
           LGPASS=.TRUE.
*
*  CAS OU LA TRACK ARRIVE SUR LE CENTRE DE LA CELLULE
*
           IF(V0(I).EQ.1) THEN
             IF(ABS(XPT-X).LE.EPS1.AND.ABS(YPT-Y).LE.EPS1) THEN
               IF(K.EQ.1) THEN
                 YPT=COS2*T0(I+1)+Y0
                 XPT=COS1*T0(I+1)
                 LGPASS=.FALSE.
               ELSE
                 YPT=COS2*T0(I-1)+Y0
                 XPT=COS1*T0(I-1)
                 LGPASS=.FALSE.
               ENDIF
             ENDIF
           ENDIF
*
           IF(XPT-X.EQ.0.) THEN
             PENTE=.5*PI
           ELSE
             PENTE=ATAN((YPT-Y)/(XPT-X))
           ENDIF
           IF(XPT.GT.X) THEN
             IF(YPT.LT.Y)PENTE=2*PI-ABS(PENTE)
           ELSE
             IF(YPT.LT.Y) THEN
                PENTE=PI+ABS(PENTE)
             ELSE
                PENTE=PI-ABS(PENTE)
             ENDIF
           ENDIF
           ISV=6*(NSECT(IFV+V0(I))-1)
           JSV=ISV/2
           ANG1=2.*PI/FLOAT(ISV)
           ANG2=PI/FLOAT(JSV)
           ATAN0=0.
           DEL0=PENTE
           ISEC=0
           DO 231 IS1=1,ISV
              ATAN1=FLOAT(IS1)*ANG1
              DEL1=ABS(PENTE-ATAN1)
*
*  2.E-2 CORRESPAND A PEU PRES A 1 DEGRE
*
              IF(DEL0.LE.2.E-2.OR.DEL1.LE.2.E-2) THEN
                IAUX=IS1
                IF(DEL0.LE.2.E-2)IAUX=IS1-1
                IF(IAUX.EQ.0)IAUX=ISV
                IF(IAUX.EQ.JSV) THEN
                  SINR=0.
                  COSR=-1.
                ELSEIF(IAUX.EQ.ISV) THEN
                  SINR=0.
                  COSR=1.
                ELSE
                  SINR=SIN(FLOAT(IAUX)*ANG1)
                  COSR=COS(FLOAT(IAUX)*ANG1)
                ENDIF
                YPH=-(XPT-X)*SINR+(YPT-Y)*COSR
                YPF=0.
                IF(K.EQ.1) THEN
                  YPH2=-(COS1*T0(KUP)-X)*SINR+(COS2*T0(KUP)+Y0-Y)*COSR
                  IF(YPH.GT.YPF) THEN
                    IF(YPH2.GT.YPF) THEN
                      ISEC=IAUX+1
                    ELSE
                      IF(MOD(ISV,IAUX).EQ.0) THEN
                        IF(ISV/IAUX.NE.4) THEN
                          JAUX=IAUX
                          IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                          TGA=TAN(FLOAT(JAUX)*ANG2)
                          TTHYP=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                        ELSE
                          TTHYP=X*COS1I
                        ENDIF
                      ELSE
                        JAUX=IAUX
                        IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                        TGA=TAN(FLOAT(JAUX)*ANG2)
                        TTHYP=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                      ENDIF
                      ISEC=IAUX
                      IF(TTHYP.GT.T0MIN) THEN
                      TT0=ABS(TTHYP-T0(I))
                      IF(TT0.LE.EPS1) THEN
                        ISEC=IAUX
                        TT1=ABS(TTHYP-T0(I+1))
                        JVOUT1=IAUX+1
                        IF(TT1.LE.EPS1) THEN
                          IF(TT1.LT.TT0) THEN
                             ISEC=IAUX+1
                             JVOUT1=0
                             JVOUT2=IAUX
                             IF(JVOUT2.GT.ISV) JVOUT2=JVOUT2-ISV
                          ENDIF
                        ENDIF
                        IF(JVOUT1.GT.ISV) JVOUT1=JVOUT1-ISV
                      ELSE
                        ISEC=IAUX+1
                      ENDIF
                      ENDIF
                    ENDIF
                  ELSEIF(YPH.LT.YPF) THEN
                    IF(YPH2.LT.YPF) THEN
                      ISEC=IAUX
                    ELSE
                      IF(MOD(ISV,IAUX).EQ.0) THEN
                        IF(ISV/IAUX.NE.4) THEN
                          JAUX=IAUX
                          IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                          TGA=TAN(FLOAT(JAUX)*ANG2)
                          TTHYP=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                        ELSE
                          TTHYP=X*COS1I
                        ENDIF
                      ELSE
                        JAUX=IAUX
                        IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                        TGA=TAN(FLOAT(JAUX)*ANG2)
                        TTHYP=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                      ENDIF
                      ISEC=IAUX+1
                      IF(TTHYP.GT.T0MIN)THEN
                      TT0=ABS(T0(I)-TTHYP)
                      IF(TT0.LE.EPS1) THEN
                        ISEC=IAUX+1
                        TT1=ABS(T0(I+1)-TTHYP)
                        JVOUT1=IAUX
                        IF(TT1.LE.EPS1) THEN
                          IF(TT1.LE.TT0) THEN
                             ISEC=IAUX
                             JVOUT1=0
                             JVOUT2=IAUX+1
                             IF(JVOUT2.GT.ISV) JVOUT2=JVOUT2-ISV
                          ENDIF
                        ENDIF
                        IF(JVOUT1.GT.ISV) JVOUT1=JVOUT1-ISV
                      ELSE
                        ISEC=IAUX
                      ENDIF
                      ENDIF
                    ENDIF
                  ELSE
                    ISEC=IAUX
                    IF(YPH2.GT.YPF) ISEC=IAUX+1
                  ENDIF
                ELSE
               YPH2=-(COS1*T0(KDOWN)-X)*SINR+(COS2*T0(KDOWN)+Y0-Y)*COSR
                  IF(YPH.GT.YPF) THEN
                    IF(YPH2.GT.YPF) THEN
                      ISEC=IAUX+1
                    ELSE
                      IF(MOD(ISV,IAUX).EQ.0) THEN
                        IF(ISV/IAUX.NE.4) THEN
                          JAUX=IAUX
                          IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                          TGA=TAN(FLOAT(JAUX)*ANG2)
                          TTHYP=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                        ELSE
                          TTHYP=X*COS1I
                        ENDIF
                      ELSE
                        JAUX=IAUX
                        IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                        TGA=TAN(FLOAT(JAUX)*ANG2)
                        TTHYP=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                      ENDIF
                      ISEC=IAUX
                      IF(TTHYP.LE.T0MAX) THEN
                      IF(LGPASS) THEN
                      IF(NCYL(ICELL).GT.1) THEN
                        XAUX=X-COS1*TTHYP
                        YAUX=Y-COS2*TTHYP-Y0
                        DIST=XAUX*XAUX+YAUX*YAUX
                        ISEC=IAUX
                        IF(V0(I).GT.NCYL(ICELL)) THEN
                          IF(DIST.LE.RAUX(NCYL(ICELL))) GOTO 331
                        ELSEIF(V0(I).EQ.1) THEN
                          IF(DIST.GE.RAUX(1)) GOTO 331
                        ELSE
                          IF(DIST.GE.RAUX(V0(I))) GOTO 331
                          IF(DIST.LE.RAUX(V0(I)-1)) GOTO 331
                        ENDIF
                      ENDIF
                      ENDIF
                      TT0=ABS(T0(I)-TTHYP)
                      IF(TT0.LE.EPS1) THEN
                        ISEC=IAUX
                        TT1=ABS(T0(I-1)-TTHYP)
                        JVOUT2=IAUX+1
                        IF(TT1.LE.EPS1) THEN
                          IF(TT1.LT.TT0) THEN
                            ISEC=IAUX+1
                            JVOUT2=0
                            JVOUT1=IAUX
                            IF(JVOUT1.GT.ISV) JVOUT1=JVOUT1-ISV
                          ENDIF
                        ENDIF
                        IF(JVOUT2.GT.ISV) JVOUT2=JVOUT2-ISV
                      ELSE
                        ISEC=IAUX+1
                      ENDIF
                      ENDIF
                    ENDIF
                  ELSEIF(YPH.LT.YPF) THEN
                    IF(YPH2.LT.YPF) THEN
                      ISEC=IAUX
                    ELSE
                      IF(MOD(ISV,IAUX).EQ.0) THEN
                        IF(ISV/IAUX.NE.4) THEN
                          JAUX=IAUX
                          IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                          TGA=TAN(FLOAT(JAUX)*ANG2)
                          TTHYP=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                        ELSE
                          TTHYP=X*COS1I
                        ENDIF
                      ELSE
                        JAUX=IAUX
                        IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                        TGA=TAN(FLOAT(JAUX)*ANG2)
                        TTHYP=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                      ENDIF
                      ISEC=IAUX+1
                      IF(TTHYP.LT.T0MAX) THEN
                      IF(LGPASS) THEN
                      IF(NCYL(ICELL).GT.1) THEN
                        XAUX=X-COS1*TTHYP
                        YAUX=Y-COS2*TTHYP-Y0
                        DIST=XAUX*XAUX+YAUX*YAUX
                        ISEC=IAUX+1
                        IF(V0(I).GT.NCYL(ICELL)) THEN
                          IF(DIST.LE.RAUX(NCYL(ICELL))) GOTO 331
                        ELSEIF(V0(I).EQ.1) THEN
                          IF(DIST.GE.RAUX(1)) GOTO 331
                        ELSE
                          IF(DIST.GE.RAUX(V0(I))) GOTO 331
                          IF(DIST.LE.RAUX(V0(I)-1)) GOTO 331
                        ENDIF
                      ENDIF
                      ENDIF
                      TT0=ABS(T0(I)-TTHYP)
                      IF(TT0.LE.EPS1) THEN
                        ISEC=IAUX+1
                        TT1=ABS(T0(I-1)-TTHYP)
                        JVOUT2=IAUX
                        IF(TT1.LE.EPS1) THEN
                          IF(TT1.LE.TT0) THEN
                            ISEC=IAUX
                            JVOUT2=0
                            JVOUT1=IAUX+1
                            IF(JVOUT1.GT.ISV) JVOUT1=JVOUT1-ISV
                          ENDIF
                        ENDIF
                        IF(JVOUT2.GT.ISV) JVOUT2=JVOUT2-ISV
                      ELSE
                        ISEC=IAUX
                      ENDIF
                      ENDIF
                    ENDIF
                  ELSE
                    ISEC=IAUX
                    IF(YPH2.GT.YPF) ISEC=IAUX+1
                  ENDIF
                ENDIF
  331           CONTINUE
                IF(ISEC.GT.ISV)ISEC=ISEC-ISV
                GOTO 232
              ELSEIF(PENTE.LT.ATAN1) THEN
                ISEC=IS1
                IAUX=IS1-1
                IF(IAUX.EQ.0)IAUX=ISV
                IF(MOD(ISV,IAUX).EQ.0) THEN
                  IF(ISV/IAUX.NE.4) THEN
                    JAUX=IAUX
                    IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                    TGA=TAN(FLOAT(JAUX)*ANG2)
                    TTHYP0=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                  ELSE
                    TTHYP0=X*COS1I
                  ENDIF
                ELSE
                  JAUX=IAUX
                  IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                  TGA=TAN(FLOAT(JAUX)*ANG2)
                  TTHYP0=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                ENDIF
                IF(ABS(T0(I)-TTHYP0).LE.EPS1) THEN
                  IF(IAUX.EQ.JSV) THEN
                    SINR=0.
                    COSR=-1.
                  ELSEIF(IAUX.EQ.ISV) THEN
                    SINR=0.
                    COSR=1.
                  ELSE
                    SINR=SIN(FLOAT(IAUX)*ANG1)
                    COSR=COS(FLOAT(IAUX)*ANG1)
                  ENDIF
                  IF(K.EQ.1) THEN
                   YPH2=-(COS1*T0(KUP)-X)*SINR+(COS2*T0(KUP)+Y0-Y)*COSR
                  ELSE
               YPH2=-(COS1*T0(KDOWN)-X)*SINR+(COS2*T0(KDOWN)+Y0-Y)*COSR
                  ENDIF
                  IF(YPH2.LE.0.)ISEC=IS1-1
                  GOTO 232
                ENDIF
                IAUX=IS1
                IF(MOD(ISV,IAUX).EQ.0) THEN
                  IF(ISV/IAUX.NE.4) THEN
                    JAUX=IAUX
                    IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                    TGA=TAN(FLOAT(JAUX)*ANG2)
                    TTHYP1=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                  ELSE
                    TTHYP1=X*COS1I
                  ENDIF
                ELSE
                  JAUX=IAUX
                  IF(IAUX.GT.JSV)JAUX=IAUX-JSV
                  TGA=TAN(FLOAT(JAUX)*ANG2)
                  TTHYP1=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
                ENDIF
                IF(ABS(T0(I)-TTHYP1).LE.EPS1) THEN
                  IF(IAUX.EQ.JSV) THEN
                    SINR=0.
                    COSR=-1.
                  ELSEIF(IAUX.EQ.ISV) THEN
                    SINR=0.
                    COSR=1.
                  ELSE
                    SINR=SIN(FLOAT(IAUX)*ANG1)
                    COSR=COS(FLOAT(IAUX)*ANG1)
                  ENDIF
                  IF(K.EQ.1) THEN
                   YPH2=-(COS1*T0(KUP)-X)*SINR+(COS2*T0(KUP)+Y0-Y)*COSR
                  ELSE
               YPH2=-(COS1*T0(KDOWN)-X)*SINR+(COS2*T0(KDOWN)+Y0-Y)*COSR
                  ENDIF
                  IF(YPH2.GT.0.)ISEC=IS1+1
                  IF(ISEC.GT.ISV)ISEC=ISEC-ISV
                  GOTO 232
                ENDIF
                GOTO 232
              ENDIF
              ATAN0=ATAN1
              DEL0=DEL1
  231      CONTINUE
           CALL XABORT('TRKHEX: ALGORITHME FAILURE -D')
  232    CONTINUE
       ELSE
         ISEC=1
       ENDIF
       VSEC(K)=IVS+ISEC
       IF(JVOUT1.NE.0) THEN
         IVSEC1=IVS+JVOUT1
         JVOUT1=0
       ENDIF
       IF(JVOUT2.NE.0) THEN
         IVSEC2=IVS+JVOUT2
         JVOUT2=0
       ENDIF
  125 CONTINUE
      KMAX=K
      JCC0=IFV
*
*  CAS OU IL N'Y A PAS D'INTERSECTION AVEC LES SECTEURS
*
      IF(KMAX.EQ.2) THEN
      IF(VSEC(1).EQ.VSEC(KMAX)) GOTO 527
      ENDIF
*
      DO 227 JC=IV01,IV02
      JCC=JCC0+JC
      IF(NSECT(JCC).GT.1) THEN
*
*             INTERSECTION DE LA DIRECTION AVEC LES SECTEURS
*
        JSEC=3*(NSECT(JCC)-1)
        JSEC2=JSEC+JSEC
        ANG0=PI/FLOAT(JSEC)
        IVS0=VSEC(1)
        IVS1=VSEC(1)-1
        IF(IVS1.EQ.0)IVS1=JSEC
        DO 225 J=1,JSEC
           IF(MOD(JSEC,J).EQ.0) THEN
             IF(JSEC/J.NE.2) THEN
               TGA=TAN(FLOAT(J)*ANG0)
               THYP=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
             ELSE
               THYP=X*COS1I
             ENDIF
           ELSE
             TGA=TAN(FLOAT(J)*ANG0)
             THYP=(-X*TGA+Y-Y0)/(COS2-COS1*TGA)
           ENDIF
           IF(THYP.LE.T0MIN.OR.THYP.GE.T0MAX) GOTO 223
*
* LES PARCOURS NULS AU CENTRE DE LA CELLULE NE SONT PAS ELEMINES
*
*          IF(LGDIM) THEN
           IZ0=1
           IF(IFACE1.EQ.7) THEN
             IF(ABS(THYP-T0MIN).LE.EPS1) GOTO 223
             IZ0=2
           ENDIF
           IF(ABS(THYP-T0MAX).LE.EPS1) GOTO 223
           DO 171 IZZ=IZ0,KMAX-1
              IF(ABS(THYP-TSEC(IZZ)).LE.EPS1) THEN
                IF(JC.EQ.1) THEN
                  IF(J.EQ.IVS0.OR.J.EQ.IVS1) THEN
                    IF(VSEC(1).NE.VSEC(2)) GOTO 114
                  ENDIF
                ENDIF
                GOTO 223
              ENDIF
  171      CONTINUE
  114      CONTINUE
*
*  CE TESTE EST NECESSAIRE SI TOUTES LES ZONES NE SONT PAS SECTORISEES
*  ON SUPPOSE ICI QU'ON A AFFAIRE A DES CYLINDRES CONCENTRIQUES
*
           IF(NCYL(ICELL).GT.1) THEN
             XAUX=X-COS1*THYP
             YAUX=Y-COS2*THYP-Y0
             DIST=XAUX*XAUX+YAUX*YAUX
             IF(JC.GT.NCYL(ICELL)) THEN
               IF(DIST.LE.RAUX(NCYL(ICELL))) THEN
                 IF(KMAX.EQ.2) THEN
                 IF(ABS(DIST-RAUX(NCYL(ICELL))).GT.EPS5) THEN
                   IF(VSEC(1).EQ.VSEC(2)) GOTO 223
                 ENDIF
                 ELSE
                   GOTO 223
                 ENDIF
               ENDIF
             ELSEIF(JC.EQ.1) THEN
               IF(DIST.GE.RAUX(1)) GOTO 223
             ELSE
               IF(DIST.GE.RAUX(JC)) GOTO 223
               IF(DIST.LE.RAUX(JC-1)) GOTO 223
             ENDIF
           ENDIF
           LMIN=2
           LMAX=KMAX
           LSTEP=1
           LL=0
           L0=2
           IF(COS3.LT.0.) THEN
             LMAX=1
             LMIN=KMAX-1
             LSTEP=-1
             LL=1
             L0=1
           ENDIF
           DO 220 K=LMIN,LMAX,LSTEP
              IF(THYP.LT.TSEC(K)) THEN
                IF(K.EQ.L0)THEN
                  IF(VSEC(1).EQ.VSEC(2)) GOTO 223
                ENDIF
                DO 241 L=KMAX,K+LL,-1
                   TSEC(L+1)=TSEC(L)
                   VSEC(L+1)=VSEC(L)
  241           CONTINUE
                TSEC(K+LL)=THYP
                KMAX=KMAX+1
                IZONE=1
                DO 210 N=KDOWN,KUP
                   IF(THYP.LE.T0(N)) THEN
                     IZONE=V0(N)
                     GO TO 211
                   ENDIF
  210           CONTINUE
                CALL XABORT('TRKHEX: IMPOSSIBLE CASE ')
  211           XT=COS1*THYP
                YT=COS2*THYP+Y0
*
* TRAITEMENT DU PROBLEME DU CENTRE DE LA CELLULE SECTORISEE
*
                IF(JC.EQ.1) THEN
                  IF(ABS(XT-X).LE.EPS1.AND.ABS(YT-Y).LE.EPS1) THEN
                    IF(K.EQ.2) THEN
                      VSEC(2)=VSEC(1)
                      GOTO 224
                    ENDIF
                    XT=COS1*T0MIN
                    YT=COS2*T0MIN+Y0
                    IF(ABS(XT-X).LE.EPS1.AND.ABS(YT-Y).LE.EPS1) THEN
                      CALL XABORT('TRKHEX: TRACK PARALLEL TO OZ 2')
                    ENDIF
                  ENDIF
                ENDIF
*
                IAUX=J
                IF(J.EQ.JSEC) THEN
                  SINR=0.
                  COSR=-1.
                  IF(XT.GT.X) THEN
                    IAUX=JSEC2
                    COSR=1.
                  ENDIF
                ELSE
                  IF(YT.LT.Y) IAUX=J+JSEC
                  SINR=SIN(IAUX*ANG0)
                  COSR=COS(IAUX*ANG0)
                ENDIF
                YRH=0.
                YRF=-(XDR(IFACE2)-X)*SINR+(YDR(IFACE2)-Y)*COSR
                IF(YRF.GT.YRH)THEN
                   ISECTR=IAUX
                ELSE
                   ISECTR=IAUX+1
                ENDIF
                IF(ISECTR.GT.JSEC2)ISECTR=ISECTR-JSEC2
                ISS=0
                DO 240 ISC=1,IZONE-1
                   IF(NSECT(IFV+ISC).GT.1) THEN
                     ISS=ISS+6*(NSECT(IFV+ISC)-1)
                   ELSE
                     ISS=ISS+1
                   ENDIF
  240           CONTINUE
                IVAUX=ISS+ISECTR
                VSEC(K+LL)=IVAUX
                GOTO 223
              ENDIF
  220      CONTINUE
  223      CONTINUE
  225   CONTINUE
      ENDIF
  224 CONTINUE
  227 CONTINUE
  527 CONTINUE
      IF(VSEC(1).NE.VSEC(2)) THEN
         IF(ABS(TSEC(1)-TSEC(2)).LE.EPS1) THEN
         IF(VSEC(1).EQ.VSEC(3))THEN
           DO 528 ISSX=2,KMAX-1
              VSEC(ISSX)=VSEC(ISSX+1)
              TSEC(ISSX)=TSEC(ISSX+1)
  528      CONTINUE
           KMAX=KMAX-1
           GOTO 529
         ENDIF
         ENDIF
      ENDIF
  529 CONTINUE
*
*  CALCUL DES VOLUMES DE SORTIES
*
      LSV=NVOL(ICELL)-1
      DO 228 JS=1,NCYL(ICELL)+1
         IF(NSECT(IFV+JS).GT.1) THEN
            LSV=LSV+6*(NSECT(IFV+JS)-1)
         ELSE
            LSV=LSV+1
         ENDIF
  228 CONTINUE
      IFV=NVOL(ICELL)-1
*
*       PARCOURS OPTIQUES ET VOLUMES ASSOCIES
*
      IF(DIRUP) THEN
        I1=1
        I2=KMAX-1
        I3=1
        KPOP=IPOP+1
      ELSE
        I2=1
        I1=KMAX-1
        I3=-1
        KPOP=IPOP+KMAX-1
      ENDIF
      DO 230 I=I1,I2,I3
         IPOP=IPOP+1
         POP(IPOP)=ABS(TSEC(I+1)-TSEC(I))
         IVOL=VSEC(I+1)+IFV
         IF(IVOL.GT.LSV) IVOL=IVOL-LSV+FVOL(ICELL)-1
         MAT(IPOP)=IVOL
  230 CONTINUE
      IMAT2=MAT(IPOP)
 739  CONTINUE
      IF(START) THEN
        IFDOWN=IFACE1
        ICDOWN=ICELL
        IMATD=MAT(KPOP)
        KFACED=KFACE1
        KCELD=KCEL1
        LGOUTD=LGOUT1
        MFACD=MFAC1
        MFACU=MFAC2
        IFACEU=IFACE2
        KPLD=KPL
        DO 229 IXX=1,MFACD
          FACEMD(IXX)=FACEM1(IXX)
 229    CONTINUE
        DO 329 IXX=1,MFACU
          FACEMU(IXX)=FACEM2(IXX)
 329    CONTINUE
      ENDIF
*
*       RECHERCHE DES PROCHAINS VOLUMES A TRACKER
*
      IF(DIRUP) THEN
        START=.FALSE.
        IVOUT=ICELL
        IMAT=IMAT2
        IF(IFACE2.LE.6) THEN
          ICELLX=ICELL
          IF(ICELL.LE.STAIRS(1)) THEN
             ICELL=VOISIN(IFACE2,ICELL)
             IF(ICELL.GT.NCEL2) THEN
               IFACE=IFACE2
               GOTO 35
             ENDIF
             GOTO 31
          ENDIF
          KPL=INT(AINT(REAL(ICELL)/REAL(STAIRS(1))))
          IF(MOD(ICELL,STAIRS(1)).NE.0)KPL=KPL+1
          ICELL2=ICELL-STAIRS(KPL-1)
          ICELLX=ICELL2
          ICELL2=VOISIN(IFACE2,ICELL2)
          IF(ICELL2.GT.NCEL2) THEN
            IFACE=IFACE2
            GOTO 35
          ENDIF
            ICELL=STAIRS(KPL-1)+ICELL2
        ELSE
           ICELLX=ICELL
           IF(.NOT.LGPG1) THEN
             IFACE=IFACE2
             GOTO 35
           ENDIF
           IF(ICELL.LE.STAIRS(1)) THEN
             ICELL=ICELL+STAIRS(1)
             GOTO 31
           ENDIF
           KPL=INT(AINT(REAL(ICELL)/REAL(STAIRS(1))))
           IF(MOD(ICELL,STAIRS(1)).NE.0)KPL=KPL+1
           ICELLX=ICELL-STAIRS(KPL-1)
           ICELL=ICELLX+STAIRS(KPL)
        ENDIF
   31   IF(ICELL.GT.NCEL) THEN
          IFACE=IFACE2
          GOTO 35
        ENDIF
        KPER=KPER+1
        MAT2(KPER)=ICELL
        IF(KPER-1.GE.3) THEN
          IF(MAT2(KPER-1).EQ.ICELL) THEN
            GOTO 631
          ELSEIF(MAT2(KPER-2).EQ.ICELL) THEN
            GOTO 631
          ELSEIF(MAT2(KPER-3).EQ.ICELL) THEN
            GOTO 631
          ENDIF
          GOTO 632
  631     CONTINUE
        IF(LGPAS0) GOTO 56
        ISURC=1
        IETAG=0
        ISURM=6
        SURC(1)=1
        SURC(2)=2
        SURC(3)=3
        SURC(4)=4
        SURC(5)=5
        SURC(6)=6
        IF(LGDIM) THEN
          IF(ICELL.LT.STAIRS(1)) THEN
            IVOLC=ICELL
            I=VOISIN(1,IVOLC)
          ELSE
            KPL=INT(AINT(REAL(ICELL)/REAL(STAIRS(1))))
            IF(MOD(ICELL,STAIRS(1)).NE.0)KPL=KPL+1
            IVOLC=ICELL-STAIRS(KPL-1)
            IF(KPL.GT.2) THEN
            IETAG=STAIRS(KPL-2)
            ENDIF
            I=VOISIN(1,IVOLC)
          ENDIF
        ELSE
          IVOLC=ICELL
          I=VOISIN(1,IVOLC)
        ENDIF
        LGPAS0=.TRUE.
        IF(I.GT.NCEL2) THEN
          I=NCEL+10
          GOTO 56
        ENDIF
        I=I+IETAG
        LGPER=.TRUE.
        GOTO 42
        ENDIF
  632   CONTINUE
        X=REMESH(ICELL)
        Y=REMESH(NCEL+ICELL)
        LGPAS0=.FALSE.
        GOTO 10
      ENDIF
*
   30 IVOUT=ICELL
      IMAT=IMATD
      START=.TRUE.
      IF(IFDOWN.LE.6) THEN
          IF(ICELL.LE.STAIRS(1)) THEN
            ICELLD=ICELL
            ICELL=VOISIN(IFDOWN,ICELL)
            IF(ICELL.GT.NCEL2) THEN
              IFACE=IFDOWN
              IVOUT2=IVOUT
              LGSTOR=.TRUE.
              GOTO 34
            ENDIF
            GOTO 33
          ENDIF
          KPL=INT(AINT(REAL(ICELL)/REAL(STAIRS(1))))
          IF(MOD(ICELL,STAIRS(1)).NE.0) KPL=KPL+1
          ICELL2=ICELL-STAIRS(KPL-1)
          ICELLD=ICELL2
          ICELL2=VOISIN(IFDOWN,ICELL2)
          IF(ICELL2.GT.NCEL2) THEN
            IFACE=IFDOWN
            IVOUT2=IVOUT
            LGSTOR=.TRUE.
            GOTO 34
          ENDIF
            ICELL=STAIRS(KPL-1)+ICELL2
          GOTO 33
      ELSE
          ICELLD=ICELL
          IF(ICELL.LE.STAIRS(1)) THEN
            IFACE=IFDOWN
            IVOUT2=IVOUT
            LGSTOR=.TRUE.
            GOTO 34
          ENDIF
          KPL=INT(AINT(REAL(ICELL)/REAL(STAIRS(1))))
          IF(MOD(ICELL,STAIRS(1)).NE.0) KPL=KPL+1
          ICELL=ICELL-STAIRS(KPL-1)
          ICELLD=ICELL
          IF(KPL.GT.2)ICELL=ICELL+STAIRS(KPL-2)
      ENDIF
   33 IF(ICELL.GT.NCEL) THEN
        IFACE=IFDOWN
        IVOUT2=IVOUT
        LGSTOR=.TRUE.
        GOTO 35
      ENDIF
      KPER=KPER+1
      MAT2(KPER)=ICELL
        IF(KPER-1.GE.3) THEN
          IF(MAT2(KPER-1).EQ.ICELL) THEN
            GOTO 635
          ELSEIF(MAT2(KPER-2).EQ.ICELL) THEN
            GOTO 635
          ELSEIF(MAT2(KPER-3).EQ.ICELL) THEN
            GOTO 635
          ENDIF
          GOTO 637
  635     CONTINUE
        IF(LGPAS0) GOTO 56
        ISURC=1
        IETAG=0
        ISURM=6
        SURC(1)=1
        SURC(2)=2
        SURC(3)=3
        SURC(4)=4
        SURC(5)=5
        SURC(6)=6
        IF(LGDIM) THEN
          IF(ICELL.LT.STAIRS(1)) THEN
            IVOLC=ICELL
            I=VOISIN(1,IVOLC)
          ELSE
            KPC=INT(AINT(REAL(ICELL)/REAL(STAIRS(1))))
            IF(MOD(ICELL,STAIRS(1)).NE.0)KPC=KPC+1
            IVOLC=ICELL-STAIRS(KPC-1)
            IF(KPC.EQ.IPLANZ) THEN
            IETAG=STAIRS(KPC-1)
            ELSE
            IETAG=STAIRS(KPC)
            ENDIF
            I=VOISIN(1,IVOLC)
          ENDIF
        ELSE
          IVOLC=ICELL
          I=VOISIN(1,IVOLC)
        ENDIF
        LGPAS0=.TRUE.
        IF(I.GT.NCEL2) THEN
          I=NCEL+10
          GOTO 56
        ENDIF
        I=I+IETAG
        LGPER=.TRUE.
        GOTO 42
        ENDIF
  637 IF(ICELL+NCEL.GT.MESH) 
     > CALL XABORT('TRKHEX: ALGORITHME FAILURE -E')
      X=REMESH(ICELL)
      Y=REMESH(NCEL+ICELL)
      LGDIR=.FALSE.
      LGPAS0=.FALSE.
      GOTO 10
*
   34 CONTINUE
*
*       RECHERCHE DE LA FACE DE SORTIE
*
   35 CONTINUE
      IMAT0=IMAT
      NFAC=1
      NFXX=1
*
      IF(LGMAT3) THEN
        DO 666 IAC=1,KPER3
          KPER=KPER+1
          MAT2(KPER)=MAT3(IAC)
  666   CONTINUE
      ENDIF
      LGMAT3=.FALSE.
      KPER3=0
      PASSD=.FALSE.
      PASSU=.FALSE.
      IF(DIRUP) THEN
        IYY=0
        IF(LGFAC) THEN
        DO 335 IXX=1,MFAC2
          IF(FACEM2(IXX).LE.6) THEN
            ICXX=VOISIN(FACEM2(IXX),ICELLX)
            IF(ICXX.GT.NCEL2) THEN
              IYY=IYY+1
              FACES(IYY)=FACEM2(IXX)
            ENDIF
          ELSE
            IF(LGPG1) THEN
              IF(IVOUT.LE.STAIRS(IPLANZ-1)) THEN
                PASSU=.TRUE.
                GOTO 435
              ENDIF
            ENDIF
            IYY=IYY+1
            FACES(IYY)=FACEM2(IXX)
          ENDIF
  435   CONTINUE
  335   CONTINUE
        ENDIF
        NFACX=IYY+1
      ELSE
        IYY=0
        DO 336 IXX=1,MFACD
          IF(FACEMD(IXX).LE.6) THEN
            ICXX=VOISIN(FACEMD(IXX),ICELLD)
            IF(ICXX.GT.NCEL2) THEN
              IYY=IYY+1
              FACES(IYY)=FACEMD(IXX)
            ENDIF
          ELSE
            IF(IVOUT.LE.STAIRS(1))THEN
              IYY=IYY+1
              FACES(IYY)=FACEMD(IXX)
              PASSD=.TRUE.
            ENDIF
          ENDIF
  336   CONTINUE
        NFACX=IYY+1
*
      ENDIF
      ISURF1=0
      ISURF2=0
      IF(IVOUT.LE.STAIRS(1)) THEN
        ISTAIR=0
        IFACST=-1
        IF(LGDIM) THEN
           ISS=0
           DO 300 IW=0,NCYL(NCEL2)
             IF(NSECT(FVOL(NCEL2)+IW).GT.1)THEN
                ISS=6*(NSECT(FVOL(NCEL2)+IW)-1)+ISS
             ELSE
                ISS=1+ISS
             ENDIF
  300      CONTINUE
           IFACST=FACB(NCEL2)+ISS-1
        ENDIF
        IMAX=1
        GOTO 37
      ENDIF
      I=INT(AINT(REAL(IVOUT)/REAL(STAIRS(1))))
      IF(MOD(IVOUT,STAIRS(1)).NE.0)I=I+1
      ISTAIR=STAIRS(I-1)
      IFACST=FACST(I-1)-1
      IMAX=I
   37 KFACE=IVOUT-IFONC(NCOUR,0)-ISTAIR
      IMAT=IMAT0
      IFOUT2=2*KFACE+2+IFACST
      IF(IFACE.EQ.7) THEN
        IMAM=1
        IF(LPOP.NE.IPOP)IMAM=IPOP
        IFOUT2=MAT(IMAM)
      ELSEIF(IFACE.EQ.8) THEN
        IPP=IVOUT
        IF(LGPG1)IPP=IVOUT-STAIRS(IPLANZ-1)
        IFOUT2=FACB(NCEL2+IPP)+VSEC(KMAX)
      ELSE
      IF(NCOUR.EQ.1) THEN
        IFOUT2=IFACE
        IF(LGDIM)IFOUT2=IFOUT2+IFACST+1
      ELSEIF(IVOUT.LE.(IFONC(NCOUR,1)+ISTAIR)) THEN
        IF(IFACE.EQ.2) THEN
          IFOUT2=IFOUT2+1
        ELSEIF(IFACE.EQ.6) THEN
          IFOUT2=FACST(IMAX)
        ELSEIF(IFACE.EQ.3) THEN
          IFOUT2=IFOUT2+2
        ENDIF
      ELSEIF(IVOUT.LE.(IFONC(NCOUR,2)+ISTAIR)) THEN
        IFOUT2=2*KFACE+3+IFACST
        IF(IFACE.EQ.3) THEN
          IFOUT2=IFOUT2+1
        ELSEIF(IFACE.EQ.4) THEN
          IFOUT2=IFOUT2+2
        ENDIF
      ELSEIF(IVOUT.LE.(IFONC(NCOUR,3)+ISTAIR)) THEN
        IFOUT2=2*KFACE+4+IFACST
        IF(IFACE.EQ.4) THEN
          IFOUT2=IFOUT2+1
        ELSEIF(IFACE.EQ.5) THEN
          IFOUT2=IFOUT2+2
        ENDIF
      ELSEIF(IVOUT.LE.(IFONC(NCOUR,4)+ISTAIR)) THEN
        IFOUT2=2*KFACE+5+IFACST
        IF(IFACE.EQ.5) THEN
          IFOUT2=IFOUT2+1
        ELSEIF(IFACE.EQ.6) THEN
          IFOUT2=IFOUT2+2
        ENDIF
      ELSEIF(IVOUT.LE.(IFONC(NCOUR,5)+ISTAIR)) THEN
        IFOUT2=2*KFACE+6+IFACST
        IF(IFACE.EQ.6) THEN
          IFOUT2=IFOUT2+1
        ELSEIF(IFACE.EQ.1) THEN
          IFOUT2=IFOUT2+2
        ENDIF
      ELSEIF(IVOUT.LE.STAIRS(IMAX)) THEN
        IFOUT2=2*KFACE+7+IFACST
        IF(IFACE.EQ.1)IFOUT2=IFOUT2+1
      ELSE
        CALL XABORT('TRKHEX: ALGORITHME FAILURE -F')
      ENDIF
*
*  PRISE EN COMPTE DES SECTEURS
*
      KSECT=NSECT(FVOL(IVOUT)+NCYL(IVOUT))
      IF(KSECT.GT.2) THEN
        IXY=0
        DO 371 IXX=0,NCYL(IVOUT)-1
           IF(NSECT(FVOL(IVOUT)+IXX).GT.1) THEN
             IXY=IXY+6*(NSECT(FVOL(IVOUT)+IXX)-1)
           ELSE
             IXY=IXY+1
           ENDIF
  371   CONTINUE
        JMIN=(IFACE-1)*(KSECT-1)+IXY+NVOL(IVOUT)-1
        JMAX=KSECT-1+JMIN
        IF(IMAT.GT.JMAX) THEN
          IF(IFACE.EQ.1) THEN
            IMAT=NVOL(IVOUT)+IXY
          ELSE
            IMAT=IMAT-1
          ENDIF
        ELSEIF(IMAT.LE.JMIN) THEN
          IF(IFACE.EQ.6) THEN
            IMAT=NVOL(IVOUT)+IXY-1+6*(KSECT-1)
          ELSE
            IMAT=IMAT+1
          ENDIF
        ENDIF
        IMAT=IMAT-(IFACE-1)*(KSECT-1)-IXY-NVOL(IVOUT)+1
        IFOUT2=SURB(IFOUT2)+IMAT
      ELSEIF(SECTOR) THEN
        IFOUT2=SURB(IFOUT2)+1
      ENDIF
      ENDIF
*------ DISTRIBUTION DES NEUTRONS SUR DES FACES DES PLANS SUPERIEURS
*       OU INFERIEURS AUTRE QUE LE PLAN CONSIDERE
      IF(LGPG1) THEN
        IF(DIRUP) THEN
          IF(PASSU) THEN
            IF(IVOUT.LE.STAIRS(1)) THEN
              LFOUT2=FACST(1)+IFOUT2-STAIRS(1)
              LVOUT=STAIRS(1)+IVOUT
              GOTO 477
            ELSEIF(IVOUT.LE.STAIRS(IPLANZ-1)) THEN
              LFOUT2=FACST(IMAX)+IFOUT2-FACST(IMAX-1)
              LVOUT=STAIRS(IMAX)+IVOUT-STAIRS(IMAX-1)
              GOTO 477
            ENDIF
          ENDIF
          GOTO 478
        ELSE
          IF(PASSD) THEN
            IF(IVOUT.GT.STAIRS(1)) THEN
              IF(IVOUT.GT.STAIRS(2)) THEN
                LFOUT2=IFOUT2-FACST(IMAX-1)+FACST(IMAX-2)
                LVOUT=IVOUT-STAIRS(IMAX-1)+STAIRS(IMAX-2)
                GOTO 477
              ELSE
                LFOUT2=IFOUT2-FACST(1)+STAIRS(1)
                LVOUT=IVOUT-STAIRS(1)
                GOTO 477
              ENDIF
            ENDIF
          ENDIF
          GOTO 478
        ENDIF
 477    IF(SECTOR) THEN
          MSECT=NSECT(FVOL(LVOUT)+NCYL(LVOUT))
          IF(MSECT.EQ.KSECT) THEN
            LFOUT2=SURB(LFOUT2)+IMAT
          ELSEIF(MSECT.LT.KSECT) THEN
            LFOUT2=SURB(LFOUT2)+NINT(REAL(MSECT*IMAT)/REAL(KSECT))
          ELSE
            KFXX=NINT(REAL(MSECT*IMAT)/REAL(KSECT))
            IF(KFXX.EQ.1) THEN
              LFOUT2=SURB(LFOUT2)+1
              SURFX(NFXX)=LFOUT2+1
              NFXX=NFXX+1
            ELSEIF(KFXX.EQ.MSECT)THEN
              LFOUT2=SURB(LFOUT2)+KFXX
              SURFX(NFXX)=LFOUT2-1
              NFXX=NFXX+1
            ELSE
              LFOUT2=SURB(LFOUT2)+KFXX
              SURFX(NFXX)=LFOUT2-1
              SURFX(NFXX+1)=LFOUT2+1
              NFXX=NFXX+2
            ENDIF
          ENDIF
          SURFX(NFXX)=LFOUT2
          NFXX=NFXX+1
        ENDIF
      ENDIF
*
 478  SURFX(NFXX)=IFOUT2
      NFXX=NFXX+1
      NFAC=NFAC+1
      IF(NFAC.LE.NFACX) THEN
        IFACE=FACES(NFAC-1)
        GOTO 37
      ENDIF
*
      IF(DIRUP) THEN
        IFOUT1=IFOUT2
        IVOUT1=IVOUT
        DIRUP=.FALSE.
        ICELL=ICDOWN
        ICELL0=ICELL
        KPL=INT(AINT(REAL(ICELL)/REAL(STAIRS(1))))
        IF(MOD(ICELL,STAIRS(1)).NE.0) KPL=KPL+1
        LGDIR=.TRUE.
        LPOP=IPOP
        NSURF1=NFXX-1
        DO 377 IXX=1,NSURF1
          SURF1(IXX)=SURFX(IXX)
          IF(SURF1(IXX).GT.NSOUT) CALL XABORT('TRKHEX: SURF1 OVERFLOW.')
 377    CONTINUE
        GOTO 30
      ELSE
        NSURF2=NFAC-1
        DO 378 IXX=1,NSURF2
          SURF2(IXX)=SURFX(IXX)
          IF(SURF2(IXX).GT.NSOUT) CALL XABORT('TRKHEX: SURF2 OVERFLOW.')
 378    CONTINUE
      ENDIF
*
      IF(LGSTOR) THEN
        IF(LGDIR) THEN
          KOUT=IVOUT1
          IVOUT1=IVOUT
          IVOUT=KOUT
          LGDIR=.FALSE.
        ENDIF
        GOTO 38
      ENDIF
*
*       STOCKAGE DES PARCOURS OPTIQUES
*
   38 CONTINUE
      WEIGHT=ABS(COS1)*POIDS
      AUX1=1./DBLE(NSURF)
      ISXX=0
      DO 439 ISS=NSURF1+1,NSURF
         ISXX=ISXX+1
         IF(ISXX.GT.NSURF1)ISXX=1
         SURF1(ISS)=SURF1(ISXX)
 439  CONTINUE
      ISXX=0
      DO 442 ISS=NSURF2+1,NSURF
        ISXX=ISXX+1
        IF(ISXX.GT.NSURF2)ISXX=1
        SURF2(ISS)=SURF2(ISXX)
 442  CONTINUE
      IF(IPRT .GE. 500) THEN
        WRITE(IOUT,6000) ILINE,WEIGHT,IPOP,Y0,Z0,-SURF1(NSURF),
     >                   -SURF2(NSURF)
        WRITE(IOUT,6001) (MAT(K),POP(K),K=LPOP,1,-1),
     >                   (MAT(K),POP(K),K=LPOP+1,IPOP)
      ENDIF
      WRITE(IFILE) 1,IPOP+NSURF+NSURF,WEIGHT,IANGL,
     +             (-SURF1(KS),KS=1,NSURF),(MAT(K),K=LPOP,1,-1),
     +             (MAT(K),K=LPOP+1,IPOP),(-SURF2(KS),KS=1,NSURF),
     +             (AUX1,KS=1,NSURF),(POP(K),K=LPOP,1,-1),
     +             (POP(K),K=LPOP+1,IPOP),(AUX1,KS=1,NSURF)
*
*  INITIALIZATION
*
      DO 383 IS=1,5
         FACEM(IS)=0
         SURF1(IS)=0
         SURF2(IS)=0
  383 CONTINUE
*
  384 CONTINUE
      Y0=Y00
      Z0=Z00
*
*       POSSIBILITE DE LA CONTINUITE DE LA TRACK SUR D'AUTRES VOLUMES
*       APRES AVOIR INTERSECTEE UN VOLUME PERIPHERIQUE
*
      LGPER=.FALSE.
      IF(NCEL.GT.1) THEN
      LGPER=.TRUE.
*
*-- RECHERCHE DES CELLULES DANS CORN OU LA TRACK PEUT ABOUTIR
*
      IX=1
      I=CORN(IX)
      LGDEB=.TRUE.
      GO TO 401
      ENDIF
*
*       DEPLACEMENT HORIZENTAL PUIS VERTICAL DANS LE PLAN YOZ
*
   40 Y0=Y0+PASY

      ILINE=ILINE+1
      LGOUT1=.FALSE.
      LGOUTD=.FALSE.
      LGOUTU=.FALSE.
      LGDEB=.TRUE.
      LGPAS0=.FALSE.
      IF(Y0.GT.Y0MAX) GOTO 53
      IF(LGDIM) THEN
  711   SSQ=R2AZ-SCOS1*Y0*Y0
        IF(SSQ.GT.0.) THEN
          SSS=SQRT(SSQ)
          Y0COS=Y0*COS2
          XZTMAX=Z0+(SSS-Y0COS)*COAZ
          XZTMIN=Z0+(-SSS-Y0COS)*COAZ
          ZTMIN=MIN(XZTMAX,XZTMIN)
          ZTMAX=MAX(XZTMAX,XZTMIN)
          IF(ZTMAX.LT.ZMIN) THEN
            Y0=Y0+PASY
            ILINE=ILINE+1
            IF(Y0.GT.Y0MAX) GOTO 53
            GOTO 711
          ELSEIF(ZTMIN.GT.ZMAX) THEN
            Y0=Y0+PASY
            ILINE=ILINE+1
            IF(Y0.GT.Y0MAX) GOTO 53
            GOTO 711
          ENDIF
        ENDIF
      ENDIF
*CC   I=1
      IX=1
      I=CORN(1)
      KPER=0
  401 Y00=Y0
      IPOP=0
      START=.TRUE.
      DIRUP= .TRUE.
      LGPAS0=.FALSE.
*
*        RECHERCHE DE LA PREMIERE CELLULE INTERSECTEE
*
   42 CONTINUE
   50 CONTINUE
         IF(LGPER) THEN
           DO 421 IP=1,KPER
           IF(I.EQ.MAT2(IP))GOTO 56
  421      CONTINUE
         ENDIF
         X=REMESH(I)
         Y=REMESH(NCEL+I)
         IF(LGDIM) THEN
           Z2=REMESH(NCELZ+I)
           Z1=0.
           IF(LGPG1) THEN
           IF(I.GT.STAIRS(1))THEN
             KPL=INT(AINT(REAL(I)/REAL(STAIRS(1))))
             IF(MOD(I,STAIRS(1)).NE.0) KPL=KPL+1
             ICC=I-STAIRS(KPL-1)
             IF(KPL.GT.2)ICC=ICC+STAIRS(KPL-2)
             Z1=REMESH(NCELZ+ICC)
           ENDIF
           ENDIF
         ENDIF
         YY=Y-Y0
         TERM1=SQRT3*(X+A)
         TERM2=SQRT3*(X-A)
         T(1)=(TERM1+YY)*DIV1
         T(2)=(YY+ACOS6)*COS2I
         T(3)=(YY-TERM2)*DIV2
         T(4)=(YY+TERM2)*DIV1
         T(5)=(YY-ACOS6)*COS2I
         T(6)=(YY-TERM1)*DIV2
         XDR(1)=COS1*T(1)
         YDR(1)=SQRT3*(X-XDR(1)+A)+Y
         XDR(2)=COS1*T(2)
         YDR(2)=Y+ACOS6
         XDR(3)=COS1*T(3)
         YDR(3)=SQRT3*(-X+A+XDR(3))+Y
         XDR(4)=COS1*T(4)
         YDR(4)=SQRT3*(X-XDR(4)-A)+Y
         XDR(5)=COS1*T(5)
         YDR(5)=Y-ACOS6
         XDR(6)=COS1*T(6)
         YDR(6)=SQRT3*(-X+XDR(6)-A)+Y
         IF(LGDIM) THEN
           T(7)=(Z1-Z0)/COS3
           T(8)=(Z2-Z0)/COS3
           ZDR(1)=COS3*T(1)+Z0
           ZDR(2)=COS3*T(2)+Z0
           ZDR(3)=COS3*T(3)+Z0
           ZDR(4)=COS3*T(4)+Z0
           ZDR(5)=COS3*T(5)+Z0
           ZDR(6)=COS3*T(6)+Z0
           XDR(7)=COS1*T(7)
           YDR(7)=COS2*T(7)+Y0
           ZDR(7)=Z1
           XDR(8)=COS1*T(8)
           YDR(8)=COS2*T(8)+Y0
           ZDR(8)=Z2
         ENDIF
         DO 45 J=1,LFACE
            YDROIT=YDR(J)
            XDROIT=XDR(J)
            IF(LGDIM) THEN
              ZDROIT=ZDR(J)
              IF(ZDROIT.LE.Z2.AND.ZDROIT.GE.Z1) GOTO 64
              GOTO 45
            ENDIF
   64       CONTINUE
*----
*  MODIFIED FOR PRECISION ON LINUX
*----
            YT1=YDROIT-YDR(5)
            IF(ABS(YT1) .LT. EPST) YT1=0.0D0
            YT2=YDR(2)-YDROIT
            IF(ABS(YT2) .LT. EPST) YT2=0.0D0
            YT3=YDROIT-SQRT3*(X-XDROIT-A)-Y
            IF(ABS(YT3) .LT. EPST) YT3=0.0D0 
            YT4=SQRT3*(-X+XDROIT+A)+Y-YDROIT
            IF(ABS(YT4) .LT. EPST) YT4=0.0D0
            YT5=YDROIT+SQRT3*(X-XDROIT+A)-Y
            IF(ABS(YT5) .LT. EPST) YT5=0.0D0
            YT6=SQRT3*(X-XDROIT+A)+Y-YDROIT
            IF(ABS(YT6) .LT. EPST) YT6=0.0D0
            XT1=XDROIT-X+A
            IF(ABS(XT1) .LT. EPST) XT1=0.0D0
            XT2=X+A-XDROIT
            IF(ABS(XT2) .LT. EPST) XT2=0.0D0
            XT3=X-A*0.5D0-XDROIT
            IF(ABS(XT3) .LT. EPST) XT3=0.0D0
            XT4=XDROIT-X-A*0.5D0
            IF(ABS(XT4) .LT. EPST) XT4=0.0D0
              IF(YT1 .GE. DZ0 .AND. YT2 .GE. DZ0) THEN
                IF(XT1 .GE. DZ0 .AND. XT2 .GE. DZ0) THEN
                  IF(XT3 .GT. DZ0)THEN
                    IF(YT4 .GE. DZ0 .AND. YT3 .GE. DZ0)THEN
                      ICELL=I
                      GOTO 15
                    ENDIF
                    GOTO 45
                  ENDIF
                  IF(XT4 .GT. DZ0)THEN
                    IF(YT6 .GE. DZ0 .AND. YT5 .GE. DZ0)THEN
                      ICELL=I
                      GOTO 15
                    ENDIF
                    GOTO 45
                  ENDIF
                  ICELL=I
                  GOTO 15
                ENDIF
              ENDIF
   45    CONTINUE
   56 CONTINUE
      IF(.NOT.LGPAS0) THEN
        IX=IX+1
        IF(IX.LE.ICOR) THEN
          I=CORN(IX)
          GOTO 50
        ENDIF
      ELSE
 57     ISURC=ISURC+1
        IF(ISURC.LE.ISURM) THEN
          I=VOISIN(SURC(ISURC),IVOLC)
          IF(I.GT.NCEL2) GOTO 57
          I=I+IETAG
          ICELL0=IVOLC+IETAG
          GOTO 50
        ENDIF
        IF(LGDIM) THEN
        IF(DIRUP) THEN
          LGPASU=.TRUE.
          IF(IETAG.LT.STAIRS(IPLANZ)) THEN
            IETAG=IETAG+STAIRS(1)
            ISURC=0
            IF(LGPG1) THEN
             I=IVOLC+IETAG
             IF(I.LE.NCEL)GOTO 50
            ENDIF
            GOTO 57
          ENDIF
          IF(LGPASU) THEN
             GOTO 384
          ENDIF
        ELSE
          LGPASD=.TRUE.
          IF(IETAG.GT.STAIRS(1)) THEN
            IETAG=IETAG-STAIRS(1)
            ISURC=0
            IF(LGPG1) THEN
             I=IVOLC+IETAG
             IF(I.LE.NCEL)GOTO 50
            ENDIF
            GOTO 57
          ENDIF
          IF(LGPASD) THEN
             GOTO 384
          ENDIF
        ENDIF
        ENDIF
      ENDIF
      IF(LGPER) THEN
        LGPER=.FALSE.
        GOTO 40
      ENDIF
*
   58 LGPAS0=.FALSE.
      Y0=Y0+PASY 
      ILINE=ILINE+1
      IF(Y0.GT.Y0MAX) GOTO 53
      IF(LGDIM) THEN
  811 SSQ=R2AZ-SCOS1*Y0*Y0
      IF(SSQ.GT.0.) THEN
      SSS=SQRT(SSQ)
      Y0COS=Y0*COS2
      XZTMAX=Z0+(SSS-Y0COS)*COAZ
      XZTMIN=Z0+(-SSS-Y0COS)*COAZ
      ZTMAX=MAX(XZTMAX,XZTMIN)
      ZTMIN=MIN(XZTMAX,XZTMIN)
      IF(ZTMAX.LT.ZMIN) THEN
        Y0=Y0+PASY
        ILINE=ILINE+1
        IF(Y0.GT.Y0MAX) GOTO 53
        GOTO 811
      ELSEIF(ZTMIN.GT.ZMAX) THEN
        Y0=Y0+PASY
        ILINE=ILINE+1
        IF(Y0.GT.Y0MAX) GOTO 53
        GOTO 811
      ENDIF
      ENDIF
      ENDIF
      Y00=Y0
      IX=1
      I=CORN(1)
      KPER=0
      LGDEB=.TRUE.
      GOTO 42
   53 IF(LGDIM) THEN
          KPER=0
          IPOP=0
          START=.TRUE.
          DIRUP= .TRUE.
        LGDEB=.TRUE.
        LGPAS0=.FALSE.
        Z0=Z0+PASZ
        Z00=Z0
        IF(Z0.LT.Z0MAX) GOTO 777
      ENDIF
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT(' INTEGRATION LINE ',I10,5X,'WEIGHT =',1P,E15.6,5X,
     >       'NUMBER OF SEGMENTS =',I10/
     >       ' Y0  = ',E15.6,10X,'Z0 = ',E15.6,5X,
     >       ' SURFI=',I10,5X,'SURFF=',I10)
 6001 FORMAT(1P,(I6,E15.6))
      END
