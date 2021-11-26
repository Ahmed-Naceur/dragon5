*DECK MESHST
      SUBROUTINE MESHST(IPTRK,IPGEOM,REMESH,FVOL,STAIRS,FACST,NCEL,
     >      IPLANZ,ISTATE,NCYL,NSECT,NCPHY,VOLSUR,MATALB,SIDE,NCOUR,
     >      NSMIN,NSMAX,NS,FACB,NVOL,SURB,VSYM,SSYM,IHEX,LXI,NV,MCODE,
     >      SURL,IPLANI,VLAT,ZMIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* General numbering of hexagonal assembly.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): M. Ouisloumen
*
*Parameters: input
* MCODE   =1 for Z- symmetry; =2 for Z+ symmetry; =0 otherwise.
*
*Parameters: input/output
* REMESH  coordinates of geometry. 
* FVOL    first zone number.
* NVOL    first volume number.
* STAIRS  see TRKHEX.
* FACST   see TRKHEX.
* NCEL    number of cells.
* NSECT   number of sectors in region.
* NCYL    number of cylinders in cell.
* NCPHY   number of physical cells.
* VOLSUR  volumes and external surfaces
* MATALB  material albedo vector.
* FACB    first face number.
* SURB    first cell number.
* VSYM    initial geometry volumes.
* SSYM    initial geometry surfaces.
* IPTRK   undefined.
* IPGEOM  undefined.
* IPLANZ  undefined.
* ISTATE  undefined.
* SIDE    undefined.
* NCOUR   undefined.
* NSMIN   undefined.
* NSMAX   undefined.
* NS      undefined.
* IHEX    undefined.
* LXI     undefined.
* NV      undefined.
* SURL    undefined.
* IPLANI  undefined.
* VLAT    undefined.
* ZMIN    undefined.
*
*-----------------------------------------------------------------------
*
      USE         GANLIB
      IMPLICIT    NONE
      INTEGER     NSTATE
      REAL        PI,SQRT3
      PARAMETER  (NSTATE=40)
      PARAMETER  (PI=3.141592653589793,SQRT3=1.732050807568877)
      TYPE(C_PTR) IPTRK,IPGEOM
      INTEGER     NCEL,IPLANZ,NCPHY,NCOUR,IHEX,LXI,NV,MCODE,IPLANI
      INTEGER     FVOL(NCEL),STAIRS(IPLANZ),FACST(IPLANZ),
     >            ISTATE(NSTATE),MATALB(*),NVOL(*),NSMIN,NSMAX,NS,
     >            NSECT(*),NCYL(NCEL),FACB(*),SURB(*),
     >            VSYM(*),SSYM(*),SURL(*),VLAT(*)
      REAL        REMESH(*),VOLSUR(*),SIDE,ZMIN
*----
* LOCAL VARIABLES
*----
      CHARACTER   TEXT12*12
      LOGICAL     LGCELL,L1CELL,LGSS,LGTURN,LG3D,LGPASS,LSPLIT,LGMERG
      REAL        POIDSH
      INTEGER     IFONC,IFCOUR,N,L
      INTEGER     NEIGHB
      INTEGER     LLL,KSS,LX,LZ,ISTAT9,NDIM,LZZ,ILENT,ITP,ICEL2,
     >            I,ICEL0,NCOUS,ICOU,KCEL,LCEL,ILENZ,IAZ,K,IAUX,
     >            ISURFP,ISURF,ISAUX,JSUR,ISF,IPOINT,KREG,KSURF1,
     >            KSURF2,NSURF1,NSURF2,ISUR,ITC,ILENR,MDR,IOFS,
     >            NRAY,ILENSP,MADD,J,ISEC,LSECT,MIXF,
     >            IFVOL,NZONE,MVOL,JVOL,IVOL,ISURB,JSURB,IST,IY,KZAUX,
     >            NAUX,MIXX,JSP,KBB,JJ,IV,IS,ILENS,ITS,IC,KMIX,NMIX,
     >            IXZPL,KXI,NXI,LFF,LZF,ISTO,ISS1,IMERMX,ISZ,
     >            IFRT,IMER0,IMER1,ISX,IMERG,JMERG,MMCYL,IZON,IFF,
     >            ILENM,ISS,ICEL
      INTEGER     JSECT,LL,IZZ,KY,IMMS,IBSS,IFR0,IXX,ICCOIN,KCCOIN,
     >            MROT,NCPER,ISURSY,MSMAX,ISSXX,IPPZ,IPP1,IXP,JCCOIN,
     >            JP,MXI,IP,ISURB6,ITT9,ISURSX,KSURBX,K1,ICX,ICY,ICZ,
     >            IDEBX,NBASE,ICELC0,ICELC1,NCC,MCYL,NZONE1,IA,
     >            NSECTO,IV1,JVT,IVTURN,ITURN,KAUX,KVOL1,KVOL2,
     >            KSECTX,NZZ,ICXX,ICELC,ICLIM,IXX0,ISS0,JSS0,M,ISS2,
     >            IJSUR,KKB,IXV,LFROT,IVFIN,IVORIN,IVMIN,IVSYM,NZON,
     >            IVAUX,IVMAX,IVOR,IVV,IVSYM0,IVV1,KVV,KVOR,ISY,MSAUX,
     >            IVLMAX,NSAUX,IPPX,ISYAUX,LSMAX,LSPLZM,IDDX,KX2,KSECT
      INTEGER     KX1,KXP,IDEB,IYAUX,II,KSAUX,IVLAT,JTX,IX1,JX,IFR,IW,
     >            JW,MSUR,MCPHY,ISYX,LSPLZP
      REAL        VEX,X,Y,XM,YM,XP,YP,YPP,ZAUX,ZBUX,Z,SURF,R1,R2,PAS,
     >            VOLUM1,VOLUM2,R,VOLUME,VOLUMS,VCYL,SAUX,XTAN,SAUX1
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: VTURN,PHTURN,VOL1,ISECT,MXX,
     + IMX,IMX2,ISSS,IMIX,IPSECT,ICC,ICELL,KNUM,NUMG,ITRN,ISPZ,KVOL,IBB,
     + ITT,IAA
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISXY
      REAL, ALLOCATABLE, DIMENSION(:) :: ZZZ,RAY,KRAY,RR
      TYPE PP
        REAL, POINTER, DIMENSION(:) :: R
        INTEGER, POINTER, DIMENSION(:) :: I1,I2
      END TYPE PP
      TYPE(PP), ALLOCATABLE, DIMENSION(:) :: IOF
*----
*  DATA
*----
      INTEGER     ROT(12)
      SAVE        ROT
      DATA        ROT /2,3,4,5,6,1,6,1,2,3,4,5/
*
*  STATEMENT FUNCTIONS
*----
      IFONC(N,L)= 2+(N-1)*(L+3*(N-2))
      IFCOUR(N)=NINT( (4.+SQRT(1.+4.*FLOAT(N-1)/3.)
     +                 +SQRT(1.+4.*FLOAT(N-2)/3.))*.25)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(VOL1(NCPHY),VTURN(NCPHY),PHTURN(NCPHY))
*
      LZ=ISTATE(5)
      LX=ISTATE(3)
      ISTAT9=ISTATE(9)
      LGCELL=.FALSE.
      L1CELL=.FALSE.
      LGSS=.FALSE.
      LGTURN=.FALSE.
      LG3D=.FALSE.
      IF(ISTATE(8).EQ.1) LGCELL=.TRUE.
      NDIM=2
      LZZ=1
      L1CELL=(LX.EQ.1)
      IF(LZ.GT.0) THEN
        NDIM=3
        LG3D=.TRUE.
        LZZ=LZ
        L1CELL=(LX*LZ.EQ.1)
      ENDIF
      VEX=1.5*SQRT(3.)*SIDE*SIDE
      CALL LCMLEN(IPGEOM,'TURN',ILENT,ITP)
      IF(ILENT.GT.0) LGTURN=.TRUE.
      NCOUR=1
      IF(LX.GT.1)NCOUR=IFCOUR(LX)
*
*              COORDONNEES DES CENTRES DES HEXAGONES
*
      X=0.
      Y=0.
      XM=0.
      YM=0.
      YP=0.
      XP=0.
      ICEL=1
   18 REMESH(ICEL)=X
      REMESH(NCEL+ICEL)=Y
      ICEL2=NEIGHB(ICEL,2,9,LX,POIDSH)
      IF(ICEL2.GT.LX) GOTO 20
      ICEL=ICEL2
      Y=Y+SIDE*SQRT3
      GO TO 18
   20 ICEL=1
      Y=0
   21 ICEL2=NEIGHB(ICEL,5,9,LX,POIDSH)
      IF(ICEL2.GT.LX) GOTO 22
      ICEL=ICEL2
      REMESH(ICEL)=X
      Y=Y-SIDE*SQRT3
      REMESH(NCEL+ICEL)=Y
      GOTO 21
   22 CONTINUE
      DO 30 I=2,NCOUR
         LGPASS=.FALSE.
         XP=XP+1.5*SIDE
         YP=YP+.5*SIDE*SQRT3
         ICEL=IFONC(I,0)
         ICEL0=ICEL
         X=XP
         Y=YP
         YPP=YP
   23    REMESH(ICEL)=X
         REMESH(NCEL+ICEL)=Y
         ICEL2=NEIGHB(ICEL,2,9,LX,POIDSH)
         IF(ICEL2.GT.LX) GOTO 25
         Y=Y+SIDE*SQRT3
         ICEL=ICEL2
         GO TO 23
   25    ICEL=ICEL0
         Y=YP
   26    ICEL2=NEIGHB(ICEL,5,9,LX,POIDSH)
         IF(ICEL2.GT.LX) GOTO 27
         Y=Y-SIDE*SQRT3
         ICEL=ICEL2
         REMESH(ICEL)=X
         REMESH(NCEL+ICEL)=Y
         GOTO 26
   27    IF(LGPASS) GOTO 28
         XM=XM-1.5*SIDE
         YM=YM+.5*SIDE*SQRT3
         YPP=YM
         ICEL=IFONC(I,2)
         ICEL0=ICEL
         X=XM
         Y=YM
         LGPASS=.TRUE.
         GO TO 23
   28    CONTINUE
   30 CONTINUE
*
*         POUR EVITER DES EFFETS DE DIFFERENCE
*
      NCOUS=7
      DO 31 ICOU=3,NCOUR,2
         KCEL=IFONC(ICOU,0)+6*ICOU-NCOUS
         LCEL=KCEL-3*(ICOU-1)
         REMESH(NCEL+KCEL)=0.
         REMESH(NCEL+LCEL)=0.
         NCOUS=NCOUS+1
   31 CONTINUE
      STAIRS(1)=LX
      IF(LG3D) THEN
        CALL LCMSIX(IPTRK,'DATA_DUP',1)
        CALL LCMLEN(IPTRK,'MESHZ',ILENZ,ITP)
        IF(ILENZ.NE.LZZ+1)CALL XABORT('MESHST: MISSING DIMENSION MESHZ')
        ALLOCATE(ZZZ(ILENZ))
        CALL LCMGET(IPTRK,'MESHZ',ZZZ)
        CALL LCMSIX(IPTRK,' ',2)
*
*--- TEST SUR L'ORDONANCE DES PLANS DES CELLULES SUIVANT L'AXE Z
*
        ZAUX=ZZZ(1)
        ZMIN=ZAUX
        DO 29 IAZ=1,ILENZ-1
          ZBUX=ZZZ(IAZ+1)
        IF(ZAUX.GE.ZBUX)CALL XABORT('MESHST: WRONG CELLS PLAN POSITION')
          ZAUX=ZBUX
   29   CONTINUE
        Z=ZZZ(2)
        DO 32 K=1,LX
          REMESH(2*NCEL+K)=Z
   32   CONTINUE
        IAUX=0
        DO 35 I=2,LZ
          STAIRS(I)=I*LX
          Z=ZZZ(I+1)
          DO 33 K=1,LX
            IAUX=IAUX+1
            REMESH(LX+IAUX)=REMESH(K)
            REMESH(NCEL+LX+IAUX)=REMESH(NCEL+K)
            REMESH(2*NCEL+LX+IAUX)=Z
   33     CONTINUE
   35   CONTINUE
      ENDIF
*
*           CALCUL DES SURFACES EXTERNES ET AFFECTATION DES ALBEDOS
*
      ISURFP=6
      IF(NCOUR.GT.1)ISURFP=6*(3+2*(NCOUR-2))
      ISS=1
      ISURF=NS
      ISAUX=ISURF
      IF(LG3D) THEN
        ISAUX=ISURF-NSMAX
        ISS=NSMIN+1
        DO 70 I=1,NSMIN
           MATALB(I)=-6
   70   CONTINUE
        DO 71 I=1,NSMAX
           MATALB(ISAUX+I)=-5
   71   CONTINUE
      ENDIF
      DO 75 I=ISS,ISAUX
         MATALB(I)=-1
   75 CONTINUE
      JSUR=NS-NSMAX+1
      DO 80 K=1,LZZ
         ISF=ISS+ISURFP-1
         SURF=SIDE
         FACST(K)=ISF
         ISS=ISF+1
   80 CONTINUE
      IPOINT=2*NCEL
      IF(LG3D) THEN
         IPOINT=3*NCEL
      ENDIF
*
      KREG=ISURF+1
      VOLSUR(KREG)=0.
      MATALB(KREG)=0
      KSURF1=0
      KSURF2=ISAUX
      NSURF1=NSMIN+1
      NSURF2=NS+1
      ISUR=0
      Z=0.0
      KSECT=0
      IF(L1CELL) THEN
*
*            CAS D'UNE SEULE CELLULE
*
        CALL LCMLEN(IPGEOM,'SECTOR',ILENS,ITC)
        CALL LCMLEN(IPGEOM,'RADIUS',ILENR,ITP)
        MDR=ILENR
        IF(ILENR.EQ.0) MDR=1
        ALLOCATE(ISECT(MDR),MXX(MDR))
        CALL LCMLEN(IPGEOM,'MIX',ILENM,ITP)
        IF(ILENS.GT.0) THEN
          IF(ILENS.NE.MDR) CALL XABORT('MESHST: INCONSISTENT LENGTHS(1'
     +    //').')
          CALL LCMGET(IPGEOM,'SECTOR',ISECT)
          K=0
          DO J=1,MDR
            K=K+6*(ISECT(J)-1)
          ENDDO
          IF(K.NE.ILENM) THEN
            CALL XABORT('MESHST: SECTOR-DEFINED MIX INDICES EXPECTED(1'
     +      //').')
          ENDIF
          ALLOCATE(IMX2(ILENM))
          CALL LCMGET(IPGEOM,'MIX',IMX2)
          IOFS=0
          DO J=1,MDR
            IOFS=IOFS+1
            MXX(J)=IMX2(IOFS)
            DO K=2,6*(ISECT(J)-1)
              IOFS=IOFS+1
              IF(IMX2(IOFS).NE.MXX(J)) THEN
                CALL XABORT('MESHST: SECTOR-DEPENDENT MIX NOT IM'
     +          //'PLEMENTED(1).')
              ENDIF
            ENDDO
          ENDDO
          DEALLOCATE(IMX2)
        ELSE
          IF(ILENM.NE.MDR) CALL XABORT('MESHST: INCONSISTENT LENGTHS(2'
     +    //').')
          CALL LCMGET(IPGEOM,'MIX',MXX)
          CALL XDISET(ISECT,MDR,1)
        ENDIF
        IF(ILENR.GT.0) THEN
          ALLOCATE(RAY(ILENR))
          CALL LCMGET(IPGEOM,'RADIUS',RAY)
          CALL LCMLEN(IPGEOM,'SPLITR',ILENSP,ITP)
          IF(ILENSP.GT.0) THEN
            IF(ITP.NE.1.OR.ILENSP.NE.ILENR-1)
     +         CALL XABORT('MESHST: '//
     +                     'MISSING TYPE OR DIMENSION OF SPLITR')
            ALLOCATE(ISSS(ILENSP))
            CALL LCMGET(IPGEOM,'SPLITR',ISSS)
            NRAY=0
            DO 101 I=1,ILENSP
               NRAY=NRAY+ABS(ISSS(I))
  101       CONTINUE
            ALLOCATE(KRAY(NRAY+1),IMIX(NRAY+1),IPSECT(NRAY+1))
            R1=0.0
            MADD=-1
            KRAY(1)=0.0
            DO 103 J=1,ILENSP
               ISEC=ISSS(J)
               LSECT=ISECT(J)
               MIXF=MXX(J)
               IF(ISEC.EQ.0)CALL XABORT('MESHST: INVALID SPLITR')
               IF(ISEC.GT.0) THEN
                 R2=RAY(J+1)
                 PAS=(R2-R1)/FLOAT(ISEC)
                 DO 111 K=1,ISEC
                   MADD=MADD+1
                   KRAY(MADD+2)=R1+PAS*FLOAT(K)
                   IMIX(MADD+1)=MIXF
                   IPSECT(MADD+1)=LSECT
  111            CONTINUE
                 R1=R2
               ELSE
                 R2=RAY(J+1)
                 R1=R1**2
                 R2=R2**2
                 PAS=(R2-R1)/FLOAT(-ISEC)
                 DO 112 K=1,-ISEC
                   MADD=MADD+1
                   KRAY(MADD+2)=SQRT(R1+PAS*FLOAT(K))
                   IMIX(MADD+1)=MIXF
                   IPSECT(MADD+1)=LSECT
  112            CONTINUE
                 R1=SQRT(R2)
               ENDIF
  103       CONTINUE
            IMIX(NRAY+1)=MXX(ILENR)
            IPSECT(NRAY+1)=ISECT(ILENR)
            DEALLOCATE(MXX,ISSS,ISECT)
            ALLOCATE(ISECT(NRAY+1))
            DO 104 J=1,NRAY+1
              ISECT(J)=IPSECT(J)
  104       CONTINUE
            DEALLOCATE(IPSECT)
          ELSE
            NRAY=ILENR-1
            ALLOCATE(IMIX(NRAY+1),KRAY(NRAY+1))
            DO 105 J=1,NRAY+1
              IMIX(J)=MXX(J)
              KRAY(J)=RAY(J)
  105       CONTINUE
          ENDIF
          DEALLOCATE(MXX,RAY)
        ELSE
          ALLOCATE(IMIX(MDR))
          DO 106 J=1,MDR
            IMIX(J)=MXX(J)
  106     CONTINUE
          DEALLOCATE(MXX)
          NRAY=0
        ENDIF
        IAUX=2
        IF(LG3D)IAUX=3*NCEL
        IFVOL=0
        NZONE=NRAY+1
        MVOL=0
        JVOL=1
        IVOL=0
        ISURB=NSMIN
        JSURB=0
        DO 119 I=1,NCEL
           NVOL(I)=JVOL
           IF(LG3D) THEN
             Z=REMESH(2*NCEL+I)
             IF(I.GT.STAIRS(1)) THEN
             IST=1
             DO 129 IY=2,LZZ
             IF(I.LE.STAIRS(IY)) THEN
               IST=IY-1
               GOTO 139
             ENDIF
  129        CONTINUE
  139        KZAUX=I-STAIRS(IST)
             IF(IST.GT.1)KZAUX=KZAUX+STAIRS(IST-1)
             Z=Z-REMESH(2*NCEL+KZAUX)
             ENDIF
           ENDIF
           NCYL(I)=NRAY
           FVOL(I)=IFVOL+1
           IFVOL=IFVOL+NZONE
           DO 117 J=1,NRAY
              IAUX=IAUX+1
              REMESH(IAUX)=KRAY(J+1)
  117      CONTINUE
           VOLUM1=0.
           IF(LG3D) THEN
             IF(IPLANZ.GT.1) THEN
               IF(I.LE.STAIRS(1)) THEN
                 FACB(I)=KSURF1
               ELSEIF(I.GT.STAIRS(IPLANZ-1)) THEN
                 ISUR=ISUR+1
                 FACB(LX+ISUR)=KSURF2
               ENDIF
             ELSE
               FACB(I)=KSURF1
               FACB(LX+I)=KSURF2
             ENDIF
           ENDIF
           NAUX=0
           DO 120 J=1,NRAY+1
             MVOL=MVOL+1
             KSECT=ISECT(J)
             NSECT(MVOL)=KSECT
             NAUX=1
             IF(KSECT.GT.1) NAUX=6*(KSECT-1)
             IVOL=IVOL+NAUX
             MIXX=IMIX(J)
             KREG=KREG+1
             IF(J.EQ.NZONE) THEN
                VOLUM2=VEX
             ELSE
                R=KRAY(J+1)
                VOLUM2=PI*R*R
             ENDIF
             VOLUME=(VOLUM2-VOLUM1)/REAL(NAUX)
             VOLUMS=0.0
             IF(LG3D) THEN
                VOLUMS=VOLUME*.25
*
*     SURFACES SUPERIEURES ET INFERIEURES (SELON L'AXE Z)
*
                IF(IPLANZ.GT.1) THEN
                  IF(I.LE.STAIRS(1)) THEN
                    KSURF1=KSURF1+1
                    NSURF1=NSURF1-1
                    VOLSUR(NSURF1)=VOLUMS
                  ELSEIF(I.GT.STAIRS(IPLANZ-1)) THEN
                    KSURF2=KSURF2+1
                    NSURF2=NSURF2-1
                    VOLSUR(NSURF2)=VOLUMS
                  ENDIF
                ELSE
                  KSURF1=KSURF1+1
                  NSURF1=NSURF1-1
                  VOLSUR(NSURF1)=VOLUMS
                  KSURF2=KSURF2+1
                  NSURF2=NSURF2-1
                  VOLSUR(NSURF2)=VOLUMS
                ENDIF
                VOLUME=VOLUME*Z
             ENDIF
             VOLUM1=VOLUM2
             VOLSUR(KREG)=VOLUME
             MATALB(KREG)=MIXX
             DO 121 K=2,NAUX
                KREG=KREG+1
                VOLSUR(KREG)=VOLUME
                MATALB(KREG)=MIXX
                IF(LG3D) THEN
                  IF(IPLANZ.GT.1) THEN
                    IF(I.LE.STAIRS(1)) THEN
                      KSURF1=KSURF1+1
                      NSURF1=NSURF1-1
                      VOLSUR(NSURF1)=VOLUMS
                    ELSEIF(I.GT.STAIRS(IPLANZ-1)) THEN
                      KSURF2=KSURF2+1
                      NSURF2=NSURF2-1
                      VOLSUR(NSURF2)=VOLUMS
                    ENDIF
                  ELSE
                    KSURF1=KSURF1+1
                    NSURF1=NSURF1-1
                    VOLSUR(NSURF1)=VOLUMS
                    KSURF2=KSURF2+1
                    NSURF2=NSURF2-1
                    VOLSUR(NSURF2)=VOLUMS
                  ENDIF
                ENDIF
  121        CONTINUE
             JVOL=IVOL+1
  120      CONTINUE
           SURF=SIDE
           JSP=6
           KBB=1
           IF(KSECT.GT.1) THEN
             SURF=SIDE/(KSECT-1)
             JSP=NAUX
             KBB=KSECT-1
           ENDIF
           IF(LG3D)SURF=SURF*Z
           DO 118 JJ=1,JSP
              JSUR=JSUR-1
              VOLSUR(JSUR)=.25*SURF
  118      CONTINUE
           DO 122 JJ=1,6
              JSURB=JSURB+1
              SURB(JSURB)=ISURB
              ISURB=ISURB+KBB
  122      CONTINUE
  119   CONTINUE
        DEALLOCATE(IMIX,ISECT)
        IF(LG3D) DEALLOCATE(ZZZ)
        IF(NRAY.GT.0) DEALLOCATE(KRAY)
        DO 715 IV=1,NV
         VSYM(IV)=IV
  715   CONTINUE
        DO 716 IS=1,NS
         SSYM(IS)=IS
  716   CONTINUE
        GO TO 800
      ENDIF
*         CAS D'UN ASSEMBLAGE DE CELLULE
*
      IPP1=0
      MROT=0
      MSMAX=0
      ALLOCATE(ISXY(6,MAX(IPLANZ,IPLANI)))
      IF(LGCELL) THEN
*
*               SPLITING DES RAYONS ET DEFINITION DES SECTEURS POUR
*               LES CELLULES GENERATRICES
*
         ALLOCATE(IOF(ISTAT9),ICC(ISTAT9),ICELL(3*ISTAT9))
         CALL LCMGET(IPGEOM,'CELL',ICELL)
         DO 10 I=1,ISTAT9
            WRITE(TEXT12(1:4),'(A4)') ICELL(3*I-2)
            WRITE(TEXT12(5:8),'(A4)') ICELL(3*I-1)
            WRITE(TEXT12(9:12),'(A4)') ICELL(3*I)
            CALL LCMSIX(IPGEOM,TEXT12,1)
            CALL LCMLEN(IPGEOM,'RADIUS',ILENR,ITP)
            CALL LCMLEN(IPGEOM,'SECTOR',ILENS,ITS)
            IF(ILENR.GT.0) THEN
              ALLOCATE(ISECT(ILENR),IMX(ILENR))
              CALL LCMLEN(IPGEOM,'MIX',ILENM,ITP)
              IF(ILENS.GT.0) THEN
                IF(ILENS.NE.ILENR) CALL XABORT('MESHST: INCONSISTENT L'
     +          //'ENGTHS(3).')
                CALL LCMGET(IPGEOM,'SECTOR',ISECT)
                K=0
                DO J=1,ILENR
                  K=K+6*(ISECT(J)-1)
                ENDDO
                IF(K.NE.ILENM) THEN
                  CALL XABORT('MESHST: SECTOR-DEFINED MIX INDICES EXPE'
     +            //'CTED(2).')
                ENDIF
                ALLOCATE(IMX2(ILENM))
                CALL LCMGET(IPGEOM,'MIX',IMX2)
                IOFS=0
                DO J=1,ILENR
                  IOFS=IOFS+1
                  IMX(J)=IMX2(IOFS)
                  DO K=2,6*(ISECT(J)-1)
                    IOFS=IOFS+1
                    IF(IMX2(IOFS).NE.IMX(J)) THEN
                      CALL XABORT('MESHST: SECTOR-DEPENDENT MIX NOT IM'
     +                //'PLEMENTED(2).')
                    ENDIF
                  ENDDO
                ENDDO
                DEALLOCATE(IMX2)
              ELSE
                IF(ILENM.NE.ILENR) CALL XABORT('MESHST: INCONSISTENT L'
     +          //'ENGTHS(4).')
                CALL XDISET(ISECT,ILENR,1)
                CALL LCMGET(IPGEOM,'MIX',IMX)
              ENDIF
              CALL LCMLEN(IPGEOM,'SPLITR',ILENSP,ITP)
              IF(ILENSP.GT.0) THEN
                IF(ITP.NE.1.OR.ILENSP.NE.ILENR-1)
     +            CALL XABORT('MESHST: '//
     +                        'MISSING TYPE OR DIMENSION OF SPLITR')
                ALLOCATE(RR(ILENR))
                CALL LCMGET(IPGEOM,'RADIUS',RR)
                ALLOCATE(ISSS(ILENSP))
                CALL LCMGET(IPGEOM,'SPLITR',ISSS)
                NRAY=1
                DO 5 J=1,ILENSP
                 NRAY=NRAY+ABS(ISSS(J))
    5           CONTINUE
                ALLOCATE(IOF(I)%R(NRAY),IOF(I)%I1(NRAY),IOF(I)%I2(NRAY))
                KMIX=0
                DO 7 J=1,ILENR-1
                   NMIX=IMX(J)
                   ISEC=ISSS(J)
                   LSECT=ISECT(J)
                   IF(ISEC.EQ.0)CALL XABORT('MESHST: INVALID SPLITR')
                   DO 6 L=1,ABS(ISEC)
                      IOF(I)%I1(KMIX+1)=NMIX
                      IOF(I)%I2(KMIX+1)=LSECT
                      KMIX=KMIX+1
    6              CONTINUE
    7           CONTINUE
                IOF(I)%I1(KMIX+1)=IMX(ILENR)
                IOF(I)%I2(KMIX+1)=ISECT(ILENR)
                IOF(I)%R(1)=0.0
                R1=0.
                IAUX=0
                DO 109 K=1,ILENSP
                  ISEC=ISSS(K)
                  IF(ISEC.EQ.0)CALLXABORT('MESHST: INVALID SPLITR')
                  IF(ISEC.GT.0) THEN
                   R2=RR(K+1)
                   PAS=(R2-R1)/FLOAT(ISEC)
                   DO 8 L=1,ISEC
                    IAUX=IAUX+1
                    IOF(I)%R(IAUX+1)=R1+PAS*REAL(L)
    8              CONTINUE
                   R1=R2
                  ELSE
                   R2=RR(K+1)
                   R1=R1**2
                   R2=R2**2
                   PAS=(R2-R1)/FLOAT(-ISEC)
                   DO 108 L=1,-ISEC
                    IAUX=IAUX+1
                    IOF(I)%R(IAUX+1)=SQRT(R1+PAS*REAL(L))
  108              CONTINUE
                   R1=SQRT(R2)
                  ENDIF
  109           CONTINUE
                DEALLOCATE(ISSS,RR)
                ICC(I)=NRAY-1
              ELSE
                ICC(I)=ILENR-1
                ALLOCATE(IOF(I)%R(ILENR))
                CALL LCMGET(IPGEOM,'RADIUS',IOF(I)%R)
                ALLOCATE(IOF(I)%I1(ILENR),IOF(I)%I2(ILENR))
                DO 9 J=1,ILENR
                   IOF(I)%I1(J)=IMX(J)
                   IOF(I)%I2(J)=ISECT(J)
    9           CONTINUE
              ENDIF
              DEALLOCATE(IMX,ISECT)
             ELSE
              ICC(I)=0
              ALLOCATE(IOF(I)%I1(1),IOF(I)%I2(1))
              CALL LCMLEN(IPGEOM,'MIX',ILENM,ITP)
              IF(ILENS.GT.0) THEN
                CALL LCMGET(IPGEOM,'SECTOR',IOF(I)%I2(1))
                IF(6*(IOF(I)%I2(1)-1).NE.ILENM) THEN
                  CALL XABORT('MESHST: SECTOR-DEFINED MIX INDICES EXPE'
     +            //'CTED(3).')
                ENDIF
                ALLOCATE(IMX2(ILENM))
                CALL LCMGET(IPGEOM,'MIX',IMX2)
                DO K=2,6*(IOF(I)%I2(1)-1)
                  IF(IMX2(K).NE.IMX2(1)) THEN
                    CALL XABORT('MESHST: SECTOR-DEPENDENT MIX NOT IMPL'
     +              //'EMENTED(3).')
                  ENDIF
                ENDDO
                IOF(I)%I1(1)=IMX2(1)
                DEALLOCATE(IMX2)
              ELSE
                IF(ILENM.NE.1) CALL XABORT('MESHST: INCONSISTENT LENGT'
     +          //'HS(5).')
                CALL LCMGET(IPGEOM,'MIX',IOF(I)%I1(1))
                IOF(I)%I2(1)=1
              ENDIF
             ENDIF
             IF(ILENS.EQ.0) THEN
               IF(LGTURN)CALL XABORT('MESHST: SECTOR MUST BE DEFINED '//
     +                               'OR CANCEL TURN ')
             ENDIF
             CALL LCMSIX(IPGEOM,' ',2)
   10    CONTINUE
*--- NUMEROTATION DES VOLUMES DE LA SYMETRIE D'ENTREE
         IXZPL=IPLANZ
         IF(MCODE.GT.0)IXZPL=IPLANI
         KXI=LXI*IXZPL
         NXI=KXI
         IF(IHEX.LE.9) THEN
           CALL LCMSIX(IPTRK,'DATA_DUP',1)
           CALL LCMLEN(IPTRK,'SPLITZD',LFF,KSS)
           CALL LCMSIX(IPTRK,' ',2)
           LSPLIT=.FALSE.
           LZF=IXZPL
           IF(LFF.GT.0) THEN
             LZF=LFF
             NXI=LXI*LZF
             LSPLIT=.TRUE.
           ENDIF
           ALLOCATE(KNUM(2*NXI),NUMG(NXI))
           CALL LCMLEN(IPGEOM,'MERGE',LLL,KSS)
           LGMERG=.FALSE.
           IF(LLL.GT.0) THEN
             CALL LCMGET(IPGEOM,'MERGE',KNUM)
             LGMERG=.TRUE.
             CALL LCMGET(IPGEOM,'MIX',NUMG)
             DO I=1,NXI
               NUMG(I)=-NUMG(I)
             ENDDO
           ELSE
             CALL LCMGET(IPGEOM,'MIX',KNUM)
             DO I=1,NXI
               KNUM(I)=-KNUM(I)
               NUMG(I)=KNUM(I)
             ENDDO
           ENDIF
           IF(LGTURN) THEN
             CALL LCMGET(IPGEOM,'TURN',KNUM(NXI+1))
           ELSE
             CALL XDISET(KNUM(NXI+1),NXI,1)
             LGTURN=.TRUE.
           ENDIF
           IVOL=0
           ISTO=0
           ALLOCATE(ITRN(KXI*3),ISPZ(LZF))
           IF(LSPLIT) THEN
              CALL LCMSIX(IPTRK,'DATA_DUP',1)
              CALL LCMGET(IPTRK,'SPLITZD',ISPZ)
              CALL LCMSIX(IPTRK,' ',2)
           ELSE
              CALL XDISET(ISPZ,LZF,1)
           ENDIF
           ALLOCATE(KVOL(KXI))
           ISS=-1
           ISS1=-1
           IMERMX=0
           DO 182 ISZ=1,LZF
           IFRT=ISS1
           IMER0=KNUM(LXI*ISZ)
           IMER1=IMER0
           DO 482 ISX=1,LXI
             ISS=ISS+1
             ISS1=ISS1+1
             IMERG=KNUM(ISS1+1)
             JMERG=NUMG(ISS1+1)
             MMCYL=ICC(JMERG)
             IF(LSPLIT) THEN
               IMER0=MAX(IMER0,IMERG)
               IMER1=MIN(IMER1,IMERG)
             ENDIF
             IZON=0
             IZZ=0
             LL=0
             DO 171 IFF=1,MMCYL+1
               JSECT=IOF(JMERG)%I2(IFF)
               IF(JSECT.GT.1) THEN
                 LL=JSECT-1
                 IZZ=6*LL
               ELSE
                 LL=1
                 IZZ=1
               ENDIF
               IZON=IZON+IZZ
  171        CONTINUE
             IMERG=IMERG+IMERMX
             IF(LGMERG) THEN
               DO 281 KY=1,ISTO
                 IF(IMERG.EQ.ITRN(KXI+KY)) THEN
                   IMMS=ITRN(2*KXI+KY)
                   IBSS=KVOL(IMMS+1)
                   KVOL(ISS+1)=IBSS
                   ITRN(ISS+1)=ITRN(IMMS+1)
                   GOTO 381
                 ENDIF
  281          CONTINUE
             ENDIF
             ITRN(KXI+ISTO+1)=IMERG
             ITRN(2*KXI+ISTO+1)=ISS
             ISTO=ISTO+1
             KVOL(ISS+1)=IVOL+1
             IVOL=IVOL+IZON
             ITRN(ISS+1)=IVOL-IZZ+LL*ROT(KNUM(NXI+ISS1+1))
  381      CONTINUE
  482      CONTINUE
           IF(LSPLIT) THEN
             IFR0=IFRT
             DO 582 IXX=1,ISPZ(ISZ)-1
             IMERMX=IMERMX+IMER0-IMER1+1
             IFRT=IFR0
           DO 682 ISX=1,LXI
             ISS=ISS+1
             IFRT=IFRT+1
             IMERG=KNUM(IFRT+1)
             JMERG=NUMG(IFRT+1)
             MMCYL=ICC(JMERG)
             IZON=0
             IZZ=0
             LL=0
             DO 671 IFF=1,MMCYL+1
               JSECT=IOF(JMERG)%I2(IFF)
               IF(JSECT.GT.1) THEN
                 LL=JSECT-1
                 IZZ=6*LL
               ELSE
                 LL=1
                 IZZ=1
               ENDIF
               IZON=IZON+IZZ
  671        CONTINUE
             IMERG=IMERG+IMERMX
             IF(LGMERG) THEN
               DO 681 KY=1,ISTO
                 IF(IMERG.EQ.ITRN(KXI+KY)) THEN
                   IMMS=ITRN(2*KXI+KY)
                   IBSS=KVOL(IMMS+1)
                   KVOL(ISS+1)=IBSS
                   ITRN(ISS+1)=ITRN(IMMS+1)
                   GOTO 781
                 ENDIF
  681          CONTINUE
             ENDIF
             ITRN(KXI+ISTO+1)=IMERG
             ITRN(2*KXI+ISTO+1)=ISS
             ISTO=ISTO+1
             KVOL(ISS+1)=IVOL+1
             IVOL=IVOL+IZON
             ITRN(ISS+1)=IVOL-IZZ+LL*ROT(KNUM(NXI+IFRT+1))
  781      CONTINUE
  682      CONTINUE
  582      CONTINUE
           ENDIF
  182      CONTINUE
*--- RECHERCHE DU NOMBRE DE SURFACES QUE PRESENTE LA SYMETRIE
           ICCOIN=-1
           JCCOIN=-1
           KCCOIN=-1
           MROT=6
           NCPER=0
           IF(IHEX.EQ.1) THEN
             NCPER=NINT(REAL(NCOUR)/2.)
           ELSEIF(IHEX.EQ.2) THEN
             NCPER=NCOUR
             MROT=3
           ELSEIF(IHEX.EQ.3) THEN
             NCPER=2*NINT(REAL(NCOUR)/2.)-1
             ICCOIN=LXI-NINT(REAL(NCPER)/2.)+1
             MROT=3
           ELSEIF(IHEX.EQ.4) THEN
             NCPER=3*NINT(REAL(NCOUR)/2.)-2
             IF(MOD(NCOUR,2).EQ.0)NCPER=NCPER+1
             ICCOIN=LXI-NCPER+NINT(REAL(NCOUR)/2.)
             MROT=2
           ELSEIF(IHEX.EQ.5) THEN
             NCPER=1
             IF(NCOUR.GT.1)NCPER=2*(NCOUR-1)
             ICCOIN=LXI-NCPER+NCOUR-1
             KCCOIN=LXI
             MROT=3
           ELSEIF(IHEX.EQ.6) THEN
             NCPER=1
             IF(NCOUR.GT.1)NCPER=4+3*(NCOUR-2)-1
             KCCOIN=LXI-NCPER+NCOUR-1
             ICCOIN=KCCOIN+NCOUR-1
             JCCOIN=ICCOIN+NCOUR-1
             MROT=2
           ELSEIF(IHEX.EQ.7) THEN
             NCPER=1
             IF(NCOUR.GT.1)NCPER=3*NCOUR-2
             ICCOIN=LXI-NCOUR+1
             KCCOIN=ICCOIN-NCOUR+1
             MROT=1
           ELSEIF(IHEX.EQ.8) THEN
             NCPER=1
             IF(NCOUR.GT.1)NCPER=2*(NCOUR+NINT(REAL(NCOUR)/2.))-3
             JCCOIN=LXI-NINT(REAL(NCOUR)/2.)+1
             KCCOIN=JCCOIN-NCOUR+1
             ICCOIN=KCCOIN-NCOUR+1
             MROT=1
           ELSEIF(IHEX.EQ.9) THEN
             MROT=1
             NCPER=1
             IF(NCOUR.GT.1) NCPER=6*(NCOUR-1)
             ICCOIN=IFONC(NCOUR,0)
             JCCOIN=IFONC(NCOUR,1)
             KCCOIN=IFONC(NCOUR,2)
           ELSE
             CALL XABORT('MESHST: INVALID TYPE OF SYMETRIE ')
           ENDIF
           ISURSY=0
           ISAUX=1
           MSMAX=0
           ISSXX=0
           IPPZ=LZF
           IPP1=IPLANZ
           IF(MCODE.GT.0) IPP1=IPLANI
           CALL XDISET(ISXY,6*IPP1,0)
           IXP=0
           DO 500 JP=1,IPPZ
           MXI=LXI*(JP-1)
           DO 501 IP=IXP+1,IXP+ISPZ(JP)
           ISURB6=0
           ITT9=NUMG(LXI+MXI)
           ISURSX=IOF(ITT9)%I2(ICC(ITT9)+1)-1
           KSURBX=ISURSX
           IF(ISURSX.EQ.0) THEN
             ISURSX=1
             KSURBX=1
           ENDIF
           JSECT=0
           DO 183 K1=LXI+MXI,LXI+MXI-NCPER+1,-1
             K=NUMG(K1)
             JSECT=IOF(K)%I2(ICC(K)+1)
             ICX=2
             IF(IHEX.GT.2) THEN
               IF(K1.EQ.ICCOIN+MXI)ICX=3
               IF(K1.EQ.JCCOIN+MXI)ICX=3
               IF(K1.EQ.KCCOIN+MXI)ICX=3
               IF(IHEX.EQ.9) THEN
                 IF(K1.EQ.MXI+IFONC(NCOUR,3))ICX=3
                 IF(K1.EQ.MXI+IFONC(NCOUR,4))ICX=3
                 IF(K1.EQ.MXI+IFONC(NCOUR,5))ICX=3
               ENDIF
               IF(K1.GE.ICCOIN+MXI) THEN
                 IF(IHEX.LE.4) THEN
                   IF(K1.EQ.MXI+LXI.AND.MOD(NCOUR,2).NE.0) THEN
                     ICY=1
                     IF(IHEX.EQ.4)ICY=2
                     ISURB6=ISURB6+ICY*(JSECT-1)
                     IF(JSECT.EQ.1)ISURB6=ISURB6+ICY
                     IF(IHEX.EQ.3) ICX=1
                   ELSE
                     ISURB6=ISURB6+2*(JSECT-1)
                     IF(JSECT.EQ.1)ISURB6=ISURB6+2
                   ENDIF
                 ELSE
                   ICZ=ICX
                   IF(K1.EQ.ICCOIN+MXI)ICZ=2
                   ISURB6=ISURB6+ICZ*(JSECT-1)
                   IF(JSECT.EQ.1)ISURB6=ISURB6+ICZ
                 ENDIF
               ENDIF
             ENDIF
             IF(JSECT.GT.1) THEN
               ISURSY=ISURSY+ICX*(JSECT-1)
             ELSE
               ISURSY=ISURSY+ICX
             ENDIF
  183      CONTINUE
           IDEBX=2
           IF(JSECT.GT.1)IDEBX=JSECT
           IF(IHEX.EQ.1.OR.IHEX.EQ.3.OR.IHEX.EQ.4.OR.IHEX.EQ.8) THEN
           IF(MOD(NCOUR,2).NE.0) THEN
             IF(JSECT.GT.1) THEN
               ISURSY=ISURSY-JSECT+1
             ELSE
               ISURSY=ISURSY-1
             ENDIF
             IF(IHEX.EQ.8) THEN
               ISURSY=ISURSY-ISURSX
               ISURB6=ISURB6-ISURSX
             ENDIF
           ENDIF
           ENDIF
           ISXY(1,IP)=ISURSY
           ISXY(2,IP)=ISURSX
           ISXY(3,IP)=ISURB6
           ISXY(4,IP)=IDEBX+ISAUX-1
           ISXY(5,IP)=KSURBX
           ISXY(6,IP)=ISAUX
           ISAUX=ISURSY+1
           MSMAX=MSMAX+ISURSY-ISSXX
           ISSXX=ISURSY
  501      CONTINUE
           IXP=IXP+ISPZ(JP)
  500      CONTINUE
           DEALLOCATE(KNUM,NUMG,ISPZ)
         ENDIF
         DEALLOCATE(ICELL)
         NBASE=NCEL+NCEL
         IF(LGTURN)NBASE=NBASE+NCEL
         ALLOCATE(IBB(NBASE))
         CALL LCMSIX(IPTRK,'DATA_DUP',1)
         CALL LCMGET(IPTRK,'GENERATING',IBB)
         CALL LCMGET(IPTRK,'MERGE',IBB(NCEL+1))
         IF(LGTURN)CALL LCMGET(IPTRK,'TURN',IBB(2*NCEL+1))
         CALL LCMSIX(IPTRK,' ',2)
         CALL XDISET(VOL1,NCPHY,0)
         CALL XDISET(VTURN,NCPHY,0)
         IAUX=IPOINT
         IFVOL=1
         MVOL=0
         ICELC0=IFONC(NCOUR,0)
         ICELC1=IFONC(NCOUR,1)-ICELC0
         JVOL=1
         IVOL=0
         ISURB=NSMIN
         JSURB=0
         DO 100 I=1,NCEL
            NVOL(I)=JVOL
            IF(LG3D) THEN
              DO 40 J=1,LZZ
                 NCC=J*LX
                 IF(I.LE.NCC) THEN
                   Z=ZZZ(J)
                   GO TO 42
                 ENDIF
   40         CONTINUE
            ENDIF
   42       CONTINUE
            MCYL=ICC(IBB(I))
            IF(MCYL.GT.0) THEN
*
*               STORAGE DES COORDONNEES DES CYLINDRES
*        ==> ATTENTION:  LES AXES DES CYLINDRES SONT SELON Z
*
              DO 43 J=1,MCYL
                 IAUX=IAUX+1
                 REMESH(IAUX)=IOF(IBB(I))%R(J+1)
   43         CONTINUE
            ENDIF
            NCYL(I)=MCYL
            NZONE1=1+MCYL
            NZONE=0
            DO 60 IA=1,NZONE1
              NZONE=NZONE+IOF(IBB(I))%I2(IA)
   60       CONTINUE
            NSECTO=IOF(IBB(I))%I2(NZONE1)
            MCPHY=IBB(NCEL+I)
            IV1=VOL1(MCPHY)
            IF(IV1.GT.0) THEN
              IF(LGTURN) THEN
                JVT=PHTURN(MCPHY)
                IVTURN=VTURN(MCPHY)
                ITURN=IBB(2*NCEL+I)
                KAUX=NSECTO-1
                IF(ITURN.LE.6) THEN
                  IF(ITURN.NE.6)KAUX=KAUX*(ITURN+1)
                ELSEIF(ITURN.GE.9) THEN
                  KAUX=(ITURN-8)*KAUX
                ELSEIF(ITURN.EQ.7) THEN
                  KAUX=5*KAUX
                ENDIF
              ENDIF
            ELSE
              KVOL1=1
              IF(MCPHY.GT.1)KVOL1=VOL1(MCPHY-1)
              KVOL2=KVOL1+NZONE-1
              IF(LGTURN) THEN
                ITURN=IBB(2*NCEL+I)
                PHTURN(MCPHY)=ITURN
                IF(ITURN.LE.6) THEN
                  IF(ITURN.EQ.6) THEN
                    VTURN(MCPHY)=KVOL2+NSECTO-1
                  ELSE
                    VTURN(MCPHY)=KVOL2+(NSECTO-1)*(ITURN+1)
                  ENDIF
                ELSEIF(ITURN.GE.8) THEN
                  VTURN(MCPHY)=KVOL2+(NSECTO-1)*(ITURN-8)+1
                ELSE
                  VTURN(MCPHY)=KVOL2+(NSECTO-1)*5+1
                ENDIF
              ENDIF
              VOL1(MCPHY)=KVOL2+1
            ENDIF
            VOLUM1=0.
            IF(LG3D) THEN
              IF(IPLANZ.GT.1) THEN
                IF(I.LE.STAIRS(1)) THEN
                  FACB(I)=KSURF1
                ELSEIF(I.GT.STAIRS(IPLANZ-1)) THEN
                  ISUR=ISUR+1
                  FACB(LX+ISUR)=KSURF2
                ENDIF
              ELSE
                FACB(I)=KSURF1
                FACB(LX+I)=KSURF2
              ENDIF
              Z=REMESH(2*NCEL+I)
              IF(I.GT.STAIRS(1)) THEN
                IST=1
                DO 329 IY=2,LZZ
                  IF(I.LE.STAIRS(IY)) THEN
                    IST=IY-1
                    GOTO 339
                  ENDIF
  329           CONTINUE
  339           KZAUX=I-STAIRS(IST)
                IF(IST.GT.1)KZAUX=KZAUX+STAIRS(IST-1)
                Z=Z-REMESH(2*NCEL+KZAUX)
              ENDIF
            ENDIF
            KSECTX=IOF(IBB(I))%I2(NZONE1)
            NZZ=NZONE1
            IF(KSECTX.GT.3)NZZ=NZONE1-1
            DO 45 J=1,NZZ
               MIXX=IOF(IBB(I))%I1(J)
               KSECT=IOF(IBB(I))%I2(J)
               NAUX=1
               IF(KSECT.GT.1) NAUX=6*(KSECT-1)
               IVOL=IVOL+NAUX
               MVOL=MVOL+1
               NSECT(MVOL)=KSECT
               KREG=KREG+1
               IF(J.EQ.NZONE1) THEN
                 VOLUM2=VEX
               ELSE
                 R=IOF(IBB(I))%R(J+1)
                 VOLUM2=PI*R*R
               ENDIF
               VOLUMS=0.0
               VOLUME=(VOLUM2-VOLUM1)/REAL(NAUX)
               IF(LG3D) THEN
                 VOLUMS=VOLUME*.25
                 VOLUME=VOLUME*Z
*
*     SURFACES SUPERIEURES ET INFERIEURES (SELON L'AXE Z)
*
                 IF(IPLANZ.GT.1) THEN
                   IF(I.LE.STAIRS(1)) THEN
                     KSURF1=KSURF1+1
                     NSURF1=NSURF1-1
                     VOLSUR(NSURF1)=VOLUMS
                   ELSEIF(I.GT.STAIRS(IPLANZ-1)) THEN
                     KSURF2=KSURF2+1
                     NSURF2=NSURF2-1
                     VOLSUR(NSURF2)=VOLUMS
                   ENDIF
                 ELSE
                   KSURF1=KSURF1+1
                   NSURF1=NSURF1-1
                   VOLSUR(NSURF1)=VOLUMS
                   KSURF2=KSURF2+1
                   NSURF2=NSURF2-1
                   VOLSUR(NSURF2)=VOLUMS
                 ENDIF
               ENDIF
               VOLUM1=VOLUM2
               VOLSUR(KREG)=VOLUME
               MATALB(KREG)=MIXX
               DO 44 K=2,NAUX
                  KREG=KREG+1
                  VOLSUR(KREG)=VOLUME
                  MATALB(KREG)=MIXX
                  IF(LG3D) THEN
                  IF(IPLANZ.GT.1) THEN
                    IF(I.LE.STAIRS(1)) THEN
                      KSURF1=KSURF1+1
                      NSURF1=NSURF1-1
                      VOLSUR(NSURF1)=VOLUMS
                    ELSEIF(I.GT.STAIRS(IPLANZ-1)) THEN
                      KSURF2=KSURF2+1
                      NSURF2=NSURF2-1
                      VOLSUR(NSURF2)=VOLUMS
                    ENDIF
                  ELSE
                    KSURF1=KSURF1+1
                    NSURF1=NSURF1-1
                    VOLSUR(NSURF1)=VOLUMS
                    KSURF2=KSURF2+1
                    NSURF2=NSURF2-1
                    VOLSUR(NSURF2)=VOLUMS
                  ENDIF
                  ENDIF
   44          CONTINUE
   45       CONTINUE
            ICXX=LX*INT(AINT(REAL(I/(LX+1))))
            ICELC=ICELC0+ICXX
            ICLIM=LX+ICXX
            IF(KSECTX.GT.3) THEN
*
* TRAITEMENT DU VOLUME A BORDURE HEXAGONALE DANS LE CAS OU KSECT>3
*
             KSECT=KSECTX
             MIXX=IOF(IBB(I))%I1(NZONE1+1)
             NAUX=6*(KSECT-1)
             IVOL=IVOL+NAUX
             MVOL=MVOL+1
             NSECT(MVOL)=KSECT
             VCYL=VOLUM1/REAL(NAUX)
             SAUX=0.
             IXX0=KREG
             ISS0=NSURF1
             JSS0=NSURF2
             DO 445 M=1,KSECT-1
               IXX=IXX0
               ISS1=ISS0
               ISS2=JSS0
               XTAN=TAN(REAL(M)*PI/(3*(KSECT-1)))
               SAUX1=XTAN/(1.+XTAN/SQRT3)
               VOLUME=.5*SIDE*SIDE*(SAUX1-SAUX)
               VOLUME=VOLUME-VCYL
               VOLUMS=VOLUME
               IF(LG3D) VOLUME=VOLUME*Z
               DO 444  K=1,6
                 IXX=IXX+1
                 KREG=KREG+1
                 VOLSUR(IXX)=VOLUME
                 MATALB(IXX)=MIXX
                 IF(LG3D) THEN
                 IF(IPLANZ.GT.1) THEN
                   IF(I.LE.STAIRS(1)) THEN
                     KSURF1=KSURF1+1
                     NSURF1=NSURF1-1
                     ISS1=ISS1-1
                     VOLSUR(ISS1)=VOLUMS
                   ELSEIF(I.GT.STAIRS(IPLANZ-1)) THEN
                     KSURF2=KSURF2+1
                     NSURF2=NSURF2-1
                     ISS2=ISS2-1
                     VOLSUR(ISS2)=VOLUMS
                   ENDIF
                 ELSE
                   KSURF1=KSURF1+1
                   NSURF1=NSURF1-1
                   ISS1=ISS1-1
                   VOLSUR(ISS1)=VOLUMS
                   KSURF2=KSURF2+1
                   NSURF2=NSURF2-1
                   ISS2=ISS2-1
                   VOLSUR(ISS2)=VOLUMS
                 ENDIF
                 ISS1=ISS1-KSECT+2
                 ISS2=ISS2-KSECT+2
                 ENDIF
                 IXX=IXX+KSECT-2
 444           CONTINUE
               IF(I.GE.ICELC.AND.I.LE.ICLIM) THEN
                 SURF=.5*SIDE*(SAUX1-SAUX)/SQRT3
                 IF(LG3D)SURF=SURF*Z
                 JSUR=JSUR-1
                 VOLSUR(JSUR)=SURF
                 IJSUR=JSUR-KSECT+1
                 VOLSUR(IJSUR)=SURF
                 IF(MOD(I-ICELC,ICELC1).EQ.0) THEN
                   IJSUR=IJSUR-KSECT+1
                   VOLSUR(IJSUR)=SURF
                 ENDIF
               ENDIF
               IXX0=IXX0+1
               ISS0=ISS0-1
               JSS0=JSS0-1
               SAUX=SAUX1
 445         CONTINUE
             IF(I.GE.ICELC.AND.I.LE.ICLIM) THEN
               JSUR=JSUR-KSECT+1
               JSURB=JSURB+2
               SURB(JSURB-1)=ISURB
               KKB=KSECT-1
               ISURB=ISURB+KKB+KKB
               SURB(JSURB)=ISURB-KKB
               IF(MOD(I-ICELC,ICELC1).EQ.0) THEN
                 JSUR=JSUR-KSECT+1
                 JSURB=JSURB+1
                 SURB(JSURB)=ISURB
                 ISURB=ISURB+KKB
               ENDIF
             ENDIF
            ELSE
            IF(I.GE.ICELC.AND.I.LE.ICLIM) THEN
              SURF=SIDE*.25
              IF(LG3D)SURF=SURF*Z
              KKB=1
              IF(KSECT.GT.1) THEN
                KKB=KSECT-1
                SURF=SURF/KKB
              ENDIF
              DO 99 IXX=1,KKB
              JSUR=JSUR-2
              VOLSUR(JSUR)=SURF
              VOLSUR(JSUR+1)=SURF
              IF(MOD(I-ICELC,ICELC1).EQ.0) THEN
                JSUR=JSUR-1
                VOLSUR(JSUR)=SURF
              ENDIF
   99         CONTINUE
              JSURB=JSURB+2
              SURB(JSURB-1)=ISURB
              ISURB=ISURB+KKB+KKB
              SURB(JSURB)=ISURB-KKB
              IF(MOD(I-ICELC,ICELC1).EQ.0) THEN
                JSURB=JSURB+1
                SURB(JSURB)=ISURB
                ISURB=ISURB+KKB
              ENDIF
            ENDIF
            ENDIF
            JVOL=IVOL+1
            FVOL(I)=IFVOL
            IFVOL=IFVOL+NZONE1
  100    CONTINUE
         DO 110 I=1,ISTAT9
            DEALLOCATE(IOF(I)%I1,IOF(I)%I2)
            IF(ICC(I).GT.0) DEALLOCATE(IOF(I)%R)
  110    CONTINUE
         DEALLOCATE(IOF,ICC)
      ELSE
*
*           CAS DE CELLULE HOMOGENES
*
         ALLOCATE(IBB(2*NCEL))
         CALL LCMSIX(IPTRK,'DATA_DUP',1)
         CALL LCMGET(IPTRK,'MERGE',IBB)
         CALL LCMGET(IPTRK,'MIX',IBB(NCEL+1))
         CALL LCMSIX(IPTRK,' ',2)
         CALL LCMLEN(IPGEOM,'SECTOR',ILENS,ITP)
         ALLOCATE(ICC(NCEL))
         IF(ILENS.GT.0) THEN
            CALL XABORT('MESHST: SECTORS NOT IMPLEMENTED.')
         ELSE
            CALL XDISET(ICC,NCEL,1)
         ENDIF
         IF(LGTURN) THEN
           CALL LCMSIX(IPTRK,'DATA_DUP',1)
           ALLOCATE(ITT(NCEL))
           CALL LCMGET(IPTRK,'TURN',ITT)
           CALL LCMSIX(IPTRK,' ',2)
         ENDIF
         CALL XDISET(VOL1,NCPHY,0)
         IFVOL=1
         MVOL=0
         ICELC0=IFONC(NCOUR,0)
         ICELC1=IFONC(NCOUR,1)-ICELC0
         JVOL=1
         IVOL=0
         ISURB=NSMIN
         JSURB=0
         DO 200 I=1,NCEL
            NVOL(I)=JVOL
            NSECTO=ICC(I)
            ICXX=LX*INT(AINT(REAL(I/(LX+1))))
            ICELC=ICELC0+ICXX
            ICLIM=LX+ICXX
            NSECT(I)=NSECTO
            NCYL(I)=0
            NZONE=1
            IF(NSECTO.GT.1)NZONE=6*(NSECTO-1)
            IVOL=IVOL+NZONE
            JVOL=IVOL+1
            MCPHY=IBB(I)
            IV1=VOL1(MCPHY)
            IF(IV1.GT.0) THEN
              IF(LGTURN) THEN
                JVT=PHTURN(MCPHY)
                IVTURN=VTURN(MCPHY)
                KAUX=NSECTO-1
                ITURN=ITT(I)
                IF(ITURN.LE.6) THEN
                  IF(ITURN.NE.6)KAUX=KAUX*(ITURN+1)
                ELSEIF(ITURN.GE.9) THEN
                  KAUX=(ITURN-8)*KAUX
                ELSEIF(ITURN.EQ.7) THEN
                  KAUX=5*KAUX
                ENDIF
              ENDIF
            ELSE
              KVOL1=1
              IF(MCPHY.GT.1)KVOL1=VOL1(MCPHY-1)
              KVOL2=KVOL1+NZONE-1
              IF(LGTURN) THEN
                ITURN=ITT(I)
                PHTURN(MCPHY)=ITURN
                IF(ITURN.LE.6) THEN
                  IF(ITURN.EQ.6) THEN
                    VTURN(MCPHY)=KVOL2+NSECTO-1
                  ELSE
                    VTURN(MCPHY)=KVOL2+(NSECTO-1)*(ITURN+1)
                  ENDIF
                ELSEIF(ITURN.GE.8) THEN
                  VTURN(MCPHY)=KVOL2+(NSECTO-1)*(ITURN-8)+1
                ELSE
                  VTURN(MCPHY)=KVOL2+(NSECTO-1)*5+1
                ENDIF
              ENDIF
              VOL1(MCPHY)=KVOL2+1
            ENDIF
            IF(NSECTO.LE.3) THEN
            IF(I.GE.ICELC.AND.I.LE.ICLIM) THEN
              SURF=SIDE*.25
              IF(LG3D)SURF=SURF*Z
              KKB=1
              IF(NSECTO.GT.1) THEN
                KKB=NSECTO-1
                SURF=SURF/(NSECTO-1)
              ENDIF
              DO 556 IXX=1,KKB
              JSUR=JSUR-2
              VOLSUR(JSUR)=SURF
              VOLSUR(JSUR+1)=SURF
              IF(MOD(I-ICELC,ICELC1).EQ.0) THEN
                JSUR=JSUR-1
                VOLSUR(JSUR)=SURF
              ENDIF
  556         CONTINUE
              JSURB=JSURB+2
              SURB(JSURB-1)=ISURB
              ISURB=ISURB+KKB+KKB
              SURB(JSURB)=ISURB-KKB
              IF(MOD(I-ICELC,ICELC1).EQ.0) THEN
                JSURB=JSURB+1
                SURB(JSURB)=ISURB
                ISURB=ISURB+KKB
              ENDIF
            ENDIF
            VOLUMS=0.0
            VOLUME=VEX/NZONE
            IF(LG3D) THEN
              Z=REMESH(2*NCEL+I)
              IF(I.GT.STAIRS(1)) THEN
                IST=1
                DO 429 IY=2,LZZ
                  IF(I.LE.STAIRS(IY)) THEN
                    IST=IY-1
                    GOTO 439
                  ENDIF
  429           CONTINUE
  439           KZAUX=I-STAIRS(IST)
                IF(IST.GT.1)KZAUX=KZAUX+STAIRS(IST-1)
                Z=Z-REMESH(2*NCEL+KZAUX)
              ENDIF
              VOLUMS=VOLUME*.25
              VOLUME=VOLUME*Z
*
*     SURFACES SUPERIEURES ET INFERIEURES (SELON L'AXE Z)
*
              IF(IPLANZ.GT.1) THEN
                IF(I.LE.STAIRS(1)) THEN
                   FACB(I)=KSURF1
                   KSURF1=KSURF1+1
                   NSURF1=NSURF1-1
                   VOLSUR(NSURF1)=VOLUMS
                 ELSEIF(I.GT.STAIRS(IPLANZ-1)) THEN
                   ISUR=ISUR+1
                   FACB(LX+ISUR)=KSURF2
                   KSURF2=KSURF2+1
                   NSURF2=NSURF2-1
                   VOLSUR(NSURF2)=VOLUMS
                 ENDIF
               ELSE
                 FACB(I)=KSURF1
                 FACB(LX+I)=KSURF2
                 KSURF1=KSURF1+1
                 NSURF1=NSURF1-1
                 VOLSUR(NSURF1)=VOLUMS
                 KSURF2=KSURF2+1
                 NSURF2=NSURF2-1
                 VOLSUR(NSURF2)=VOLUMS
               ENDIF
             ENDIF
             MIXX=IBB(NCEL+I)
             DO 65 J=1,NZONE
               KREG=KREG+1
               VOLSUR(KREG)=VOLUME
               MATALB(KREG)=MIXX
               IF(LG3D) THEN
               IF(IPLANZ.GT.1) THEN
                 IF(I.LE.STAIRS(1)) THEN
                   KSURF1=KSURF1+1
                   NSURF1=NSURF1-1
                   VOLSUR(NSURF1)=VOLUMS
                 ELSEIF(I.GT.STAIRS(IPLANZ-1)) THEN
                   KSURF2=KSURF2+1
                   NSURF2=NSURF2-1
                   VOLSUR(NSURF2)=VOLUMS
                 ENDIF
               ELSE
                 KSURF1=KSURF1+1
                 NSURF1=NSURF1-1
                 VOLSUR(NSURF1)=VOLUMS
                 KSURF2=KSURF2+1
                 NSURF2=NSURF2-1
                 VOLSUR(NSURF2)=VOLUMS
               ENDIF
               ENDIF
   65        CONTINUE
             ELSE
*
*  TRAITEMENT DU CAS OU NSECTO>3
*
             KSECT=NSECTO
             MIXX=IBB(NCEL+I)
             NAUX=6*(KSECT-1)
             IVOL=IVOL+NAUX
             MVOL=MVOL+1
             NSECT(MVOL)=KSECT
             SAUX=0.
             IXX0=KREG
             ISS0=NSURF1
             JSS0=NSURF2
             DO 555 M=1,KSECT-1
               IXX=IXX0
               ISS1=ISS0
               ISS2=JSS0
               XTAN=TAN(REAL(M)*PI/(3*(KSECT-1)))
               SAUX1=XTAN/(1.+XTAN/SQRT3)
               VOLUME=.5*SIDE*SIDE*(SAUX1-SAUX)
               VOLUMS=VOLUME
               IF(LG3D) VOLUME=VOLUME*Z
               DO 554  K=1,6
                 IXX=IXX+1
                 KREG=KREG+1
                 VOLSUR(IXX)=VOLUME
                 MATALB(IXX)=MIXX
                 IF(LG3D) THEN
                 IF(IPLANZ.GT.1) THEN
                   IF(I.LE.STAIRS(1)) THEN
                     KSURF1=KSURF1+1
                     NSURF1=NSURF1-1
                     ISS1=ISS1-1
                     VOLSUR(ISS1)=VOLUMS
                   ELSEIF(I.GT.STAIRS(IPLANZ-1)) THEN
                     KSURF2=KSURF2+1
                     NSURF2=NSURF2-1
                     ISS2=ISS2-1
                     VOLSUR(ISS2)=VOLUMS
                   ENDIF
                 ELSE
                   KSURF1=KSURF1+1
                   NSURF1=NSURF1-1
                   ISS1=ISS1-1
                   VOLSUR(ISS1)=VOLUMS
                   KSURF2=KSURF2+1
                   NSURF2=NSURF2-1
                   ISS2=ISS2-1
                   VOLSUR(ISS2)=VOLUMS
                 ENDIF
                 ISS1=ISS1-KSECT+2
                 ISS2=ISS2-KSECT+2
                 ENDIF
                 IXX=IXX+KSECT-2
 554           CONTINUE
               IF(I.GE.ICELC.AND.I.LE.ICLIM) THEN
                 SURF=.5*SIDE*(SAUX1-SAUX)/SQRT3
                 IF(LG3D)SURF=SURF*Z
                 JSUR=JSUR-1
                 VOLSUR(JSUR)=SURF
                 IJSUR=JSUR-KSECT+1
                 VOLSUR(IJSUR)=SURF
                 IF(MOD(I-ICELC,ICELC1).EQ.0) THEN
                   IJSUR=IJSUR-KSECT+1
                   VOLSUR(IJSUR)=SURF
                 ENDIF
               ENDIF
               IXX0=IXX0+1
               ISS0=ISS0-1
               JSS0=JSS0-1
               SAUX=SAUX1
 555         CONTINUE
             IF(I.GE.ICELC.AND.I.LE.ICLIM) THEN
               JSUR=JSUR-KSECT+1
               JSURB=JSURB+2
               SURB(JSURB-1)=ISURB
               KKB=KSECT-1
               ISURB=ISURB+KKB+KKB
               SURB(JSURB)=ISURB-KKB
               IF(MOD(I-ICELC,ICELC1).EQ.0) THEN
                 JSUR=JSUR-KSECT+1
                 JSURB=JSURB+1
                 SURB(JSURB)=ISURB
                 ISURB=ISURB+KKB
               ENDIF
             ENDIF
             ENDIF
             FVOL(I)=IFVOL
             IFVOL=IFVOL+1
  200    CONTINUE
         DEALLOCATE(ICC)
         IF(LGTURN) DEALLOCATE(ITT)
      ENDIF
      IF(LG3D) DEALLOCATE(ZZZ)
*
*--- CONSTRUCTION DU VECTEUR VSYM QUI CONTIENT POUR CHAQUE VOLUME DE
*--- DE L'ASSEMBLAGE SON EQUIVALENT DANS LA SYMETRIE D'ENTREE
*
      IF(LGCELL.AND.(IHEX.LE.9)) THEN
        ALLOCATE(IAA(NCEL))
        CALL LCMSIX(IPTRK,'DATA_DUP',1)
        CALL LCMGET(IPTRK,'GENER0',IAA)
        CALL LCMSIX(IPTRK,' ',2)
        IXV=0
        DO 300 IC=1,NCEL
           MCYL=NCYL(IC)
           LFROT=ROT(IBB(2*NCEL+IC))
           IVFIN=KVOL(IAA(IC))
           IVORIN=ITRN(IAA(IC))
           IVMIN=IVFIN
           IVSYM=NVOL(IC)
           IXX=IXV
           NZON=0
           DO 310 I=1,MCYL
             IXX=IXX+1
             LSECT=NSECT(IXX)
             IF(LSECT.GT.1) THEN
               NZON=NZON+6*(LSECT-1)
             ELSE
               NZON=NZON+1
             ENDIF
  310      CONTINUE
           LSECT=NSECT(IXX+1)
           IF(LSECT.GT.1) THEN
             IVAUX=INT(REAL(IVORIN-IVFIN-NZON+1)/REAL(LSECT-1))
           ELSE
             IVAUX=IVORIN-IVFIN-NZON+1
           ENDIF
           DO 320 K=1,MCYL+1
              IXV=IXV+1
              LSECT=NSECT(IXV)
              IF(LSECT.EQ.1) THEN
                VSYM(IVSYM)=IVMIN
                IVMAX=IVMIN
              ELSE
                IVOR=IVMIN+IVAUX*(LSECT-1)-1
                IVV=LFROT-1
                IVSYM0=IVSYM
                IVSYM=IVSYM+6*(LSECT-1)-1
                IVMAX=IVMIN+6*(LSECT-1)-1
                IF(IBB(2*NCEL+IC).LE.6) THEN
                  IVV1=LFROT*(LSECT-1)
                  IVV=(5-IVV)*(LSECT-1)
                ELSE
                  IVV1=IVV*(LSECT-1)+1
                  IVV=(6-IVV)*(LSECT-1)-1
                ENDIF
                IF(IBB(2*NCEL+IC).LE.6) THEN
                  KVV=IVSYM-IVV
                  KVOR=IVOR
                  DO 315 L=KVV,IVSYM
                    VSYM(L)=KVOR
                    KVOR=KVOR+1
                    IF(KVOR.GT.IVMAX)KVOR=IVMIN
  315             CONTINUE
                  DO 316 L=IVSYM0,KVV-1
                    VSYM(L)=KVOR
                    KVOR=KVOR+1
                    IF(KVOR.GT.IVMAX)KVOR=IVMIN
  316             CONTINUE
                ELSE
                  KVV=IVSYM0+IVV1-1
                  KVOR=IVOR
                  DO 317 L=KVV,IVSYM0,-1
                    VSYM(L)=KVOR
                    KVOR=KVOR+1
                    IF(KVOR.GT.IVMAX)KVOR=IVMIN
  317             CONTINUE
                  DO 318 L=IVSYM,KVV+1,-1
                    VSYM(L)=KVOR
                    KVOR=KVOR+1
                    IF(KVOR.GT.IVMAX)KVOR=IVMIN
  318             CONTINUE
                ENDIF
              ENDIF
              IVMIN=IVMAX+1
              IVSYM=IVSYM+1
  320      CONTINUE
  300   CONTINUE
        DEALLOCATE(IAA)
*--- CONSTRUCTION DU VECTEUR SSYM QUI CONTIENT POUR CHAQUE SURFACE
*--- DE L'ASSEMBLAGE SON EQUIVALENT DANS LA SYMETRIE D'ENTREE
        ISY=0
        MSAUX=NS
        ISYX=NSMIN+1
        IVLMAX=0
        IF(LG3D) THEN
          MSAUX=(NS-NSMIN-NSMAX)/IPLANZ
*--- CAS DE LA SYMETRIE Z- (MCODE=1) SSYM EST REMPLIT PLUS LOIN
          IF(MCODE.NE.1) THEN
          IF(IHEX.NE.9) THEN
          DO 301 ISUR=1,NSMIN
             SSYM(ISUR)=VSYM(ISUR)
             IVLMAX=MAX(IVLMAX,VSYM(ISUR))
  301     CONTINUE
          ELSE
          DO 311 ISUR=1,NSMIN
             SSYM(ISUR)=ISUR
  311     CONTINUE
          IVLMAX=NSMIN
          ENDIF
          ENDIF
          ISY=NSMIN
        ENDIF
        IPPZ=IPP1
        NSAUX=0
        IF(MCODE.EQ.1) THEN
          DO 302 IPPX=IPPZ+1,IPLANZ
            NSAUX=NSAUX+(SURL(IPPX)-SURL(IPPX-1))
  302     CONTINUE
          ISYX=ISYX+NSAUX
        ENDIF
        MSAUX=0
        ISYAUX=ISYX
        LSMAX=NSMIN+NSAUX
        IF(MCODE.EQ.0) LSMAX=NS
        LSPLZM=2*IPPZ-IPLANZ
        DO 600 IP=1,IPPZ
         IF(LG3D) THEN
          IF(MCODE.EQ.1) THEN
            LSMAX=LSMAX+(SURL(IPPZ-LSPLZM+IP)-SURL(IPPZ-LSPLZM-1+IP))
            ISYAUX=ISYAUX+MSAUX
            MSAUX=SURL(IPPZ-LSPLZM+IP)-SURL(IPPZ-LSPLZM-1+IP)
          ELSE
            LSMAX=NSMIN+SURL(IP)
            ISYAUX=ISYX+MSAUX
            MSAUX=SURL(IP)
          ENDIF
         ENDIF
         ISY=ISYAUX-1
        ISURSY=ISXY(1,IP)
        ISURSX=ISXY(2,IP)
        ISURB6=ISXY(3,IP)
        IDEBX=ISXY(4,IP)
        KSURBX=ISXY(5,IP)
        ISAUX=ISXY(6,IP)
        IDDX=ISURSX
        KX2=ISURSY
        KX1=ISURSY-ISURSX+1
        KXP=1
        IF(IHEX.GE.3) THEN
           ISY=ISURB6+ISY
           IF(IHEX.EQ.3.OR.IHEX.EQ.8) THEN
             IDDX=0
             KX2=KX1
             KX1=ISURSY
             KXP=-1
           ENDIF
           IF(IHEX.EQ.4.OR.IHEX.EQ.7)ISY=ISY-KSURBX
        ENDIF
        IDEB=ISAUX
        IF(IHEX.EQ.2.OR.IHEX.EQ.7)IDEB=IDEBX
        IF(IHEX.EQ.5.OR.IHEX.EQ.6) THEN
          IDDX=-ISURSY
          ISURSX=0
          ISURSY=0
          KX1=0
          KX2=-1
        ENDIF
        DO 400 I=1,MROT
        DO 349 K=KX1,KX2,KXP
          ISY=ISY+1
          IF(IHEX.GE.3) THEN
            IF(ISY.GT.LSMAX) ISY=ISYAUX
          ENDIF
          SSYM(ISY)=K+IVLMAX
  349   CONTINUE
        IF(IHEX.LT.9) THEN
        DO 350 K=ISURSY-ISURSX,ISAUX,-1
          ISY=ISY+1
          IF(IHEX.GE.3) THEN
            IF(ISY.GT.LSMAX) ISY=ISYAUX
          ENDIF
          SSYM(ISY)=K+IVLMAX
  350   CONTINUE
        ENDIF
        DO 351 K=IDEB,ISURSY-IDDX
          ISY=ISY+1
          IF(IHEX.GE.3) THEN
            IF(ISY.GT.LSMAX) ISY=ISYAUX
          ENDIF
          SSYM(ISY)=K+IVLMAX
  351   CONTINUE
  400   CONTINUE
  600   CONTINUE
*
*---SURFACE LATERALES SUPERIEURES DANS LE CAS 3D
*
        IF(LG3D) THEN
        IYAUX=0
        IF(MCODE.LE.1) THEN
        IF(IPLANZ.GT.1) THEN
          IYAUX=IVLMAX
          DO 599 II=NSMIN+1,NV-NSMAX
            IYAUX=MAX(VSYM(II),IYAUX)
  599     CONTINUE
        ENDIF
        IF(IHEX.NE.9) THEN
          IAUX=NV-NSMAX
          DO 601 II=NS-NSMAX+1,NS
           IAUX=IAUX+1
           SSYM(II)=IVLMAX+MSMAX+VSYM(IAUX)-IYAUX
  601     CONTINUE
        ELSE
          IAUX=0
          IF(MCODE.GT.0) THEN
           KSAUX=SURL(IPLANI)
          ELSE
           KSAUX=SSYM(NS-NSMAX)
          ENDIF
          DO 611 II=NS-NSMAX+1,NS
           IAUX=IAUX+1
           SSYM(II)=IAUX+KSAUX
  611     CONTINUE
        ENDIF
        ENDIF
*--- CAS DE SYMETRIE Z-
        IF(MCODE.EQ.1) THEN
          IF(IHEX.NE.9) THEN
            DO 602 IV=1,NSMAX
             SSYM(IV)=IVLMAX+MSMAX+VSYM(IV)-IYAUX
  602       CONTINUE
          ELSE
            IVLAT=1
            IAUX=0
            JTX=NS-NSMAX
            DO 612 IX1=1,NCOUR
              DO 613 JX=JTX+VLAT(NCOUR+IX1),JTX+IVLAT,-1
                IAUX=IAUX+1
                SSYM(IAUX)=SSYM(JX)
  613         CONTINUE
            DO 614 JX=JTX+VLAT(IX1),JTX+VLAT(NCOUR+IX1)+1,-1
                IAUX=IAUX+1
                SSYM(IAUX)=SSYM(JX)
  614       CONTINUE
            IVLAT=VLAT(IX1)+1
  612       CONTINUE
          ENDIF
          IFR=NSMIN
          LSPLZM=2*IPLANI-IPLANZ
          DO 605 IW=IPLANZ,IPPZ+1,-1
            DO 603 JW=SURL(IPLANZ+IW),SURL(IW-1)+1,-1
               IFR=IFR+1
               SSYM(IFR)=SSYM(NSMIN+JW)
  603       CONTINUE
            DO 604 JW=SURL(IW),SURL(IPLANZ+IW)+1,-1
               IFR=IFR+1
               SSYM(IFR)=SSYM(NSMIN+JW)
  604       CONTINUE
  605     CONTINUE
        ELSEIF(MCODE.EQ.2) THEN
          IFR=NS-NSMAX+1
          MSUR=0
          LSPLZP=2*IPLANI-IPLANZ
          DO 705 IW=1,IPPZ-LSPLZP
            DO 704 JW=SURL(IPLANZ+IW)+1,SURL(IW)
               IFR=IFR-1
               SSYM(IFR)=SSYM(NSMIN+JW)
  704       CONTINUE
            DO 703 JW=MSUR+1,SURL(IPLANZ+IW)
               IFR=IFR-1
               SSYM(IFR)=SSYM(NSMIN+JW)
  703       CONTINUE
            MSUR=SURL(IW)
  705     CONTINUE
          IF(IHEX.NE.9) THEN
            IAUX=NV-NSMIN
            DO 706 II=NS-NSMAX+1,NS
             IAUX=IAUX+1
             SSYM(II)=VSYM(IAUX)
  706       CONTINUE
          ELSE
            IVLAT=1
            IAUX=NS-NSMAX
            DO 712 IX1=1,NCOUR
              DO 713 JX=VLAT(NCOUR+IX1),IVLAT,-1
                IAUX=IAUX+1
                SSYM(IAUX)=SSYM(JX)
  713         CONTINUE
            DO 714 JX=VLAT(IX1),VLAT(NCOUR+IX1)+1,-1
                IAUX=IAUX+1
                SSYM(IAUX)=SSYM(JX)
  714       CONTINUE
            IVLAT=VLAT(IX1)+1
  712       CONTINUE
          ENDIF
        ENDIF
        ENDIF
        DEALLOCATE(IBB,KVOL,ITRN)
      ELSE
        DO 717 ISUR=1,NS
          SSYM(ISUR)=ISUR
  717   CONTINUE
        DO 718 IVOL=1,NV
          VSYM(IVOL)=IVOL
  718   CONTINUE
      ENDIF
      DEALLOCATE(ISXY)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  800 DEALLOCATE(PHTURN,VTURN,VOL1)
      RETURN
      END
