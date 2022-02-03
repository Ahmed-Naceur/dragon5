*DECK NXTEGI
      SUBROUTINE NXTEGI(IPGEO ,IPRINT,ITYPG ,MAXMSH,NMIX  ,NM    ,
     >                  MAXMSS,NMS   ,NREG  ,NREGS ,NSUR  ,NSURS ,
     >                  MIX   ,ISPLT ,DAMESH,
     >                  RMESH ,MIXC  )
*
*-----------------------------------------------------------------------
*
*Purpose:
* To extract cell or pin geometry information.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPGEO   pointer to the reference geometry data structure.
* IPRINT  intermediate printing level for output.
* ITYPG   geometry type.
* MAXMSH  maximum number of elements in MESH array.
* NMIX    number of elements in MIX array.
* NM      mesh size in all directions ($X$, $Y$, $Z$ and $R$).
*
*Parameters: output
* MAXMSS  maximum number of elements in MESH array after split.
* NMS     mesh size in all directions ($X$, $Y$, $Z$ and $R$)
*         after split.
* NREG    number of regions.
* NREGS   number of regions after split.
* NSUR    number of surfaces.
* NSURS   number of surfaces after split.
* MIX     final mixture description for geometry  (including HMIX).
* ISPLT   final split desctiption for geometry.
* DAMESH  final mesh description for geometry.
*
*Parameters: temporary storage
* RMESH   temporary vector for reading cell mesh array.
* MIXC    temporary mixture for cell rotation  (including HMIX).
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPGEO
      INTEGER          IPRINT,ITYPG,MAXMSH,NMIX,NM(4),MAXMSS,
     >                 NMS(4),NREG,NREGS,NSUR,NSURS
      INTEGER          MIX(NMIX,2),ISPLT(MAXMSH,4)
      DOUBLE PRECISION DAMESH(-1:MAXMSH,4)
      REAL             RMESH(0:MAXMSH)
      INTEGER          MIXC(NMIX,2,2)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTEGI')
      DOUBLE PRECISION DZERO,DONE
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0)
      DOUBLE PRECISION DSQ3O2
      PARAMETER       (DSQ3O2=0.86602540378444D0)
*----
*  Local variables
*----
      INTEGER          IDIR,IR,IX,IY,IZ,IMTN,IMTO
      INTEGER          ILCMLN,ILCMTY,IMRGLN,IMRGTY
      INTEGER          NX,NY,NZ,NR,NRM,NXS,NYS,NZS,NRS,NRMS,
     >                 NMREAD(4)
      CHARACTER        NAMREC*12,NAMMRG*12
      REAL             OFFCEN(3)
      REAL             SIDE,SIDET
      DOUBLE PRECISION DSIDE,DSIDET
*----
*  Data
*----
      CHARACTER        CDIR(4)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      NRM=0
      NRMS=0
      NX=MAX(1,NM(1))
      NY=MAX(1,NM(2))
      NZ=MAX(1,NM(3))
      NR=MAX(1,NM(4))
*----
*  Read geometry information
*  1- Cartesian MESH
*----
      IF(ITYPG .EQ.  8 .OR. ITYPG .EQ.  9 .OR.
     >   ITYPG .EQ. 12 .OR. ITYPG .EQ. 13 .OR.
     >   ITYPG .EQ. 26 .OR. ITYPG .EQ. 27) THEN
*----
*  Hexagons
*----
        IDIR=1
        NAMREC='SIDE        '
        CALL LCMGET(IPGEO,NAMREC,SIDE)
        DSIDE=DBLE(SIDE)*DSQ3O2
        IF(ITYPG .EQ. 12 .OR. ITYPG .EQ. 13 .OR.
     >     ITYPG .EQ. 26 .OR. ITYPG .EQ. 27) THEN
          NAMREC='SIDET       '
          CALL LCMGET(IPGEO,NAMREC,SIDET)
          DSIDET=DBLE(SIDET)*DSQ3O2
          DAMESH(0,IDIR)=-DSIDE
          DAMESH(1,IDIR)=-DSIDET*(NM(1)-1)
          DO IX=2,2*NM(1)-1
            DAMESH(IX,IDIR)=DAMESH(IX-1,IDIR)+DSIDET
          ENDDO
          DAMESH(2*NM(1),IDIR)=DSIDE
        ELSE
          DAMESH(-1,IDIR)=2*DSIDE
          DAMESH(0,IDIR)=-DSIDE
          DAMESH(1,IDIR)=DSIDE
        ENDIF
        IDIR=2
        DO IX=-1,2*NM(1)
          DAMESH(IX,IDIR)=DAMESH(IX,IDIR-1)
        ENDDO
        IDIR=3
        NMREAD(IDIR)=0
        IF(NM(IDIR) .GT. 0) THEN
          NAMREC='MESH'//CDIR(IDIR)//'       '
          CALL LCMLEN(IPGEO,NAMREC,NMREAD(IDIR),ILCMTY)
          IF(NMREAD(IDIR) .EQ. NM(IDIR)+1) THEN
            CALL LCMGET(IPGEO,NAMREC,RMESH)
            DO IX=0,NM(IDIR)
              DAMESH(IX,IDIR)=DBLE(RMESH(IX))
            ENDDO
            CALL XDISET(ISPLT(1,IDIR),NM(IDIR),1)
          ENDIF
        ELSE
          DAMESH(0,IDIR)=DZERO
          IX=1
          DAMESH(IX,IDIR)=DONE
        ENDIF
      ELSE
*----
*  Parallepiped
*----
        DO IDIR=1,3
          NMREAD(IDIR)=0
          IF(NM(IDIR) .GT. 0) THEN
            NAMREC='MESH'//CDIR(IDIR)//'       '
            CALL LCMLEN(IPGEO,NAMREC,NMREAD(IDIR),ILCMTY)
            IF(NMREAD(IDIR) .EQ. NM(IDIR)+1) THEN
              CALL LCMGET(IPGEO,NAMREC,RMESH)
              DO IX=0,NM(IDIR)
                DAMESH(IX,IDIR)=DBLE(RMESH(IX))
              ENDDO
              CALL XDISET(ISPLT(1,IDIR),NM(IDIR),1)
            ENDIF
          ELSE
            DAMESH(0,IDIR)=DZERO
            IX=1
            DAMESH(IX,IDIR)=DONE
          ENDIF
        ENDDO
      ENDIF
*----
*  2- Read cell OFFCENTER and store in position -1 of DAMESH
*----
      CALL XDRSET(OFFCEN,3,0.0)
      CALL LCMLEN(IPGEO,'OFFCENTER   ',ILCMLN,ILCMTY)
      IF(ILCMLN .GT. 0)
     >  CALL LCMGET(IPGEO,'OFFCENTER   ',OFFCEN)
      DO IDIR=1,3
        DAMESH(-1,IDIR)=DBLE(OFFCEN(IDIR))
      ENDDO
*----
*  3- Radial mesh
*----
      IDIR=4
      NMREAD(IDIR)=0
      NAMREC='RADIUS      '
      CALL LCMLEN(IPGEO,NAMREC,NMREAD(IDIR),ILCMTY)
      IF(NMREAD(IDIR) .EQ. NM(IDIR)+1) THEN
        CALL LCMGET(IPGEO,NAMREC,RMESH)
        DO IX=0,NM(IDIR)
          DAMESH(IX,IDIR)=DBLE(RMESH(IX))
        ENDDO
        CALL XDISET(ISPLT(1,IDIR),NM(IDIR),1)
      ENDIF
*----
*  4- Cartesian, radial and hexagonal split
*----
      DO IDIR=1,4
        IF(NM(IDIR) .GT. 0) THEN
          NAMREC='SPLIT'//CDIR(IDIR)//'      '
          CALL XDISET(ISPLT(1,IDIR),NM(IDIR),1)
          CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
          IF(ILCMLN .EQ. NM(IDIR))
     >      CALL LCMGET(IPGEO,NAMREC,ISPLT(1,IDIR))
          NMS(IDIR)=0
          DO IR=1,NM(IDIR)
            NMS(IDIR)=NMS(IDIR)+ABS(ISPLT(IR,IDIR))
          ENDDO
        ELSE
          NMS(IDIR)=NM(IDIR)
        ENDIF
      ENDDO
      IF(ITYPG .EQ. 12 .OR. ITYPG .EQ. 13 .OR.
     >   ITYPG .EQ. 26 .OR. ITYPG .EQ. 27) THEN
        NAMREC='SPLITH      '
        IDIR=1
        CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
        CALL XDISET(ISPLT(1,IDIR),1,1)
        IF(ILCMLN .EQ. NM(IDIR))
     >  CALL LCMGET(IPGEO,NAMREC,ISPLT(1,IDIR))
        NMS(IDIR)=NM(IDIR)*ABS(ISPLT(1,IDIR))
      ENDIF
      NXS=MAX(1,NMS(1))
      NYS=MAX(1,NMS(2))
      NZS=MAX(1,NMS(3))
      NRS=MAX(1,NMS(4))
      IF(IPRINT .GE. 100) THEN
        IF(ITYPG .EQ.  8 .OR. ITYPG .EQ.  9 .OR.
     >     ITYPG .EQ. 12 .OR. ITYPG .EQ. 13 .OR.
     >     ITYPG .EQ. 26 .OR. ITYPG .EQ. 27) THEN
          WRITE(IOUT,6015) ITYPG,NM(1),NMREAD(3),NMIX
          IF(NM(1) .GT. 0) THEN
            WRITE(IOUT,6011) 'MESHH ='
            WRITE(IOUT,6012) (DAMESH(IX,1),IX=-1,2*NM(1))
            WRITE(IOUT,6011) 'SPLTH ='
            WRITE(IOUT,6013) (ISPLT(1,1),IX=1,2*NM(1))
          ENDIF
          IF(NMREAD(3) .GT. 0) THEN
            WRITE(IOUT,6011) 'MESHZ ='
            WRITE(IOUT,6012) (DAMESH(IX,3),IX=-1,NM(3))
            WRITE(IOUT,6011) 'SPLTZ ='
            WRITE(IOUT,6013) (ISPLT(IX,3),IX=1,NM(3))
          ENDIF
        ELSE
          WRITE(IOUT,6010) ITYPG,(CDIR(IDIR),NM(IDIR),IDIR=1,4),NMIX
          DO IDIR=1,3
            IF(NMREAD(IDIR) .GT. 0) THEN
              WRITE(IOUT,6011) 'MESH'//CDIR(IDIR)//' ='
              WRITE(IOUT,6012) (DAMESH(IX,IDIR),IX=-1,NM(IDIR))
              WRITE(IOUT,6011) 'SPLT'//CDIR(IDIR)//' ='
              WRITE(IOUT,6013) (ISPLT(IX,IDIR),IX=1,NM(IDIR))
            ENDIF
          ENDDO
        ENDIF
        IDIR=4
        IF(NMREAD(IDIR) .GT. 0) THEN
          WRITE(IOUT,6011) 'RADIUS='
          WRITE(IOUT,6012) (DAMESH(IX,IDIR),IX=0,NM(IDIR))
          WRITE(IOUT,6011) 'SPLT'//CDIR(IDIR)//' ='
          WRITE(IOUT,6013) (ISPLT(IX,IDIR),IX=1,NM(IDIR))
        ENDIF
      ENDIF
*----
*  5- Get MIX
*----
      NAMREC='MIX         '
      CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
      IF(ILCMLN .LT. 0 .OR. ILCMLN .GT. NMIX) CALL XABORT(NAMSBR//
     >': Size of MIX vector is invalid')
      NAMMRG='HMIX        '
      CALL LCMLEN(IPGEO,NAMMRG,IMRGLN,IMRGTY)
      IF(IMRGLN .LE. 0 ) THEN
        NAMMRG=NAMREC
      ELSE IF(IMRGLN .NE. ILCMLN) THEN
        NAMMRG=NAMREC
        WRITE(IOUT,8000) NAMSBR
      ENDIF
      IF(ILCMLN .GT. 0) THEN
        IF     (ITYPG .EQ.  3 ) THEN
*----
*  TUBE
*----
          NRM=NR
          NSUR=NZ
          NRMS=NRS
          NSURS=NZS
          IF(ILCMLN .LT. NY*NX*NRM) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,NY*NX*NRM
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE IF(ITYPG .EQ.  6 ) THEN
*----
*  TUBEZ
*----
          NRM=NR
          NSUR=2*(NRM*NX*NY)+NZ
          NRMS=NRS
          NSURS=2*(NRMS*NXS*NYS)+NZS
          IF(ILCMLN .LT. NZ*NY*NX*NRM) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,NZ*NY*NX*NRM
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE IF(ITYPG .EQ.  5) THEN
*----
*  CAR2D
*----
          NRM=1
          NSUR=2*(NY+NX)
          NRMS=1
          NSURS=2*(NYS+NXS)
          IF(ILCMLN .LT. NY*NX) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,NY*NX
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE IF(ITYPG .EQ.  7 ) THEN
*----
*  CAR3D
*----
          NRM=1
          NSUR=2*(NX*NY+NY*NZ+NZ*NX)
          NRMS=1
          NSURS=2*(NXS*NYS+NYS*NZS+NZS*NXS)
          IF(ILCMLN .LT. NZ*NY*NX) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,NZ*NY*NX
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE IF(ITYPG .EQ. 20) THEN
*----
*  CARCEL
*----
          NRM=NR+1
          NSUR=2*(NY+NX)
          NRMS=NRS+1
          NSURS=2*(NYS+NXS)
          IF(ILCMLN .LT. NY*NX*NRM) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,NY*NX*NRM
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE IF(ITYPG .EQ. 23 ) THEN
*----
*  CARCELZ
*----
          NRM=NR+1
          NSUR=2*(NX*NY*NRM+NY*NZ+NZ*NX)
          NRMS=NRS+1
          NSURS=2*(NXS*NYS*NRMS+NYS*NZS+NZS*NXS)
          IF(ILCMLN .LT. NZ*NY*NX*NRM) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,NZ*NY*NX*NRM
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE IF(ITYPG .EQ. 10 .OR. ITYPG .EQ. 21) THEN
*----
*  TUBEX and CARCELX
*----
          IF(ITYPG .EQ.21) THEN
            NRM=NR+1
            NSUR=2*(NX*NY+NY*NZ*NRM+NZ*NX)
            NRMS=NRS+1
            NSURS=2*(NXS*NYS+NYS*NZS*NRMS+NZS*NXS)
          ELSE
            NRM=NR
            NSUR=NX+2*NY*NZ*NRM
            NRMS=NRS
            NSURS=NXS+2*NYS*NZS*NRMS
          ENDIF
*----
*  For CARCELX reorder mixtures from $(R,Y,Z,X)$ to $(R,X,Y,Z)$
*----
          IF(ILCMLN .LT. NZ*NY*NX*NRM) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,NZ*NY*NX*NRM
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIXC(1,1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIXC(1,2,1))
          IMTN=0
          DO 20 IZ=1,NZ
            DO 21 IY=1,NY
              DO 22 IX=1,NX
                DO 23 IR=1,NRM
                  IMTN=IMTN+1
                  IMTO=(IX-1)*NY*NZ*NRM
     >                +(IZ-1)*NY*NRM
     >                +(IY-1)*NRM+IR
                  MIX(IMTN,1)=MIXC(IMTO,1,1)
                  MIX(IMTN,2)=MIXC(IMTO,2,1)
 23             CONTINUE
 22           CONTINUE
 21         CONTINUE
 20       CONTINUE
        ELSE IF(ITYPG .EQ. 11 .OR. ITYPG .EQ. 22) THEN
*----
*  TUBEY and CARCELY
*----
          IF(ITYPG .EQ.22) THEN
            NRM=NR+1
            NSUR=2*(NX*NY+NY*NZ+NZ*NX*NRM)
            NRMS=NRS+1
            NSURS=2*(NXS*NYS+NYS*NZS+NZS*NXS*NRMS)
          ELSE
            NRM=NR
            NSUR=NY+2*NZ*NX*NRM
            NRMS=NRS
            NSURS=NYS+2*NZS*NXS*NRMS
          ENDIF
*----
*  For CARCELX reorder mixtures from $(R,Z,X,Y)$ to $(R,X,Y,Z)$
*----
          IF(ILCMLN .LT. NZ*NY*NX*NRM) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,NZ*NY*NX*NRM
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIXC(1,1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIXC(1,2,1))
          IMTN=0
          DO 30 IZ=1,NZ
            DO 31 IY=1,NY
              DO 32 IX=1,NX
                DO 33 IR=1,NRM
                  IMTN=IMTN+1
                  IMTO=(IY-1)*NZ*NX*NRM
     >                +(IX-1)*NZ*NRM
     >                +(IZ-1)*NRM+IR
                  MIX(IMTN,1)=MIXC(IMTO,1,1)
                  MIX(IMTN,2)=MIXC(IMTO,2,1)
 33             CONTINUE
 32           CONTINUE
 31         CONTINUE
 30       CONTINUE
        ELSE IF(ITYPG .EQ.  8) THEN
*----
*  HEX
*----
          NRM=1
          NSUR=6
          NRMS=1
          NSURS=6
          IF(ILCMLN .LT. 1) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,1
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE IF(ITYPG .EQ.  9) THEN
*----
*  HEXZ
*----
          NRM=1
          NSUR=6*NZ+2
          NRMS=1
          NSURS=6*NX+2
          IF(ILCMLN .LT. NZ) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,NZ
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
        ELSE IF(ITYPG .EQ. 12) THEN
*----
*  HEXT
*----
          NRM=1
          NSUR=6*(2*NX-1)
          NRMS=1
          NSURS=6*(2*NXS-1)
          IF(ILCMLN .LT. 6*NX*NX) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,6*NX*NX
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE IF(ITYPG .EQ. 13) THEN
*----
*  HEXTZ
*----
          NRM=1
          NSUR=6*(2*NX-1)*NZ+12*NX*NX
          NRMS=1
          NSURS=6*(2*NXS-1)*NZS+12*NXS*NXS
          IF(ILCMLN .LT. 6*NX*NX*NZ) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,6*NX*NX*NZ
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE IF(ITYPG .EQ. 26) THEN
*----
*  HEXTCEL
*----
          NRM=NR+1
          NSUR=6*(2*NX-1)
          NRMS=NRS+1
          NSURS=6*(2*NXS-1)
          IF(ILCMLN .LT. 6*NX*NX*NRM) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,6*NX*NX*NRM
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE IF(ITYPG .EQ. 27) THEN
*----
*  HEXTCELZ
*----
          NRM=NR+1
          NSUR=6*(2*NX-1)*NZ+12*NX*NX*NRM
          NRMS=NRS+1
          NSURS=6*(2*NXS-1)*NZS+12*NXS*NXS*NRMS
          IF(ILCMLN .LT. 6*NX*NX*NZ*NRM) THEN
            WRITE(IOUT,9000) NAMSBR,ILCMLN,6*NX*NX*NZ*NRM
            CALL XABORT(NAMSBR//
     >    ': Invalid number of mixtures provided')
          ENDIF
          CALL LCMGET(IPGEO,NAMREC,MIX(1,1))
          CALL LCMGET(IPGEO,NAMMRG,MIX(1,2))
        ELSE
          CALL XABORT(NAMSBR//
     >    ': Geometry type invalid for cell or pin')
        ENDIF
      ENDIF
      IF(ITYPG .EQ. 12 .OR. ITYPG .EQ. 13 .OR.
     >   ITYPG .EQ. 26 .OR. ITYPG .EQ. 27) THEN
        NREG=6*NX*NX*NZ*NRM
        NREGS=6*NXS*NXS*NZS*NRMS
        MAXMSS=MAX(NRMS,2*(NXS+1),2*(NYS+1),NZS,MAXMSH)+1
      ELSE
        NREG=NRM*NX*NY*NZ
        NREGS=NRMS*NXS*NYS*NZS
        MAXMSS=MAX(NRMS,NXS,NYS,NZS,MAXMSH)+1
      ENDIF
*----
*  Print mesh if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6011) 'MIX   ='
        WRITE(IOUT,6013) (MIX(IX,1),IX=1,NMIX)
        WRITE(IOUT,6011) 'HMIX  ='
        WRITE(IOUT,6013) (MIX(IX,2),IX=1,NMIX)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(1X,'Geometry type =',I10/
     >       1X,'Original mesh dimensions ='/
     >       4(1X,A1,'=',1X,I8)/
     >       1X,'Number of regions =',i8)
 6011 FORMAT(1X,A7)
 6012 FORMAT(5F15.9)
 6013 FORMAT(5I15)
 6015 FORMAT(1X,'Geometry type =',I10/
     >       1X,'Original hexagonal mesh dimensions =',I10,/
     >       1X,'Original z mesh dimensions         =',I10,/
     >       1X,'Number of regions =',i8)
 8000 FORMAT(' ***** Warning in ',A6,' *****'/
     >       ' HMIX not compatible with MIX '/
     >       ' HMIX mixture are replaced by MIX mixtures' )
 9000 FORMAT(' ***** Error in ',A6,' *****'/
     >       ' Number of mixtures provided = ',I10/
     >       ' Number of mixtures required = ',I10)
      END
