*DECK PSPTHR
      SUBROUTINE PSPTHR(IPTRK ,ISPSP ,IPRINT,ICEL  ,NDIM  ,NFREG ,
     >                  MAXMSH,MXGREG,MAXPIN,KPSP  ,COLREG,FACT  ,
     >                  CELLPO,IDREG ,ITPIN ,DCMESH,DRAPIN,
     >                  NBPTS  ,POSTRI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To draw an hexagonal cell with triangular mesh according 
* to its explicit position in the assembly.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal.
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure in
*         update or creation mode.
* ISPSP   POSTSCRIPT file index.
* IPRINT  print level.
* ICEL    cell number.
* NDIM    problem dimensions.
* NFREG   number of regions.
* MAXMSH  maximum number of elements in MESH array.
* MXGREG  maximum number of region for any geometry.
* MAXPIN  maximum number of pins in a cell.
* KPSP    PSP plot options.
* FACT    scale factor for drawing.
* CELLPO  global cell position in space.
* COLREG  region color.
*
*Parameters: temporary storage
* IDREG   local region identifier.
* ITPIN   pin type identifier.
* DCMESH  meshing vector for geometries.
* DRAPIN  pin position identifier.
* NBPTS   number of corners for regions in first sector.
* POSTRI  positions of corners for regions in first sector.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK,ISPSP
      INTEGER          IPRINT,ICEL,NDIM,NFREG,MAXMSH,MXGREG,MAXPIN
      INTEGER          KPSP(7)
      REAL             COLREG(4,NFREG)
      DOUBLE PRECISION FACT,CELLPO(2,2)
      INTEGER          IDREG(MXGREG),ITPIN(3,MAXPIN)
      DOUBLE PRECISION DCMESH(-1:MAXMSH,4),DRAPIN(-1:4,MAXPIN)
      INTEGER          NBPTS(MXGREG)
      DOUBLE PRECISION POSTRI(2,4,MXGREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='PSPTHR')
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
      REAL             WLINE
      PARAMETER       (WLINE=0.002)
      DOUBLE PRECISION DZERO,DONE,DTWO,DHALF,DSQ3O2
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0,
     >                 DHALF=0.5D0,DSQ3O2=0.86602540378444D0)
      INTEGER          MAXDIM
      PARAMETER       (MAXDIM=4)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
*----
*  Local variables
*----
      INTEGER          ILEV,ICONT,ICOL,KFS,KFR,KSS,KSR,NPTS,NSEG,NINT
      INTEGER          IEDIMC(NSTATE),IEDIMP(NSTATE)
      INTEGER          NX,MESHC(4),NREGC,NTPIN,NRTP
      INTEGER          IDIR,ILXY,IREG,IR,IY,IX,ILOC,IORDER(16)
      INTEGER          JR,NAR,ITREG,IAREG,JREG
      REAL             WLFAC,XYPOS(2,4),CENTER(2),RCIRC,RADANG(2,16)
      CHARACTER        NAMCEL*9,NAMREC*12
      REAL             COLWHI(4),CENTEP(2),CENTED(2),CENTEB(2)
      INTEGER          IWCOL,IPIN,ISEG,NBRP,MESHP(4)
      DOUBLE PRECISION ROTAX,COSDIR(3)
      INTEGER          ISECT,ITRI,KREG,NCPT
      DOUBLE PRECISION SIDET,SIDEL,SIDEH,POSCXD,POSCYD,POSCX,POSCY,
     >                 SIDEU,SIDEUH,DTMPX,DTMPY
*----
*  Data
*----
      CHARACTER        CDIR(MAXDIM)*1
      SAVE             CDIR
      CHARACTER        CLEV(2)*1
      SAVE             CLEV
      DATA             CDIR /'X','Y','Z','R'/
      DATA             CLEV /'C','P'/
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      IF(NDIM .GE. 3) CALL XABORT(NAMSBR//
     >': PSP cannot treat 3D geometries')
      PI=XDRCST('Pi',' ')
      ILEV=1
      IWCOL=1
      CALL XDRSET(COLWHI,4,1.0)
*----
*  PSP print control
*----
      WLFAC=1.0
      ICONT=KPSP(1)
      ICOL=KPSP(2)
      IF(KPSP(3) .EQ. 1) WLFAC=2.5
      KFS=KPSP(4)
      KFR=KPSP(5)
      KSS=KPSP(6)
      KSR=KPSP(7)
      NPTS=3
      NINT=16
*----
*  Read cell information
*----
      WRITE(NAMCEL,'(A1,I8.8)') CLEV(ILEV),ICEL
      NAMREC=NAMCEL//'DIM'
      CALL XDISET(IEDIMC,NSTATE,0)
      CALL LCMGET(IPTRK,NAMREC,IEDIMC)
      NX=IEDIMC(3)
      NAR=IEDIMC(2)
      MESHC(1)=NX
      MESHC(2)=NX
      MESHC(3)=IEDIMC(5)
      MESHC(4)=NAR
      NREGC=IEDIMC(8)
      NTPIN=IEDIMC(18)
      NAMREC=NAMCEL//'RID'
      CALL LCMGET(IPTRK,NAMREC,IDREG)
      DO IDIR=1,4
        NAMREC=NAMCEL//'SM'//CDIR(IDIR)
        IF(MESHC(IDIR) .GT. 0) THEN
          CALL LCMGET(IPTRK,NAMREC,DCMESH(-1,IDIR))
        ENDIF
      ENDDO
      IF(NTPIN .GT .0) THEN
        NAMREC=NAMCEL//'PIN'
        CALL LCMGET(IPTRK,NAMREC,DRAPIN)
        NAMREC=NAMCEL//'PNT'
        CALL LCMGET(IPTRK,NAMREC,ITPIN)
      ENDIF
      CENTER(1)=REAL(FACT*(DCMESH(-1,1)+CELLPO(1,1)))
      CENTER(2)=REAL(FACT*(DCMESH(-1,2)+CELLPO(2,1)))
      CENTEB(1)=-CENTER(1)
      CENTEB(2)=-CENTER(2)
      NRTP=NX**2
      SIDEU=DCMESH(1,1)-DCMESH(0,1)
      SIDET=SIDEU
      IF(NX .GT. 1) THEN
        SIDET=DCMESH(2,1)-DCMESH(1,1)
      ENDIF
      SIDEL=SIDET/DSQ3O2
      SIDEH=DHALF*SIDEL
*----
*  Find triangle corners
*  1- First sectors
*     All crown
*----
      IREG=0
      POSCXD=DZERO
      POSCYD=SIDEH
      DO IX=1,NX
        POSCX=POSCXD
        POSCY=POSCYD
*----
*  right triangles
*----
        DO IR=1,IX-1
          IREG=IREG+1
          POSTRI(1,1,IREG)=POSCX+SIDET
          POSTRI(2,1,IREG)=POSCY
          POSTRI(1,2,IREG)=POSCXD
          POSTRI(2,2,IREG)=POSCY+SIDEH
          POSTRI(1,3,IREG)=POSCXD
          POSTRI(2,3,IREG)=POSCY-SIDEH
          NBPTS(IREG)=3
          POSCY=POSCY+SIDEL
        ENDDO
*----
*  Loop over Left triangles on the line
*----
        POSCYD=POSCYD-SIDEH
        POSCX=POSCXD
        POSCY=POSCYD
        DO IR=1,IX
          IREG=IREG+1
          POSTRI(1,1,IREG)=POSCX
          POSTRI(2,1,IREG)=POSCY
          POSTRI(1,2,IREG)=POSCX+SIDET
          POSTRI(2,2,IREG)=POSCY-SIDEH
          POSTRI(1,3,IREG)=POSCX+SIDET
          POSTRI(2,3,IREG)=POSCY+SIDEH
          POSCY=POSCY+SIDEL
          NBPTS(IREG)=3
        ENDDO
        POSCXD=POSCXD+SIDET
      ENDDO
*----
*  Replace right triangles in last crown by 
*  Isosceles trapezoid when required
*----
      IF(SIDEU .LT. SIDET) THEN
        SIDEU=SIDET-SIDEU
        SIDEUH=DHALF*SIDEU/DSQ3O2
        IX=NX
        IREG=IREG-2*IX+1
*----
*  Cut left triangle in last crown when required
*----
        DO IR=1,IX-1
          IREG=IREG+1
          POSTRI(1,4,IREG)=POSTRI(1,1,IREG)-SIDEU
          POSTRI(2,4,IREG)=POSTRI(2,1,IREG)-SIDEUH
          POSTRI(1,1,IREG)=POSTRI(1,1,IREG)-SIDEU
          POSTRI(2,1,IREG)=POSTRI(2,1,IREG)+SIDEUH
          NBPTS(IREG)=4
        ENDDO
*----
*  Loop over Left triangles on the line
*----
        DO IR=1,IX
          IREG=IREG+1
          POSTRI(1,2,IREG)=POSTRI(1,2,IREG)-SIDEU
          POSTRI(2,2,IREG)=POSTRI(2,2,IREG)+SIDEUH
          POSTRI(1,3,IREG)=POSTRI(1,3,IREG)-SIDEU
          POSTRI(2,3,IREG)=POSTRI(2,3,IREG)-SIDEUH
        ENDDO
      ENDIF
*----
*  General location of regions in first sector with respect to
*  cell center are now known
*  Draw all regions in first sector
*----
      IREG=0
      ISECT=1
      DO IR=1,NRTP
        IREG=IREG+1
        ITREG=IREG*(NAR+1)
        NCPT=NBPTS(IR) 
        KREG=ABS(IDREG(ITREG))
        DO ITRI=1,NCPT
          XYPOS(1,ITRI)=REAL(FACT*(POSTRI(1,ITRI,IR)+CELLPO(1,1)))
          XYPOS(2,ITRI)=REAL(FACT*(POSTRI(2,ITRI,IR)+CELLPO(2,1)))
        ENDDO
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6010) ISECT,IR,IREG,ITREG,KREG,NCPT
          WRITE(IOUT,6011) (XYPOS(1,ITRI),XYPOS(2,ITRI),ITRI=1,NCPT)
        ENDIF
        IF(KREG .NE. 0) THEN
*----
*  Color and trace result
*----
          CALL PSDREG(ISPSP,NCPT,XYPOS)
          IF(ICOL. GT. 0) THEN
            CALL PSFILL(ISPSP,ICOL,COLREG(1,KREG),KFS,KFR)
          ENDIF
          IF(ICONT.EQ.1) THEN
            CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
          ENDIF
        ENDIF
        DO JR=NAR,1,-1
          IAREG=(NAR+1)*(IREG-1)+JR
          JREG=ABS(IDREG(IAREG))
          IF(IPRINT .GE. 100) THEN
            WRITE(IOUT,6020) JR,IAREG,JREG,FACT*DCMESH(JR,4)
          ENDIF
*----
*  Annular region
*----
          IF(JREG .NE. 0) THEN
            RCIRC=REAL(FACT*DCMESH(JR,4))
*----
*  Move cursor to center of annulus
*----
            CALL PSMOVE(ISPSP,CENTER,-3)
*----
*  Color and trace result
*----
            CALL PSPRAI(NINT,NPTS,XYPOS,CENTER,RCIRC,
     >                  NSEG,IORDER,RADANG)
            CALL PSDRAI(ISPSP,NSEG,IORDER,CENTER,RADANG)
            IF(ICOL. GT. 0) THEN
              CALL PSFILL(ISPSP,ICOL,COLREG(1,JREG),KFS,KFR)
            ENDIF
            IF(ICONT.EQ.1) THEN
              CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
            ENDIF
*----
*  Return cursor to original position
*----
            CALL PSMOVE(ISPSP,CENTEB,-3)
          ENDIF
        ENDDO
      ENDDO
*----
*  2) Rotate regions for remaining sectors
*     and draw
*----
      DO ISECT=2,6
        DO IR=1,NRTP
          IREG=IREG+1
          ITREG=IREG*(NAR+1)
          NCPT=NBPTS(IR) 
          KREG=ABS(IDREG(ITREG))
          DO ITRI=1,NCPT
*----
*  rotation by pi/3 from previous sector
*----
            DTMPX=DHALF*POSTRI(1,ITRI,IR)
     >           -DSQ3O2*POSTRI(2,ITRI,IR)
            DTMPY=DSQ3O2*POSTRI(1,ITRI,IR)
     >           +DHALF*POSTRI(2,ITRI,IR)
            POSTRI(1,ITRI,IR)=DTMPX
            POSTRI(2,ITRI,IR)=DTMPY
            XYPOS(1,ITRI)=REAL(FACT*(POSTRI(1,ITRI,IR)+CELLPO(1,1)))
            XYPOS(2,ITRI)=REAL(FACT*(POSTRI(2,ITRI,IR)+CELLPO(2,1)))
          ENDDO
          IF(IPRINT .GE. 100) THEN
            WRITE(IOUT,6010) ISECT,IR,IREG,ITREG,KREG,NCPT
            WRITE(IOUT,6011) (XYPOS(1,ITRI),XYPOS(2,ITRI),ITRI=1,NCPT)
          ENDIF
          IF(KREG .NE. 0) THEN
*----
*  Color and trace result
*----
            CALL PSDREG(ISPSP,NCPT,XYPOS)
            IF(ICOL. GT. 0) THEN
              CALL PSFILL(ISPSP,ICOL,COLREG(1,KREG),KFS,KFR)
            ENDIF
            IF(ICONT.EQ.1) THEN
              CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
            ENDIF
          ENDIF
          DO JR=NAR,1,-1
            IAREG=(NAR+1)*(IREG-1)+JR
            JREG=ABS(IDREG(IAREG))
            IF(IPRINT .GE. 100) THEN
              WRITE(IOUT,6020) JR,IAREG,JREG,DCMESH(JR,4)
            ENDIF
*----
*  Annular region
*----
            IF(JREG .NE. 0) THEN
              RCIRC=REAL(FACT*DCMESH(JR,4))
*----
*  Move cursor to center of annulus
*----
              CALL PSMOVE(ISPSP,CENTER,-3)
*----
*  Color and trace result
*----
              CALL PSPRAI(NINT,NPTS,XYPOS,CENTER,RCIRC,
     >                    NSEG,IORDER,RADANG)
              CALL PSDRAI(ISPSP,NSEG,IORDER,CENTER,RADANG)
              IF(ICOL. GT. 0) THEN
                CALL PSFILL(ISPSP,ICOL,COLREG(1,JREG),KFS,KFR)
              ENDIF
              IF(ICONT.EQ.1) THEN
                CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
              ENDIF
*----
*  Return cursor to original position
*----
              CALL PSMOVE(ISPSP,CENTEB,-3)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
*----
* Pins
*----
      IF(NTPIN .GT. 0) THEN
        NPTS=4
        ILEV=2
        MESHP(1)=1
        MESHP(2)=1
        MESHP(3)=1
        MESHP(4)=1
        CENTEP(1)=0.0
        CENTEP(2)=0.0
        DO IPIN=1,NTPIN
*----
*  Locate pin position
*----
          COSDIR(1)=DRAPIN(0,IPIN)*COS(DRAPIN(-1,IPIN))
          COSDIR(2)=DRAPIN(0,IPIN)*SIN(DRAPIN(-1,IPIN))
          CENTED(1)=CENTER(1)+REAL(FACT*COSDIR(1))
          CENTED(2)=CENTER(2)+REAL(FACT*COSDIR(2))
          RCIRC=REAL(FACT*DRAPIN(4,IPIN))
          ROTAX=PI/DTWO-DRAPIN(-1,IPIN)
*----
*  Move cursor to center of pin
*----
          CALL PSMOVE(ISPSP,CENTED,-3)
*----
*  Read pin information
*----
          WRITE(NAMCEL,'(A1,I8.8)') CLEV(ILEV),ITPIN(2,IPIN)
          NAMREC=NAMCEL//'DIM'
          CALL XDISET(IEDIMP,NSTATE,0)
          CALL LCMGET(IPTRK,NAMREC,IEDIMP)
          MESHP(1)=IEDIMP(3)
          MESHP(2)=IEDIMP(4)
          MESHP(3)=IEDIMP(5)
          MESHP(4)=IEDIMP(2)
          NBRP=MESHP(4)
          NAMREC=NAMCEL//'RID'
          CALL LCMGET(IPTRK,NAMREC,IDREG)
          DO IDIR=1,4
            NAMREC=NAMCEL//'SM'//CDIR(IDIR)
            IF(MESHP(IDIR) .GT. 0) THEN
              CALL LCMGET(IPTRK,NAMREC,DCMESH(-1,IDIR))
            ENDIF
          ENDDO
          DO IY=MESHP(2),1,-1
            DO IX=MESHP(1),1,-1
              ILXY=((IY-1)*MESHP(1)+(IX-1))*NBRP
*----
*  Cartesian region
*----
              XYPOS(1,1)=REAL(FACT*(DCMESH(IX-1,1)))
              XYPOS(2,1)=REAL(FACT*(DCMESH(IY-1,2)))
              XYPOS(1,2)=REAL(FACT*(DCMESH(IX,1)))
              XYPOS(2,2)=XYPOS(2,1)
              XYPOS(1,3)=XYPOS(1,2)
              XYPOS(2,3)=REAL(FACT*(DCMESH(IY,2)))
              XYPOS(1,4)=XYPOS(1,1)
              XYPOS(2,4)=XYPOS(2,3)
              DO IR=MESHP(4),1,-1
                ILOC=ILXY+IR
                IREG=ABS(IDREG(ILOC))
                RCIRC=REAL(FACT*DCMESH(IR,4))
*----
*  Annular pin regions
*----
                IF(IREG .NE. 0) THEN
                  CALL PSPRAI(NINT,NPTS,XYPOS,CENTEP,RCIRC,
     >                        NSEG,IORDER,RADANG)
*----
*  Rotate pins intersection points
*----
                  DO ISEG=1,NSEG
                    RADANG(2,ISEG)=RADANG(2,ISEG)-REAL(ROTAX)
                  ENDDO
*----
*  Color and trace result
*----
                  CALL PSDRAI(ISPSP,NSEG,IORDER,CENTEP,RADANG)
                  IF(ICOL. GT. 0) THEN
                    CALL PSFILL(ISPSP,ICOL,COLREG(1,IREG),KFS,KFR)
                  ENDIF
                  IF(ICONT.EQ.1) THEN
                    CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO
*----
*  Return cursor to original position
*----
          CENTED(1)=-CENTED(1)
          CENTED(2)=-CENTED(2)
          CALL PSMOVE(ISPSP,CENTED,-3)
        ENDDO
      ENDIF
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(' Sector ',I5,3X,'Triangle ',I5,5X,'T region ',I5,5X,
     >'G region ',I5,5X,
     >'Region number',I5,5X,'Number of corners ',I5)
 6011 FORMAT(2F20.7)
 6020 FORMAT(' Annulus ',I5,3X,'A Region ',I5,5X,
     >'Region number',I5,5X,'radius ',F20.10)
      END
