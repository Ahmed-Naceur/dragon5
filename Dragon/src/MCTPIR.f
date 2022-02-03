*DECK MCTPIR
      SUBROUTINE MCTPIR(IPTRK,IPRINT,NDIM,MAXMSH,MXGSUR,MXGREG,NTPIN,
     1           NBIND,INDX,ITPIN,DRAPIN,DCMESH,MESHP,NSURP,NREGP,
     2           PINCEN,INDEX,IDREG,POSC,IDIRP,INPIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Search if point is within a pin. If so, load pin contents and change
* coordinates to pin local ones.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* IPRINT  print level.
* NDIM    problem dimensions.
* MAXMSH  maximum number of elements in MESH array.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* NTPIN   number of pins within the cell.
* NBIND   first dimension of INDX.
* ITPIN   integer pin descriptor.
* DRAPIN  double pin descriptor.

*Parameters: input/output
* INDX    position index in the geometry structure.
* DCMESH  cell and pin (if point is in pin) meshing.
* POSC    local cell/pin (if point is in pin, cell otherwise)
*         coordinates of the point.
*
*Parameters: output
* MESHP   pin meshes size (if point is in pin).
* NSURP   number of surfaces for the pin (if point is in pin).
* NREGP   number of regions for the pin (if point is in pin).
* PINCEN  cell coordinates of the pin center (if point is in pin).
* INDEX   pin  index vector (if point is in pin).
* IDREG   pin region id array (if point is in pin).
* IDIRP   pin orientation (if point is in pin).
* INPIN   logical flag: true (if point is within a pin).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER IPRINT,NDIM,MAXMSH,MXGSUR,MXGREG,NTPIN,NBIND,
     1 INDX(NBIND,2),ITPIN(3,NTPIN),IDIRP,MESHP(4),NSURP,NREGP,
     2 INDEX(5,-MXGSUR:MXGREG),IDREG(NREGP)
      DOUBLE PRECISION DCMESH(-1:MAXMSH,4,2),DRAPIN(-1:4,NTPIN),
     1 PINCEN(3),POSC(4)
      LOGICAL INPIN
*----
*  LOCAL COORDINATES
*----
      INTEGER IDIR,ITPINO,IPIN,IDG1,IDG2
      DOUBLE PRECISION POSO(3),HALFL,R,PI,XDRCST,ANGP,COSP,SINP
      CHARACTER NAMPIN*9,NAMREC*12
      INTEGER INDOS(2,3)
      CHARACTER RIDNAM*3
      PARAMETER (RIDNAM='RIC')
      DATA INDOS / 2,3,
     1             3,1,
     2             1,2 /
*----
*  SCAN PINS TO FIND IF THE POINT IS LOCATED WITHIN ONE OF THEN 
*----
      PI=XDRCST('Pi',' ')
*     CHANGE COORDINATES TO TAKE INTO ACCOUNT OFFCEN FOR CYLINDERS AND PINS
      DO IDIR=1,NDIM
         POSO(IDIR)=POSC(IDIR)-DCMESH(-1,IDIR,1)
      ENDDO
      ITPINO=INDX(6,1)
*     SCAN PINS
      DO 10 IPIN=1,NTPIN
         IDIRP=ABS(ITPIN(3,IPIN))
         IF (NDIM.EQ.3) THEN
*        VERIFY THE POSITION ALONG THE PIN AXIS
            HALFL=0.5D0*DRAPIN(IDIRP,IPIN)
            IF ((POSO(IDIRP).LT.-HALFL).OR.
     1          (POSO(IDIRP).GT.HALFL)) GOTO 10
         ENDIF
*        VERIFY RADIAL POSITION
*        IDG1 is first direction of plane perpendicular
*                to main direction ($Y, $Z$ or $X$).
*        IDG2 is second direction of plane perpendicular
*                to main direction ($Z$, $X$ or $Y$).
*        IDIRP is main cylinder direction ($X$, $Y$ or $Z$)
*                for 2D case IDIRP=3
         IDG1=INDOS(1,IDIRP)
         IDG2=INDOS(2,IDIRP)
*        pin center         
         PINCEN(IDG1)=DRAPIN(0,IPIN)*COS(DRAPIN(-1,IPIN))
         PINCEN(IDG2)=DRAPIN(0,IPIN)*SIN(DRAPIN(-1,IPIN))
         PINCEN(IDIRP)=0.D0
*        distance with respect to this center
         R=((POSO(IDG1)-PINCEN(IDG1))**2+(POSO(IDG2)-PINCEN(IDG2))**2)
         R=SQRT(R)
         IF (R.LT.DRAPIN(4,IPIN)) THEN
            INDX(5,1)=IPIN
            INDX(6,1)=ITPIN(2,IPIN)
            GOTO 20
         ENDIF
 10   CONTINUE
*     NO
      INPIN=.FALSE.
      RETURN
*     YES
 20   CONTINUE
      INPIN=.TRUE.
      POSC(4)=R
      IF (IPRINT.GT.4) THEN
         WRITE(6,*) 'PIN: TYPE:',INDX(6,1),' INDEX IN CELL:',INDX(5,1)
      ENDIF
      IF (INDX(6,1).EQ.ITPINO) THEN
         IF (INDX(7,2).LE.0) THEN
            WRITE(NAMPIN,'(A1,I8.8)') 'P',ITPINO
            NAMREC=NAMPIN//RIDNAM
            CALL LCMGET(IPTRK,NAMREC,IDREG)
         ENDIF
      ELSE
         CALL MCTLDP(IPTRK,IPRINT,MAXMSH,MXGSUR,MXGREG,INDX(6,1),
     1        RIDNAM,MESHP,NSURP,NREGP,DCMESH(-1,1,2),INDEX,IDREG)
      ENDIF
*     Change coordinates to (origin=pin center)
      DO IDIR=1,NDIM
         POSO(IDIR)=POSO(IDIR)-PINCEN(IDIR)
      ENDDO
*     Rotate geometry by (Pi/2-alpha) for pin at alpha.
      ANGP=0.5D0*PI-DRAPIN(-1,INDX(5,1))
      COSP=COS(ANGP)
      SINP=SIN(ANGP)
      POSC(IDG1)=POSO(IDG1)*COSP-POSO(IDG2)*SINP
      POSC(IDG2)=POSO(IDG1)*SINP+POSO(IDG2)*COSP
      POSC(IDIRP)=POSO(IDIRP)
*
      RETURN
      END
