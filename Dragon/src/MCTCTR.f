*DECK MCTCTR
      SUBROUTINE MCTCTR(IPTRK,IPRINT,NDIM,MAXMSH,NUCELL,MXGSUR,MXGREG,
     1                  MAXPIN,IUNFLD,DGMESH,NBIND,ODIR,POS,IREG,INDX,
     2                  IDIRC,MESHC,NSURC,NREGC,NTPIN,CELLPO,PINCEN,
     3                  INDEX,IDREG,DCMESH,ITPIN,DRAPIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find region index from global coordinates of a point within the
* geometry.
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
* NUCELL  number of cell after unfolding in $X$, $Y$ and
*         $Z$ directions.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* MAXPIN  maximum number of pins in a cell.
* IUNFLD  description of unfolded geometry.
* DGMESH  meshing vector for global geometry.
* NBIND   first dimension of INDX.
* ODIR    search (octant) direction. 
* POS     spatial coordinates of the point to locate.
*
*Parameters: input/output
* IREG    region index.
* INDX    position index in the geometry structure.
*
*Parameters: scratch
* IDIRC   undefined.
* MESHC   undefined.
* NSURC   undefined.
* NREGC   undefined.
* NTPIN   undefined.
* CELLPO  undefined.
* PINCEN  undefined.
* INDEX   undefined.
* IDREG   undefined.
* DCMESH  undefined.
* ITPIN   undefined.
* DRAPIN  undefined.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER IPRINT,NDIM,MAXMSH,NUCELL(3),MXGSUR,MXGREG,MAXPIN,
     1 IUNFLD(2,NUCELL(1),NUCELL(2),NUCELL(3)),NBIND,ODIR(3),IREG,
     2 INDX(NBIND,0:2),IDIRC(2),MESHC(4,2),NSURC(2),NREGC(2),NTPIN,
     3 INDEX(5,-MXGSUR:MXGREG,2),IDREG(MXGREG,2),ITPIN(3,MAXPIN)
      DOUBLE PRECISION DGMESH(-1:MAXMSH,4),POS(3),CELLPO(3,2),PINCEN(3),
     1 DCMESH(-1:MAXMSH,4,2),DRAPIN(-1:4,MAXPIN)
*----
*  LOCAL VARIABLES
*----
      INTEGER JJ,ICELO,ICEL,ITRNO,ITRN,ILEV,ODIRC(3),IDG1,IDG2,IDIR
      DOUBLE PRECISION POSC(4)
      CHARACTER NAMCEL*9,NAMREC*12
      LOGICAL INPIN
      INTEGER INDOS(2,3)
      CHARACTER RIDNAM*3
      PARAMETER (RIDNAM='RIC')
      DATA INDOS / 2,3,
     1             3,1,
     2             1,2 /
*----
*  FIND LOCATION WITHIN THE GLOBAL ASSEMBLY (level 0)
*---- 
      ICELO=INDX(5,0)
      ITRNO=INDX(6,0)
      DO IDIR=1,NDIM
         IF (ODIR(IDIR).EQ.1) THEN
            INDX(IDIR,0)=1
            DO WHILE (POS(IDIR).GT.DGMESH(INDX(IDIR,0),IDIR))
               INDX(IDIR,0)=INDX(IDIR,0)+1
            ENDDO
         ELSE     
            INDX(IDIR,0)=NUCELL(IDIR)-1
            DO WHILE (POS(IDIR).LE.DGMESH(INDX(IDIR,0),IDIR))
               INDX(IDIR,0)=INDX(IDIR,0)-1
            ENDDO
            INDX(IDIR,0)=INDX(IDIR,0)+1
         ENDIF
      ENDDO
      DO IDIR=NDIM+1,3
         INDX(IDIR,0)=1
      ENDDO
      ICEL=IUNFLD(1,INDX(1,0),INDX(2,0),INDX(3,0))
      ITRN=IUNFLD(2,INDX(1,0),INDX(2,0),INDX(3,0))
      INDX(5,0)=ICEL
      INDX(6,0)=ITRN
      IF (IPRINT.GT.4) THEN
         WRITE(6,*) '**** FIND REGION WITHIN THE GEOMETRY ****'
         IF (IPRINT.GT.5) THEN
            WRITE(6,*) '0 - GLOBAL MESH:'
            WRITE(6,'(3(I3,1X,1H(,F8.6,1H<,F8.6,1H<,F8.6,1H),2X))')
     1      (INDX(JJ,0),DGMESH(INDX(JJ,0)-1,JJ),POS(JJ),
     2                  DGMESH(INDX(JJ,0),JJ),JJ=1,NDIM)
         ENDIF
         WRITE(6,*) 'CELL:',ICEL,' ROTATION:',ITRN
      ENDIF
*----
*  ENTER THE CORRESPONDING CELL GEOMETRY
*----
      ILEV=1
      IF (ICEL.EQ.ICELO) THEN
         IF (INDX(7,ILEV).LE.0) THEN
            WRITE(NAMCEL,'(A1,I8.8)') 'C',ICEL
            NAMREC=NAMCEL//RIDNAM
            CALL LCMGET(IPTRK,NAMREC,IDREG(1,ILEV))
         ENDIF
         IF (ITRN.NE.ITRNO) THEN
            CELLPO(1,2)=DGMESH(INDX(1,0),1)
            CELLPO(1,1)=DGMESH(INDX(1,0)-1,1)
            CELLPO(2,2)=DGMESH(INDX(2,0),2)
            CELLPO(2,1)=DGMESH(INDX(2,0)-1,2)
            IF(NDIM.EQ.3) THEN
               CELLPO(3,2)=DGMESH(INDX(3,0),3)
               CELLPO(3,1)=DGMESH(INDX(3,0)-1,3) 
            ENDIF      
         ENDIF
      ELSE
      CALL MCTLDC(IPTRK,IPRINT,NDIM,MAXMSH,MXGSUR,MXGREG,NBIND,ICEL,
     1     INDX(1,0),RIDNAM,DGMESH,IDIRC(ILEV),MESHC(1,ILEV),
     2     NSURC(ILEV),NREGC(ILEV),NTPIN,CELLPO,DCMESH(-1,1,ILEV),
     3     INDEX(1,-MXGSUR,ILEV),ITPIN,DRAPIN,IDREG(1,ILEV))
      ENDIF
*----
*  AFTER THE CALL TO MCTCCC THE COORDINATES ARE IN TERMS OF CELL COORDINATES
*----
      CALL MCTCCC(NDIM,ITRN,CELLPO,ODIR,POS,ODIRC,POSC)
*----
*  FIND IF IT IS WITHIN A PIN
*----
      IF (NTPIN.GT.0) THEN
         CALL MCTPIR(IPTRK,IPRINT,NDIM,MAXMSH,MXGSUR,MXGREG,NTPIN,
     1        NBIND,INDX(1,1),ITPIN,DRAPIN,DCMESH,MESHC(1,2),NSURC(2),
     2        NREGC(2),PINCEN,INDEX(1,-MXGSUR,2),IDREG(1,2),POSC,
     3        IDIRC(2),INPIN)
      ELSE
         INPIN=.FALSE.
      ENDIF
*----
*  FIND LOCATION WITHIN CELL OR PIN
*----
      IF (INPIN) ILEV=2
*     Find location in Cartesian mesh
      DO IDIR=1,NDIM
         JJ=1
         DO WHILE (POSC(IDIR).GT.DCMESH(JJ,IDIR,ILEV))
            JJ=JJ+1
         ENDDO
         INDX(IDIR,ILEV)=JJ
      ENDDO
      DO IDIR=NDIM+1,3
         INDX(IDIR,ILEV)=1
      ENDDO
      IF (IDIRC(ILEV).GT.0) THEN
*     Find location in cylindrical mesh
         IF (ILEV.EQ.1) THEN
*           calculate distance to the cell center
*           (already calculated for a pin in MCTPIR) 
            IDG1=INDOS(1,IDIRC(ILEV))
            IDG2=INDOS(2,IDIRC(ILEV))
            POSC(4)=((POSC(IDG1)-DCMESH(-1,IDG1,ILEV))**2
     1              +(POSC(IDG2)-DCMESH(-1,IDG2,ILEV))**2)
            POSC(4)=SQRT(POSC(4))
         ENDIF
         JJ=MESHC(4,ILEV)
         IF (POSC(4).GT.DCMESH(JJ,4,ILEV)) THEN
            INDX(4,ILEV)=0
         ELSE
            JJ=JJ-1
            DO WHILE(POSC(4).LT.DCMESH(JJ,4,ILEV)) 
               JJ=JJ-1
            ENDDO
            INDX(4,ILEV)=JJ+1
         ENDIF
      ELSE
         INDX(4,ILEV)=0
      ENDIF
      IF (IPRINT.GT.5) THEN
         WRITE(6,*) ILEV,'- LOCAL MESH:'
         WRITE(6,'(3(I3,1X,1H(,F8.6,1H<,F8.6,1H<,F8.6,1H),2X))')
     1        (INDX(JJ,ILEV),DCMESH(INDX(JJ,ILEV)-1,JJ,ILEV),POSC(JJ),
     2                       DCMESH(INDX(JJ,ILEV),JJ,ILEV),JJ=1,NDIM)
      ENDIF
*     Find location of this element in the index
      IREG=1
      DO WHILE ((INDX(1,ILEV).NE.INDEX(1,IREG,ILEV)).OR.
     1          (INDX(2,ILEV).NE.INDEX(2,IREG,ILEV)).OR.
     2          (INDX(3,ILEV).NE.INDEX(3,IREG,ILEV)).OR.
     3          (INDX(4,ILEV).NE.INDEX(4,IREG,ILEV)))
         IREG=IREG+1
         IF (IREG.GT.NREGC(ILEV)) THEN
            WRITE(6,*) (INDX(JJ,ILEV),JJ=1,4)
            CALL XABORT('MCTCTR: INDEXES DO NOT MATCH')
         ENDIF
      ENDDO
*
      IREG=ABS(IDREG(IREG,ILEV))
      DO JJ=1,2
         INDX(7,JJ)=0
      ENDDO
      INDX(7,ILEV)=IREG
      IF (IPRINT.GT.4) THEN
         WRITE(6,*) 'REGION:',IREG
         WRITE(6,*) '*** FOUND REGION WITHIN THE GEOMETRY ***' 
      ENDIF
      IF (IREG.LE.0) THEN
         WRITE(6,*) INDX(1,ILEV),INDX(2,ILEV),INDX(3,ILEV),INDX(4,ILEV)
         CALL XABORT('MCTCTR: INVALID REGION FOUND.')
      ENDIF
*
      RETURN
      END
