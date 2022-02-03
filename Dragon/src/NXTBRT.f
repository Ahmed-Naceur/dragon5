*DECK NXTBRT
      SUBROUTINE NXTBRT(IPTRK ,IPRINT,NDIM  ,ITYPBC,ISAXIS,NBOCEL,
     >                  MAXMSP,MAXPIN,NFSUR ,MXGSUR,MXGREG,IDFRT ,
     >                  MATRT)
*
*----------
*
*Purpose:
* To built the surface reflection/transmission coupling
* array.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* IPRINT  print level.
* NDIM    problem dimensions.
* ITYPBC  type of boundary conditions where:
*         =0 for geometry with Cartesian boundaries;
*         =1 for geometry with annular boundary;
*         =2 for geometry with hexagonal boundary.
* ISAXIS  symmetry vector for each direction.
* NBOCEL  number of cells in original geometry.
* MAXMSP  maximum number of elements in MESH array.
* MAXPIN  maximum number of pins in clusters.
* IDFRT   identify reflection/transmission faces.
* NFSUR   final number of surfaces.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
*
*Parameters: output
* MATRT   reflection/transmission surface coupling array.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IPRINT,NDIM,ITYPBC,ISAXIS(3),NBOCEL,
     >                 MAXMSP,MAXPIN,IDFRT(8,NBOCEL),
     >                 NFSUR,MXGSUR,MXGREG
      INTEGER          MATRT(NFSUR)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTBRT')
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
*----
*  Local variables
*----
      INTEGER          ISV,IDT,IDIR,ICEL,IGEN(2),ILEV,IG,NBSD,NBST,
     >                 NR1,NS1,NUNK1,IG1,ICL1,IPIN1,IFPIN1,ILPIN1,
     >                 NR2,NS2,NUNK2,IG2,ICL2,IPIN2,IFPIN2,ILPIN2,
     >                 MXRUNK,IDO
      INTEGER          IEDIMC(NSTATE,2),IEDIMP(NSTATE,2)
      CHARACTER        NAMREC*12
      INTEGER          ILCMLN,ILCMTY
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ID1,ID2
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IX1,IX2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SV1,SV2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DAMESH,DRAPIN
*----
*  Data
*----
      CHARACTER        CDIR(1:4)*1,CLEV(2)*1
      SAVE             CDIR,CLEV
      DATA             CDIR /'X','Y','Z','R'/
      DATA             CLEV /'C','P'/
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. -100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      MXRUNK=MXGSUR+MXGREG+1
*----
*  Scratch storage allocation
*   DRAPIN  temporary vector for storing global pin positions.
*   DAMESH  temporary vector for storing global mesh array.
*----
      ALLOCATE(ID1(MXGSUR),ID2(MXGSUR))
      ALLOCATE(IX1(5,MXRUNK),IX2(5,MXRUNK))
      ALLOCATE(SV1(MXRUNK),SV2(MXRUNK))
      ALLOCATE(DAMESH(-1:MAXMSP,4,2),DRAPIN(-1:4,MAXPIN,2))
*----
*  Initialize MATRT assuming all surfaces are reflective
*----
      DO ISV=1,NFSUR
        MATRT(ISV)=ISV
      ENDDO
*----
*  X, Y, and Z translation
*  Scan over cells and locate those with X- surface boundary
*  Find X+ cell from which neutrons are generated
*----
      DO IDT=1,3
        IDO=2*IDT-1
        IF(ISAXIS(IDT) .EQ. 3) THEN
          DO ICEL=1,NBOCEL
            IGEN(1)=ICEL
            ILEV=1
            IGEN(2)=IDFRT(IDO,ICEL)
            IF(IGEN(2) .GT. 0) THEN
*----
*  Cells are identified:
*  Extract dimensioning vectors.
*----
              CALL XDISET(IEDIMC,NSTATE*2,0)
              DO IG=1,2
                WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(IG),'DIM'
                CALL LCMGET(IPTRK,NAMREC,IEDIMC(1,IG))
*----
*  Read meshes
*----
                IF(ITYPBC .EQ. 2) THEN
*----
*  Hexagons
*----
                  IDIR=1
                  WRITE(NAMREC,'(A1,I8.8,A3)')
     >            CLEV(ILEV),IGEN(IG),'SM'//CDIR(IDIR)
                  CALL LCMLEN(IPTRK,NAMREC,ILCMLN,ILCMTY)
                  IF(ILCMLN .GT. 0)
     >            CALL LCMGET(IPTRK,NAMREC,DAMESH(-1,IDIR,IG))
                  IDIR=3
                  WRITE(NAMREC,'(A1,I8.8,A3)')
     >            CLEV(ILEV),IGEN(IG),'SM'//CDIR(IDIR)
                  CALL LCMLEN(IPTRK,NAMREC,ILCMLN,ILCMTY)
                  IF(ILCMLN .GT. 0) THEN
                    CALL LCMGET(IPTRK,NAMREC,DAMESH(-1,IDIR,IG))
                  ELSE
                    CALL XDDSET(DAMESH(-1,IDIR,IG),MAXMSP+2,0.0D0)
                  ENDIF
                ELSE
*----
*  Cartesian, annluar or spherical
*----
                  DO IDIR=1,4
                    WRITE(NAMREC,'(A1,I8.8,A3)')
     >              CLEV(ILEV),IGEN(IG),'SM'//CDIR(IDIR)
                    CALL LCMLEN(IPTRK,NAMREC,ILCMLN,ILCMTY)
                    IF(ILCMLN .GT. 0) THEN
                      CALL LCMGET(IPTRK,NAMREC,DAMESH(-1,IDIR,IG))
                    ELSE
                      CALL XDDSET(DAMESH(-1,IDIR,IG),MAXMSP+2,0.0D0)
                    ENDIF
                  ENDDO
                ENDIF
                WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(IG),'PIN'
                CALL LCMLEN(IPTRK,NAMREC,ILCMLN,ILCMTY)
                IF(ILCMLN .GT. 0) THEN
                  CALL LCMGET(IPTRK,NAMREC,DRAPIN(-1,1,IG))
                ELSE
                  CALL XDDSET(DRAPIN(-1,1,IG),6,0.0D0)
                ENDIF
              ENDDO
*----
*  Find maximum surfaces and regions and retreive
*  MESH, DRAPIN, INDXSR, IDSUR and SURVOL
*----
              NR1=IEDIMC(8,1)
              NS1=IEDIMC(9,1)
              NUNK1=NR1+NS1+1
              WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(1),'VSE'
              CALL LCMGET(IPTRK,NAMREC,SV1)
              WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(1),'VSI'
              CALL LCMGET(IPTRK,NAMREC,IX1)
              WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(1),'SID'
              CALL LCMGET(IPTRK,NAMREC,ID1)
              NR2=IEDIMC(8,1)
              NS2=IEDIMC(9,1)
              NUNK2=NR2+NS2+1
              WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(2),'VSE'
              CALL LCMGET(IPTRK,NAMREC,SV2)
              WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(2),'VSI'
              CALL LCMGET(IPTRK,NAMREC,IX2)
              WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(2),'SID'
              CALL LCMGET(IPTRK,NAMREC,ID2)
*----
*  Find equivalent translated surface
*----
              IF(ITYPBC .EQ. 2) THEN
                CALL NXTETH(IPRINT,IDT   ,ILEV  ,NFSUR ,MAXMSP,
     >                      NS1   ,NR1   ,NS2   ,NR2   ,IEDIMC,DAMESH,
     >                      IX1,ID1,SV1,IX2,ID2 ,SV2,
     >                      MATRT ,NBSD  ,NBST  )
              ELSE
                CALL NXTETS(IPRINT,IDT   ,ILEV  ,NFSUR ,MAXMSP,
     >                      NS1   ,NR1   ,NS2   ,NR2   ,IEDIMC,DAMESH,
     >                      IX1,ID1,SV1,IX2,ID2 ,SV2,
     >                      MATRT ,NBSD  ,NBST  )
              ENDIF
*----
*  For EACH pin in first geometry, find if a pin at an equivalent position
*  in second geometry can be found.
*----
*----
*  Start correction 2010/11/10
*  Pin analysis not required in 3 dimensions
          IF(NDIM .EQ. 3) THEN
*  Start correction 2010/11/10
*----
              ILEV=2
              CALL XDISET(IEDIMP,NSTATE*2,0)
              IG1=1
              IG2=2
              IGEN(IG1)=IEDIMC(17,IG1)-1
              DO ICL1=1,IEDIMC(16,IG1)
                IGEN(IG1)=IGEN(IG1)+1
                IG=IG1
                WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(IG),'DIM'
                CALL LCMGET(IPTRK,NAMREC,IEDIMP(1,IG))
                IFPIN1=IEDIMP(16,IG)
                ILPIN1=IFPIN1+IEDIMP(17,IG)-1
                DO IDIR=1,4
                  WRITE(NAMREC,'(A1,I8.8,A3)')
     >            CLEV(ILEV),IGEN(IG),'SM'//CDIR(IDIR)
                  CALL LCMLEN(IPTRK,NAMREC,ILCMLN,ILCMTY)
                  IF(ILCMLN .GT. 0) THEN
                    CALL LCMGET(IPTRK,NAMREC,DAMESH(-1,IDIR,IG))
                  ELSE
                    CALL XDDSET(DAMESH(-1,IDIR,IG),MAXMSP+2,0.0D0)
                  ENDIF
                ENDDO
                NR1=IEDIMP(8,1)
                NS1=IEDIMP(9,1)
                NUNK1=NR1+NS1+1
                WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(1),'VSE'
                CALL LCMGET(IPTRK,NAMREC,SV1)
                WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(1),'VSI'
                CALL LCMGET(IPTRK,NAMREC,IX1)
                WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(1),'SID'
                CALL LCMGET(IPTRK,NAMREC,ID1)
                IGEN(IG2)=IEDIMC(17,IG2)-1
                DO ICL2=1,IEDIMC(16,IG2)
                  IGEN(IG2)=IGEN(IG2)+1
                  IG=IG2
                  WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(IG),'DIM'
                  CALL LCMGET(IPTRK,NAMREC,IEDIMP(1,IG))
                  IFPIN2=IEDIMP(16,IG)
                  ILPIN2=IFPIN2+IEDIMP(17,IG)-1
                  DO IDIR=1,4
                    WRITE(NAMREC,'(A1,I8.8,A3)')
     >              CLEV(ILEV),IGEN(IG),'SM'//CDIR(IDIR)
                    CALL LCMLEN(IPTRK,NAMREC,ILCMLN,ILCMTY)
                    IF(ILCMLN .GT. 0) THEN
                      CALL LCMGET(IPTRK,NAMREC,DAMESH(-1,IDIR,IG))
                    ELSE
                      CALL XDDSET(DAMESH(-1,IDIR,IG),MAXMSP+2,0.0D0)
                    ENDIF
                  ENDDO
                  NR2=IEDIMP(8,IG2)
                  NS2=IEDIMP(9,IG2)
                  NUNK2=NR2+NS2+1
                  WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(2),'VSE'
                  CALL LCMGET(IPTRK,NAMREC,SV2)
                  WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(2),'VSI'
                  CALL LCMGET(IPTRK,NAMREC,IX2)
                  WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEN(2),'SID'
                  CALL LCMGET(IPTRK,NAMREC,ID2)
*----
*  Find equivalent translated surface
*----
                  CALL NXTETS(IPRINT,IDT   ,ILEV  ,NFSUR ,MAXMSP,
     >                        NS1   ,NR1   ,NS2   ,NR2   ,IEDIMP,DAMESH,
     >                        IX1,ID1,SV1,IX2,ID2 ,SV2,
     >                        MATRT ,NBSD  ,NBST  )
*----
*  No IDT directed face for direct cluster
*  go to next direct cluster
*----
                  IF(NBSD .EQ. 0) GO TO 105
                  IF(NBSD .EQ. NBST) THEN
*----
*  Test if pin position locations are adequate
*----
                    DO IPIN1=IFPIN1,ILPIN1
                      DO IPIN2=IFPIN2,ILPIN2
                        IF(DRAPIN(-1,IPIN1,1) .EQ. DRAPIN(-1,IPIN2,2)
     >               .AND. DRAPIN( 0,IPIN1,1) .EQ. DRAPIN( 0,IPIN2,2)
     >               .AND. DRAPIN( 4,IPIN1,1) .EQ. DRAPIN( 4,IPIN2,2)
     >                     ) THEN
*----
*   Pin positions are identical, select next pin
*----
                          GO TO 125
                        ENDIF
                      ENDDO
*----
*  Pin positions are not compatible
*  go to next translated cluster
*----
                      GO TO 115
 125                  CONTINUE
                    ENDDO
*----
*   Translation surfaces found here go to next direct cluster
*----
                    GO TO 105
                  ENDIF
 115              CONTINUE
                ENDDO
*----
*  Translated surfaces for directed pin not found
*  send warning signal and continue
*----
                WRITE(IOUT,9000) ICEL,ICL1
 105            CONTINUE
              ENDDO
*----
*  Start correction 2010/11/10
*  Pin analysis not required in 3 dimensions
          ENDIF
*  End correction 2010/11/10
*----
            ENDIF
          ENDDO
        ENDIF
      ENDDO
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(DRAPIN,DAMESH)
      DEALLOCATE(SV2,SV1)
      DEALLOCATE(IX2,IX1)
      DEALLOCATE(ID2,ID1)
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 9000 FORMAT(' ***** Warning ***** '/
     >       '       Translated surface for CELL ',I5,1X,
     >       ' and PIN :',I5,' is absent')
      END
