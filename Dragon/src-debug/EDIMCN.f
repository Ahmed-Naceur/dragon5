*DECK EDIMCN
      SUBROUTINE EDIMCN(IPTRK ,IPRINT,NDIM  ,NUCELL,NBUCEL,MAXREG,
     >                  NFREG ,NFSUR ,NNC   ,NREGIO,NMERGE,IMERGE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read NXT geometry and generate merging index.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
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
* IPRINT  print level.
* NDIM    dimension of the problem.
* NUCELL  number of cell after unfolding in
*         $X$, $Y$ and $Z$ directions.
* NBUCEL  number of cells in unfolded geometry.
* MAXREG  maximum number of region for any geometry.
* NFREG   final number of regions.
* NFSUR   final number of surfaces.
* NNC     number of saved cells.
* NREGIO  number of regions.
*
*Parameters: output
* NMERGE  final number of merged regions.
* IMERGE  merged region index.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IPRINT,NDIM,NUCELL(3),NBUCEL,MAXREG,NFREG,NFSUR,
     >                 NNC,NREGIO,NMERGE,IMERGE(NREGIO)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='EDIMCN')
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYMRG,ICMRG,IDREG
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ICIS,IUNFLD
*----
*  Local variables
*----
      INTEGER          IEDIMC(NSTATE)
      CHARACTER        NAMREC*12,NAMCEL*9,NAMPIN*9
      INTEGER          NX,NY,NZ,NXY,IX,IY,IZ,ICELL,ICELT,ITRN,ILEV,
     >                 NREGC,IFPIN,ILPIN,IR,IREG,IREGM,IPIN,NBRP,KCIS,
     >                 ICS,ITYLCM
*----
*  Data
*----
      CHARACTER        CLEV(2)*1
      SAVE             CLEV
      DATA             CLEV /'C','P'/
*----
*  Scratch storage allocation
*   KEYMRG  merge region array
*   ICIS    internal cell symmetry
*   IUNFLD  description of unfolded geometry
*   IDREG   region identification array
*   ICMRG   cell material array
*----
      ALLOCATE(KEYMRG(-NFSUR:NFREG),ICMRG(NBUCEL),ICIS(4,NNC),
     > IUNFLD(2,NBUCEL),IDREG(MAXREG))
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Initialise some arrays
*----
      CALL LCMGET(IPTRK,'KEYMRG      ',KEYMRG)
      CALL XDISET(ICMRG,NBUCEL,0)
      CALL XDISET(IMERGE,NREGIO,0)
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,*) 'MAXREG =',MAXREG
        WRITE(IOUT,*) 'KEYMRG=',-NFSUR,NFREG
        WRITE(IOUT,'(17I6)') (KEYMRG(IR),IR=-NFSUR,NFREG)
      ENDIF
*----
*  Read global mesh for geometry
*  and determine graphics size
*----
      CALL LCMGET(IPTRK,'G00000001CIS',ICIS)
      CALL LCMGET(IPTRK,'G00000001CUF',IUNFLD)
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,*) 'IUNFLD=',NBUCEL
        WRITE(IOUT,'(2I6)') (IUNFLD(1,IR),IUNFLD(2,IR),IR=1,NBUCEL)
      ENDIF
      NX=NUCELL(1)
      NY=NUCELL(2)
      NZ=MAX(NUCELL(3),1)
      NXY=NX*NY
      NMERGE=0
*----
*  Scan over $Z$ directions
*----
      DO IZ=1,NZ
*----
*  Scan over $Y$ directions
*----
        DO IY=1,NY
*----
*  Scan over $X$ directions
*----
          DO IX=1,NX
            ICELL=NXY*(IZ-1)+NX*(IY-1)+IX
            ICELT=IUNFLD(1,ICELL)
            ITRN=IUNFLD(2,ICELL)
*----
*  If cell not already merged create new merged mixture
*  and associate cell regions to this mixture
*----
            IF(IPRINT .GE. 100) THEN
              WRITE(IOUT,'(A6,6(1X,I8))') 'CELL  ',
     >        IX,IY,IZ,ICELL,ICELT,ITRN
            ENDIF
            IF(ITRN .EQ.1) THEN
              IF(ICMRG(ICELT) .NE. 0) GO TO 100
              NMERGE=NMERGE+1
              ICMRG(ICELT)=NMERGE
*----
*  Read cell info
*----
              ILEV=1
              WRITE(NAMCEL,'(A1,I8.8)') CLEV(ILEV),ICELT
              NAMREC=NAMCEL//'DIM'
              CALL XDISET(IEDIMC,NSTATE,0)
              CALL LCMGET(IPTRK,NAMREC,IEDIMC)
              NREGC=IEDIMC(8)
              IF(NREGC .GT. MAXREG) CALL XABORT(NAMSBR//': MAXREG for '
     >        //'main geometry not coherent with NREGC for cells')
              IFPIN=IEDIMC(17)
              ILPIN=IFPIN+IEDIMC(16)-1
              NAMREC=NAMCEL//'RID'
              CALL LCMGET(IPTRK,NAMREC,IDREG)
              IF(IPRINT .GE. 100) THEN
                WRITE(IOUT,*) NAMREC//'=',NREGC,IFPIN,ILPIN
                WRITE(IOUT,'(17I6)') (IDREG(IR),IR=1,NREGC)
              ENDIF
              KCIS=0
              DO ICS=1,4
                IF(ICIS(ICS,ICELT) .NE. 0) KCIS=1 
              ENDDO
              DO IR=1,NREGC
                IREG=IDREG(IR)
                IF(IREG .GT. 0) THEN
                  IREGM=KEYMRG(IREG)
                  IF(IMERGE(IREGM) .EQ. 0) THEN
                    IMERGE(IREGM)=NMERGE
                  ELSE IF(IMERGE(IREGM) .NE. NMERGE) THEN
                    WRITE(IOUT,9000) NAMSBR,ICELL,ICELT,
     >                               IREG,IREGM,IMERGE(IREGM)
                    CALL XABORT(NAMSBR//
     >              ': Problem in cells for merge by cell') 
                  ENDIF
                ELSE IF(IREG .LT. 0) THEN
                  IF(KCIS .NE. 1) THEN
                    WRITE(IOUT,9002) NAMSBR,ICELL,ICELT,IREG,IREGM
                    CALL XABORT(NAMSBR//
     >            ': Negative region number for cell without symmetry')
                  ENDIF                    
                ENDIF
              ENDDO
*----
*  Read pin info
*----
              ILEV=2
              DO IPIN=IFPIN,ILPIN
                WRITE(NAMPIN,'(A1,I8.8)') CLEV(ILEV),IPIN
                NAMREC=NAMPIN//'RID'
                CALL LCMLEN(IPTRK,NAMREC,NBRP,ITYLCM)
                IF(NBRP .GT. MAXREG) CALL XABORT(NAMSBR//': MAXREG for'
     >          //' main geometry not coherent with NBRP for pins')
                CALL LCMGET(IPTRK,NAMREC,IDREG)
                DO IR=1,NBRP
                  IREG=ABS(IDREG(IR))
                  IF(IREG .NE. 0) THEN
                    IREGM=KEYMRG(IREG)
                    IF(IMERGE(IREGM) .EQ. 0) THEN
                      IMERGE(IREGM)=NMERGE
                    ELSE IF(IMERGE(IREGM) .NE. NMERGE) THEN
                      WRITE(IOUT,9001) NAMSBR,IPIN,ICELL,ICELT,
     >                               IREG,IREGM,IMERGE(IREGM)
                      CALL XABORT(NAMSBR//
     >                ': Problem in pins for merge by cell')
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
 100          CONTINUE
            ENDIF
          ENDDO
        ENDDO
      ENDDO
*----
*  Verify if all cells analysed
*----
      DO ICELL=1,NX*NY*NZ
        ICELT=IUNFLD(1,ICELL)
        IF(ICMRG(ICELT) .EQ. 0) THEN
          WRITE(IOUT,*) 'Merge Error',ICELL,ICELT
          CALL XABORT(NAMSBR//': Some cells not merged')
        ENDIF
      ENDDO
*----
*  print routine closing  header if required
*----
      IF(IPRINT .GE. 10) THEN
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6010)
          DO IZ=1,NZ
*----
*  Scan over $Y$ directions
*----
            IF(NDIM .EQ. 3) THEN
              WRITE(IOUT,6011) IZ
            ENDIF
            WRITE(IOUT,6012) (IX,IX=1,NX)
            WRITE(IOUT,6013) ('------',IX=1,NX)
            DO IY=NY,1,-1
*----
*  Scan over $X$ directions
*----
              WRITE(IOUT,6014) IY,(ICMRG(IUNFLD(1,ICELL)),
     >        ICELL=NXY*(IZ-1)+NX*(IY-1)+1,NXY*(IZ-1)+NX*IY)
            ENDDO
          ENDDO
          WRITE(IOUT,6020)
          WRITE(IOUT,6021) (IMERGE(IREGM),IREGM=1,NREGIO)
        ENDIF
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(IDREG,IUNFLD,ICIS,ICMRG,KEYMRG)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('Material homogenisation indices for cells'//
     >       ' --  Unfolded geometry')
 6011 FORMAT('Plan Z =',5x,I5)
 6012 FORMAT('    Y | X=',100(1X,I5))
 6013 FORMAT('-----------',100(A6))
 6014 FORMAT(I6,'    |',100(1X,I5))
 6020 FORMAT('Merging Index :')
 6021 FORMAT(12(1X,I5))
 9000 FORMAT(' Error in ',A6,' virtual cell ',I5,
     >       ' (real cell=',I5,') analysis'/3I10)
 9001 FORMAT(' Error in ',A6,' pin ',I5,' virtual cell ',I5,
     >       ' (real cell=',I5,') analysis'/3I10)
 9002 FORMAT(' Internal symmetries problem in ',A6,' virtual cell ',I5,
     >       ' (real cell=',I5,') analysis'/3I10)
      END
