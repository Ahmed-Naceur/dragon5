*DECK EDIMRG
      SUBROUTINE EDIMRG(IPTRK ,IPMRG ,IPRINT,GEONAM,ITM   ,NREGIO,
     >                  NMERGE,IMERGE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find merge vector by mixtures.
*
*Copyright:
* Copyright (C) 2001 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* IPMRG   merge geometry (ITM=-1) or 
*         tracking (ITM=1) data structure.
* IPRINT  print level.                     
* GEONAM  geometry name.
* ITM     type of merge:
*         =-1 merge from second geometry;
*         = 0 merge from calculation tracking;
*         = 1 merge from second tracking.
* NREGIO  number of regions.                     
*
*Parameters: output
* NMERGE  final number of merged regions.
* IMERGE  merged region index.
*
*-----------------------------------------------------------------------
*
*
*---------------------------  EDIMRG  ---------------------------------
*
*  1- PROGRAMME STATISTICS:
*      NAME     : EDIMRG
*      USE      : FIND MERGE VECTOR BY MIXTURES IN IPMRG
*                 WITH IPTRK
*      MODIFIED : 2001/10/30 (G.M)
*      AUTHOR   : G.MARLEAU
*
*  2- ROUTINE PARAMETERS: 
*      IPTRK    : CALCULATION TRACKING DATA STRUCTURE    
*                 ***> INTEGER IPTRKI
*      IPMRG    : MERGE GEOMETRY (ITM=-1) OR 
*                 TRACKING (ITM=1) DATA STRUCTURE        I
*                 ***> INTEGER IPMRG
*      IPRINT   : PRINT LEVEL                     
*                 ***> INTEGER IPRINT
*     GEONAM : GEOMETRY NAME
*              ***> CHARACTER*12 GEONAM
*      ITM      : TYPE OF MERGE                          I
*                 ITM = -1 -> FROM SECOND GEOMETRY
*                 ITM = 0  -> FROM CALCULATION TRACKING
*                 ITM = 1  -> FROM SECOND TRACKING
*                 ***> INTEGER ITM
*      NREGIO   : NUMBER OF REGIONS                     
*                 ***> INTEGER NREGIO
*      NMERGE   : FINAL NUMBER OF MERGED REGIONS         
*                 ***> INTEGER NMERGE
*      IMERGE   : MERGED REGIONS POSITION                
*                 ***> INTEGER IMERGE(NREGIO)
*
*---------------------------   EDIMRG  --------------------------------
*
      USE         GANLIB
      IMPLICIT    NONE
      INTEGER     IOUT,NSTATE
      CHARACTER   NAMSBR*6
      PARAMETER  (IOUT=6,NSTATE=40,NAMSBR='EDIMRG')
*----
*  ROUTINE PARAMETERS
*----
      TYPE(C_PTR) IPTRK
      TYPE(C_PTR) IPMRG
      INTEGER     IPRINT
      CHARACTER   GEONAM*12
      INTEGER     ITM
      INTEGER     NREGIO
      INTEGER     NMERGE
      INTEGER     IMERGE(NREGIO)
*----
*  LOCAL PARAMETERS
*----
      INTEGER     ISTATE(NSTATE)
      TYPE(C_PTR) IPTRK2
      INTEGER     IMODT2,IMEDT2,ICLST2,IPRIN2
      INTEGER     ITYPEG,ITGEO
      CHARACTER   NAMTR2*12
      CHARACTER   HSIGN*12
      INTEGER     NV,NS,NSOUT,NREG,NUNK,ICODE(6)
      REAL        EXTKOP(NSTATE)
      INTEGER     ITROP,MAXMIX,IREG,ISYMM
      INTEGER     IUEXP,KDROPN,KDRCLS,IRC
      INTEGER     ITYPM
      LOGICAL     LASS,LDRASS
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYMRG,MATALB,MATMRG,MERT
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLSUR,VOLMRG
*----
*  LOCAL NAME OF TEMPORARY TRACKING FILE
*  WHEN IPGEO IS PROVIDED
*----
      NAMTR2='EDIMRGIPTRK2'
      IMODT2=0
      IMEDT2=1
      IPRIN2=0
      ICLST2=2
      ITYPM=1
      ITROP=0
      IF(ITM .EQ. -1) THEN
        LASS=LDRASS(IPMRG,IPRINT)
        CALL LCMOP(IPTRK2,NAMTR2,IMODT2,IMEDT2,IPRIN2)
        HSIGN='L_TRACK     '
        CALL LCMPTC(IPTRK2,'SIGNATURE',12,1,HSIGN)
        HSIGN='EXCELL      '
        CALL LCMPTC(IPTRK2,'TRACK-TYPE',12,1,HSIGN)
*----
*  ANALYZE GEOMETRY
*----
        CALL XDISET(ISTATE,NSTATE,0)
        CALL LCMGET(IPMRG,'STATE-VECTOR',ISTATE)
        ITYPEG= ISTATE(1)
        IF(ITYPEG .EQ. 3 .OR. ITYPEG .EQ. 6 ) THEN
          ITGEO= 1
        ELSE IF(ITYPEG .EQ.  8 .OR. ITYPEG .EQ.  9 .OR.
     >          ITYPEG .EQ. 24 .OR. ITYPEG .EQ. 25      ) THEN
          ITGEO= 2
        ELSE IF(ITYPEG .EQ.  5 .OR. ITYPEG .EQ.  7 .OR.
     >          ITYPEG .EQ. 20 .OR. ITYPEG .EQ. 21 .OR.
     >          ITYPEG .EQ. 22 .OR. ITYPEG .EQ. 23      ) THEN
          ITGEO= 3
        ELSE
          ITGEO= 0
        ENDIF
        IF(ISTATE(13) .GE. 1) THEN
*----
*  CLUSTER GEOMETRY
*----
          ISYMM=1
*c          CALL AXGXCW(IPMRG ,IPTRK2,IPRINT,GEONAM,ISYMM )
          CALL XABORT('EDIMRG: NOT IMPLEMENTED(1):'//GEONAM)
          ITROP=3
        ELSE IF(ITGEO .EQ. 2 ) THEN
*----
*  HEXAGONAL 2D GEOMETRIES
*----
*          CALL AXGXHX(IPMRG ,IPTRK2,IPRINT,GEONAM)
          ITROP=2
        ELSE IF(ITGEO .EQ. 3 ) THEN
*----
*  CARTESIAN 2D/3D ASSEMBLIES
*  CALL XELPRP TO GET GEOMETRY DIMENSIONING INFORMATION
*----
*c          CALL AXGXEL(IPMRG ,IPTRK2,IPRINT,GEONAM)
          CALL XABORT('EDIMRG: NOT IMPLEMENTED(2):'//GEONAM)
          ITROP=1
        ELSE
          CALL XABORT(NAMSBR//': INVALID TYPE OF GEOMETRY')
        ENDIF
        CALL LCMGET(IPTRK2,'ICODE       ',ICODE)
        CALL LCMSIX(IPTRK2,'EXCELL      ',1)
        CALL XDISET(ISTATE,NSTATE,0)
        CALL LCMGET(IPTRK2,'STATE-VECTOR',ISTATE)
        NV=ISTATE(3)
        NS=ISTATE(2)
        NUNK=NV+NS+1
        ALLOCATE(KEYMRG(NUNK),MATALB(NUNK),VOLSUR(NUNK))
        CALL LCMGET(IPTRK2,'KEYMRG      ',KEYMRG)
        CALL LCMGET(IPTRK2,'MATALB      ',MATALB)
        CALL LCMGET(IPTRK2,'VOLSUR      ',VOLSUR)
        CALL LCMSIX(IPTRK2,'EXCELL      ',2)
        ALLOCATE(MATMRG(NUNK),VOLMRG(NUNK))
        CALL XELCMP(NS    ,NV    ,VOLSUR,MATALB,KEYMRG,NSOUT ,
     >              NREG  ,VOLMRG,MATMRG,ITGEO ,ICODE )
        MAXMIX=0
        DO 100 IREG=1,NREG
          KEYMRG(IREG+NSOUT+1)= IREG
          MAXMIX=MAX(MAXMIX,MATMRG(IREG+NSOUT+1))
 100    CONTINUE
        CALL LCMPUT(IPTRK2,'MATCOD',NREG,1,MATMRG(NSOUT+2))
        CALL LCMPUT(IPTRK2,'VOLUME',NREG,2,VOLMRG(NSOUT+2))
        CALL LCMPUT(IPTRK2,'KEYFLX',NREG,1,KEYMRG(NSOUT+2))
        CALL XDRSET(EXTKOP,NSTATE,0.0)
        CALL LCMPUT(IPTRK2,'EXCELTRACKOP',NSTATE,2,EXTKOP)
        CALL XDISET(ISTATE,NSTATE,0)
        ISTATE(1)=NREG
        ISTATE(2)=NREG
        ISTATE(4)=MAXMIX
        ISTATE(5)=NSOUT
        ISTATE(7)=ITROP
        ISTATE(8)=-1
        CALL LCMPUT(IPTRK2,'STATE-VECTOR',NSTATE,1,ISTATE)
        DEALLOCATE(VOLSUR,MATALB,KEYMRG)
        DEALLOCATE(VOLMRG,MATMRG)
*----
*  IF IPRINT >= 10
*  EXPORT TEMPORARY TRACKING FILE
*----
        IF(IPRINT .GE. 10) THEN
          IUEXP=KDROPN('EDIMRGEXP',0,3,0)
          CALL LCMEXP(IPTRK2,IPRINT,IUEXP,2,1)
          IRC=KDRCLS(IUEXP,1)
        ENDIF
      ELSE
        IPTRK2=IPMRG
      ENDIF
*----
*  DESTROY TEMPORARY TRACKING FILE
*  WHEN IPGEO IS PROVIDED
*----
      ALLOCATE(MERT(NREGIO+1))
      CALL XDISET(MERT,NREGIO+1,1)
*C      CALL MRGTRK(IPTRK ,IPTRK2,IPRINT,ITYPM ,NREGIO, MERT)
        CALL LCMLIB(IPTRK)
      CALL XABORT('EDIMRG: NOT IMPLEMENTED(3)')
      NMERGE=0
      DO 110 IREG=1,NREGIO
        IMERGE(IREG)=MERT(IREG+1)
        NMERGE=MAX(NMERGE,IMERGE(IREG))
 110  CONTINUE
      DEALLOCATE(MERT)
      IF(ITM .EQ. -1) THEN
        CALL LCMCL(IPTRK2,ICLST2)
      ENDIF
      RETURN
      END
