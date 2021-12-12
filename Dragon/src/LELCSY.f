*DECK LELCSY
      FUNCTION LELCSY(IPGEOM,IPRT   ,GEONAM,GEOCV ,NXYZ  ,IGT   ,
     >                GMESH ,ISPLG ,TMESH ,ISPLT )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Verify if geometry satisfies an intrinsic symmetry.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal.
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and G. Marleau
*
*Parameters: input
* IPGEOM  pointer to the reference geometry data structure.
* IPRT    intermediate printing level for output.
* GEONAM  name of the reference geometry.
* GEOCV   name of the geometry to analyze if different from 
*         reference geometry.
* NXYZ    maximum mesh size in directions $x$, $y$ and $z$.
* IGT     geometry turn number.
*
*Parameters: scratch
* GMESH   general mesh description for geometry to analyze.
* ISPLG   general split desctiption for geometry to analyze.
* TMESH   temporary mesh description for geometry comparison.
* ISPLT   temporary split desctiption for geometry comparison.
* LELCSY  result of geometry testing (true or false).
*
*-----------------------------------------------------------------------
*
      USE                GANLIB
      IMPLICIT           NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL            LELCSY
      TYPE(C_PTR)        IPGEOM
      INTEGER            IPRT,NXYZ
      CHARACTER*12       GEONAM,GEOCV 
      INTEGER            IGT(2)
      INTEGER            ISPLG(0:NXYZ-1)
      INTEGER            ISPLT(3,3,0:NXYZ-1)
      REAL               GMESH(0:NXYZ)
      DOUBLE PRECISION   TMESH(3,3,0:NXYZ)
*----
*  LOCAL VARIABLES
*----
      INTEGER            IOUT,NSTATE
      CHARACTER          NAMSBR*6
      PARAMETER         (IOUT=6,NSTATE=40,NAMSBR='LELCSY')
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, TARGET, DIMENSION(:) :: IMIX,IMIX1,IM1,IM2
      INTEGER, POINTER, DIMENSION(:) :: IMIX2
*----
*  LOCAL PARAMETERS
*----
      INTEGER            ISTATE(NSTATE),IKT(2)
      INTEGER            NR,NX,NY,NZ,NK,
     >                   NTM(2,3)
      REAL               OFFCEN(3),OFFTR(3,2)
      INTEGER            ILCMLN,ILCMTY
      INTEGER            IG,IDIR,ITMI,IDMI,
     >                   IR,IX,IY,IZ
      DOUBLE PRECISION   DDM(3),DDD 
      INTEGER            ITYPG,IMTN,IMTO
      LELCSY=.TRUE. 
      IF(IPRT .GE. 10) THEN 
        IF(GEOCV .EQ. ' ') THEN
          WRITE(IOUT,6000) NAMSBR,GEONAM,IGT(1),IGT(2)
        ELSE
          WRITE(IOUT,6000) NAMSBR,GEOCV,IGT(1),IGT(2)
        ENDIF
      ENDIF
*----
*  READ GEOMETRY INFORMATION
*----
      IF(GEOCV .NE. '            ') THEN  
        CALL LCMSIX(IPGEOM,GEOCV,1) 
      ENDIF
*----
*  STATE-VECTOR
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      ITYPG= ISTATE(1)
      NR=ISTATE(2)
      NX=ISTATE(3)
      NY=ISTATE(4)
      NZ=ISTATE(5)
      NK=ISTATE(6)
      CALL XDDSET(TMESH,3*3*(NXYZ+1),0.0D0)
      CALL XDISET(ISPLT,3*3*NXYZ,1)
*----
*  MESHX AND SPLITX
*----
      IDIR=1
      CALL XDRSET(GMESH,(NXYZ+1),0.0)
      CALL XDISET(ISPLG,NXYZ,1)
      CALL LCMLEN(IPGEOM,'MESHX',ILCMLN,ILCMTY)
      IF(ILCMLN .EQ. NX+1)
     >  CALL LCMGET(IPGEOM,'MESHX',GMESH)
      CALL LCMLEN(IPGEOM,'SPLITX',ILCMLN,ILCMTY)
      IF(ILCMLN .EQ. NX)
     >  CALL LCMGET(IPGEOM,'SPLITX',ISPLG)
      DDM(IDIR)=0.0D0 
      DO 10 IX=0,NX-1
        TMESH(3,IDIR,IX)=DBLE(GMESH(IX))
        DDM(IDIR)=MAX(DDM(IDIR),ABS(TMESH(3,IDIR,IX)))
        ISPLT(3,IDIR,IX)=ISPLG(IX)
 10   CONTINUE 
      TMESH(3,IDIR,NX)=DBLE(GMESH(NX))
      DDM(IDIR)=MAX(DDM(IDIR),ABS(TMESH(3,IDIR,NX)))
*----
*  MESHY AND SPLITY
*----
      IDIR=2
      CALL XDRSET(GMESH,(NXYZ+1),0.0)
      CALL XDISET(ISPLG,NXYZ,1)
      CALL LCMLEN(IPGEOM,'MESHY',ILCMLN,ILCMTY)
      IF(ILCMLN .EQ. NY+1)
     >  CALL LCMGET(IPGEOM,'MESHY',GMESH) 
      CALL LCMLEN(IPGEOM,'SPLITY',ILCMLN,ILCMTY)
      IF(ILCMLN .EQ. NY)
     >  CALL LCMGET(IPGEOM,'SPLITY',ISPLG) 
      DDM(IDIR)=0.0D0 
      DO 11 IY=0,NY-1
        TMESH(3,IDIR,IY)=DBLE(GMESH(IY))
        DDM(IDIR)=MAX(DDM(IDIR),ABS(TMESH(3,IDIR,IX)))
        ISPLT(3,IDIR,IY)=ISPLG(IY)
 11   CONTINUE 
      TMESH(3,IDIR,NY)=DBLE(GMESH(NY))
      DDM(IDIR)=MAX(DDM(IDIR),ABS(TMESH(3,IDIR,NX)))
*----
*  MESHZ AND SPLITZ
*----
      IDIR=3
      DDM(IDIR)=0.0D0 
      IF(NZ .GT. 0) THEN
        CALL XDRSET(GMESH,(NXYZ+1),0.0)
        CALL XDISET(ISPLG,NXYZ,1)
        CALL LCMLEN(IPGEOM,'MESHZ',ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. NZ+1) 
     >    CALL LCMGET(IPGEOM,'MESHZ',GMESH)
        CALL LCMLEN(IPGEOM,'SPLITZ',ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. NZ) 
     >    CALL LCMGET(IPGEOM,'SPLITZ',ISPLG)
        DO 12 IZ=0,NZ-1
          TMESH(3,IDIR,IZ)=DBLE(GMESH(IZ))
          DDM(IDIR)=MAX(DDM(IDIR),ABS(TMESH(3,IDIR,IX)))
          ISPLT(3,IDIR,IZ)=ISPLG(IZ)
 12     CONTINUE 
        TMESH(3,IDIR,NZ)=DBLE(GMESH(NZ))
        DDM(IDIR)=MAX(DDM(IDIR),ABS(TMESH(3,IDIR,NX)))
      ELSE
        NZ=1
      ENDIF
C-- MIX
      ALLOCATE(IMIX(NK))    
      CALL XDISET(IMIX,NK,0)
      CALL LCMLEN(IPGEOM,'MIX',ILCMLN,ILCMTY)
      IF(ITYPG .EQ. 21) THEN
C----
C  FOR CARCELX REORDER MIXTURE
C----
        ALLOCATE(IMIX1(NK))    
        IF(ILCMLN .GT. 0 .AND. ILCMLN .LE. NK)
     >    CALL LCMGET(IPGEOM,'MIX',IMIX1)
        IMTN=0
        DO 20 IZ=1,NZ
          DO 21 IY=1,NY
            DO 22 IX=1,NX
              DO 23 IR=1,NR+1
                IMTN=IMTN+1 
                IMTO=(IX-1)*NZ*NY*(NR+1)+(IZ-1)*NY*(NR+1)+
     >               (IY-1)*(NR+1)+IR
                IMIX(IMTN)=IMIX1(IMTO)
 23           CONTINUE
 22         CONTINUE
 21       CONTINUE
 20     CONTINUE
        DEALLOCATE(IMIX1)    
      ELSE IF(ITYPG .EQ. 22) THEN
C----
C  FOR CARCELY REORDER MIXTURE
C----
        ALLOCATE(IMIX1(NK))    
        IF(ILCMLN .GT. 0 .AND. ILCMLN .LE. NK)
     >    CALL LCMGET(IPGEOM,'MIX',IMIX1)
        IMTN=0
        DO 30 IZ=1,NZ
          DO 31 IY=1,NY
            DO 32 IX=1,NX
              DO 33 IR=1,NR+1
                IMTN=IMTN+1 
                IMTO=(IY-1)*NZ*NX*(NR+1)+(IX-1)*NZ*(NR+1)+
     >               (IZ-1)*(NR+1)+IR
                IMIX(IMTN)=IMIX1(IMTO)
 33           CONTINUE
 32         CONTINUE
 31       CONTINUE
 30     CONTINUE
        DEALLOCATE(IMIX1)    
      ELSE
        IF(ILCMLN .GT. 0 .AND. ILCMLN .LE. NK)
     >    CALL LCMGET(IPGEOM,'MIX',IMIX) 
      ENDIF
C- OFF CENTER 
      CALL XDRSET(OFFCEN,3,0.0)
      CALL LCMLEN(IPGEOM,'OFFCENTER',ILCMLN,ILCMTY)
      IF(ILCMLN .EQ. 3)
     >  CALL LCMGET(IPGEOM,'OFFCENTER',OFFCEN)
C----
C  PRINT INITIAL MESH
C----
      IF(IPRT .GE. 20) THEN 
        WRITE(IOUT,6001) NX,NY,NZ,NR,NK
        WRITE(IOUT,6002) 'MESHX ='
        WRITE(IOUT,6003) (TMESH(3,1,IX),IX=0,NX)
        WRITE(IOUT,6002) 'SPLTX ='
        WRITE(IOUT,6004) (ISPLT(3,1,IX),IX=0,NX-1)
        WRITE(IOUT,6002) 'MESHY ='
        WRITE(IOUT,6003) (TMESH(3,2,IY),IY=0,NY)
        WRITE(IOUT,6002) 'SPLTY ='
        WRITE(IOUT,6004) (ISPLT(3,2,IY),IY=0,NY-1)
        WRITE(IOUT,6002) 'MESHZ ='
        WRITE(IOUT,6003) (TMESH(3,3,IZ),IZ=0,NZ)
        WRITE(IOUT,6002) 'SPLTZ ='
        WRITE(IOUT,6004) (ISPLT(3,3,IZ),IZ=0,NZ-1)
        WRITE(IOUT,6002) 'MIXT  ='
        WRITE(IOUT,6004) (IMIX(IX),IX=1,NK)
        WRITE(IOUT,6002) 'OFFC  ='
        WRITE(IOUT,6003) (OFFCEN(IX),IX=1,3)
      ENDIF
C----
C  TURN GEOMETRY WITH IGT
C---- 
      DO 1000 IG=1,2
        IF(IGT(IG) .GT. 12 ) THEN
          IKT(IG)=12-IGT(IG)
        ELSE
          IKT(IG)=IGT(IG)
        ENDIF
        IF(IG.EQ.1) THEN
          ALLOCATE(IM1(2*NK))
          IMIX2=>IM1
        ELSE IF(IG.EQ.2) THEN
          ALLOCATE(IM2(2*NK))
          IMIX2=>IM2
        ENDIF
        CALL XDISET(IMIX2,NK,0)
        IF(IKT(IG) .LT. 0) THEN 
          OFFTR(3,IG)=-OFFCEN(3)
        ELSE
          OFFTR(3,IG)=OFFCEN(3)
        ENDIF
        IF     (ABS(IKT(IG)) .EQ. 1) THEN
          NTM(IG,1)=NX
          NTM(IG,2)=NY
          NTM(IG,3)=NZ
          DO 100 IX=0,NX-1 
            TMESH(IG,1,IX)=TMESH(3,1,IX+1)-TMESH(3,1,IX)
            ISPLT(IG,1,IX)=ISPLT(3,1,IX)
 100      CONTINUE
          DO 110 IY=0,NY-1
            TMESH(IG,2,IY)=TMESH(3,2,IY+1)-TMESH(3,2,IY)
            ISPLT(IG,2,IY)=ISPLT(3,2,IY)
 110      CONTINUE
          OFFTR(1,IG)=OFFCEN(1)
          OFFTR(2,IG)=OFFCEN(2)
          IF(IKT(IG) .LT. 0) THEN
            DO 120 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,NZ-IZ)-TMESH(3,3,NZ-IZ-1)
              ISPLT(IG,3,IZ)=ISPLT(3,3,NZ-IZ-1)
              ITMI=IZ*NX*NY*(NR+1)
              IDMI=(NZ-IZ-1)*NX*NY*(NR+1)
              DO 121 IY=0,NY-1
                DO 122 IX=0,NX-1
                  DO 123 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 123              CONTINUE
 122            CONTINUE
 121          CONTINUE
 120        CONTINUE
          ELSE
            DO 130 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,IZ+1)-TMESH(3,3,IZ)
              ISPLT(IG,3,IZ)=ISPLT(3,3,IZ)
              ITMI=IZ*NX*NY*(NR+1)
              IDMI=ITMI
              DO 131 IY=0,NY-1
                DO 132 IX=0,NX-1
                  DO 133 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 133              CONTINUE
 132            CONTINUE
 131          CONTINUE
 130        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT(IG)) .EQ. 2) THEN 
C----
C  ROTATION OF PI/2 
C----
          NTM(IG,1)=NY
          NTM(IG,2)=NX
          NTM(IG,3)=NZ
          DDD=DDM(2)
          DDM(2)=DDM(1)
          DDM(1)=DDD
          DO 200 IX=0,NY-1 
            TMESH(IG,1,IX)=TMESH(3,2,IX+1)-TMESH(3,2,IX)
            ISPLT(IG,1,IX)=ISPLT(3,2,IX)
 200      CONTINUE
          DO 210 IY=0,NX-1
            TMESH(IG,2,IY)=TMESH(3,1,NX-IY)-TMESH(3,1,NX-IY-1)
            ISPLT(IG,2,IY)=ISPLT(3,1,NX-IY-1)
 210      CONTINUE
          OFFTR(1,IG)=OFFCEN(2)
          OFFTR(2,IG)=-OFFCEN(1)
          IF(IKT(IG) .LT. 0) THEN
            DO 220 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,NZ-IZ)-TMESH(3,3,NZ-IZ-1)
              ISPLT(IG,3,IZ)=ISPLT(3,3,NZ-IZ-1)
              DO 221 IY=0,NX-1
                DO 222 IX=0,NY-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NY*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=(NZ-IZ-1)*NX*NY*(NR+1)+(NY-IX-1)*NX*(NR+1)+
     >                 IY*(NR+1)
                  DO 223 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 223              CONTINUE
 222            CONTINUE
 221          CONTINUE
 220        CONTINUE
          ELSE
            DO 230 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,IZ+1)-TMESH(3,3,IZ)
              ISPLT(IG,3,IZ)=ISPLT(3,3,IZ)
              DO 231 IY=0,NX-1
                DO 232 IX=0,NY-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NY*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=IZ*NX*NY*(NR+1)+IX*NX*(NR+1)+
     >                 (NX-IY-1)*(NR+1)
                  DO 233 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 233              CONTINUE
 232            CONTINUE
 231          CONTINUE
 230        CONTINUE 
          ENDIF
        ELSE IF(ABS(IKT(IG)) .EQ. 3) THEN
C----
C  ROTATION OF PI 
C----
          NTM(IG,1)=NX
          NTM(IG,2)=NY
          NTM(IG,3)=NZ
          DO 300 IX=0,NX-1 
            TMESH(IG,1,IX)=TMESH(3,1,NX-IX)-TMESH(3,1,NX-IX-1)
            ISPLT(IG,1,IX)=ISPLT(3,1,NX-IX-1)
 300      CONTINUE
          DO 310 IY=0,NY-1
            TMESH(IG,2,IY)=TMESH(3,2,NY-IY)-TMESH(3,2,NY-IY-1)
            ISPLT(IG,2,IY)=ISPLT(3,2,NY-IY-1)
 310      CONTINUE
          OFFTR(1,IG)=-OFFCEN(1)
          OFFTR(2,IG)=-OFFCEN(2)
          IF(IKT(IG) .LT. 0) THEN
            DO 320 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,NZ-IZ)-TMESH(3,3,NZ-IZ-1)
              ISPLT(IG,3,IZ)=ISPLT(3,3,NZ-IZ-1)
              DO 321 IY=0,NY-1
                DO 322 IX=0,NX-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NX*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=(NZ-IZ-1)*NX*NY*(NR+1)+(NY-IY-1)*NX*(NR+1)+
     >                 (NX-IX-1)*(NR+1)
                  DO 323 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 323              CONTINUE
 322            CONTINUE
 321          CONTINUE
 320        CONTINUE
          ELSE
            DO 330 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,IZ+1)-TMESH(3,3,IZ)
              ISPLT(IG,3,IZ)=ISPLT(3,3,IZ)
              DO 331 IY=0,NY-1
                DO 332 IX=0,NX-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NX*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=IZ*NX*NY*(NR+1)+(NY-IY-1)*NX*(NR+1)+
     >                 (NX-IX-1)*(NR+1)
                  DO 333 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 333              CONTINUE
 332            CONTINUE
 331          CONTINUE
 330        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT(IG)) .EQ. 4) THEN
C----
C  ROTATION OF 3*PI/2 
C----
          NTM(IG,1)=NY
          NTM(IG,2)=NX
          NTM(IG,3)=NZ
          DDD=DDM(2)
          DDM(2)=DDM(1)
          DDM(1)=DDD
          DO 400 IX=0,NY-1 
            TMESH(IG,1,IX)=TMESH(3,2,NY-IX)-TMESH(3,2,NY-IX-1)
            ISPLT(IG,1,IX)=ISPLT(3,2,NY-IX-1)
 400      CONTINUE
          DO 410 IY=0,NX-1
            TMESH(IG,2,IY)=TMESH(3,1,IY+1)-TMESH(3,1,IY)
            ISPLT(IG,2,IY)=ISPLT(3,1,IY)
 410      CONTINUE
          OFFTR(1,IG)=-OFFCEN(2)
          OFFTR(2,IG)=OFFCEN(1)
          IF(IKT(IG) .LT. 0) THEN
            DO 420 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,NZ-IZ)-TMESH(3,3,NZ-IZ-1)
              ISPLT(IG,3,IZ)=ISPLT(3,3,NZ-IZ-1)
              DO 421 IY=0,NX-1
                DO 422 IX=0,NY-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NY*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=(NZ-IZ-1)*NX*NY*(NR+1)+(NY-IX-1)*NX*(NR+1)+
     >                 IY*(NR+1)
                  DO 423 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 423              CONTINUE
 422            CONTINUE
 421          CONTINUE
 420        CONTINUE
          ELSE
            DO 430 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,IZ+1)-TMESH(3,3,IZ)
              ISPLT(IG,3,IZ)=ISPLT(3,3,IZ)
              DO 431 IY=0,NX-1
                DO 432 IX=0,NY-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NY*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=IZ*NX*NY*(NR+1)+(NY-IX-1)*NX*(NR+1)+
     >                 IY*(NR+1)
                  DO 433 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 433              CONTINUE
 432            CONTINUE
 431          CONTINUE
 430        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT(IG)) .EQ. 5) THEN
C----
C  REFLECTION WITH RESPECT TO AXIS  // TO Y
C----
          NTM(IG,1)=NX
          NTM(IG,2)=NY
          NTM(IG,3)=NZ
          DO 500 IX=0,NX-1 
            TMESH(IG,1,IX)=TMESH(3,1,NX-IX)-TMESH(3,1,NX-IX-1)
            ISPLT(IG,1,IX)=ISPLT(3,1,NX-IX-1)
 500      CONTINUE
          DO 510 IY=0,NY-1
            TMESH(IG,2,IY)=TMESH(3,2,IY+1)-TMESH(3,2,IY)
            ISPLT(IG,2,IY)=ISPLT(3,2,IY)
 510      CONTINUE
          OFFTR(1,IG)=-OFFCEN(1)
          OFFTR(2,IG)=OFFCEN(2)
          IF(IKT(IG) .LT. 0) THEN
            DO 520 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,NZ-IZ)-TMESH(3,3,NZ-IZ-1)
              ISPLT(IG,3,IZ)=ISPLT(3,3,NZ-IZ-1)
              DO 521 IY=0,NY-1
                DO 522 IX=0,NX-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NX*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=(NZ-IZ-1)*NX*NY*(NR+1)+IY*NX*(NR+1)+
     >                 (NX-IX-1)*(NR+1)
                  DO 523 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 523              CONTINUE
 522            CONTINUE
 521          CONTINUE
 520        CONTINUE
          ELSE
            DO 530 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,IZ+1)-TMESH(3,3,IZ)
              ISPLT(IG,3,IZ)=ISPLT(3,3,IZ)
              DO 531 IY=0,NY-1
                DO 532 IX=0,NX-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NX*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=IZ*NX*NY*(NR+1)+IY*NX*(NR+1)+
     >                 (NX-IX-1)*(NR+1)
                  DO 533 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 533              CONTINUE
 532            CONTINUE
 531          CONTINUE
 530        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT(IG)) .EQ. 6) THEN
C----
C  ROTATION OF PI/2 FOLLOWED BY 
C  REFLECTION WITH RESPECT TO AXIS  // TO Y
C----
          NTM(IG,1)=NY
          NTM(IG,2)=NX
          NTM(IG,3)=NZ
          DDD=DDM(2)
          DDM(2)=DDM(1)
          DDM(1)=DDD
          DO 600 IX=0,NY-1 
            TMESH(IG,1,IX)=TMESH(3,2,IX+1)-TMESH(3,2,IX)
            ISPLT(IG,1,IX)=ISPLT(3,2,IX)
 600      CONTINUE
          DO 610 IY=0,NX-1
            TMESH(IG,2,IY)=TMESH(3,1,IY+1)-TMESH(3,1,IY)
            ISPLT(IG,2,IY)=ISPLT(3,1,IY)
 610      CONTINUE
          OFFTR(1,IG)=OFFCEN(2)
          OFFTR(2,IG)=OFFCEN(1)
          IF(IKT(IG) .LT. 0) THEN
            DO 620 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,NZ-IZ)-TMESH(3,3,NZ-IZ-1)
              ISPLT(IG,3,IZ)=ISPLT(3,3,NZ-IZ-1)
              DO 621 IY=0,NX-1
                DO 622 IX=0,NY-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NY*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=(NZ-IZ-1)*NX*NY*(NR+1)+IX*NX*(NR+1)+
     >                 IY*(NR+1)
                  DO 623 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 623              CONTINUE
 622            CONTINUE
 621          CONTINUE
 620        CONTINUE
          ELSE
            DO 630 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,IZ+1)-TMESH(3,3,IZ)
              ISPLT(IG,3,IZ)=ISPLT(3,3,IZ)
              DO 631 IY=0,NX-1
                DO 632 IX=0,NY-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NY*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=IZ*NX*NY*(NR+1)+IX*NX*(NR+1)+
     >                 IY*(NR+1)
                  DO 633 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 633              CONTINUE
 632            CONTINUE
 631          CONTINUE
 630        CONTINUE 
          ENDIF
        ELSE IF(ABS(IKT(IG)) .EQ. 7) THEN
C----
C  REFLECTION WITH RESPECT TO AXIS // TO X
C----
          NTM(IG,1)=NX
          NTM(IG,2)=NY
          NTM(IG,3)=NZ
          DO 700 IX=0,NX-1 
            TMESH(IG,1,IX)=TMESH(3,1,IX+1)-TMESH(3,1,IX)
            ISPLT(IG,1,IX)=ISPLT(3,1,IX)
 700      CONTINUE
          DO 710 IY=0,NY-1
            TMESH(IG,2,IY)=TMESH(3,2,NY-IY)-TMESH(3,2,NY-IY-1)
            ISPLT(IG,2,IY)=ISPLT(3,2,NY-IY-1)
 710      CONTINUE
          OFFTR(1,IG)=OFFCEN(1)
          OFFTR(2,IG)=-OFFCEN(2)
          IF(IKT(IG) .LT. 0) THEN
            DO 720 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,NZ-IZ)-TMESH(3,3,NZ-IZ-1)
              ISPLT(IG,3,IZ)=ISPLT(3,3,NZ-IZ-1)
              DO 721 IY=0,NY-1
                ITMI=IZ*NX*NY*(NR+1)+IY*NX*(NR+1)
                IDMI=(NZ-IZ-1)*NX*NY*(NR+1)+(NY-IY-1)*NX*(NR+1)
                DO 722 IX=0,NX-1
                  DO 723 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 723              CONTINUE
 722            CONTINUE
 721          CONTINUE
 720        CONTINUE
          ELSE
            DO 730 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,IZ+1)-TMESH(3,3,IZ)
              ISPLT(IG,3,IZ)=ISPLT(3,3,IZ)
              DO 731 IY=0,NY-1
                ITMI=IZ*NX*NY*(NR+1)+IY*NX*(NR+1)
                IDMI=IZ*NX*NY*(NR+1)+(NY-IY-1)*NX*(NR+1)
                DO 732 IX=0,NX-1
                  DO 733 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 733              CONTINUE
 732            CONTINUE
 731          CONTINUE
 730        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT(IG)) .EQ. 8) THEN
C----
C  ROTATION OF PI/2 FOLLOWED BY
C  REFLECTION WITH RESPECT TO AXIS // TO X
C----
          NTM(IG,1)=NY
          NTM(IG,2)=NX
          NTM(IG,3)=NZ
          DDD=DDM(2)
          DDM(2)=DDM(1)
          DDM(1)=DDD
          DO 800 IX=0,NY-1 
            TMESH(IG,1,IX)=TMESH(3,2,NY-IX)-TMESH(3,2,NY-IX-1)
            ISPLT(IG,1,IX)=ISPLT(3,2,NY-IX-1)
 800      CONTINUE
          DO 810 IY=0,NX-1
            TMESH(IG,2,IY)=TMESH(3,1,NX-IY)-TMESH(3,1,NX-IY-1)
            ISPLT(IG,2,IY)=ISPLT(3,1,NX-IY-1)
 810      CONTINUE
          OFFTR(1,IG)=-OFFCEN(2)
          OFFTR(2,IG)=-OFFCEN(1)
          IF(IKT(IG) .LT. 0) THEN
            DO 820 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,NZ-IZ)-TMESH(3,3,NZ-IZ-1)
              ISPLT(IG,3,IZ)=ISPLT(3,3,NZ-IZ-1)
              DO 821 IY=0,NX-1
                DO 822 IX=0,NY-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NY*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=(NZ-IZ-1)*NX*NY*(NR+1)+(NY-IX-1)*NX*(NR+1)+
     >                 (NX-IY-1)*(NR+1)
                  DO 823 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 823              CONTINUE
 822            CONTINUE
 821          CONTINUE
 820        CONTINUE
          ELSE
            DO 830 IZ=0,NZ-1
              TMESH(IG,3,IZ)=TMESH(3,3,IZ+1)-TMESH(3,3,IZ)
              ISPLT(IG,3,IZ)=ISPLT(3,3,IZ)
              DO 831 IY=0,NX-1
                DO 832 IX=0,NY-1
                  ITMI=IZ*NX*NY*(NR+1)+IY*NY*(NR+1)+
     >                 IX*(NR+1)
                  IDMI=IZ*NX*NY*(NR+1)+(NY-IX-1)*NX*(NR+1)+
     >                 (NX-IY-1)*(NR+1)
                  DO 833 IR=0,NR
                    IMIX2(ITMI+1)=IMIX(IDMI+1)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 833              CONTINUE
 832            CONTINUE
 831          CONTINUE
 830        CONTINUE 
          ENDIF
        ENDIF
C----
C  PRINT TURNED MESH
C----
        IF(IPRT .GE. 20) THEN 
          WRITE(IOUT,6010) IKT(IG)
          WRITE(IOUT,6002) 'MESHX ='
          WRITE(IOUT,6003) (TMESH(IG,1,IX),IX=0,NX)
          WRITE(IOUT,6002) 'SPLTX ='
          WRITE(IOUT,6004) (ISPLT(IG,1,IX),IX=0,NX-1)
          WRITE(IOUT,6002) 'MESHY ='
          WRITE(IOUT,6003) (TMESH(IG,2,IY),IY=0,NY)
          WRITE(IOUT,6002) 'SPLTY ='
          WRITE(IOUT,6004) (ISPLT(IG,2,IY),IY=0,NY-1)
          WRITE(IOUT,6002) 'MESHZ ='
          WRITE(IOUT,6003) (TMESH(IG,3,IZ),IZ=0,NZ)
          WRITE(IOUT,6002) 'SPLTZ ='
          WRITE(IOUT,6004) (ISPLT(IG,3,IZ),IZ=0,NZ-1)
          WRITE(IOUT,6002) 'MIXT  ='
          WRITE(IOUT,6004) (IMIX2(IX),IX=1,NK)
          WRITE(IOUT,6002) 'OFFC  ='
          WRITE(IOUT,6003) (OFFTR(IX,IG),IX=1,3) 
        ENDIF
 1000 CONTINUE
C----
C  COMPARE GEOMETRY
C  1- MESH AND SPLIT IN X, Y AND Z
C  2- MIXTURES
C  3- OFFCENTER
C----
      DO 900 IDIR=1,3 
        IF(NTM(1,IDIR) .EQ. NTM(2,IDIR)) THEN
          DO 910 IX=0,NTM(1,IDIR)-1
            DDD=ABS(TMESH(2,IDIR,IX)-TMESH(1,IDIR,IX))
            IF(DDD .GT. 1.0D-6*ABS(DDM(IDIR))       .OR.
     >         ISPLT(2,IDIR,IX) .NE. ISPLT(1,IDIR,IX) ) THEN
               WRITE(IOUT,6020) IDIR,IX,
     >           ISPLT(1,IDIR,IX),ISPLT(2,IDIR,IX),
     >           TMESH(1,IDIR,IX),TMESH(2,IDIR,IX),DDD 
              LELCSY=.FALSE. 
              GO TO 995
            ENDIF
 910      CONTINUE
        ELSE
          LELCSY=.FALSE.
          GO TO 995
        ENDIF
 900  CONTINUE
      DO 920 IX=1,NK
        IF(IM1(IX) .NE. IM2(IX) ) THEN
          LELCSY=.FALSE.
          WRITE(IOUT,6021) IX,IM1(IX),IM2(IX)
          GO TO 995
        ENDIF
 920  CONTINUE
      IF(OFFTR(1,1) .NE. OFFTR(1,2) .OR.
     >   OFFTR(2,1) .NE. OFFTR(2,2) .OR.   
     >   OFFTR(3,1) .NE. OFFTR(3,2) ) THEN
        LELCSY=.FALSE.
        GO TO 995
      ENDIF
 995  CONTINUE
C----
C  RELEASE MEMORY
C----
      DEALLOCATE(IM2,IM1,IMIX)
      IF(GEOCV .NE. ' ') THEN  
        CALL LCMSIX(IPGEOM,GEOCV,2) 
      ENDIF
C----
C  RETURN 
C----
      RETURN
C----
C  FORMATS
C----
 6000 FORMAT(1X,A6,'-- ANALYZING :',A12,2I10)
 6001 FORMAT(1X,' DIMENSIONS =',5I10/1X,' ORIGINAL MESH ')
 6002 FORMAT(1X,A7)
 6003 FORMAT(5F15.9)
 6004 FORMAT(5I15) 
 6010 FORMAT(1X,' GEOMETRY AFTER TURN = ',I10)
 6020 FORMAT(1X,'ERROR FOR DIRECTION ',2I10/
     >       1X,'SPLIT = ',2I10/
     >       1X,'MESH  = ',3F15.9)
 6021 FORMAT(1X,'ERROR FOR MIXTURE ',I10,2I10)
      END
