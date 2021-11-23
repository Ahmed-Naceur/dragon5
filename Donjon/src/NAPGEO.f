*DECK NAPGEO
      SUBROUTINE NAPGEO(IPGNW,IPGOD,IPCPO,NSTATE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Split geometry from homogeneous to heterogeneous assemblies
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal.
*
*Author(s): 
* R. Chambon
*
*Parameters: input/output
* IPGNW   LCM object address of heterogeneous assembly Geometry.
* IPGOD   LCM object address of homogeneous assembly Geometry.
* IPCPO   LCM object address of Multicompo.
* NSTATE  length of the state vector
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NSTATE
      TYPE(C_PTR)  IPGNW,IPGOD,IPCPO
*----
*  LOCAL VARIABLES
*----
      INTEGER IOUT,MAXLIN
      REAL REPS
      PARAMETER (REPS=1.0E-4,IOUT=6,MAXLIN=50)
      TYPE(C_PTR) JPGEO,KPGEO,JPCPO,KPCPO
      INTEGER INDIC,NITMA,LENGTH
      CHARACTER TEXT*12
      REAL FLOT
      DOUBLE PRECISION DFLOT
      INTEGER ISTATE(NSTATE),IMPX,NCODE(6),ICODE(6),ITYPGP,STYPP,
     1  KCHAR(3)
      REAL ZCODE(6)
      INTEGER NXP,NYP,NREGP,NMIXP,NCOMLI,NXD,NYD,NZD,NREGD,NMIXD,NXF,
     1  NYF,NZF,NREGF,NMIXF,NXA,NYA,NMIXA,NXPTMP,NYPTMP,NMIXD2,NASS
      CHARACTER DIRHET*12,COMMEN(MAXLIN)*80,HMSG*131
      INTEGER I,J,K,L,IP,JP,JF,IFBEG,JFBEG,LMIX,IASS,IZ,IZT,JM,JN,IM
      LOGICAL LSPX,LSPY,LSPZ,LPOS,LMGEO
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIXA,NBAX,AZONE,IBAX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MIXP
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: MIXD,MIXF
      REAL, ALLOCATABLE, DIMENSION(:) :: MXP,MYP,MXD,MYD,MZD,MXF,
     1  MYF
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SXP,SYP,SXD,SYD,SZD,
     1  SXF,SYF,AXD,AYD

      IMPX=0
      LSPX=.FALSE.
      LSPY=.FALSE.
      LSPZ=.FALSE.
      LMGEO=.FALSE.
C Read mandatory inputs
      CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
      IF(INDIC.NE.3) CALL XABORT('NAPGEO: character data expected.')
      IF(TEXT.EQ.'EDIT') THEN
        CALL REDGET(INDIC,IMPX,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.1) CALL XABORT('NAPGEO: integer data expected.')
        CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.3) CALL XABORT('NAPGEO: character data expected.')
      ENDIF
      IF(TEXT.NE.'DIRGEO') CALL XABORT('NAPGEO: ''DIRGEO'' '//
     1  'expected.')
      CALL REDGET(INDIC,NITMA,FLOT,DIRHET,DFLOT)
      IF(INDIC.NE.3) CALL XABORT('NAPGEO: character data expected.')
      CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.EQ.'MACGEO') THEN
        LMGEO=.TRUE.
        CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
      ENDIF
      IF(TEXT.NE.'MIXASS') THEN
        CALL XABORT('NAPGEO: ''MIXASS'' expected.')
      ENDIF
      CALL REDGET(INDIC,NMIXA,FLOT,DIRHET,DFLOT)
      IF(INDIC.NE.1)CALL XABORT('@NAPGEO: integer data expected.')
      ALLOCATE(MIXA(NMIXA*2))
      DO I=1,NMIXA
        CALL REDGET(INDIC,MIXA(I),FLOT,DIRHET,DFLOT)
        IF(INDIC.NE.1)CALL XABORT('@NAPGEO: integer data expected.')
      ENDDO
      CALL LCMSIX(IPCPO,DIRHET,1)

*----
*  RECOVER TABLE-OF-CONTENT INFORMATION FOR THE COMPO.
*----
      CALL LCMGET(IPCPO,'STATE-VECTOR',ISTATE)
      NCOMLI=ISTATE(10)
      CALL LCMGTC(IPCPO,'COMMENT',80,NCOMLI,COMMEN)
      IF(IMPX.GT.10)WRITE(IOUT,'(1X,A)') (COMMEN(I),I=1,NCOMLI)
* Get Geometry from calculation #1
* for pin by pin geometries
      IF(LMGEO) THEN
        JPCPO=LCMGID(IPCPO,'MIXTURES')
        KPCPO=LCMGIL(JPCPO,1)
        JPGEO=LCMGID(KPCPO,'CALCULATIONS')
        KPGEO=LCMGIL(JPGEO,1)
        CALL LCMSIX(KPGEO,'MACROLIB    ',1)
        CALL LCMSIX(KPGEO,'GFF         ',1)
        CALL LCMSIX(KPGEO,'GFF-GEOM    ',1)
      ELSE
* for heterogeneous geometries
        JPGEO=LCMGID(IPCPO,'GEOMETRIES')
        KPGEO=LCMGIL(JPGEO,1)
      ENDIF
C get dimension in geometry from L_MULTICOMPO
      if(impx.ge.100)write(6,*) 'debug: get ISTATE multicompo Geometry'
      CALL LCMGET(KPGEO,'STATE-VECTOR',ISTATE)
      ITYPGP=ISTATE(1)
      STYPP=ISTATE(11)
      IF(ITYPGP.NE.5) CALL XABORT('NAPGEO: CAR2D geometry type '
     1 //'expected in L_MULTICOMPO.')
      IF(STYPP.NE.0) CALL XABORT('NAPGEO: No split in geometry '
     1 //'expected in L_MULTICOMPO.')
      NXP=ISTATE(3)
      NYP=ISTATE(4)
      NREGP=ISTATE(6)
      NMIXP=ISTATE(7)
      ALLOCATE(MXP(NXP+1),MYP(NYP+1))
      ALLOCATE(SXP(NXP),SYP(NYP))
      ALLOCATE(MIXP(NXP,NYP))
      CALL LCMGET(KPGEO,'MESHX',MXP)
      CALL LCMGET(KPGEO,'MESHY',MYP)
      CALL LCMGET(KPGEO,'MIX',MIXP)
      CALL XDISET(SXP,NXP,1)
      CALL XDISET(SYP,NYP,1)
C get dimension in homogeneous assembly core Geometry
      if(impx.ge.100)write(6,*) 'debug: ISTATE homog. ass. core Geom.'
      CALL LCMGET(IPGOD,'SIGNATURE',KCHAR)
      CALL LCMGET(IPGOD,'STATE-VECTOR',ISTATE)
      ITYPGP=ISTATE(1)
      NXD=ISTATE(3)
      NYD=ISTATE(4)
      NZD=ISTATE(5)
      NREGD=ISTATE(6)
      NMIXD=ISTATE(7)
      NMIXD2=NMIXD
      ALLOCATE(MXD(NXD+1),MYD(NYD+1))
      ALLOCATE(MZD(NZD+1))
      ALLOCATE(SXD(NXD),SYD(NYD))
      ALLOCATE(SZD(NZD))
      ALLOCATE(AXD(NXD),AYD(NYD))
      ALLOCATE(MIXD(NXD,NYD,NZD))
      CALL LCMGET(IPGOD,'MESHX',MXD)
      CALL LCMGET(IPGOD,'MESHY',MYD)
      CALL LCMGET(IPGOD,'MESHZ',MZD)
      CALL LCMGET(IPGOD,'MIX',MIXD)
      CALL XDISET(AXD,NXD,0)
      CALL XDISET(AYD,NYD,0)
      CALL LCMLEN(IPGOD,'SPLITX',LENGTH,INDIC)
      IF(LENGTH.NE.0) THEN
        CALL LCMGET(IPGOD,'SPLITX',SXD)
        LSPX=.TRUE.
      ELSE
        CALL XDISET(SXD,NXD,1)
      ENDIF
      CALL LCMLEN(IPGOD,'SPLITY',LENGTH,INDIC)
      IF(LENGTH.NE.0) THEN
        CALL LCMGET(IPGOD,'SPLITY',SYD)
        LSPY=.TRUE.
      ELSE
        CALL XDISET(SYD,NYD,1)
      ENDIF
      CALL LCMLEN(IPGOD,'SPLITZ',LENGTH,INDIC)
      IF(LENGTH.NE.0) THEN
        CALL LCMGET(IPGOD,'SPLITZ',SZD)
        LSPZ=.TRUE.
      ELSE
        CALL XDISET(SZD,NZD,1)
      ENDIF
      CALL LCMGET(IPGOD,'NCODE',NCODE)
      CALL LCMGET(IPGOD,'ZCODE',ZCODE)
      CALL LCMGET(IPGOD,'ICODE',ICODE)
C get assembly mixture in homogeneous core geometry
      if(impx.ge.100)write(6,*) 'debug: get assembly mixture'
      DO 40 K=1,NZD
      DO 30 J=1,NYD
      DO 20 I=1,NXD
      DO 10 L=1,NMIXA
         IF(MIXA(L).EQ.MIXD(I,J,K)) THEN
           AXD(I)=1
           AYD(J)=1
           GOTO 20
         ENDIF
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
      if(impx.ge.5) then
        write(6,*) 'Original mesh corresponding to assemblies'
        write(6,*) 'X direction: AXD(1 : NXD)=',(AXD(I),I=1,NXD)
        write(6,*) 'Y direction: AYD(1 : NYD)=',(AYD(I),I=1,NYD)
      endif
C specify splitting in heterogeneous assembly geometry
   50 CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
      IF(INDIC.NE.3) CALL XABORT('@NAPGEO: character data expected 1.')
      IF(TEXT.EQ.'SPLITX-ASS') THEN
        DO I=1,NXP
          CALL REDGET(INDIC,SXP(I),FLOT,TEXT,DFLOT)
          IF(INDIC.NE.1) THEN
            WRITE(HMSG,*) '@NAPGEO: integer number expected'
     1  //' for SPLITX-ASS: ',I,'out of ',NXP
            CALL XABORT(HMSG)
          ENDIF
        ENDDO
        LSPX=.TRUE.
        GOTO 50
      ELSEIF(TEXT.EQ.'SPLITY-ASS') THEN
        DO I=1,NYP
          CALL REDGET(INDIC,SYP(I),FLOT,TEXT,DFLOT)
          IF(INDIC.NE.1) THEN
            WRITE(HMSG,*) '@NAPGEO: integer number expected'
     1  //' for SPLITY-ASS: ',I,'out of ',NYP
            CALL XABORT(HMSG)
          ENDIF
        ENDDO
        LSPY=.TRUE.
        GOTO 50
C read final ';'
      ELSEIF(TEXT.EQ.'MAX-MIX-GEO') THEN
        CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.1) CALL XABORT('@NAPGEO: integer data expected.')
        NMIXD2=MAX(NMIXD,NITMA)
        GOTO 50
C read final ';'
      ELSEIF(TEXT.EQ.';') THEN
        if(impx.ge.5) then
          write(6,*) 'Splitting within assemblies:'
          write(6,*) 'SXP',(SXP(I),I=1,NXP)
          write(6,*) 'SYP',(SYP(I),I=1,NYP)
        endif
        GOTO 60
      ELSE
        CALL XABORT('@NAPGEO: '//TEXT//' WRONG KEYWORD')
      ENDIF
C compute new dimension
C     get number of x and y original mesh to split
   60 NXA=0
      DO I=1,NXD
        NXA=NXA+AXD(I)
      ENDDO
      NYA=0
      DO I=1,NYD
        NYA=NYA+AYD(I)
      ENDDO
C     compute new dimension
      NXF=NXD+NXA*(NXP-1)
      NYF=NYD+NYA*(NYP-1)
      NZF=NZD
      NREGF=NXF*NYF*NZF
C allocate new geometry dimensions
      ALLOCATE(MXF(NXF+1),MYF(NYF+1))
      ALLOCATE(MIXF(NXF,NYF,NZF))
      ALLOCATE(SXF(NXF),SYF(NYF))
C Compute new x/y mesh and new x/y split
      if(impx.ge.100)write(6,*) 'debug: Compute new x/y mesh and split'
      J=1
      MXF(J)=MXD(1)
      DO I=2,NXD+1
        IF(AXD(I-1).EQ.0) THEN
          J=J+1
          MXF(J)=MXD(I)
          SXF(J-1)=SXD(I-1)
        ELSE
          DO K=2,NXP+1
            J=J+1
            MXF(J)=MXD(I-1)+MXP(K)-MXP(1)
            SXF(J-1)=SXP(K-1)
          ENDDO
        ENDIF
      ENDDO
      J=1
      MYF(J)=MYD(1)
      DO I=2,NYD+1
        IF(AYD(I-1).EQ.0) THEN
          J=J+1
          MYF(J)=MYD(I)
          SYF(J-1)=SYD(I-1)
        ELSE
          DO K=2,NYP+1
            J=J+1
            MYF(J)=MYD(I-1)+MYP(K)-MYP(1)
            SYF(J-1)=SYP(K-1)
          ENDDO
        ENDIF
      ENDDO
      IF(MXF(NXF+1).NE.MXD(NXD+1)) CALL XABORT('@NAPGEO: OLD and NEW'
     1 //' X MESH do not match.')
      IF(MYF(NYF+1).NE.MYD(NYD+1)) CALL XABORT('@NAPGEO: OLD and NEW'
     1 //' Y MESH do not match.')
C Compute new mixture
      if(impx.ge.100)write(6,*) 'debug: Compute new mixture'
      NMIXF=NMIXD2+NMIXA*NMIXP
      CALL XDISET(MIXF,NREGF,-1)
      DO 100 L=1,NMIXA
         MIXA(L+NMIXA)=NMIXD2+(L-1)*NMIXP+1 
  100 CONTINUE
      DO 240 K=1,NZD
      JFBEG=0
      DO 230 J=1,NYD
      IFBEG=0
      DO 220 I=1,NXD
        LMIX=0
        DO 110 L=1,NMIXA
          IF(MIXA(L).EQ.MIXD(I,J,K)) THEN
            LMIX=L
          ENDIF
  110   CONTINUE
        IF((AXD(I).EQ.1).AND.(AYD(J).EQ.1).AND.(LMIX.NE.0)) THEN
          DO 130 JP=1,NYP
          JF=JFBEG+JP
          DO 120 IP=1,NXP
            MIXF(IFBEG+IP,JF,K)=MIXA(LMIX+NMIXA)-1+MIXP(IP,JP)
  120     CONTINUE
  130     CONTINUE
          IFBEG=IFBEG+NXP
        ELSE
          NXPTMP=1
          IF(AXD(I).EQ.1) NXPTMP=NXP
          NYPTMP=1
          IF(AYD(J).EQ.1) NYPTMP=NYP
          DO 150 JP=1,NYPTMP
          JF=JFBEG+JP
          DO 140 IP=1,NXPTMP
            MIXF(IFBEG+IP,JF,K)=MIXD(I,J,K)
  140     CONTINUE
  150     CONTINUE
          IFBEG=IFBEG+NXPTMP
        ENDIF
  220 CONTINUE
      NYPTMP=1
      IF(AYD(J).EQ.1) NYPTMP=NYP
      JFBEG=JFBEG+NYPTMP
  230 CONTINUE
  240 CONTINUE
      
C Compute A-ZONE
      if(impx.ge.100)write(6,*) 'debug: Compute A-ZONE'
      IASS=0
      ALLOCATE(NBAX(NYD))
      ALLOCATE(IBAX(NYD))
      DO 340 J=1,NYD
      NBAX(J)=0
      IBAX(J)=0
      DO 330 I=1,NXD
      DO 320 K=1,NZD
        DO 310 L=1,NMIXA
          IF(MIXA(L).EQ.MIXD(I,J,K)) THEN
            IASS=IASS+1
            NBAX(J)=NBAX(J)+1
            IF(IBAX(J).EQ.0) IBAX(J)=I
            GOTO 330
          ENDIF
  310   CONTINUE
  320 CONTINUE
  330 CONTINUE
  340 CONTINUE
      NASS=IASS
*
      ALLOCATE(AZONE(NASS*NXP*NYP))
      IZ=0
      IASS=0
      DO 370 J=1,NYD
      DO 360 I=1,NBAX(J)
      IASS=IASS+1
      DO 365 JP=1,NYP
      DO 355 IP=1,NXP
      IZT=IZ+(JP-1)*NXP*NBAX(J)+(I-1)*NXP+IP
      AZONE(IZT)=IASS
  355 CONTINUE
  365 CONTINUE
  360 CONTINUE
      IZ=IZ+NBAX(J)*NXP*NYP
  370 CONTINUE
*      
      if(impx.ge.5)then
        write(6,*) 'New mixtures:'
        do K=1,NZF
          write(6,*) 'plane #',K
          do J=1,NYF
            write(6,*) (MIXF(I,J,K),I=1,NXF)
          enddo
        enddo

        write(6,*) 'Assembly zones:'
        IZ=0
        do J=1,NYD
          do K=1,NYP
            write(6,*) (AZONE(I),I=IZ+(K-1)*NBAX(J)*NXP+1,
     1      IZ+K*NBAX(J)*NXP)
          enddo
          IZ=IZ+NBAX(J)*NXP*NYP
        enddo
      endif
C Verify new mixture
      DO K=1,NZF
      DO J=1,NYF
      DO I=1,NXF
      IF(MIXF(I,J,K).EQ.-1) CALL XABORT('@NAPGEO: new '
     1 //'geometry mixture not assigned')
      ENDDO
      ENDDO
      ENDDO

      DEALLOCATE(MXD,MYD)
      DEALLOCATE(SXD,SYD)
      DEALLOCATE(MIXD)

      DEALLOCATE(MXP,MYP)
      DEALLOCATE(SXP,SYP)
      DEALLOCATE(MIXP)
C Compute relative position of assembly in original geometry
      if(impx.ge.100)write(6,*) 'debug: Compute relative position'
      JM=0
      JN=0
      LPOS=.TRUE.
      DO J=1,NYD
        IF(NBAX(J).NE.0) THEN
          JN=JN+1
          IF(LPOS) THEN
            JM=J
            LPOS=.FALSE.
          ENDIF
        ENDIF
      ENDDO

      IM=10000000
      DO J=JM,JM+JN-1
        IM=MIN(IBAX(J),IM)
      ENDDO
      DO J=1,NYD
        IBAX(J)=IBAX(J)-IM+1
      ENDDO

C Save heterogeneous core geometry
      if(impx.ge.100)write(6,*) 'debug: Save heter. core geometry'
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=ITYPGP
      ISTATE(3)=NXF
      ISTATE(4)=NYF
      ISTATE(5)=NZF
      ISTATE(6)=NREGF
      ISTATE(7)=NMIXF
*      ISTATE(39)=NMIXA
*      ISTATE(40)=NMIXP
      IF(LSPX .OR. LSPY .OR. LSPZ) ISTATE(11)=1
      CALL LCMPUT(IPGNW,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPGNW,'MESHX',NXF+1,2,MXF)
      CALL LCMPUT(IPGNW,'MESHY',NYF+1,2,MYF)
      CALL LCMPUT(IPGNW,'MESHZ',NZF+1,2,MZD)
      IF(LSPX) CALL LCMPUT(IPGNW,'SPLITX',NXF,1,SXF)
      IF(LSPY) CALL LCMPUT(IPGNW,'SPLITY',NYF,1,SYF)
      IF(LSPZ) CALL LCMPUT(IPGNW,'SPLITZ',NZF,1,SZD)
      CALL LCMPUT(IPGNW,'MIX',NREGF,1,MIXF)
      CALL LCMPUT(IPGNW,'NCODE',6,1,NCODE)
      CALL LCMPUT(IPGNW,'ICODE',6,1,ICODE)
      CALL LCMPUT(IPGNW,'ZCODE',6,2,ZCODE)
      CALL LCMPUT(IPGNW,'MIX-ASBLY',2*NMIXA,1,MIXA)
      CALL LCMPUT(IPGNW,'SIGNATURE',3,3,KCHAR)
      CALL LCMPUT(IPGNW,'A-ZONE',NASS*NXP*NYP,1,AZONE)
      CALL LCMPUT(IPGNW,'A-NX',JN,1,NBAX(JM))
      CALL LCMPUT(IPGNW,'A-IBX',JN,1,IBAX(JM))
      CALL LCMPUT(IPGNW,'A-NMIXP',1,1,NMIXP)
! 
      if(impx.ge.100)write(6,*) 'debug: beging deallacate'

      DEALLOCATE(IBAX)
      DEALLOCATE(NBAX)
      DEALLOCATE(AZONE)
      DEALLOCATE(MXF,MYF)
      DEALLOCATE(MZD)
      DEALLOCATE(SXF,SYF)
      DEALLOCATE(SZD)
      DEALLOCATE(MIXF)
      DEALLOCATE(MIXA)


      RETURN
      END
