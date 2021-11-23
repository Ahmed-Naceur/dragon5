*DECK D2PADF
      SUBROUTINE D2PADF (IPDAT,IPRINT,NG,NMIL, ADF, NSF, DIFC,CURRN,
     1 SRFLX,ZAFLX,RPAR,IPAR,ADF_T,STAIDX,NVAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* CALL to GET_SAP_ADF to recover ADF information.
*
*Author(s): 
* J. Taforeau
*
*Parameters: 
* IPDAT    
* IPRINT   
* NG       
* NMIL     
* ADF      
* NSF      
* DIFC     
* CURRN    
* SRFLX    
* ZAFLX    
* RPAR     
* IPAR     
* ADF_T    
* STAIDX   
* NVAR     
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT,IPTH,KPTH
      INTEGER IPRINT,NG,NMIL,NSF,IPAR(3,NSF),NVAR
      REAL ADF(NSF,NG,10),DIFC(NG),CURRN(NSF,NG,2),SRFLX(NSF,NG),
     1 ZAFLX(NMIL,NG)
      DOUBLE PRECISION RPAR(6,NSF)
      CHARACTER*3 ADF_T
      INTEGER STAIDX(NVAR)   ! Index of current branch state values
*----
*  LOCAL VARIABLES
*----
      REAL SIDE,APOTHEM,VOLUME
      INTEGER :: NSD = 4
      INTEGER TOS(-1:1,-1:1)
      REAL SIGF(4)
      INTEGER DX, DY, SOT, IAXIS
      INTEGER  NAXIS,NPAIR(2),CAXIS(4),PAXIS(0:1,2),TRF_I(2,4)
      INTEGER  ICELL(2),NSURF(2)
      INTEGER IND, NZ, NC, IPAIR, IA, IP, NSURFAC,P, TR,NS,NGRP
      REAL*8 :: J_NET,J_PLUS,J_MINOS,FI_HET,TRANSV_CURR,FI_HOMOG,FAVE
      REAL*8 :: J_SUMM
      REAL :: B2_VECT(NMIL,NG), DIFF_C(NMIL,NG) ! B2 and D vectors
      REAL :: APOTH(NMIL,4)
      LOGICAL :: HASSYM(2,NMIL)
      INTEGER INTCORR(0:1,1,2)
      REAL CURR_INFO(1:(NMIL+1),NG,NSF,9)

      IF(NMIL > 1) CALL XABORT ('@D2P: MORE THAN 1 MIXTRURE ')
      IF(NSF .NE. NSD) CALL XABORT('@D2PADF: NUMBER OF SURFACE NE 4')

      SIDE= REAL(MAXVAL(RPAR(5,:)))
      APOTHEM= SIDE/2.0
      VOLUME= NSF*SIDE*APOTHEM/2.0
      CURR_INFO= 0.0
      ! TOS is the interface number corresponding to the cell
      ! to the right of the equation number (interface)
      TOS= 0
      TOS( 0, 1)= 4 !DX=0 DY>0 west
      TOS( 0,-1)= 2 !DX=0 DY<0 east
      TOS( 1, 0)= 1 !DX>0 DY=0 north
      TOS(-1, 0)= 3 !DX<0 DY=0 south

      SIGF(1)= 1.
      SIGF(2)= 1.
      SIGF(3)= 1.
      SIGF(4)= 1.

      !deltas in sense counterclokwise around the geometry
      !AXIS 1 DX>0 DY=0
      !AXIS 2 DX=0 DY>0
      NPAIR= 0
      NAXIS= 2

      INTCORR= 0
      !AXIS 1
      INTCORR(0,1,1)= 1
      INTCORR(1,1,1)= 3
      NPAIR(1)= 1
      !AXIS 2
      INTCORR(0,1,2)= 2
      INTCORR(1,1,2)= 4
      NPAIR(2)= 1

      !axis not crossing the surface
      CAXIS(1)= 1
      CAXIS(2)= 2
      CAXIS(3)= CAXIS(1)
      CAXIS(4)= CAXIS(2)
      !axis crossing a surface
      PAXIS(0,1)= 2
      PAXIS(1,1)= 4
      PAXIS(0,2)= 1
      PAXIS(1,2)= 3

      HASSYM= .FALSE.
      ! coefficient related to the transversal component of the J+.
      ! each surface has its 2 transversal components
      ! first surface
      TRF_I(1,1)= 2
      TRF_I(2,1)= 4

      ! 2-nd surface
      TRF_I(1,2)= 1
      TRF_I(2,2)= 3

      ! 3-th surface
      TRF_I(1,3)= 2
      TRF_I(2,3)= 4

      ! 4-th surface
      TRF_I(1,4)= 1
      TRF_I(2,4)= 3

      ADF=0.0
      SOT=0

      CURR_INFO= 0.0 !this is needed to know where to apply simmetries

      DO NS= 1,NSF

         ICELL(1)= IPAR(2,NS)
         ICELL(2)= IPAR(3,NS)

         IF(RPAR(3,NS).LT.-1.E-3) THEN
          DX = -1
         ELSEIF(RPAR(3,NS).GT.1.E-3) THEN
          DX = 1
         ELSE
          DX = 0
         ENDIF

         IF(RPAR(4,NS).LT.-1.E-3) THEN
          DY = -1
         ELSEIF(RPAR(4,NS).GT.1.E-3) THEN
          DY = 1
         ELSE
          DY = 0
         ENDIF
         ! check for the boundary regions

         IF(ICELL(1).LE.0) THEN
          ICELL(1)= NMIL+1
!           WRITE (*,*) 'BORDER TO THE RIGHT! MESH CH  ', ICELL(1)
         ENDIF

         IF(ICELL(2).LE.0) THEN
           ICELL(2)= NMIL+1
!           WRITE (*,*) 'BORDER TO THE LEFT! MESH CH   ', ICELL(2)
         ENDIF
         ! equations at the boundary:
         ! mesh on the left indicator of the surface ------------
         IF(TOS(DX,DY).EQ.1) SOT= 3
         IF(TOS(DX,DY).EQ.2) SOT= 4
         IF(TOS(DX,DY).EQ.3) SOT= 1
         IF(TOS(DX,DY).EQ.4) SOT= 2
         !
         !-------------------------------------------------------
         ! loop for the values of the J+-, J, FI
         DO NGRP= 1,NG
          ! J+
          CURR_INFO(ICELL(1),NGRP,TOS(DX,DY),1)=
     >     CURRN(NS,NGRP,2)/REAL(RPAR(5,NS))
          CURR_INFO(ICELL(2),NGRP,SOT,1)=
     >    CURRN(NS,NGRP,1)/REAL(RPAR(5,NS))
          ! J-
          CURR_INFO(ICELL(1),NGRP,TOS(DX,DY),2)=
     >    CURRN(NS,NGRP,1)/REAL(RPAR(5,NS))
          CURR_INFO(ICELL(2),NGRP,SOT       ,2)=
     >    CURRN(NS,NGRP,2)/REAL(RPAR(5,NS))

          ! J
          CURR_INFO(ICELL(1),NGRP,TOS(DX,DY),3)=
     >    CURR_INFO(ICELL(1),NGRP,TOS(DX,DY),1) -
     >    CURR_INFO(ICELL(1),NGRP,TOS(DX,DY),2)

          CURR_INFO(ICELL(2),NGRP,SOT       ,3)=
     >    CURR_INFO(ICELL(2),NGRP,SOT       ,1) -
     >    CURR_INFO(ICELL(2),NGRP,SOT       ,2)
          ! F-surf(het)
          IF(ICELL(1).EQ.(NMIL+1)) THEN
             IF(HASSYM(CAXIS(SOT),ICELL(2))) THEN
                CURR_INFO(ICELL(2),NGRP,SOT,4) = 0.0
             ELSE
                CURR_INFO(ICELL(2),NGRP,SOT,4) = SRFLX(NS,NGRP)
     >          / REAL(RPAR(5,NS))
             ENDIF
          ELSEIF(ICELL(2).EQ.(NMIL+1)) THEN
             IF(HASSYM(CAXIS(TOS(DX,DY)),ICELL(1))) THEN
                CURR_INFO(ICELL(1),NGRP,TOS(DX,DY),4) = 0.0
             ELSE
                CURR_INFO(ICELL(1),NGRP,TOS(DX,DY),4) =
     >           SRFLX(NS,NGRP)/REAL(RPAR(5,NS))
             ENDIF
          ELSE ! both cells are real
            CURR_INFO(ICELL(1),NGRP,TOS(DX,DY),4) =
     >       SRFLX(NS,NGRP)/REAL(RPAR(5,NS))
            CURR_INFO(ICELL(2),NGRP,SOT       ,4) =
     >       SRFLX(NS,NGRP)/REAL(RPAR(5,NS))
          ENDIF
         ! side dimension
          CURR_INFO(ICELL(1),NGRP,TOS(DX,DY),9)= REAL(RPAR(5,NS))
          CURR_INFO(ICELL(2),NGRP,SOT       ,9)= REAL(RPAR(5,NS))

          NSURF(1)= TOS(DX,DY)
          NSURF(2)= SOT
          DO IND= 1,2
           IF(ICELL(IND) < (NMIL+1)) THEN
            NZ= ICELL(IND)
            ! FI
            CURR_INFO(NZ,NGRP,:,5)=ZAFLX(NZ,NGRP)
           ENDIF
          ENDDO
         ENDDO
      ENDDO ! NS

      DO NC= 1,NMIL
         DO IAXIS= 1,NAXIS
          IF(HASSYM(IAXIS,NC)) THEN
           DO IPAIR= 1,NPAIR(IAXIS)
            ! put current value in the interface in front of it
            IF(CURR_INFO(NC,1,INTCORR(0,IPAIR,IAXIS),4).NE.0.) THEN
               CURR_INFO(NC,:,INTCORR(1,IPAIR,IAXIS),1:9)=
     >         CURR_INFO(NC,:,INTCORR(0,IPAIR,IAXIS),1:9)
            ELSEIF(CURR_INFO(NC,1,INTCORR(1,IPAIR,IAXIS),4).NE.0.)
     >       THEN
               CURR_INFO(NC,:,INTCORR(0,IPAIR,IAXIS),1:9)=
     >         CURR_INFO(NC,:,INTCORR(1,IPAIR,IAXIS),1:9)
            ENDIF
           ENDDO
          ENDIF
         ENDDO

         ! now put the possible zero dimension values
         DO IA= 1,NAXIS
          DO IP= 1,NPAIR(IA)
           ! put current value in the interface in front of it
           IF(CURR_INFO(NC,1,INTCORR(0,IP,IA),9) .NE. 0.) THEN
           ELSE
            CURR_INFO(NC,:,INTCORR(0,IP,IA),1:9) =
     >       CURR_INFO(NC,:,INTCORR(1,IP,IA),1:9)
           ENDIF
           IF(CURR_INFO(NC,1,INTCORR(1,IP,IA),9) .NE. 0.)THEN
           ELSE
            CURR_INFO(NC,:,INTCORR(1,IP,IA),1:9) =
     >       CURR_INFO(NC,:,INTCORR(0,IP,IA),1:9)
           ENDIF
          ENDDO
         ENDDO
      ENDDO ! NC

!-------------------------------------------------------
      DO NC= 1,NMIL
         DO NSURFAC= 1,NSD
          DO NGRP= 1,NG

           DIFF_C(NC,NGRP)= DIFC(NGRP)
           J_PLUS = CURR_INFO(NC,NGRP,NSURFAC,1)
           J_MINOS= CURR_INFO(NC,NGRP,NSURFAC,2)
           J_NET  = CURR_INFO(NC,NGRP,NSURFAC,3)
           FI_HET = CURR_INFO(NC,NGRP,NSURFAC,4)
           FAVE   = CURR_INFO(NC,NGRP,NSURFAC,5)

           APOTH(NC,NSURFAC)=
     >     CURR_INFO(NC,NGRP,PAXIS(0,CAXIS(NSURFAC)),9)/2.0
           CURR_INFO(NC,NGRP,NSURFAC,8)= APOTH(NC,NSURFAC)
           FI_HOMOG = SIGF(NSURFAC)*J_NET * APOTH(NC,NSURFAC)
     >     / DIFF_C(NC,NGRP) + FAVE
           ! FG:
           CURR_INFO(NC,NGRP,NSURFAC,6)= REAL(FI_HET / FI_HOMOG)
           ! FS:
           CURR_INFO(NC,NGRP,NSURFAC,7)= REAL(2. *
     >                ( J_PLUS + J_MINOS ) / FI_HOMOG)

          ENDDO !NGRP
         ENDDO !NSURFAC
      ENDDO !NC
      !
      ! B2 loop:
      !
      DO NCELL= 1,NMIL
         DO NGRP= 1,NG
          J_SUMM = SUM(CURR_INFO(NCELL,NGRP,:,3))

          B2_VECT(NCELL,NGRP)= REAL(J_SUMM / ( DIFF_C(NCELL,NGRP)
     >    * CURR_INFO(NCELL,NGRP,1,5) ))
         ENDDO
      ENDDO

      DO NCELL= 1,NMIL
         DO NGRP= 1,NG
          DO NSURFAC= 1,NSD
           ! TRANSVERSAL CURRENTS SUMMATION
           TRANSV_CURR= 0.
           DO TR= 1,2
            TRANSV_CURR= TRANSV_CURR +
     >      CURR_INFO(NCELL,NGRP,TRF_I(TR,NSURFAC),3)
           ENDDO
            ! no need to be stored !!!!
            ! CURR_INFO(NCELL,NGRP,NSURFAC,8)= TRANSV_CURR
          ENDDO
         ENDDO
      ENDDO
      ! store new IDF in the corresponding module to be used in
      ! writenemtab
      DO NCELL= 1,NMIL
        DO NGRP= 1,NG
          ! B2XS(K,NCELL,NGRP)=B2_VECT(NCELL,NGRP)
          DO NSURFAC= 1,NSD
            DO P=1,9
              ! 1 ->  J+
              ! 2 ->  J-
              ! 3 ->  J
              ! 4 ->  F-surf
              ! 5 ->  F-ave
              ! 6 ->  GET_IDF
              ! 7 ->  SEL_IDF
              ! 8 ->  apotheme
              ! 9 ->  side length
              ADF(NSURFAC,NGRP,P)=CURR_INFO(NCELL,NGRP,NSURFAC,P)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF(IPRINT > 1) THEN
        WRITE(6,*) "*** RECOVER ASSEMBLY DISCONTINUITY FACTOR ***"
        IF(ADF_T.EQ.'GET') WRITE(6,*) "ADF TYPE : GET "
        IF(ADF_T.EQ.'SEL') WRITE(6,*) "ADF TYPE : SELENGUT "
       DO NGRP=1, NG
        WRITE(6,*) "GROUP             :",NGRP
         IF(ADF_T.EQ.'GET') WRITE(6,*)"ADF(N/E/S/W) :",ADF (:,NGRP,6)
         IF(ADF_T.EQ.'SEL') WRITE(6,*)"ADF(N/E/S/W) :",ADF (:,NGRP,7)
       ENDDO
      ENDIF

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      IPTH=LCMGID(IPDAT,'CROSS_SECT')
      KPTH=LCMDIL(IPTH,STAIDX(NVAR))
      CALL LCMSIX(KPTH,'MACROLIB_XS',1)
      IF(ADF_T.EQ.'GET') THEN
       CALL LCMPUT(KPTH,'ADF',NSF*NG,2,ADF(:,:,6))
      ELSEIF(ADF_T.EQ.'SEL') THEN
       CALL LCMPUT(KPTH,'ADF',NSF*NG,2,ADF(:,:,7))
      ELSE
       CALL XABORT('@D2PADF: UNKNOW ADF TYPE'//ADF_T//'.')
      ENDIF
      END
