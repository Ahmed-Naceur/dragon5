*DECK LIBND3
      SUBROUTINE LIBND3(NGF,NGFR,NGRO,NBDIL,SN,SB,DILUS,DELTA,NF,XA,XS,
     > XF,XN,GAR1,SCAT,GAR2,WT0)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolate dilution-tabulated information, perform 
* Livolant-Jeanpierre normalization and compute self-shielded 
* cross sections at specific dilution.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: input
* NGF     number of fast groups without self-shielding.
* NGFR    number of fast and resonance groups.
* NGRO    number of energy groups.
* NBDIL   number of dilutions.
* SN      actual dilution of the nuclide.
* SB      dilution of the nuclide, as used in Livolant-Jeanpierre
*         normalization.
* DILUS   tabulation points in dilution.
* DELTA   lethargy widths.
* NF      flag set to 3 if fission information is present.
* XA      tabulated absorption effective reaction rates.
* XS      tabulated scattering effective reaction rates.
* XF      tabulated nu*fission effective reaction rates.
* XN      tabulated NJOY flux.
* GAR1    infinite dilution cross sections.
* SCAT    infinite dilution P0 differential scattering cross sections.
*
*Parameters: output
* GAR2    interpolated self-shielded cross sections.
* WT0     NJOY flux.
*
*Reference:
* Copyright (C) from NDAS Atomic Energy of Canada Limited utility (2006)
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  Subroutine arguments
*----
      INTEGER NGF,NGFR,NGRO,NBDIL,NF
      REAL SN(NGRO),SB(NGRO),DILUS(NBDIL),DELTA(NGRO),
     1 XA(NGFR-NGF,NBDIL),XS(NGFR-NGF,NBDIL),XF(NGFR-NGF,NBDIL),
     2 XN(NGFR-NGF,NBDIL),GAR1(NGRO,6),SCAT(NGRO,NGRO),GAR2(NGRO,6),
     3 WT0(NGRO)
*----
*  Local variables
*----
      REAL WW,ZNGAR,SSFACT,AUX,XN3
      INTEGER I,IG,IG2
      CHARACTER HSMG*131
      LOGICAL LCUBIC
      PARAMETER(LCUBIC=.TRUE.)
      REAL, ALLOCATABLE, DIMENSION(:) :: TERPD,DD,XA2,XS2,XF2,WK
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LINF
*----
*  Scratch storage allocation
*----
      ALLOCATE(TERPD(NBDIL),LINF(NGRO))
*----
*  Dilution interpolation
*----
      ALLOCATE(DD(NBDIL))
      DO I=1,NBDIL
        DD(I)=LOG10(DILUS(I))
      ENDDO
      ALLOCATE(XA2(NGRO),XS2(NGRO),XF2(NGRO))
      DO IG=1,NGF
        XA2(IG)=GAR1(IG,2)
        XS2(IG)=GAR1(IG,5)
        XF2(IG)=GAR1(IG,4)
        WT0(IG)=1.0
      ENDDO
      ALLOCATE(WK(3*NBDIL))
      DO IG=NGF+1,NGFR
        CALL XDRSET(TERPD,NBDIL,0.0)
        LINF(IG)=SN(IG).EQ.DILUS(NBDIL)
        DO I=1,NBDIL
          IF(ABS(SN(IG)-DILUS(I)).LE.1.0E-5*ABS(SN(IG))) THEN
            TERPD(I)=1.0
            GO TO 10
          ENDIF
        ENDDO
        IF((NBDIL.EQ.1).OR.(SN(IG).GE.DILUS(NBDIL))) THEN
*         No interpolation above infinite dilution
          TERPD(NBDIL)=1.0
        ELSE IF((NBDIL.EQ.2).OR.(SN(IG).GE.DILUS(NBDIL-1))) THEN
*         One over SN interpolation near infinite dilution
          LINF(IG)=.TRUE.
          TERPD(NBDIL-1)=DILUS(NBDIL-1)/SN(IG)
          TERPD(NBDIL)=1.0-DILUS(NBDIL-1)/SN(IG)
        ELSE
*         Perform Ceschino cubic interpolation
          CALL ALTERP(LCUBIC,NBDIL-1,DD,LOG10(SN(IG)),.FALSE.,
     >    TERPD,WK)
        ENDIF
   10   XA2(IG)=0.0
        XS2(IG)=0.0
        XF2(IG)=0.0
        WT0(IG)=0.0
        DO I=1,NBDIL
          WW=TERPD(I)
          IF(ABS(WW).GT.1.0E-6) THEN
            XA2(IG)=XA2(IG)+WW*XA(IG-NGF,I)
            XS2(IG)=XS2(IG)+WW*XS(IG-NGF,I)
            IF(NF.EQ.3) XF2(IG)=XF2(IG)+WW*XF(IG-NGF,I)
            WT0(IG)=WT0(IG)+WW*XN(IG-NGF,I)
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(WK,DD)
*----
*  Livolant-Jeanpierre normalization
*----
      DO IG=NGF+1,NGFR
        IF(SB(IG).NE.SN(IG)) THEN
          ZNGAR=-(XA2(IG)+XS2(IG))
          DO IG2=1,IG
            SSFACT=XS2(IG2)/GAR1(IG2,5)
            ZNGAR=ZNGAR+SCAT(IG,IG2)*SSFACT*DELTA(IG2)/DELTA(IG)
          ENDDO
          XN3=(WT0(IG)-1.0)*SN(IG)
          IF(LINF(IG)) THEN
*           Use an interpolated value near infinite dilution
            AUX=(DILUS(NBDIL-1)/SB(IG))**2
            XN3=AUX*XN3+(1.0-AUX)*ZNGAR
            XN3=1.0+XN3/SB(IG)
          ELSE
            XN3=1.0+XN3/SB(IG)
          ENDIF
          IF((XN3.LE.0.0).OR.(XN3.GT.2.0)) THEN
            WRITE (HSMG,100) XN3,IG,SB(IG),SN(IG)
            CALL XABORT(HSMG)
          ELSE IF(XN3.GT.1.2) THEN
            WRITE (HSMG,100) XN3,IG,SB(IG),SN(IG)
            WRITE(6,'(1X,A)') HSMG
          ENDIF
          WT0(IG)=XN3
        ENDIF
      ENDDO
*----
*  Divide effective reaction rates by NJOY flux for obtaining
*  self-shielded cross sections
*----
      DO I=1,6
        DO IG=1,NGRO
          GAR2(IG,I)=GAR1(IG,2)
        ENDDO
      ENDDO
      DO IG=NGF+1,NGFR
*       Absorption xs
        GAR2(IG,2)=XA2(IG)/WT0(IG)
*       P0 scattering xs
        GAR2(IG,5)=XS2(IG)/WT0(IG)
*       nu*fission xs
        GAR2(IG,4)=XF2(IG)/WT0(IG)
*       Transport-corrected total xs
        GAR2(IG,1)=GAR1(IG,1)*(GAR2(IG,2)+GAR2(IG,5))/(GAR1(IG,2)+
     >  GAR1(IG,5))
*       Fission xs
        IF((NF.EQ.3).AND.(XF2(IG).EQ.0.0)) THEN
          GAR2(IG,3)=GAR1(IG,3)
        ELSE IF(NF.EQ.3) THEN
          GAR2(IG,3)=GAR1(IG,3)*GAR2(IG,4)/GAR1(IG,4)
        ENDIF
*       P1 scattering xs
        GAR2(IG,6)=GAR1(IG,6)*GAR2(IG,5)/GAR1(IG,5)
      ENDDO
      DEALLOCATE(XF2,XS2,XA2)
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(LINF,TERPD)
      RETURN
*
  100 FORMAT(37H LIBND3: Invalid value of NJOY flux (,1P,E11.3,
     1 10H) in group,I4,11H. Dilution=,E11.3,2H (,E11.3,2H).)
      END
