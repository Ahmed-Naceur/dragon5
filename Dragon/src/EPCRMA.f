*DECK EPCRMA
      SUBROUTINE EPCRMA(IPMIC,IPRINT,NGR,NXS,NCV,NMIXT,NIFISS,IMIX,
     >                  ISOF,ITOTL,ISCAT,NAMDXS,DENSI,ICOV,COV,DRVAR,
     >                  XSREC,XSMAC,VAR,VAROLD,RST)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Add contribution of random error distribution to MACROLIB.
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPMIC   pointer to microlib.
* IPRINT  print level.
* NGR     number of groups.
* NXS     number of cross section types on EPC.
* NCV     maximum dimension of symmetrized covariance matrix.
* NMIXT   number of MIXTURES on MICROLIB.
* NIFISS  number of fissiles isotopes on MICROLIB.
* IMIX    mixture containing isotope.
* ISOF    fissile isotope number on MACROLIB.
* ITOTL   position of total XS in NAMDXS.
* ISCAT   position of scattering XS in NAMDXS.
* NAMDXS  names of XS.
* DENSI   isotope density.
* COV     variance and covariance matrix.
* ICOV    variance and covariance index.
* DRVAR   random distribution (width=1.).
* XSREC   microscopic vector XS.
*
*Parameters: input/output
* XSMAC   macroscopic vector XS contribution.
* VAR     multigroup variance.
* VAROLD  multigroup variance from previous correction.
* RST     ratio of scattering over total XS.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPMIC
      INTEGER          IPRINT,NGR,NXS,NCV,NMIXT,NIFISS,IMIX,ISOF,ITOTL,
     >                 ISCAT
      CHARACTER*6      NAMDXS(NXS)
      REAL             DENSI
      INTEGER          ICOV(NGR,NXS)
      REAL             COV(NCV,NXS)
      DOUBLE PRECISION DRVAR(NGR,NXS)
      REAL             XSREC(NGR,NXS),
     >                 XSMAC(NGR,NXS,NMIXT,NIFISS),VAR(NGR),
     >                 VAROLD(NGR),RST(NGR)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='EPCRMA')
*----
*  Local variables
*----
      INTEGER          IPRTL,IXS,IGR,ILOC,IDR,ILCMLN,ILCMTY
      INTEGER          KTOPT,ILOCS
      DOUBLE PRECISION CTOTL,CSCAT,CABS
*----
*  Write header
*----
      IPRTL=IPRINT
      IF(IPRTL .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Two options :
*  KTOPT = 0 : Scattering and total are independent
*  KTOPT = 1 : RST is scattering over total XS ratio
*  KTOPT = 2 : RST is absorption XS
*----
      KTOPT=2
      IF(KTOPT .EQ. 1) THEN
        DO IGR=1,NGR
          RST(IGR)=XSREC(IGR,ISCAT)/XSREC(IGR,ITOTL)
        ENDDO
      ELSE IF(KTOPT .EQ. 2) THEN
        DO IGR=1,NGR
          RST(IGR)=XSREC(IGR,ITOTL)-XSREC(IGR,ISCAT)
*          write(6,'(1X,I5,1P,3E15.7)')
*     >    IGR,XSREC(IGR,ITOTL),XSREC(IGR,ISCAT),RST(IGR)
        ENDDO
      ELSE
        CALL XDRSET(RST,NGR,0.0)
      ENDIF
*----
*  Process all vector cross section types
*----
      DO IXS=1,NXS
        ILOC=1
        ILOCS=1
        IF(NAMDXS(IXS) .EQ. 'NUSIGF') THEN
*        write(6,*) NAMDXS(IXS)
          IF(ISOF .GT. 0) THEN
            CALL LCMLEN(IPMIC,NAMDXS(IXS),ILCMLN,ILCMTY)
            IF(ILCMLN .EQ. NGR) THEN
              CALL LCMGET(IPMIC,NAMDXS(IXS),VAROLD)
            ELSE
              CALL XDRSET(VAROLD,NGR,0.0)
            ENDIF
            IDR=ISOF
            DO IGR=1,NGR
              VAR(IGR)=REAL(COV(ILOC,IXS)*XSREC(IGR,IXS)*DRVAR(IGR,IXS))
              IF(ABS(VAR(IGR))/ABS(XSREC(IGR,IXS)) .GT. 1.0) THEN
                WRITE(IOUT,*) 'VAR ',NAMDXS(IXS),IMIX,ISOF,IGR,
     >         VAR(IGR),XSREC(IGR,IXS),COV(ILOC,IXS),DRVAR(IGR,IXS)
              ELSE
                XSMAC(IGR,IXS,IMIX,IDR)=XSMAC(IGR,IXS,IMIX,IDR)+
     >                            DENSI*(VAR(IGR)-VAROLD(IGR))
              ENDIF
              ILOC=ILOC+ICOV(IGR,IXS)
            ENDDO
            IF(IPRTL .GE. 10) THEN
              WRITE(IOUT,*) 'VAR ',NAMDXS(IXS),IMIX,ISOF
              WRITE(IOUT,6010) (VAR(IGR),IGR=1,NGR)
              WRITE(IOUT,*) 'VAROLD ',NAMDXS(IXS),IMIX,ISOF
              WRITE(IOUT,6010) (VAROLD(IGR),IGR=1,NGR)
            ENDIF
            CALL LCMPUT(IPMIC,NAMDXS(IXS),NGR,2,VAR)
          ENDIF
        ELSE IF(IXS .EQ. ITOTL) THEN
*        write(6,*) NAMDXS(IXS)
          IDR=1
          CALL LCMLEN(IPMIC,NAMDXS(IXS),ILCMLN,ILCMTY)
          IF(ILCMLN .EQ. NGR) THEN
            CALL LCMGET(IPMIC,NAMDXS(IXS),VAROLD)
          ELSE
            CALL XDRSET(VAROLD,NGR,0.0)
          ENDIF
          IF(KTOPT .EQ. 1) THEN
*----
*  Assume the perturbation term comes from scattering with a contribution
*  proportional to RST and from total with a contribution
*  proportional to (1-RST)
*----
            DO IGR=1,NGR
              CTOTL=COV(ILOC,ITOTL)*XSREC(IGR,ITOTL)*DRVAR(IGR,ITOTL)
              CSCAT=COV(ILOCS,ISCAT)*XSREC(IGR,ISCAT)*DRVAR(IGR,ISCAT)
              VAR(IGR)=REAL((1.-RST(IGR))*CTOTL+RST(IGR)*CSCAT)
*          write(6,'(1X,I5,1P,6E15.7)')
*     >    IGR,XSREC(IGR,ITOTL),XSREC(IGR,ISCAT),RST(IGR),
*     >    VAR(IGR),CTOTL,CSCAT
              IF(ABS(VAR(IGR))/ABS(XSREC(IGR,IXS)) .GT. 1.0) THEN
                WRITE(IOUT,*) 'VAR ',NAMDXS(IXS),IMIX,ISOF,IGR,
     >          VAR(IGR),XSREC(IGR,IXS),COV(ILOC,IXS),DRVAR(IGR,IXS)
              ELSE
                XSMAC(IGR,IXS,IMIX,IDR)=XSMAC(IGR,IXS,IMIX,IDR)+
     >                            DENSI*(VAR(IGR)-VAROLD(IGR))
              ENDIF
              ILOC=ILOC+ICOV(IGR,IXS)
              ILOCS=ILOCS+ICOV(IGR,ISCAT)
            ENDDO
          ELSE IF(KTOPT .EQ. 2) THEN
*----
*  Assume the perturbation term comes from scattering and absorption
*----
            DO IGR=1,NGR
              CABS=COV(ILOC,ITOTL)*RST(IGR)*DRVAR(IGR,ITOTL)
              CTOTL=COV(ILOC,ITOTL)*XSREC(IGR,ITOTL)*DRVAR(IGR,ITOTL)
              CSCAT=COV(ILOCS,ISCAT)*XSREC(IGR,ISCAT)*DRVAR(IGR,ISCAT)
              VAR(IGR)=REAL(CABS+CSCAT)
*          write(6,'(1X,I5,1P,7E15.7)')
*     >    IGR,XSREC(IGR,ITOTL),XSREC(IGR,ISCAT),RST(IGR),
*     >    VAR(IGR),CABS,CSCAT,CTOTL
              IF(ABS(VAR(IGR))/ABS(XSREC(IGR,IXS)) .GT. 1.0) THEN
                WRITE(IOUT,*) 'VAR ',NAMDXS(IXS),IMIX,ISOF,IGR,
     >          VAR(IGR),XSREC(IGR,IXS),COV(ILOC,IXS),DRVAR(IGR,IXS)
              ELSE
                XSMAC(IGR,IXS,IMIX,IDR)=XSMAC(IGR,IXS,IMIX,IDR)+
     >                            DENSI*(VAR(IGR)-VAROLD(IGR))
              ENDIF
              ILOC=ILOC+ICOV(IGR,IXS)
              ILOCS=ILOCS+ICOV(IGR,ISCAT)
            ENDDO
          ELSE
*----
*  Assume total and scattering are independant
*----
            DO IGR=1,NGR
              CTOTL=COV(ILOC,ITOTL)*XSREC(IGR,ITOTL)*DRVAR(IGR,ITOTL)
              VAR(IGR)=REAL(CTOTL)
*          write(6,'(1X,I5,1P,6E15.7)')
*     >    IGR,XSREC(IGR,ITOTL),XSREC(IGR,ISCAT),RST(IGR),
*     >    VAR(IGR)
              IF(ABS(VAR(IGR))/ABS(XSREC(IGR,IXS)) .GT. 1.0) THEN
                WRITE(IOUT,*) 'VAR ',NAMDXS(IXS),IMIX,ISOF,IGR,
     >          VAR(IGR),XSREC(IGR,IXS),COV(ILOC,IXS),DRVAR(IGR,IXS)
              ELSE
                XSMAC(IGR,IXS,IMIX,IDR)=XSMAC(IGR,IXS,IMIX,IDR)+
     >                            DENSI*(VAR(IGR)-VAROLD(IGR))
              ENDIF
              ILOC=ILOC+ICOV(IGR,IXS)
            ENDDO
          ENDIF
          IF(IPRTL .GE. 10) THEN
            WRITE(IOUT,*) 'VAR ',NAMDXS(IXS),IMIX
            WRITE(IOUT,6010) (VAR(IGR),IGR=1,NGR)
            WRITE(IOUT,*) 'VAROLD ',NAMDXS(IXS),IMIX,ISOF
            WRITE(IOUT,6010) (VAROLD(IGR),IGR=1,NGR)
          ENDIF
          CALL LCMPUT(IPMIC,NAMDXS(IXS),NGR,2,VAR)
        ELSE
*        write(6,*) NAMDXS(IXS)
          IDR=1
          CALL LCMLEN(IPMIC,NAMDXS(IXS),ILCMLN,ILCMTY)
          IF(ILCMLN .EQ. NGR) THEN
            CALL LCMGET(IPMIC,NAMDXS(IXS),VAROLD)
          ELSE
            CALL XDRSET(VAROLD,NGR,0.0)
          ENDIF
          DO IGR=1,NGR
            VAR(IGR)=REAL(COV(ILOC,IXS)*XSREC(IGR,IXS)*DRVAR(IGR,IXS))
*          write(6,'(1X,I5,1P,6E15.7)')
*     >    IGR,XSREC(IGR,ITOTL),XSREC(IGR,ISCAT),XSREC(IGR,IXS),VAR(IGR)
            IF(ABS(VAR(IGR))/ABS(XSREC(IGR,IXS)) .GT. 1.0) THEN
              WRITE(IOUT,*) 'VAR ',NAMDXS(IXS),IMIX,ISOF,IGR,
     >        VAR(IGR),XSREC(IGR,IXS),COV(ILOC,IXS),DRVAR(IGR,IXS)
            ELSE
              XSMAC(IGR,IXS,IMIX,IDR)=XSMAC(IGR,IXS,IMIX,IDR)+
     >                          DENSI*(VAR(IGR)-VAROLD(IGR))
            ENDIF
            ILOC=ILOC+ICOV(IGR,IXS)
          ENDDO
          IF(IPRTL .GE. 10) THEN
            WRITE(IOUT,*) 'VAR ',NAMDXS(IXS),IMIX
            WRITE(IOUT,6010) (VAR(IGR),IGR=1,NGR)
            WRITE(IOUT,*) 'VAROLD ',NAMDXS(IXS),IMIX,ISOF
            WRITE(IOUT,6010) (VAROLD(IGR),IGR=1,NGR)
          ENDIF
          CALL LCMPUT(IPMIC,NAMDXS(IXS),NGR,2,VAR)
        ENDIF
      ENDDO
*----
*  Print if required
*----
      IF(IPRTL .GE. 10) THEN
        DO IXS=1,NXS
          IF(NAMDXS(IXS) .EQ. 'NUSIGF') THEN
            IF(ISOF .GT. 0) THEN
              WRITE(IOUT,*) 'MAC ',NAMDXS(IXS),IMIX,ISOF
              WRITE(IOUT,6010) (XSMAC(IGR,IXS,IMIX,ISOF),IGR=1,NGR)
            ENDIF
          ELSE
            WRITE(IOUT,*) 'MAC ',NAMDXS(IXS),IMIX
            WRITE(IOUT,6010) (XSMAC(IGR,IXS,IMIX,1),IGR=1,NGR)
          ENDIF
        ENDDO
      ENDIF
*----
*  Write header and return
*----
      IF(IPRTL .GE. 10) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(1P,10E12.4,:,' ...')
      END
