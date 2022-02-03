*DECK COMGFF
      SUBROUTINE COMGFF(MPCPO,IPEDI2,FNORM,NGFF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the group form factor information from an edition object.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* MPCPO   pointer to a microlib directory of the multicompo.
* IPEDI2  pointer to the edition object containing group form factor
*         information (L_EDIT signature).
* FNORM   flux normalization factor.
* NGFF    number of form factors per energy group (set to -1 if not
*         initialized).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) MPCPO,IPEDI2
      REAL FNORM
      INTEGER NGFF
*----
*  LOCAL PARAMETERS
*----
      TYPE(C_PTR) JPEDI2,KPEDI2
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLUME
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLUX,HFACT,SIGF
*----
*  RECOVER GFF INFO FROM THE ROOT OF THE EDITION OBJECT
*----
      CALL LCMGTC(IPEDI2,'SIGNATURE',12,1,TEXT12)
      IF(TEXT12.NE.'L_EDIT') THEN
         CALL XABORT('COMGFF: SIGNATURE OF OBJECT IS '//TEXT12//
     1   '. L_EDIT EXPECTED.')
      ENDIF
      CALL LCMGET(IPEDI2,'STATE-VECTOR',ISTATE)
      IF(NGFF.EQ.-1) THEN
         NGFF=ISTATE(1)
      ELSE IF(NGFF.NE.ISTATE(1)) THEN
         CALL XABORT('COMGFF: INVALID NUMBER OF FORM FACTORS IN EDITIO'
     1   //'N OBJECT.')
      ENDIF
      IF(ISTATE(20).EQ.0) CALL XABORT('COMGFF: MISSING MACRO-GEOMETRY '
     1 //'IN EDITION OBJECT.')
      CALL LCMSIX(MPCPO,'MACROLIB',1)
      CALL LCMSIX(MPCPO,'GFF',1)
*----
*  RECOVER THE MACRO-GEOMETRY
*----
      CALL LCMSIX(IPEDI2,'MACRO-GEOM',1)
      CALL LCMSIX(MPCPO,'GFF-GEOM',1)
      CALL LCMEQU(IPEDI2,MPCPO)
      CALL LCMSIX(MPCPO,' ',2)
      CALL LCMSIX(IPEDI2,' ',2)
*----
*  RECOVER GFF INFO FROM THE LAST-EDIT DIRECTORY IN THE EDITION OBJECT
*----
      CALL LCMGTC(IPEDI2,'LAST-EDIT',12,1,TEXT12)
      CALL LCMSIX(IPEDI2,TEXT12,1)
         CALL LCMSIX(IPEDI2,'MACROLIB',1)
         CALL LCMGET(IPEDI2,'STATE-VECTOR',ISTATE)
         NG=ISTATE(1)
         ALLOCATE(VOLUME(NGFF),FLUX(NGFF,NG),HFACT(NGFF,NG),
     1   SIGF(NGFF,NG))
         IF(NGFF.NE.ISTATE(2)) THEN
            CALL XABORT('COMGFF: INVALID NUMBER OF FORM FACTORS IN MAC'
     1      //'ROLIB OF THE EDTION OBJECT.')
         ENDIF
         CALL LCMGET(IPEDI2,'VOLUME',VOLUME)
         JPEDI2=LCMGID(IPEDI2,'GROUP')
         DO IG=1,NG
            KPEDI2=LCMGIL(JPEDI2,IG)
            CALL LCMGET(KPEDI2,'FLUX-INTG',FLUX(1,IG))
            DO IBM=1,NGFF
              FLUX(IBM,IG)=FNORM*FLUX(IBM,IG)/VOLUME(IBM)
            ENDDO
            CALL LCMGET(KPEDI2,'H-FACTOR',HFACT(1,IG))
            CALL LCMGET(KPEDI2,'NFTOT',SIGF(1,IG))
         ENDDO
         CALL LCMPUT(MPCPO,'VOLUME',NGFF,2,VOLUME)
         CALL LCMPUT(MPCPO,'NWT0',NGFF*NG,2,FLUX)
         CALL LCMPUT(MPCPO,'H-FACTOR',NGFF*NG,2,HFACT)
         CALL LCMPUT(MPCPO,'NFTOT',NGFF*NG,2,SIGF)
         DEALLOCATE(SIGF,HFACT,FLUX,VOLUME)
         CALL LCMSIX(IPEDI2,' ',2)
      CALL LCMSIX(IPEDI2,' ',2)
      CALL LCMSIX(MPCPO,' ',2)
*----
*  SET STATE-VECTOR INDEX FOR MICROLIB IN MULTICOMPO
*----
      CALL LCMLEN(MPCPO,'STATE-VECTOR',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
         CALL LCMGET(MPCPO,'STATE-VECTOR',ISTATE)
         ISTATE(16)=NGFF
         CALL LCMPUT(MPCPO,'STATE-VECTOR',NSTATE,1,ISTATE)
      ENDIF
      CALL LCMSIX(MPCPO,' ',2)
      RETURN
      END
