*DECK SCRFND
      SUBROUTINE SCRFND(MAXISO,NBISOI,NBISO,INAME,IBM,HRESID,HUSE,HNAME,
     > IMIX,JSO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the isotope index of an isotope in the microlib.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* MAXISO  maximum number of isotopes in the microlib.
* NBISOI  initial number of isotopes in the microlib.
* NBISO   exact number of isotopes in the microlib.
* INAME   name of an isotope.
* IBM     mixture in which the isotope is present.
* HRESID  character*8 name of the residual isotope in the Saphyb.
* HUSE    alias names of microlib isotopes.
* HNAME   reference name of microlib isotopes.
* IMIX    full-core mixture belonging to each isotope.
*
*Parameters: output
* NBISO   exact number of isotopes in the microlib.
* HUSE    names of microlib isotopes.
* HNAME   reference name of microlib isotopes.
* IMIX    full-core mixture belonging to each isotope.
* JSO     position of isotope INAME in isotope list.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXISO,NBISOI,NBISO,INAME(2),IBM,HUSE(3,MAXISO),
     > HNAME(3,MAXISO),IMIX(MAXISO),JSO
      CHARACTER HRESID*8
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT4*4,TEXT8*8
      INTEGER IBLANK, I0, ISO, ISAVE
      INTEGER IHRES(2)
      SAVE IBLANK,IHRES,ISAVE
      DATA TEXT4,TEXT8/'    ','*MAC*RES'/
      DATA ISAVE/0/
*
      IF(ISAVE.EQ.0) THEN
        READ(TEXT4,'(A4)') IBLANK
        READ(TEXT8,'(2A4)') IHRES(1),IHRES(2)
        ISAVE=1   
      ENDIF
*
      JSO=0
      DO ISO=1,NBISOI
        IF(IMIX(ISO).NE.IBM) CYCLE
        IF((INAME(1).EQ.HUSE(1,ISO)).AND.(INAME(2).EQ.HUSE(2,ISO))) THEN
          JSO=ISO
          RETURN
        ENDIF
      ENDDO
      NBISO=NBISO+1
      IF(NBISO.GT.MAXISO) CALL XABORT('SCRFND: MAXISO OVERFLOW.')
      JSO=NBISO
      HUSE(1,JSO)=INAME(1)
      HUSE(2,JSO)=INAME(2)
      HUSE(3,JSO)=IBLANK
      IF((INAME(1).EQ.IHRES(1)).AND.(INAME(2).EQ.IHRES(2))) THEN
        READ(HRESID,'(2A4)') (HNAME(I0,JSO),I0=1,2)
      ELSE
        HNAME(1,JSO)=INAME(1)
        HNAME(2,JSO)=INAME(2)
      ENDIF
      HNAME(3,JSO)=IBLANK
      IMIX(JSO)=IBM
      RETURN
      END
