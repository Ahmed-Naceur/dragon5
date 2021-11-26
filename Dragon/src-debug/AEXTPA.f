*DECK AEXTPA
      SUBROUTINE AEXTPA(NOMFIC,ISFICH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Determination of a SAPHYR archive file characteristics.
* Component of a FORTRAN-77 emulator of the SAPHYR archive system.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NOMFIC  name of SAPHYR archive file.
*
*Parameters: output
* ISFICH  file characteristics with:
*         ISFICH(1) is the address of the table of content;
*         ISFICH(2) is the number of archive objects on file;
*         ISFICH(3) is the direct access record length in words.
*
*-----------------------------------------------------------------------
*
      INTEGER    ISFICH(3)
      CHARACTER  NOMFIC*(*),HSMG*131
*
      IULFIC = KDROPN(NOMFIC,2,4,1)
      IF(IULFIC.LE.0) THEN
         WRITE(HSMG,'(33HAEXTPA: KDROPN FAILURE WITH CODE=,I3)') IULFIC
         CALL XABORT(HSMG)
      ENDIF
      ISTATE = 5
      I2 = 3
*
   40 READ(IULFIC,REC=I2,ERR=50,IOSTAT=IOS) MOTLU
      IF(IOS.NE.0) GO TO 50
*
      ISTATE = ISTATE + 1
      IF(ISTATE .EQ. 8) THEN
        ISFICH(3) = MOTLU
        I2 = 4
        GO TO 40
      ELSEIF(ISTATE .EQ. 9) THEN
        ISFICH(2) = MOTLU
      ELSEIF(ISTATE .EQ. 7) THEN
        I2 = MOTLU + 7
        GO TO 40
      ELSEIF(ISTATE .EQ. 6) THEN
        ISFICH(1) = MOTLU
        I2 = MOTLU + 3
        GO TO 40
      ENDIF
*
   50 IER = KDRCLS(IULFIC,1)
      IF(IER.LT.0) THEN
         WRITE(HSMG,'(33HAEXTPA: KDRCLS FAILURE WITH CODE=,I3)') IER
         CALL XABORT(HSMG)
      ENDIF
      RETURN
      END
