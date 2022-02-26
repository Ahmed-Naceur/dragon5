*DECK PSPFIL
      SUBROUTINE PSPFIL(ISPSP,JSPSP,NAMPSP,NPAGE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* PSP file analysis.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input/output
* ISPSP   PSP file unit.
* JSPSP   PSP file mode:
*         = 0 new;
*         = 1 update.
* NAMPSP  PSP file name.
* NPAGE   page number.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6,PROGNM*6
      PARAMETER       (IOUT=6,NAMSBR='PSPFIL',PROGNM='DRAGON')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          JSPSP,ISPSP,NPAGE
      CHARACTER        NAMPSP*12
*----
*  LOCAL VARIABLES
*----
      INTEGER          IRL,IDR,ILINE,IPF,IPN
      CHARACTER        CMDSTR*132,CFMT*16
      REAL             XYPOS(2)
      NPAGE=0
      IF(JSPSP .EQ. 1) THEN
*----
*  TEST IF ADEQUATE DRAGON PS FILE TYPE
*----
        DO 100 IRL=1,3
          READ(ISPSP,'(A132)') CMDSTR
 100    CONTINUE
        READ(ISPSP,'(A132)') CMDSTR
        IDR=INDEX(CMDSTR,PROGNM)
        IF(IDR .EQ. 0) CALL XABORT(NAMSBR//
     >    ': NOT A DRAGON GENERATED POSTSCRIPT FILE')
        ILINE=0
        IPF=1
*----
*  LOCATE LAST PAGE NUMBER
*----
 110    CONTINUE
          READ(ISPSP,'(A132)',END=115) CMDSTR
          IPN=INDEX(CMDSTR,'%%Page')
          IF(IPN .NE. 0) THEN
            IPN=INDEX(CMDSTR,' ')
            IPF=INDEX(CMDSTR(IPN+1:132),' ')-1
            CFMT=' '
            WRITE(CFMT,'(2H(I,I1,1H))') IPF
            READ(CMDSTR(IPN+1:IPN+IPF),CFMT) NPAGE
          ENDIF
          GO TO 110
 115    CONTINUE
        BACKSPACE ISPSP
*----
*  SET NEXT PAGE NUMBER AND PREPARE FOR OUTPUT
*----
        NPAGE=NPAGE+1
        XYPOS(1)=0.5
        XYPOS(2)=0.5
        CALL PSPAGE(ISPSP,NPAGE,XYPOS)
      ELSE
        CALL PSHEAD(ISPSP,NAMPSP,PROGNM)
        NPAGE=NPAGE+1
        XYPOS(1)=0.5
        XYPOS(2)=0.5
        CALL PSPAGE(ISPSP,NPAGE,XYPOS)
      ENDIF
      RETURN
      END
