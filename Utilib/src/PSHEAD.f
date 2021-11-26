*DECK PSHEAD
      SUBROUTINE PSHEAD(ISPSP,NAMPSP,PROGNM)
C
C---------------------------  PSHEAD  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSHEAD
C      USE      : SET POSTSCRIPT HEADER
C                 REPLACES PART OF PSPLOT ROUTINE PSINIT
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      NAMPSP   : PSP FILE NAME                          C*12
C      PROGNM   : PAGE PROGRAM NAME                      C*6
C---------------------------   PSHEAD  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP
      CHARACTER        NAMPSP*12,PROGNM*6
C----
C  LOCAL VARIABLES
C----
      REAL             CONVER
      CHARACTER        NAMSBR*6
      PARAMETER       (CONVER=72.0,NAMSBR='PSHEAD')
      CHARACTER        CMDSTR*132
C----
C  PREPARE HEADER
C----
      CMDSTR='%!PS-Adobe-2.0 EPSF-2.0'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR= '%%Title: '//NAMPSP
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR= '%%CreationDate: 1999/03/29'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR= '%%Created with: PSPLOT PostScript Plotting Package'//
     >        ' in '//PROGNM
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR= '%%Reference: Kevin E. Kohler '//
     >        '<kevin@ocean.nova.edu> '//
     >        '- DRAGON implementation'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='%%EndComments'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/inch {72 mul} bind def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Ah {moveto lineto lineto stroke} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Ar {moveto 2 copy lineto 4 -2 roll'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     moveto lineto lineto stroke } def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/arcit {S /A2 exch def /A1 exch def /Rad exch def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        /Yc exch def /Xc exch def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        Xc Rad A1 cos mul add Yc Rad A1 sin mul add'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        moveto newpath'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        Xc Yc Rad A1 A2 arc stroke} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/C {/Rad exch def /Yc exch def /Xc exch def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         Xc Yc Rad 0 360 arc closepath'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        } def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/c0sf {closepath 0 setgray fill} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/cf {closepath fill} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Cs {closepath stroke} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Cln {newpath 3 1 roll'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='      moveto {lineto} repeat clip newpath'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     } def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Cs {closepath stroke} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Fb {newpath moveto '
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR=' Dx 0 rlineto 0 Dy rlineto Dx neg 0 rlineto closepath'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR=' fill } def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Fbn { newpath 3 1 roll moveto {lineto} repeat'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='       closepath fill } def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Fbnc { newpath 3 1 roll moveto'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='       {lineto} repeat closepath fill } def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/L /lineto load def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Lend {/Strlen exch stringwidth pop def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Lendi {/Strlen exch stringwidth pop 1.5 mul def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Lends {/Strlen exch stringwidth pop 1.1 mul def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Lenssd '//
     >       '{/Strlenss exch stringwidth pop 3 mul 4 div def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/LSM {2 copy lineto stroke moveto} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/lsm {Xp Yp lineto stroke mover} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/M /moveto load def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/mover {Xp Yp moveto} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Np {newpath} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/S /stroke load def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Scrgb {setrgbcolor} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Scmyk {setcmykcolor} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Schsb {sethsbcolor} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Sgray {setgray} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/FScrgb {setrgbcolor fill} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/FScmyk {setcmykcolor fill} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/FSchsb {sethsbcolor fill} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/FSgray {setgray fill} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Setf {Curfnt exch scalefont setfont} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/SM {stroke moveto} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/sm {stroke mover} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR=' '
      WRITE(CMDSTR,'(6H/Slw {,f7.4,22H mul setlinewidth} def)') CONVER
      CALL PSCPUT(ISPSP,CMDSTR)
      WRITE(CMDSTR,'(7H/SSlw {,f7.4,29H mul setlinewidth stroke} def)')
     >                         CONVER
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Slw0 {.24 setlinewidth} bind def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/SSlw0 {.24 setlinewidth stroke} bind def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR= '%Line Breaking Procedure'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/TurnLineFL'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='   { /T exch def /spacewidth space stringwidth pop def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     /currentw 0 def /wordspace_count 0 def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     /restart 0 def  /remainder T def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     {remainder space search'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='       {/nextword exch def pop'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        /remainder exch def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        /nextwordwidth nextword stringwidth pop def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        currentw nextwordwidth add lw gt'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        {T restart wordspace_count restart sub'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         getinterval showline'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         /restart wordspace_count def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         /currentw nextwordwidth spacewidth add def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        }'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        {/currentw currentw nextwordwidth add'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         spacewidth add def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        } '
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        ifelse'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        /wordspace_count wordspace_count'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='        nextword length add 1 add def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='       }'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='       {pop exit}'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='       ifelse'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     } loop'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     /lrem remainder stringwidth pop def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     currentw lrem add lw gt'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     {T restart wordspace_count restart sub '
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='      getinterval showline remainder showline}'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     {/lastchar T length def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='      T restart lastchar restart sub getinterval '
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='      lm y moveto show}'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='     ifelse'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='   } def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR=' /parms {/y exch def /lm exch def /rm exch def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         /leading exch def /pointsize exch def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         /lw rm lm sub def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         findfont pointsize scalefont setfont '
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         /showline {lm y moveto show'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         /y y leading sub def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='         lm y moveto } def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Xposd {/Xpos exch def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Xposjd  '//
     >       '{/Xpos exch Xpos exch Strlen mul sub def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/xydef {/Xp Xpos def /Yp Ypos def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='%/Xypd {/Yp exch def /Xp exch def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Xypos0d {/Xpos0 Xpres def /Ypos0 Ypres def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Xyprset {dup '//
     >       '/Xpres exch cos Strlen mul Xpos add def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='              '//
     >       '/Ypres exch sin Strlen mul Ypos add def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Xyprset0 {dup '//
     >       '/Xpres exch cos Strlen mul Xpos0 add def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='               '//
     >       '/Ypres exch sin Strlen mul Ypos0 add def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Yposd {/Ypos exch def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/Yposjd  '//
     >       '{/Ypos exch Ypos exch Strlen mul sub def} def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='%%EndProlog'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='/space ( ) def'
      CALL PSCPUT(ISPSP,CMDSTR)
      RETURN
      END
