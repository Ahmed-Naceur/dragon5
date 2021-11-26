!
!-----------------------------------------------------------------------
!
!Purpose:
! Generate a Postscript representation of the surfacic geometry.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s):
! G. Civario (CS-SI)
!
!Comments:
! Ce fichier est derive de la bibliotheque PSPLOT de Kevin E. Kohler,
! developpee au Nova Southeastern University Oceanographic Center en Floride.
! Le code inital a simplement ete encapsule dans un module, pour eviter les
! eventuels conflics de noms, et ampute de toutes les routines non utilisees
! ici.
! Quelques modifications mineures ont par ailleurs ete faites, pour permettre
! la production d'un fichier eps, et assurer le centrage de la figure.
!
!-----------------------------------------------------------------------
!
module derivedPSPLOT
  !!implicit none

  logical,save :: g_psp_isEpsFile !true if file is eps
  real ,save   :: g_psp_bBoxXmin,g_psp_bBoxYmin,g_psp_bBoxXmax,g_psp_bBoxYmax

contains

  subroutine line(fx,fy,tx,ty)
    implicit none
    double precision,intent(in)    :: fx,fy,tx,ty

    call PLOT(real(fx),real(fy),3)
    call PLOT(real(tx),real(ty),2)
  end subroutine line


  subroutine arc(xc,yc,rad,ang1,ang2)
    character*132 cmdstr
    common/plt1/cmdstr
    common/cnvcom/conver
    radi=rad*conver
    xci=xc*conver
    yci=yc*conver
    cmdstr=' '
    write(cmdstr,'(f10.2,'' '',f10.2,'' '',f10.2,'' '',2f10.2,'' arcit'')')&
         & xci,yci,radi,ang1,ang2
    call filler
  end subroutine arc

  subroutine filler
    !nfild is the last position filled in the compressed aaa buffer array
    !work is a work array used to load array aaa
    character*132 cmdstr,cmdc(132)*1
    common/plt1/cmdstr
    common/outcom/iunit
    logical ispace
    !equivalence (cmdstr,cmdc(1))
    ibslash=92
    lc=0
    lcc=lenstr(cmdstr,132)
    ispace=.false.

    !itot is running total of left/right parentheses in text string
    !if itot=0 then we are not in text mode, i.e. left=right

    itot=0
    icclst=-999
    do l=1,lcc
       icc=ichar(cmdstr(l:l))
       if(icc.eq.32.and.ispace.and.itot.eq.0) cycle !Don't place 2 or more
       ! spaces together if
       ! not in text mode
       if(icc.ge.32.and.icc.le.127) then
          lc=lc+1
          cmdc(lc)=cmdstr(l:l)
       endif
       if(icc.eq.32) then
          ispace=.true.
       else
          ispace=.false.
       endif
       if(icc.eq.40.and.icclst.ne.ibslash)itot=itot+1
       if(icc.eq.41.and.icclst.ne.ibslash)itot=itot-1
       icclst=icc
    end do

    !Write cmdstr to output file
    write(iunit,'(132a1)')(cmdc(ii),ii=1,lc)
    return
  end subroutine filler

  subroutine psinit(fileNbr,portrait)
    !initializes plot for hp plotter
    integer fileNbr
    logical first,portrait,prtrt
    character*132 cmdstr,curfnt
    character*80 fileout
    character tim*10,dat*8,zn*5
    integer*4 val(8)
    character*1 timer(8),dater(9)
    equivalence(timer(1),tim),(dater(1),dat)
    common/conre1/ioffp,spval
    common/plt1/cmdstr
    common/cnvcom/conver
    common/plt2/fac
    common/io/fileout,inew
    common/kkplot/szrat
    common/chpcom/ientry,prtrt
    common/fntcom/curfnt,ifntsz,nfont
    common/outcom/iunit
    common/pagcom/npage
    !data ioffp,spval/0,0.0/
    ioffp=0
    spval=0.0

    !Set conversion factor (conver=72. for inches, conver=72./25.4 for mm, etc.)
    !conver
    conver=72.

    npage=1

    prtrt=portrait

    first=.true.
    pi=4.*abs(atan(1.))


    !Use default name unless newdev has already been called (inew=999).
    if(inew.eq.0) then
       fileout='psplot.ps'
       inew=1
    else if(inew.eq.999)then
       inew=1
    endif

    !output file opened
    iunit=fileNbr

    cmdstr='%!PS-Adobe-3.0'
    if (g_psp_isEpsFile) cmdstr='%!PS-Adobe-3.0 EPSF-3.0'
    call filler

    cmdstr= '%%Title: '//fileout(1:lenstr(fileout,80))

    call filler
    call date_and_time(dat,tim,zn,val)
    if(timer(1).eq.' ')timer(1)='0'
    if(dater(1).eq.' ')dater(1)='0'
    cmdstr= '%%CreationDate: '//DAT//' '//TIM
    call filler

    cmdstr= '%%Creator: PSPLOT PostScript Plotting Package'
    call filler

    if (g_psp_isEpsFile) then
       cmdstr='%%BoundingBox: 0 0 595 842'
       call filler
    end if

    cmdstr= '%%Pages: (atend)'
    call filler

    cmdstr='%%EndComments'
    call filler

    cmdstr='%Library Creator: Kevin E. Kohler <kevin@ocean.nova.edu>'
    call filler

    cmdstr='%%BeginProlog'
    call filler

    cmdstr='/inch {72 mul} bind def'
    call filler

    cmdstr='/Ah {moveto lineto lineto stroke} def'
    call filler

    cmdstr='/Ar {moveto 2 copy lineto 4 -2 roll'
    call filler
    cmdstr='     moveto lineto lineto stroke } def'
    call filler

    cmdstr='/arcit {S /A2 exch def /A1 exch def /Rad exch def'
    call filler
    cmdstr='        /Yc exch def /Xc exch def'
    call filler
    cmdstr='        Xc Rad A1 cos mul add Yc Rad A1 sin mul add'
    call filler
    cmdstr='        moveto newpath'
    call filler
    cmdstr='        Xc Yc Rad A1 A2 arc stroke} def'
    call filler

    cmdstr='/C {/Rad exch def /Yc exch def /Xc exch def'
    call filler
    cmdstr='         Xc Yc Rad 0 360 arc closepath'    !closepath needed to
    !avoid notch
    call filler                                        !with fat line width
    cmdstr='        } def'
    call filler

    cmdstr='/c0sf {closepath 0 setgray fill} def'
    call filler

    cmdstr='/cf {closepath fill} def'
    call filler

    cmdstr='/Cs {closepath stroke} def'
    call filler

    cmdstr='/Cln {newpath 3 1 roll'
    call filler
    cmdstr='      moveto {lineto} repeat clip newpath'
    call filler
    cmdstr='     } def'
    call filler

    cmdstr='/Cs {closepath stroke} def'
    call filler

    cmdstr='/Fb {newpath moveto '
    call filler
    cmdstr=' Dx 0 rlineto 0 Dy rlineto Dx neg 0 rlineto closepath'
    call filler
    cmdstr=' fill } def'
    call filler

    cmdstr='/Fbn { newpath 3 1 roll moveto {lineto} repeat'
    call filler
    cmdstr='       closepath fill } def'
    call filler

    cmdstr='/Fbnc { newpath 3 1 roll moveto'
    call filler
    cmdstr='       {lineto} repeat closepath fill } def'
    call filler

    cmdstr='/L /lineto load def'
    call filler

    cmdstr='/Lend {/Strlen exch stringwidth pop def} def'
    call filler

    !Define stringlength slightly increased for integrand placement
    cmdstr='/Lendi {/Strlen exch stringwidth pop 1.5 mul def} def'
    call filler

    !Define stringlength slightly increased for summation placement
    cmdstr='/Lends {/Strlen exch stringwidth pop 1.1 mul def} def'
    call filler

    cmdstr='/Lenssd {/Strlenss exch stringwidth pop 3 mul 4 div def} def'
    call filler

    cmdstr='/LSM {2 copy lineto stroke moveto} def'
    call filler

    cmdstr='/lsm {Xp Yp lineto stroke mover} def'
    call filler

    cmdstr='/M /moveto load def'
    call filler

    cmdstr='/mover {Xp Yp moveto} def'
    call filler

    cmdstr='/Np {newpath} def'
    call filler

    cmdstr='/S /stroke load def'
    call filler

    cmdstr='/Sc {setrgbcolor} def'
    call filler

    cmdstr='/Sg {setgray} def'
    call filler

    cmdstr='/Setf {Curfnt exch scalefont setfont} def'
    call filler

    cmdstr='/SM {stroke moveto} def'
    call filler

    cmdstr='/sm {stroke mover} def'
    call filler

    write(cmdstr,'(''/Slw {'',f7.4,'' mul setlinewidth} def'')') conver
    call filler

    cmdstr='/Slw0 {.24 setlinewidth} bind def'  !Minimum line width 300 dpi
    call filler

    !Add this for fun
    cmdstr= '%Line Breaking Procedure'
    call filler

    cmdstr='/TurnLineFL'
    call filler
    cmdstr='   { /T exch def /spacewidth space stringwidth pop def'
    call filler
    cmdstr='     /currentw 0 def /wordspace_count 0 def'
    call filler
    cmdstr='     /restart 0 def  /remainder T def'
    call filler
    cmdstr='     {remainder space search'
    call filler
    cmdstr='       {/nextword exch def pop'
    call filler
    cmdstr='        /remainder exch def'
    call filler
    cmdstr='        /nextwordwidth nextword stringwidth pop def'
    call filler
    cmdstr='        currentw nextwordwidth add lw gt'
    call filler
    cmdstr='        {T restart wordspace_count restart sub'
    call filler
    cmdstr='         getinterval showline'
    call filler
    cmdstr='         /restart wordspace_count def'
    call filler
    cmdstr='         /currentw nextwordwidth spacewidth add def'
    call filler
    cmdstr='        }'
    call filler
    cmdstr='        {/currentw currentw nextwordwidth add'
    call filler
    cmdstr='         spacewidth add def'
    call filler
    cmdstr='        } '
    call filler
    cmdstr='        ifelse'
    call filler
    cmdstr='        /wordspace_count wordspace_count'
    call filler
    cmdstr='        nextword length add 1 add def'
    call filler
    cmdstr='       }'
    call filler
    cmdstr='       {pop exit}'
    call filler
    cmdstr='       ifelse'
    call filler
    cmdstr='     } loop'
    call filler
    cmdstr='     /lrem remainder stringwidth pop def'
    call filler
    cmdstr='     currentw lrem add lw gt'
    call filler
    cmdstr='     {T restart wordspace_count restart sub '
    call filler
    cmdstr='      getinterval showline remainder showline}'
    call filler
    cmdstr='     {/lastchar T length def'
    call filler
    cmdstr='      T restart lastchar restart sub getinterval '
    call filler
    cmdstr='      lm y moveto show}'
    call filler
    cmdstr='     ifelse'
    call filler
    cmdstr='   } def'
    call filler

    cmdstr=' /parms {/y exch def /lm exch def /rm exch def'
    call filler
    cmdstr='         /leading exch def /pointsize exch def'
    call filler
    cmdstr='         /lw rm lm sub def'
    call filler
    cmdstr='         findfont pointsize scalefont setfont '
    call filler
    cmdstr='         /showline {lm y moveto show'
    call filler
    cmdstr='         /y y leading sub def} def'
    call filler
    cmdstr='         lm y moveto } def'
    call filler

    cmdstr='/Xposd {/Xpos exch def} def'
    call filler

    cmdstr='/Xposjd  {/Xpos exch Xpos exch Strlen mul sub def} def'
    call filler

    cmdstr='/xydef {/Xp Xpos def /Yp Ypos def} def'
    call filler

    cmdstr='%/Xypd {/Yp exch def /Xp exch def} def'
    call filler

    cmdstr='/Xypos0d {/Xpos0 Xpres def /Ypos0 Ypres def} def'
    call filler

    cmdstr='/Xyprset {dup /Xpres exch cos Strlen mul Xpos add def'
    call filler
    cmdstr='              /Ypres exch sin Strlen mul Ypos add def} def'
    call filler

    cmdstr='/Xyprset0 {dup /Xpres exch cos Strlen mul Xpos0 add def'
    call filler
    cmdstr='               /Ypres exch sin Strlen mul Ypos0 add def} def'
    call filler

    cmdstr='/Yposd {/Ypos exch def} def'
    call filler

    cmdstr='/Yposjd  {/Ypos exch Ypos exch Strlen mul sub def} def'
    call filler

    cmdstr='/space ( ) def'
    call filler

    cmdstr='%%EndProlog'
    call filler

    cmdstr='%%Page: 1 1'
    call filler

    !Szrat is the ratio of width to height of characters. Determined empirically.
    szrat=.6
    !Set initial font to helvetica, 12 point
    ifntsz=12
    call setfnt(20)
    !Set factor to 1 for initialization, reset later if chopit called
    fac=1.
    call factor(fac)

    fact=min(595./(g_psp_bBoxXmax-g_psp_bBoxXmin), &
         842./(g_psp_bBoxYmax-g_psp_bBoxYmin))/72.
    write(cmdstr,'(2f8.3,a)') fact,fact,' scale'
    call filler

    write(cmdstr,'(2f10.2,a)') -g_psp_bBoxXmin*72., &
         -g_psp_bBoxYmin*72., &
         ' translate'
    call filler

    !Set initial lineweight to 0
    call setlw(0.)
    !Set initial grayscale to 0
    call setgry(0.)
    !Set initial rgb colors to black(0)
    call setcolr(0.,0.,0.)

    xsh=0.
    ysh=0.
    call plot(xsh,ysh,-3)
  end subroutine psinit

  subroutine setcolr(red,green,blue)
    !this routines sets the current color
    !red, green blue are the saturation ratios between 0 and 1
    character*132 cmdstr
    common/plt1/cmdstr
    common/colrcom/cred,cgreen,cblue,cgry

    r=red
    r=amin1(1.,r)
    r=amax1(0.,r)
    g=green
    g=amin1(1.,g)
    g=amax1(0.,g)
    b=blue
    b=amin1(1.,b)
    b=amax1(0.,b)

    cmdstr=' '
    write(cmdstr,'(3F7.3,'' Sc'')')r,g,b
    call filler
    cred=r
    cgreen=g
    cblue=b
  end subroutine setcolr

  subroutine setfnt(numfnt)
    !This routines changes the typeface of the current font
    character*132 cmdstr,scrc
    character*132 curfnt
    common/fntcom/curfnt,ifntsz,nfont
    character*40 fntnam(35)
    common/plt1/cmdstr
    data fntnam/'AvantGarde-Book','AvantGarde-BookOblique',&
         &'AvantGarde-Demi','AvantGarde-DemiOblique','Bookman-Demi',&
         &'Bookman-DemiItalic','Bookman-Light','Bookman-LightItalic',&
         &'Courier-Bold','Courier-BoldOblique','Courier-Oblique', 'Courier',&
         &'Helvetica-Bold','Helvetica-BoldOblique', 'Helvetica-Narrow-Bold',&
         &'Helvetica-Narrow-BoldOblique', 'Helvetica-Narrow-Oblique',&
         &'Helvetica-Narrow', 'Helvetica-Oblique','Helvetica',&
         &'NewCenturySchlbk-Bold','NewCenturySchlbk-BoldItalic',&
         &'NewCenturySchlbk-Italic','NewCenturySchlbk-Roman',&
         &'Palatino-Bold','Palatino-BoldItalic','Palatino-Italic',&
         &'Palatino-Roman','Symbol','Times-Bold','Times-BoldItalic',&
         &'Times-Italic','Times-Roman','ZapfChancery-MediumItalic',&
         &'ZapfDingbats'/


    nfont=numfnt
    if(numfnt.lt.1.or.numfnt.gt.35) then
       print *,'Invalid font number encountered in **setfnt**'
       print *,'Using Helvetica default'
       nfont=20
    endif
    scrc=fntnam(nfont)
    cmdstr='/Curfnt /'//scrc(1:lenstr(scrc,132))//' findfont def'
    call filler
    write(cmdstr,'(i3,'' Setf'')')ifntsz
    call filler
  end subroutine setfnt

  subroutine setgry(gry)
    !This routines sets the current gray level
    !Gry is set to be between 0 and 1
    character*132 cmdstr
    common/plt1/cmdstr
    common/colrcom/cred,cgreen,cblue,cgry

    g=gry
    g=amin1(1.,g)
    g=amax1(0.,g)

    cmdstr=' '
    write(cmdstr,'(F7.3,'' Sg'')')g
    call filler
    cgry=g
  end subroutine setgry

  subroutine setlw(rlwi)
    !this routines sets the current linewidth
    !rlwi is linewidth in inches
    character*132 cmdstr
    common/plt1/cmdstr
    common/lcom/curwid

    if(abs(rlwi).lt.1.e-5) then  !0
       cmdstr='Slw0'
    else
       cmdstr=' '
       write(cmdstr,'(F7.3,'' Slw'')')rlwi
    endif
    call filler
    curwid=rlwi
  end subroutine setlw

  subroutine factor(facc)
    common/plt2/fac
    character*132 cmdstr
    common/plt1/cmdstr

    !Unscale previous scaling
    recipx=1./fac
    recipy=1./fac
    write(cmdstr,'(2f7.3,a)')recipx,recipy,' scale'
    call filler
    fac=facc
    write(cmdstr(:14),'(2f7.3)')fac,fac
    call filler
  end subroutine factor

  subroutine plot(xcall,ycall,ip)
    character*132 cmdstr
    character*80 scr
    common/plt1/cmdstr
    common/cnvcom/conver
    common/outcom/iunit
    common/pagcom/npage

    ipp=iabs(ip)

    if(ip.eq.999) then   !Terminate plot session.
       cmdstr='stroke showpage'
       call filler

       cmdstr='%%Trailer'
       call filler

       write(scr,'(i6)')npage
       call blkstp(scr,80,scr,nch)
       cmdstr='%%Pages: '//scr(1:nch)
       call filler

       cmdstr='%%EOF'
       call filler

       return
    endif

    !Moving pen
    if(ipp.eq.3) then    !Stroke to paint previous path, then moveto
       write(cmdstr,'(2f10.2,'' SM'')')xcall*conver,ycall*conver
    else                 !Lineto
       write(cmdstr,'(2f10.2,'' LSM'')')xcall*conver,ycall*conver
    endif
    call filler

    !Reset origin if ip.lt.0
    if(ip.lt.0) then
       write(cmdstr,'(2f10.2,'' translate'')')xcall*conver,ycall*conver
       call filler
       ipen=ipp
    endif

  end subroutine plot

  subroutine plotnd
    call plot(0.,0.,999)
  end subroutine plotnd

  subroutine circle(xc,yc,rad,fill)
    character*132 cmdstr,scrc
    common/plt1/cmdstr
    common/cnvcom/conver
    logical fill
    xci=xc*conver
    yci=yc*conver
    radi=rad*conver
    scrc=' '
    write(scrc,'(f10.2,'' '',f10.2,'' '',f10.2,'' C'')') xci,yci,radi
    if(fill) then
       cmdstr='Np '//scrc(1:lenstr(scrc,132))//' fill'
    else
       cmdstr='Np '//scrc(1:lenstr(scrc,132))//' stroke'
    endif
    call filler
  end subroutine circle

  function lenstr(string,ls) 
    !This routine finds actual length of string by eliminating trailing blanks
    character*(*) string

    do i=ls,1,-1
       is=i
       if(string(i:i).ne.char(32)) goto 10
    enddo
    is=0
10  lenstr=is
  end function lenstr

  subroutine blkstp(ch,ndim,a,leng)
    !character*1 ch(ndim),a(ndim)
    character(len=*) ch,a
    !Strip out blanks only (leave in esc, etc.)
    i=1
    leng=0
10  continue
    if(ichar(ch(i:i)).ne.32)then
       leng=leng+1
       a(leng:leng)=ch(i:i)
    endif

    if(i.eq.ndim) then
       !Blankfill remainder of output array
       do l=leng+1,ndim
          a(l:l)=' '
       enddo
       return
    endif

    i=i+1
    goto 10
  end subroutine blkstp

  subroutine keknum(xp,yp,size,fpn,ang,ndec,mjus)
    !Just assume that user really wants kekflt.
    call kekflt(xp,yp,size,fpn,ang,ndec,mjus)
  end subroutine keknum

  subroutine kekflt(xp,yp,size,fpn,ang,ndec,mjus)
    dimension ichrnum(20)

    fnum=fpn
    !Get number in character form
    call numsym(fnum,ndec,ichrnum,ndigit,.false.)
    call keksym(xp,yp,size,ichrnum,ang,ndigit,mjus)
  end subroutine kekflt

  subroutine keksym(xp,yp,size,ltitle1,ang,nchar1,mjus)
    character*132 cmdstr
    character*132 curfnt,scrc
    common/fntcom/curfnt,ifntsz,nfont
    character*80 titlec,titleb
    character*1 bslash
    dimension ltitle(20),ltitle1(20)
    equivalence(ititle,ltitle(1))
    common/plt1/cmdstr
    common/cnvcom/conver
    common/kkplot/szrat

    !Stroke previous paths before this write
    cmdstr='S'
    call filler

    bslash=char(92)

    pi=4*abs(atan(1.))

    if(nchar1.eq.-999) then !octal code
       do n=1,20
          ltitle(n)=ltitle1(n)
       enddo
       nchar=1
    else
       nchar=nchar1
       write(titlec,'(20a4)')ltitle1
       if(iabs(nchar).lt.80)titlec(nchar+1:80)=' '
       read(titlec,'(20a4)')ltitle
    endif

    !Choose proper font height, using current font
    mchar=iabs(nchar)
    !Set character size
    iht=max(1,int(size*conver/.6))     !.6 FACTOR IS EMPIRICAL

    if(iht.ne.ifntsz) then
       cmdstr=' '
       write(cmdstr,'(I3,'' Setf'')')iht
       call filler
       ifntsz=iht
    endif

    if(nchar1.eq.-999) then        !Octal code
       write(titlec,'(A1,I10)')bslash,ititle
       call blkstp(titlec,80,titlec,numc)
    else
       write(titlec,'(20a4)')ltitle
       !Check if titlec contains ( or ) or \.  These characters must be treated
       !specially by preceding them with a "\".  Do this to ( and ) even though
       !they might be balanced, i.e. () within a string, which can be treated
       !normally.

       titleb=titlec
       numc=0
       do m=1,mchar
          if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')' .or. &
               & titleb(m:m).eq.bslash) then
             numc=numc+1
             titlec(numc:numc)=bslash
          endif
          numc=numc+1
          titlec(numc:numc)=titleb(m:m)
       enddo
    endif

    mchar=numc
    xpos=xp
    ypos=yp
    if(nchar.lt.0) then
       njus=0
    else
       njus=mjus
    endif
    rsize=size
    !Character space height is 2.0 x char height
    !Character space width is 1.5 x char width
    !Actual string length is (nc-1)*1.5*char width + char width
    strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat

    if(xpos.eq.999.) then
       cmdstr='/Xpos Xpres def'
       njus=0
    else
       cmdstr=' '
       write(cmdstr,'(f10.2,'' Xposd'')')xp*conver
    endif
    call filler

    if(ypos.eq.999.) then
       cmdstr='/Ypos Ypres def'
       njus=0
    else
       cmdstr=' '
       write(cmdstr,'(f10.2,'' Yposd'')')yp*conver
    endif
    call filler

    if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
       print 110, njus
110    format(1x,'incorrect justification code ',i5,'found in ',&
            &'KEKSYM, zero used')
       njus=0
    endif
    !Strlen has already been "factored" by the choice of font height
    !Since it will eventually be factored again, we must divide by
    !factor now.
    cmdstr='('//titlec(1:mchar)//') Lend'
    arg=ang*4.*abs(atan(1.))/180.
    xarg=cos(arg)*njus/2.
    yarg=sin(arg)*njus/2.

    if(xarg.ne.0.) then
       scrc=' '
       write(scrc,'(f7.3,'' Xposjd'')')xarg
       cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
    endif

    if(yarg.ne.0.) then
       scrc=' '
       write(scrc,'(f7.3,'' Yposjd'')')yarg
       cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
    endif

    if (nchar.eq.-1) then    !centered symbol
       high=rsize
       wide=high*szrat

       xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
       ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

       xarg=high/2.*sin(arg)-wide/2.*cos(arg)
       yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

       if(xarg.ne.0.) then
          scrc=' '
          write(scrc, '(''/Xpos Xpos'',f10.2,'' add def'')')xarg*conver
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
       endif

       if(yarg.ne.0.) then
          scrc=' '
          write(scrc,'(''/Ypos Ypos'',f10.2,'' add def'')')yarg*conver
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
       endif
    endif

    call filler

    !Move pen to proper coordinates
    scrc='xydef mover'
    cmdstr=scrc(1:lenstr(scrc,132))

    !Set angle
    if(ang.ne.0.) then
       scrc=' '
       write(scrc,'(F7.1,'' rotate'')') ang
       cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
    endif

    scrc='('//titlec(1:mchar)//') show'
    cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

    !Reset angle
    if(ang.ne.0.) then
       scrc=' '
       write(scrc,'(F7.1,'' rotate'')') -ang
       cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
    endif
    call filler

    !Start next char at .5 char width away
    argdeg=arg*180./pi
    cmdstr=' '
    write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
    call filler
  end subroutine keksym

  subroutine numsym(fpn,ndec,itext,nchar,eform) 
    !eform:  true for exponential format
    !false for floating pt format
    logical eform
    dimension itext(20)
    character*10 ifrmt
    character*80 a
    a=' '

    !Check if ndec is valid
    if(ndec.gt.15) then
       print 100, ndec
100    format(1x,'In call to numsym, ndec gt 15 ',i20,' program',' abandoned')
       stop
    else if(eform.and.ndec.lt.0) then
       print 110, ndec
110    format(1x,'in numsym, exponential format specified with ','ndec= '&
            &,i3,' program abandoned')
       stop
    endif

    if(eform) then
       write(ifrmt,'(''(1pe16.'',i2,'')'')')ndec
    else if(ndec.lt.0) then
       ifrmt='(f16.1)'
    else
       write(ifrmt,'(''(f16.'',i2,'')'')')ndec
    endif

    write(a,ifrmt)fpn

    !Strip off all blanks in a
    call blkstp(a,80,a,nchar)
    ipos=index(a,'.')
    if(eform) then
       if(ndec.eq.0) then
          !Delete characters between '.' and 'E'
          ie=index(a,'E')
          do n=ie,nchar
             nind=ipos+n-ie+1
             a(nind:nind)=a(n:n)
          enddo
          nchar=nchar-(ie-ipos-1)
       endif
    else if(ndec.lt.0) then
       nchar=ipos-1
    else if(ndec.eq.0) then
       nchar=ipos
    endif
    read(a,'(20a4)')itext
  end subroutine numsym

end module derivedPSPLOT
