!
!-----------------------------------------------------------------------
!
!Purpose:
! Definition of parameter types used in G2S: module.
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
! De plus, quelques fonctions utilitaires sont aussi definies, telles que :
!  - isEqual : teste l'egalite a epsilon pres de deux doubles, et les egalise
!              et si leur valeur est proche, les egalise a la moyenne des deux
!  - isEqualConst : teste l'egalite a epsilon pres de deux doubles
!  - isEqualAngl : teste l'egalite a epsilon pres de deux angles en les
!                  normalisant si besoin
!  - isEqualConstAngl : teste l'egalite a epsilon pres de deux angles normaux
!  - angleNormal : retourne la valeur de l'angle centree entre -Pi et Pi
!  - distance : calcule la distance entre un point et une droite, et calcule
!               les coordonnees de la projection du point sur la droite
!  - longVect : calcule la norme d'un vecteur
!  - calculeAngle : calcule l'angle d'un vecteur par rapport a Ox
!  - estColli :  teste la colinearite de deux vecteurs
!  - estAligne : teste l'alignement de trois points
!  - estAGauche : teste si un point est a gauche d'une droite
!  - estAGaucheStrict : teste si un point est strictement a gauche d'une droite
!  - estADroite : teste si un point est a droite d'une droite
!  - estADroiteStrict : teste si un point est strictement a droite d'une droite
!  - isPi : teste l'egalite avec Pi d'un angle
!  - estPaire : teste la parite d'un entier
!  - isAngleInInterval : teste l'appartenance d'un angle a un domaine angulaire
!
!-----------------------------------------------------------------------
!
module constUtiles
  implicit none

  double precision,parameter  :: epsilon  = 1.d-4
  double precision,parameter  :: muleps1  = 5.d0
  double precision,parameter  :: muleps2  = 2.d0
  double precision,parameter  :: pi_c     = 3.14159265358979d0
  double precision,parameter  :: pi_2_c   = 1.57079632679490d0
  double precision,parameter  :: pi_3_c   = 1.04719755119660d0
  double precision,parameter  :: dpi_c    = 6.28318530717959d0  !=2.*pi
  double precision,parameter  :: rad2deg  = 5.72957795130823d1  !=180./pi
  double precision,parameter  :: deg2rad  = 1.74532925199433d-2 !=pi/180.
  double precision,parameter  :: infinity = 1.d99
  double precision,parameter  :: sqrt3_2d = 8.66025403784439d-1 !=sqrt(3)/2
  character*5,parameter       :: formatr = 'e18.7'
  character*5,parameter       :: formati = 'i6'
  real,parameter              :: gSALeps = 1e-4

  ! Numerical Constants
  double precision, parameter :: dp_0 = 0.0d0
  double precision, parameter :: dp_1 = 1.0d0
  double precision, parameter :: dp_05 = 0.5d0

contains

  logical function isEqual(a,b)
    double precision,intent(inout) :: a,b
    isEqual = (abs(b-a)<epsilon)
    if (isEqual) then ; a = (a+b)*0.5d0 ; b = a ; end if
  end function isEqual

  logical function isEqualAngl(a,b)
    double precision,intent(inout) :: a,b
    a=angleNormal(a) ; b=angleNormal(b) ; isEqualAngl = isEqual(a,b) 
!CS-IV    isEqualAngl = isEqual(a,b) .or. isEqualConst(a+dpi_c,b) &
!CS-IV         .or. isEqualConst(a,b+dpi_c)
  end function isEqualAngl

  logical function isEqualAnglConst(a,b)
    double precision,intent(inout) :: a,b
    a=angleNormal(a) ; b=angleNormal(b) ; isEqualAnglConst = isEqualConst(a,b) 
  end function isEqualAnglConst

  ! CS-IV : fonction a supprimer, elle n'apporte rien
  function isEqualConstAngl(a,b)
    double precision,intent(in) :: a,b
    logical                     :: isEqualConstAngl
    double precision :: aa, bb
    aa = a ; bb = b
    isEqualConstAngl = isEqualAngl(aa,bb)
  end function isEqualConstAngl

  double precision function angleNormal(a)
    double precision,intent(in) :: a
    angleNormal = modulo(a,dpi_c)
    if (angleNormal>pi_c) angleNormal = angleNormal - dpi_c
  end function angleNormal

  logical function isEqualConst(a,b)
    double precision,intent(in) :: a,b
    isEqualConst = (abs(b-a)<epsilon)
  end function isEqualConst

  function distance(ptx,pty,dtox,dtoy,dtex,dtey,prjx,prjy)
    double precision,intent(in)  :: ptx,pty,dtox,dtoy,dtex,dtey
    double precision,intent(out) :: prjx,prjy
    double precision             :: distance
    
    double precision :: A,B,C,invn
    
    !retourne |D(pt)| avec D:ax+by+c=0 ou a**2+b**2=1
    !et le projete de pt sur D prj
    !(la droite D est (O,E) ie origine,extremite)
    ! alors a=A/n , b=B/n et c=C/n avec A=dtoy-dtey , B=dtex-dtox ,
    ! C=dtox*dtey-dtex*dtoy et n=SQRT(A**2+B**2)
    A=dtoy-dtey ; B=dtex-dtox ; C=dtox*dtey-dtex*dtoy
    invn = 1.d0/sqrt(A*A+B*B)
    A=A*invn ; B=B*invn ; C=C*invn
    distance = abs(A*ptx+B*pty+C)
    prjx = B*B*ptx-A*B*pty-A*C
    prjy = A*A*pty-A*B*ptx-B*C
  end function distance

  function longVect(vx,vy)
    double precision,intent(in) :: vx,vy
    double precision            :: longVect

    longVect = sqrt(vx*vx+vy*vy)
  end function longVect

!!$  double precision function getAngle(u,v)
!!$    ! CS-IV : Calculate the angle between the vector and the Ox axis without 
!!$    ! seeking to bring this value to a usual characteristic value.
!!$    double precision, intent(in) :: u, v
!!$    double precision :: norm
!!$    norme = sqrt(u*u+v*v)
!!$    if (IsEqualConst(norme, dble_zero)) then
!!$       getAngle = dble_zero
!!$    else
!!$       getAngle = acos(u/norm)
!!$    endif
!!$  end function getAngle
 
!!$  double precision function getAngleVect(u1,u2,v1,v2)
!!$    !  CS-IV : Calculate the angle between the two vectors without seeking 
!!$    ! to bring this value to a usual characteristic value.
!!$    double precision, intent(in) :: u1, u2, v1, v2
!!$    double precision:: normU, normV, scalProd
!!$    
!!$    normU = sqrt(u1*u1+u2*u2)
!!$    normV = sqrt(v1*v1+v2*v2)
!!$    if (isEqualConst(normU,dble_zero).or.isEqualConst(normV,dble_zero)) then
!!$       getAngleVect = dble_zero
!!$    else
!!$       scalProd = u1*v1+v2*v2
!!$       getAngleVect = acos(scalProd/(normU*normV))
!!$    end if
!!$
!!$  end function getAngleVect
    
 function calculeAngle(cx,cy,px,py)
    double precision,intent(in) :: cx,cy,px,py
    double precision            :: calculeAngle

    double precision :: u,v

    u=px-cx ; v=py-cy
    if (isEqualConst(u,0.d0).and.isEqualConst(v,0.d0)) then
       calculeAngle=0.d0
       return
    end if
    calculeAngle = acos(u/sqrt(u*u+v*v))
    if (isEqualConst(calculeAngle,0.d0)) then
       calculeAngle=0.d0
       return
    else if (isEqualConst(calculeAngle,pi_2_c)) then
       calculeAngle=pi_2_c
    else if (isEqualConst(calculeAngle,pi_c)) then
       calculeAngle=pi_c
       return
    end if
    if (py<cy) calculeAngle = -calculeAngle
  end function calculeAngle

  function estColi(ax,ay,bx,by)
    double precision,intent(in) :: ax,ay,bx,by
    logical                     :: estColi
    double precision :: la, lb, anx, any, bnx, bny

    la = longvect(ax,ay)
    lb = longvect(bx,by)
    if (isequalconst(la,0.0d0).or.isequalconst(lb,0.0d0)) then
       estColi=.true.
    else
       anx = ax/la ; any = ay/la
       bnx = bx/lb ; bny = by/lb
       estColi = isEqualConst(anx*bny,bnx*any)
    endif
  end function estColi

  function pointsAlignes(ax,ay,bx,by,cx,cy)
    double precision,intent(in) :: ax,ay,bx,by,cx,cy
    logical                     :: pointsAlignes

    PointsAlignes=(estColi(bx-ax,by-ay,cx-ax,cy-ay) &
         & .and. estColi(bx-ax,by-ay,cx-bx,cy-by))
  end function pointsAlignes

  function estAGauche(ox,oy,ex,ey,ptx,pty)
    double precision,intent(in) :: ox,oy,ex,ey,ptx,pty
    logical                     :: estAGauche

    double precision  :: tmp

    tmp = sin(calculeAngle(ox,oy,ptx,pty)-calculeAngle(ox,oy,ex,ey))
    estAGauche=(tmp>0.d0).or.pointsAlignes(ox,oy,ex,ey,ptx,pty)
  end function estAGauche

  function estAGaucheStrict(ox,oy,ex,ey,ptx,pty)
    double precision,intent(in) :: ox,oy,ex,ey,ptx,pty
    logical                     :: estAGaucheStrict

    double precision :: tmp

    tmp = sin(calculeAngle(ox,oy,ptx,pty)-calculeAngle(ox,oy,ex,ey))
    estAGaucheStrict=(tmp>0.d0).and.(.not.pointsAlignes(ox,oy,ex,ey,ptx,pty))
  end function estAGaucheStrict

  function estADroite(ox,oy,ex,ey,ptx,pty)
    double precision,intent(in) :: ox,oy,ex,ey,ptx,pty
    logical                     :: estADroite

    double precision :: tmp

    tmp = sin(calculeAngle(ox,oy,ptx,pty)-calculeAngle(ox,oy,ex,ey))
    estADroite=(tmp<0.d0).or.pointsAlignes(ox,oy,ex,ey,ptx,pty)
  end function estADroite

  function estADroiteStrict(ox,oy,ex,ey,ptx,pty)
    double precision,intent(in) :: ox,oy,ex,ey,ptx,pty
    logical                     :: estADroiteStrict

    double precision :: tmp

    tmp = sin(calculeAngle(ox,oy,ptx,pty)-calculeAngle(ox,oy,ex,ey))
    estADroiteStrict=(tmp<0.d0).and.(.not.pointsAlignes(ox,oy,ex,ey,ptx,pty))
  end function estADroiteStrict

  function isPi(angl)
    double precision,intent(in) :: angl
    logical                     :: isPi

    isPi=isEqualConst(abs(angl),pi_c)
  end function isPi

  logical function estPaire(n)
    integer,intent(in) :: n
    estPaire = ((n/2)*2==n)
  end function estPaire

  function isAngleInInterval(a,o,e)
    double precision,intent(in) :: a,o,e
    integer                     :: isAngleInInterval
    !dit si un angle est sur sur l'interval [o,e]:
    ! 0 -> pas dessus,
    ! 1 -> c'est l'origine,
    ! 2 -> entre les deux,
    ! 3 -> c'est l'extremite
    double precision :: aa,oo,ee
    aa = angleNormal(a) ; oo = angleNormal(o) ; ee = angleNormal(e)
    if (isEqualConstAngl(oo,aa)) then
       isAngleInInterval = 1
    else if (isEqualConstAngl(ee,aa)) then
       isAngleInInterval = 3 
    else if (oo<ee) then
       if((oo<aa).and.(aa<ee)) then
          isAngleInInterval = 2
       else
          isAngleInInterval = 0
       end if
    else
       if((ee<aa).and.(aa<oo)) then
          isAngleInInterval = 0
       else
          isAngleInInterval = 2
       end if
    end if
  end function isAngleInInterval

end module constUtiles
