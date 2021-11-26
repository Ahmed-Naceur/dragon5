!
!-----------------------------------------------------------------------
!
!Purpose:
! Definition of parameter types used in Dragon.
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
!-----------------------------------------------------------------------
!
module constType
  implicit none

  !type of geometry (ST(1))
  integer,parameter  :: G_Virtual=0
  integer,parameter  :: G_Homoge=1
  integer,parameter  :: G_Car1d=2
  integer,parameter  :: G_Tube=3
  integer,parameter  :: G_Car2d=5
  integer,parameter  :: G_Hex=8
  integer,parameter  :: G_Tri=13
  integer,parameter  :: G_Carcel=20
  integer,parameter  :: G_Hexcel=24
  integer,parameter  :: G_Group=30

  !type of sectorisation (ST(14))
  integer,parameter  :: S_X_tot=-1
  integer,parameter  :: S_not=0
  integer,parameter  :: S_T_tot=1  !"T" = "+"
  integer,parameter  :: S_TX_tot=2
  integer,parameter  :: S_TXS_tot=3
  integer,parameter  :: S_WM_tot=4

  !type of hexagonal symetry (IHEX)
  integer,parameter  :: H_S30=1
  integer,parameter  :: H_SA60=2
  integer,parameter  :: H_SB60=3
  integer,parameter  :: H_S90=4
  integer,parameter  :: H_R120=5
  integer,parameter  :: H_R180=6
  integer,parameter  :: H_SA180=7
  integer,parameter  :: H_SB180=8
  integer,parameter  :: H_Complete=9

  !type of triangular symetry (ITRI)
  integer,parameter  :: T_S30=1
  integer,parameter  :: T_SA60=2
  integer,parameter  :: T_ST60=3
  integer,parameter  :: T_Complete=4

  !type of boundary condition (NCODE)
  integer,parameter  :: B_NotUsed=0
  integer,parameter  :: B_Void=1
  integer,parameter  :: B_Refl=2
  integer,parameter  :: B_Diag=3
  integer,parameter  :: B_Tran=4
  integer,parameter  :: B_Syme=5
  integer,parameter  :: B_Albe=6
  integer,parameter  :: B_Zero=7
  integer,parameter  :: B_Pi_2=8
  integer,parameter  :: B_Pi=9
  integer,parameter  :: B_Ssym=10

  !type of rotation (TURN)
  integer,parameter  :: R_A=1
  integer,parameter  :: R_B=2
  integer,parameter  :: R_C=3
  integer,parameter  :: R_D=4
  integer,parameter  :: R_E=5
  integer,parameter  :: R_F=6
  integer,parameter  :: R_G=7
  integer,parameter  :: R_H=8  !max for cartesian
  integer,parameter  :: R_I=9
  integer,parameter  :: R_J=10
  integer,parameter  :: R_K=11
  integer,parameter  :: R_L=12 !max for hex and tri

  ! 
  type t_point
     double precision :: x
     double precision :: y
  end type t_point
end module constType
