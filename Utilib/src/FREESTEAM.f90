!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for freesteam (light water).
!
!Copyright:
! Copyright (C) 2012 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): A. Hebert
!
!-----------------------------------------------------------------------
!
subroutine THMSAP(p, t)
! return the saturation pressure (Pa) as a function of the temperature (K)
   use, intrinsic :: iso_c_binding
   real :: p, t
   real(c_double) :: td
   interface 
      real(c_double) function region4_psat_T (td) bind(c, name='freesteam_region4_psat_T')
      use, intrinsic :: iso_c_binding
      real(c_double), value :: td
      end function region4_psat_T
   end interface
   td=t
   p = real(region4_psat_T (td))
end subroutine THMSAP
!
subroutine THMSAT(p, t)
   ! return the saturation temperature (K) as a function of the pressure (Pa)
   use, intrinsic :: iso_c_binding
   real :: p, t
   real(c_double) :: pd
   interface 
      real(c_double) function region4_Tsat_p (pd) bind(c, name='freesteam_region4_Tsat_p')
      use, intrinsic :: iso_c_binding
      real(c_double), value :: pd
      end function region4_tsat_p
   end interface
   pd=p
   t = real(region4_tsat_p (pd))
end subroutine THMSAT
!
subroutine THMPT(p, t, rho, h, zk, zmu, cp)
   ! return the remaining thermohydraulics parameters as a function of the pressure (Pa)
   ! and temperature (K)
   use, intrinsic :: iso_c_binding
   real :: p, t, rho, h, zk, zmu, cp
   real(c_double) :: pd, td, rhod, hd, zkd, zmud, cpd
   interface 
      subroutine free_pT (pd, td, rhod, hd, zkd, zmud, cpd) bind(c, name='free_pT')
      use, intrinsic :: iso_c_binding
      real(c_double) :: pd, td, rhod, hd, zkd, zmud, cpd
      end subroutine free_pT
   end interface
   pd=p
   td=t
   call free_pT(pd, td, rhod, hd, zkd, zmud, cpd)
   rho=real(rhod)
   h=real(hd)
   zk=real(zkd)
   zmu=real(zmud)
   cp=real(cpd)
end subroutine THMPT
!
subroutine THMTX(t, x, rho, h, zk, zmu, cp)
   ! return the remaining thermohydraulics parameters as a function of the temperature (K)
   ! and quality
   use, intrinsic :: iso_c_binding
   real :: t, x, rho, h, zmu
   real(c_double) :: td, xd, rhod, hd, zkd, zmud, cpd
   interface 
      subroutine free_Tx (td, xd, rhod, hd, zkd, zmud, cpd) bind(c, name='free_Tx')
      use, intrinsic :: iso_c_binding
      real(c_double) :: td, xd, rhod, hd, zkd, zmud, cpd
      end subroutine free_Tx
   end interface
   td=t
   xd=x
   call free_Tx(td, xd, rhod, hd, zkd, zmud, cpd)
   rho=real(rhod)
   h=real(hd)
   zk=real(zkd)
   zmu=real(zmud)
   cp=real(cpd)
end subroutine THMTX
