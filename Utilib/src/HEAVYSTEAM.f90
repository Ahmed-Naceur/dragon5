!
!-----------------------------------------------------------------------
!
!Purpose:
! Heavy water properties
!
!Copyright:
! Copyright (C) 2018 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modIFy it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): G. Marleau and A. Hebert
!
!-----------------------------------------------------------------------
!
subroutine THMHSP(p, t)
! return the saturation pressure (Pa) as a function of the temperature (K)
! Ref: Ji. Zhang, January 20, 98
   real :: p, t
   double precision :: tcd,zd,pd
   character hsmg*131
   !
   tcd=dble(t-273.16)
   pd=0.0D0
   IF(tcd .LT. 90.5D0 .OR. tcd.GT.370.74D0) THEN
     WRITE(hsmg,*) 'THMHSP: T =',tcd,'C exceeds the valid temperature range', &
     & ' for automatic pressure evaluation (90.5<T<370.74)'
     call XABORT(hsmg)
   ELSEIF (tcd .GT. 320.0D0) THEN
     zd = (tcd - 320.0D0)/46.0D0
     pd = 11.442428D6 + (7.0289472D6 + 1.7422423D6*zd)*zd &
     + 0.2384117D6*zd**3.680D0
   ELSEIF (tcd .GT. 252.0D0) THEN
     zd = (tcd - 252.0D0)/68.0D0
     pd = 4.1333143D6 + (4.75435D6 + 2.124772D6*zd)*zd &
     + 0.4299918D6*zd**3.2250D0
   ELSEIF (tcd .GT. 179.0D0) THEN
     zd = (tcd - 179.0D0)/73.0D0
     pd = 0.9682855D6 + (1.6569005D6 + 1.1322316D6*zd)*zd &
     + 0.3758966D6*zd**3.1460D0
   ELSEIF (tcd .GT. 127.0D0) THEN
     zd = (tcd - 127.0D0)/52.0D0
     pd =0.2385152D6 + (0.384338812D6 + 0.2559459D6*zd)*zd &
     + 0.089485567D6*zd**3.1740D0
   ELSE
     zd = (tcd - 90.50D0)/36.50D0
     pd =0.06736076D6 + (0.09519366D6 + 0.05709235D6*zd)*zd &
     + 0.01886846D6*zd**3.2010D0
   ENDIF
   p =real(pd)
end subroutine THMHSP
!
subroutine THMHST(p, t)
! return the saturation temperature (K) as a function of the pressure (Pa)
! Ref: Ji. Zhang, January 20, 98
   real :: p, t
   double precision :: pd,z,td
   character hsmg*131
   !
   pd=p
   td=0.0d0
   IF (0.0741494D6 .GT. pd .OR. pd .GT. 21.2082144D6) THEN
     WRITE(hsmg,*) 'THMHST: P =',pd,'Pa exceeds the valid pressure range', &
     & ' for temperature evaluation (0.0741494<P<21.2082144) MPa'
     call XABORT(hsmg)
   ELSEIF (pd .GT. 6.9829216D6) THEN
       z = (pd - 6.9829216D6)/14.2252928D6
       td = 2.85D2 + (136.2916149D0 - 517.5749172D0*z)*z &
        + 465.2833023D0*z**2.0509579D0        
   ELSEIF (pd .GT. 2.2074482D6) THEN
       z = (pd - 2.2074482D6)/4.7754734D6
       td = 2.175D2 + (110.6079115D0 - 414.5685862D0*z)*z &
        + 371.4606747D0*z**2.0575065D0        
   ELSEIF (pd .GT. 0.6891798D6) THEN                      
       z = (pd - 0.6891798D6)/1.5182684D6
       td = 1.65D2 + (87.5356881D0 - 423.3068041D0*z)*z &
        + 388.271116D0*z**2.0455901D0        
   ELSEIF (pd .GT. 0.2241012D6) THEN                      
       z = (pd - 0.2241012D6)/0.4650786D6
       td = 1.25D2 + (66.1987487D0 - 225.285045D0*z)*z &
        + 199.0862962D0*z**2.0653628D0        
   ELSE                       
       z = (pd - 0.0741494D6)/0.1499518D6
       td = 9.3D1 + (53.0754191D0 - 165.0407938D0*z)*z &
        + 143.9653747D0*z**2.0723743D0
   ENDIF
   t=REAL(td)+273.16
end subroutine THMHST
!
subroutine THMHPT(p, t, rho, h, zk, zmu, cp)
! return the remaining thermohydraulics parameters as a function of the pressure (Pa)
! and temperature (K)
! Ref: Ji. Zhang, January 20, 98
   use, intrinsic :: iso_c_binding
   implicit real*8(a-h,o-z)
   real :: p, t, rho, h, zk, zmu, cp, ps
   interface 
      subroutine free_pT (pd, td, rhod, hd, zkd, zmud, cpd) bind(c, name='free_pT')
      use, intrinsic :: iso_c_binding
      real(c_double) :: pd, td, rhod, hd, zkd, zmud, cpd
      end subroutine free_pT
   end interface
   !
   tcd=dble(t-273.16)
   pd=dble(p)
   call THMHSP(ps, t)
   psd=dble(ps)
   !
   ! compute density in kg/m3
   drho=0.0d0
   IF (90.50D0 .GT. tcd .OR. tcd .GT. 369.0D0) THEN
     call XABORT('THMHPT: exceed the valid temperature range')
   ELSEIF (tcd .LT. 307.50D0) THEN
     zd = (tcd - 90.50D0)/217.0D0
     drho = 1.07065471D3 - (0.1611572D3 + 0.12188254D3*zd)*zd &
       - 0.2106382D2*zd**6.830
   ELSEIF (tcd .LT. 350.0D0) THEN
     zd = (tcd - 307.50D0)/42.50D0
     drho = 0.7665511D3 - (0.1074816D3 + 0.2506107D2*zd)*zd &
       - 0.1042465D2*zd**4.640
   ELSEIF (tcd .LT. 364.0D0) THEN
     zd = (tcd - 350.0D0)/14.0D0
     drho = 0.6235838D3 - (0.0678503D3 + 0.01551348D3*zd)*zd &
       - 0.7259508D1*zd**4.470
   ELSE   
     zd = (tcd - 364.0D0)/4.50D0
     drho = 0.5329606D3 - (0.4286569D2 + 0.1064107D2*zd)*zd &
       - 0.4761284D1*zd**5
   ENDIF
   A  = (0.91725D0+0.61D-7*(drho-825.0D0)**2)*(drho+103.65D0)-drho
   B  = 0.12D-6*psd-0.3D0
   DT = 1.0D6*(370.74D0-tcd)
   RR = (B+DT/(22.13D6 - psd))/(B+DT/(pd - psd))
   rho  = REAL((drho+(A*RR)))
   !
   ! compute specific enthalpy, in J/kg
   IF (90.5D0 .GT. tcd .OR. tcd .GT. 358.5D0) THEN
     call XABORT('THMHPT: exceed the valid temperature range(1).')
   ELSEIF (p .GT. 22.131475E6) THEN
     call XABORT('THMHPT: exceed the valid pressure range(1).')
   ELSE
     H0=0.0d0
     IF (tcd .LT. 257.5D0) THEN
        Z = (tcd - 90.5D0)/167.D0
        H0 = 365.515323D0 + (696.60209D0 - 9.3229851D0*Z)*Z &
         + 30.7092221D0*Z**3.535D0        
     ELSEIF (tcd .LT. 340.D0) THEN
        Z = (tcd - 257.5D0)/82.5D0
        H0 = 1083.50185D0 + (388.54695D0 + 40.6830346D0*Z)*Z &
         + 21.4058154D0*Z**4.79D0        
     ELSEIF (tcd .LT. 364.D0) THEN                      
        Z = (tcd - 340.D0)/24.D0
        H0 = 1534.137695D0 + (166.53D0 + 27.8460448D0*Z)*Z &
        + 15.9652752D0*Z**5.67D0        
     ELSE                       
        Z = (tcd - 364.D0)/5.D0
        H0 = 1744.47897D0 + (65.11525D0 + 15.3342016D0*Z)*Z &
        + 8.2093984D0*Z**5.74D0
     ENDIF
     A = (0.9769D0 - 0.695D-6*(tcd - 199.4D0)**2)*(H0+27.93D0) - H0
     B = 0.13D-1*tcd - 2.D0
     DT = 1.D06*(370.74D0 - tcd)
     HH = (B + DT/(22.131475D06 - psd))/(B+DT/(pd - psd))
     h = REAL(H0 + A*HH)*1.0E3
   ENDIF
   !
   ! compute specific heat capacity, in J/kg.K
   IF (90.5 .GT. tcd .OR. tcd .GT. 354.5) THEN
     call XABORT('THMHPT: exceed the valid temperature range(2).')
   ELSEIF (p .GT. 20.06E6) THEN
     call XABORT('THMHPT: exceed the valid pressure range(2).')
   ELSE
     CP0=0.0d0
     IF (tcd .LT. 216.D0) THEN
        Z = (tcd - 90.5D0)/125.5D0
        CP0 = 0.2398399D0 + (0.7260874D-2 - 0.8715753D-2*Z)*Z &
         - 0.1067229D-1*Z**2.43D0        
        CP0 = 1.D0/CP0
     ELSEIF (tcd .LT. 289.D0) THEN
        Z = (tcd - 216.D0)/73.D0
        CP0 = 4.3914986D0 + (0.4050098D0 + 0.2650923D0*Z)*Z &
         + 0.152426D0*Z**4.26D0        
     ELSEIF (tcd .LT. 334.D0) THEN                      
        Z = (tcd - 289.D0)/45.D0
        CP0 = 0.1917903D0 - (0.3592811D-1 + 0.1396726D-1*Z)*Z &
         - 0.5187553D-2*Z**4.07D0        
        CP0 = 1.D0/CP0
     ELSEIF (tcd .LT. 357.D0) THEN                      
        Z = (tcd - 334.D0)/23.D0
        CP0 = 0.1367074D0 - (0.4343216D-1 + 0.1360508D-1*Z)*Z &
         - 0.5536208D-2*Z**3.82D0        
        CP0 = 1.D0/CP0
     ELSE                       
        Z = (tcd - 357.D0)/9.5D0
        CP0 = 0.07413399D0 - (0.3791352D-1 + 0.6539416D-2*Z)*Z &
         - 0.2756629D-2*Z**2.58D0
        CP0 = 1.D0/CP0
     ENDIF
     A = -0.19878D0 + (1.521D0 - 0.393D0*CP0)*CP0
     B = -0.293594D0 + (0.45876D0 + 0.57448D-02*CP0)*CP0
     DT = 1.D06*(370.74D0 - tcd)
     CC = B + DT/(pd - psd)
     cp = REAL(CP0 + A/CC)*1.0E3
   ENDIF
   !
   ! use thermal conductivity and dynamic viscosity of light water
   td=dble(t)
   call free_pT(pd, td, rhod, hd, zkd, zmud, cpd)
   zk=real(zkd)
   zmu=real(zmud)
end subroutine THMHPT
!
subroutine THMHTX(t, x, rho, h, zk, zmu, cp)
! return the remaining thermohydraulics parameters as a function of the temperature (K)
! and quality
! Ref: Ji. Zhang, January 20, 98
   use, intrinsic :: iso_c_binding
   implicit real*8(a-h,o-z)
   real :: t, x, rho, h, zk, zmu, cp
   interface 
      subroutine free_Tx (td, xd, rhod, hd, zkd, zmud, cpd) bind(c, name='free_Tx')
      use, intrinsic :: iso_c_binding
      real(c_double) :: td, xd, rhod, hd, zkd, zmud, cpd
      end subroutine free_Tx
   end interface
   !
   tcd=dble(t-273.16)
   RO = 0.0D0
   H0 = 0.0D0
   CP0 = 0.0D0
   IF (x.EQ.0.0) THEN
      ! saturated liquid
      if (90.5D0 .GT. tcd .OR. tcd .GT. 367.D0) then
         call XABORT('THMHTX: the valid range of temperature is exceeded(1).')
      ENDIF
      !
      ! compute density in kg/m3
      IF (tcd .LT. 307.5D0) then
         Z = (tcd - 90.5D0)/217.D0
         RO = 1.07065471D3 - (0.1611572D3 + 0.12188254D3*Z)*Z &
               - 0.2106382D2*Z**6.83D0              
      ELSEIF (tcd .LT. 350.D0) then
         Z = (tcd - 307.5D0)/42.5D0
         RO = 0.7665511D3 - (0.1074816D3 + 0.2506107D2*Z)*Z &
               - 0.1042465D2*Z**4.64D0              
      ELSEIF (tcd .LT. 364.D0) then                                        
         Z = (tcd - 350.D0)/14.D0
         RO = 0.6235838D3 - (0.0678503D3 + 0.01551348D3*Z)*Z &
               - 0.7259508D1*Z**4.47D0              
      ELSE                                         
         Z = (tcd - 364.D0)/4.5D0
         RO = 0.5329606D3 - (0.4286569D2 + 0.1064107D2*Z)*Z &
               - 0.4761284D1*Z**5
      ENDIF
      !
      ! compute specific enthalpy, in J/kg
      IF (tcd .LT. 257.5D0) then
         Z = (tcd - 90.5D0)/167.D0
         H0 = 365.515323D0 + (696.60209D0 - 9.3229851D0*Z)*Z &
               + 30.7092221D0*Z**3.535D0              
      ELSEIF (tcd .LT. 340.D0) then
         Z = (tcd - 257.5D0)/82.5D0
         H0 = 1083.50185D0 + (388.54695D0 + 40.6830346D0*Z)*Z &
               + 21.4058154D0*Z**4.79D0              
      ELSEIF (tcd .LT. 364.D0) then                                        
         Z = (tcd - 340.D0)/24.D0
         H0 = 1534.137695D0 + (166.53D0 + 27.8460448D0*Z)*Z &
              + 15.9652752D0*Z**5.67D0              
      ELSE                                         
         Z = (tcd - 364.D0)/5.D0
         H0 = 1744.47897D0 + (65.11525D0 + 15.3342016D0*Z)*Z &
              + 8.2093984D0*Z**5.74D0
      ENDIF
      !
      ! compute specific heat capacity at constant pressure, in J/kg/K
      IF (tcd .LT. 216.D0) then
         Z = (tcd - 90.5D0)/125.5D0
         CP0 = 0.2398399D0 + (0.7260874D-2 - 0.8715753D-2*Z)*Z &
               - 0.1067229D-1*Z**2.43D0              
         CP0 = 1.D0/CP0
      ELSEIF (tcd .LT. 289.D0) then
         Z = (tcd - 216.D0)/73.D0
         CP0 = 4.3914986D0 + (0.4050098D0 + 0.2650923D0*Z)*Z &
               + 0.152426D0*Z**4.26D0              
      ELSEIF (tcd .LT. 334.D0) then                                        
         Z = (tcd - 289.D0)/45.D0
         CP0 = 0.1917903D0 - (0.3592811D-1 + 0.1396726D-1*Z)*Z &
               - 0.5187553D-2*Z**4.07D0              
         CP0 = 1.D0/CP0
      ELSEIF (tcd .LT. 357.D0) then                                        
         Z = (tcd - 334.D0)/23.D0
         CP0 = 0.1367074D0 - (0.4343216D-1 + 0.1360508D-1*Z)*Z &
               - 0.5536208D-2*Z**3.82D0              
         CP0 = 1.D0/CP0
      ELSE                                         
         Z = (tcd - 357.D0)/9.5D0
         CP0 = 0.07413399D0 - (0.3791352D-1 + 0.6539416D-2*Z)*Z &
               - 0.2756629D-2*Z**2.58D0
         CP0 = 1.D0/CP0
      ENDIF
  ELSEIF (x.EQ.1.0) THEN
      ! saturated steam
      IF (90.5D0 .GT. tcd .OR. tcd .GT. 367.D0) then
         call XABORT('THMHTX: the valid range of temperature is exceeded(2).')
      ENDIF
      !
      ! compute density in kg/m3
      IF (tcd .GT. 350.D0) then
         Z = (tcd - 350.D0)/16.D0
         RO = 7.5017113D0 - (2.7448665D0 - 0.4712976D-1*Z)*Z &
               - 0.1461081D0*Z**5.45D0              
         RO = 1.D3/RO
      ELSEIF (tcd .GT. 288.D0) then
         Z = (tcd - 288.D0)/62.D0
         RO = 23.2557129D0 - (24.2703096D0 - 12.7121893D0*Z)*Z &
               - 4.1958813D0*Z**2.81D0              
         RO = 1.D3/RO
      ELSEIF (tcd .GT. 221.D0) then                                        
         Z = (tcd - 221.D0)/67.D0
         RO = 0.01319205D3 + (0.1687167D2 + 0.9533432D1*Z)*Z &
               + 0.3402844D1*Z**3.69D0              
      ELSEIF (tcd .GT. 147.5D0) then                                        
         Z = (tcd - 147.5D0)/73.5D0
         RO = 0.2595185D1 + (0.4966204D1 + 0.3980025D1*Z)*Z &
               + 0.1650632D1*Z**3.382D0              
      ELSE                                         
         Z = (tcd - 90.5D0)/57.D0
         RO = 0.451542D0 + (0.9337486D0 + 0.8179223D0*Z)*Z &
               + 0.3919721D0*Z**3.27D0
      ENDIF
      !
      ! compute specific enthalpy, in J/kg
      IF (tcd .LT. 259.D0) then
         Z = (tcd - 90.5D0)/168.5D0
         H0 = 2465.02D0 + (257.038325D0 - 69.298269D0*Z)*Z &
              - 59.547036D0*Z**3.39D0              
      ELSEIF (tcd .LT. 333.D0) then
         Z = (tcd - 259.D0)/74.D0
         H0 = 2593.21302D0 - (36.63666D0 + 70.6712028D0*Z)*Z &
              - 26.1578672D0*Z**4.41D0              
      ELSEIF (tcd .LT. 359.D0) then                                        
         Z = (tcd - 333.D0)/26.D0
         H0 = 2459.74729D0 - (103.06374D0 + 40.5743371D0*Z)*Z &
              - 17.2902029D0*Z**4.76D0              
      ELSE                                         
         Z = (tcd - 359.D0)/8.5D0
         H0 = 2298.81901D0 - (87.129505D0 + 27.7163491D0*Z)*Z &
              - 14.0347359D0*Z**5.2D0
      ENDIF
      !
      ! compute specific heat capacity at constant pressure, in J/kg/K
      IF (tcd .LT. 208.D0) then
         Z = (tcd - 90.5D0)/117.5D0
         CP0 = 1.8689755D0 + (0.3394869D0 + 0.2728998D0*Z)*Z &
               + 0.2686182D0*Z**3.507D0              
      ELSEIF (tcd .LT. 270.D0) then
         Z = (tcd - 208.D0)/62.D0
         CP0 = 2.7499805D0 + (0.9642085D0 + 0.4429098D0*Z)*Z &
               + 0.141205D0*Z**3.945D0              
      ELSEIF (tcd .LT. 339.D0) then                                        
         Z = (tcd - 270.D0)/69.D0
         CP0 = 0.2326499D0 - (0.1449925D0 - 0.19935345D-2*Z)*Z &
               - 0.441272D-2*Z**3.96D0              
         CP0 = 1.D0/CP0
      ELSEIF (tcd .LT. 363.D0) then                                        
         Z = (tcd - 339.D0)/24.D0
         CP0 = 0.08523823D0 - (0.5512341D-1 + 0.3706342D-2*Z)*Z &
               - 0.2012056D-2*Z**4.26D0              
      ELSE                                         
         Z = (tcd - 363.D0)/4.D0
         CP0 = 0.02439642D0 - (0.1185124D-1 + 0.4619397D-3*Z)*Z &
               - 0.1003777D-3*Z**2.4D0
         CP0 = 1.D0/CP0
      ENDIF
   ELSE
      CALL XABORT('THMHTX: quality = 0.0 or 1.0 expected.')
   ENDIF
   rho = REAL(RO)
   h = REAL(H0)*1.0E3
   cp = REAL(CP0)*1.0E3
   !
   ! use thermal conductivity and dynamic viscosity of light water
   td=dble(t)
   xd=dble(x)
   call free_Tx(td, xd, rhod, hd, zkd, zmud, cpd)
   zk=real(zkd)
   zmu=real(zmud)
end subroutine THMHTX
