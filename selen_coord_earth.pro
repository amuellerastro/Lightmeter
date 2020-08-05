pro selen_coord_earth, jd, dist_em, geolong, geolat, F, M, Mprime, D, E, l_dum, b_dum

;=====================================================================
;calculate the heliocentric coordinates "l" and "b" of the moon (p376)
;=====================================================================


  delta = dist_em	;Earth-Moon distance
  lambda = geolong
  beta = geolat


; calculate l_prime and b_prime

  nutate_mod, jd, nut_long, nut_obliq, Omega
  Omega = Omega/!dtor
  deltapsi = nut_long/3600.0d0

  W = lambda - deltapsi - Omega
  cirrange, W

  I_lun = ten(1.,32.,32.7)	;inclination of the mean lunar equator to the ecliptic

  dum1 = ((sin(W*!dtor)*cos(beta*!dtor)*cos(I_lun*!dtor)) - (sin(beta*!dtor)*sin(I_lun*!dtor)))
  dum2 = (cos(W*!dtor)*cos(beta*!dtor))
  A = atan(dum1, dum2)

  l_prime = (A/!dtor) - (F/!dtor)
  b_prime = asin(-(sin(W*!dtor)*cos(beta*!dtor)*sin(I_lun*!dtor)) - (sin(beta*!dtor)*cos(I_lun*!dtor)))
  b_prime = b_prime/!dtor


; calculate l_prime2 and b_prime2

  M_prime = Mprime
  rho = (-0.02752d*cos(M_prime))-(0.02245d*sin(F))+(0.00684d*cos(M_prime-(2.d*F)))-(0.00293d*cos(2.d*F))-(0.00085d*cos((2.d*F)-(2.d*D))) $
        -(0.00054d*cos(M_prime-(2.d*D)))-(0.00020d*sin(M_prime+F))-(0.00020d*cos(M_prime+(2.d*F)))-(0.00020d*cos(M_prime-F)) $
        +(0.00014d*cos(M_prime+(2.d*F)-(2.d*D)))

  sigma = (-0.02816d*sin(M_prime))+(0.02244d*cos(F))-(0.00682d*sin(M_prime-(2.d*F)))-(0.00279d*sin(2.d*F))-(0.00083d*sin((2.d*F)-(2.d*D))) $
          +(0.00069d*sin(M_prime-(2.d*D)))+(0.00040d*cos(M_prime+F))-(0.00025d*sin(2.d*M_prime))-(0.00023d*sin(M_prime+(2.d*F))) $
          +(0.0002d*cos(M_prime-F))+(0.00019d*sin(M_prime-F))+(0.00013d*sin(M_prime+(2.d*F)-(2.d*D)))-(0.00010d*cos(M_prime-(3.d*F)))

  K1 = 119.75d0+(131.849d0*(jd-2451545.0d0)/36525.0d0)
  K2 = 72.56d0+(20.186d0*(jd-2451545.0d0)/36525.0d0)
  tau = (0.0252d*E*sin(M))+(0.00473d*sin((2.d*M_prime)-(2.d*F)))-(0.00467d*sin(M_prime))+(0.00396d*sin(K1*!dtor)) $
        +(0.00276d*sin((2.d*M_prime)-(2.d*D)))+(0.00196d*sin(Omega*!dtor))-(0.00183d*cos(M_prime-F))+(0.00115d*sin(M_prime-(2.d*D))) $
        -(0.00096d*sin(M_prime-D))+(0.00046d*sin((2.d*F)-(2.d*D)))-(0.00039d*sin(M_prime-F))-(0.00032d*sin(M_prime-M-D)) $
        +(0.00027d*sin((2.d*M_prime)-M-(2.d*D)))+(0.00023d*sin(K2*!dtor))-(0.00014d*sin(2.d*D))+(0.00014d*cos((2.d*M_prime)-(2.d*F))) $
        -(0.00012d*sin(M_prime-(2.d*F)))-(0.00012d*sin(2.d*M_prime))+(0.00011d*sin((2.d*M_prime)-(2.d*M)-(2.d*D)))

  l_prime2 = -(tau*!dtor) + (((rho*!dtor)*cos(A*!dtor))+((sigma*!dtor)*sin(A*!dtor)))*tan(b_prime*!dtor)
  b_prime2 = ((sigma*!dtor)*cos(A*!dtor))-((rho*!dtor)*sin(A*!dtor))
  l_prime2 = l_prime2/!dtor
  b_prime2 = b_prime2/!dtor


; calculate l and b

  l_dum = l_prime + l_prime2

  idx = where(l_dum lt -180.0d0)
  if (idx[0] ne -1) then l_dum[idx] = 360.0d0 + l_dum[idx]
  ;if (l_dum lt -180.0d0) then l_dum = 360.0d0 + l_dum

  b_dum = b_prime + b_prime2

  ;print, rho, sigma, tau
  ;print, l, b

return
end