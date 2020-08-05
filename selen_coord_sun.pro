pro selen_coord_sun, jd, dis, geolong, geolat, app_long, dist,$
                     F, M, Mprime, D, E, lambda_0_dum, lambda_dum, beta_dum, $
                     l0_dum, b0_dum, R_dum

;==================================================
;selenographic coordinates "l0" and "b0" of the Sun
;==================================================


;calculate the heliocentric coordinates of the moon (p376)
  delta_dum = dis	;Earth-Moon distance
  lambda_dum = geolong
  beta_dum = geolat

  R_dum = dist*149597871.464d0	;Astronomical Almanac, 2010, K6
  lambda_0_dum = app_long	;apparent geocentric longitude of the sun

  lambda_h = lambda_0_dum + 180.0d + (delta_dum/R_dum)*57.296*cos(beta_dum*!dtor)*sin((lambda_0_dum - lambda_dum)*!dtor) 
  lambda_h = lambda_h mod 360.0d
  beta_h = (delta_dum/R_dum)*beta_dum


; calculate l0_prime and b0_prime

  nutate_mod, jd, nut_long, nut_obliq, Omega
  Omega = Omega/!dtor
  deltapsi = nut_long/3600.0d0

  W = lambda_h - deltapsi - Omega
  cirrange, W

  I_lun = ten(1.d0,32.d0,32.7d0)	;inclination of the mean lunar equator to the ecliptic

  dum1 = ((sin(W*!dtor)*cos(beta_h*!dtor)*cos(I_lun*!dtor)) - (sin(beta_h*!dtor)*sin(I_lun*!dtor)))
  dum2 = (cos(W*!dtor)*cos(beta_h*!dtor))
  A = atan(dum1, dum2)

  l0_prime = (A/!dtor) - (F/!dtor)
  b0_prime = asin(-(sin(W*!dtor)*cos(beta_h*!dtor)*sin(I_lun*!dtor)) - (sin(beta_h*!dtor)*cos(I_lun*!dtor)))
  b0_prime = b0_prime/!dtor


; calculate l0_prime2 and b0_prime2

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

  l0_prime2 = -(tau*!dtor) + (((rho*!dtor)*cos(A*!dtor))+((sigma*!dtor)*sin(A*!dtor)))*tan(b0_prime*!dtor)
  b0_prime2 = ((sigma*!dtor)*cos(A*!dtor))-((rho*!dtor)*sin(A*!dtor))
  l0_prime2 = l0_prime2/!dtor
  b0_prime2 = b0_prime2/!dtor

  l0_dum = l0_prime + l0_prime2

  idx = where(l0_dum lt -180.d0)
  if (idx[0] ne -1) then l0_dum[idx] = 360.0d0 + l0_dum[idx]
;  if (l0_dum lt -180.0d0) then l0_dum = 360.0d0 + l0_dum

  b0_dum = b0_prime + b0_prime2

; Colongitude of the Sun
;   if (l0 lt 90.) then begin
;      co=90.-l0
;   endif else begin
;     co = 450.-l0
;   endelse

;  print, l0, b0;, co

return
end