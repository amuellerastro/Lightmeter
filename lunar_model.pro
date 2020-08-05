; Model based on Kieffer&Stone - 2005 - The Spectral Irradiance of the Moon, AJ, 129,2887

; function lightmeter_lunar_model_fit_func, x, p
; 
;   common data0, counts_moon_cut, moon_irr_tot_cut, temp_moon_cut
; 
;   fit = p[2]*(p[1]*(p[0]*exp((counts_moon_cut*(1.d0+3.d-3*temp_moon_cut))/p[0])-1.d0)+(counts_moon_cut*(1.d0+3.d-3*temp_moon_cut)))
; 
;   return, total((fit-moon_irr_tot_cut)^2./moon_irr_tot_cut)
;   ;return, total((counts_moon-fit)^2.d0 * abs(dummy_err))
; 
; end



pro lunar_model, pathbase, pathlunfit, longdeg, latdeg, altitude, jd_day, et, jd, hour, longmed, dist_EarthSun, t18set, t18rise, $
		counts, temp_kelvin, temp_celsius, binsize, hrbin, rvarcounts, h_moon, h_moon_orig, k_moon, g, hr_moon, counts_moon, moon_irr_tot, temp_moon, check

;   common data0, counts_moon_cut, moon_irr_tot_cut, temp_moon_cut


;path contains 'hard coded data'
  path = pathbase+'fix/'

  app_long = longmed
  dist_ES = dist_EarthSun

;compute some basic ephemerides of the moon

;compute LST, UT assumed!
  ct2lst, lst, longdeg, 1., double(jd)

;compute RA, DEC [deg] of moon
  print, '   computation of lunar ephemerides...'
  moonpos_mod, et, ra_moon, dec_moon, dist_EarthMoon, geolong, geolat, F, M, Mprime, D, E

;computation of h_moon
  eq2hor_mod, ra_moon, dec_moon, et, h_moon, az_moon, ha_moon, LON=longdeg, LAT=latdeg, precess_=1, nutate_=1, refract_=1, aberration_=1, altitude=altitude, temperature=temp_kelvin

  h_moon_orig = h_moon

;compute zenith distance zd of the moon
  zd = 90.d0 - h_moon

;sort for AIRMASS
;do not consider an airmass lower than XX because atmospheric models becoming bad for low altitudes
;if moon is not visible in during night (h_sun<-18deg) and don't reach AM<2.9 skip lunar model computation and return to main level
  idx = where(abs(zd) le 70.)

  if (idx[0] eq -1) then begin

    check = 'no'
    hr_moon = hour
    h_moon = h_moon_orig
    counts_moon = counts
    temp_moon = temp_celsius
    moon_irr_tot = 0.

    et_moon = et

    ;computation lunar phase angle and illumination fraction
      selen_coord_sun, et_moon, dist_EarthMoon, geolong, geolat, app_long, dist_ES, F, M, Mprime, D, E, $
		      lambda_0_tmp, lambda_tmp, beta_tmp, l0_tmp, b0_tmp, R_tmp

      lambda_0 = lambda_0_tmp
      lambda = lambda_tmp
      beta = beta_tmp
      l0 = l0_tmp
      b0 = b0_tmp
      delta = dist_EarthMoon
      R = R_tmp

      ;calculate absolute phase angle "g"
      psi = acos(cos(beta*!dtor)*cos((lambda*!dtor)-(lambda_0*!dtor)))
      g_dum1 = (R*sin(psi))
      g_dum2 = (delta-(R*cos(psi)))
      g = (atan(g_dum1, g_dum2))


      ;calculate illuminated fraction of the Moon
      k_moon = (1.d0+cos(g))/2.d0

    return

  endif

  zd = zd[idx]	;[deg]
  jd_moon = jd[idx]
  et_moon = et[idx]
  counts_moon = counts[idx]
  temp_moon = temp_celsius[idx]
  hr_moon = hour[idx]
  ;ha_moon = ha_moon[idx]
  h_moon = h_moon[idx]
  ra_moon = ra_moon[idx]
  dec_moon = dec_moon[idx]
  dist_EarthMoon = dist_EarthMoon[idx]
  geolong = geolong[idx]
  geolat = geolat[idx]
  F = F[idx]
  M = M[idx]
  Mprime = Mprime[idx]
  D = D[idx]
  E = E[idx]
  app_long = app_long[idx]
  dist_ES = dist_ES[idx]


  ;do check of visibility during night
  timesort_moon = hr_moon
;   tdummy = dblarr(n_elements(timesort_moon))
;   for j=0L,n_elements(timesort_moon)-1 do begin
;     if (timesort_moon[j] ge 12.) then begin
;       tdummy[j] = timesort_moon[j]-12.
;     endif
;     if (timesort_moon[j] lt 12.) then begin
;       tdummy[j] = timesort_moon[j]+12.
;     endif
;   endfor
;   timesort_moon = tdummy


;sort for VISIBILITY during NIGHT
  idx_check = where(timesort_moon ge t18set and timesort_moon le t18rise)

  if (idx_check[0] ne -1) then begin

    check = 'yes'	;i.e. moon is visible during dark time under given conditions, i.e. AM

    zd = zd[idx_check]	;[deg]
    jd_moon = jd_moon[idx_check]
    et_moon = et_moon[idx_check]
    counts_moon = counts_moon[idx_check]
    temp_moon = temp_moon[idx_check]
    hr_moon = hr_moon[idx_check]
    ;ha_moon = ha_moon[idx_check]
    h_moon = h_moon[idx_check]
    ra_moon = ra_moon[idx_check]
    dec_moon = dec_moon[idx_check]
    dist_EarthMoon = dist_EarthMoon[idx_check]
    geolong = geolong[idx_check]
    geolat = geolat[idx_check]
    F = F[idx_check]
    M = M[idx_check]
    Mprime = Mprime[idx_check]
    D = D[idx_check]
    E = E[idx_check]
    app_long = app_long[idx_check]
    dist_ES = dist_ES[idx_check]


  ;cut out areas having clouds
    idx_clean = where(rvarcounts lt 0.02d0)
    if (idx_clean[0] ne -1) then begin

      hrbin_clean = hrbin[idx_clean]
      rvarcounts_clean = rvarcounts[idx_clean]

    endif else begin

      hrbin_clean = hrbin
      rvarcounts_clean = rvarcounts

    endelse


    ;create time array in order to match hr_moon with the non cloudy parts of the night
    tleft = hrbin_clean - binsize[0]/(2.d0*3600.d0)
    tright = hrbin_clean + binsize[0]/(2.d0*3600.d0)

    hr_tmp = dblarr(n_elements(hr_moon)) & jd_tmp = hr_tmp & et_tmp = hr_tmp & zd_tmp = hr_tmp & counts_moon_tmp = hr_tmp
    h_moon_tmp = hr_tmp & ra_moon_tmp = hr_tmp & dec_moon_tmp = hr_tmp & temp_moon_tmp = hr_tmp & app_long_tmp = hr_tmp & dist_ES_tmp = hr_tmp
    dist_EarthMoon_tmp = hr_tmp & geolong_tmp = hr_tmp & geolat_tmp = hr_tmp & F_tmp = hr_tmp 
    M_tmp = hr_tmp & Mprime_tmp = hr_tmp & D_tmp = hr_tmp & E_tmp = hr_tmp	;ha_moon_tmp = hr_tmp
    for i=0L,n_elements(tleft)-1 do begin

      idx = where(hr_moon ge tleft[i] and hr_moon le tright[i])
;      idx = where(tleft[i] le hr_moon and tright[i] ge hr_moon)

      if (idx[0] ne -1) then begin

	idx0 = where(hr_tmp eq 0.)
	hr_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = hr_moon[idx]
	jd_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = jd_moon[idx]
	et_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = et_moon[idx]
	zd_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = zd[idx]
	counts_moon_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = counts_moon[idx]
	temp_moon_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = temp_moon[idx]
	;ha_moon_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = ha_moon[idx]
	h_moon_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = h_moon[idx]
	ra_moon_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = ra_moon[idx]
	dec_moon_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = dec_moon[idx]
	dist_EarthMoon_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = dist_EarthMoon[idx]
	geolong_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = geolong[idx]
	geolat_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = geolat[idx]
	F_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = F[idx]
	M_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = M[idx]
	Mprime_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = Mprime[idx]
	D_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = D[idx]
	E_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = E[idx]
	app_long_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = app_long[idx]
	dist_ES_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = dist_ES[idx]

      endif

    endfor

    idxval = where(hr_tmp ne 0.)
    if (idxval[0] ne -1) then begin

      hr_tmp = hr_tmp[idxval]
      jd_tmp = jd_tmp[idxval]
      et_tmp = et_tmp[idxval]
      zd_tmp = zd_tmp[idxval]
      counts_moon_tmp = counts_moon_tmp[idxval]
      temp_moon_tmp = temp_moon_tmp[idxval]
      ;ha_moon_tmp = ha_moon_tmp[idxval]
      h_moon_tmp = h_moon_tmp[idxval]
      ra_moon_tmp = ra_moon_tmp[idxval]
      dec_moon_tmp = dec_moon_tmp[idxval]
      dist_EarthMoon_tmp = dist_EarthMoon_tmp[idxval]
      geolong_tmp = geolong_tmp[idxval]
      geolat_tmp = geolat_tmp[idxval]
      F_tmp = F_tmp[idxval]
      M_tmp = M_tmp[idxval]
      Mprime_tmp = Mprime_tmp[idxval]
      D_tmp = D_tmp[idxval]
      E_tmp = E_tmp[idxval]
      app_long_tmp = app_long_tmp[idxval]
      dist_ES_tmp = dist_ES_tmp[idxval]



      jd_moon = jd_tmp
      et_moon = et_tmp
      hr_moon = hr_tmp
      zd = zd_tmp
      counts_moon = counts_moon_tmp
      temp_moon = temp_moon_tmp
      ;ha_moon = ha_moon_tmp
      h_moon = h_moon_tmp
      ra_moon = ra_moon_tmp
      dec_moon = dec_moon_tmp
      dist_EarthMoon = dist_EarthMoon_tmp
      geolong = geolong_tmp
      geolat = geolat_tmp
      F = F_tmp
      M = M_tmp
      Mprime = Mprime_tmp
      D = D_tmp
      E = E_tmp
      app_long = app_long_tmp
      dist_ES = dist_ES_tmp


    endif else begin

      check = 'no'
      hr_moon = hour
      h_moon = h_moon_orig
      counts_moon = counts
      temp_moon = temp_celsius
      moon_irr_tot = 0.

      ;computation lunar phase angle and illumination fraction
	selen_coord_sun, et_moon, dist_EarthMoon, geolong, geolat, app_long, dist_ES, F, M, Mprime, D, E, $
			lambda_0_tmp, lambda_tmp, beta_tmp, l0_tmp, b0_tmp, R_tmp

	lambda_0 = lambda_0_tmp
	lambda = lambda_tmp
	beta = beta_tmp
	l0 = l0_tmp
	b0 = b0_tmp
	delta = dist_EarthMoon
	R = R_tmp

	;calculate absolute phase angle "g"
	psi = acos(cos(beta*!dtor)*cos((lambda*!dtor)-(lambda_0*!dtor)))
	g_dum1 = (R*sin(psi))
	g_dum2 = (delta-(R*cos(psi)))
	g = (atan(g_dum1, g_dum2))


	;calculate illuminated fraction of the Moon
	k_moon = (1.d0+cos(g))/2.d0

      return

    endelse

    ;computation of selenographic coordinates, lunar phase angle, and illumination

      print, '   computation of selenographic coordinates...'

      ;calculate selenographic coordinates "l0" and "b0" of the Sun; Reference: Meeus - Astronomical Algorithms, chapter 53
      selen_coord_sun, et_moon, dist_EarthMoon, geolong, geolat, app_long, dist_ES, F, M, Mprime, D, E, $
		      lambda_0_tmp, lambda_tmp, beta_tmp, l0_tmp, b0_tmp, R_tmp

      lambda_0 = lambda_0_tmp
      lambda = lambda_tmp
      beta = beta_tmp
      l0 = l0_tmp
      b0 = b0_tmp
      delta = dist_EarthMoon
      R = R_tmp

      ;calculate selenographic coordinates "l" and "b" of the Earth
      selen_coord_earth, et_moon, dist_EarthMoon, geolong, geolat, F, M, Mprime, D, E, l_tmp, b_tmp

      l = l_tmp
      b = b_tmp

      ;calculate absolute phase angle "g"
      psi = acos(cos(beta*!dtor)*cos((lambda*!dtor)-(lambda_0*!dtor)))
      g_dum1 = (R*sin(psi))
      g_dum2 = (delta-(R*cos(psi)))
      g = (atan(g_dum1, g_dum2))


      ;calculate illuminated fraction of the Moon
      k_moon = (1.d0+cos(g))/2.d0


;apply ROLO lunar model and compute irradiance

      lnA = dblarr(n_elements(jd_moon),18) & A = lnA
      moon_irr_0 = lnA & moon_irr = moon_irr_0

    ; read in of the Model Coefficients, v311g
      readcol, path+'Kieffer2005_table4_modelcoefficients.dat',wavelength_Kieffer,a0,a1,a2,a3,b1,b2,b3,d1,d2,d3,format='d,d,d,d,d,d,d,d,d,d,d', /silent


    ;constants from Kieffer et al. 2005
      c1 = 0.00034115d0
      c2 = -0.0013425d0
      c3 = 0.00095906d0
      c4 = 0.00066229d0
      p1 = 4.06054d0
      p2 = 12.8802d0
      p3 = -30.5858d0
      p4 = 16.7498d0


    ; disk-equivalent reflectance A (depends on wavelength)
    ; A is computed for 350.0 - 774.8 nm


      for i=0L,n_elements(wavelength_Kieffer)-1 do begin
	;Kieffer changes dimension of the constants faster than I cange my socks...2nd expression should be OK
; 	lnA[*,i] = a0[i]*g^0. + a1[i]*g^1. + a2[i]*g^2. + a3[i]*g^3. + b1[i]*(l0*!dtor)^1. $
; 		+ b2[i]*(l0*!dtor)^3. + b3[i]*(l0*!dtor)^5. + c1*(b*!dtor) + c2*(l*!dtor) $
; 		+ c3*(l0*!dtor)*(b*!dtor) + c4*(l0*!dtor)*(l*!dtor) + d1[i]*exp(-g/p1) + d2[i]*exp(-g/p2) $
; 		+ d3[i]*cos((g-p3)/p4)

	lnA[*,i] = a0[i]*g^0. + a1[i]*g^1. + a2[i]*g^2. + a3[i]*g^3. $
		+ b1[i]*(l0*!dtor)^1. + b2[i]*(l0*!dtor)^3. + b3[i]*(l0*!dtor)^5. $
		+ c1*b + c2*l + c3*l0*(b*!dtor) + c4*l0*(l*!dtor) $
		+ d1[i]*exp(-(g/!dtor)/p1) + d2[i]*exp(-(g/!dtor)/p2) + d3[i]*cos((((g/!dtor)-p3)/p4)*!dtor)

      endfor

      A = exp(lnA[*,*])


    ;read irradiance data from the sun ASTM E490 00a
      readcol, path+'ASTM_E490_00a.dat', wavelength_sun_irr, sun_irr, format='d,d', /silent	;[nm],W/m^2/nm]

    ;compute irradiance received by the moon
      r_moon = 1737.40d0	;astronomical almanac 2010, K7
      ;omega_moon = !DPI*r_moon*r_moon/(delta*delta)	; calculate solid angle of the Moon
      omega_moon = 2.d0*!DPI*(1.d0-cos(atan(r_moon/delta)))	; calculate solid angle of the Moon


      for i=0L,n_elements(wavelength_Kieffer)-1 do begin

	moon_irr_0[*,i] = A[*,i]*sun_irr[i]*omega_moon/(!DPI)

      endfor


    ;compute extinction "A_prime" after Hayes and Latham 1975
      extinction, altitude, wavelength_Kieffer, A_prime

    ;cimpute airmass
      AM = airmass_rozenberg(zd)

      readcol, path+'Lightmeter_SpecResp.dat', lam, val, format='d,d', /silent ;spectral response of solar cell

      ;moon_irr_lux = moon_irr
      for i=0L,n_elements(wavelength_Kieffer)-1 do begin
	moon_irr[*,i] = moon_irr_0[*,i]*val[i]/(100.0d0^0.2d0^(A_prime[i]*AM))	;same as moon_irr_0[i]*val[i]*10.^(-(A_prime[i]*AM)/2.5d0)
      endfor

      ;compute total irradiance from moon
	moon_irr_tot = dblarr(n_elements(jd_moon))


  ;the routines tsum and int_tabulated produce different results. values of int_tabulated are smaller
      for i=0L,n_elements(jd_moon)-1 do begin
	moon_irr_tot[i] = tsum(wavelength_Kieffer, moon_irr[i,*])
      endfor

      ;constrain for geometrical projection
      moon_irr_tot = moon_irr_tot*cos(zd*!dtor)


  endif else begin

    check = 'no'
    hr_moon = hour
    h_moon = h_moon_orig
    counts_moon = counts
    temp_moon = temp_celsius
    moon_irr_tot = 0.

    print, '   Moon is not visible during this night or does not reach AM<2.9.'

    ;computation lunar phase angle and illumination fraction
      selen_coord_sun, et_moon, dist_EarthMoon, geolong, geolat, app_long, dist_ES, F, M, Mprime, D, E, $
		      lambda_0_tmp, lambda_tmp, beta_tmp, l0_tmp, b0_tmp, R_tmp

      lambda_0 = lambda_0_tmp
      lambda = lambda_tmp
      beta = beta_tmp
      l0 = l0_tmp
      b0 = b0_tmp
      delta = dist_EarthMoon
      R = R_tmp

      ;calculate absolute phase angle "g"
      psi = acos(cos(beta*!dtor)*cos((lambda*!dtor)-(lambda_0*!dtor)))
      g_dum1 = (R*sin(psi))
      g_dum2 = (delta-(R*cos(psi)))
      g = (atan(g_dum1, g_dum2))


      ;calculate illuminated fraction of the Moon
      k_moon = (1.d0+cos(g))/2.d0


    return

  endelse


end