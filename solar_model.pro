pro solar_model, pathbase, altitude, t0set, t0rise, time_orig, counts_orig, temp_celsius_orig, $
		 timesort_synth=timesort_synth, h_sun_synth, sun_binsize, sun_binsize10, sun_hrbin, sun_rvarcounts, $
		 time_sun, counts_sun, temp_sun, h_sun, sun_irr_tot

;path contains 'hard coded data'
  path = pathbase+'fix/'


;rename
  time = time_orig
  counts = counts_orig
  temp_celsius = temp_celsius_orig


;select day time data

  ;lower limit for sun altitude for solar model
  hlow_sun = 20.d0

  ;cut out night and twilight
  idxday = where(time le t0set or time ge t0rise)
  time_sun = time[idxday]
  counts_sun = counts[idxday]
  temp_sun = temp_celsius[idxday]

if keyword_set(timesort_synth) then begin

;extract corresponding altitude of sun
  ;synthetic data have to be synchronized with measurd data
  ;necessary because of limited machine precision
  tmp1 = round(time_sun*1.d12, /L64)
  tmp2 = round(timesort_synth*1.d12, /L64)

  match2, tmp1, tmp2, idx1, idx2

  h_sun = h_sun_synth[idx1]

endif else begin

  h_sun = h_sun_synth

endelse


  idxh = where(h_sun ge hlow_sun)
  if (idxh[0] ne -1) then begin

    time_sun = time_sun[idxh]
    counts_sun = counts_sun[idxh]
    temp_sun = temp_sun[idxh]
    h_sun = h_sun[idxh]

  endif


; zenith distance of sun
  zd_sun = 90.d0-h_sun


; read irradiance data from the sun ASTM E490 00a
  readcol, path+'ASTM_E490_00a.dat', wavelength_sun_irr, sun_model, format='d,d', /silent	;[nm],W/m^2/nm]

; read in of the Model Coefficients, v311g, but only need wavelength
  readcol, path+'Kieffer2005_table4_modelcoefficients.dat',wavelength_Kieffer,a0,a1,a2,a3,b1,b2,b3,d1,d2,d3,format='d,d,d,d,d,d,d,d,d,d,d', /silent

; read in Lightmeter spectral response
  readcol, path+'Lightmeter_SpecResp.dat', lam, val, format='d,d', /silent ;spectral response of solar cell

; compute extinction "A_prime" after Hayes and Latham 1975
  extinction, altitude, wavelength_Kieffer, A_prime

; compute Airmass
  AM_sun = airmass_rozenberg(zd_sun)

; compute total irradiance of sun
  sun_irr = dblarr(n_elements(time_sun),18)
  sun_irr_tot = dblarr(n_elements(time_sun))
  for i=0L,n_elements(wavelength_Kieffer)-1 do begin
    sun_irr[*,i] = sun_model[i]*val[i]/(100.0d0^0.2d0^(A_prime[i]*AM_sun))	;same as sun_model[i]*val[i]*10.^(-(A_prime[i]*AM_sun)/2.5d0)
  endfor

  for i=0L,n_elements(time_sun)-1 do begin
    sun_irr_tot[i] = tsum(wavelength_Kieffer, sun_irr[i,*])
  endfor

  sun_irr_tot = sun_irr_tot*cos(zd_sun*!dtor)
  ;680 W/m^2 maximum (pure model w/o extinction and spectral sensitivity)


;SUN: consider only 'clear' parts
;select data under good conditions using LOSSAM computation

    idx_clean = where(sun_rvarcounts lt 0.01d0)
    if (idx_clean[0] ne -1) then begin

      hrbin_clean = sun_hrbin[idx_clean]
      rvarcounts_clean = sun_rvarcounts[idx_clean]

    endif else begin

      hrbin_clean = sun_hrbin
      rvarcounts_clean = sun_rvarcounts

    endelse

    tleft = hrbin_clean - sun_binsize[0]/(2.d0*3600.d0)
    tright = hrbin_clean + sun_binsize[0]/(2.d0*3600.d0)

    time_sun_tmp = dblarr(n_elements(time_sun)) & temp_sun_tmp = time_sun_tmp
    counts_sun_tmp = time_sun_tmp & h_sun_tmp = time_sun_tmp & sun_irr_tot_tmp = time_sun_tmp

    for i=0L,n_elements(tleft)-1 do begin

      idx = where(time_sun ge tleft[i] and time_sun le tright[i])

      if (idx[0] ne -1) then begin

	idx0 = where(time_sun_tmp eq 0.)
	time_sun_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = time_sun[idx]
	temp_sun_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = temp_sun[idx]
	counts_sun_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = counts_sun[idx]
	h_sun_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = h_sun[idx]
	sun_irr_tot_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = sun_irr_tot[idx]

      endif

    endfor

    idxval = where(time_sun_tmp ne 0.)
    if (idxval[0] ne -1) then begin

      time_sun_tmp = time_sun_tmp[idxval]
      temp_sun_tmp = temp_sun_tmp[idxval]
      counts_sun_tmp = counts_sun_tmp[idxval]
      h_sun_tmp = h_sun_tmp[idxval]
      sun_irr_tot_tmp = sun_irr_tot_tmp[idxval]

      time_sun = time_sun_tmp
      temp_sun = temp_sun_tmp
      counts_sun = counts_sun_tmp
      h_sun = h_sun_tmp
      sun_irr_tot = sun_irr_tot_tmp

    endif


end
;======================================================================================
;NOT USED SIMPLE SOLAR MODEL
; 
; 
; ; airmass
;   AM = 1.d0/(cos(zd_sun*!dtor) + (0.025d0*exp(-11.0d0*cos(zd_sun*!dtor))))
; 
; ; estimate Pressure based on altitude, using U.S. Standard Atmosphere formula.
;   pressure = (1010.d0*(1.d0-6.5d0/288000.d0*altitude)^5.255d0)*0.1d0
; 
; ; Rayleigh and permanent gas scattering and absorption transmission coefficient
;   TrTg = 1.021d0-0.084d0*(AM*((949.d0*pressure*1.d-5)+0.051d0))^(-0.5d0)
; 
; ; read in average dew point
;   readcol, path+loc1+'_DewPoint.txt', month_dt, dewT_all, format='a,d', skipline=3, /silent
; 
;   if (month_st eq '01' or month_st eq '02' or month_st eq '03') then dewT = dewT_all[0]
;   if (month_st eq '04' or month_st eq '05' or month_st eq '06') then dewT = dewT_all[1]
;   if (month_st eq '07' or month_st eq '08' or month_st eq '09') then dewT = dewT_all[2]
;   if (month_st eq '10' or month_st eq '11' or month_st eq '12') then dewT = dewT_all[3]
; 
;   ;convert into Fahrenheit
;   dewT = dewT*1.8d0 + 32.d0
; 
;   ;read in "lam", an empirically derived constant based on season and latitude (Table 2.1 in Smith 1966)
;   ;these valus are only valid for northern hemisphere but using the average value for the moment
;   ;readcol, 'fix/Smith1966.txt', ....
;   ;const = 2.61d0
; 
;   readcol, path+'Smith1966.txt', z1, z2, winter, spring, summer, fall, ann_ave, format='d,d,d,d,d,d,d,d', /silent
; 
;   idxz = where(abs(latdeg) le z2)
;   idxz = idxz[0]
; 
;   ;northern hemisphere
;   if (latdeg ge 0.) then begin
; 
;     if (month_st eq '01' or month_st eq '11' or month_st eq '12') then const = winter[idxz]
;     if (month_st eq '02' or month_st eq '03' or month_st eq '04') then const = spring[idxz]
;     if (month_st eq '05' or month_st eq '06' or month_st eq '07') then const = summer[idxz]
;     if (month_st eq '08' or month_st eq '09' or month_st eq '10') then const = fall[idxz]
; 
;   endif
; 
;   ;southern hemisphere
;   if (latdeg lt 0.) then begin
; 
;     if (month_st eq '01' or month_st eq '11' or month_st eq '12') then const = summer[idxz]
;     if (month_st eq '02' or month_st eq '03' or month_st eq '04') then const = fall[idxz]
;     if (month_st eq '05' or month_st eq '06' or month_st eq '07') then const = winter[idxz]
;     if (month_st eq '08' or month_st eq '09' or month_st eq '10') then const = spring[idxz]
; 
;   endif
; 
; 
; ; precipitable water vapor
;   u = exp(0.1133d0-alog(const+1.d0)+(0.0393d0*dewT))
; 
; ; water vapor transmission coefficient
;   Tw = 1.d0-0.077d0*(u*AM)^0.3
; 
; ; aerosol transmission coefficient
;   Ta = 0.935d0^AM
; 
;   sun_irr_tot = 1367.d0*cos(zd_sun*!dtor)*TrTg*Tw*Ta

;======================================================================================
