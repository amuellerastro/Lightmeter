pro night_log, pathlog, loc1, loc2, day_st, month_st, year_st, t18set, t18rise, jd_orig, hour_orig, timesort_orig, irradiance_orig, temp_celsius_orig, $
		temp_kelvin, h_moon_orig, moon_binsize, moon_hrbin_orig, moon_rvarcounts, illfrac_moon


;length of night (dark time) sun < -18deg
  st_nightlength = sigfig(t18rise-t18set,4)	;string

;time of astronomical twilight in UT
  if (t18set lt 12.) then st_t18set = sigfig(t18set+12.,5) else st_t18set = sigfig(t18set-12.,5)
  if (t18rise lt 12.) then st_t18rise = sigfig(t18rise+12.,5) else st_t18rise = sigfig(t18rise-12.,5)

;moon illumination fraction
  st_illfrac_moon = sigfig(mean(illfrac_moon),3)



;consider only astronomical night, i.e., sun < -18 deg
  idx = where(timesort_orig ge t18set and timesort_orig le t18rise)
    jd = jd_orig[idx]
    hour = hour_orig[idx]
    timesort = timesort_orig[idx]
    irradiance = irradiance_orig[idx]
    temp_celsius = temp_celsius_orig[idx]
    h_moon = h_moon_orig[idx]


;extract lossam data
  moon_hrbin = moon_hrbin_orig
;   idx1 = where(moon_hrbin ge 12.)
;   idx2 = where(moon_hrbin lt 12.)
;   if (idx1[0] ne -1) then moon_hrbin[idx1] = moon_hrbin[idx1]-12.
;   if (idx2[0] ne -1) then moon_hrbin[idx2] = moon_hrbin[idx2]+12.


  tleft = moon_hrbin - moon_binsize[0]/(2.d0*3600.d0)
  tright = moon_hrbin + moon_binsize[0]/(2.d0*3600.d0)

  time_tmp = dblarr(n_elements(hour)) & lossam_tmp = time_tmp
  for i=0L,n_elements(tleft)-1 do begin

    idx = where(timesort ge tleft[i] and timesort le tright[i])

    if (idx[0] ne -1) then begin

	idx0 = where(time_tmp eq 0.)
	time_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = timesort[idx]
	lossam_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = moon_rvarcounts[i]

      endif

    endfor

    idxval = where(time_tmp ne 0.)
    if (idxval[0] ne -1) then begin

      time_tmp = time_tmp[idxval]
      lossam_tmp = lossam_tmp[idxval]

      time_lossam = time_tmp
      lossam = lossam_tmp

    endif


;sinchronize with lossam data
  idx = where(timesort ge time_lossam[0] and timesort le time_lossam[n_elements(time_lossam)-1])
    jd = jd[idx]
    hour = hour[idx]
    timesort = timesort[idx]
    irradiance = irradiance[idx]
    temp_celsius = temp_celsius[idx]
    h_moon = h_moon[idx]


openw, lun, pathlog+loc2+'_'+loc1+'_NightLog_'+strcompress(year_st,/rem)+strcompress(month_st,/rem)+strcompress(day_st,/rem)+'.txt', /get_lun

  printf, lun, day_st+'.'+month_st+'.'+year_st
  printf, lun, ''
  printf, lun, 'Night length (H_sun < -18 deg): '+st_nightlength+' hr'
  printf, lun, 'Illumintaed fraction of the Moon: '+st_illfrac_moon
  printf, lun, ''
  printf, lun, ''


  printf, lun, 'JD                    UT          Irradiance     Rel.RMS   Temp.  H_moon'
  printf, lun, '[days]                [hr]        [W/m^2]                  [degC] [deg]'
  for i=0L, n_elements(timesort)-1 do begin

    printf, lun, jd[i], hour[i], irradiance[i], lossam[i], temp_celsius[i], h_moon[i], format='(f18.10,f12.5,e15.5,f10.4,f8.1,f8.2)'

  endfor


close, lun
free_lun, lun


end