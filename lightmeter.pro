;own routines
@airmass_rozenberg.pro
@ave_night_sky_brightness.pro
@determine_time_offset.pro
@extinction.pro
@fit_global_irradiance_model.pro
@irradiance_statistic.pro
@lunar_model.pro
@lightmeter_lossamlike_moon.pro
@lightmeter_lossamlike_sun.pro
@night_log.pro
@median_fit_params.pro
@plot_lightcurve_single_night.pro
@plot_lightmeter_results,pro
@selen_coord_earth.pro
@selen_coord_sun.pro
@solar_model.pro
@sort_data_paranal.pro
@sort_data_armazones.pro
@sun_times.pro
;modified external routines
@co_aberration.pro
@co_nutate.pro
@co_refract_mod.pro
@eq2hor_mod.pro
@moonpos_mod.pro
@nutate_mod.pro
@robust_mean_mod.pro
@sunpos_mod.pro
;@twilight_model.pro
;external routines
@avg.pro
@bsort.pro
@caldat.pro
@cirrange.pro
@closest.pro
@convert_to_type.pro
@ct2lst.pro
@exponent.pro
@fsc_color.pro
@gettok.pro
@hadec2altaz.pro
@histoplot.pro
@interpol.pro
@isarray.pro
@julday.pro
@legend.pro
@logticks_exp.pro
@match2.pro
@mean.pro
@moment.pro
@mpfitfun.pro
@mpfit.pro
@multiplot.pro
@nutate.pro
@oploterror.pro
@plothist.pro
@poly.pro
@precess.pro
@premat.pro
@readcol.pro
@recpol.pro
@remchar.pro
@showsym.pro
@sigfig.pro
@sixlin.pro
@stddev.pro
@strnumber.pro
@strsplit.pro
@sunpos.pro
@ten.pro
@ts_diff.pro
@tsum.pro
@uniq.pro


pro lightmeter


;select location
  readcol, 'fix/LOCATION.txt', loc2, loc1, long1, long2, long3, lat1, lat2, lat3, altitude, format='a,a,d,d,d,d,d,d,d', /silent

  nlocs = indgen(n_elements(loc1))+1
  print, ''
  print, ''
  for i=0,n_elements(loc1)-1 do begin
    print, nlocs[i], '  ', loc1[i]
  endfor

  print, ''
  read, 'Select Location by Number: ',num
  num = num-1

  loc1 = loc1[num]
  loc2 = loc2[num]
  long1 = long1[num]
  long2 = long2[num]
  long3 = long3[num]
  lat1 = lat1[num]
  lat2 = lat2[num]
  lat3 = lat3[num]
  altitude = altitude[num]

  longdeg = ten(long1,long2,long3)
  latdeg = ten(lat1,lat2,lat3)


;sort raw data or reduce sorted data?
  quest = ''
  print, ''
  read, 'Sort Raw Data (s) / Reduce Data (r): ', quest

  if (quest ne 's' and quest ne 'r') then begin

    print, 'Enter correct option. Stop.'
    return

  endif


if (quest eq 's') then begin

  if (loc1 eq 'Armazones') then sort_data_armazones
  if (loc1 eq 'Paranal') then sort_data_paranal

endif


if (quest eq 'r') then begin


  ;read in setup file
    readcol, 'fix/SETUP.txt', dum, var, format='a,a', skipline=2, /silent

    pathbase = var[0]

    pathsort = pathbase+'data_sorted/'+loc1+'/'	;directory containing the data of each night, if processed it get moved to $pathproc

    pathproc = pathsort+'processed/'	;directory containing the processed nights
    file_mkdir, pathproc

    pathresult = pathbase+loc1+'/'	;directory containing plots and analysis
    file_mkdir, pathresult

    pathlightcurve = pathresult+'Lightcurve/'
    file_mkdir, pathlightcurve

    pathlossamlike = pathresult+'LOSSAMlike/'
    file_mkdir, pathlossamlike

    pathmodelfit = pathresult+'IrradianceModelFit/'
    file_mkdir, pathmodelfit

    pathnsb = pathresult+'NightSkyBrightness/'
    file_mkdir, pathnsb

    pathstat = pathresult+'IrradianceStatistic/'
    file_mkdir, pathstat

    pathlog = pathresult+'NightLog/'
    file_mkdir, pathlog


  ;read in time correction
    readcol, 'fix/SETUP.txt', dum, dT, format='a,d', skipline=7, /silent


  ;read in starting values for the model fits
    readcol, 'fix/SETUP.txt', loc1dum, irrad_case, p1dum, p2dum, p3dum, p4dum, format='a,a,d,d,d,d', skipline=9, /silent
    idx = where(loc1 eq loc1dum)
    if (idx[0] eq -1) then begin

      print, 'Check setup file. Aborting here.'
      return

    endif else begin

      loc1dum = loc1dum[idx]
      irrad_case = irrad_case[idx]
      p1dum = p1dum[idx]
      p2dum = p2dum[idx]
      p3dum = p3dum[idx]
      p4dum = p4dum[idx]

      start_params_sun = [p1dum[0],p2dum[0],p3dum[0],p4dum[0]]
      start_params_sunmoon = [p1dum[1],p2dum[1],p3dum[1],p4dum[1]]

    endelse


  ;search files
    file = file_search(pathsort+loc2+'_'+loc1+'_*.txt', count=nfiles)
    print, ''
    print, 'processing '+strtrim(string(nfiles),2)+' nights...'
    print, ''

  ;flag for time offset computation 0 - off / 1 - on
    to = 0
    print, ''
    if (to eq 0) then print, 'time offset computation disabled'
    if (to eq 1) then print, 'time offset computation enabled'

  for i=0L,nfiles-1 do begin

;time = systime(/sec)

    print, ''
    print, 'processing night '+strtrim(string(i+1),2)+' of '+strtrim(string(nfiles),2)

    ;read in data of each night
    print, 'read in Lightmeter data...'
    nlines = file_lines(file[i])
    openr,lun,file[i],/get_lun
      rows = nlines
      data = dblarr(3, rows)  
      READF, lun, data
      tmp_jd = data(0,*)
      tmp_temp = data(1,*)
      tmp_counts = data(2,*)
      jd_orig = REFORM(tmp_jd)
      temp_celsius_orig = REFORM(tmp_temp)
      counts_orig = REFORM(tmp_counts)
    close, lun
    free_lun, lun

    ;release some memory
    tmp_jd = 0 & tmp_temp = 0 & tmp_counts = 0

    temp_kelvin_orig = temp_celsius_orig+273.15d0

  ;due to incomplete datasets the date of the night has to be extracted from the file name!
    pos1 = strpos(file[i], '_', /reverse_search)
    pos2 = strpos(file[i], '.', /reverse_search)
    tmp_date = strmid(file[i],pos1+1, pos2-pos1-1)
    year_st = strmid(tmp_date, 0, 4)
    month_st = strmid(tmp_date, 4, 2)
    day_st = strmid(tmp_date, 6, 2)

  ;JD at 12:00 a.m. UT which gets processed
    jd_day = long(julday(double(month_st), double(day_st), double(year_st), 12., 0., 0.))



;simple check if data for all day are available (12:00:00 p.m. UT to 11:59:59 a.m. UT) -- MIGHT BE NOT VERY ROBUST
    ;compute median sampling rate in seconds
    sampling = abs(median(ts_diff(jd_orig,1))*86400.d0)

    ;one file has 24hrs coverage, check the number of data points giving the sampling rate
    n_theory = 86400.d0/sampling	;total number of points theoretically recorded

    ;set the threshold of missing data to 10%
    n_theory_missing = 0.9d0*n_theory


    if (n_elements(jd_orig) ge n_theory_missing) then begin


;**************************
;*SUN RISING/SETTING TIMES*
;**************************

    ;compute civil (sun -6deg), nautical (sun -12deg) and astronomical (sun -18deg) twilight, sun set and rise (50" below horizon)
      print, 'computing rising and setting times of the sun...'
      sun_times, longdeg, latdeg, altitude, jd_day, jd_orig, temp_kelvin_orig, dt, t0set, t0rise, t6set, t6rise, t12set, t12rise, t18set, t18rise, timesort_synth, h_sun_synth


;*************
;*TIME OFFSET*
;*************

    if (to eq 1) then begin

    ;get gregorian date from JD
      caldat, double(jd_orig), month, day, year, hh, mm, ss
      hour_orig = hh+mm/60.d0+double(floor(ss))/3600.d0

    ;sort time in a proper way
      timesort_orig = hour_orig
      idx1_orig = where(timesort_orig ge 12.)
      idx2_orig = where(timesort_orig lt 12.)
      if (idx1_orig[0] ne -1) then timesort_orig[idx1_orig] = timesort_orig[idx1_orig]-12.
      if (idx2_orig[0] ne -1) then timesort_orig[idx2_orig] = timesort_orig[idx2_orig]+12.


    ;lossam like data for day time for time offsets
      print, 'computing lossam like data for day time for time offset computation...'
      lightmeter_lossamlike_sun, t0set, t0rise, timesort_orig, counts_orig, $
				sun_binsize, sun_binsize10, sun_hrbin, sun_rvarcounts

    ;solar model
      print, 'computing solar model for time offset computation...'
      solar_model, pathbase, altitude, t0set, t0rise, timesort_orig, counts_orig, temp_celsius_orig, $
		   timesort_synth=timesort_synth, h_sun_synth, sun_binsize, sun_binsize10, sun_hrbin, sun_rvarcounts, $
		   time_sun_tmp, counts_sun_tmp, temp_sun_tmp, h_sun_tmp, sun_irr_tot_tmp	;already selected for h>20deg
 
    ;check for time shifts
      print, 'check for time shifts...will take a while...'
      determine_time_offset, time_sun_tmp, counts_sun_tmp, temp_sun_tmp, h_sun_tmp, sun_irr_tot_tmp, start_params_sun, $
			    t_offset

      jd = jd_orig + t_offset

      ;account for date jumps due to offset correction
      idx = where(jd ge jd_orig[0] and jd le ceil(jd_orig[n_elements(jd_orig)-1]))
      jd = jd[idx]
      counts = counts_orig[idx]
      temp_celsius = temp_celsius_orig[idx]
      temp_kelvin = temp_kelvin_orig[idx]

    endif else begin

      jd = jd_orig
      counts = counts_orig
      temp_celsius = temp_celsius_orig
      temp_kelvin = temp_kelvin_orig

      t_offset = 0.	;dummy

    endelse


    ;get gregorian date from JD
      caldat, double(jd), month, day, year, hh, mm, ss
      hour = hh+mm/60.d0+double(floor(ss))/3600.d0

    ;sort time in a proper way
      timesort = hour
      idx1 = where(timesort ge 12.)
      idx2 = where(timesort lt 12.)
      if (idx1[0] ne -1) then timesort[idx1] = timesort[idx1]-12.
      if (idx2[0] ne -1) then timesort[idx2] = timesort[idx2]+12.

      print, 'computing solar epehemerides...'
      et = jd + (dt[0]/86400.d0)
      sunpos_mod, et, ra_sun, dec_sun, longmed, dist_ES
      eq2hor_mod, ra_sun, dec_sun, et, h_sun_all, az_sun, ha_sun, LON=longdeg, LAT=latdeg, precess_=1, nutate_=1, $
		  refract_=1, aberration_=1, altitude=altitude, temperature=temp_kelvin

;*************************
;*GLOBAL IRRADIANCE MODEL*
;*************************

      ;lossam like data for day time for time offsets
	print, 'computing lossam like data for day time for global model...'
	lightmeter_lossamlike_sun, t0set, t0rise, timesort, counts, altsun=h_sun, $
				  sun_binsize, sun_binsize10, sun_hrbin, sun_rvarcounts

      ;solar model
	print, 'computing solar model for global model...'
	solar_model, pathbase, altitude, t0set, t0rise, timesort, counts, temp_celsius, h_sun_all, sun_binsize, sun_binsize10, sun_hrbin, sun_rvarcounts, $
		     time_sun, counts_sun, temp_sun, h_sun, sun_irr_tot

;       ;twilight model
; 	print, 'computing twilight model for global model...'
; 	twilight_model, pathbase, t0set, t0rise, t18set, t18rise, timesort, counts, temp_celsius, h_sun_all, $
; 		     time_twi, counts_twi, temp_twi, h_sun_twi, twi_irr_tot

      ;create plot and ASCII file of lightmeter data handled like LOSSAM
	print, 'computing lossam like data of night...'
	lightmeter_lossamlike_moon, pathlossamlike, loc1, loc2, hour, day_st, month_st, year_st, counts, t18set, t18rise, $
				    moon_binsize, moon_hrbin, moon_rvarcounts

      ;lunar model
	print, 'computing lunar model...'
	lunar_model, pathbase, pathlunfit, longdeg, latdeg, altitude, jd_day, et, jd, timesort, longmed, dist_ES, $
		    t18set, t18rise, counts, temp_kelvin, temp_celsius, moon_binsize, moon_hrbin, moon_rvarcounts, $
		    h_moon, h_moon_orig, illfrac_moon, phaseangle_moon, time_moon, counts_moon, moon_irr_tot, temp_moon, check

      ;global (Sun+Moon) irradiance Model
	print, 'fit global irradiance model...will take a while'
	fit_global_irradiance_model, loc1, time_sun, counts_sun, temp_sun, sun_irr_tot, $
				      time_moon, counts_moon, temp_moon, moon_irr_tot, check, $;time_twi, counts_twi, temp_twi, h_sun_twi, twi_irr_tot, $
				      start_params_sun, start_params_sunmoon, $
				      parinfo, fit_params, fit_status, fit_rms, fit_bestnorm	;bestnorm = the value of the summed squared residuals for the returned parameter values.


      ;write out results file of fit parameters
	openw, lun, pathmodelfit+'IrradianceModelFit.txt', width=2500, /get_lun, /append
		  ;  JD      Lunar phase angle             IllumFrac       Moon used?         fitted parameters
	  printf, lun, jd_day, mean(phaseangle_moon/!dtor),mean(illfrac_moon), check, fit_params[0], fit_params[1], $
			;                              RMS      SSR           MPfitStatus  time offset
			fit_params[2], fit_params[3], fit_rms, fit_bestnorm, fit_status, t_offset*86400.d0, $
			format='(i7, f14.8, f13.8, a5, 4e18.9, f13.8, e18.9, i4, f10.1)'

	close, lun
	free_lun, lun


;*************
;*SOME OUTPUT*
;*************

      ;check if fit parameters already exist, if not, use estimates from the setup file
      fitfile = file_search(pathmodelfit+'IrradianceModelFit.txt', count=nfile)

      if (nfile eq 1) then begin

	median_fit_params, pathmodelfit, start_params_sunmoon, a, delta_a, b, delta_b, c, delta_c, d, delta_d
	irradiance = c*(b*(a*exp((counts*(1.d0+d*temp_celsius))/a)-1.d0)+(counts*(1.d0+d*temp_celsius)))

      endif else begin

	a = start_params_sunmoon[0]
	b = start_params_sunmoon[1]
	c = start_params_sunmoon[2]
	d = start_params_sunmoon[3]

	irradiance = c*(b*(a*exp((counts*(1.d0+d*temp_celsius))/a)-1.d0)+(counts*(1.d0+d*temp_celsius)))

      endelse

    ;create plot of light curve of one night
      print, 'creating plot of lightcurve...'
      plot_lightcurve_single_night, loc1, loc2, pathlightcurve, day_st, month_st, year_st, timesort, irradiance, $
				    t0rise, t0set, t6rise, t6set, t12rise, t12set, t18rise, t18set

    ;average night sky brightness (without moon) ans h_sun < -18deg
      print, 'average night sky brightness...'	;lossam file has to be read in
      ave_night_sky_brightness, pathnsb, jd_day, timesort, irradiance, t18set, t18rise, h_moon_orig, moon_binsize, moon_hrbin, moon_rvarcounts

    ;counts statistic (histogram)
      print, 'irradiance statistics...'
      irradiance_statistic, loc1, loc2, pathstat, timesort, irradiance, jd_day, t18rise, t18set, moon_binsize, moon_hrbin, moon_rvarcounts

    ;night log
      print, 'writing night log...'
      night_log, pathlog, loc1, loc2, day_st, month_st, year_st, t18set, t18rise, jd, hour, timesort, irradiance, temp_celsius, $
		 temp_kelvin, h_moon_orig, moon_binsize, moon_hrbin, moon_rvarcounts, illfrac_moon



; print, 'took ', systime(/sec)-time, ' sec'

    endif else begin	;condition for sufficient number of data points

      print, 'No complete data set available.'

    ;compute civil (sun -6deg), nautical (sun -12deg) and astronomical (sun -18deg) twilight, sun set and rise (50" below horizon)
      print, 'computing rising and setting times of the sun...'
      sun_times, longdeg, latdeg, altitude, jd_day, jd_orig, temp_kelvin_orig, dt, t0set, t0rise, t6set, t6rise, t12set, t12rise, t18set, t18rise, timesort_synth, h_sun_synth


;*************
;*TIME OFFSET*
;*************

    if (to eq 1) then begin

    ;get gregorian date from JD
      caldat, double(jd_orig), month, day, year, hh, mm, ss
      hour_orig = hh+mm/60.d0+double(floor(ss))/3600.d0

    ;sort time in a proper way
      timesort_orig = hour_orig
      idx1_orig = where(timesort_orig ge 12.)
      idx2_orig = where(timesort_orig lt 12.)
      if (idx1_orig[0] ne -1) then timesort_orig[idx1_orig] = timesort_orig[idx1_orig]-12.
      if (idx2_orig[0] ne -1) then timesort_orig[idx2_orig] = timesort_orig[idx2_orig]+12.


    ;lossam like data for day time for time offsets
      print, 'computing lossam like data for day time for time offset computation...'
      lightmeter_lossamlike_sun, t0set, t0rise, timesort_orig, counts_orig, $
				sun_binsize, sun_binsize10, sun_hrbin, sun_rvarcounts

    ;solar model
      print, 'computing solar model for time offset computation...'
      solar_model, pathbase, altitude, t0set, t0rise, timesort_orig, counts_orig, temp_celsius_orig, $
		   timesort_synth=timesort_synth, h_sun_synth, sun_binsize, sun_binsize10, sun_hrbin, sun_rvarcounts, $
		   time_sun_tmp, counts_sun_tmp, temp_sun_tmp, h_sun_tmp, sun_irr_tot_tmp	;already selected for h>20deg
 
    ;check for time shifts
      print, 'check for time shifts...will take a while...'
      determine_time_offset, time_sun_tmp, counts_sun_tmp, temp_sun_tmp, h_sun_tmp, sun_irr_tot_tmp, start_params_sun, $
			    t_offset

      jd = jd_orig + t_offset

      ;account for date jumps due to offset correction
      idx = where(jd ge jd_orig[0] and jd le ceil(jd_orig[n_elements(jd_orig)-1]))
      jd = jd[idx]
      counts = counts_orig[idx]
      temp_celsius = temp_celsius_orig[idx]
      temp_kelvin = temp_kelvin_orig[idx]

    endif else begin

      jd = jd_orig
      counts = counts_orig
      temp_celsius = temp_celsius_orig
      temp_kelvin = temp_kelvin_orig

      t_offset = 0.	;dummy

    endelse


    ;get gregorian date from JD
      caldat, double(jd), month, day, year, hh, mm, ss
      hour = hh+mm/60.d0+double(floor(ss))/3600.d0

    ;sort time in a proper way
      timesort = hour
      idx1 = where(timesort ge 12.)
      idx2 = where(timesort lt 12.)
      if (idx1[0] ne -1) then timesort[idx1] = timesort[idx1]-12.
      if (idx2[0] ne -1) then timesort[idx2] = timesort[idx2]+12.

      idxdata = where(timesort ge t18set and timesort le t18rise)
      if (idxdata[0] ne -1) then begin

	print, 'computing solar epehemerides...'
	et = jd + (dt[0]/86400.d0)
	sunpos_mod, et, ra_sun, dec_sun, longmed, dist_ES
	eq2hor_mod, ra_sun, dec_sun, et, h_sun_all, az_sun, ha_sun, LON=longdeg, LAT=latdeg, precess_=1, nutate_=1, $
		    refract_=1, aberration_=1, altitude=altitude, temperature=temp_kelvin

	;create plot and ASCII file of lightmeter data handled like LOSSAM
	  print, 'computing lossam like data of night...'
	  lightmeter_lossamlike_moon, pathlossamlike, loc1, loc2, hour, day_st, month_st, year_st, counts, t18set, t18rise, $
				      moon_binsize, moon_hrbin, moon_rvarcounts

	;lunar model
	  print, 'computing lunar model...'
	  lunar_model, pathbase, pathlunfit, longdeg, latdeg, altitude, jd_day, et, jd, timesort, longmed, dist_ES, $
		      t18set, t18rise, counts, temp_kelvin, temp_celsius, moon_binsize, moon_hrbin, moon_rvarcounts, $
		      h_moon, h_moon_orig, illfrac_moon, phaseangle_moon, time_moon, counts_moon, moon_irr_tot, temp_moon, check


  ;*************
  ;*SOME OUTPUT*
  ;*************

	;check if fit parameters already exist, if not, use estimates from the setup file
	fitfile = file_search(pathmodelfit+'IrradianceModelFit.txt', count=nfile)

	if (nfile eq 1) then begin

	  median_fit_params, pathmodelfit, start_params_sunmoon, a, delta_a, b, delta_b, c, delta_c, d, delta_d
	  irradiance = c*(b*(a*exp((counts*(1.d0+d*temp_celsius))/a)-1.d0)+(counts*(1.d0+d*temp_celsius)))

	endif else begin

	  a = start_params_sunmoon[0]
	  b = start_params_sunmoon[1]
	  c = start_params_sunmoon[2]
	  d = start_params_sunmoon[3]

	  irradiance = c*(b*(a*exp((counts*(1.d0+d*temp_celsius))/a)-1.d0)+(counts*(1.d0+d*temp_celsius)))

	endelse

      ;create plot of light curve of one night
	print, 'creating plot of lightcurve...'
	plot_lightcurve_single_night, loc1, loc2, pathlightcurve, day_st, month_st, year_st, timesort, irradiance, $
				      t0rise, t0set, t6rise, t6set, t12rise, t12set, t18rise, t18set

      ;average night sky brightness (without moon) ans h_sun < -18deg
	print, 'average night sky brightness...'	;lossam file has to be read in
	ave_night_sky_brightness, pathnsb, jd_day, timesort, irradiance, t18set, t18rise, h_moon_orig, moon_binsize, moon_hrbin, moon_rvarcounts

      ;counts statistic (histogram)
	print, 'irradiance statistics...'
	irradiance_statistic, loc1, loc2, pathstat, timesort, irradiance, jd_day, t18rise, t18set, moon_binsize, moon_hrbin, moon_rvarcounts

      ;night log
	print, 'writing night log...'
	night_log, pathlog, loc1, loc2, day_st, month_st, year_st, t18set, t18rise, jd, hour, timesort, irradiance, temp_celsius, $
		  temp_kelvin, h_moon_orig, moon_binsize, moon_hrbin, moon_rvarcounts, illfrac_moon


  ;       openw, lun, pathresult+'not_processed.txt', width=1400, /get_lun, /append
  ; 	printf, lun, loc2+'_'+loc1+'_'+year_st+month_st+day_st
  ;       close, lun
  ;       free_lun, lun

      endif

    endelse

    spawn, 'mv '+file[i]+' '+pathproc

  endfor	;loop over nights


;plot several results
  print, ''
  print, 'creating several plots...'
  print, ''
  plot_lightmeter_results, loc1, loc2, pathlossamlike, pathlunfit, pathnsb, pathstat


endif	;reduction / sorting of data



end