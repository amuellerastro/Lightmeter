pro twilight_model, pathbase, t0set, t0rise, t18set, t18rise, time_orig, counts_orig, temp_celsius_orig, h_sun_all, $
		     time_twi, counts_twi, temp_twi, h_sun_twi, twi_irr_tot

;path contains 'hard coded data'
  path = pathbase+'fix/'

;upper and lower sun altitude for twilight
  hlow = -16.d0
  hup = -6.d0

;rename
  time = time_orig
  counts = counts_orig
  temp_celsius = temp_celsius_orig
  h_sun = h_sun_all

;read in twilight model
  file = path+'twilight_brightness.txt'
  nlines = file_lines(file)
  openr,lun,file,/get_lun
    rows = nlines
    data = dblarr(3, rows)  
    READF, lun, data
    tmp_h = data(0,*)
    tmp_wm2 = data(1,*)
    tmp_lux = data(2,*)
    h_sun_model_orig = REFORM(tmp_h)
    wm2_model_orig = REFORM(tmp_wm2)
    lux_model_orig = REFORM(tmp_lux)
  close, lun
  free_lun, lun

  ;release some memory
  tmp_h = 0 & tmp_wm2 = 0 & tmp_lux = 0

;******************
;*twilight evening*
;******************

  idx = where(time ge t0set and time le t18set)
  ;idx = where(time ge t18rise and time le t0rise)
  counts_eve = counts[idx]
  h_sun_eve = h_sun[idx]
  time_eve = time[idx]
  temp_eve = temp_celsius[idx]

  idx = where(h_sun_eve le hup and h_sun_eve ge hlow)
  counts_eve = counts_eve[idx]
  h_sun_eve = h_sun_eve[idx]
  time_eve = time_eve[idx]
  temp_eve = temp_eve[idx]

;sort h_sun_model_orig in decreasing order, only necessary for evening twilight
  idxsort = bsort(h_sun_model_orig, /reverse)
  h_sun_model = h_sun_model_orig[idxsort]
  wm2_model = wm2_model_orig[idxsort]
  lux_model = lux_model_orig[idxsort]


;select model altitude depending on 'measured' ones, i.e., at given time of measurement
  ;accounts for time gaps
  idx = dblarr(n_elements(h_sun_eve))
  for i=0L,n_elements(h_sun_eve)-1 do begin
    idx[i] = closest(h_sun_eve[i], h_sun_model)
  endfor

  h_sun_model_eve = h_sun_model[idx]
  wm2_model_eve = wm2_model[idx]
  lux_model_eve = lux_model[idx]


;******************
;*twilight morning*
;******************

  idx = where(time ge t18rise and time le t0rise)
  counts_mor = counts[idx]
  h_sun_mor = h_sun[idx]
  time_mor = time[idx]
  temp_mor = temp_celsius[idx]

  idx = where(h_sun_mor ge hlow and h_sun_mor le hup)
  counts_mor = counts_mor[idx]
  h_sun_mor = h_sun_mor[idx]
  time_mor = time_mor[idx]
  temp_mor = temp_mor[idx]

  h_sun_model = h_sun_model_orig
  wm2_model = wm2_model_orig
  lux_model = lux_model_orig

;select model altitude depending on 'measured' ones, i.e., at given time of measurement
  ;accounts for time gaps
  idx = dblarr(n_elements(h_sun_mor))
  for i=0L,n_elements(h_sun_mor)-1 do begin
    idx[i] = closest(h_sun_mor[i], h_sun_model)
  endfor

  h_sun_model_mor = h_sun_model[idx]
  wm2_model_mor = wm2_model[idx]
  lux_model_mor = lux_model[idx]

  time_twi = [time_eve, time_mor]
  counts_twi = [counts_eve, counts_mor]
  temp_twi = [temp_eve, temp_mor]
  h_sun_twi = [h_sun_eve, h_sun_mor]
  twi_irr_tot = [wm2_model_eve, wm2_model_mor]

end