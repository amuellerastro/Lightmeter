pro sun_times, longdeg, latdeg, altitude, jd_day, jd, temp_kelvin, dt, t0set, t0rise, t6set, t6rise, t12set, t12rise, t18set, t18rise, $
		timesort_synth, h_sun_synth


;generate own time array (1sec resolution) in order to be independent of measurement gaps, all parameters have therefore the '_synth' because its not related to the measurements
  jd_synth = dindgen(86400.)/86400.d0 + jd_day
  ;caldat, double(jd_synth), month_synth, day_synth, year_synth, hh_synth, mm_synth, ss_synth
  ;hour_synth = hh_synth+mm_synth/60.d0+ss_synth/3600.d0
  hour_synth = dindgen(86400)/3600.d0

;sort time array
;   timesort_synth = hour_synth
;   idx1 = where(timesort_synth ge 12.)
;   idx2 = where(timesort_synth lt 12.)
;   if (idx1[0] ne -1) then timesort_synth[idx1] = timesort_synth[idx1]-12.
;   if (idx2[0] ne -1) then timesort_synth[idx2] = timesort_synth[idx2]+12.
  timesort_synth = hour_synth


  et_synth = jd_synth + (dt[0]/86400.d0)

  sunpos_mod, et_synth, ra_sun_synth, dec_sun_synth, longmed_synth, dist_ES_synth

  temp_kelvin_synth = interpol(temp_kelvin, jd, jd_synth, /spline)

  ;somtimes it happens that some elements are "NaN" after interpolation, don't know where this comes from 
  nan = finite(temp_kelvin_synth)
  idxnan = where(nan ne 1)
  if (idxnan[0] ne -1) then begin

    for i=0L,n_elements(idxnan)-1 do begin

      temp_kelvin_synth[idxnan[i]] = temp_kelvin_synth[idxnan[i]-1]

    endfor

  endif


  ;edges of temp-kelvin_synth have to be handled separately if data are missing, interpolation sucks there
  ;left side
  if (jd[0]-jd_synth[0] gt 1.d0/86400.d0) then begin

    idx = closest(jd[0], jd_synth, /upper)
    temp_kelvin_synth[0:idx-1] = temp_kelvin[0]

  endif

  ;right side
  if (jd_synth[n_elements(jd_synth)-1]-jd[n_elements(jd)-1] gt 1.d0/86400.d0) then begin

    idx = closest(jd[n_elements(jd)-1], jd_synth, /upper)
    temp_kelvin_synth[idx:n_elements(temp_kelvin_synth)-1] = temp_kelvin[n_elements(temp_kelvin)-1]

  endif


  eq2hor_mod, ra_sun_synth, dec_sun_synth, et_synth, h_sun_synth, az_sun_synth, ha_sun_synth, LON=longdeg, $
		  LAT=latdeg, precess_=1, nutate_=1, refract_=1, aberration_=1, altitude=altitude, temperature=temp_kelvin_synth

;MIGHT BE NOT VERY ROBUST FOR 'EXTREME' LONGITUDES

  ;sun set/rise: h="0" deg
  ;offset due to observers altitude
  h_offset0 = -(5.d0/6.d0)-(0.347d*sqrt(altitude/1000.d0))
  sun_set_0deg = hour_synth[closest(h_offset0, h_sun_synth[0:floor(n_elements(jd_synth)/2)])]
  sun_rise_0deg = hour_synth[floor(n_elements(jd_synth)/2)+closest(h_offset0, h_sun_synth[floor(n_elements(jd_synth)/2):*])]

;   if (sun_set_0deg ge 12.) then t0setdum = sun_set_0deg-12.
;   if (sun_set_0deg lt 12.) then t0setdum = sun_set_0deg+12.
;   t0set = t0setdum[0]
; 
;   if (sun_rise_0deg ge 12.) then t0risedum = sun_rise_0deg-12.
;   if (sun_rise_0deg lt 12.) then t0risedum = sun_rise_0deg+12.
;   t0rise = t0risedum[0]
  t0set = sun_set_0deg[0]
  t0rise = sun_rise_0deg[0]


  ;h_offset = -0.347d*sqrt(altitude/1000.d0) -- not sure if I need to use this correction for other elevations too
  h_offset = 0.d0

; use index to check if offset is bigger than XX deg than use mean sun coordinates!
; um wert zu bekommen: chck wieviel zeit diff ist bei bestimmter hoehe manuell

  ;civil twilight
  sun_set_6deg = hour_synth[closest(-6.d0+h_offset, h_sun_synth[0:floor(n_elements(jd_synth)/2)])]
  sun_rise_6deg = hour_synth[floor(n_elements(jd_synth)/2)+closest(-6.d0+h_offset, h_sun_synth[floor(n_elements(jd_synth)/2):*])]

;   if (sun_set_6deg ge 12.) then t6setdum = sun_set_6deg-12.
;   if (sun_set_6deg lt 12.) then t6setdum = sun_set_6deg+12.
;   t6set = t6setdum[0]
; 
;   if (sun_rise_6deg ge 12.) then t6risedum = sun_rise_6deg-12.
;   if (sun_rise_6deg lt 12.) then t6risedum = sun_rise_6deg+12.
;   t6rise = t6risedum[0]
  t6set = sun_set_6deg[0]
  t6rise = sun_rise_6deg[0]

  ;nautical twilight
  sun_set_12deg = hour_synth[closest(-12.d0+h_offset, h_sun_synth[0:floor(n_elements(jd_synth)/1.3d0)])]
  sun_rise_12deg = hour_synth[floor(n_elements(jd_synth)/1.3d0)+closest(-12.d0+h_offset, h_sun_synth[floor(n_elements(jd_synth)/1.3d0):*])]

;   if (sun_set_12deg ge 12.) then t12setdum = sun_set_12deg-12.
;   if (sun_set_12deg lt 12.) then t12setdum = sun_set_12deg+12.
;   t12set = t12setdum[0]
; 
;   if (sun_rise_12deg ge 12.) then t12risedum = sun_rise_12deg-12.
;   if (sun_rise_12deg lt 12.) then t12risedum = sun_rise_12deg+12.
;   t12rise = t12risedum[0]
  t12set = sun_set_12deg[0]
  t12rise = sun_rise_12deg[0]

  ;astronomical twilight
  sun_set_18deg = hour_synth[closest(-18.d0+h_offset, h_sun_synth[0:floor(n_elements(jd_synth)/1.3d0)])]
  sun_rise_18deg = hour_synth[floor(n_elements(jd_synth)/1.3d0)+closest(-18.d0+h_offset, h_sun_synth[floor(n_elements(jd_synth)/1.3d0):*])]

;   if (sun_set_18deg ge 12.) then t18setdum = sun_set_18deg-12.
;   if (sun_set_18deg lt 12.) then t18setdum = sun_set_18deg+12.
;   t18set = t18setdum[0]
; 
;   if (sun_rise_18deg ge 12.) then t18risedum = sun_rise_18deg-12.
;   if (sun_rise_18deg lt 12.) then t18risedum = sun_rise_18deg+12.
;   t18rise = t18risedum[0]
  t18set = sun_set_18deg[0]
  t18rise = sun_rise_18deg[0]

end