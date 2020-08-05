pro ave_night_sky_brightness, pathnsb, jd_day, time, counts_wm2, t18set, t18rise, hm, binsize, hrbin, rvarcounts


  idxnight = where(time ge t18set and time le t18rise)

  if (idxnight[0] ne -1) then begin

    hr_dark = time[idxnight]
    counts_wm2_dark = counts_wm2[idxnight]
    hm_dark = hm[idxnight]

  ;consider only times where moon is below -Xdeg
    idxmoon = where(hm_dark le -5.)

    if (idxmoon[0] ne -1) then begin

      hr_nomoon = hr_dark[idxmoon]
      counts_wm2_nomoon = counts_wm2_dark[idxmoon]

      ;remove the first/last 2 hours of night because of zodiacal light
      idxzod = where(hr_nomoon ge (t18set+2.d0) and hr_nomoon le (t18rise-2.d0))
      if (idxzod[0] ne -1) then begin

	hr_nomoon_zod = hr_nomoon[idxzod]
	counts_wm2_nomoon_zod = counts_wm2_nomoon[idxzod]

	mean_counts_wm2_zod = mean(counts_wm2_nomoon_zod)
	median_counts_wm2_zod = median(counts_wm2_nomoon_zod, /even)
	stddev_counts_wm2_zod = stddev(counts_wm2_nomoon_zod)

	openw, lun, pathnsb+'MedianNightSkyBrightness_Wm2.txt', width=1400, /get_lun, /append
	  printf, lun, jd_day, median_counts_wm2_zod, stddev_counts_wm2_zod, format='(i7, 2e16.7)'
	close, lun
	free_lun, lun

;ONLY DARK TIME WITHOUT CLOUDS
	;cut out areas having clouds
	idx_clean = where(rvarcounts lt 0.02d0)
	if (idx_clean[0] ne -1) then begin

	  hrbin_clean = hrbin[idx_clean]
	  rvarcounts_clean = rvarcounts[idx_clean]

	endif else begin

	  return  ;i.e. if whole night is cloudy do nothing
	;     hrbin_clean = hrbin
	;     rvarcounts_clean = rvarcounts

	endelse

	;create time array in order to match hr_nomoon_zod with the non cloudy parts of the night
	  tleft = hrbin_clean - binsize[0]/(2.d0*3600.d0)
	  tright = hrbin_clean + binsize[0]/(2.d0*3600.d0)

	counts_wm2_tmp = dblarr(n_elements(counts_wm2_nomoon_zod)) & time_tmp = counts_wm2_tmp

	for i=0L,n_elements(tleft)-1 do begin

	  idx = where(hr_nomoon_zod ge tleft[i] and hr_nomoon_zod le tright[i])

	  if (idx[0] ne -1) then begin

	    idx0 = where(counts_wm2_tmp eq 0.)
	    counts_wm2_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = counts_wm2_nomoon_zod[idx]
	    time_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = hr_nomoon_zod[idx]

	  endif

	endfor

	idxval = where(counts_wm2_tmp ne 0.)
	counts_wm2_clean = counts_wm2_tmp[idxval]
	time_clean = time_tmp[idxval]

	mean_counts_wm2_clean = mean(counts_wm2_clean)
	median_counts_wm2_clean = median(counts_wm2_clean, /even)
	stddev_counts_wm2_clean = stddev(counts_wm2_clean)


	openw, lun, pathnsb+'MedianNightSkyBrightness_Wm2_NoClouds.txt', width=1400, /get_lun, /append
	  printf, lun, jd_day, median_counts_wm2_clean, stddev_counts_wm2_clean, format='(i7, 2e16.7)'
	close, lun
	free_lun, lun


      endif

    endif

  endif

end