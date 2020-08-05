pro lightmeter_lossamlike_sun, t0set, t0rise, time_tmp, counts_tmp, altsun=altsun, binsize, binsize10, hrbin_sun, rvarcounts_sun

; binsize = 10.d0
; binsize10 = 600.d0

;read in setup file
  readcol, 'fix/SETUP.txt', dum, binsize, binsize10, format='a,d,d', skipline=5, /silent

  binsize = binsize[0]
  binsize10 = binsize10[0]

;SUN BEFORE NIGHT
  counts = counts_tmp
  hr = time_tmp
;  h_sun = h_sun_tmp


;altitude selection

  if keyword_set(altsun) then begin

    idx_sun = where(altsun ge 20.d0)
    counts = counts[idx_sun]
    hr = hr[idx_sun]

  endif

;time selection
  idx = where(hr le t0set)

  if (idx[0] ne -1 and n_elements(idx) gt binsize10) then begin

    check1 = 'yes'
    hr = hr[idx]
    counts = counts[idx]

    diff_all = abs(ts_diff(hr, 1))*3600.d0	;difference to the next time stamp
    ave_samp = median(diff_all)	;mean sampling rate
    ave_samp = round(median(diff_all[where(diff_all lt 1.5*ave_samp)]))	;just to make it more reliable

    ;if there are gaps in time fill them
    idx_tmp = where(diff_all gt 1.5*ave_samp)
    if (idx_tmp[0] ne -1) then begin

      nts = (hr[n_elements(hr)-1] - hr[0])*3600.d0*ave_samp
      time = dblarr(nts+2)	;new time array, number of elements is the same if there would be no time gaps
      counts_new = time

      ;new time base if there would be no gaps in time
      time[*] = (dindgen(nts+2)*ave_samp/3600.d0)+hr[0]
      idxval = dblarr(n_elements(hr))
      for i=0L,n_elements(hr)-1 do begin

	idxval[i] = where(round(time*3600.d0) eq round(hr[i]*3600.d0))	;necessary because of internal precision
	counts_new[idxval[i]] = counts[i]

      endfor

    ;fill time gaps with counts from previous measurement
      for i=0L,n_elements(time)-1 do begin

	idx0 = where(counts_new eq 0.d0)
	if (idx0[0] ne -1) then counts_new[idx0] = counts_new[idx0-1] else break

      endfor


    ;find out where are the first, e.g. 10 min (binsize10) of a night
    ;starting computation from that point
      diff = time-time[0]
      index10 = closest(binsize10/3600.d0, diff, value=value, /upper)	;starting index
      nbins = fix(((time[n_elements(time)-1]-time[index10[0]])/binsize)*3600.d0)

      counts = counts_new
      hr = time

    ;array definition
      ;time bin belongs to 'binsize', i.e. small bin
      hrbin = dblarr(nbins) & hrdum = hrbin & tmp_hr = hr[index10:n_elements(hr)-1] & countsbin10_rms = hrbin & countsbin10_mean = hrbin



    endif else begin
    ;find out where are the first, e.g. 10 min (binsize10) of a night
    ;starting computation from that point
      diff = hr-hr[0]
      index10 = closest(binsize10/3600.d0, diff, value=value, /upper)	;starting index
      nbins = fix(((hr[n_elements(hr)-1]-hr[index10[0]])/binsize)*3600.d0)

    ;array definition
      ;time bin belongs to 'binsize', i.e. small bin
      hrbin = dblarr(nbins) & hrdum = hrbin &  tmp_hr = hr[index10:n_elements(hr)-1] & countsbin10_rms = hrbin & countsbin10_mean = hrbin

    endelse

  ;first element needs to be handled separately
    idx_1st = closest(diff[index10]-binsize/3600.d0, diff[0:index10], /decide)
    hrbin[0] = mean(hr[idx_1st:index10])
    countsbin10_mean[0] = mean(counts[0:index10])
    countsbin10_rms[0] = stddev(counts[0:index10])

    for i=1L,nbins[0]-1 do begin

    ;small bin
      diff = tmp_hr - tmp_hr[0]
      idx = closest(binsize/3600.d0, diff, /decide)

      if (idx[0] eq 0) then begin

	print, ''
	print, 'Binsize is smaller than sampling. Can be caused by data gaps. Stop.'
	print, ''
	stop

      endif

      hrbin[i] = mean(tmp_hr[0:idx-1])
      hrdum[i] = tmp_hr[idx-1]

      tmp_hr = tmp_hr[idx:n_elements(tmp_hr)-1]

    ;big bin
      idx1 = where(hrdum[i] eq hr)
      tdum = hr[idx1]-(binsize10/3600.d0)
      idx2 = closest(tdum, hr, /decide)

      countsbin10_mean[i] = mean(counts[idx2:idx1])
      countsbin10_rms[i] = stddev(counts[idx2:idx1])

    endfor

    rvarcounts_1 = countsbin10_rms/countsbin10_mean
    hrbin_1 = hrbin

  endif else begin

    check1 = 'no'

  endelse



;SUN AFTER NIGHT

  ;check2 = 'yes'
  counts = counts_tmp
  hr = time_tmp
  idx = where(hr ge t0rise)

  if (idx[0] ne -1 and n_elements(idx) gt binsize10) then begin

    check2 = 'yes'
    hr = hr[idx]
    counts = counts[idx]

    diff_all = abs(ts_diff(hr, 1))*3600.d0	;difference to the next time stamp
    ave_samp = median(diff_all)	;mean sampling rate
    ave_samp = round(median(diff_all[where(diff_all lt 1.5*ave_samp)]))	;just to make it more reliable

    ;if there are gaps in time fill them
    idx_tmp = where(diff_all gt 1.5*ave_samp)
    if (idx_tmp[0] ne -1) then begin

      nts = (hr[n_elements(hr)-1] - hr[0])*3600.d0*ave_samp
      time = dblarr(nts+2)	;new time array, number of elements is the same if there would be no time gaps
      counts_new = time

      ;new time base if there would be no gaps in time
      time[*] = (dindgen(nts+2)*ave_samp/3600.d0)+hr[0]
      idxval = dblarr(n_elements(hr))
      for i=0L,n_elements(hr)-1 do begin

	idxval[i] = where(round(time*3600.d0) eq round(hr[i]*3600.d0))	;necessary because of internal precision
	counts_new[idxval[i]] = counts[i]

      endfor

    ;fill time gaps with counts from previous measurement
      for i=0L,n_elements(time)-1 do begin

	idx0 = where(counts_new eq 0.d0)
	if (idx0[0] ne -1) then counts_new[idx0] = counts_new[idx0-1] else break

      endfor


    ;find out where are the first, e.g. 10 min (binsize10) of a night
    ;starting computation from that point
      diff = time-time[0]
      index10 = closest(binsize10/3600.d0, diff, value=value, /upper)	;starting index
      nbins = fix(((time[n_elements(time)-1]-time[index10[0]])/binsize)*3600.d0)

      counts = counts_new
      hr = time

    ;array definition
      ;time bin belongs to 'binsize', i.e. small bin
      hrbin = dblarr(nbins) & hrdum = hrbin & tmp_hr = hr[index10:n_elements(hr)-1] & countsbin10_rms = hrbin & countsbin10_mean = hrbin



    endif else begin
    ;find out where are the first, e.g. 10 min (binsize10) of a night
    ;starting computation from that point
      diff = hr-hr[0]
      index10 = closest(binsize10/3600.d0, diff, value=value, /upper)	;starting index
      nbins = fix(((hr[n_elements(hr)-1]-hr[index10[0]])/binsize)*3600.d0)

    ;array definition
      ;time bin belongs to 'binsize', i.e. small bin
      hrbin = dblarr(nbins) & hrdum = hrbin &  tmp_hr = hr[index10:n_elements(hr)-1] & countsbin10_rms = hrbin & countsbin10_mean = hrbin

    endelse

  ;first element needs to be handled separately
    idx_1st = closest(diff[index10]-binsize/3600.d0, diff[0:index10], /decide)
    hrbin[0] = mean(hr[idx_1st:index10])
    countsbin10_mean[0] = mean(counts[0:index10])
    countsbin10_rms[0] = stddev(counts[0:index10])

    for i=1L,nbins[0]-1 do begin

    ;small bin
      diff = tmp_hr - tmp_hr[0]
      idx = closest(binsize/3600.d0, diff, /decide)

      if (idx[0] eq 0) then begin

	print, ''
	print, 'Binsize is smaller than sampling. Can be caused by data gaps. Stop.'
	print, ''
	stop

      endif

      hrbin[i] = mean(tmp_hr[0:idx-1])
      hrdum[i] = tmp_hr[idx-1]

      tmp_hr = tmp_hr[idx:n_elements(tmp_hr)-1]

    ;big bin
      idx1 = where(hrdum[i] eq hr)
      tdum = hr[idx1]-(binsize10/3600.d0)
      idx2 = closest(tdum, hr, /decide)

      countsbin10_mean[i] = mean(counts[idx2:idx1])
      countsbin10_rms[i] = stddev(counts[idx2:idx1])

    endfor

    rvarcounts_2 = countsbin10_rms/countsbin10_mean
    hrbin_2 = hrbin

  endif else begin

    check2 = 'no'

  endelse


;merge both datasets
  if (check1 ne 'no' and check2 ne 'no') then begin

    rvarcounts_sun = [rvarcounts_1, rvarcounts_2]
    hrbin_sun = [hrbin_1, hrbin_2]

  endif

  if (check1 eq 'no' and check2 ne 'no') then begin

    rvarcounts_sun = [rvarcounts_2]
    hrbin_sun = [hrbin_2]

  endif

  if (check1 ne 'no' and check2 eq 'no') then begin

    rvarcounts_sun = [rvarcounts_1]
    hrbin_sun = [hrbin_1]

  endif


end