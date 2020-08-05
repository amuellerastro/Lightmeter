pro irradiance_statistic, loc1, loc2, pathstat, timesort, counts, jd_day, t18rise, t18set, binsize, hrbin, rvarcounts

;counts is irradiance


;ALL COUNTS, i.e. WITH CLOUDS

;only consider counts during dark time
  idx = where(timesort ge t18set and timesort le t18rise)

  counts_dark = counts[idx]

;sort counts
  counts_sort = counts_dark[sort(counts_dark)]

;LOGARITHMIC COMPUTATION
  ;assumption: irradiance ranging from 1.d-5 to 1d-2
  bs = 0.01d0	;log bin size
  min_count = alog10(1.0d-6)
  max_count = alog10(6.0d-4)
  oom = abs(max_count-min_count)	;orders of magnitude

; 1.3d-5  -- -4.88
; 5.2d-4 -- -3.29

win = min_count+dindgen(oom/bs+1)*bs

totnum = dblarr(n_elements(win)-1)	;number of counts per window, e.g. bs=5000 -> win[0]=0 to 4999 counts but for plotting its 2500
bin = totnum


for i=0L,n_elements(win)-2 do begin
  bin[i] = win[i]+(win[i+1]-win[i])/2.d0
endfor
;last element
  bin[n_elements(bin)-1] = win[n_elements(win)-2]+(win[n_elements(win)-1]-win[n_elements(win)-2])/2.d0

counts_dark = alog10(counts_sort)

for i=0L,n_elements(win)-2 do begin

  idx = where(counts_dark ge win[i] and counts_dark lt win[i+1])

  if (idx[0] eq -1) then begin

    totnum[i] = 0.

  endif else begin

    totnum[i] = n_elements(idx)

  endelse

endfor

;write out
  openw, lun, pathstat+'IrradianceStatistic.txt' , width=5000, /get_lun, /append
    printf, lun, jd_day, round(totnum)
  close, lun
  free_lun, lun

;write out properties, done every time but okay...
  openw, lun, pathstat+'Setup_IrradianceStatistic.txt', /get_lun

    printf, lun, 'Binsize: ', bs
;     printf, lun, 'Maximum value: ', round(max_count)
;     printf, lun, ''

    for i=0,n_elements(bin)-1 do begin
      printf, lun, bin[i]
    endfor

  close, lun
  free_lun, lun



;SELECTED COUNTS, i.e. NO CLOUDS

;only consider counts during dark time
  idx = where(timesort ge t18set and timesort le t18rise)
  counts_dark = counts[idx]
  time_dark = timesort[idx]

;cut out areas having clouds
  idx_clean = where(rvarcounts lt 0.02d0)
  if (idx_clean[0] ne -1) then begin

    hrbin_clean = hrbin[idx_clean]
    rvarcounts_clean = rvarcounts[idx_clean]

  endif else begin

    return	;i.e. if whole night is cloudy do nothing
;     hrbin_clean = hrbin
;     rvarcounts_clean = rvarcounts

  endelse

;create time array in order to match time_dark with the non cloudy parts of the night
  tleft = hrbin_clean - binsize[0]/(2.d0*3600.d0)
  tright = hrbin_clean + binsize[0]/(2.d0*3600.d0)

counts_tmp = dblarr(n_elements(counts_dark)) & time_tmp = counts_tmp

for i=0L,n_elements(tleft)-1 do begin

  idx = where(time_dark ge tleft[i] and time_dark le tright[i])

  if (idx[0] ne -1) then begin

    idx0 = where(counts_tmp eq 0.)
    counts_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = counts_dark[idx]
    time_tmp[idx0[0]:idx0[0]+n_elements(idx)-1] = time_dark[idx]

  endif

endfor

idxval = where(counts_tmp ne 0.)
counts_clean = counts_tmp[idxval]
time_clean = time_tmp[idxval]

totnum = dblarr(n_elements(win)-1)
counts_clean = alog10(counts_clean)
for i=0L,n_elements(win)-2 do begin

  idx = where(counts_clean ge win[i] and counts_clean lt win[i+1])

  if (idx[0] eq -1) then begin

    totnum[i] = 0.

  endif else begin

    totnum[i] = n_elements(idx)

  endelse

endfor

;write out
  openw, lun, pathstat+'IrradianceStatistic_NoClouds.txt' , width=5000, /get_lun, /append
    printf, lun, jd_day, round(totnum)
  close, lun
  free_lun, lun


end