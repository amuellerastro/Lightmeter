;assuming that DIT is even number, i.e. no 1.5sec but 1 or 2 sec

pro lightmeter_lossamlike_moon, path, loc1, loc2, hr_orig, file_day, file_month, file_year, counts_tmp, t18set, t18rise, binsize, hrbin, rvarcounts ;, check, hr_moon, counts_model


; binsize = 60.d0	;60 seconds
; binsize10 = 120.d0	;120

;read in setup file
  readcol, 'fix/SETUP.txt', dum, binsize, binsize10, format='a,d,d', skipline=4, /silent

binsize = binsize[0]
binsize10 = binsize10[0]

;shift time by 12hr to avoid the 24->0 jump
timeplot = hr_orig
idx1 = where(timeplot ge 12.)
idx2 = where(timeplot lt 12.)
if (idx1[0] ne -1) then timeplot[idx1] = timeplot[idx1]-12.
if (idx2[0] ne -1) then timeplot[idx2] = timeplot[idx2]+12.
hr = timeplot

; if (check eq 'yes') then begin
;   timeplot_moon = hr_moon
;   idx1 = where(timeplot_moon ge 12.)
;   idx2 = where(timeplot_moon lt 12.)
;   if (idx1[0] ne -1) then timeplo_moont[idx1] = timeplot_moon[idx1]-12.
;   if (idx2[0] ne -1) then timeplot_moon[idx2] = timeplot_moon[idx2]+12.
;   hr_moon = timeplot_moon
; endif


;considering only night time, i.e. sun -18 deg below horizon
  counts = counts_tmp
  idx = where(hr le t18rise and hr ge t18set)
  hr = hr[idx]
  counts = counts[idx]

;========================================================================
;substract model counts from measured counts caused by lunar illumination
; if (check eq 'yes') then begin
; 
;   idx1 = where(hr_moon[0] eq hr)
;   idx2 = where(hr_moon[n_elements(hr_moon)-1] eq hr)
; 
;   counts[idx1:idx2] = counts[idx1:idx2]-counts_model
; 
;   counts = counts+abs(min(counts))
; 
; endif
;========================================================================


;check for time gaps
;   t_tmp = hr;[index10:n_elements(hr)-1]
;   c_tmp = counts;[index10:n_elements(hr)-1]


  diff_all = abs(ts_diff(hr, 1))*3600.d0	;difference to the next time stamp
  ave_samp = median(diff_all)	;mean sampling rate
  ave_samp = round(median(diff_all[where(diff_all lt 1.5*ave_samp)]))	;just to make it more reliable

;=============================================================================================================
;NEEDS to be optimized concerning speed

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
    ;fix value can be a problem for huge nbins
    nbins = fix(((time[n_elements(time)-1]-time[index10[0]])/binsize)*3600.d0)


    counts = counts_new
    hr = time

  ;array definition
    ;time bin belongs to 'binsize', i.e. small bin
    hrbin = dblarr(nbins) & hrdum = hrbin & tmp_hr = hr[index10:n_elements(hr)-1] & countsbin10_rms = hrbin & countsbin10_mean = hrbin

;   ;if there are gaps in time fill them
;   idx = where(diff_all gt 1.5*ave_samp)
;   if (idx[0] ne -1) then begin
; 
;     nts = (t_tmp[n_elements(t_tmp)-1] - t_tmp[0])*3600.d0*ave_samp
;     time = dblarr(nts+2)	;new time array, number of elements is the same if there would be no time gaps
;     counts_new = time
; 
;     ;first element
;     time[0:idx[0]] = t_tmp[0:idx[0]]
;     counts_new[0:idx[0]] = c_tmp[0:idx[0]]
;     tspan = round((t_tmp[idx[0]+1] - t_tmp[idx[0]])*3600.d0)
;     arr_tspan = dindgen(tspan-1)+1
;     time[idx[0]+1:idx[0]+tspan-1] = arr_tspan/ave_samp/3600.d0 + time[idx[0]]
;     counts_new[idx[0]+1:idx[0]+tspan-1] = c_tmp[idx[0]]
; 
; 
;     for i=1L,n_elements(idx)-1 do begin
; 
;       ;original data
;       idx0 = where(time eq 0.)
;       time[idx0[0]:idx0[0]+(idx[i]-idx[i-1])-1] = t_tmp[idx[i-1]+1:idx[i]]
;       counts_new[idx0[0]:idx0[0]+(idx[i]-idx[i-1])-1] = c_tmp[idx[i-1]+1:idx[i]]
; 
;       ;dummy data
;       tspan = round((t_tmp[idx[i]+1] - t_tmp[idx[i]])*3600.d0)
;       arr_tspan = dindgen(tspan-1)+1
;       idx0 = where(time eq 0.)
;       time[idx0[0]:idx0[0]+tspan-2] = arr_tspan/ave_samp/3600.d0 + time[idx0[0]-1]
;       counts_new[idx0[0]:idx0[0]+tspan-2] = c_tmp[idx[i]-1]	;fill the gaps with the last measured value
; 
;     endfor
; 
;     ;last element
;     idx0 = where(time eq 0.)
;     time[idx0[0]:n_elements(time)-1] = t_tmp[idx[n_elements(idx)-1]+1:n_elements(t_tmp)-1]
;     counts_new[idx0[0]:n_elements(counts_new)-1] = c_tmp[idx[n_elements(idx)-1]+1:n_elements(c_tmp)-1]
; 
;   ;find out where are the first, e.g. 10 min (binsize10) of a night
;   ;starting computation from that point
;     diff = time-time[0]
;     index10 = closest(binsize10/3600.d0, diff, value=value, /upper)	;starting index
;     nbins = fix(((time[n_elements(time)-1]-time[index10[0]])/binsize)*3600.d0)
; 
;     counts = counts_new
;     hr = time
; 
;   ;array definition
;     ;time bin belongs to 'binsize', i.e. small bin
;     hrbin = dblarr(nbins) & countsbin_rms = hrbin & countsbin10_mean = hrbin & hrdum = hrbin
;     tmp_hr = hr[index10:n_elements(hr)-1] & tmp_counts = counts_new[index10:n_elements(counts_new)-1] 
;     tmp_hr10 = tmp_hr & tmp_counts10 = counts_new[index10:n_elements(counts_new)-1]

;=============================================================================================================

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

rvarcounts = countsbin10_rms/countsbin10_mean


;OUTPUT
;write PS file

  ;Calculate the aspect ratio of display window.
  aspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE 
  ;Now calculate the correct size on a Portrait page.
  xsize = 8.0
  ysize = xsize * aspectRatio
  IF ysize GT 10.5 THEN BEGIN
    ysize = 10.5
    xsize = ysize / aspectRatio
  ENDIF
  ;Calculate the offsets, so the output window is not off the page.
  xoffset = (8.5 - xsize) / 2.0
  yoffset = (11.0 - ysize) / 2.0

  set_plot, 'ps'
  device, filename=path+loc2+'_'+loc1+'_'+strcompress(file_year,/rem)+strcompress(file_month,/rem)+strcompress(file_day,/rem)+'_LOSSAMlike.ps',/color,XSIZE=20, YSIZE=11.6, XOffset=xoffset, YOffset=yoffset

    !P.Font=1
    !p.thick=4
    !x.thick=3
    !y.thick=3

  plot, hrbin, rvarcounts,$
	yr=[0,0.1],yst=1,xtitle='Time UT [hr]',ytitle='Rel. Flux RMS',$
	background=fsc_color('white'),Color=fsc_color('black'),/nodata,$
	title=loc1+' '+strcompress(file_day,/rem)+'.'+strcompress(file_month,/rem)+'.'+strcompress(file_year,/rem),$
	xr=[0,24],xst=1,xtickinterval=3,xminor=3, yminor=2, charsize=1.3, xtickname=['12','15','18','21','00','03','06','09','12']
  oplot, hrbin, rvarcounts,color=fsc_color('blue')
  plots, [t18set,t18set],[0,0.1],/data, linestyle=2, color=fsc_color('orange red')
  plots, [t18rise,t18rise],[0,0.1],/data, linestyle=2, color=fsc_color('orange red')
  plots, [t18set,t18rise],[0.02,0.02], linestyle=2,/data,color=fsc_color('gray')
  plots, [t18set,t18rise],[0.05,0.05], linestyle=2,/data,color=fsc_color('gray')

device,/close
set_plot,'x'
!p.font=0
!p.thick=1
!x.thick=1
!y.thick=1

;spawn, 'gv '+path+loc2+'_'+loc1+'_'+strcompress(file_year,/rem)+strcompress(file_month,/rem)+strcompress(file_day,/rem)+'_likeLOSSAM.ps'
;stop

;convert hrbin back to 'normal' time scale
hrout = hrbin
idx1 = where(hrout ge 12.)
idx2 = where(hrout lt 12.)
if (idx1[0] ne -1) then hrout[idx1] = hrout[idx1]-12.
if (idx2[0] ne -1) then hrout[idx2] = hrout[idx2]+12.


;write data in an ASCII file
  openw, lun, path+loc2+'_'+loc1+'_'+strcompress(file_year,/rem)+strcompress(file_month,/rem)+strcompress(file_day,/rem)+'_LOSSAMlike.txt', /get_lun
    printf, lun, 'Time UT [hr] RelFlux rms'
    for i=0L,n_elements(hrout)-1 do begin
      printf, lun, hrout[i], rvarcounts[i], format='(f9.5,f15.10)'
    endfor

  close, lun
  free_lun, lun


end