@showsym.pro

;the value of 683.002 for W/m^2 to lux from http://en.wikipedia.org/wiki/Candela

pro plot_lightmeter_results, loc1, loc2, pathll, pathlunfit, pathnsb, pathstat


; Calculate the aspect ratio of display window.
aspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE 
; Now calculate the correct size on a Portrait page.
xsize = 8.0
ysize = xsize * aspectRatio
if ysize gt 10.5 then begin
  ysize = 10.5
  xsize = ysize / aspectRatio
endif
; Calculate the offsets, so the output window is not off the page.
xoffset = (8.5 - xsize) / 2.0
yoffset = (11.0 - ysize) / 2.0

!P.Font=1
!p.thick=4
!x.thick=3
!y.thick=3



;=========================================================================================================================



; ;CONVERSION FACTOR
;   filelunfit = file_search(pathlunfit+'LunarModelFit_linearity.txt', count=nlunfit)
;   if (nlunfit eq 1) then begin
;   ;     print, 'No file or multiple files found. Stop.'
;   ;     stop
;   ;   endif
; 
;     readcol, filelunfit[0], jdlunfit, moonill, confact1, confact2, fitres, flag, format='d,d,d,d,d,a', /silent
;     ;consider only fitresults where convfact was not at the limit
;     idx = where(flag eq '+')
;     jdlunfit = jdlunfit[idx]
;     moonill = moonill[idx]
;     confact1 = confact1[idx]
;     confact2 = confact2[idx]
;     fitres = fitres[idx]
; 
;     jdlunfit = jdlunfit - 2450000.d0
; 
;     idxrms = where(fitres lt 0.05)
;     jdlunfit = jdlunfit[idxrms]
;     moonill = moonill[idxrms]
;     confact1 = confact1[idxrms]
;     confact2 = confact2[idxrms]
;     fitres = fitres[idxrms]
; 
;   ;   med_confact2 = median(confact2)
;   ;   sd_confact2 = stddev(confact2)
; 
;   ;sigma clipping ;reason see A01.02.2010 where cloudy parts are set as clear due to insuffiecient lossam like computation
;     med_confact2 = robust_mean_mod(confact2, 3.0, sigma, num_rej, GoodInd=GoodInd2)
;     sigma_confact2 = stddev(confact2[GoodInd2])
; 
;     jdlunfit = jdlunfit[GoodInd2]
;     moonill = moonill[GoodInd2]
;     confact1 = confact1[GoodInd2]
;     confact2 = confact2[GoodInd2]
;     fitres = fitres[GoodInd2]
; 
;     print, ''
;     print, 'Confact2 [W/m^2/counts]: ', med_confact2, sigma_confact2
; 
;     ;fit trend with linear function
;     sixlin, jdlunfit, confact2, a, siga, b, sigb
; 
;   ;plot confact stuff including the avergae value+error + linear fit showing trend of solar cell degeneration
;     set_plot, 'ps'
;     device, filename=pathlunfit+loc2+'_'+loc1+'_Confact2_Time.ps',/color,XSIZE=20, YSIZE=14, XOffset=xoffset, YOffset=yoffset
; 
;       plot, jdlunfit, confact2, /nodata, /ynozero, title=loc1, $
; 	  xtitle='Time [JD-2450000]', ytitle='Conversion Factor [W/m!U2!N/counts]', charsize=1.6, $
; 	  xminor=2, yminor=2, color=fsc_color('black')
;       oplot, jdlunfit, confact2, psym=sym(1), symsize=0.5, color=fsc_color('blue')
; 
;       ;plot linear fit
;       x = !x.crange
;       y = a[0] + b[0]*x
;       oplot, x, y, linestyle=2, thick=4, color=fsc_color('dark gray')
;       tmp1 = sigfig(med_confact2, 4, /sci)
;       tmp2 = sigfig(sigma_confact2, 3, /sci)
;       legend, ['ConFact = '+strcompress(tmp1,/rem)+' +/- '+strcompress(tmp2,/rem)+' W/m!U2!N/counts'], box=0, margin=0, /left, charsize=1.4
; 
;     device,/close
;     set_plot,'x'
;     ;spawn, 'gv '+pathlunfit+loc2+'_'+loc1+'_Confact2_Time.ps'
; 
;   endif


;=========================================================================================================================



;MEDIAN NIGHT SKY BRIGHTNESS
  filensb = file_search(pathnsb+'MedianNightSkyBrightness_Wm2.txt', count=nnsb)
  if (nnsb eq 1) then begin
  ;     print, 'No file or multiple files found. Stop.'
  ;     stop
  ;   endif

    readcol, filensb[0], jdnsb, avecounts, sdcounts, format='d,d,d', /silent
    jdnsb = jdnsb - 2450000.d0

    tmp_tspan = strtrim(string(jdnsb[n_elements(jdnsb)-1]-jdnsb[0]+1),2)
    pos = strpos(tmp_tspan, '.')
    tspan = strmid(tmp_tspan,0,pos)

    nnights = n_elements(jdnsb)

    ;get range of y-axis, is necessary because of outliers
    up = max(avecounts+sdcounts)
    low = min(avecounts-sdcounts)

    set_plot, 'ps'
    device, filename=pathnsb+loc2+'_'+loc1+'_MedianNightSkyBrightness_Wm2.ps',/color,XSIZE=20, YSIZE=14, XOffset=xoffset, YOffset=yoffset

      plot, jdnsb, avecounts, /nodata, title=loc1, $
	  xtitle='Time [JD - 2450000 days]', ytitle='Median Sky Irradiance Dark Time [W/m!U2!N]', charsize=1.6, $
	  xminor=2, yminor=2, color=fsc_color('black'), /ynozero, yr=[low-0.5d-6, up+0.5d-6], yst=1
      oploterror, jdnsb, avecounts, sdcounts, psym=sym(1), color=fsc_color('blue'), errcolor=fsc_color('blue')
      plots, [min(jdnsb),max(jdnsb)],[mean(avecounts),mean(avecounts)],/data, linestyle=2, color=fsc_color('gray'), thick=4
      legend, ['T!Dspan!N = '+strcompress(tspan,/rem)+' days', 'n!Dnights!N = '+strcompress(nnights,/rem)], box=0, margin=0, /bottom, charsize=1.6

    device,/close
    set_plot,'x'
    ;spawn, 'gv '+pathnsb+loc2+'_'+loc1+'_MedianNightSkyBrightness_Wm2.ps'

  endif

;=========================================================================================================================

;MEDIAN NIGHT SKY BRIGHTNESS - NO CLOUDS
  filensb = file_search(pathnsb+'MedianNightSkyBrightness_Wm2_NoClouds.txt', count=nnsb)
  if (nnsb eq 1) then begin
  ;     print, 'No file or multiple files found. Stop.'
  ;     stop
  ;   endif

    readcol, filensb[0], jdnsb, avecounts, sdcounts, format='d,d,d', /silent
    jdnsb = jdnsb - 2450000.d0

    tmp_tspan = strtrim(string(jdnsb[n_elements(jdnsb)-1]-jdnsb[0]+1),2)
    pos = strpos(tmp_tspan, '.')
    tspan = strmid(tmp_tspan,0,pos)

    nnights = n_elements(jdnsb)

    ;get range of y-axis, is necessary because of outliers
    up = max(avecounts+sdcounts)
    low = min(avecounts-sdcounts)

    set_plot, 'ps'
    device, filename=pathnsb+loc2+'_'+loc1+'_MedianNightSkyBrightness_Wm2_NoClouds.ps',/color,XSIZE=20, YSIZE=14, XOffset=xoffset, YOffset=yoffset

      plot, jdnsb, avecounts, /nodata, title=loc1, $
	  xtitle='Time [JD - 2450000 days]', ytitle='Median Sky Irradiance / Dark Time [W/m!U2!N]', charsize=1.6, $
	  xminor=2, yminor=2, color=fsc_color('black'), /ynozero, yr=[low-0.5d-6, up+0.5d-6], yst=1
      oploterror, jdnsb, avecounts, sdcounts, psym=sym(1), color=fsc_color('blue'), errcolor=fsc_color('blue')
      plots, [min(jdnsb),max(jdnsb)],[mean(avecounts),mean(avecounts)],/data, linestyle=2, color=fsc_color('gray'), thick=4
      legend, ['T!Dspan!N = '+strcompress(tspan,/rem)+' days', 'n!Dnights!N = '+strcompress(nnights,/rem)], box=0, margin=0, /bottom, charsize=1.6

    device,/close
    set_plot,'x'
    ;spawn, 'gv '+pathnsb+loc2+'_'+loc1+'_MedianNightSkyBrightness_Wm2_NoClouds.ps'

  endif

;=========================================================================================================================

;COUNTS STATISTIC
  filestat = file_search(pathstat+'IrradianceStatistic.txt', count=nstat)
  if (nstat eq 1) then begin
;     print, 'No file or multiple files found. Stop.'
;     stop
;   endif
    rows = file_lines(filestat[0])

    ;find number of cols
    openr, lun, filestat, /get_Lun
      line = ""
      readf, lun, line
    close, lun
    free_lun, lun
    cols = N_Elements(StrSplit(line, /RegEx, /Extract))


    openr, lun, filestat[0], /get_lun
      data = dblarr(cols, rows)  
      READF, lun, data
    close, lun
    free_lun, lun


    jdstat = reform(data[0,*]) - 2450000.d0

    sz = size(data)
    tmp_data = dblarr(sz[2],sz[1]-1)
    bindata = dblarr(sz[1]-1)
    for i=0L,sz[1]-2 do begin
      tmp_data[*,i] = reform(data[i+1,*])
      bindata[i] = total(tmp_data[*,i])
    endfor

    ;read in setup file
    filesetup = file_search(pathstat+'Setup_IrradianceStatistic.txt', count=nsetup)
    if (nsetup gt 1 or nsetup eq 0) then begin
      print, 'No file or multiple files found. Stop.'
      stop
    endif

    readcol, filesetup[0], dum, bs, format='a,d', /silent
    readcol, filesetup[0], bin, format='d', skipline=1, /silent

    idxval = where(bindata ne 0.)
    bin = bin[idxval]
    bindata = bindata[idxval]

    winl = dblarr(n_elements(bin)) & winr = winl
    for i=0L,n_elements(bin)-1 do begin
      winl[i] = bin[i]-bs/2.d0
      winr[i] = bin[i]+bs/2.d0
    endfor

    tmp_tspan = strtrim(string(jdstat[n_elements(jdstat)-1]-jdstat[0]+1),2)
    pos = strpos(tmp_tspan, '.')
    tspan = strmid(tmp_tspan,0,pos)

    nnights = n_elements(jdstat)
    ;counts
    set_plot, 'ps'
    device, filename=pathstat+loc2+'_'+loc1+'_IrradianceStatistic.ps',/color,XSIZE=20, YSIZE=14, XOffset=xoffset, YOffset=yoffset

      plot, bin, bindata, /nodata, title=loc1, xr=[-5.4,-3], xst=1, $
	  xtitle='log(Irradiance) [W/m!U2!N]', ytitle='Number of Measurements', charsize=1.6, $
	  xminor=2, yminor=2, color=fsc_color('black'), ytickformat='exponent'
  ;    oplot, bin, bindata, psym=sym(1), color=fsc_color('blue')
      for i=0L,n_elements(bin)-1 do begin
	polyfill, [winl[i],winl[i],winr[i],winr[i]], [0,bindata[i],bindata[i],0], color=fsc_color('blue')
      endfor
      legend, ['T!Dspan!N = '+strcompress(tspan,/rem)+' days', 'n!Dnights!N = '+strcompress(nnights,/rem)], box=0, margin=0, /right, charsize=1.6

    device,/close
    set_plot,'x'
;    spawn, 'gv '+pathstat+loc2+'_'+loc1+'_IrradianceStatistic.ps'

  endif

;=========================================================================================================================



;COUNTS STATISTIC - NO CLOUDS
  filestat = file_search(pathstat+'IrradianceStatistic_NoClouds.txt', count=nstat)
  if (nstat eq 1) then begin
  ;     print, 'No file or multiple files found. Stop.'
  ;     stop
  ;   endif


    openr, lun, filestat, /get_Lun
      line = ""
      readf, lun, line
    close, lun
    free_lun, lun
    cols = N_Elements(StrSplit(line, /RegEx, /Extract))


    rows = file_lines(filestat[0])
    openr, lun, filestat[0], /get_lun
      data = dblarr(cols, rows)  
      READF, lun, data
    close, lun
    free_lun, lun


    jdstat = reform(data[0,*]) - 2450000.d0

    sz = size(data)
    tmp_data = dblarr(sz[2],sz[1]-1)
    bindata = dblarr(sz[1]-1)
    for i=0L,sz[1]-2 do begin
      tmp_data[*,i] = reform(data[i+1,*])
      bindata[i] = total(tmp_data[*,i])
    endfor

    ;read in setup file
    filesetup = file_search(pathstat+'Setup_IrradianceStatistic.txt', count=nsetup)
    if (nsetup gt 1 or nsetup eq 0) then begin
      print, 'No file or multiple files found. Stop.'
      stop
    endif

    readcol, filesetup[0], dum, bs, format='a,d', /silent
    readcol, filesetup[0], bin, format='d', skipline=1, /silent

    idxval = where(bindata ne 0.)
    bin = bin[idxval]
    bindata = bindata[idxval]

    winl = dblarr(n_elements(bin)) & winr = winl
    for i=0L,n_elements(bin)-1 do begin
      winl[i] = bin[i]-bs/2.d0
      winr[i] = bin[i]+bs/2.d0
    endfor

    tmp_tspan = strtrim(string(jdstat[n_elements(jdstat)-1]-jdstat[0]+1),2)
    pos = strpos(tmp_tspan, '.')
    tspan = strmid(tmp_tspan,0,pos)

    nnights = n_elements(jdstat)
    ;counts
    set_plot, 'ps'
    device, filename=pathstat+loc2+'_'+loc1+'_IrradianceStatistic_NoClouds.ps',/color,XSIZE=20, YSIZE=14, XOffset=xoffset, YOffset=yoffset

      plot, bin, bindata, /nodata, title=loc1, xr=[-5.4,-3], xst=1, $
	  xtitle='log(Irradiance) [W/m!U2!N]', ytitle='Number of Measurements', charsize=1.6, $
	  xminor=2, yminor=2, color=fsc_color('black'), ytickformat='exponent'
  ;    oplot, bin, bindata, psym=sym(1), color=fsc_color('blue')
      for i=0L,n_elements(bin)-1 do begin
	polyfill, [winl[i],winl[i],winr[i],winr[i]], [0,bindata[i],bindata[i],0], color=fsc_color('blue')
      endfor
      legend, ['T!Dspan!N = '+strcompress(tspan,/rem)+' days', 'n!Dnights!N = '+strcompress(nnights,/rem)], box=0, margin=0, /right, charsize=1.6

    device,/close
    set_plot,'x'
    ;spawn, 'gv '+pathstat+loc2+'_'+loc1+'_IrradianceStatistic_NoClouds.ps'


  endif

;=========================================================================================================================


;CLOUDY NIGHTS STATISTIC
  filell = file_search(pathll+loc2+'_'+loc1+'*_LOSSAMlike.txt', count=nll)
  cloud = dblarr(nll) & jdll = cloud	;array containing the cloud coverage per night in percent (0-100%)

  for i=0L,nll-1 do begin

    readcol, filell[i], time, relcounts, format='d,d', /silent

    ;extract date (this could be retrieved from 'counts statistic', too)
      pos1 = strpos(filell[i], '_L', /reverse_search)
      dum = strmid(filell[i],0,pos1)
      pos2 = strpos(dum, '_', /reverse_search)
      date = strmid(dum,pos2+1,8)
      yr = double(strmid(date,0,4))
      mn = double(strmid(date,4,2))
      dy = double(strmid(date,6,2))
      jdll[i] = (julday(mn, dy, yr, 12., 0., 0.))-2450000.d0

    idxval = where(relcounts ge 0.02)
    if (idxval[0] eq -1) then begin

      cloud[i] = 0.

    endif else begin

      cloud[i] = 100.d0*double(n_elements(idxval))/double(n_elements(relcounts))

    endelse

  endfor

  set_plot, 'ps'
  device, filename=pathstat+loc2+'_'+loc1+'_Clouds_Time.ps',/color,XSIZE=25, YSIZE=10, XOffset=xoffset, YOffset=yoffset

    plot, jdll, cloud, /nodata, /ynozero, title=loc1, $
	xtitle='Time [JD-2450000]', ytitle='Cloud Appearence per Night  [%]', charsize=1.4, $
	xminor=2, yminor=2, color=fsc_color('black')
    oplot, jdll, cloud, psym=sym(1), symsize=0.5, color=fsc_color('blue')
  device,/close
  set_plot,'x'
  ;spawn, 'gv '+pathstat+'Clouds_Time.ps'


  set_plot, 'ps'
  device, filename=pathstat+loc2+'_'+loc1+'_Clouds_Distribution.ps',/color,XSIZE=20, YSIZE=14, XOffset=xoffset, YOffset=yoffset

    cloudbs = 5.	;binsize

    ;the next 2 commands are necessary to define the axis because one cannot set x-range with histoplot ...
    plothist, cloud, xhist, yhist, bin=cloudbs, /noplot	;to get histogram data
    ;to set the plot and axis
    plot, xhist, yhist, /nodata, xr=[0,100], xst=1, charsize=1.6, $
	title=loc1, xtitle='Cloud Appearence per Night  [%]', ytitle='Number of Nights', xminor=2, yminor=2, $
	color=fsc_color('black')
    ;plot histogram
    histoplot, cloud, /FILL, POLYCOLOR='blue', binsize=cloudbs, datacolorname='black', /oplot

    legend, ['n!Dnights!N = '+strcompress(n_elements(jdll),/rem)], box=0, margin=0, /right, charsize=1.6

  device,/close
  set_plot,'x'
  ;spawn, 'gv '+pathstat+loc2+'_'+loc1+'_Clouds_Distribution.ps'


!p.font=0
!p.thick=1
!x.thick=1
!y.thick=1


end