function fit_day_func, x, p

  common data_model, irrad_model_sel, temp_sel

  fit = p[2]*(p[1]*(p[0]*exp((x*(1.d0+p[3]*temp_sel))/p[0])-1.d0)+(x*(1.d0+p[3]*temp_sel)))

TRY TO RETURN FIT istead of this, number of iterations, no limits, output of fits ith RMS
  return, total(abs(fit-irrad_model_sel)/irrad_model_sel)


  ;return, total(fit-irrad_model_sel)^2./irrad_model_sel
  ;return, total(abs(fit-irrad_model_sel))

end


pro fit_global_irradiance_model, loc1, time_sun_orig, counts_sun_orig, temp_sun_orig, sun_irr_tot_orig, $
      time_moon_orig, counts_moon_orig, temp_moon_orig, moon_irr_tot_orig, check, $; time_twi_orig, counts_twi_orig, temp_twi_orig, h_sun_twi_orig, twi_irr_tot_orig, $
      start_params_sun, start_params_sunmoon, $
      pi, fit_params, status, rmsfit, bestnorm

  common data_model, irrad_model_sel, temp_sel


;read in fit parameter limits
  readcol, 'fix/fit_parameter_limits.txt', dum1, dum2, parlim, format='a,a,d', skipline=2, /silent


  if (check eq 'yes') then begin

  ;in order to speed up fit, select onlt evry, e.g., 10th data point depending on the original number of elements
    xx = 10

    ;MOON
    if (n_elements(counts_moon_orig) gt 10000) then begin

      counts_moon_tmp = dblarr(ceil(n_elements(counts_moon_orig)/double(xx)))
      time_moon_tmp = counts_moon_tmp
      temp_moon_tmp = counts_moon_tmp
      moon_irr_tot_tmp = counts_moon_tmp

      for i=0L,n_elements(counts_moon_orig)-1, xx do begin

	counts_moon_tmp[i/xx] = counts_moon_orig[i]
	time_moon_tmp[i/xx] = time_moon_orig[i]
	temp_moon_tmp[i/xx] = temp_moon_orig[i]
	moon_irr_tot_tmp[i/xx] = moon_irr_tot_orig[i]

      endfor

      counts_moon_sel = counts_moon_tmp
      time_moon_sel = time_moon_tmp
      temp_moon_sel = temp_moon_tmp
      moon_irr_tot_sel = moon_irr_tot_tmp

    endif else begin

      counts_moon_sel = counts_moon_orig
      time_moon_sel = time_moon_orig
      temp_moon_sel = temp_moon_orig
      moon_irr_tot_sel = moon_irr_tot_orig

    endelse


    ;SUN
    if (n_elements(counts_sun_orig) gt 10000) then begin

      counts_sun_tmp = dblarr(ceil(n_elements(counts_sun_orig)/double(xx)))
      time_sun_tmp = counts_sun_tmp
      temp_sun_tmp = counts_sun_tmp
      sun_irr_tot_tmp = counts_sun_tmp

      for i=0L,n_elements(counts_sun_orig)-1, xx do begin

	counts_sun_tmp[i/xx] = counts_sun_orig[i]
	time_sun_tmp[i/xx] = time_sun_orig[i]
	temp_sun_tmp[i/xx] = temp_sun_orig[i]
	sun_irr_tot_tmp[i/xx] = sun_irr_tot_orig[i]

      endfor

      counts_sun_sel = counts_sun_tmp
      time_sun_sel = time_sun_tmp
      temp_sun_sel = temp_sun_tmp
      sun_irr_tot_sel = sun_irr_tot_tmp

    endif else begin

      counts_sun_sel = counts_sun_orig
      time_sun_sel = time_sun_orig
      temp_sun_sel = temp_sun_orig
      sun_irr_tot_sel = sun_irr_tot_orig

    endelse


    counts_sel = [counts_sun_sel, counts_moon_sel]
    time_sel = [time_sun_sel, time_moon_sel]
    temp_sel = [temp_sun_sel, temp_moon_sel]
    irrad_model_sel = [sun_irr_tot_sel, moon_irr_tot_sel]


    ;make sure time is increasing
    idxsort = bsort(time_sel)
    time_sel = time_sel[idxsort]
    counts_sel = counts_sel[idxsort]
    temp_sel = temp_sel[idxsort]
    irrad_model_sel = irrad_model_sel[idxsort]


  ;do fitting
    dummy_err = dblarr(n_elements(counts_sel))
    dummy_err[*] = 1.d0

    start_val = start_params_sunmoon

    pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
    pi[0].limited(0) = parlim[0]
    pi[0].limits(0) = parlim[1]
    pi[0].limited(1) = parlim[2]
    pi[0].limits(1) = parlim[3]
    pi[1].limited(0) = parlim[4]
    pi[1].limits(0) = parlim[5]
    pi[1].limited(1) = parlim[6]
    pi[1].limits(1) = parlim[7]
    pi[2].limited(0) = parlim[8]
    pi[2].limits(0) = parlim[9]
    pi[2].limited(1) = parlim[10]
    pi[2].limits(1) = parlim[11]
    pi[3].limited(0) = parlim[12]
    pi[3].limits(0) = parlim[13]
    pi[3].limited(1) = parlim[14]
    pi[3].limits(1) = parlim[15]


    fit_params = MPFITfun('fit_day_func', counts_sel, irrad_model_sel, dummy_err, start_val, weights = dummy_err, $
	    maxiter=50000, parinfo=pi, niter=niter, status=status, bestnorm=bestnorm, /quiet);, ftol=1.d-15, gtol=1.d-15, xtol=1d-15)

    fit_result = fit_params[2]*(fit_params[1]*(fit_params[0]*exp((counts_sel*(1.d0+fit_params[3]*temp_sel))/fit_params[0])-1.d0)+(counts_sel*(1.d0+fit_params[3]*temp_sel)))

    rel_res = (fit_result-irrad_model_sel)/irrad_model_sel	;relative residuals

    rmsfit = sqrt(total((fit_result - irrad_model_sel)^2.d0,/double)/double(n_elements(irrad_model_sel)))

; ;screen output
;     window, 0
;     !p.multi=[0,1,2]
; 
;     plot, time_sel, fit_result, psym=3, /ylog;, xr=[0,24], xst=1
;       oplot, time_sel, irrad_model_sel, psym=3, color=fsc_color('red')
;     plot, time_sel, rel_res, psym=3, xr=[!x.crange(0),!x.crange(1)], xst=1, yr=[-0.3,0.3], yst=1, /nodata
;       oplot, time_sel, rel_res, psym=3, color=fsc_color('blue')
;       plots, [!x.crange(0),!x.crange(1)], [0.0,0.0], linestyle=2, color=fsc_color('gray')
;       plots, [!x.crange(0),!x.crange(1)], [0.1,0.1], linestyle=1
;       plots, [!x.crange(0),!x.crange(1)], [-0.1,-0.1], linestyle=1
;       plots, [!x.crange(0),!x.crange(1)], [0.2,0.2], linestyle=1
;       plots, [!x.crange(0),!x.crange(1)], [-0.2,-0.2], linestyle=1
; 
;     !p.multi=[0,1,0]
; 
; ;graphical output
;     ; Calculate the aspect ratio of display window.
;     aspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE 
;     ; Now calculate the correct size on a Portrait page.
;     xsize = 8.0
;     ysize = xsize * aspectRatio
;     if ysize gt 10.5 then begin
;       ysize = 10.5
;       xsize = ysize / aspectRatio
;     endif
;     ; Calculate the offsets, so the output window is not off the page.
;     xoffset = (8.5 - xsize) / 2.0
;     yoffset = (11.0 - ysize) / 2.0
; 
;     !P.Font=1
;     !p.thick=4
;     !x.thick=3
;     !y.thick=3
; 
;     set_plot, 'ps'
;     device, filename=loc1+'_20100329.ps',/color,XSIZE=20, YSIZE=16, XOffset=xoffset, YOffset=yoffset
; 
;     !p.multi=[0,1,2]
; 
;       plot, time_sel, fit_result, xr=[0,24], xst=1, yr=[1.d-4, 1.d3], yst=1, ytickformat='logticks_exp', title=loc1+' 29.03.2010', $
; 	    /nodata, /ylog, /ynozero, ytitle='Irradiance [W/m!U2!N]', charsize=1.4, pos=[0.16,0.41,0.97,0.94], $
; 	    background=fsc_color('white'), Color=fsc_color('black'), xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '], xtickinterval=3, xminor=3
;       oplot, time_sel, fit_result, color=fsc_color('blue'), thick=5, psym=3
;       oplot, time_sel, irrad_model_sel, color=fsc_color('red'), thick=2, psym=3
; 
; 
;     plot, time_sel, rel_res, psym=3, xr=[0,24], xst=1, yr=[-0.2,0.2], yst=1, /nodata, charsize=1.4, yminor=2, $
; 	  background=fsc_color('white'), Color=fsc_color('black'), xtitle='Time UT [hr]', ytitle='Rel. RMS', $
; 	  xtickname=['12','15','18','21','00','03','06','09','12'], xtickinterval=3, xminor=3, pos=[0.16,0.1,0.97,0.4]
;       oplot, time_sel, rel_res, psym=3, color=fsc_color('blue')
;       plots, [!x.crange(0),!x.crange(1)], [0.0,0.0], linestyle=2, color=fsc_color('dark gray')
;       plots, [!x.crange(0),!x.crange(1)], [0.1,0.1], linestyle=1, color=fsc_color('gray')
;       plots, [!x.crange(0),!x.crange(1)], [-0.1,-0.1], linestyle=1, color=fsc_color('gray')
; ;       plots, [!x.crange(0),!x.crange(1)], [0.2,0.2], linestyle=1
; ;       plots, [!x.crange(0),!x.crange(1)], [-0.2,-0.2], linestyle=1
; 
; 
;   device,/close
;   set_plot,'x'
;   spawn, 'gv '+loc1+'_20100329.ps'
; 
; 
; stop

  endif else begin	;only fit sun


  ;in order to speed up fit, select onlt evry, e.g., 10th data point depending on the original number of elements
    xx = 10

    ;SUN
    if (n_elements(counts_sun_orig) gt 10000) then begin

      counts_sun_tmp = dblarr(ceil(n_elements(counts_sun_orig)/double(xx)))
      time_sun_tmp = counts_sun_tmp
      temp_sun_tmp = counts_sun_tmp
      sun_irr_tot_tmp = counts_sun_tmp

      for i=0L,n_elements(counts_sun_orig)-1, xx do begin

	counts_sun_tmp[i/xx] = counts_sun_orig[i]
	time_sun_tmp[i/xx] = time_sun_orig[i]
	temp_sun_tmp[i/xx] = temp_sun_orig[i]
	sun_irr_tot_tmp[i/xx] = sun_irr_tot_orig[i]

      endfor

      counts_sun_sel = counts_sun_tmp
      time_sun_sel = time_sun_tmp
      temp_sun_sel = temp_sun_tmp
      sun_irr_tot_sel = sun_irr_tot_tmp

    endif else begin

      counts_sun_sel = counts_sun_orig
      time_sun_sel = time_sun_orig
      temp_sun_sel = temp_sun_orig
      sun_irr_tot_sel = sun_irr_tot_orig

    endelse


    counts_sel = [counts_sun_sel]
    time_sel = [time_sun_sel]
    temp_sel = [temp_sun_sel]
    irrad_model_sel = [sun_irr_tot_sel]


;     counts_sel = [counts_sun_sel, counts_twi_orig]
;     time_sel = [time_sun_sel, time_twi_orig]
;     temp_sel = [temp_sun_sel, temp_twi_orig]
;     irrad_model_sel = [sun_irr_tot_sel, twi_irr_tot_orig]

    ;make sure time is increasing
    idxsort = bsort(time_sel)
    time_sel = time_sel[idxsort]
    counts_sel = counts_sel[idxsort]
    temp_sel = temp_sel[idxsort]
    irrad_model_sel = irrad_model_sel[idxsort]


  ;do fitting
    dummy_err = dblarr(n_elements(counts_sel))
    dummy_err[*] = 1.d0

    start_val = start_params_sun

    pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
    pi[0].limited(0) = parlim[0]
    pi[0].limits(0) = parlim[1]
    pi[0].limited(1) = parlim[2]
    pi[0].limits(1) = parlim[3]
    pi[1].limited(0) = parlim[4]
    pi[1].limits(0) = parlim[5]
    pi[1].limited(1) = parlim[6]
    pi[1].limits(1) = parlim[7]
    pi[2].limited(0) = parlim[8]
    pi[2].limits(0) = parlim[9]
    pi[2].limited(1) = parlim[10]
    pi[2].limits(1) = parlim[11]
    pi[3].limited(0) = parlim[12]
    pi[3].limits(0) = parlim[13]
    pi[3].limited(1) = parlim[14]
    pi[3].limits(1) = parlim[15]


    fit_params = MPFITfun('fit_day_func', counts_sel, irrad_model_sel, dummy_err, start_val, weights = dummy_err, $
	    maxiter=5000, parinfo=pi, niter=niter, status=status, bestnorm=bestnorm, /quiet);, ftol=1.d-15, gtol=1.d-15, xtol=1d-15)

    fit_result = fit_params[2]*(fit_params[1]*(fit_params[0]*exp((counts_sel*(1.d0+fit_params[3]*temp_sel))/fit_params[0])-1.d0)+(counts_sel*(1.d0+fit_params[3]*temp_sel)))

    rel_res = (fit_result-irrad_model_sel)/irrad_model_sel	;relative residuals

    rmsfit = sqrt(total((fit_result - irrad_model_sel)^2.d0,/double)/double(n_elements(irrad_model_sel)))


;     window, 1
;     !p.multi=[0,1,2]
; 
;     plot, time_sel, fit_result, psym=3, /ylog;, xr=[0,24], xst=1
;       oplot, time_sel, irrad_model_sel, psym=3, color=fsc_color('red')
;     plot, time_sel, rel_res, psym=3, xr=[!x.crange(0),!x.crange(1)], xst=1, yr=[-0.3,0.3], yst=1, /nodata
;       oplot, time_sel, rel_res, psym=3, color=fsc_color('blue')
;       plots, [!x.crange(0),!x.crange(1)], [0.0,0.0], linestyle=2, color=fsc_color('gray')
;       plots, [!x.crange(0),!x.crange(1)], [0.1,0.1], linestyle=1
;       plots, [!x.crange(0),!x.crange(1)], [-0.1,-0.1], linestyle=1
;       plots, [!x.crange(0),!x.crange(1)], [0.2,0.2], linestyle=1
;       plots, [!x.crange(0),!x.crange(1)], [-0.2,-0.2], linestyle=1
; 
;     !p.multi=[0,1,0]

  endelse


end

