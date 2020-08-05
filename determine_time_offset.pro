function fit_sun_func, x, p

  common data_sun, sun_irr_tot_sel, temp_sun_sel

  fit = p[2]*(p[1]*(p[0]*exp((x*(1.d0+p[3]*temp_sun_sel))/p[0])-1.d0)+(x*(1.d0+p[3]*temp_sun_sel)))

  return, total(abs(fit-sun_irr_tot_sel)/sun_irr_tot_sel)

  ;return, total((fit-sun_irr_tot_sel)^2./sun_irr_tot_sel)
  ;return, total((fit-sun_irr_tot)/sun_irr_tot)

end


pro determine_time_offset, time_sun, counts_sun, temp_sun, h_sun, sun_irr_tot, start_params_sun, t_offset

  common data_sun, sun_irr_tot_sel, temp_sun_sel


;only use every 10th data point to speed up
;   pointsel = 1	;only use every XX point
; 
;   counts_sun_tmp = dblarr(ceil(n_elements(counts_sun)/double(pointsel)))
;   time_sun_tmp = counts_sun_tmp
;   temp_sun_tmp = counts_sun_tmp
;   h_sun_tmp = counts_sun_tmp
;   sun_irr_tot_tmp = counts_sun_tmp
; 
;   for i=0L,n_elements(counts_sun)-1, pointsel do begin
; 
;     counts_sun_tmp[i/pointsel] = counts_sun[i]
;     time_sun_tmp[i/pointsel] = time_sun[i]
;     temp_sun_tmp[i/pointsel] = temp_sun[i]
;     h_sun_tmp[i/pointsel] = h_sun[i]
;     sun_irr_tot_tmp[i/pointsel] = sun_irr_tot[i]
; 
;   endfor
; 
;   counts_sun_sel = counts_sun_tmp
;   time_sun_sel = time_sun_tmp
;   temp_sun_sel = temp_sun_tmp
;   h_sun_sel = h_sun_tmp
;   sun_irr_tot_sel = sun_irr_tot_tmp
  counts_sun_sel = counts_sun
  time_sun_sel = time_sun
  temp_sun_sel = temp_sun
  h_sun_sel = h_sun
  sun_irr_tot_sel = sun_irr_tot


;do fitting

  counts_sun_sel_smooth = smooth(counts_sun_sel, 11, /edge_truncate)

  dummy_err = dblarr(n_elements(counts_sun_sel))
  dummy_err[*] = 1.d0

  readcol, 'fix/fit_parameter_limits.txt', dum1, dum2, parlim, format='a,a,d', skipline=2, /silent

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


  params = MPFITfun('fit_sun_func', counts_sun_sel_smooth, sun_irr_tot_sel, dummy_err, start_val, weights = dummy_err, $
	  maxiter=500, parinfo=pi, niter=niter, status=status, bestnorm=bestnorm, /quiet)

  fit_result = params[2]*(params[1]*(params[0]*exp((counts_sun_sel*(1.d0+params[3]*temp_sun_sel))/params[0])-1.d0)+(counts_sun_sel*(1.d0+params[3]*temp_sun_sel)))

  rel_res = (fit_result-sun_irr_tot_sel)/sun_irr_tot_sel	;relative residuals

;   window, 0
;   !p.multi=[0,1,2]
; 
;   plot, time_sun_sel, fit_result, psym=3
;     oplot, time_sun_sel, sun_irr_tot_sel, psym=3, color=fsc_color('red')
;   plot, time_sun_sel, rel_res, psym=3, xr=[!x.crange(0),!x.crange(1)], xst=1, yr=[-0.3,0.3], yst=1, /nodata
;     oplot, time_sun_sel, rel_res, psym=3, color=rgb(0,0,255)
;     plots, [!x.crange(0),!x.crange(1)], [0.0,0.0], linestyle=2, color=fsc_color('gray')
;     plots, [!x.crange(0),!x.crange(1)], [0.1,0.1], linestyle=1
;     plots, [!x.crange(0),!x.crange(1)], [-0.1,-0.1], linestyle=1
;     plots, [!x.crange(0),!x.crange(1)], [0.2,0.2], linestyle=1
;     plots, [!x.crange(0),!x.crange(1)], [-0.2,-0.2], linestyle=1
; 
;   !p.multi = [0,1,0]


;check for time shifts, sampling rate of 1 second assumed!

  win = 1200	;allow +/-20min shift
  offset = dindgen(2*win+1)-win
  rms = dblarr(n_elements(offset))

  for i=0L,n_elements(offset)-1 do begin

    time_shift = time_sun_sel+offset[i]/3600.d0

    ;adjust arrays to compensate for time shift
    if (time_sun_sel[0] ge time_shift[0]) then begin

      idx1 = closest(time_sun_sel[0],time_shift)
      idx2 = closest(time_shift[n_elements(fit_result)-1],time_sun_sel)

      fit_result_cut = fit_result[idx1[0]:n_elements(fit_result)-1]
      model_cut = sun_irr_tot_sel[0:idx2[0]]

      rms[i] = sqrt(total((fit_result_cut - model_cut)^2.d0,/double)/double(n_elements(model_cut)))

    endif else begin

      idx1 = closest(time_shift[0],time_sun_sel)
      idx2 = closest(time_sun_sel[n_elements(fit_result)-1],time_shift)

      fit_result_cut = fit_result[0:idx2[0]]
      model_cut = sun_irr_tot_sel[idx1[0]:n_elements(fit_result)-1]

      rms[i] = sqrt(total((fit_result_cut - model_cut)^2.d0,/double)/double(n_elements(model_cut)))

    endelse

  endfor

  dum = min(rms, idxmin)

  t_offset = (offset[idxmin])/86400.d0

  if (t_offset gt 7.d-4) then begin

    t_offset = (offset[idxmin])/86400.d0	;days
    print, '   time offset: '+sigfig(t_offset*1440.d0,4)+' minutes'

  endif else begin

    ;if best fit value is at the edges than stop
    if (idxmin eq n_elements(offset)-1 or idxmin eq 0) then begin

      print, '   Limits of considered time offset (20 min) reached! Stop.'
      stop

    endif else begin

      t_offset = 0.
      print, '   time offset: '+sigfig(t_offset*1440.d0,4)+' minutes'

    endelse

  endelse


end
