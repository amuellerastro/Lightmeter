pro median_fit_params, pathmodelfit, start_params_sunmoon, a, delta_a, b, delta_b, c, delta_c, d, delta_d
;pro median_fit_params, a, delta_a, b, delta_b, c, delta_c, d, delta_d

;pathmodelfit = '/home/amueller/LightPollution/Lightmeter/Paranal/IrradianceModelFit/'
;pathmodelfit = '/home/amueller/LightPollution/Lightmeter/Armazones/IrradianceModelFit/'


;read in fit results, its a global file and shouldn't changed
  readcol, pathmodelfit+'IrradianceModelFit.txt', jd_day, phaseangle_moon, illfrac_moon, moon_included, p1, p2, p3, p4, rmsfit, ssr, fit_status, t_offset, $
	    format='d,d,d,a,d,d,d,d,d,d,i,d',/silent

;make selection for computation of mean parameters
  idxmoon = where(moon_included eq 'yes')
  ;if (idxmoon[0] ne -1) then begin
  if (n_elements(idxmoon) gt 400) then begin

    jd_day = jd_day[idxmoon]
    phaseangle_moon = phaseangle_moon[idxmoon]
    illfrac_moon = illfrac_moon[idxmoon]
    moon_included = moon_included[idxmoon]
    p1 = p1[idxmoon]
    p2 = p2[idxmoon]
    p3 = p3[idxmoon]
    p4 = p4[idxmoon]
    rmsfit = rmsfit[idxmoon]
    ssr = ssr[idxmoon]
    fit_status = fit_status[idxmoon]
    t_offset = t_offset[idxmoon]

  endif else begin

    ;if no fit values present, than use estimates from the setup file
    a = start_params_sunmoon[0]
    b = start_params_sunmoon[1]
    c = start_params_sunmoon[2]
    d = start_params_sunmoon[3]

    delta_a = 0.
    delta_b = 0.
    delta_c = 0.
    delta_d = 0.

    return
  endelse


  readcol, 'fix/fit_parameter_limits.txt', dum1, dum2, parlim, format='a,a,d', skipline=2, /silent
  lima = [parlim[1], parlim[3]]
  limb = [parlim[5], parlim[7]]
  limc = [parlim[9], parlim[11]]
  limd = [parlim[13], parlim[15]]

  idxa = where(p1 ne lima[0] or p1 ne lima[1])
  if (idxa[0] ne -1) then begin

    jd_day = jd_day[idxa]
    phaseangle_moon = phaseangle_moon[idxa]
    illfrac_moon = illfrac_moon[idxa]
    moon_included = moon_included[idxa]
    p1 = p1[idxa]
    p2 = p2[idxa]
    p3 = p3[idxa]
    p4 = p4[idxa]
    rmsfit = rmsfit[idxa]
    ssr = ssr[idxa]
    fit_status = fit_status[idxa]
    t_offset = t_offset[idxa]

  endif

  idxb = where(p2 ne limb[0] or p2 ne limb[1])
  if (idxb[0] ne -1) then begin

    jd_day = jd_day[idxb]
    phaseangle_moon = phaseangle_moon[idxb]
    illfrac_moon = illfrac_moon[idxb]
    moon_included = moon_included[idxb]
    p1 = p1[idxb]
    p2 = p2[idxb]
    p3 = p3[idxb]
    p4 = p4[idxb]
    rmsfit = rmsfit[idxb]
    ssr = ssr[idxb]
    fit_status = fit_status[idxb]
    t_offset = t_offset[idxb]

  endif

  idxc = where(p3 ne limc[0] or p3 ne limc[1])
  if (idxc[0] ne -1) then begin

    jd_day = jd_day[idxc]
    phaseangle_moon = phaseangle_moon[idxc]
    illfrac_moon = illfrac_moon[idxc]
    moon_included = moon_included[idxc]
    p1 = p1[idxc]
    p2 = p2[idxc]
    p3 = p3[idxc]
    p4 = p4[idxc]
    rmsfit = rmsfit[idxc]
    ssr = ssr[idxc]
    fit_status = fit_status[idxc]
    t_offset = t_offset[idxc]

  endif

  idxd = where(p4 ne limd[0] or p4 ne limd[1])
  if (idxd[0] ne -1) then begin

    jd_day = jd_day[idxd]
    phaseangle_moon = phaseangle_moon[idxd]
    illfrac_moon = illfrac_moon[idxd]
    moon_included = moon_included[idxd]
    p1 = p1[idxd]
    p2 = p2[idxd]
    p3 = p3[idxd]
    p4 = p4[idxd]
    rmsfit = rmsfit[idxd]
    ssr = ssr[idxd]
    fit_status = fit_status[idxd]
    t_offset = t_offset[idxd]

  endif

  idxfitstatus = where(fit_status gt 0)
  if (idxfitstatus[0] ne -1) then begin

    jd_day = jd_day[idxfitstatus]
    phaseangle_moon = phaseangle_moon[idxfitstatus]
    illfrac_moon = illfrac_moon[idxfitstatus]
    moon_included = moon_included[idxfitstatus]
    p1 = p1[idxfitstatus]
    p2 = p2[idxfitstatus]
    p3 = p3[idxfitstatus]
    p4 = p4[idxfitstatus]
    rmsfit = rmsfit[idxfitstatus]
    ssr = ssr[idxfitstatus]
    fit_status = fit_status[idxfitstatus]
    t_offset = t_offset[idxfitstatus]

  endif

  idxrms = where(rmsfit lt median(rmsfit))
  if (idxrms[0] ne -1) then begin

    jd_day = jd_day[idxrms]
    phaseangle_moon = phaseangle_moon[idxrms]
    illfrac_moon = illfrac_moon[idxrms]
    moon_included = moon_included[idxrms]
    p1 = p1[idxrms]
    p2 = p2[idxrms]
    p3 = p3[idxrms]
    p4 = p4[idxrms]
    rmsfit = rmsfit[idxrms]
    ssr = ssr[idxrms]
    fit_status = fit_status[idxrms]
    t_offset = t_offset[idxrms]

  endif


;   idxssr = where(ssr lt median(ssr))
;   if (idxssr[0] ne -1) then begin
; 
;     jd_day = jd_day[idxssr]
;     phaseangle_moon = phaseangle_moon[idxssr]
;     illfrac_moon = illfrac_moon[idxssr]
;     moon_included = moon_included[idxssr]
;     p1 = p1[idxssr]
;     p2 = p2[idxssr]
;     p3 = p3[idxssr]
;     p4 = p4[idxssr]
;     rmsfit = rmsfit[idxssr]
;     ssr = ssr[idxssr]
;     fit_status = fit_status[idxssr]
;     t_offset = t_offset[idxssr]
; 
;   endif

  n = n_elements(p1)
  if (n lt 40) then begin

    ;if no fit values present, than use estimates from the setup file
    a = start_params_sunmoon[0]
    b = start_params_sunmoon[1]
    c = start_params_sunmoon[2]
    d = start_params_sunmoon[3]

    delta_a = 0.
    delta_b = 0.
    delta_c = 0.
    delta_d = 0.

    return

  endif

  a = median(p1)
  delta_a = stddev(p1)
  b = median(p2)
  delta_b = stddev(p2)
  c = median(p3)
  delta_c = stddev(p3)
  d = median(p4)
  delta_d = stddev(p4)

; ;y = c*(b*(a*exp((x*(1.d0+d*T))/a)-1.d0)+(x*(1.d0+d*T)))
; ;partial derivatives
;   dyda = b*c*exp(x*(1.d0+d*T)/a)*(1.d0 - (x*(1.d0+d*T))/a)
;   dydb = c*(a*exp(x*(1.d0+d*T)/a) - 1.d0)
;   dydc = b*(a*exp(x*(1.d0+d*T)/a) - 1.d0) + x*(1.d0+d*T)
;   dydd = c*x*T*(b*exp(x*(1.d0+d*T)/a) + 1.d0)
;   ;dydx = c*(1.d0 + d*T)*(b*exp(x*(1.d0+d*T)/a) + 1.d0)
; 
;   delta_y = abs(dyda)*delta_a + abs(dydb)*delta_b + abs(dydc)*delta_c + abs(dydd)*delta_d	;groessfehler
;   rel_delta_y = sqrt((delta_a/a)^2. + (delta_b/b)^2. + (delta_c/c)^2. + (delta_d/d)^2.)	;relative fehler


end