
pro plot_lightcurve_single_night, loc1, loc2, path, file_day, file_month, file_year, hr, counts, $
				  t0rise, t0set, t6rise, t6set, t12rise, t12set, t18rise, t18set


timeplot = hr


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
  device, filename=path+loc2+'_'+loc1+'_'+strcompress(file_year,/rem)+strcompress(file_month,/rem)+strcompress(file_day,/rem)+'.ps',/color,XSIZE=20, YSIZE=15, XOffset=xoffset, YOffset=yoffset

    !P.Font=1
    !p.thick=4
    !x.thick=3
    !y.thick=3

    !p.multi=[0,1,2]

    multiplot, mXtitle = 'Time UT [hr]', mYtitle = 'Irradiance [W/m!U2!N]', mxTitSize=1.4, mxTitOffset=0.3, $
	    myTitSize=1.4, myTitOffset=0, gap=0.005, $;, ytickformat='exponent', $;, /square
	    mtitle=loc1+' '+strcompress(file_day,/rem)+'.'+strcompress(file_month,/rem)+'.'+strcompress(file_year,/rem), mTitSize=1.4, mTitOffset=-.3

      plot, timeplot, counts, xr=[0,24], xst=1, ytickformat='logticks_exp', charsize=1.4, $
	    /nodata, /ylog, /ynozero, yr=[1.d0,1.d3], yst=1, $
	    background=fsc_color('white'), Color=fsc_color('black'), xtickinterval=3, xminor=3
      oplot, timeplot, counts, color=fsc_color('blue')
      plots, [t0set,t0set],[1.d0,1.d3],/data, linestyle=2, color=fsc_color('yellow')
      plots, [t0rise,t0rise],[1.d0,1.d3],/data, linestyle=2, color=fsc_color('yellow')


    multiplot

      plot, timeplot, counts,xr=[0,24],xst=1, ytickformat='logticks_exp',$
	    /nodata, /ylog, /ynozero, yr=[5.d-6,1.D-3], yst=1, $	  
	    background=fsc_color('white'), Color=fsc_color('black'),xtickinterval=3, xminor=3, $
	    xtickname=['12','15','18','21','00','03','06','09','12'], charsize=1.4
      oplot, timeplot, counts, color=fsc_color('blue')
      plots, [t18set,t18set],[5.d-6,1.d-3],/data, linestyle=2, color=fsc_color('orange red')
      plots, [t18rise,t18rise],[5.d-6,1.d-3],/data, linestyle=2, color=fsc_color('orange red')
      plots, [t12set,t12set],[5.d-6,1.d-3],/data, linestyle=2, color=fsc_color('orange')
      plots, [t12rise,t12rise],[5.d-6,1.d-3],/data, linestyle=2, color=fsc_color('orange')


  device,/close
  set_plot,'x'

  multiplot, /reset

  !p.font=0
  !p.thick=1
  !x.thick=1
  !y.thick=1

  !p.multi = [0,1,0]

;spawn, 'gv '+path+loc2+'_'+loc1+'_'+strcompress(file_year,/rem)+strcompress(file_month,/rem)+strcompress(file_day,/rem)+'.ps'

; 
;   set_plot, 'ps'
;   device, filename=path+loc2+'_'+loc1+'_'+strcompress(file_year,/rem)+strcompress(file_month,/rem)+strcompress(file_day,/rem)+'.ps',/color,XSIZE=20, YSIZE=11.6, XOffset=xoffset, YOffset=yoffset
; 
;     !P.Font=1
;     !p.thick=4
;     !x.thick=3
;     !y.thick=3
; 
;     plot, timeplot, counts,xr=[0,24],xst=1,$
; 	  /nodata,/YLOG,/ynozero, ytitle='Irradiance [W/m!U2!N]',yr=[1.d-6,1.D3],yst=1,$
; 	  title=loc1+' '+strcompress(file_day,/rem)+'.'+strcompress(file_month,/rem)+'.'+strcompress(file_year,/rem),$
; 	  background=fsc_color('white'),Color=fsc_color('black'),xtitle='Time UT [hr]',xtickinterval=3,xminor=3,$
; 	  xtickname=['12','15','18','21','00','03','06','09','12'], charsize=1.3
; 
;     oplot, timeplot, counts, color=fsc_color('blue')
; 
;     plots, [t18set,t18set],[1.d-6,1.D3],/data, linestyle=2, color=fsc_color('orange red')
;     plots, [t18rise,t18rise],[1.d-6,1.D3],/data, linestyle=2, color=fsc_color('orange red')
;     plots, [t12set,t12set],[1.d-6,1.D3],/data, linestyle=2, color=fsc_color('orange')
;     plots, [t12rise,t12rise],[1.d-6,1.D3],/data, linestyle=2, color=fsc_color('orange')
; ;     plots, [t6set,t6set],[1.d-6,1.D3],/data, linestyle=2, color=fsc_color('yellow')
; ;     plots, [t6rise,t6rise],[1.d-6,1.D3],/data, linestyle=2, color=fsc_color('yellow')
;     plots, [t0set,t0set],[1.d-6,1.D3],/data, linestyle=2, color=fsc_color('yellow')
;     plots, [t0rise,t0rise],[1.d-6,1.D3],/data, linestyle=2, color=fsc_color('yellow')
; 
;   device,/close
;   set_plot,'x'
;   !p.font=0
;   !p.thick=1
;   !x.thick=1
;   !y.thick=1

;spawn, 'gv '+path+loc2+'_'+loc1+'_'+strcompress(file_year,/rem)+strcompress(file_month,/rem)+strcompress(file_day,/rem)+'.ps'

end