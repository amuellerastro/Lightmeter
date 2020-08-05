pro sort_data_paranal

;choose dataset
  loc1 = 'Paranal'
  loc2 = 'CL'

;read in of setup file
  readcol, 'fix/SETUP.txt', dum, var, format='a,a', skipline=2
  pathbase = var[0]

  pathraw = pathbase+'data_original/'+loc1+'/'
  pathsort = pathbase+'data_sorted/'+loc1+'/'

  file_mkdir, pathsort
  ;pathraw already exist because there we copied already the raw data
  file_mkdir, pathraw+'old/'
  pathold = pathraw+'old/'	;in case it does not exist

;search and selection of raw file
  file = file_search(pathraw+'*.csv', count=nfiles)

  if (nfiles eq 0) then begin
    print, 'No files found. Stop.'
    return
  endif

  print, 'Following files were found: '

  for i=0L,nfiles-1 do begin
    print, uint(i+1), ' ', file[i]
  endfor

for i=0L,nfiles-1 do begin

  ;read in selected file, can be very large files, usage of readcol
    print, 'Processing File '+strcompress(i+1,/rem)+' of '+strcompress(nfiles,/rem)
    readcol, file[i], day, month, year, hh, mm, ss, t1, t2, dum1, counts, dum2, format='i,i,i,i,i,i,d,d,a,d,a', delimiter='.;:,', /silent

  ;remove last line because new file is started if there was an error and last line is broken (most of the time)
    nlines = file_lines(file[i])

    if (nlines gt 1) then begin

      day = day[0:nlines-2]
      month = month[0:nlines-2]
      year = year[0:nlines-2]
      hh = hh[0:nlines-2]
      mm = mm[0:nlines-2]
      ss = ss[0:nlines-2]
      t1 = t1[0:nlines-2]
      t2 = t2[0:nlines-2]
      counts = counts[0:nlines-2]


    ;merge temperature
      temp = t1+t2/10.0

    ;convert day,month,year,hh,mm,ss into strings (needed for, e.g. '2' -> '02')
      year_st = strtrim(string(year),2)
      month_st = strtrim(string(month),2)
      day_st = strtrim(string(day),2)
    ;   hh_st = strtrim(string(hh),2)
    ;   mm_st = strtrim(string(mm),2)
    ;   ss_st = strtrim(string(ss),2)

      ;get length of the strings
      lmonth = strlen(month_st)
      lday = strlen(day_st)
    ;   lhh = strlen(hh_st)
    ;   lmm = strlen(mm_st)
    ;   lss = strlen(ss_st)

      idxmonth = where(lmonth ne 2)
      if (idxmonth[0] ne -1) then month_st[idxmonth] = '0'+month_st[idxmonth]

      idxday = where(lday ne 2)
      if (idxday[0] ne -1) then day_st[idxday] = '0'+day_st[idxday]

    ;   idxhh = where(lhh ne 2)
    ;   if (idxhh[0] ne -1) then hh_st[idxhh] = '0'+hh_st[idxhh]
    ; 
    ;   idxmm = where(lmm ne 2)
    ;   if (idxmm[0] ne -1) then mm_st[idxmm] = '0'+mm_st[idxmm]
    ; 
    ;   idxss = where(lss ne 2)
    ;   if (idxss[0] ne -1) then ss_st[idxss] = '0'+ss_st[idxss]

    ;release memory
      dum1 = 0 & dum2 = 0
      lmonth = 0 & lday = 0 & idxmonth = 0 & idxday = 0
      ;lhh = 0 & lmm = 0 & lss = 0 & idxhh = 0 & idxmm = 0 & idxss = 0

    ;convert every date/time into JD
      jd = julday(month, day, year, hh, mm, ss)

    ;extract each day with data - a day starts at 12:00am UT (thats why 'floor' is used)
      jdround = floor(jd)
      idx = uniq(jdround)

      for j=0L,n_elements(idx)-1 do begin

	idxval = where(jdround eq jdround[idx(j)])


	;if record starts after midnight one has to substract a day to get coorect date of the night, and this could change year, month, day
	;I know, code looks crazy
	if (hh[idxval(0)] lt 12.) then begin

;print, '2nd half of night, modify day'

	  tmp_day = uint(day_st[idxval(0)])-1
	  tmp_day = strtrim(string(tmp_day),2)
	  day_st[idxval(0)] = tmp_day

	  if (strlen(tmp_day) ne 2) then begin

;print, 'put additional zero to day'

	    tmp_day = '0'+tmp_day
	    day_st[idxval(0)] = tmp_day

	    ;in case of a jump to the previous day
	    if (uint(tmp_day) lt 1) then begin

;print, 'handling of jumps in month'

	      ;feb -> jan
	      if (month[idxval(0)] eq 2) then begin
		month_st[idxval(0)] = '01'
		day_st[idxval(0)] = '31'
	      endif

	      ;mar -> feb
	      if (month[idxval(0)] eq 3) then begin
		month_st[idxval(0)] = '02'
		if (year[idxval(0)] ne 2010 or year[idxval(0)] ne 2016 or year[idxval(0)] ne 2020) then $
		  day_st[idxval(0)] = '28' else day_st[idxval(0)] = '29'
	      endif

	      ;apr -> mar
	      if (month[idxval(0)] eq 4) then begin
		month_st[idxval(0)] = '03'
		day_st[idxval(0)] = '31'
	      endif

	      ;may -> apr
	      if (month[idxval(0)] eq 5) then begin
		month_st[idxval(0)] = '04'
		day_st[idxval(0)] = '30'
	      endif

	      ;jun -> may
	      if (month[idxval(0)] eq 6) then begin
		month_st[idxval(0)] = '05'
		day_st[idxval(0)] = '31'
	      endif

	      ;jul -> jun
	      if (month[idxval(0)] eq 7) then begin
		month_st[idxval(0)] = '06'
		day_st[idxval(0)] = '30'
	      endif

	      ;aug -> jul
	      if (month[idxval(0)] eq 8) then begin
		month_st[idxval(0)] = '07'
		day_st[idxval(0)] = '31'
	      endif

	      ;sep -> aug
	      if (month[idxval(0)] eq 9) then begin
		month_st[idxval(0)] = '08'
		day_st[idxval(0)] = '31'
	      endif

	      ;oct -> sep
	      if (month[idxval(0)] eq 10) then begin
		month_st[idxval(0)] = '09'
		day_st[idxval(0)] = '30'
	      endif

	      ;nov -> oct
	      if (month[idxval(0)] eq 11) then begin
		month_st[idxval(0)] = '10'
		day_st[idxval(0)] = '31'
	      endif

	      ;dec -> nov
	      if (month[idxval(0)] eq 12) then begin
		month_st[idxval(0)] = '11'
		day_st[idxval(0)] = '30'
	      endif

	      ;jan -> dec
	      if (month[idxval(0)] eq 1) then begin
		month_st[idxval(0)] = '12'
		day_st[idxval(0)] = '31'
		year_st[idxval(0)] = strtrim(string(year[idxval(0)]-1),2)
	      endif


	    endif

	  endif

	endif

	;write data out
	;output file has the structure: country_location_YYYYMMDD.txt
	openw, lun, pathsort+loc2+'_'+loc1+'_'+strcompress(year_st[idxval(0)],/rem)+strcompress(month_st[idxval(0)],/rem)+strcompress(day_st[idxval(0)],/rem)+'.txt', width=1400, /append, /get_lun
    ;    printf, lun, '                JD    Temp.     Counts'
	for k=0L,n_elements(idxval)-1 do begin

	  ;printf, lun, jd[idxval(j)], day_st[idxval(j)], month_st[idxval(j)], year_st[idxval(j)], temp[idxval(j)], counts[idxval(j)], format='(f18.10,2a7,a9,f9.1,i11.0)'
	  printf, lun, jd[idxval(k)], temp[idxval(k)], counts[idxval(k)], format='(f18.10,f9.1,i11.0)'

	endfor

	close, lun
	free_lun, lun

      endfor

    endif

  ;move the log file into a separate folder
  spawn, 'mv '+file[i]+' '+pathold

endfor

end