
begin
  level1 = 850
  level2 = 850

; set parameters the filter
  npd = 4       ; number of timesteps per day of input data
  dt = 1./npd   ;days per time step
  t1 = 10.  ;days  (low frequency cutoff, expressed in time domain)
  t2 = 2.   ;days  (high frequency cutoff, expressed in time domain)
  ihp = 2   ;1=high pass 2=bp
  fca = dt/t1
  fcb = dt/t2
  sigma = 1.0                           ;Lanczos sigma

; number of days for buffer on either side
   bufDays = 5
;   
   nWgt    = npd*bufDays*2 + 1                        
   wgt     = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
   swgt = sum(wgt)

; the era data are in daily files
   baseDir    =  "/tamrind/era/3D/"
   baseDirOut = "/tamrind/era/filtered/2to10/"
   filtprefix = "2-10DayBP"
   
   variable = "U"
   variableOut = "UBP"
   
   
   timeUnitsEra = "hours since 1-1-1 00:00:0.0"
   
   do iy = 2015,2015
     
; need the units in days since the files are daily  
     units     = "days since 1900-01-01 00:00:00"
     dt        = 24./npd
     startDate = cd_inv_calendar(iy,06,01,0,0,0,units,0)
     endDate   = cd_inv_calendar(iy,09,30,18,0,0,units,0)
     nDaysFilt =  tointeger(endDate-startDate) + 1

     date1 = startDate  - bufDays 
     date2 = startDate  + bufDays
     ndays = tointeger(date2-date1+1)
     print(ndays) ; ndays = nfiles
     
; create a list of filenames   
     dateArrayCoded := fspan(date1,date2,tointeger(ndays))
     dateArrayCoded@units = units
     dateArray := cd_calendar(dateArrayCoded,-2)
     eraFiles =  baseDir + dateArray +".nc"
     
     in      = addfiles(eraFiles, "r")
     ListSetType (in, "cat")             ; concatenate or "merge" (default)
     dat := in[:]->$variable$(:,{level1:level2},:,:)  
     printVarSummary(dat)
     print ( "read data")
     dat!0    = "time" 
     dat!1    = "lev"     
     dat!2    = "lat"     
     dat!3    = "lon"       
     print(dat&time + " " + cd_calendar(dat&time,-3) )
     print (nWgt  + " " +ndays ) 
     
     FDAT = dat(0:3,:,:,:)
     printVarSummary(FDAT)
     
     do ih = 0,3
       FDAT(ih,:,:,:) = dim_sum_n(dat(ih:nWgt+ih-1,:,:,:)*conform(dat(ih:nWgt+ih-1,:,:,:),wgt,0),0 )
       print ( dat&time( (nWgt-1)/2 + ih ) + " " + cd_calendar( dat&time( (nWgt-1)/2 + ih ), -3 ) + " " + dat( (nWgt-1)/2 + ih, 0, {10.}, {10.} ) )
       FDAT&time(ih) = dat&time((nWgt-1)/2 + ih )
       print (FDAT&time(ih) + " " + cd_calendar( FDAT&time(ih) , -3 )  + " " + FDAT(ih, 0, {10.}, {10.} ) )
     end do

     datFilt = FDAT


; Having done day 0 now lets move ahead to read the next day's data
     dateCoded = dateArrayCoded(ndays-1)

     print (dateArrayCoded + " " + cd_calendar(dateArrayCoded,-2) )

     do idays = 1,nDaysFilt-1
       dateCoded = dateCoded + 1.0
       ;print ( dateCoded + " " + cd_calendar(dateCoded,-2) )
       ;print ( "----------")
       ;print (dat&time + " " + cd_calendar(dat&time,-3) + " " +  dat(:,0, {10.}, {10.} ) )
        ; print ( "----------")



       dat(0:(ndays-1)*npd -1,:,:,:)  = dat (npd:(ndays*npd-1),:,:,:)
       print (dat&time + " " + cd_calendar(dat&time,-3) + " " +  dat(:,0, {10.}, {10.} ) )

       exit




     end do
     

  end do

  



end

