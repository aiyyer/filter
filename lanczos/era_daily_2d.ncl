;---------------------------------------------------------------------------
; AEW Diagnostics Package
; NCSU Tropical Dynamics Group
;
;
;
; program to create high pass filtered fields by reading year by year
; currently configured for daily era interim files (one file per day) on
; our server.
;
; Reads the 2D fields in erai daily data files
;
; Aiyyer 06/04/2019
;
; updated 07/12/2019 to fix time units
;---------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

  baseDir ="/tamrind/eradaily/2D_derived/"
  baseDirOut = "/tamrind/eradaily/filtered/2to10/"

  variable    = "ColAvgQ"
  variableOut = "ColAvgQBP"


; set parameters the filter
  dt = 1.     ;days per time step
  t1 = 10.    ;days  (low frequency cutoff, expressed in time domain)
  t2 = 2.    ;days  (high frequency cutoff, expressed in time domain)
  ihp = 2     ;1=high pass 2=bp
  fca = dt/t1
  fcb = dt/t2
  sigma   = 1.0                           ;Lanczos sigma
  nWgt    = 121                           
  wgt     = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )

  
  do iy = 1998,2014

; define timeunits for ease of calculation
; this is not the actual timeunits of daily era data in our server!
    units = "days since 1900-01-01 00:00:00";

    startDate = cd_inv_calendar(iy,06,01,0,0,0,units,0)
    endDate   = cd_inv_calendar(iy,09,30,0,0,0,units,0)
    
; number of days to be processed
    NDates = endDate - startDate + 1
    print (startDate + " " + cd_calendar(startDate,-2))
    print (endDate + " " + cd_calendar(endDate,-2))
    print (NDates)
    
    date1 = startDate  - (nWgt-1)/2
    date2 =    endDate + (nWgt+1)/2
    date1@units = units
    date2@units = units
    ndays = tointeger(date2-date1+1)
    print (date1 + " "  + cd_calendar(date1,-2))
    print (date2 + " "  + cd_calendar(date2,-2))
    print ( "ndays , nWgt = " + ndays + " " + nWgt)
    dateArrayCoded := fspan(date1,date2,tointeger(ndays))
    dateArrayCoded@units = units

; the raw era files are arranged daily. we now createthe file list
    dateArray := cd_calendar(dateArrayCoded,-2)
    ;print(dateArray) 
    fils := baseDir+dateArray+".nc"
; fils is the filelist. next open and read them
;
    in   := addfiles (fils, "r")
    ListSetType (in, "join")             ; concatenate or "merge" (default)
    dat := in[:]->$variable$ 
    dat!0    = "time" 
    dat!1    = "lat"     
    dat!2    = "lon" 


; the era daily raw files have units of "hours since 1-1-1 00:00:0.0"
;    datTimeUnits = dat@time_units
; for now delete these:
 
    delete (dat@time_units)
    delete (dat@time)
    print ( "read data")
    
    
    dat&time = dateArrayCoded   ; cannot be string
    
; apply the filter    
    FDAT     :=  wgt_runave_n(dat,  wgt, 0, 0 )
    print("filtered data")
    copy_VarMeta(dat,FDAT)
    
    
; subset the FDAT to retain only days of filtered output desired by user
    FDAT := FDAT({startDate:endDate},:,:) 
;   printVarSummary(FDAT&time)
;   print( cd_calendar(FDAT&time,-2))
    setfileoption("nc","format","netcdf4") 
    outfile = baseDirOut+"2-10DayBP_daily" + iy + ".nc"     ; change file name to reflect variable
    if ( fileexists(outfile) ) then    
      fout    = addfile (outfile, "w" )
    else
      fout    = addfile (outfile, "c" )
    end if
    fout->$variableOut$=FDAT   ; change variable name 
    print("wrote data")
    
  end do
end
