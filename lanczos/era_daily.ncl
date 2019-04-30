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
;
; Aiyyer 12/2018
;
; 
;---------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

  baseDir ="/mnt/tamrind/eradaily/3D/"
  baseDirOut = "/typhoon/projects/data/filtered/"

  variable = "V"
  variableOut = "VBP"


; set parameters the filter
  dt = 1.     ;days per time step
  t1 = 12.      ;days  (low frequency cutoff, expressed in time domain)
  t2 = 2.   ;days  (high frequency cutoff, expressed in time domain)
  ihp = 2     ;1=high pass 2=bp
  fca = dt/t1
  fcb = dt/t2
  sigma   = 1.0                           ;Lanczos sigma
  nWgt    = 121                           
  wgt     = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )

  
  do iy = 1998,2013
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
    dateArray := cd_calendar(dateArrayCoded,-2)
    ;print(dateArray) 
    fils := baseDir+dateArray+".nc"
    in   := addfiles (fils, "r")
    ListSetType (in, "join")             ; concatenate or "merge" (default)
    dat := in[:]->$variable$(:,:,:,:)  
    delete (dat@time_units)
    delete (dat@time)
    print ( "read data")
    dat!0    = "time" 
    dat&time = dateArrayCoded   ; cannot be string
    dat!1    = "lev"     
    dat!2    = "lat"     
    dat!3    = "lon"    
    
    FDAT     :=  wgt_runave_n(dat,  wgt, 0, 0 )
    print("filtered data")
    copy_VarMeta(dat,FDAT)

    FDAT := FDAT({startDate:endDate},:,:,:) 
    printVarSummary(FDAT)

    ;print( cd_calendar(FDAT&time,-2))
    setfileoption("nc","format","netcdf4") 

    outfile = baseDirOut+"2-12DayBP_daily" + iy + ".nc"     ; change file name to reflect variable
    if ( fileexists(outfile) ) then    
      fout    = addfile (outfile, "w" )
    else
      fout    = addfile (outfile, "c" )
    end if
    fout->$variableOut$=FDAT   ; change variable name 
    print("wrote data")
    
  end do
end
