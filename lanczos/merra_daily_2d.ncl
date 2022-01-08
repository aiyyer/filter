;---------------------------------------------------------------------------
; AEW Diagnostics Package
; NCSU Tropical Dynamics Group
;
;
;
; Reads merra sfc flux data and averages to get daily data and then
; calculates the band=pass filtered fields
;
;---------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

 baseDir = "/typhoon/merra/sfcFlux/"
 baseDirOut = "/lotus/projects/data/nasa/merra/"
 inVar = "EFLUX"
 outVar = inVar+"2to10"

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

  
  do iy = 2021,2021

; define timeunits for ease of calculation
; use the same units as merra-2
    units = "days since 1900-01-01 00:00:00";

    startDate = cd_inv_calendar(iy,07,01,0,0,0,units,0)
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
    dateArrayCoded := fspan(date1,date2,tointeger(ndays))
    dateArrayCoded@units = units
    print ( "ndays , nWgt = " + ndays + " " + nWgt)

    date = date1
    firstTime = True

    do iday = 0,ndays-1
      print(iday + " " + date + " " +  cd_calendar(date,-2) )
      fileName = systemfunc("ls " + baseDir + "*" + cd_calendar(date,-2) + "*.nc4" )
      print(fileName)
      in   := addfile (fileName, "r")

; read and average the data over the 24 1-hourly samples and get daily average
      dat = dim_avg_n_Wrap(in->$inVar$,0) ; for sfc variables
      if ( firstTime) then
        dims = dimsizes(dat)
        ny = dims(0)
        nx = dims(1)
        datFilt := new( (/ndays,ny,nx/), typeof(dat) )
        printVarSummary(datFilt)
        firstTime = False
      end if      
      datFilt(iday,:,:) = dat
      date = date + 1.0
    end do

    dat!0    = "lat"     
    dat!1    = "lon" 
    
; apply the filter    
    datFilt     :=  wgt_runave_n(datFilt,  wgt, 0, 0 )
    print("filtered data")
    
    datFilt!1    = "lat"     
    datFilt!2    = "lon" 
    datFilt!0    = "time"     

    datFilt&time =  dateArrayCoded 
    datFilt&time@units = units
    datFilt&lat = dat&lat
    datFilt&lon = dat&lon
    datFilt@standard_name = dat@standard_name
    ;printVarSummary(datFilt)
    datFilt := datFilt({startDate:endDate},:,:)
    printVarSummary(datFilt)
    print(datFilt&time + " " + cd_calendar(datFilt&time,-3))


    setfileoption("nc","format","netcdf4") 
    outfile = baseDirOut+"2-10DayBP_daily_flux" + iy + ".nc"     ; change file name to reflect variable
    if ( fileexists(outfile) ) then    
      fout    = addfile (outfile, "w" )
    else
      fout    = addfile (outfile, "c" )
    end if
    fout->$outVar$=datFilt   ; change variable name 
    print("wrote data")
    
  end do
end
