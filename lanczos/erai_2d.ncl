;---------------------------------------------------------------------------
; AEW Diagnostics Package
; NCSU Tropical Dynamics Group
;
;
;
; program to create band pass filtered fields by reading year by year
; currently configured for era interim files (one file per hour) on
; our server.
; cannot cross over years (needs update to do so). 
;
; Aiyyer 10/2015
;        updated 11/2015 to the faster netcdf write
;        updated 04/2019 for 2d fields see separate code for 3d fields
;---------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin


;
; set parameters the filter
;
  dt = .25 ;days per time step
  t1 = 12  ;days  (low frequency cutoff, expressed in time domain)
  t2 = 2   ;days  (high frequency cutoff, expressed in time domain)
  ihp     = 2                             ; band pass
  fca = dt/t1
  fcb = dt/t2
  sigma   = 1.0                           ; Lanczos sigma
  nWgt    = 175                           ; loose 115 points each end                            
  wgt     = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
;
; 
;
  icount = 0


  
  do iy = 2006,2006 ;1980,2013
   units = "hours since 1900-01-01 00:00:00";
   startDate = cd_inv_calendar(iy,01,01,0,0,0,units,0)
   endDate   = cd_inv_calendar(iy,12,31,0,0,0,units,0)
   print (startDate + " " + cd_calendar(startDate,-2))
   print (endDate + " " + cd_calendar(endDate,-2))

   date1 =    startDate  - (nWgt-1)/(2*dt)
   date2 =    endDate    + (nWgt+1)/(2*dt)
   date1@units = units
   date2@units = units
   
   print (date1 + " "  + cd_calendar(date1,-2))
   print (date2 + " "  + cd_calendar(date2,-2))
 

   exit
   





    print ( "Year = " + iy )
    sDate = iy*1000000 +  MMDDHH_start
    eDate = iy*1000000 +  MMDDHH_end
    allTimes = yyyymmddhh_time(iy,iy,6, "integer")
    TIME = allTimes({sDate:eDate})

;************************************************************
; for u and v
;    fils = "/typhoon/era/" + TIME + ".uv.grib"    

; for z and T
    fils = "/typhoon/era/" + TIME + ".geo.grib"    
;
;   fils = "/typhoon/era/" + TIME + ".omega.grib"    

;************************************************************

    in      = addfiles (fils, "r")
    ListSetType (in, "join")             ; concatenate or "merge" (default)

; make changes to variable name here
;   dat = in[:]->V_GDS4_ISBL(:,{50:1000},{-20:45},:)     
;   dat = in[:]->U_GDS4_ISBL(:,{50:1000},{-20:45},:)   
;   dat = in[:]->Z_GDS4_ISBL(:,{50:1000},{-20:45},:)  
;   dat = in[:]->W_GDS4_ISBL(:,{50:1000},{-20:45},:)


; temperature    
    dat = in[:]->T_GDS4_ISBL(:,{50:1000},{-20:45},:)  

    dat!0    = "time"  
    dat&time = TIME   ; cannot be string
    dat!1    = "lev"     
    dat!2    = "lat"     
    dat!3    = "lon"    
    
    FDAT      =  wgt_runave_n_Wrap ( dat,  wgt, 0, 0 )
    delete(dat)
 

; just retain a subset of the data since the ends of the filtered data are missing (due to weights)    
    sDate = iy*1000000 + MMDDHH_fiter_start
    eDate = iy*1000000 + MMDDHH_filter_end
    WBP = FDAT({sDate:eDate},:,:,:) 
    delete(FDAT)
    printVarSummary(WBP)    
    
; get the arrays for lev, lat and lon
    time= WBP&time
    lev = WBP&lev
    lat = WBP&lat
    lon = WBP&lon
    print ( time + " " + WBP(:,{650.},{10.},{345.} ) )
;
;    
;    printVarSummary(WBP)
    ntim  = dimsizes(time)                 ; get dimension sizes  
    klev  = dimsizes(lev)                                               
    nlat  = dimsizes(lat)  
    nlon  = dimsizes(lon)  
;    
;    
    setfileoption("nc","format","netcdf4") 
    outfile = "T2T10Day" + iy + ".nc"
    system( "rm " + outfile )

    fout    = addfile (outfile, "c" )

;===================================================================
; explicitly declare file definition mode. Improve efficiency.
;===================================================================
    setfileoption(fout,"DefineMode",True)
;===================================================================
; create global attributes of the file
;===================================================================
    fAtt               = True            ; assign file attributes
    fAtt@title         = "NCL Efficient Approach to netCDF Creation"  
    fAtt@source_file   =  fils(0)
    fAtt@Conventions   = "None"   
    fAtt@creation_date = systemfunc ("date")        
    fileattdef( fout, fAtt )            ; copy file attributes  
;===================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;===================================================================
    dimNames = (/"time", "lev", "lat", "lon"/)  
    dimSizes = (/ ntim,   klev,  nlat,  nlon/) 
    dimUnlim = (/ False , False, False, False/)   

    filedimdef(fout,dimNames,dimSizes,dimUnlim)
    
    filevardef(fout, "time" ,typeof(time),getvardims(time)) 
    filevardef(fout, "lev"  ,typeof(lev), getvardims(lev) )                           
    filevardef(fout, "lat"  ,typeof(lat), getvardims(lat))                          
    filevardef(fout, "lon"  ,typeof(lon), getvardims(lon))   
    
    filevardef(fout, "WBP"    ,typeof(WBP), getvardims(WBP))    
    
;===================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;===================================================================
    setfileoption(fout,"DefineMode",False)

    fout->TBP     = (/WBP/)
    fout->time   = (/time/)     
    fout->lev    = (/lev/)
    fout->lat    = (/lat/)
    fout->lon    = (/lon/) 
    
    delete(allTimes)
    delete(time)
    delete (WBP)
    delete(fils)
  end do
  
  
end



