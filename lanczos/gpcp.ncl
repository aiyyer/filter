;---------------------------------------------------------------------------
;
; NCSU Tropical Dynamics Group
;
;
;
; program to create band pass filtered fields 
; Aiyyer 10/2015
;        updated 11/2015 to the faster netcdf write
;        updated 11/2018 to work with trmmm data
;
;
;---------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
;
  baseDir = "/tamrind/gpcp/data/Daily_Data/"
 ;
; set parameters the filter
;
  dt = 1.0 ;.25 ;days per time step
  t1 = 10  ;days  (low frequency cutoff, expressed in time domain)
  t2 = 2   ;days  (high frequency cutoff, expressed in time domain)
;
;
;
  ihp     = 2                             ; band pass
  fca = dt/t1
  fcb = dt/t2
  sigma    = 1.0                           ; Lanczos sigma
  nWgt     = 151                           ; loose 120 points each end                            
  wgt      = filwgts_lanczos(nWgt, ihp, fca, fcb, sigma )
;  FDAT     =  wgt_runave_n_Wrap(rain,wgt,0,0 )
;  delete(rain)
  
  first = True
  units = "days since 1970-01-01 00:00:00 0:00"

  do iy = 1998,2015
    startDate = cd_inv_calendar(iy,01,01,0,0,0,units,0)
    endDate   = cd_inv_calendar(iy,12,31,0,0,0,units,0)
    date1 =    startDate  - (nWgt-1)/(2*dt)
    date2 =    endDate    + (nWgt+1)/(2*dt)
    date1@units = units
    date2@units = units
    
    ndays = tointeger(date2-date1+1)
    print (date1 + " "  + cd_calendar(date1,-2))
    print (date2 + " "  + cd_calendar(date2,-2))
    print ( "ndays , nWgt = " + ndays + " " + nWgt)

    dateArrayCoded := fspan(date1,date2,tointeger(ndays))
    dateArrayCoded@units = units
    dateTemp := cd_calendar(dateArrayCoded,0)
    id := tointeger(dateTemp(:,2))
    im := tointeger(dateTemp(:,1))
    yy := tointeger(dateTemp(:,0))
    fPre :="gpcp_daily_cdr_v13rA1_y"+yy+"_m"
    yyyymm := yy + "" + sprinti("%0.2i",im)
    fils :=baseDir + yyyymm  + "/" +  fPre + sprinti("%0.2i",im) + "_d" + sprinti("%0.2i",id)+".nc"

    in   := addfiles (fils, "r")
    ListSetType (in, "join")             ; concatenate or "merge" (default)
    rain := in[:]->precip
    if (first) then
      fName=fils(0)
      in   := addfile (fName, "r")
      lat = in->latitude
      lon = in->longitude
      first = False
    end if

    rain!0    = "time"
    rain&time = dateArrayCoded
    rain!1    = "lat"
    rain!2    = "lon"
    rain&lat = lat
    rain&lon = lon
    ;printVarSummary(rain)
    rain := rain(:,{-35:35},:)
    ;printVarSummary(rain)

; some data checks
    if ( any(ismissing(rain)) ) then
      print ("data missing")
      exit
    end if
    
    numM = num(rain.lt. 0.)
    print (iy +" "  + numM)
    if (numM .gt. 0 ) then
      print ("Unexpected < 0 rain value/s found")
      print ("Abort Code")
      exit
    end if      
    
    FDAT   :=  wgt_runave_n_Wrap(rain,wgt,0,0 )
    dimens  =  dimsizes(FDAT)
    nt      =  dimens(0)


    rainBP :=  FDAT({startDate:endDate},:,:)

    fAtt               = True            ; assign file attributes
    fAtt@title         = "2-10 day BP filtered gpcp precip"
    fAtt@source_file   =  fils(0)
    fAtt@creation_date = systemfunc ("date")        
    

    ;print ( cd_calendar(rainBP&time,-3) + " " + rainBP(:,{0.},{10.}))
    ;exit

    setfileoption("nc","format","netcdf4") 
    outfile = "gpcp_2to10day_" + iy + ".nc"
    system( "rm " + outfile )
    fout    = addfile (outfile, "c" )
    fileattdef( fout, fAtt )            ; copy file attributes  
    fout->rainBP=rainBP




   ;exit
  end do














 ; dimens = dimsizes(FDAT)
 ; nt = dimens(0)
 ; FDAT := lonFlip (FDAT(nWgt:nt-nWgt,:,:))
 ; printVarSummary(FDAT)


 ; outfile = dir + "trmm3b42_daily_2to12day_filtered.nc"

 ; system( "rm " + outfile )  
 ; fout    = addfile (outfile, "c" )

  
 ; fout->rain2to12=FDAT
  
end



