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




; trmm file name

;----------------------------------------------
;coarse 1-day aggregated trmm data 
  dir  = "/share01/trmm/carl_3b42/"
  fil  = dir + "trmm3b42_daily.nc"
  in   = addfile(fil, "r")
  rain = in->raindaily
  printVarSummary(rain)
;----------------------------------------------
;
;
; set parameters the filter
;
  dt = 1.0 ;.25 ;days per time step
  t1 = 12  ;days  (low frequency cutoff, expressed in time domain)
  t2 = 2   ;days  (high frequency cutoff, expressed in time domain)
;
;
;
  ihp     = 2                             ; band pass
  fca = dt/t1
  fcb = dt/t2
  sigma    = 1.0                           ; Lanczos sigma
  nWgt     = 121                           ; loose 120 points each end                            
  wgt      = filwgts_lanczos(nWgt, ihp, fca, fcb, sigma )
  FDAT     =  wgt_runave_n_Wrap(rain,wgt,0,0 )
  delete(rain)

  dimens = dimsizes(FDAT)
  nt = dimens(0)
  FDAT := lonFlip (FDAT(nWgt:nt-nWgt,:,:))
  printVarSummary(FDAT)


  outfile = dir + "trmm3b42_daily_2to12day_filtered.nc"

  system( "rm " + outfile )  
  fout    = addfile (outfile, "c" )

  
  fout->rain2to12=FDAT
  
end



