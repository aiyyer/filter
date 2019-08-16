; a program to filter gridded data in space and time
;;
;
;
;
;
;---------------------------------------------------------------------------
load "./diagnostics_cam.ncl"


begin

  plev = 850


; set region to work on. Choose from the list below
; regions = (/africa,india,nwpac,global,equatorial/)

  region = "equatorial"
  pdfName = region + "_ST"

  region = "africa"
  pdfName = region + "_ST"


  taperx = True  ; taper along the x direction except if global data
  
  if (region .eq. "africa" ) then
    lonW =  -52 ; integer
    lonE =   55 ; integer  
; lat bounds
    latS = 5.
    latN = 20.
    FlipLon = True
  end if
  
  if (region .eq. "india" ) then
    lonW =  50  ; integer
    lonE = 150  ; integer
; lat bounds
    latS = 10.
    latN = 25.
    FlipLon = False
  end if

  if (region .eq. "nwpac" ) then
    latS = 10.
    latN = 25.
    lonW = 90   ; integer
    lonE = 190  ; integer
    FlipLon = False
  end if


; global
  if (region .eq. "global" ) then
    lonW = 0  ;-40 ; integer
    lonE = 360 ; 40 ; integer
; lat bounds
    latS = 5.
    latN = 25.
    FlipLon = False
    taperx = False
  end if
  
  if (region .eq. "equatorial" ) then
    lonW = 0  
    lonE = 360
; lat bounds
    latS = -15.
    latN =  15.
    FlipLon = False
    taperx = False
  end if


  
  Title = region


  ; set taper percent for both time and x directions
  taperPercent = 0.1

  

  sTime = 0615 ;050100
  eTime = 1015 ;113018

  spd = 4 ; 4 samples per day for 6 hourly data
  ystart = 1980
  nyears = 24

  samplewindow = 60*spd
  overlap      = 30*spd
 
  frequency      = fspan(-1*samplewindow/2,samplewindow/2,samplewindow+1)/(samplewindow/spd)
  frequency!0    = "frequency"
  frequency&frequency =  frequency
  
  itS   = 0
  count = 0
  firstTime = True

; set variable to operate on  
  ivar = 0


  dataDir = "/tamrind/era/3D/"

  do iy = ystart,ystart+nyears-1


    sDate = iy*10000 + sTime
    eDate = iy*10000 + eTime
    allTimes := yyyymmdd_time(iy,iy,"integer")
    TIME     := allTimes({sDate:eDate})
    ;print (TIME)
; now define all files 
    fils :=  dataDir + TIME +   ".nc"
    in   := addfiles (fils, "r")
    ListSetType (in, "cat")             ; concatenate or "merge" (default)
; now read the data 
    if (FlipLon) then
      dat := lonFlip(in[:]->V(:,{plev},{latS:latN},:))
    else
      dat := in[:]->V(:,{plev},{latS:latN},:)
    end if
    ;printVarSummary(dat)
    

; the array dat has all the data for the months of interest
; now we will 
    
  

; iterate 3 times in 61 day intervals and 30 day overlap
    spd = 4
    iB  = 0
 
    do itimes = 0,2
      iE  = iB + 61*spd - 1
      print ( iB + " " + iE)
      dat2 := dat(iB:iE,:,{lonW:lonE}) 
      ; taper in time
      dat2 = taper_n(dat2,taperPercent,0,0)

      ; taper in longitude only if not global data along x
      if (taperx) then
        dat2 = taper_n(dat2,taperPercent,0,2)
      end if

      dims = dimsizes(dat2)
      nt   = dims(0)
      ny   = dims(1)
      nx   = dims(2)      
      ndays = nt/spd
     ; printVarSummary(dat2)

; create the array to hold the transform data

      datx     := new( (/2,nx,nt/),      typeof(dat) )
      POWER    := new( (/ny,nx+1,nt+1/), typeof(dat) )
      


      do j = 0,ny-1 ; loop over latitudes     
; first the space transform in longitude direction      
        do l = 0,nt-1 ; loop over times
          datx(:,:,l) = cfftf (dat2(l,j,:), 0.0, 0)    ; imaginary part set to 0.0
        end do
        datx = datx/nx        
; next the transform in time for each longitudinal point      
        do i = 0,nx-1 ; loop over longs
          datx(:,i,:) = cfftf (datx(0,i,:),datx(1,i,:), 0)   
        end do
        datx = datx/nt
        ; split into westward and eastward components and calculate the power
        POWER(j,:,:) = resolveWavesHayashi( datx, ndays, spd )
      end do


      if (firstTime) then
        print ( "Creating array to hold average PS")
        POWERAve = POWER
        firstTime = False
      end if
      
      POWERAve = (POWERAve*count + POWER )/ (count + 1)
      count = count + 1
      
      delete(POWER)
      delete(datx)
      iB = iE -30*spd

      print ( "iy count = " + iy + " "  + count ) 
    end do

   
  end do
  

  POWER = POWERAve
  delete(POWERAve)
  printVarSummary(POWER)


; write out the data
;  Define the output file name and write some data
  fileName = "stPower_" + region +  ".nc"
  system( "rm " + fileName )
  outFile = addfile( fileName, "c" )


; write the timeseries to output file
  outFile->POWER=POWER

  
  
    
; be sure to redefine wavenumber so that it is relative to the full earth circumference
; wavenumber is defined to be an integer. So this is an approximate conversion
  wavenumber = POWER&wave
  wavenumber = wavenumber*360/(lonE-lonW)
  POWER&wave  = wavenumber
  delete(wavenumber)

  printVarSummary(POWER)
    
;-------------------------------------------------------------------
; plotting parameters freq and wavenumbers to plot
; much of the following is lifted from NCL's wkSpaceTime code
;------------------------------------------------------------------

  frqfftwin      = POWER&freq
  frqfftwin&freq = POWER&freq
  
  wavep1         = POWER&wave
  wavep1&wave    = POWER&wave
  
  minfrq4plt =  0.0
  maxfrq4plt =  0.5
  minwav4plt = -30.
  maxwav4plt =  30.
  
  minfrq     = minfrq4plt
  maxfrq     = min((/maxfrq4plt,max(frqfftwin)/))
  
  fillVal    = 1e20           ; miscellaneous
  
  
  freq       = frqfftwin({freq|minfrq:maxfrq})
  wavenumber = wavep1({wave|minwav4plt:maxwav4plt})
  NWVN       = dimsizes(wavenumber)         ; number of wavenumbers
  

  ; sum the power over all latitudes
  PAVE = dim_sum_n_Wrap(POWER,0)
  printVarSummary(PAVE)
; set the mean to missing
  PAVE(:,{0.0}) = (/PAVE@_FillValue /)  
  
  minwav4smth = -30 ;-27
  maxwav4smth =  30 ; 27

  do wv=minwav4smth,maxwav4smth
    wk_smooth121( PAVE({wv},nt/2+1:nt-1) )
  end do    


  
; create the background spectrum  
  BACK = dim_sum_n_Wrap(POWER,0)
  BACK@long_name = "Background Spectrum"  
  BACK@_FillValue        = fillVal    
  BACK(wave|:,freq|nt/2) = fillVal
  BACK@_FillValue        = fillVal
  
; smooth to get the background  
  do tt = nt/2+1,nt
    smthlen = maxwav4smth-minwav4smth+1
    if (frqfftwin(tt).lt.0.1) then
     do i = 1,5
        wk_smooth121( BACK(freq|tt,{wave|minwav4smth:maxwav4smth}) )
      end do
    end if
    if (frqfftwin(tt).ge.0.1.and.frqfftwin(tt).lt.0.2) then
      do i = 1,10
        wk_smooth121( BACK(freq|tt,{wave|minwav4smth:maxwav4smth}) )
      end do
    end if
    if (frqfftwin(tt).ge.0.2.and.frqfftwin(tt).lt.0.3) then
      do i = 1,20
         wk_smooth121( BACK(freq|tt,{wave|minwav4smth:maxwav4smth}) )
       end do
     end if     
     if (frqfftwin(tt).ge.0.3) then
       do i = 1,40
         wk_smooth121(BACK(freq|tt,{wave|minwav4smth:maxwav4smth}))
       end do
     end if
   end do
   
   do nw = minwav4smth,maxwav4smth ; smth frequency up to .8 cycles per day
     pt8cpd  = min((/closest_val(.8,frqfftwin),dimsizes(frqfftwin)-1/))
     smthlen = pt8cpd-(nt/2+1)+1
     do i = 1,10
       wk_smooth121( BACK({nw},nt/2+1:pt8cpd) )
     end do
   end do    
   
   BACK(wave|:,freq|nt/2) = fillVal
   BACK@long_name = "background spec"
   
   wks  = gsn_open_wks("pdf",pdfName)
   gsn_define_colormap(wks,"MPL_YlOrRd")
   
   res = True
   res@gsnSpreadColors      = True                  ; use full colormap
   tmFontHgtF            = 0.015     ; not sure why
   tiFontHgtF            = 0.018
   lbFontHgtF            = 0.015
   txFontHgtF            = 0.013
   
   res@gsnFrame          = False
   res@gsnMaximize       = True
   res@gsnPaperOrientation = "portrait"
   
   res@gsnLeftString     = "Westward"
   res@gsnRightString    = "Eastward"
   
                                ;res@lbBoxMinorExtentF = 0.18
   res@lbLabelFontHeightF= lbFontHgtF
   res@lbOrientation     = "vertical"
   
   res@cnFillOn          = True
   res@cnLinesOn     = False
   
   
   res@tmYLMode          = "Explicit"
   res@tmYLValues        = fspan(minfrq,maxfrq,9)
   res@tmYLLabels        = fspan(minfrq,maxfrq,9)
   res@tmYLMinorValues   = fspan(minfrq,maxfrq,17)
   
   res@tmYLLabelFontHeightF = tmFontHgtF
   res@tmXBLabelFontHeightF = tmFontHgtF
   
   res@tiXAxisString     = "Zonal Wave Number"
   res@tiXAxisFontHeightF= tiFontHgtF
   
   res@tiYAxisString     = "Frequency (cpd)"
   res@tiYAxisFontHeightF= res@tiXAxisFontHeightF
   
  
;-----------------------------------------------------------------------------   
   ; plot the log of the raw spectrum
   raw       = log(PAVE({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt}))
   raw!0     = "freq"
   raw&freq  =  freq 
   raw!1     = "wave"
   raw&wave  =  wavenumber
   raw@long_name = "Log of Raw Power"
   res@gsnCenterString = "Raw"
   plot = gsn_csm_contour(wks,raw,res)
   frame(wks)
;-----------------------------------------------------------------------------
; plot the log of the background spectrum
   backg       = log(BACK({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt}))
   backg!0     = "freq"
   backg&freq  =  freq 
   backg!1     = "wave"
   backg&wave  =  wavenumber
   backg@long_name = "Log of Background Power"
   res@gsnCenterString = "Background"

   plot = gsn_csm_contour(wks,backg,res)
   frame(wks)

;-----------------------------------------------------------------------------

; Now plot the signal relative to the background   
     
     raw    = PAVE({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
     backg  = BACK({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
     signal = raw
     signal = raw/backg

     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMinLevelValF       = 1.0               ; set min contour level
     res@cnMaxLevelValF       = 2.0                 ; set max contour level
     res@cnLevelSpacingF      =  .1               ; set contour spacing
     res@cnFillOn             = True               ; turn on color for contours
     ;res@cnLinesOn            = True              ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line label
     
     res@gsnCenterString = Title
     plot = gsn_csm_contour(wks,signal,res)
     
end