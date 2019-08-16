;--------------------------------------------------------------------------
;Program to read TRMM data and create power spectrum in (w,k) space
;for zonally limited regions. Hardcoded for JJAS months
;
;
;Much of the code is borrowed from NCL's wkSpaceTime function
;
;NCSU Tropical Dynamics
;
;---------------------------------------------------------------------------


begin
  
  wks  = gsn_open_wks("pdf","spaceTimeSpectra")
  gsn_define_colormap(wks,"sunshine_9lev")

  ;pltColorMap = "amwg_blueyellowred"
  ;gsn_define_colormap(wks,pltColorMap)


  resP                     = True                ; modify the panel plot
  resP@gsnPanelLabelBar    = True                ; add common colorbar
  resP@lbLabelFontHeightF  = 0.007               ; make labels smaller

  res = True
  res@gsnSpreadColors      = True                  ; use full colormap
  tmFontHgtF            = 0.015     ; not sure why
  tiFontHgtF            = 0.018
  lbFontHgtF            = 0.015
  txFontHgtF            = 0.013
    
  res@gsnDraw             = False           ; don't draw
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
  
  
  
  
  res@tmYLLabelFontHeightF = tmFontHgtF
  res@tmXBLabelFontHeightF = tmFontHgtF
   
  res@tiXAxisString     = "Zonal Wave Number"
  res@tiXAxisFontHeightF= tiFontHgtF
   
  res@tiYAxisString     = "Frequency (cpd)"
  res@tiYAxisFontHeightF= res@tiXAxisFontHeightF
  
  
  
; set region to work on. Choose from the list below
  regions = (/"africa","india","nwpac","equatorial"/)
  

  iregion = 0
  region = regions(iregion)
  fileName = "stPower_" + region +  ".nc"
  inFile = addfile( fileName, "r" )
  Title = regions(iregion)
  print (fileName)
  
; write the timeseries to output file
  POWER :=inFile->POWER
  
  if (region .eq. "africa" ) then
    lonW =  -50 ; integer
    lonE =   50 ; integer  
  end if
  if (region .eq. "india" ) then
    lonW =  50  ; integer
    lonE = 150  ; integer
  end if
  if (region .eq. "nwpac" ) then
    lonW = 90   ; integer
    lonE = 190  ; integer
  end if
  if (region .eq. "equatorial" ) then
    lonW = 0  ;-40 ; integer
    lonE = 360 ; 40 ; integer
  end if
  
  dims = dimsizes(POWER)
  nt = dims(2) - 1 ; -1 since power has dimension ny,nx+1,nt+1, 
  
  
  
;-------------------------------------------------------------------
; plotting parameters freq and wavenumbers to plot
; much of the following is lifted from NCL's wkSpaceTime code
;------------------------------------------------------------------
  
  frqfftwin      := POWER&freq
  frqfftwin&freq := POWER&freq
  
  wavep1         := POWER&wave
  wavep1&wave    := POWER&wave
  
  
  
  minfrq4plt =  0.0
  maxfrq4plt =  .8
  minwav4plt = -30.
  maxwav4plt =  30.
  
  minfrq     = minfrq4plt
  maxfrq     = min((/maxfrq4plt,max(frqfftwin)/))
  
  fillVal    = 1e20           ; miscellaneous
  
  
  freq       := frqfftwin({freq|minfrq:maxfrq})
  wavenumber := wavep1({wave|minwav4plt:maxwav4plt})
  NWVN       = dimsizes(wavenumber)         ; number of wavenumbers
  
  
                                ; sum the power over all latitudes
  PAVE := dim_sum_n_Wrap(POWER,0)
  printVarSummary(PAVE)
; set the mean to missing
  PAVE(:,{0.0}) = (/PAVE@_FillValue /)  
  
  minwav4smth = -30 ;-27
  maxwav4smth =  30 ; 27
  
  do wv=minwav4smth,maxwav4smth
    wk_smooth121( PAVE({wv},nt/2+1:nt-1) )
  end do    
  
  
  
; create the background spectrum  
  BACK := dim_sum_n_Wrap(POWER,0)
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
  
  
  res@tmYLMode          = "Explicit"
  res@tmYLValues        = fspan(minfrq,maxfrq,9)
  res@tmYLLabels        = fspan(minfrq,maxfrq,9)
  res@tmYLMinorValues   = fspan(minfrq,maxfrq,17)   
  
; Now plot the signal relative to the background   
  
  
; be sure to redefine wavenumber so that it is relative to the full earth circumference
; wavenumber is defined to be an integer. So this is an approximate conversion
  
  
  if ( (region .ne. "equatorial") .and. (region .ne. "global") ) then
    printVarSummary(PAVE)
    wavenumber := PAVE&wave
    wavenumber = wavenumber*360/(lonE-lonW)
    PAVE&wave  = wavenumber
    BACK&wave  = wavenumber
    
    delete(wavenumber)
  end if
  
  raw    := PAVE({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
  backg  := BACK({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
  signal := raw
  signal = raw/backg
  
  
  
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 1.0               ; set min contour level
  res@cnMaxLevelValF       = 3.0                 ; set max contour level
  res@cnLevelSpacingF      =  .2               ; set contour spacing
  
  
                                ;res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
                                ;res@cnLevels = (/0.8,0.9 ,1.0,1.1,1.15,1.2,1.25 \
                                ; ,1.3,1.35,1.4,1.45,1.5,1.6,18,2./)
  
  res@cnFillOn             = True               ; turn on color for contours
  res@cnLinesOn            = False              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off contour line label
; res@lbLabelBarOn = False
  
  res@gsnCenterString = Title
  plot = gsn_csm_contour(wks,signal,res)
  draw(plot)

    
end

