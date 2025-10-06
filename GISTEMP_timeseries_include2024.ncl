; GISTEMP_timeseries_include2023.ncl
; last updated April 16, 2025
; coded by: Jennifer.E.Kay@colorado.edu

;********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

begin

;; ######## BEGIN USER MODS

    ;plot_type = "x11"  ;; x11 testing
    plot_type = "ps" ;;; save graphs to postscript files
    plot_name="GISTEMP_plot_global_arctic_timeseries_through2024"
    makedata="yes" ;;; first time through - set this to "yes".  afterwards, set to "no"
    fname="GISTEMP_timeseries_through2024.nc"

;; ########## END USER MODS 

;; #########################################################################
;; ###  BEGIN PROGRAM, NO MORE USER MODS UNLESS YOU ARE CHANGING PROGRAM ###
;; #########################################################################

if makedata.eq."yes" then

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; READ IN THE GISTEMP OBSERVATIONS, MAKE TIMESERIES ;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Observations were downloaded March 2022 
;;; Data from http://www.esrl.noaa.gov/psd/data/gridded/data.gistemp.html

   ;;; read in the data - anomaly surface temperature with respect to 1951-1980
   giss_obs = addfile("air.2x2.1200.mon.anom.comb.nc","r")
   obs_TS = giss_obs->air
   lat_obs = giss_obs->lat
   lon_obs = giss_obs->lon
   time_foo = giss_obs->time

   ;;; convert the time to a usable year unit and find the years covered in the dataset
   time_decode=calendar_decode2(time_foo,0)
   minyear=min(time_decode(:,0))
   maxyear=max(time_decode(:,0))-1  ;;note: I put the minus 1 here because not all of the last year (2025) is available
   endidx=floattointeger((maxyear-minyear+1)*12-1)
   years_decode=time_decode(0:endidx,0)
   years_obs=fspan(minyear,maxyear,floattointeger((maxyear-minyear+1))) ;;; years from website
   print("minyear: "+minyear)
   print("maxyear: "+maxyear)

   ;; select the data only for full years and create a data array
   foo=ndtooned(obs_TS(0:endidx,:,:))
   foo1=onedtond(foo,(/floattointeger((maxyear-minyear+1)),12,90,180/))
   TS_ann_anom_obs=dim_avg_n(foo1,1)
   delete(foo)
   delete(foo1)
   TS_ann_anom_obs!0="year"
   TS_ann_anom_obs!1="lat"
   TS_ann_anom_obs!2="lon"
   TS_ann_anom_obs&lat=lat_obs
   TS_ann_anom_obs&lon=lon_obs
   print("min TS anomaly (deg C): "+min(TS_ann_anom_obs))
   print("max TS anomaly (deg C): "+max(TS_ann_anom_obs))

   ;; calculate the cosine weight for the observations
   if (typeof(lat_obs) .eq. "double") then
    wgt_obs = doubletofloat(NormCosWgtGlobe(lat_obs))
   else
    wgt_obs = NormCosWgtGlobe(lat_obs)
   end if
   wgt_obs!0="lat"
   wgt_obs&lat=lat_obs

   ;; calculate the annual mean global surface temperature anomaly
   gts_ann_anom_obs=dim_avg(dim_avg_wgt_n(TS_ann_anom_obs,wgt_obs,1,1))
   gts_ann_anom_obs!0="year"
   gts_ann_anom_obs&year=years_obs
   print("min global TS anomaly (deg C): "+min(gts_ann_anom_obs))  ;;; should be around -0.4 based on http://data.giss.nasa.gov/gistemp/graphs_v3/
   print("max global TS anomaly (deg C): "+max(gts_ann_anom_obs))  ;;; should be around 0.8 based on http://data.giss.nasa.gov/gistemp/graphs_v3/

   ;; calulate the annual mean Arctic surface temperature anomaly
   foo=TS_ann_anom_obs(:,{70:90},:)
   arcts_70to90_ann_anom_obs=dim_avg(dim_avg_wgt_n(foo,wgt_obs({70:90}),1,1))
   arcts_70to90_ann_anom_obs!0="year"
   arcts_70to90_ann_anom_obs&year=years_obs
   delete(foo)
   print("min arctic TS (70-90) anomaly (deg C):"+min(arcts_70to90_ann_anom_obs))
   print("max arctic TS (70-90) anomaly (deg C):"+max(arcts_70to90_ann_anom_obs))

   ;; calulate the annual mean Arctic surface temperature anomaly
   foo=TS_ann_anom_obs(:,{60:90},:)
   arcts_60to90_ann_anom_obs=dim_avg(dim_avg_wgt_n(foo,wgt_obs({60:90}),1,1))
   arcts_60to90_ann_anom_obs!0="year"
   arcts_60to90_ann_anom_obs&year=years_obs
   delete(foo)
   print("min arctic (60-90) TS anomaly (deg C):"+min(arcts_60to90_ann_anom_obs))
   print("max arctic (60-90) TS anomaly (deg C):"+max(arcts_60to90_ann_anom_obs))

   ;; calulate the annual mean Antarctic (60-90 S) surface temperature anomaly
   foo=TS_ann_anom_obs(:,{-90:-60},:)
   antts_60to90_ann_anom_obs=dim_avg(dim_avg_wgt_n(foo,wgt_obs({-90:-60}),1,1))
   antts_60to90_ann_anom_obs!0="year"
   antts_60to90_ann_anom_obs&year=years_obs
   delete(foo)
   print("min antarctic (60-90) TS anomaly (deg C):"+min(antts_60to90_ann_anom_obs))
   print("max antarctic (60-90) TS anomaly (deg C):"+max(antts_60to90_ann_anom_obs))

   ;; calulate the annual mean Antarctic (60-90 S) surface temperature anomaly
   foo=TS_ann_anom_obs(:,{-75:-55},:)
   antts_55to75_ann_anom_obs=dim_avg(dim_avg_wgt_n(foo,wgt_obs({-75:-55}),1,1))
   antts_55to75_ann_anom_obs!0="year"
   antts_55to75_ann_anom_obs&year=years_obs
   delete(foo)
   print("min antarctic (55-75) TS anomaly (deg C):"+min(antts_55to75_ann_anom_obs))
   print("max antarctic (55-75) TS anomaly (deg C):"+max(antts_55to75_ann_anom_obs))


;;;; WRITE EVERYTHING TO NETCDF FILES FOR QUICK PLOTTING LATER

   print(fname)
   system("/bin/rm -f " + fname)    ; remove if exists
   fout  = addfile (fname, "c")  ; open output file
   filedimdef(fout,"case",-1,True)
   fout->gts_ann_anom_obs=gts_ann_anom_obs
   fout->arcts_70to90_ann_anom_obs=arcts_70to90_ann_anom_obs
   fout->arcts_60to90_ann_anom_obs=arcts_60to90_ann_anom_obs
   fout->antts_60to90_ann_anom_obs=antts_60to90_ann_anom_obs
   fout->antts_55to75_ann_anom_obs=antts_55to75_ann_anom_obs
   fout->minyear=minyear
   fout->maxyear=maxyear
   ;attributes:
   fAtt = True
   fAtt@source = "GISTEMP_timeseries.ncl"
   fAtt@creation_date = systemfunc("date")
   fAtt@case = "GISTEMP"
   fileattdef(fout,fAtt)

;;; save to plotting variable names
   plotdata_gts_TS_obs=gts_ann_anom_obs
   plotdata_ats_TS_obs=arcts_70to90_ann_anom_obs

print("wrote the file")

else

   fplotdata=addfile(fname,"r")
   plotdata_gts_TS_obs=fplotdata->gts_ann_anom_obs
   plotdata_ats_TS_obs=fplotdata->ats_ann_anom_obs
   minyear=fplotdata->minyear
   maxyear=fplotdata->maxyear

end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; PLOT SURFACE TEMPERATURE TIMESERIES - 1880-maxyear ;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   plot_name=plot_name
   plot_type=plot_type
   print("Making...." + plot_name)
   wks = gsn_open_wks(plot_type,plot_name)
   gsn_define_colormap(wks,"rainbow+gray")
   plot = new(2,graphic)  ;; 2 timeseries panels - one Global; one Arctic
   res  = True
   res@vpHeightF= 0.4                    ; change aspect ratio of plot and slide it over so it fits on 
   res@vpWidthF = 0.6                    ; the screen
   res@tiMainFontHeightF = 0.015
   ;res@tiXAxisString = "Year"
   res@tiYAxisString = "Global"
   res@xyDashPattern = (/0/)
   res@xyLineColor = (/"black"/)
   res@xyLineThicknessF = (/6./)
   res@gsnFrame = False
   res@gsnDraw = False
   res@trYMinF    = -2.0
   res@trYMaxF    = 2.0
   res@tmYLMode          = "Explicit"              ; explicit labels
   res@tmYLValues        = (/-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0/)          ; location of labels 
   res@tmYLLabels        = (/-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0/)         ; labels themselves
   res@trXMinF    = minyear
   res@trXMaxF    = maxyear

   res1  = res
   res1@tiYAxisString = "Arctic (red) Global (black)"
   res1@trYMinF    = -4.0
   res1@trYMaxF    = 4.0
   res1@xyLineColor = (/"red"/)
   res1@tmYLValues        = (/-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0/)          ; location of labels 
   res1@tmYLLabels        = (/-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0/)         ; labels themselves

   plot(0) = gsn_csm_xy(wks,plotdata_gts_TS_obs&year,plotdata_gts_TS_obs,res) 
   plot(1) = gsn_csm_xy(wks,plotdata_ats_TS_obs&year,plotdata_ats_TS_obs,res1) 
   reso=res
   reso@xyLineColor=(/"black"/)
   plot1 = gsn_csm_xy(wks,plotdata_gts_TS_obs&year,plotdata_gts_TS_obs,reso) 
   overlay(plot(1),plot1)


   resP                     = True                ; modify the panel plot
   resP@txString            = "Surface Temperature Anomaly, deg C (1951-1980 base)"
   resP@lbLabelFontHeightF  = 0.007               ; make labels smaller
   gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot



;;;; WRITE OUT A CSV FILE WITH THE DATA

;;print(plotdata_gts_TS_obs&year)
;;plotdata_gts_TS_obs
;;plotdata_ats_TS_obs

year=plotdata_gts_TS_obs&year
printVarSummary(year)

 alist = [/year,plotdata_gts_TS_obs,plotdata_ats_TS_obs/]

 header = (/"--------------------------------", \
            "Annual Mean Surface Temperature Anomaly Observations", \
            "Column #1 = Global Surface, units = deg C", \
            "Column #2 = Arctic Surface Temperature (70-90 N), units = deg C", \
            "Source: GISTEMP, Hansen et al. 2010; anomaly calculated with respect to 1951-1980 base period", \
            "Data were downloaded April 2025 from http://www.esrl.noaa.gov/psd/data/gridded/data.gistemp.html",\
            "Data munged by Jennifer.E.Kay@colorado.edu using GISTEMP_timeseries_through2024.ncl, last update April 2025", \
            "--------------------------------"/)

 hlist = [/header/]

 fname = "GISTEMP_Arctic_Global_SurfaceTemperature_1880-2024.txt"
 write_table(fname, "w", hlist, "%s ")
 write_table(fname, "a", alist, "%16.0f,%16.2f,%16.2f")

exit



end
