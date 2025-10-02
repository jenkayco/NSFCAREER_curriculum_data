; plot_anom_tseries_include2025.ncl
; written by Jennifer.E.Kay@colorado.edu
; last modified: October 2, 2025
;; includes improved weighting from Norman Loeb implemented in 2014!

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
wks_type = "ps"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;  CERES-EBAF OBSERVATIONS - TOA  ;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fname_ceres = systemfunc("ls /glade/u/home/jenkay/small_data/CERES_EBAF/CERES_EBAF-TOA_Ed4.2.1_Subset_200003-202502.nc")
;;; I downloaded data directly from: https://asdc.larc.nasa.gov/project/CERES/CERES_EBAF_Edition4.1

f = addfile(fname_ceres,"r")
lat=f->lat
lon=f->lon
FSUTOA=f->toa_sw_all_mon
FSDTOA=f->solar_mon

year=(/2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025/)
yearstart=2000
nyrs=max(year)-yearstart+1
sizefoo=dimsizes(FSUTOA)
;print(sizefoo)
nlat=sizefoo(1)
nlon=sizefoo(2)
wgty=NormCosWgtGlobe(FSUTOA&lat)/2
wgty!0="lat"
wgty&lat=FSUTOA&lat
printVarSummary(FSUTOA)
print(nyrs)

;;; calculate the montly average net TOA fluxes
;;; pad the array so that it starts in january and ends in december
FSUTOApad=new((/nyrs*12,sizefoo(1),sizefoo(2)/),float)
FSUTOApad@_FillValue=-999
FSUTOApad=FSUTOApad@_FillValue
FSUTOApad(2:sizefoo(0)-1+2,:,:)=FSUTOA
FSUTOApad&time=ispan(0,nyrs*12-1,1)
FSUTOApad&lat@units="degrees_north"
FSUTOApad&lon@units="degrees_east"
FSUTOA_month=clmMonTLL(FSUTOApad)
FSDTOApad=new((/nyrs*12,sizefoo(1),sizefoo(2)/),float)
FSDTOApad@_FillValue=-999
FSDTOApad=FSDTOApad@_FillValue
FSDTOApad(2:sizefoo(0)-1+2,:,:)=FSDTOA
FSDTOApad&time=ispan(0,nyrs*12-1,1)
FSDTOApad&lat@units="degrees_north"
FSDTOApad&lon@units="degrees_east"
FSDTOA_month=clmMonTLL(FSDTOApad)
FLUTOApad=new((/nyrs*12,sizefoo(1),sizefoo(2)/),float)
FLUTOApad@_FillValue=-999
FLUTOApad=FSUTOApad@_FillValue
FLUTOApad(2:sizefoo(0)-1+2,:,:)=FSUTOA
FLUTOApad&time=ispan(0,nyrs*12-1,1)
FLUTOApad&lat@units="degrees_north"
FLUTOApad&lon@units="degrees_east"
FLUTOA_month=clmMonTLL(FSUTOApad)

;;; make a mask
a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
data=FSUTOA_month
lsm  = landsea_mask(lsdata,data&lat,data&lon)
lsm_ocneq1=byte2flt(lsm)
lsm_ocneq1=where(lsm.eq.0,1,0)
;; find number of longitudes covered by ocean at each latitude
cntlat=dim_sum_n_Wrap(lsm_ocneq1,1)

;;; apply mask
FSUTOA_month = mask(FSUTOA_month,(lsm.eq.1),False)
FSDTOA_month = mask(FSDTOA_month,(lsm.eq.1),False)
FLUTOA_month = mask(FLUTOA_month,(lsm.eq.1),False)
FSUTOApad = mask(FSUTOApad,(lsm.eq.1),False)
FSDTOApad = mask(FSDTOApad,(lsm.eq.1),False)
FLUTOApad = mask(FLUTOApad,(lsm.eq.1),False)

;;;; make an Arctic monthly average of the observations

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; **  find climatology **
;; calculate masked zonal mean for all months
FSUTOA_mevol_zmean=dim_avg_n_Wrap(FSUTOA_month,2)
;; calculate the weighted average from 70-90 north for each month (climatology)
FSUTOA_arctic_climo=new((/12/),float)
do j=0,11  ;;; over month
  sumnum=0.0
  sumden=0.0
  do i=160,179 ;;; over 70-90 N lat
    sumnum=sumnum + FSUTOA_mevol_zmean(j,i)*(cntlat(i)/360.)*wgty(i)
    sumden=sumden + (cntlat(i)/360.)*wgty(i)
  end do
  FSUTOA_arctic_climo(j)=sumnum/sumden
end do
;print(FSUTOA_arctic_climo)
;;; ** find timeseries ** 
;; calculate masked zonal mean for all times
FSUTOA_mevol_tevol_zmean=dim_avg_n_Wrap(FSUTOApad,2)
FSUTOA_arctic_tseries=new((/12*nyrs/),float)
do j=0,dimsizes(FSUTOA_arctic_tseries)-1 ;;; over timeseries
  sumnum=0.0
  sumden=0.0
  do i=160,179 ;;; over 70-90 N lat
    sumnum=sumnum + FSUTOA_mevol_tevol_zmean(j,i)*(cntlat(i)/360.)*wgty(i)
    sumden=sumden + (cntlat(i)/360.)*wgty(i)
  end do
  FSUTOA_arctic_tseries(j)=sumnum/sumden
end do

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; **  find climatology **
;; calculate masked zonal mean for all months
FSDTOA_mevol_zmean=dim_avg_n_Wrap(FSDTOA_month,2)
;; calculate the weighted average from 70-90 north for each month (climatology)
FSDTOA_arctic_climo=new((/12/),float)
do j=0,11  ;;; over month
  sumnum=0.0
  sumden=0.0
  do i=160,179 ;;; over 70-90 N lat
    sumnum=sumnum + FSDTOA_mevol_zmean(j,i)*(cntlat(i)/360.)*wgty(i)
    sumden=sumden + (cntlat(i)/360.)*wgty(i)
  end do
  FSDTOA_arctic_climo(j)=sumnum/sumden
end do
;print(FSDTOA_arctic_climo)
;;; ** find timeseries ** 
;; calculate masked zonal mean for all times
FSDTOA_mevol_tevol_zmean=dim_avg_n_Wrap(FSDTOApad,2)
FSDTOA_arctic_tseries=new((/12*nyrs/),float)
do j=0,dimsizes(FSDTOA_arctic_tseries)-1 ;;; over timeseries
  sumnum=0.0
  sumden=0.0
  do i=160,179 ;;; over 70-90 N lat
    sumnum=sumnum + FSDTOA_mevol_tevol_zmean(j,i)*(cntlat(i)/360.)*wgty(i)
    sumden=sumden + (cntlat(i)/360.)*wgty(i)
  end do
  FSDTOA_arctic_tseries(j)=sumnum/sumden
end do

;;; find timeseries of net TOA shortwave
FSNTOA_arctic=FSDTOA_arctic_tseries-FSUTOA_arctic_tseries
foo=onedtond(FSNTOA_arctic,(/nyrs,12/))
FSNTOA_arctic_jja=new((/4,dimsizes(year)/),float)
FSNTOA_arctic_jja(0,:)=foo(:,5)
FSNTOA_arctic_jja(1,:)=foo(:,6)
FSNTOA_arctic_jja(2,:)=foo(:,7)
FSNTOA_arctic_jja(3,:)=(foo(:,5)+foo(:,6)+foo(:,7))/3
delete(foo)

;;; print data, check it
print(FSNTOA_arctic_jja(3,:))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; **  find climatology **
;; calculate masked zonal mean for all months
FLUTOA_mevol_zmean=dim_avg_n_Wrap(FLUTOA_month,2)
;; calculate the weighted average from 70-90 north for each month (climatology)
FLUTOA_arctic_climo=new((/12/),float)
do j=0,11  ;;; over month
  sumnum=0.0
  sumden=0.0
  do i=160,179 ;;; over 70-90 N lat
    sumnum=sumnum + FLUTOA_mevol_zmean(j,i)*(cntlat(i)/360.)*wgty(i)
    sumden=sumden + (cntlat(i)/360.)*wgty(i)
  end do
  FLUTOA_arctic_climo(j)=sumnum/sumden
end do
;print(FLUTOA_arctic_climo)
;;; ** find timeseries ** 
;; calculate masked zonal mean for all times
FLUTOA_mevol_tevol_zmean=dim_avg_n_Wrap(FLUTOApad,2)
FLUTOA_arctic_tseries=new((/12*nyrs/),float)
do j=0,dimsizes(FLUTOA_arctic_tseries)-1 ;;; over timeseries
  sumnum=0.0
  sumden=0.0
  do i=160,179 ;;; over 70-90 N lat
    sumnum=sumnum + FLUTOA_mevol_tevol_zmean(j,i)*(cntlat(i)/360.)*wgty(i)
    sumden=sumden + (cntlat(i)/360.)*wgty(i)
  end do
  FLUTOA_arctic_tseries(j)=sumnum/sumden
end do


;;; ** find anomaly timeseries of net TOA shortwave**
mean=ndtooned(conform_dims((/nyrs,12/),FSUTOA_arctic_climo,1))
FSNTOA_arctic_manom=-(FSUTOA_arctic_tseries-mean) 
;;; note: negative sign in above equation converts from FSUTOA to FSNTOA.  Assuming that FSDTOA has not changed.
delete(mean)
foo=onedtond(FSNTOA_arctic_manom,(/nyrs,12/))
FSNTOA_arctic_manom_jja=new((/4,dimsizes(year)/),float)
FSNTOA_arctic_manom_jja(0,:)=foo(:,5)
FSNTOA_arctic_manom_jja(1,:)=foo(:,6)
FSNTOA_arctic_manom_jja(2,:)=foo(:,7)
FSNTOA_arctic_manom_jja(3,:)=(foo(:,5)+foo(:,6)+foo(:,7))/3
delete(foo)

;;; print data, check it
print(FSNTOA_arctic_manom_jja(3,:))

;;; find timeseries of incoming (downwelling) TOA shortwave
FSDTOA_arctic=FSDTOA_arctic_tseries
foo=onedtond(FSDTOA_arctic,(/nyrs,12/))
FSDTOA_arctic_jja=new((/4,dimsizes(year)/),float)
FSDTOA_arctic_jja(0,:)=foo(:,5)
FSDTOA_arctic_jja(1,:)=foo(:,6)
FSDTOA_arctic_jja(2,:)=foo(:,7)
FSDTOA_arctic_jja(3,:)=(foo(:,5)+foo(:,6)+foo(:,7))/3
delete(foo)

;;; print data, check it
print(FSDTOA_arctic_jja(3,:))

;;; find timeseries of outgoing longwave radiation
FLUTOA_arctic=FLUTOA_arctic_tseries
foo=onedtond(FLUTOA_arctic,(/nyrs,12/))
FLUTOA_arctic_ann=(foo(:,0)+foo(:,1)+foo(:,2)+foo(:,3)+foo(:,4)+foo(:,5)+foo(:,6)+foo(:,7)+foo(:,8)+foo(:,9)+foo(:,10)+foo(:,11))/12
FLUTOA_arctic_son=(foo(:,8)+foo(:,9)+foo(:,10))/3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; SEA ICE EXTENT/AREA OBSERVATIONS ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   nyrsobs=43+4
   lastyrobs=2021+4
   NHextent= new((/12,nyrsobs/),float)
   NHextent@_FillValue=100
   NHarea= NHextent
   months = (/"01","02","03","04","05","06","07","08","09","10","11","12"/)
   do i=0,11
      NHdata=readAsciiTable("/glade/u/home/jenkay/small_data/NSIDC/N_"+months(i)+"_extent_v3.0.csv",4,"float",1)
      ;printVarSummary(NHdata)
      ;print(NHdata(:,2))
      NHextent(i,:)=NHdata(0:nyrsobs-1,2)
      NHarea(i,:)=NHdata(0:nyrsobs-1,3)
      delete(NHdata)
    end do

   NHextent@long_name = " Arctic Ice Extent"
   NHextent@units = "million km2"
   NHextent!0="month"
   NHextent!1 = "year"
   NHextent&year = ispan(1979,lastyrobs,1)
   NHextent@_FillValue=100
   copy_VarCoords(NHextent,NHarea)
   
   ICE_month_mean=dim_avg_n_Wrap(NHextent(:,{2000:2021}),1)
   
   printVarSummary(NHextent)
   print(NHextent(7,:))
   print(ICE_month_mean)

   ;; assemble data for plotting
   year_all=(/2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025/)
   ICE_arctic_manom_jja=new((/4,dimsizes(year_all)/),float)
   ICE_arctic_manom_jja(0,:)=((NHextent(5,{2000:2025})-ICE_month_mean(5)))
   ICE_arctic_manom_jja(1,:)=((NHextent(6,{2000:2025})-ICE_month_mean(6)))
   ICE_arctic_manom_jja(2,:)=((NHextent(7,{2000:2025})-ICE_month_mean(7)))
   mean=(ICE_month_mean(5)+ICE_month_mean(6)+ICE_month_mean(7))/3
   ICE_arctic_manom_jja(3,:)=(((NHextent(5,{2000:2025})+NHextent(6,{2000:2025})+NHextent(7,{2000:2025}))/3)-mean)

   ICE_arctic_jja=new((/4,dimsizes(year_all)/),float)
   ICE_arctic_jja(0,:)=NHextent(5,{2000:2025})
   ICE_arctic_jja(1,:)=NHextent(6,{2000:2025})
   ICE_arctic_jja(2,:)=NHextent(7,{2000:2025})
   ICE_arctic_jja(3,:)=((NHextent(5,{2000:2025})+NHextent(6,{2000:2025})+NHextent(7,{2000:2025}))/3)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; FIGURE - Time series of JJA absorbed shortwave and sea ice extent

;;; plot the monthly evolution
  plot_name = "timeseries_FSNTOA_seaiceextent_through2024"

  wks = gsn_open_wks(wks_type,plot_name)
  colors = (/"white","black","white","royal blue","light sky blue",\
             "powder blue","light sea green","pale green","wheat","brown",\
             "pink","light blue","firebrick","DodgerBlue","CornflowerBlue",\
             "light steel blue","Cyan","sky blue","red","orange"/)
  gsn_define_colormap(wks, colors)  

  plot = new(2,graphic)
  res  = True
  res@vpHeightF= 0.4                    ; change aspect ratio of plot and slide it over so it fits on 
  res@vpWidthF = 1.5                    ; the screen
  res@vpXF     = 0.2
  res@gsnDraw                  = False
  res@gsnFrame                 = False
  res@tiMainFontHeightF = 0.012
  res@xyLineThicknesses = (/3, 3/)
  res@xyLineColors = (/"black","blue"/)
  res@xyMarkLineMode = (/"MarkLines","MarkLines"/)
  res@xyMarker = (/16,16/)
  res@xyDashPatterns = (/0,1/)
  res@tmXBMode          = "Explicit"                ; explicit labels
  res@trXMinF           = 2000
  res@trXMaxF           = 2025
  res@tmXBValues        = (/2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020,2022,2024/) ; location of labels 
  res@tmXBLabels        = (/2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020,2022,2024/) ; labels themselves
  res@tiYAxisString = "Absorbed solar (Wm-2)"
  res@trYMinF = 200.
  res@trYMaxF = 225.
  plot(0) = gsn_csm_xy(wks,year_all,FSNTOA_arctic_jja(3,:),res)
  res@tiYAxisString = "Sea ice extent (msqkm)"
  res@trYMinF = 7.5
  res@trYMaxF = 9.5
  plot(1) = gsn_csm_xy(wks,year_all,ICE_arctic_jja(3,:),res)

  resP                  = True                   ; modify the panel plot
  resP@txString         = "Arctic Ocean (70-90 N) Summer (JJA) Averages"  ; add common title
  resP@gsnMaximize      = True                   ; use full page
  gsn_panel(wks,plot,(/2,1/),resP)

delete(plot)
delete(res)


;;;; WRITE OUT A CSV FILE WITH THE DATA

print(year)
print(ICE_arctic_jja(3,:))
print(FSNTOA_arctic_jja(3,:))
print(FSDTOA_arctic_jja(3,:))


print(year(0:24))
print(FSDTOA_arctic_jja(3,0:24))


 alist = [/year(0:24),ICE_arctic_jja(3,0:24),FSNTOA_arctic_jja(3,0:24),FSDTOA_arctic_jja(3,0:24)/]

 header = (/"------------------------------------------------------------------------------------------------", \
            "Arctic Ocean Summer (JJA) Observations 2000-2024", \
	    "Column #1 = Sea Ice Extent, units = million square kilometers", \
	    "Source: NSIDC Sea Ice Index https://noaadata.apps.nsidc.org/NOAA/G02135/, downloaded May 2025", \
	    "Column #2 = Top-of-atmosphere absorbed shortwave radiation, units = Watts per meter squared", \
	    "Source: CERES-EBAF Ed4.2.1, downloaded from https://ceres-tool.larc.nasa.gov/ord-tool/srbavg May 2025", \
	    "Column #3 = Top-of-atmosphere incoming shortwave radiation, units = Watts per meter squared", \
	    "Source: CERES-EBAF Ed4.2.1, downloaded from https://ceres-tool.larc.nasa.gov/ord-tool/srbavg May 2025", \
	    "Data munged by Jennifer.E.Kay@colorado.edu using plot_anom_tseries_include2025.ncl, last update October 2, 2025", \
            "------------------------------------------------------------------------------------------------"/)

 hlist = [/header/]

 fname = "ArcticOceanData_2000-2024.txt"
 write_table(fname, "w", hlist, "%s ")
 write_table(fname, "a", alist, "%d,%16.2f,%16.2f,%16.2f")




print("stopping here")
exit

;############################################
;############################################
;############################################
;############################################


print("JJA TOA absorbed sw trends 2000-2024")

rc   = regline(fspan(0,21,22),FSNTOA_arctic_manom_jja(3,:))
df   = rc@nptxy-2
prob = (1 - betainc(df/(df+rc@tval^2), df/2.0, 0.5) )
print("slope: "+rc)
print("prob: "+prob)
N    = rc@nptxy                        ; convenience/clarity
acr  = esacr(FSNTOA_arctic_manom_jja(3,:),2)
if (acr(1).gt.0.0) then
   pr1     = rtest(acr(1), N, 0)
   rsiglvl = 0.10
   if (pr1.lt.rsiglvl) then
      dfa   = N*(1.0-acr(1))/(1.0+acr(1))
      proba = (1 - betainc(dfa/(dfa+rc@tval^2), dfa/2.0, 0.5) )
      print("prob (remove autocorrelation): "+proba)
   end if
end if

print("")
print("")

print("JJA ice extent trends trends 2000-2021")
rc   = regline(fspan(0,21,22),ICE_arctic_manom_jja(3,:))
df   = rc@nptxy-2
prob = (1 - betainc(df/(df+rc@tval^2), df/2.0, 0.5) )
print("slope: "+rc)
print("prob: "+prob)
N    = rc@nptxy                        ; convenience/clarity
acr  = esacr(ICE_arctic_manom_jja(3,:),2)
if (acr(1).gt.0.0) then
   pr1     = rtest(acr(1), N, 0)
   rsiglvl = 0.10
   if (pr1.lt.rsiglvl) then
      dfa   = N*(1.0-acr(1))/(1.0+acr(1))
      proba = (1 - betainc(dfa/(dfa+rc@tval^2), dfa/2.0, 0.5) )
      print("prob (remove autocorrelation): "+proba)
   end if
end if

print("")
print("")

printVarSummary(year)
printVarSummary(FSNTOA_arctic_jja)

end
