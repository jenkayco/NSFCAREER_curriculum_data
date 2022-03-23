# NSFCAREER_curriculum_data
# Last update: March 23, 2022

Short Description: This repo contains processing scripts and resulting timeseries for the Arctic Feedbacks curriculum (https://cires.colorado.edu/outreach/resources/unit/arctic-feedbacks-not-all-warming-equal) developed with Kay NSF CAREER award

Code to munge the data:
1) plot_anom_tseries_include2021.ncl -- produces Arctic Ocean Data for summer (JJA) sea ice extent and top-of-atmosphere radiation (2000-2021)
2) GISTEMP_timeseries_include2021.ncl - calculates Arctic and Global temperature anomalies from GISTEMP (1880-2021)

On Albedo calculations from these data:  Absorbed Shortwave Radiation = Net Shortwave Radiation (SWnet) = Downward Shortwave Radiation (SWdown) - Upward Shortwave Radiation (SWup).  Thus, you can calculate SWup from the fields provided: SWup=SWdown-SWnet.  Then - you can calculate albedo as SWup/SWdown.
