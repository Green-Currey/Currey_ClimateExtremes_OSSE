#!/bin/bash


CRU="/discover/nobackup/projects/SBG-DO/data/Climate/CRU/CRU407/"
tmp="$CRU/cru_ts4.07.1901.2022.tmp.dat.nc"
pre="$CRU/cru_ts4.07.1901.2022.pre.dat.nc"
out="/discover/nobackup/bcurrey/Currey_IS_ClimateExtremes_OSSE/data/"

echo -e "Select Month"
cdo selmon,7 $tmp $out/july_tmp.nc
cdo selmon,7 $pre $out/july_pre.nc

#select period and calculate the average through time
cdo timmean -selyear,1970/2020 $out/july_tmp.nc $out/tmp_long_term_avg_1970_2020.nc
cdo timmean -selyear,1970/2020 $out/july_pre.nc $out/pre_long_term_avg_1970_2020.nc

# select 2022 and subtract that from the LTA data.
cdo sub -selyear,2022 $out/july_tmp.nc  $out/tmp_long_term_avg_1970_2020.nc $out/tmp_anomaly_july_2022.nc
cdo sub -selyear,2022 $out/july_pre.nc $out/pre_long_term_avg_1970_2020.nc $out/pre_anomaly_july_2022.nc

rm $out/july_tmp.nc
rm $out/july_pre.nc

