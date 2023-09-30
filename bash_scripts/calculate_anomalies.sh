#!/bin/bash

month='Aug' #NOTE: CHANGE CDO SELMON #
CRU="/discover/nobackup/projects/SBG-DO/data/Climate/CRU/CRU407/"
tmp="$CRU/cru_ts4.07.1901.2022.tmp.dat.nc"
pre="$CRU/cru_ts4.07.1901.2022.pre.dat.nc"
out="/discover/nobackup/bcurrey/Currey_IS_ClimateExtremes_OSSE/data/"

echo -e "Select Month"
cdo selmon,8 $tmp $out/${month}_tmp.nc
cdo selmon,8 $pre $out/${month}_pre.nc

#select period and calculate the average through time
cdo timmean -selyear,1970/2020 $out/${month}_tmp.nc $out/${month}_tmp_long_term_avg_1970_2020.nc
cdo timmean -selyear,1970/2020 $out/${month}_pre.nc $out/${month}_pre_long_term_avg_1970_2020.nc

# select 2022 and subtract that from the LTA data.
cdo sub -selyear,2022 $out/${month}_tmp.nc  $out/${month}_tmp_long_term_avg_1970_2020.nc $out/tmp_anomaly_${month}_2022.nc
cdo sub -selyear,2022 $out/${month}_pre.nc $out/${month}_pre_long_term_avg_1970_2020.nc $out/pre_anomaly_${month}_2022.nc

rm $out/${month}_tmp.nc
rm $out/${month}_pre.nc

