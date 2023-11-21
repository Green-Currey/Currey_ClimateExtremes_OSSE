#!/bin/bash

month='Aug' #NOTE: CHANGE CDO SELMON #
CRU="/discover/nobackup/projects/SBG-DO/data/Climate/CRU/CRU407/"
tmp="$CRU/cru_ts4.07.1901.2022.tmp.dat.nc"
pre="$CRU/cru_ts4.07.1901.2022.pre.dat.nc"
out="/discover/nobackup/bcurrey/Currey_IS_ClimateExtremes_OSSE/data/"
ref1=1970
ref2=2010
year=2012

echo -e "Select Month"
cdo selmon,8 $tmp $out/${month}_tmp.nc
cdo selmon,8 $pre $out/${month}_pre.nc

#select period and calculate the average through time
cdo timmean -selyear,$ref1/$ref2 $out/${month}_tmp.nc $out/${month}_tmp_long_term_avg_${ref1}_${ref2}.nc
cdo timmean -selyear,$ref1/$ref2 $out/${month}_pre.nc $out/${month}_pre_long_term_avg_${ref1}_${ref2}.nc

# select 2022 and subtract that from the LTA data.
cdo sub -selyear,$year $out/${month}_tmp.nc  $out/${month}_tmp_long_term_avg_${ref1}_${ref2}.nc $out/tmp_anomaly_${month}_${year}.nc
cdo sub -selyear,$year $out/${month}_pre.nc $out/${month}_pre_long_term_avg_${ref1}_${ref2}.nc $out/pre_anomaly_${month}_${year}.nc

rm $out/${month}_tmp.nc
rm $out/${month}_pre.nc

