#!/bin/bash


mkdir temp
cdo seltimestep,1/1452 cru_ts4.07.1901.2022.tmp.dat.nc temp/first1452.nc
cdo seltimestep,1453/1464 cru_ts4.07.1901.2022.tmp.dat.nc temp/last12.nc
cdo addc,2 temp/last12.nc temp/last12_force2.nc
cdo cat temp/first1452.nc temp/last12_force2.nc cru_2deg.1901.2022.tmp.dat.nc
cdo seltimestep,-12/-1 cru_2deg.1901.2022.tmp.dat.nc cru_2deg.2022.tmp.dat.nc
cdo seltimestep,-12/-1 cru_ts4.07.1901.2022.tmp.dat.nc cru_ts4.07.2022.tmp.dat.nc

cdo fldmean -timmean cru_2deg.2022.tmp.dat.nc temp/mean2deg.nc
cdo fldmean -timmean cru_ts4.07.2022.tmp.dat.nc temp/mean.nc

mean=$(ncdump temp/mean.nc | awk '/tmp =/{getline; print $0}')
mean2deg=$(ncdump temp/mean2deg.nc | awk '/tmp =/{getline; print $0}')
echo -e "mean: $mean"
echo -e "mean 2 deg: $mean2deg"

rm temp -rf

