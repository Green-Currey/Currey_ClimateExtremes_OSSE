library(ncdf4)
library(terra)
library(tidyverse)
library(ggplot2)
library(tidyterra)
library(ggpubr)
source('~/Current Projects/SBG/LPJ/Global_trait_PLSRs/R_scripts/lpj_plsr_functions.R')


# set data path
dp <- '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/'


# Temp anomalies ----------------------------------------------------------

# read in temp anaomalies nc and quickview
tmp.anomalies <- nc_open(file.path(dp, 'tmp_anomaly_Aug_2018.nc')) %>%
    ncvar_get('tmp') %>%
    aperm(c(2,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    flip
plot(tmp.anomalies)
global(tmp.anomalies, mean, na.rm = T)
# Global average temp anomaly = 0.79 degrees

# Create region 1 box
xmin <- -125
xmax <- -105
ymin <- 35
ymax <- 55
r1 <- matrix(c(xmin, ymin, 
               xmax, ymin, 
               xmax, ymax, 
               xmin, ymax, 
               xmin, ymin), 
             ncol=2, byrow=TRUE) %>%
    vect(type="polygon")
plot(r1, add = T)
global(tmp.anomalies %>% crop(r1), mean, na.rm = T)
# R1 average temp anomaly = 2.34 degrees


# Create region 2 box
xmin <- 30
xmax <- 50
ymin <- 45
ymax <- 65
r2 <- matrix(c(xmin, ymin, 
               xmax, ymin, 
               xmax, ymax, 
               xmin, ymax, 
               xmin, ymin), 
             ncol=2, byrow=TRUE) %>%
    vect(type="polygon")
plot(r2, add = T)
global(tmp.anomalies %>% crop(r2), mean, na.rm = T)
# R2 average temp anomaly = 3.86 degrees


# Precip anomalies --------------------------------------------------------


# read in temp anaomalies nc and quickview
pre.anomalies <- nc_open(file.path(dp, 'pre_anomaly_Aug_2020.nc')) %>%
    ncvar_get('pre') %>%
    aperm(c(2,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    flip
pre.anomalies[pre.anomalies>400] <- NA
plot(pre.anomalies)
global(pre.anomalies, mean, na.rm = T)
# Global average temp anomaly = 0.79 degrees

# Create region 1 box
xmin <- -125
xmax <- -105
ymin <- 35
ymax <- 55
r3 <- matrix(c(xmin, ymin, 
               xmax, ymin, 
               xmax, ymax, 
               xmin, ymax, 
               xmin, ymin), 
             ncol=2, byrow=TRUE) %>%
    vect(type="polygon")
plot(r1, add = T)
global(tmp.anomalies %>% crop(r1), mean, na.rm = T)
# R1 average temp anomaly = 2.34 degrees


# Create region 2 box
xmin <- 30
xmax <- 50
ymin <- 45
ymax <- 65
r4 <- matrix(c(xmin, ymin, 
               xmax, ymin, 
               xmax, ymax, 
               xmin, ymax, 
               xmin, ymin), 
             ncol=2, byrow=TRUE) %>%
    vect(type="polygon")
plot(r2, add = T)
global(tmp.anomalies %>% crop(r2), mean, na.rm = T)

# Apply coeffs to LPJ data ------------------------------------------------

coeffs <- read_csv(file.path(dp, 'LPJ-PROSAIL_PLSR_coefficients_August2022.csv'))
coeffs2 <- read_csv('~/Current Projects/SBG/LPJ/Global_trait_PLSRs/LPJ-PROSAIL_PLSR_coefficients_August2022.csv')
head(coeffs2)
rfl16 <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelE_retrieved-HDR_Version021_m_2016.nc') %>%
    ncvar_get('HDR', start = c(1,1,1,8), count = c(211,-1,-1,1)) %>% 
    aperm(c(2,3,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326')
rfl16[rfl16>0.7] <- NA # remove snow
rfl16[rfl16==0] <- NA # remove zeros
plot(rfl16$lyr.1)

ldmc16 <- trait.map(rfl16, coeffs = coeffs$lma[-1], intercept = coeffs$lma[1], coeffs_wl = seq(400,2500,10))    
n16 <- trait.map(rfl16, coeffs = coeffs$n[-1], intercept = coeffs$n[1], coeffs_wl = seq(400,2500,10))    
p16 <- trait.map(rfl16, coeffs = coeffs$p[-1], intercept = coeffs$p[1], coeffs_wl = seq(400,2500,10))    
sla16 <- trait.map(rfl16, coeffs = coeffs$sla[-1], intercept = coeffs$sla[1], coeffs_wl = seq(400,2500,10))    

global(ldmc16, mean, na.rm = T)
global(ldmc16, range, na.rm = T)

global(n16, mean, na.rm = T)
global(n16, range, na.rm = T)

global(p16, mean, na.rm = T)
global(p16, range, na.rm = T)

global(sla16, mean, na.rm = T)
global(sla16, range, na.rm = T)


ldmc.r1 <- crop(ldmc16, r1)
ldmc.r2 <- crop(ldmc16, r1)



# 2020 --------------------------------------------------------------------


rfl20 <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelC_DR_Version021_m_2020.nc') %>%
    ncvar_get('DR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>% # select the month of august
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326')
rfl20[rfl20>0.7] <- NA # remove snow
rfl20[rfl20==0] <- NA # remove zeros
plot(rfl20$lyr.1)

lma20 <- trait.map(rfl20, coeffs = coeffs$lma[-1], intercept = coeffs$lma[1], coeffs_wl = seq(400,2500,10))    
n20 <- trait.map(rfl20, coeffs = coeffs$n[-1], intercept = coeffs$n[1], coeffs_wl = seq(400,2500,10))    
p20 <- trait.map(rfl20, coeffs = coeffs$p[-1], intercept = coeffs$p[1], coeffs_wl = seq(400,2500,10))    
# sla20 <- trait.map(rfl20, coeffs = coeffs$sla[-1], intercept = coeffs$sla[1], coeffs_wl = seq(400,2500,10))    

# rfl22 <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelC_DR_Version021_m_2020.nc') %>%
#     ncvar_get('DR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>%
#     aperm(c(2,1,3)) %>%
#     rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326')



# 2022 --------------------------------------------------------------------





