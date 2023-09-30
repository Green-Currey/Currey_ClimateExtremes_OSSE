source('~/R/clean.R')
library(ncdf4)
library(readr)
library(dplyr)
library(terra)
source('~/Current Projects/SBG/LPJ/Global_trait_PLSRs/R_scripts/lpj_plsr_functions.R')

fp <- '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/'
dp <- file.path(fp, 'data/')
ncdfs <- list.files(dp, pattern = '*.nc', full.names = T); ncdfs
coeffs <- read_csv(file.path(dp, 'LPJ-PROSAIL_TRY_coeffs_July.csv'))



# Baseline 2022 ----------------------------------------------------------------

r <- nc_open(ncdfs[4]) %>%
    ncvar_get('DR', start = c(1,1,1,7), count = c(-1,-1,-1,1)) %>%
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326')
r[r>0.7] <- NA # remove snow
r[r==0] <- NA # remove zeros

lma.base <- trait.map(r, coeffs = coeffs$lma[-1], intercept = coeffs$lma[1], coeffs_wl = seq(400,2500,10))    


# 2deg --------------------------------------------------------------------

r <- nc_open(ncdfs[3]) %>%
    ncvar_get('DR', start = c(1,1,1,7), count = c(-1,-1,-1,1)) %>%
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326')
r[r>0.7] <- NA # remove snow
r[r==0] <- NA # remove zeros

lma.2deg <- trait.map(r, coeffs = coeffs$lma[-1], intercept = coeffs$lma[1], coeffs_wl = seq(400,2500,10))    



library(ggplot2)
library(tidyterra)
library(ggpubr)
p1 <- ggplot() +
    geom_spatraster(data = lma.base) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), limits = c(0.25, 0.4), na.value = 'transparent') +
    theme_void(base_size = 20) +
    labs(title = 'Baseline LMA')

p2 <- ggplot() +
    geom_spatraster(data = lma.2deg) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), limits = c(0.25, 0.4), na.value = 'transparent') +
    theme_void(base_size = 20) +
    labs(title = '+2Degrees LMA')

lma.dif <- lma.2deg - lma.base
p3 <- ggplot() +
    geom_spatraster(data = lma.dif) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), limits = c(0.25, 0.4), na.value = 'transparent') +
    theme_void(base_size = 20) +
    labs(title = '+2 - Baseline')


lma.pdif <- ((lma.2deg - lma.base)/lma.base)*100
p4 <- ggplot() +
    geom_spatraster(data = lma.pdif) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), limits = c(0.25, 0.4), na.value = 'transparent') +
    theme_void(base_size = 20) +
    labs(title = '% Change')