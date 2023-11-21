source('~/R/clean.R')
source('~/R/Figure.R')
library(ncdf4)
library(terra)
library(tidyverse)
library(ggplot2)

shp <- vect('~/Current Projects/SBG/LPJ/Misc data/Shapefile/Would_boundaries_noGreenland.shp')

rfl <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelC_DR_Version021_m_2015.nc') %>%
    ncvar_get('DR', start = c(1,1,2,8), count = c(-1,-1,210,1)) %>%
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) 
rfl[rfl==0] <- NA

rtr <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelE_retrieved-HDR_Version021_m_2015.nc') %>%
    ncvar_get('HDR', start = c(2,1,1,8), count = c(210,-1,-1,1)) %>% 
    aperm(c(2,3,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) 
rtr[rtr==0] <- NA

dif <- rtr-rfl
pdif <- abs(dif)/rfl*100;
pdif[is.infinite(pdif)] <- NA; #pdif
pdif[pdif>50] <- NA; #pdif
pdif[pdif< -50] <- NA; #pdif
mean_pdif <- app(pdif, 'mean', na.rm = T)

AOT <- nc_open('~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/MERRA2-TOTSCATAU_2015.nc') %>%
    ncvar_get('TOTSCATAU', start = c(1,1,8), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp)

H2O <- nc_open('~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/MERRA2-TQV_2015.nc') %>%
    ncvar_get('TQV', start = c(1,1,8), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp)

df <- data.frame(as.data.frame(mean_pdif, na.rm = F),
                 as.data.frame(AOT, na.rm = F),
                 as.data.frame(H2O, na.rm = F)) %>%
    na.exclude
names(df) <- c('% MAE', 'AOT', 'H2O')

summary(lm(`% MAE` ~ AOT, df))
p1 <- ggplot(df) +
    aes(AOT, `% MAE`) +
    geom_point(alpha = 0.05, size = 1) +
    geom_smooth(method = 'lm', level = 0.9999, size = 1.5) +
    theme_pubr(base_size = 25) +
    coord_cartesian(ylim = c(NA, 10), xlim = c(NA, 1))+
    labs(x = bquote('Aerosol Optical Thickness @ 550 nm'), y = 'Global Average % Error'); p1

summary(lm(`% MAE` ~ H2O, df))
p2 <- ggplot(df) +
    aes(H2O/10, `% MAE`) +
    geom_point(alpha = 0.05, size = 1.5) +
    geom_smooth(method = 'lm', level = 0.9999, size = 1.5) +
    theme_pubr(base_size = 25)+
    coord_cartesian(ylim = c(NA, 10), xlim = c(NA, 6))+
    labs(x = bquote('Total Precipitable Water Vapor (g'~cm^-2*')'), y = 'Global Average % Error'); p2

figure(p1,
       path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/AOT_error.png',
       width = 10,
       save = T)

figure(p2,
       path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/H2O_error.png',
       width = 10,
       save = T)



# Exploring radiances -----------------------------------------------------

rdn <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelD_TOA-radiance_Version021_m_2015.nc') %>%
    ncvar_get('Radiance', start = c(2,1,1,8), count = c(210,-1,-1,1)) %>% 
    aperm(c(2,3,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) 
rdn[rdn==0] <- NA

h2o.df <- (H2O/10) %>% as.data.frame(xy = T)
h2o.df[h2o.df$lyr.1==max(h2o.df$lyr.1),]
h2o.df[h2o.df$lyr.1==h2o.df$lyr.1[10000],]
rdn.df <- rdn %>% as.data.frame(xy = T)
df <- data.frame(wl = c(1,2,seq(410,2500,10)), 
                 rdn_highH2O = t(rdn.df %>% filter(x == -78.25, y == -6.75)),
                 rdn_lowH2O = t(rdn.df %>% filter(x == 118.25 , y == 64.75 ))) %>%
    filter(wl>3)

ggplot(df) +
    geom_line(aes(wl, rdn_highH2O), col = 'red3', linewidth = 1) +
    geom_line(aes(wl, rdn_lowH2O), col = 'blue3', linewidth = 1) +
    theme_pubr(base_size = 25) +
    annotate('text', x = 2000, y = 12, label = 'High H2O', color = 'red3', size = 8) +
    annotate('text', x = 2000, y = 10, label = 'Low H2O', color = 'blue3', size = 8) +
    labs(y = 'Radiance')


