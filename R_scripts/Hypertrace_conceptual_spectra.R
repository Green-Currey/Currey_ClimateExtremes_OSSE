source('~/R/clean.R')
library(ncdf4)
library(ggpubr)
library(tidyverse)
library(ggplot2)
source('~/R/figure.R')

dp <- '~/Current Projects/SBG/LPJ/Reflectance_Data/version2/'


# line plots --------------------------------------------------------------



rfl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2015.nc')) %>%
    ncvar_get('DR', start = c(460,80,1,8), count = c(1,1,-1,1)) 

rdn <- nc_open(file.path(dp, 'lpj-prosail_levelD_TOA-radiance_Version021_m_2015.nc')) %>%
    ncvar_get('Radiance', start = c(1,80,460,8), count = c(-1,1,1,1))

rtr <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2015.nc')) %>%
    ncvar_get('HDR', start = c(1,80,460,8), count = c(-1,1,1,1))

df <- data.frame(wl = seq(400,2500,10), rfl = rfl, rdn = rdn, rtr = rtr)


ggplot(df) +
    geom_line(aes(wl, rfl), linewidth = 2) +
    geom_line(aes(wl, rtr), linewidth = 3, color = 'red', linetype = 'longdash') +
    theme_void()

ggplot(df) +
    geom_line(aes(wl, rdn), linewidth = 2, color = 'grey50') +
    theme_void()

# spectra plots -----------------------------------------------------------

# High error
error <- 'High'
rfl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2015.nc')) %>% # lon, lat, wl, time
    ncvar_get('DR', start = c(374,166,1,8), count = c(3,3,-1,1)) %>%
    aperm(c(3,2,1)) %>%
    as.data.frame() %>%
    mutate(wl = seq(400,2500,10)) %>%
    gather(key = "Cell", value = "Reflectance", -wl)

rtr <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2015.nc')) %>% # wl, lat, lon, time
    ncvar_get('HDR', start = c(1,166,374,8), count = c(-1,3,3,1)) %>%
    as.data.frame() %>%
    mutate(wl = seq(400,2500,10)) %>%
    gather(key = "Cell", value = "Reflectance", -wl)

# Low error
error <- 'Low'
rfl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2015.nc')) %>% # lon, lat, wl, time
    ncvar_get('DR', start = c(252,242,1,8), count = c(3,3,-1,1)) %>%
    aperm(c(3,2,1)) %>%
    as.data.frame() %>%
    mutate(wl = seq(400,2500,10)) %>%
    gather(key = "Cell", value = "Reflectance", -wl) 

rtr <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2015.nc')) %>% # wl, lat, lon, time
    ncvar_get('HDR', start = c(1,242,252,8), count = c(-1,3,3,1)) %>%
    as.data.frame() %>%
    mutate(wl = seq(400,2500,10)) %>%
    gather(key = "Cell", value = "Reflectance", -wl) 



# Differnece
dif <- data.frame(wl = rfl$wl, Cell = rtr$Cell, rfl = rfl$Reflectance, rtr = rtr$Reflectance) %>%
    mutate(difference = rtr - rfl) %>%
    group_by(wl) %>%
    summarise(mean = mean(difference, na.rm = T), sd = sd(difference, na.rm = T))

# Plot
top <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 1, color = 'grey60') +
    # geom_smooth(data = dif, aes(x = wl, y = mean), color = 'black')
    geom_line(data = dif, aes(x = wl, y = mean), color = 'black') +
    # geom_line(data = dif, aes(x = wl, y = difference, group = Cell), color = 'black') +
    geom_ribbon(data = dif, aes(x = wl, ymin = mean-sd, ymax = mean+sd), alpha = 0.3) +
    theme_pubclean(base_size = 25) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.title.x = element_blank()) +
    # coord_cartesian(ylim = c(-0., 0.1)) +
    scale_y_continuous(breaks = c(-0.05,0,0.05)) +
    labs(x = 'Wavelength', y = bquote(~Delta)); top

bottom <- ggplot() +
    geom_line(data = rfl, aes(x = wl, y = Reflectance, group = Cell), color = 'red', linewidth = .5) +
    geom_line(data = rtr, aes(x = wl, y = Reflectance, group = Cell), color = 'blue', linetype = 'longdash', linewidth = .5) +
    theme_pubr(base_size = 20) +
    coord_cartesian(ylim = c(0, 0.5))+
    scale_y_continuous(breaks = c(0.1, 0.3, 0.5))+
    labs(x = 'Wavelength'); bottom



figure(
    ggarrange(top, bottom, 
              nrow = 2, 
              heights = c(0.75,1),
              align = 'v'),
    path.name = paste0('~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/',error,'_RMSE_spectra.png'),
    height = 5, width = 8, save = T
)

