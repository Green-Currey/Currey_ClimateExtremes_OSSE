library(dplyr)
library(ncdf4)
library(readr)
library(ggpubr)
source('~/R/figure.R')

# sbg <- read_delim('~/Current Projects/SBG/LPJ/Hypertrace-LPJ-PROSAIL/data/sbg_noise_coeffs.txt', delim = '  ', col_names = c('wl', 'a', 'b', 'c', 'rmse'))[-c(1,2),]
sbg <- read_delim('~/Current Projects/SBG/LPJ/Hypertrace-LPJ-PROSAIL/data/20231107_sbg_noise.txt', delim = ' ', col_names = c('wl', 'a', 'b', 'c', 'rmse'))
wl <- seq(400,2500,10)

# this cell has low error with retrievals
L <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelD_TOA-radiance_Version021_m_2015.nc') %>%
    ncvar_get('Radiance', start = c(1,242,252,8), count = c(-1,1,1,1)) 
R <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelC_DR_Version021_m_2015.nc') %>% # lon, lat, wl, time
    ncvar_get('DR', start = c(252,242,1,8), count = c(1,1,-1,1))

# this cell has high error with retrievals
# L <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelD_TOA-radiance_Version021_m_2015.nc') %>%
#     ncvar_get('Radiance', start = c(1,166,374,8), count = c(-1,1,1,1)) 
# R <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelC_DR_Version021_m_2015.nc') %>% # lon, lat, wl, time
#     ncvar_get('DR', start = c(374,166,1,8), count = c(1,1,-1,1))
noise <- data.frame(wl, L)

noise$a <- approx(x = sbg$wl, y = sbg$a, xout = wl)$y
noise$b <- approx(x = sbg$wl, y = sbg$b, xout = wl)$y
noise$c <- approx(x = sbg$wl, y = sbg$c, xout = wl)$y
noise$rmse <- approx(x = sbg$wl, y = sbg$rmse, xout = wl)$y

# sbg2 <- noise %>% select(wl, a, b, c, rmse)
# write_delim(sbg2,'~/Current Projects/SBG/LPJ/Hypertrace-LPJ-PROSAIL/data/20231107_sbg_noise_10nm.txt', delim = ' ')

# Noise equivalent detection limit
# NeDL = a * sqrt(b + L) + c
# 
# Where:
#     
#     SNR = L / NeDL
#  or
#     NeDL = L/SNR

noise$NeDL = noise$a * sqrt(noise$b + noise$L) + noise$c
noise$wl[noise$NeDL>0.05]
# wl to exclude:
bad.wl <- c(400,1350,1370,1380,1390,1400,1830,1840,1850,1860,1870,1880,1890,1900,1910)
noise$NeDL[noise$wl %in% bad.wl] <- NA
noise$SNR = noise$L / noise$NeDL

noise %>% filter(wl < 700) %>% summarise(mean(SNR, na.rm = T)) #311 for < 700
noise %>% filter(wl >= 700) %>% summarise(mean(SNR, na.rm = T)) # 410 < 700


# Panel in Figure 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(
    ggplot(noise) +
        geom_line(aes(wl,SNR), linewidth = 1) +
        theme_pubclean(base_size = 40) +
        labs(x = 'Wavelength (nm)', y = 'SNR') +
        ggtitle('SBG Instrument Model') +
        theme(plot.title = element_text(hjust = 1)),
    path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/Instrument_SNR_plot.png', 
    save = T)



# Misc other plots --------------------------------------------------------


# noise <- read_delim('~/Current Projects/SBG/LPJ/Hypertrace-LPJ-PROSAIL/data/sbg_noise_coeffs.txt', delim = '  ')
sbg <- read_delim('~/Current Projects/SBG/LPJ/Hypertrace-LPJ-PROSAIL/data/20231107_sbg_noise.txt', delim = ' ', col_names = c('wl', 'a', 'b', 'c', 'rmse'))
wl <- seq(400,2500,10)

# this cell has low error with retrievals
Ll <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelD_TOA-radiance_Version021_m_2015.nc') %>%
    ncvar_get('Radiance', start = c(1,242,252,8), count = c(-1,1,1,1)) 
Rl <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelC_DR_Version021_m_2015.nc') %>% # lon, lat, wl, time
    ncvar_get('DR', start = c(252,242,1,8), count = c(1,1,-1,1))

# this cell has high error with retrievals
Lh <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelD_TOA-radiance_Version021_m_2015.nc') %>%
    ncvar_get('Radiance', start = c(1,166,374,8), count = c(-1,1,1,1))
Rh <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelC_DR_Version021_m_2015.nc') %>% # lon, lat, wl, time
    ncvar_get('DR', start = c(374,166,1,8), count = c(1,1,-1,1))
noise <- data.frame(wl, Ll, Rl, Lh, Rh)

noise$a <- approx(x = sbg$wl, y = sbg$a, xout = wl)$y
noise$b <- approx(x = sbg$wl, y = sbg$b, xout = wl)$y
noise$c <- approx(x = sbg$wl, y = sbg$c, xout = wl)$y
noise$rmse <- approx(x = sbg$wl, y = sbg$rmse, xout = wl)$y


noise$NeDLl = noise$a * sqrt(noise$b + noise$Ll) + noise$c
noise$NeDLh = noise$a * sqrt(noise$b + noise$Lh) + noise$c
noise$NeDLl[noise$NeDLl > 0.05] <- NA
noise$NeDLh[noise$NeDLh > 0.05] <- NA

noise$SNRl = noise$Ll / noise$NeDLl
noise$SNRh = noise$Lh / noise$NeDLh


top1 <- ggplot(noise) +
    geom_line(aes(wl,Rl), linewidth = 1) +
    theme_pubclean(base_size = 25) +
    ggtitle('Grassland') +
    labs(x = '', y = bquote('Reflectance'))

mid1 <- ggplot(noise) +
    geom_line(aes(wl,Ll), linewidth = 1) +
    theme_pubclean(base_size = 25) +
    labs(x = '', y = bquote('Radiance (W'~sr^-1~m^-2*')'))

bottom1 <- ggplot(noise) +
    geom_line(aes(wl,SNRl), linewidth = 1) +
    theme_pubclean(base_size = 25) +
    labs(x = 'Wavelength (nm)', y = 'SNR')
    

top2 <- ggplot(noise) +
    geom_line(aes(wl,Rh), linewidth = 1) +
    theme_pubclean(base_size = 25) +
    ggtitle('Tropical Forest') +
    labs(x = '', y = bquote('Reflectance'))

mid2 <- ggplot(noise) +
    geom_line(aes(wl,Lh), linewidth = 1) +
    theme_pubclean(base_size = 25) +
    labs(x = '', y = bquote('Radiance (W'~sr^-1~m^-2*')'))

bottom2 <- ggplot(noise) +
    geom_line(aes(wl,SNRh), linewidth = 1) +
    theme_pubclean(base_size = 25) +
    labs(x = 'Wavelength (nm)', y = 'SNR')

ggarrange(top1, top2,
          mid1, mid2,
          bottom1, bottom2,
          nrow = 3, ncol = 2, align = 'v')

