source('~/R/clean.R')
source('~/R/Figure.R')
library(ncdf4)
library(terra)
library(tidyverse)
library(tidyterra)
source('~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/R_scripts/global_spectra_comparsion.R')
source('~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/R_scripts/lpj_plsr_functions.R')
`%!in%` = Negate(`%in%`)


dp <- '~/Current Projects/SBG/LPJ/Reflectance_Data/version2/'
shp <- vect('~/Current Projects/SBG/LPJ/Misc data/Shapefile/Would_boundaries_noGreenland.shp')
op <- '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/'
bad.wl <- c(400,1350,1370,1380,1390,1830,1840,1850,1860,1870,1880,1890,1900,1910)
wl <- seq(400,2500,10)
good.wl <- wl[wl %!in% bad.wl]

# 2015 --------------------------------------------------------------------

rfl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2015.nc')) %>%
    ncvar_get('DR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>%
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>% 
    # subset(which(seq(400,2500,10) %in% wl)) %>%
    mask(shp) 
rfl[rfl==0] <- NA

rtr <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2015.nc')) %>%
    ncvar_get('HDR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>% 
    aperm(c(2,3,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>% 
    # subset(which(seq(400,2500,10) %in% wl)) %>%
    mask(shp) 
rtr[rtr==0] <- NA

dif <- rtr-rfl
mae <- abs(dif) # MAE
pmae <- mae/rfl*100; # %MAE
pmae[is.infinite(pmae)] <- NA; #pmae
pmae[pmae>50] <- NA; #pmae
pmae[pmae< -50] <- NA; #pmae

rfl.sd <- unlist(as.vector(global(rfl, 'sd', na.rm = T))) # Global reflectance SD (211 values)
dev.f <- mae/rfl.sd*100 # deviation fraction

mean_dif <- global(dif, 'mean', na.rm = T); sd_dif <- global(dif, 'sd', na.rm = T)
mean_mae <- global(mae, 'mean', na.rm = T); sd_mae <- global(mae, 'sd', na.rm = T)
mean_pmae <- global(pmae, 'mean', na.rm = T); sd_pmae <- global(pmae, 'sd', na.rm = T)
mean_devf <- global(dev.f, 'mean', na.rm = T); sd_devf <- global(dev.f, 'sd', na.rm = T)

df <- data.frame(wl = wl, 
                 dif = mean_dif$mean, dif_sd = sd_dif$sd, 
                 mae = mean_mae$mean, mae_sd = sd_mae$sd,
                 pmae = mean_pmae$mean, pmae_sd = sd_pmae$sd,
                 df = mean_devf$mean, df_sd = sd_devf$sd)
df[df$wl %in% bad.wl,2:9] <- NA

figure(
    ggplot(df) +
        geom_hline(yintercept = 0) +
        geom_line(aes(wl, dif)) +
        geom_ribbon(aes(x = wl, ymin = dif-dif_sd, ymax = dif+dif_sd), alpha = 0.3) +
        theme_pubr(base_size = 20) +
        ggtitle('2015 Retrievals - Reflectance')+
        labs(x = 'Wavelength (nm)', y = 'Difference'),
    path.name = paste0(op,'2015_error.png'), width = 8, height = 5, save = F
)

figure(
    ggplot(df) +
        geom_hline(yintercept = 0) +
        geom_line(aes(wl, pmae)) +
        geom_ribbon(aes(x = wl, ymin = pmae-pmae_sd, ymax = pmae+pmae_sd), alpha = 0.3) +
        theme_pubr(base_size = 20) +
        ggtitle('% MAE')+
        labs(x = 'Wavelength (nm)', y = '%'),
    path.name = paste0(op,'2015_error.png'), width = 8, height = 5, save = F
)

figure(
    ggplot(df) +
        geom_hline(yintercept = 0) +
        geom_line(aes(wl, df)) +
        geom_ribbon(aes(x = wl, ymin = df-df_sd, ymax = df+df_sd), alpha = 0.3) +
        theme_pubr(base_size = 20) +
        ggtitle('Deviation fraction (%)')+
        labs(x = 'Wavelength (nm)', y = '%'),
    path.name = paste0(op,'2015_error.png'), width = 8, height = 5, save = F
)

mean(df$pmae, na.rm = T); mean(df$pmae_sd, na.rm = T) # 3.2 +/- 3.5 %
mean(df$df, na.rm = T); mean(df$df_sd, na.rm = T) # 3.7 +/- 4.4 %


# global_spectra_comparison(rfl,
#                           rtr,
#                           array1_short_name = 'rfl15',
#                           array2_short_name = 'rtr15',
#                           nc.outpath = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/',
#                           month = 8)


# 2015 no noise ---------------------------------------------
# Reflectance vs retrieved

rfl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2015.nc')) %>%
    ncvar_get('DR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>%
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>% 
    mask(shp) 
rfl[rfl==0] <- NA

rtr <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2015_noNoise.nc')) %>%
    ncvar_get('HDR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>% 
    aperm(c(2,3,1))%>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) 
rtr[rtr==0] <- NA

dif <- rtr-rfl
mae <- abs(dif) # MAE
pmae <- mae/rfl*100; # %MAE
pmae[is.infinite(pmae)] <- NA; #pmae
pmae[pmae>50] <- NA; #pmae
pmae[pmae< -50] <- NA; #pmae

rfl.sd <- unlist(as.vector(global(rfl, 'sd', na.rm = T))) # Global reflectance SD (211 values)
dev.f <- mae/rfl.sd*100 # deviation fraction

mean_dif <- global(dif, 'mean', na.rm = T); sd_dif <- global(dif, 'sd', na.rm = T)
mean_mae <- global(mae, 'mean', na.rm = T); sd_mae <- global(mae, 'sd', na.rm = T)
mean_pmae <- global(pmae, 'mean', na.rm = T); sd_pmae <- global(pmae, 'sd', na.rm = T)
mean_devf <- global(dev.f, 'mean', na.rm = T); sd_devf <- global(dev.f, 'sd', na.rm = T)

df <- data.frame(wl = wl, 
                 dif = mean_dif$mean, dif_sd = sd_dif$sd, 
                 mae = mean_mae$mean, mae_sd = sd_mae$sd,
                 pmae = mean_pmae$mean, pmae_sd = sd_pmae$sd,
                 df = mean_devf$mean, df_sd = sd_devf$sd)
df[df$wl %in% bad.wl,2:9] <- NA

figure(
    ggplot(df) +
        geom_hline(yintercept = 0) +
        geom_line(aes(wl, dif)) +
        geom_ribbon(aes(x = wl, ymin = dif-dif_sd, ymax = dif+dif_sd), alpha = 0.3) +
        theme_pubr(base_size = 20) +
        ggtitle('2015 Retrievals - Reflectance')+
        labs(x = 'Wavelength (nm)', y = 'Difference'),
    path.name = paste0(op,'2015_error.png'), width = 8, height = 5, save = F
)

figure(
    ggplot(df) +
        geom_hline(yintercept = 0) +
        geom_line(aes(wl, df)) +
        geom_ribbon(aes(x = wl, ymin = df-df_sd, ymax = df+df_sd), alpha = 0.3) +
        theme_pubr(base_size = 20) +
        ggtitle('Deviation fraction (%)')+
        labs(x = 'Wavelength (nm)', y = '%'),
    path.name = paste0(op,'2015_error.png'), width = 8, height = 5, save = F
)

mean(df$pmae, na.rm = T); mean(df$pmae_sd, na.rm = T) # 2.8 +/- 3.3 %
mean(df$df, na.rm = T); mean(df$df_sd, na.rm = T) # 3.6 +/- 4.4 %




# rfl.m <- rfl %>% as.data.frame %>% as.matrix
# rtr.m <- rtr %>% as.data.frame %>% as.matrix
# y <- sample(dim(rtr.m)[1], 1000)
# rfl.m <- rfl.m[y,]
# rtr.m <- rtr.m[y,]
# matplot(seq(410,2500,10), t(rfl.m),
#         type = 'l', lty = 1,
#         ylab = "Simulated Reflectance",
#         xlab = "Wavelength (nm)",
#         main = 'Example Spectra subset (n = 1000)')
# matplot(seq(410,2500,10), t(rtr.m),
#         type = 'l', lty = 1,
#         ylab = "Retrieved Reflectance",
#         xlab = "Wavelength (nm)",
#         main = 'No Noise')


# 2016 --------------------------------------------------------------------

rfl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2016.nc')) %>%
    ncvar_get('DR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>%
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp)
rfl[rfl==0] <- NA

rtr <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2016.nc')) %>%
    ncvar_get('HDR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>% 
    aperm(c(2,3,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) 
rtr[rtr==0] <- NA


dif <- rtr-rfl
mae <- abs(dif) # MAE
pmae <- mae/rfl*100; # %MAE
pmae[is.infinite(pmae)] <- NA; #pmae
pmae[pmae>50] <- NA; #pmae
pmae[pmae< -50] <- NA; #pmae

rfl.sd <- unlist(as.vector(global(rfl, 'sd', na.rm = T))) # Global reflectance SD (211 values)
dev.f <- mae/rfl.sd*100 # deviation fraction

mean_dif <- global(dif, 'mean', na.rm = T); sd_dif <- global(dif, 'sd', na.rm = T)
mean_mae <- global(mae, 'mean', na.rm = T); sd_mae <- global(mae, 'sd', na.rm = T)
mean_pmae <- global(pmae, 'mean', na.rm = T); sd_pmae <- global(pmae, 'sd', na.rm = T)
mean_devf <- global(dev.f, 'mean', na.rm = T); sd_devf <- global(dev.f, 'sd', na.rm = T)

df <- data.frame(wl = wl, 
                 dif = mean_dif$mean, dif_sd = sd_dif$sd, 
                 mae = mean_mae$mean, mae_sd = sd_mae$sd,
                 pmae = mean_pmae$mean, pmae_sd = sd_pmae$sd,
                 df = mean_devf$mean, df_sd = sd_devf$sd)
df[df$wl %in% bad.wl,2:9] <- NA


figure(
    ggplot(df) +
        geom_hline(yintercept = 0) +
        geom_line(aes(wl, dif)) +
        geom_ribbon(aes(x = wl, ymin = dif-dif_sd, ymax = dif+dif_sd), alpha = 0.3) +
        theme_pubr(base_size = 20) +
        ggtitle('2016 Retrievals - Reflectance')+
        labs(x = 'Wavelength (nm)', y = 'Difference'),
    path.name = paste0(op,'2016_error.png'), width = 8, height = 5, save = F
)

figure(
    ggplot(df) +
        geom_hline(yintercept = 0) +
        geom_line(aes(wl, df)) +
        geom_ribbon(aes(x = wl, ymin = df-df_sd, ymax = df+df_sd), alpha = 0.3) +
        theme_pubr(base_size = 20) +
        ggtitle('Deviation fraction (%)')+
        labs(x = 'Wavelength (nm)', y = '%'),
    path.name = paste0(op,'2016_error.png'), width = 8, height = 5, save = F
)


mean(df$pmae, na.rm = T); mean(df$pmae_sd, na.rm = T) # 2.8 +/- 3.3 %
mean(df$df, na.rm = T); mean(df$df_sd, na.rm = T) # 3.0 +/- 3.9 %


# global_spectra_comparison(rfl,
#                           rtr,
#                           array1_short_name = 'rfl16',
#                           array2_short_name = 'rtr16',
#                           nc.outpath = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/',
#                           month = 8)

# 2016 no noise ----------

# Reflectance vs retrieved

rfl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2016.nc')) %>%
    ncvar_get('DR', start = c(1,1,2,8), count = c(-1,-1,210,1)) %>%
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) 
rfl[rfl==0] <- NA


rtr <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2016_noNoise.nc')) %>%
    ncvar_get('HDR', start = c(2,1,1,8), count = c(210,-1,-1,1)) %>% 
    aperm(c(2,3,1))%>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) 
rtr[rtr==0] <- NA


dif <- rtr-rfl
pdif <- abs(rtr-rfl)/rfl*100; pdif
pdif[is.infinite(pdif)] <- NA; #pdif
pdif[pdif>50] <- NA; #pdif
pdif[pdif< -50] <- NA; #pdif

mean_dif <- global(dif, 'mean', na.rm = T)
mean_pdif <- global(pdif, 'mean', na.rm = T)
sd_dif <- global(dif, 'sd', na.rm = T)
sd_pdif <- global(pdif, 'sd', na.rm = T)
df2 <- data.frame(wl = seq(410,2500,10), 
                 mean = mean_dif$mean, sd = sd_dif$sd, 
                 pmean = mean_pdif$mean, psd = sd_pdif$sd)

figure(
    ggplot() +
        geom_hline(yintercept = 0) +
        
        # 2016 with noise (has to be run first above)
        geom_line(data = df, aes(wl, mean)) + 
        geom_ribbon(data = df, aes(x = wl, ymin = mean-sd, ymax = mean+sd), alpha = 0.3) +
        
        # 2016 without noise
        geom_ribbon(data = df2, aes(x = wl, ymin = mean-sd, ymax = mean+sd), alpha = 0.3, fill = 'red3') +
        geom_line(data = df2, aes(wl, mean), color = 'red3', linetype = 'longdash') + 
        
        theme_pubr(base_size = 20) +
        ggtitle('2016 Retrievals - Reflectance')+
        labs(x = 'Wavelength (nm)', y = 'Difference'),
    path.name = paste0(op,'2016_error.png'), width = 8, height = 5, save = F
)

mean(df2$pmean); mean(df2$psd)

# dim(rtr)

rfl.m <- rfl %>% as.data.frame %>% as.matrix
rtr.m <- rtr %>% as.data.frame %>% as.matrix
y <- sample(dim(rtr.m)[1], 1000)
rfl.m <- rfl.m[y,]
rtr.m <- rtr.m[y,]
matplot(seq(410,2500,10), t(rfl.m),
        type = 'l', lty = 1,
        ylab = "Simulated Reflectance",
        xlab = "Wavelength (nm)",
        main = 'Example Spectra subset (n = 1000)')
matplot(seq(410,2500,10), t(rtr.m),
        type = 'l', lty = 1,
        ylab = "Retrieved Reflectance",
        xlab = "Wavelength (nm)",
        main = 'No Noise')


# global_spectra_comparison(rtr1, 
#                           rtr2, 
#                           array1_short_name = 'rtr_noise',
#                           array2_short_name = 'rtr_noNoise',
#                           nc.outpath = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/',
#                           month = 8)


# # retreived vs retrieved
# rtr1 <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2016.nc')) %>%
#     ncvar_get('HDR', start = c(2,1,1,8), count = c(210,-1,-1,1)) %>% 
#     aperm(c(2,3,1))%>%
#     rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
#     mask(shp) 
# 
# 
# rtr2 <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2016_noNoise.nc')) %>%
#     ncvar_get('HDR', start = c(2,1,1,8), count = c(210,-1,-1,1)) %>% 
#     aperm(c(2,3,1))%>%
#     rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
#     mask(shp) 
# 
# dif <- rtr2-rtr1
# pdif <- (rtr2-rtr1)/rtr1*100; pdif
# pdif[is.infinite(pdif)] <- NA; #pdif
# pdif[pdif>50] <- NA; #pdif
# pdif[pdif< -50] <- NA; #pdif
# mean_dif <- global(dif, 'mean', na.rm = T)
# mean_pdif <- global(pdif, 'mean', na.rm = T)
# sd_dif <- global(dif, 'sd', na.rm = T)
# sd_pdif <- global(pdif, 'sd', na.rm = T)
# df <- data.frame(wl = seq(410,2500,10), mean = mean_dif$mean, sd = sd_dif$sd, pmean = mean_pdif$mean, psd = sd_pdif$sd)
# 
# ggplot(df) +
#     geom_hline(yintercept = 0) +
#     geom_line(aes(wl, pmean)) +
#     geom_ribbon(aes(x = wl, ymin = pmean-psd, ymax = pmean+psd), alpha = 0.4) +
#     theme_pubr(base_size = 20) +
#     ggtitle('SBG Retrievals - Retrievals (no noise)')+
#     labs(x = 'Wavelength (nm)', y = '% Difference')
# 
# mean(df$pmean); mean(df$psd)


# 2022 --------------------------------------------------------------------

rfl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2022.nc')) %>%
    ncvar_get('DR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>%
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) 
rfl[rfl==0] <- NA

rtr <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2022.nc')) %>%
    ncvar_get('HDR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>% 
    aperm(c(2,3,1))%>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) 
rtr[rtr==0] <- NA


dif <- rtr-rfl
mae <- abs(dif) # MAE
pmae <- mae/rfl*100; # %MAE
pmae[is.infinite(pmae)] <- NA; #pmae
pmae[pmae>50] <- NA; #pmae
pmae[pmae< -50] <- NA; #pmae

rfl.sd <- unlist(as.vector(global(rfl, 'sd', na.rm = T))) # Global reflectance SD (211 values)
dev.f <- mae/rfl.sd*100 # deviation fraction

mean_dif <- global(dif, 'mean', na.rm = T); sd_dif <- global(dif, 'sd', na.rm = T)
mean_mae <- global(mae, 'mean', na.rm = T); sd_mae <- global(mae, 'sd', na.rm = T)
mean_pmae <- global(pmae, 'mean', na.rm = T); sd_pmae <- global(pmae, 'sd', na.rm = T)
mean_devf <- global(dev.f, 'mean', na.rm = T); sd_devf <- global(dev.f, 'sd', na.rm = T)

df <- data.frame(wl = wl, 
                 dif = mean_dif$mean, dif_sd = sd_dif$sd, 
                 mae = mean_mae$mean, mae_sd = sd_mae$sd,
                 pmae = mean_pmae$mean, pmae_sd = sd_pmae$sd,
                 df = mean_devf$mean, df_sd = sd_devf$sd)
df[df$wl %in% bad.wl,2:9] <- NA


figure(
    ggplot(df) +
        geom_hline(yintercept = 0) +
        geom_line(aes(wl, dif)) +
        geom_ribbon(aes(x = wl, ymin = dif-dif_sd, ymax = dif+dif_sd), alpha = 0.3) +
        theme_pubr(base_size = 20) +
        ggtitle('2016 Retrievals - Reflectance')+
        labs(x = 'Wavelength (nm)', y = 'Difference'),
    path.name = paste0(op,'2016_error.png'), width = 8, height = 5, save = F
)

figure(
    ggplot(df) +
        geom_hline(yintercept = 0) +
        geom_line(aes(wl, df)) +
        geom_ribbon(aes(x = wl, ymin = df-df_sd, ymax = df+df_sd), alpha = 0.3) +
        theme_pubr(base_size = 20) +
        ggtitle('Deviation fraction (%)')+
        labs(x = 'Wavelength (nm)', y = '%'),
    path.name = paste0(op,'2016_error.png'), width = 8, height = 5, save = F
)

mean(df$pmae, na.rm = T); mean(df$pmae_sd, na.rm = T) # 3.0 +/- 3.5 %
mean(df$df, na.rm = T); mean(df$df_sd, na.rm = T) # 3.6 +/- 4.5 %



# lineplots ---------------------------------------------------------------

wl <- seq(410,2500,10)
df <- data.frame(wl = rep(wl),
                 mean = rep(0),
                 upr = rep(0),
                 lwr = rep(0),
                 median = rep(0),
                 max = rep(0), 
                 rmse = rep(0),
                 mae = rep(0))

for (i in wl) {
    r1 <- rast(rfl[ , , which(wl==i)])
    r2 <- rast(rtr[ , , which(wl==i)])
    resids <- r2 - r1
    rmse <- sqrt(mean(resids^2, na.rm = T))#/r1 * 100
    mae <- mean(abs(resids), na.rm = T)/r1 * 100
    df$mean[df$wl==i] <- as.numeric(global(resids, mean, na.rm = T))
    # df$upr[df$wl==i] <- df$mean[df$wl==i] + as.numeric(global(resids, sd, na.rm = T))
    # df$lwr[df$wl==i] <- df$mean[df$wl==i] - as.numeric(global(resids, sd, na.rm = T))
    df$median[df$wl==i] <- as.numeric(global(resids, median, na.rm = T))
    df$max[df$wl==i] <- as.numeric(global(resids, max, na.rm = T))
    df$rmse[df$wl==i] <- as.numeric(global(rmse, mean, na.rm = T))
    df$mae[df$wl==i] <- as.numeric(global(mae, mean, na.rm = T))
}


ggplot(df) +
    geom_hline(yintercept = 0, alpha = 0.7) +
    geom_line(aes(wl, mean), linewidth = 1) +
    geom_line(aes(wl, median), linewidth = 1, color = 'red') +
    geom_line(aes(wl, rmse), linewidth = 1, color = 'grey') +
    # geom_line(aes(wl, upr), linewidth = 1, linetype = 'longdash', color = 'grey') +
    # geom_line(aes(wl, mean), linewidth = 1, linetype = 'longdash', color = 'grey') +
    theme_bw(base_size = 20) +
    labs(x = 'Wavelength', y = 'Global average difference', title = 'Retrievals - Reflectances')

ggplot(df) +
    geom_line(aes(wl, mean), linewidth = 1) +
    theme_bw(base_size = 20) +
    labs(x = 'Wavelength', y = 'Global average difference', title = 'Retrievals - Reflectances')

ggplot(df) +
    geom_line(aes(wl, max), linewidth = 1) +
    theme_bw(base_size = 20) +
    labs(x = 'Wavelength', y = 'Global average difference', title = 'Retrievals - Reflectances')

ggplot(df) +
    # geom_line(aes(wl, rmse), linewidth = 1) +
    geom_line(aes(wl, mae), linewidth = 1, color = 'red') +
    theme_bw(base_size = 20) +
    labs(x = 'Wavelength', y = 'Global Percent error', title = 'Retrievals - Reflectances')



# Geographic differences -----------------------------------------------------

dif <- nc_open('~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/Dif_rtr15_rfl15_month_8.nc')
mask_rfl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2015.nc')) %>%
    ncvar_get('DR', start = c(1,1,2,8), count = c(-1,-1,1,1)) %>%
    aperm(c(2,1)) %>%
    rast(ext = c(-180,180,-90,90), crs = 'EPSG:4326')
rfl[rfl==0] <- NA

max.r <- dif %>%
    ncvar_get(attributes(dif$var)$name, start = c(1,1,1), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>%
    mask(mask_rfl)
plot(max.r)
global(max.r, median, na.rm = T)
global(max.r, range, na.rm = T)

wl.r <- dif %>%
    ncvar_get(attributes(dif$var)$name, start = c(1,1,2), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>%
    mask(mask_rfl)
plot(wl.r)
global(wl.r, median, na.rm = T)

mean.r <- dif %>%
    ncvar_get(attributes(dif$var)$name, start = c(1,1,3), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>%
    mask(mask_rfl)
plot(mean.r)
global(mean.r, median, na.rm = T) # global = -0.0002 (<-50 - >100)
global(mean.r, range, na.rm = T)

sd.r <- dif %>%
    ncvar_get(attributes(dif$var)$name, start = c(1,1,4), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>%
    mask(mask_rfl)
plot(sd.r)
global(sd.r, median, na.rm = T) # global = 0.003
global(sd.r, range, na.rm = T)

cv.r <- dif %>%
    ncvar_get(attributes(dif$var)$name, start = c(1,1,5), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>%
    mask(mask_rfl)
plot(cv.r)
global(cv.r, median, na.rm = T) # global = 3.61
global(cv.r, range, na.rm = T)

rmse.r <- dif %>%
    ncvar_get(attributes(dif$var)$name, start = c(1,1,6), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>%
    mask(mask_rfl)
plot(rmse.r)
global(rmse.r, median, na.rm = T) # global = 0.0032% (2 - >100)
global(rmse.r, range, na.rm = T)

prmse.r <- dif %>%
    ncvar_get(attributes(dif$var)$name, start = c(1,1,7), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>%
    mask(mask_rfl)
prmse.r[prmse.r > 100] <- NA
prmse.r[prmse.r< -50] <- NA
plot(prmse.r)
global(prmse.r, median, na.rm = T) # global RMSE = 1.61% (2 - >100)
global(prmse.r, sd, na.rm = T) # sd 3.3%
global(prmse.r, range, na.rm = T)

mae.r <- dif %>%
    ncvar_get(attributes(dif$var)$name, start = c(1,1,8), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>%
    mask(mask_rfl)
plot(mae.r)
global(mae.r, median, na.rm = T) # global = 0.0025% (2 - >100)
global(mae.r, sd, na.rm = T) # global = 0.003%
global(mae.r, range, na.rm = T)

pmae.r <- dif %>%
    ncvar_get(attributes(dif$var)$name, start = c(1,1,9), count = c(-1,-1,1)) %>% 
    aperm(c(2,1)) %>%
    rast(ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>%
    mask(mask_rfl)
pmae.r[pmae.r>100] <- NA
pmae.r[pmae.r< -50] <- NA
plot(pmae.r)
global(pmae.r, median, na.rm = T) # global = 1.23%
global(pmae.r, sd, na.rm = T) # global = 3.15%
global(pmae.r, range, na.rm = T)




# Combination -------------------------------------------------------------

mask <- rmse.r<0.001


wl <- seq(410,2500,10)
df <- data.frame(wl = rep(wl),
                 mean = rep(0),
                 upr = rep(0),
                 lwr = rep(0),
                 median = rep(0),
                 max = rep(0), 
                 rmse = rep(0),
                 mae = rep(0))

for (i in wl) {
    r1 <- rast(rfl[ , , which(wl==i)], ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>% crop(ext(c(-120, -110, 35, 45))) #mask(mask, maskvalue=1, inverse=T)
    r2 <- rast(rtr[ , , which(wl==i)], ext = c(-180,180,-90,90), crs = 'EPSG:4326') %>% crop(ext(c(-120, -110, 35, 45))) #mask(mask, maskvalue=1, inverse=T)
    resids <- r2 - r1
    rmse <- sqrt(mean(resids^2, na.rm = T))#/r1 * 100
    mae <- mean(abs(resids), na.rm = T)/r1 * 100
    df$mean[df$wl==i] <- as.numeric(global(resids, mean, na.rm = T))
    # df$upr[df$wl==i] <- df$mean[df$wl==i] + as.numeric(global(resids, sd, na.rm = T))
    # df$lwr[df$wl==i] <- df$mean[df$wl==i] - as.numeric(global(resids, sd, na.rm = T))
    df$median[df$wl==i] <- as.numeric(global(resids, median, na.rm = T))
    df$max[df$wl==i] <- as.numeric(global(resids, max, na.rm = T))
    df$rmse[df$wl==i] <- as.numeric(global(rmse, mean, na.rm = T))
    df$mae[df$wl==i] <- as.numeric(global(mae, mean, na.rm = T))
}


ggplot(df) +
    geom_hline(yintercept = 0, alpha = 0.7) +
    geom_line(aes(wl, mean), linewidth = 1) +
    geom_line(aes(wl, median), linewidth = 1, color = 'red') +
    geom_line(aes(wl, rmse), linewidth = 1, color = 'grey') +
    theme_bw(base_size = 20) +
    labs(x = 'Wavelength', y = 'Global average difference', title = 'Retrievals - Reflectances')

# plotting ------------------------------------------------

p1 <- ggplot() +
    geom_spatraster(data=max.r) +
    scale_fill_gradient(low="blue", high="red", na.value = 'transparent') +
    theme_void(base_size = 15) +
    labs(title = paste0('Max = ', round(global(max.r, median, na.rm = T), 3)), fill = '')

p2 <- ggplot() +
    geom_spatraster(data=wl.r) +
    scale_fill_gradient(low="blue", high="red", na.value = 'transparent') +
    theme_void(base_size = 15) +
    labs(title = paste0('Band of max difference = ', round(global(wl.r, median, na.rm = T), 3)), fill = '')

p3 <- ggplot() +
    geom_spatraster(data=mean.r) +
    scale_fill_gradient(low="blue", high="red", na.value = 'transparent') +
    theme_void(base_size = 15) +
    labs(title = paste0('Mean = ', round(global(mean.r, median, na.rm = T), 3)), fill = '')

p4 <- ggplot() +
    geom_spatraster(data=sd.r) +
    scale_fill_gradient(low="blue", high="red", na.value = 'transparent') +
    theme_void(base_size = 15) +
    labs(title = paste0('SD = ', round(global(sd.r, median, na.rm = T), 3)), fill = '')

p5 <- ggplot() +
    geom_spatraster(data=cv.r) +
    scale_fill_gradient(low="blue", high="red", na.value = 'transparent') +
    theme_void(base_size = 15) +
    labs(title = paste0('CV = ', round(global(cv.r, median, na.rm = T), 3)), fill = '')

p6 <- ggplot() +
    geom_spatraster(data=rmse.r) +
    scale_fill_gradient(low="blue", high="red", na.value = 'transparent') +
    theme_void(base_size = 15) +
    labs(title = paste0('RMSE = ', round(global(rmse.r, median, na.rm = T), 3)), fill = '') 

p7 <- ggplot() +
    geom_spatraster(data=prmse.r) +
    scale_fill_gradient(low="blue", high="red", na.value = 'transparent') +
    theme_void(base_size = 15) +
    labs(title = paste0('% RMSE = ', round(global(prmse.r, median, na.rm = T), 3)), fill = '') 

p8 <- ggplot() +
    geom_spatraster(data=mae.r) +
    scale_fill_gradient(low="blue", high="red", na.value = 'transparent') +
    theme_void(base_size = 15) +
    labs(title = paste0('MAE = ', round(global(mae.r, median, na.rm = T), 3)), fill = '') 

p9 <- ggplot() +
    geom_spatraster(data=pmae.r) +
    scale_fill_gradient(low="cyan", high="red4", 
                        limits = c(0,5),
                        na.value = 'transparent') +
    theme_void(base_size = 15) +
    labs(title = paste0('% MAE = ', round(global(pmae.r, median, na.rm = T), 3)), fill = '') 


plot <- ggarrange(#p.MODIS,
    p1,p2,
    p3,p4,
    p6,
    p7,p8,
    p9,
    nrow = 4, ncol = 2, align='hv')

annotate_figure(
    plot, 
    top = text_grob(paste0("Comparison between Retrievals and Reflectance"), 
                    face = "bold", size = 25))

