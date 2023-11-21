source('~/R/clean.R')
library(terra)
library(ncdf4)
library(tidyverse)
source('~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/R_scripts/lpj_plsr_functions.R')
`%!in%` = Negate(`%in%`)

# paths -------------------------------------------------------------------

trait.path <- '~/Current Projects/SBG/Trait mapping/Global_trait_maps_Moreno_Martinez_2018_Version2_3km_resolution'
lpj.path <- '~/Current Projects/SBG/LPJ/Reflectance_Data/version2/'


ldmc.name <- 'LDMC_3km_v1'
n.name <- 'LNC_3km_v1'
p.name <- 'LPC_3km_v1'
sla.name <- 'SLA_3km_v1'

lpj.nc <- 'lpj-prosail_levelC_DR_Version021_m_2015.nc'
wl <- seq(400,2500,10)
bad.wl <- c(400,1350,1360,1370,1380,1390,1400,1830,1840,1850,1860,1870,1880,1890,1900,1910) # identified by SNR_plot.R

# create PLSR data.frame --------------------------------------------------
print('Reading in LPJ array')
lpj.array <- nc_open(file.path(lpj.path, lpj.nc)) %>%
    ncvar_get('DR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>%
    aperm(c(2,1,3))

lpj.r <- rast(lpj.array, ext = c(-180,180,-90,90), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>%
    subset(which(wl %!in% bad.wl))

lpj.vis <- lpj.r %>% subset(which(wl<=700))
lpj.mask <- lpj.vis>0.3
lpj.mask <- app(lpj.mask, sum, na.rm = T); 
lpj.mask[lpj.mask>1] <- 1; plot(lpj.mask)
lpj.r <- lpj.r %>% mask(lpj.mask, maskvalues = 1)
lpj.r[lpj.r==0] <- NA

cells <- vect(crds(lpj.r, na.rm = F), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
plot(lpj.r$lyr.2)


# All data ----------------------------------------------------------------

print('Extract points from TRY trait maps')
ldmc <- rast(file.path(trait.path, ldmc.name, paste0(ldmc.name,'.tif')))
ldmc[ldmc<=0] <- NA
global(ldmc, range, na.rm = T)
ldmc <- ldmc %>% terra::extract(cells, ID = F)

lnc <- rast(file.path(trait.path, n.name, paste0(n.name,'.tif')))
lnc[lnc<=0] <- NA
global(lnc, range, na.rm = T)
lnc <- lnc %>% terra::extract(cells, ID = F)

lpc <- rast(file.path(trait.path, p.name, paste0(p.name,'.tif'))) 
lpc[lpc<=0] <- NA
global(lpc, range, na.rm = T)
lpc <- lpc %>% terra::extract(cells, ID = F)

sla <- rast(file.path(trait.path, sla.name, paste0(sla.name,'.tif')))
sla[sla<=0] <- NA
global(sla, range, na.rm = T)
sla <- sla %>% terra::extract(cells, ID = F)

lpj.df <- as.data.frame(lpj.r, xy = T, na.rm = F)
plsr.data <- cbind.data.frame(lpj.df, ldmc, lnc, lpc, sla) %>% na.exclude()
names(plsr.data) <- c('x', 'y', paste0('wave',wl[wl %!in% bad.wl]), ldmc.name, n.name, p.name, sla.name)

# run PLSR 
# This is for August

ldmc.coefs <- runPLSR(plsr.data, data.var = ldmc.name, band.prefix = 'wave', train.size = 2000, plots = F,
                     jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 7,
                     wl = wl[wl %!in% bad.wl])

n.coefs <- runPLSR(plsr.data, data.var = n.name, band.prefix = 'wave', train.size = 2000, plots = F,
                   jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 7,
                   wl = wl[wl %!in% bad.wl])

p.coefs <- runPLSR(plsr.data, data.var = p.name, band.prefix = 'wave', train.size = 2000, plots = F,
                   jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 7,
                   wl = wl[wl %!in% bad.wl])

sla.coefs <- runPLSR(plsr.data, data.var = sla.name, band.prefix = 'wave', train.size = 2000, plots = F,
                     jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 7,
                     wl = wl[wl %!in% bad.wl])


coeff.df <- data.frame(coeff = c('Intercept', wl[wl %!in% bad.wl]), ldmc = ldmc.coefs, n = n.coefs, p = p.coefs, sla = sla.coefs)
write_csv(coeff.df, '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/LPJ-PROSAIL_PLSR_coefficients_masked_WL_August2015.csv')

ldmc.map <- trait.map(lpj.r, coeffs = coeff.df$ldmc[-1], intercept = coeff.df$ldmc[1], coeffs_wl = wl[wl %!in% bad.wl])
ldmc.map[ldmc.map<=0] <- NA
global(ldmc.map, median, na.rm = T)
global(ldmc.map, sd, na.rm = T)
global(ldmc.map, range, na.rm = T)

n.map <- trait.map(lpj.r,  coeffs = coeff.df$n[-1], intercept = coeff.df$n[1], coeffs_wl = wl[wl %!in% bad.wl])
n.map[n.map<=0] <- NA
global(n.map, median, na.rm = T)
global(n.map, sd, na.rm = T)
global(n.map, range, na.rm = T)

p.map <- trait.map(lpj.r, coeff.df$p[-1], coeff.df$p[1], coeffs_wl = wl[wl %!in% bad.wl])
p.map[p.map<=0] <- NA
global(p.map, median, na.rm = T)
global(p.map, sd, na.rm = T)
global(p.map, range, na.rm = T)

sla.map <- trait.map(lpj.r, coeff.df$sla[-1], coeff.df$sla[1], coeffs_wl = wl[wl %!in% bad.wl])
sla.map[sla.map<=0] <- NA
global(sla.map, median, na.rm = T)
global(sla.map, sd, na.rm = T)
global(sla.map, range, na.rm = T)




# Inner 90% of data ---------------------------------------------------------------

print('Extract points from TRY trait maps')
ldmc <- rast(file.path(trait.path, ldmc.name, paste0(ldmc.name,'.tif')))
ldmc[ldmc<=0] <- NA
global(ldmc, range, na.rm = T)
lwr.bnd <- global(ldmc, function(x) quantile(x, 0.05, na.rm=T)) %>% as.numeric()
upr.bnd <- global(ldmc, function(x) quantile(x, 0.95, na.rm=T)) %>% as.numeric()
ldmc[ldmc < lwr.bnd] <- NA
ldmc[ldmc > upr.bnd] <- NA
ldmc <- ldmc %>% terra::extract(cells, ID = F)
range(ldmc, na.rm = T)

lnc <- rast(file.path(trait.path, n.name, paste0(n.name,'.tif')))
lnc[lnc<=0] <- NA
global(lnc, range, na.rm = T)
lwr.bnd <- global(lnc, function(x) quantile(x, 0.05, na.rm=T)) %>% as.numeric()
upr.bnd <- global(lnc, function(x) quantile(x, 0.95, na.rm=T)) %>% as.numeric()
lnc[lnc < lwr.bnd] <- NA
lnc[lnc > upr.bnd] <- NA
lnc <- lnc %>% terra::extract(cells, ID = F)
range(lnc, na.rm = T)

lpc <- rast(file.path(trait.path, p.name, paste0(p.name,'.tif'))) 
lpc[lpc<=0] <- NA
global(lpc, range, na.rm = T)
lwr.bnd <- global(lpc, function(x) quantile(x, 0.05, na.rm=T)) %>% as.numeric()
upr.bnd <- global(lpc, function(x) quantile(x, 0.95, na.rm=T)) %>% as.numeric()
lpc[lpc < lwr.bnd] <- NA
lpc[lpc > upr.bnd] <- NA
lpc <- lpc %>% terra::extract(cells, ID = F)
range(lpc, na.rm = T)

sla <- rast(file.path(trait.path, sla.name, paste0(sla.name,'.tif')))
sla[sla<=0] <- NA
global(sla, range, na.rm = T)
lwr.bnd <- global(sla, function(x) quantile(x, 0.05, na.rm=T)) %>% as.numeric()
upr.bnd <- global(sla, function(x) quantile(x, 0.95, na.rm=T)) %>% as.numeric()
sla[sla < lwr.bnd] <- NA
sla[sla > upr.bnd] <- NA
sla <- sla %>% terra::extract(cells, ID = F)
range(sla, na.rm = T)


lpj.df <- as.data.frame(lpj.r, xy = T, na.rm = F)
plsr.data <- cbind.data.frame(lpj.df, ldmc, lnc, lpc, sla) %>% na.exclude()
names(plsr.data) <- c('x', 'y', paste0('wave',seq(410,2500,10)), ldmc.name, n.name, p.name, sla.name)

# run PLSR 

# This is for August
ldmc.coefs <- runPLSR(plsr.data, data.var = ldmc.name, band.prefix = 'wave', train.size = 2000, plots = F,
                      jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 7,
                      wl = seq(410, 2500, 10))

n.coefs <- runPLSR(plsr.data, data.var = n.name, band.prefix = 'wave', train.size = 2000, plots = F,
                   jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 7,
                   wl = seq(410, 2500, 10))

p.coefs <- runPLSR(plsr.data, data.var = p.name, band.prefix = 'wave', train.size = 2000, plots = F,
                   jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 7,
                   wl = seq(410, 2500, 10))

sla.coefs <- runPLSR(plsr.data, data.var = sla.name, band.prefix = 'wave', train.size = 2000, plots = F,
                     jk.test = F, jk.prop = 0.15, jk.iterations = 20, jk.comps = 7,
                     wl = seq(410, 2500, 10))


coeff.df <- data.frame(coeff = c('Intercept', seq(410,2500,10)), ldmc = ldmc.coefs, n = n.coefs, p = p.coefs, sla = sla.coefs)
write_csv(coeff.df, '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/LPJ-PROSAIL_PLSR_coefficients_inner90_August2015.csv')


ldmc.map <- trait.map(lpj.r, coeffs = coeff.df$ldmc[-1], intercept = coeff.df$ldmc[1], coeffs_wl = seq(410,2500,10))
ldmc.map[ldmc.map<=0] <- NA
global(ldmc.map, median, na.rm = T)
global(ldmc.map, sd, na.rm = T)
global(ldmc.map, range, na.rm = T)

n.map <- trait.map(lpj.r,  coeffs = coeff.df$n[-1], intercept = coeff.df$n[1], coeffs_wl = seq(410,2500,10))
n.map[n.map<=0] <- NA
global(n.map, median, na.rm = T)
global(n.map, sd, na.rm = T)
global(n.map, range, na.rm = T)

p.map <- trait.map(lpj.r, coeff.df$p[-1], coeff.df$p[1], coeffs_wl = seq(410,2500,10))
p.map[p.map<=0] <- NA
global(p.map, median, na.rm = T)
global(p.map, sd, na.rm = T)
global(p.map, range, na.rm = T)

sla.map <- trait.map(lpj.r, coeff.df$sla[-1], coeff.df$sla[1], coeffs_wl = seq(410,2500,10))
sla.map[sla.map<=0] <- NA
global(sla.map, median, na.rm = T)
global(sla.map, sd, na.rm = T)
global(sla.map, range, na.rm = T)


# plotting ----------------------------------------------------------------


library(ggplot2)
library(tidyterra)
library(ggpubr)
p1 <- ggplot() +
    geom_spatraster(data = ldmc.map) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), 
                         limits = c(0.25, 0.4),
                         na.value = 'transparent') +
    theme_void(base_size = 20) +
    labs(title = 'Estimated LDMC (g/g)')


p2 <- ggplot() +
    geom_spatraster(data = n.map) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), 
                         limits = c(15, 25),
                         na.value = 'transparent') +
    theme_void(base_size = 20) +
    labs(title = 'Estimated Leaf N (mg/g)')

p3 <- ggplot() +
    geom_spatraster(data = p.map) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), 
                         limits = c(1, 1.7),
                         na.value = 'transparent') +
    theme_void(base_size = 20) +
    labs(title = 'Estimated Leaf P (mg/g)')

p4 <- ggplot() +
    geom_spatraster(data = sla.map) +
    scale_fill_gradientn(colors = c("wheat2", "darkgreen"), 
                         limits = c(7, 20),
                         na.value = 'transparent') +
    theme_void(base_size = 20) +
    labs(title = 'Estimated SLA (mm2/mg)')


ggarrange(p1,p2,p3,p4,
          nrow = 2, ncol = 2,
          align = 'hv')

