source('~/R/Clean.R')
source('~/R/figure.R')
library(ncdf4)
library(terra)
library(tidyverse)
library(ggplot2)
library(tidyterra)
library(ggpubr)
library(rworldmap)
source('~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/R_scripts/lpj_plsr_functions.R')
`%!in%` = Negate(`%in%`)

# set data path
dp <- '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/'
map <- vect(getMap())
wl <- seq(400,2500,10)
bad.wl <- c(400,1350,1360,1370,1380,1390,1400,1830,1840,1850,1860,1870,1880,1890,1900,1910) # identified by SNR_plot.R


# Exploring Temp anomalies ----------------------------------------------------------

# read in temp anaomalies nc and quickview
tmp.anomalies <- nc_open(file.path(dp, 'Anomalies/tmp_anomaly_Aug_2022.nc')) %>%
    ncvar_get('tmp') %>%
    aperm(c(2,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    flip
plot(tmp.anomalies, range = c(-5,5), smooth = T)
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
    vect(type="polygon", crs = 'EPSG:4326')
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
    vect(type="polygon", crs = 'EPSG:4326')
plot(r2, add = T)
global(tmp.anomalies %>% crop(r2), mean, na.rm = T)
# R2 average temp anomaly = 3.86 degrees

# figure(
#     ggplot() +
#         geom_spatraster(data = tmp.anomalies) +
#         geom_spatvector(data = map, alpha = 0, color = 'grey70') +
#         geom_spatvector(data = r1, alpha = 0, linewidth = 1, color = 'black') +
#         geom_spatvector(data = r2, alpha = 0, linewidth = 1, color = 'black') +
#         scale_fill_gradient2(low = 'blue', high = 'red3', na.value = 'transparent') +
#         theme_void(base_size = 20) +
#         labs(fill = ''),
#     path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/temp_2015_regions.png',
#     height = 6, width = 12,
#     save = T)


# Exploring Precip anomalies --------------------------------------------------------


# read in temp anaomalies nc and quickview
pre.anomalies <- nc_open(file.path(dp, 'Anomalies/pre_anomaly_Aug_2016.nc')) %>%
    ncvar_get('pre') %>%
    aperm(c(2,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    flip
pre.anomalies[pre.anomalies>400] <- NA
plot(pre.anomalies, range = c(-200,200), smooth = T)
global(pre.anomalies, mean, na.rm = T)
# Global average temp anomaly = 0.79 degrees

# Create region 1 box
xmin <- -5
xmax <- 15
ymin <-  35
ymax <- 55
r3 <- matrix(c(xmin, ymin, 
               xmax, ymin, 
               xmax, ymax, 
               xmin, ymax, 
               xmin, ymin), 
             ncol=2, byrow=TRUE) %>%
    vect(type="polygon", crs = 'EPSG:4326')
plot(r3, add = T)
global(pre.anomalies %>% crop(r3), mean, na.rm = T)
# R1 average temp anomaly = 2.34 degrees


# Create region 2 box
xmin <- 160
xmax <- 180
ymin <- -32
ymax <- -52
r4 <- matrix(c(xmin, ymin, 
               xmax, ymin, 
               xmax, ymax, 
               xmin, ymax, 
               xmin, ymin), 
             ncol=2, byrow=TRUE) %>%
    vect(type="polygon", crs = 'EPSG:4326')
plot(r4, add = T)
global(pre.anomalies %>% crop(r4), mean, na.rm = T)

pre.anomalies[pre.anomalies>50] <- 50
pre.anomalies[pre.anomalies< -50] <- -50

# figure(
#     ggplot() +
#         geom_spatraster(data = pre.anomalies) +
#         geom_spatvector(data = map, alpha = 0, color = 'grey70') +
#         geom_spatvector(data = r3, alpha = 0, linewidth = 1, color = 'black') +
#         geom_spatvector(data = r4, alpha = 0, linewidth = 1, color = 'black') +
#         scale_fill_gradient2(low = 'gold', high = 'blue3', na.value = 'transparent') +
#         theme_void(base_size = 20) +
#         labs(fill = ''),
#     path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/pre_2015_regions.png',
#     height = 6, width = 12,
#     save = T)


# 2015 Retrieved vs Simulated reflectances ------------------------------------------------

coeffs <- read_csv('~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/data/LPJ-PROSAIL_PLSR_coefficients_masked_WL_August2015.csv')
head(coeffs)


rfl15 <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelC_DR_Version021_m_2015.nc') %>%
    ncvar_get('DR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>%
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    subset(which(wl %!in% bad.wl))

vis <- rfl15 %>% subset(which(wl<=700))
lpj.mask <- vis>0.3
lpj.mask <- app(lpj.mask, sum, na.rm = T)
lpj.mask[lpj.mask>1] <- 1; plot(lpj.mask)
rfl15 <- rfl15 %>% mask(lpj.mask, maskvalues = 1)
rfl15[rfl15==0] <- NA
plot(rfl15$lyr.2)


rtr15 <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelE_retrieved-HDR_Version021_m_2015.nc') %>%
    ncvar_get('HDR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>% 
    aperm(c(2,3,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    subset(which(wl %!in% bad.wl))

rtr15 <- rtr15 %>% mask(lpj.mask, maskvalues = 1)
rtr15[rtr15==0] <- NA
plot(rtr15$lyr.2)
names(rtr15) <- wl[wl %!in% bad.wl]
rtr15.r1 <- crop(rtr15, r1)
rtr15.r2 <- crop(rtr15, r2)
rtr15.r3 <- crop(rtr15, r3)
rtr15.r4 <- crop(rtr15, r4)

ldmc.rfl15 <- trait.map(rfl15, coeffs = coeffs$ldmc[-1], intercept = coeffs$ldmc[1], coeffs_wl = wl[wl %!in% bad.wl])
n.rfl15 <- trait.map(rfl15, coeffs = coeffs$n[-1], intercept = coeffs$n[1], coeffs_wl = wl[wl %!in% bad.wl])
p.rfl15 <- trait.map(rfl15, coeffs = coeffs$p[-1], intercept = coeffs$p[1], coeffs_wl = wl[wl %!in% bad.wl]) 
sla.rfl15 <- trait.map(rfl15, coeffs = coeffs$sla[-1], intercept = coeffs$sla[1], coeffs_wl = wl[wl %!in% bad.wl])

ldmc.rtr15 <- trait.map(rtr15, coeffs = coeffs$ldmc[-1], intercept = coeffs$ldmc[1], coeffs_wl = wl[wl %!in% bad.wl])
n.rtr15 <- trait.map(rtr15, coeffs = coeffs$n[-1], intercept = coeffs$n[1], coeffs_wl = wl[wl %!in% bad.wl])   
p.rtr15 <- trait.map(rtr15, coeffs = coeffs$p[-1], intercept = coeffs$p[1], coeffs_wl = wl[wl %!in% bad.wl]) 
sla.rtr15 <- trait.map(rtr15, coeffs = coeffs$sla[-1], intercept = coeffs$sla[1], coeffs_wl = wl[wl %!in% bad.wl])


# differences
rast1 <- ldmc.rfl15
rast2 <- ldmc.rtr15

rast1 <- n.rfl15
rast2 <- n.rtr15

rast1 <- p.rfl15
rast2 <- p.rtr15

rast1 <- sla.rfl15
rast2 <- sla.rtr15

dif <- rast2-rast1
mae <- abs(dif) # MAE
pmae <- mae/rast1*100; # %MAE
# pmae[is.infinite(pmae)] <- NA; #pmae
# pmae[pmae>50] <- NA; #pmae
# pmae[pmae< -50] <- NA; #pmae

global(pmae, 'mean', na.rm = T); global(pmae, 'sd', na.rm = T)

rfl.sd <- unlist(as.vector(global(rast1, 'sd', na.rm = T))) # Global reflectance SD (211 values)
dev.f <- mae/rfl.sd*100 # deviation fraction

global(dev.f, 'mean', na.rm = T); global(dev.f, 'sd', na.rm = T)


ldmc15.r1 <- crop(ldmc.rtr15, r1)
ldmc15.r2 <- crop(ldmc.rtr15, r2)
ldmc15.r3 <- crop(ldmc.rtr15, r3)
ldmc15.r4 <- crop(ldmc.rtr15, r4)

n15.r1 <- crop(n.rtr15, r1)
n15.r2 <- crop(n.rtr15, r2)
n15.r3 <- crop(n.rtr15, r3)
n15.r4 <- crop(n.rtr15, r4)

p15.r1 <- crop(p.rtr15, r1)
p15.r2 <- crop(p.rtr15, r2)
p15.r3 <- crop(p.rtr15, r3)
p15.r4 <- crop(p.rtr15, r4)

sla15.r1 <- crop(sla.rtr15, r1)
sla15.r2 <- crop(sla.rtr15, r2)
sla15.r3 <- crop(sla.rtr15, r3)
sla15.r4 <- crop(sla.rtr15, r4)

# 2016 --------------------------------------------------------------------

rtr16 <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelE_retrieved-HDR_Version021_m_2016.nc') %>%
    ncvar_get('HDR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>% 
    aperm(c(2,3,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    subset(which(wl %!in% bad.wl))

rtr16 <- rtr16 %>% mask(lpj.mask, maskvalues = 1)
rtr16[rtr16==0] <- NA
plot(rtr16$lyr.2)
names(rtr16) <- wl[wl %!in% bad.wl]
rtr16.r3 <- crop(rtr16, r3)
rtr16.r4 <- crop(rtr16, r4)


ldmc.rtr16 <- trait.map(rtr16, coeffs = coeffs$ldmc[-1], intercept = coeffs$ldmc[1], coeffs_wl = wl[wl %!in% bad.wl])
n.rtr16 <- trait.map(rtr16, coeffs = coeffs$n[-1], intercept = coeffs$n[1], coeffs_wl = wl[wl %!in% bad.wl])
p.rtr16 <- trait.map(rtr16, coeffs = coeffs$p[-1], intercept = coeffs$p[1], coeffs_wl = wl[wl %!in% bad.wl])
sla.rtr16 <- trait.map(rtr16, coeffs = coeffs$sla[-1], intercept = coeffs$sla[1], coeffs_wl = wl[wl %!in% bad.wl])


ldmc16.r3 <- crop(ldmc.rtr16, r3)
ldmc16.r4 <- crop(ldmc.rtr16, r4)

n16.r3 <- crop(n.rtr16, r3)
n16.r4 <- crop(n.rtr16, r4)

p16.r3 <- crop(p.rtr16, r3)
p16.r4 <- crop(p.rtr16, r4)

sla16.r3 <- crop(sla.rtr16, r3)
sla16.r4 <- crop(sla.rtr16, r4)


# 2022 --------------------------------------------------------------------

rtr22 <- nc_open('~/Current Projects/SBG/LPJ/Reflectance_Data/version2/lpj-prosail_levelE_retrieved-HDR_Version021_m_2022.nc') %>%
    ncvar_get('HDR', start = c(1,1,1,8), count = c(-1,-1,-1,1)) %>% 
    aperm(c(2,3,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    subset(which(wl %!in% bad.wl))

rtr22 <- rtr22%>% mask(lpj.mask, maskvalues = 1)
rtr22[rtr22==0] <- NA
plot(rtr22$lyr.2)
names(rtr22) <-  wl[wl %!in% bad.wl]
rtr22.r1 <- crop(rtr22, r1)
rtr22.r2 <- crop(rtr22, r2)


ldmc.rtr22 <- trait.map(rtr22, coeffs = coeffs$ldmc[-1], intercept = coeffs$ldmc[1], coeffs_wl = wl[wl %!in% bad.wl])
n.rtr22 <- trait.map(rtr22, coeffs = coeffs$n[-1], intercept = coeffs$n[1], coeffs_wl = wl[wl %!in% bad.wl])
p.rtr22 <- trait.map(rtr22, coeffs = coeffs$p[-1], intercept = coeffs$p[1], coeffs_wl = wl[wl %!in% bad.wl])
sla.rtr22 <- trait.map(rtr22, coeffs = coeffs$sla[-1], intercept = coeffs$sla[1], coeffs_wl = wl[wl %!in% bad.wl])


ldmc22.r1 <- crop(ldmc.rtr22, r1)
ldmc22.r2 <- crop(ldmc.rtr22, r2)

n22.r1 <- crop(n.rtr22, r1)
n22.r2 <- crop(n.rtr22, r2)

p22.r1 <- crop(p.rtr22, r1)
p22.r2 <- crop(p.rtr22, r2)

sla22.r1 <- crop(sla.rtr22, r1)
sla22.r2 <- crop(sla.rtr22, r2)



# Data frame creation, plotting, and statistics ----------------------------------------
# Bootstraps
n <- 1000

# Region 1 - 2022 and 2015
R1 <- cbind.data.frame(
    seq(ncell(ldmc15.r1)),
    ldmc15.r1 %>% as.data.frame(na.rm = F),
    ldmc22.r1 %>% as.data.frame(na.rm = F),
    
    n15.r1 %>% as.data.frame(na.rm = F),
    n22.r1 %>% as.data.frame(na.rm = F),
    
    p15.r1 %>% as.data.frame(na.rm = F),
    p22.r1 %>% as.data.frame(na.rm = F),
    
    sla15.r1 %>% as.data.frame(na.rm = F),
    sla22.r1 %>% as.data.frame(na.rm = F)
)
names(R1) <- c('ID', 'LDMC15', 'LDMC22', 'N15', 'N22', 'P15', 'P22', 'SLA15', 'SLA22')
R1 <- R1 %>% pivot_longer(-ID, names_to = 'TraitYear', values_to = 'TraitValue') %>%
    separate(TraitYear, into = c('Trait', 'Year'), sep = "(?<=\\D)(?=\\d)", convert = T) %>%
    mutate(Year = Year+2000) %>%
    # group_by(Trait) %>%
    mutate(Year = factor(Year),
           Trait = factor(Trait)) %>%
    na.exclude() 

p1 <- ggplot(R1) +
    geom_density(aes(y = TraitValue, fill = Year, color = Year), alpha = 0.3, adjust = 2) +
    facet_wrap(~ Trait, scales = 'free', nrow = 1) +
    scale_fill_manual(values = c('grey30', 'red3')) +
    scale_color_manual(values = c('grey20', 'red4')) +
    # scale_y_log10() +
    theme_pubclean(base_size = 30) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          strip.text = element_text(size = 20, face = 'bold'),
          legend.position = 'bottom')+
    labs(y = 'Trait Value');p1

df <- cbind.data.frame(n = rep(seq(n), each = 4),
                       trait = rep(c('LDMC', 'N','P','SLA'), times = n),
                       delta = rep(0),
                       pdelta = rep(0),
                       p = rep(0))

for (i in seq(n)) {
    y <- sample(dim(R1)[1], n)
    coef <- summary(lm(TraitValue ~ Year, R1 %>% filter(Trait == 'LDMC') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'LDMC'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'LDMC'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'LDMC'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R1 %>% filter(Trait == 'N') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'N'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'N'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'N'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R1 %>% filter(Trait == 'P') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'P'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'P'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'P'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R1 %>% filter(Trait == 'SLA') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'SLA'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'SLA'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'SLA'] <- coef[2,4]
}

df %>% group_by(trait) %>% summarise(mean(delta), sd(delta), mean(pdelta), sd(pdelta), median(p))

R1 %>% group_by(Trait, Year) %>%
    summarise(CV = sd(TraitValue)/mean(TraitValue)*100) %>%
    group_by(Trait) %>%
    mutate(CV_diff = (CV - lag(CV, default = first(CV)))/lag(CV, default = first(CV))*100) %>%
    filter(CV_diff != 0)




# Region 2 - 2022 and 2015
R2 <- cbind.data.frame(
    seq(ncell(ldmc15.r2)),
    ldmc15.r2 %>% as.data.frame(na.rm = F),
    ldmc22.r2 %>% as.data.frame(na.rm = F),
    
    n15.r2 %>% as.data.frame(na.rm = F),
    n22.r2 %>% as.data.frame(na.rm = F),
    
    p15.r2 %>% as.data.frame(na.rm = F),
    p22.r2 %>% as.data.frame(na.rm = F),
    
    sla15.r2 %>% as.data.frame(na.rm = F),
    sla22.r2 %>% as.data.frame(na.rm = F)
)

names(R2) <- c('ID', 'LDMC15', 'LDMC22', 'N15', 'N22', 'P15', 'P22', 'SLA15', 'SLA22')
R2 <- R2 %>% pivot_longer(-ID, names_to = 'TraitYear', values_to = 'TraitValue') %>%
    separate(TraitYear, into = c('Trait', 'Year'), sep = "(?<=\\D)(?=\\d)", convert = T) %>%
    mutate(Year = Year+2000) %>%
    mutate(Year = as.factor(Year),
           Trait = as.factor(Trait)) %>%
    na.exclude() 

p2 <- ggplot(R2) +
    geom_density(aes(y = TraitValue, fill = Year, color = Year), alpha = 0.3, adjust = 2) +
    facet_wrap(~ Trait, scales = 'free', nrow = 1) +
    scale_fill_manual(values = c('grey30', 'red3')) +
    scale_color_manual(values = c('grey20', 'red4')) +
    # scale_y_log10() +
    theme_pubclean(base_size = 30) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          strip.text = element_text(size = 20, face = 'bold'),
          legend.position = 'bottom')+
    labs(y = 'Trait Value'); p2

df <- cbind.data.frame(n = rep(seq(n), each = 4),
                       trait = rep(c('LDMC', 'N','P','SLA'), times = n),
                       delta = rep(0),
                       pdelta = rep(0),
                       p = rep(0))

for (i in seq(n)) {
    y <- sample(dim(R2)[1], n)
    coef <- summary(lm(TraitValue ~ Year, R2 %>% filter(Trait == 'LDMC') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'LDMC'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'LDMC'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'LDMC'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R2 %>% filter(Trait == 'N') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'N'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'N'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'N'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R2 %>% filter(Trait == 'P') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'P'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'P'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'P'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R2 %>% filter(Trait == 'SLA')%>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'SLA'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'SLA'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'SLA'] <- coef[2,4]
}

df %>% group_by(trait) %>% summarise(mean(delta), sd(delta), mean(pdelta), sd(pdelta), median(p))

R2 %>% group_by(Trait, Year) %>%
    summarise(CV = sd(TraitValue)/mean(TraitValue)*100) %>%
    group_by(Trait) %>%
    mutate(CV_diff = (CV - lag(CV, default = first(CV)))/lag(CV, default = first(CV))*100) %>%
    filter(CV_diff != 0)






# Region 3 - 2016 and 2015
R3 <- cbind.data.frame(
    seq(ncell(ldmc15.r3)),
    ldmc15.r3 %>% as.data.frame(na.rm = F),
    ldmc16.r3 %>% as.data.frame(na.rm = F),
    
    n15.r3 %>% as.data.frame(na.rm = F),
    n16.r3 %>% as.data.frame(na.rm = F),
    
    p15.r3 %>% as.data.frame(na.rm = F),
    p16.r3 %>% as.data.frame(na.rm = F),
    
    sla15.r3 %>% as.data.frame(na.rm = F),
    sla16.r3 %>% as.data.frame(na.rm = F)
)
names(R3) <- c('ID','LDMC15', 'LDMC16', 'N15', 'N16', 'P15', 'P16', 'SLA15', 'SLA16')

R3 <- R3 %>% pivot_longer(-ID, names_to = 'TraitYear', values_to = 'TraitValue') %>%
    separate(TraitYear, into = c('Trait', 'Year'), sep = "(?<=\\D)(?=\\d)", convert = T) %>%
    mutate(Year = Year+2000) %>%
    group_by(Trait) %>%
    mutate(Year = as.factor(Year),
           Trait = as.factor(Trait)) %>%
    na.exclude

p3 <- ggplot(R3) +
    geom_density(aes(y = TraitValue, fill = Year, color = Year), alpha = 0.3, adjust = 2) +
    facet_wrap(~ Trait, scales = 'free', nrow = 1) +
    scale_fill_manual(values = c('grey30', 'goldenrod4')) +
    scale_color_manual(values = c('grey20', 'goldenrod4')) +
    # scale_y_log10() +
    theme_pubclean(base_size = 30) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          strip.text = element_text(size = 20, face = 'bold'),
          legend.position = 'bottom')+
    labs(y = 'Trait Value'); p3

df <- cbind.data.frame(n = rep(seq(n), each = 4),
                       trait = rep(c('LDMC', 'N','P','SLA'), times = n),
                       delta = rep(0),
                       pdelta = rep(0),
                       p = rep(0))

for (i in seq(n)) {
    y <- sample(dim(R3)[1], n)
    coef <- summary(lm(TraitValue ~ Year, R3 %>% filter(Trait == 'LDMC') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'LDMC'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'LDMC'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'LDMC'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R3 %>% filter(Trait == 'N') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'N'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'N'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'N'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R3 %>% filter(Trait == 'P') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'P'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'P'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'P'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R3 %>% filter(Trait == 'SLA') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'SLA'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'SLA'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'SLA'] <- coef[2,4]
}

df %>% group_by(trait) %>% summarise(mean(delta), sd(delta), mean(pdelta), sd(pdelta), median(p))

R3 %>% group_by(Trait, Year) %>%
    summarise(CV = sd(TraitValue)/mean(TraitValue)*100) %>%
    group_by(Trait) %>%
    mutate(CV_diff = (CV - lag(CV, default = first(CV)))/lag(CV, default = first(CV))*100) %>%
    filter(CV_diff != 0)





# Region 4 - 2016 and 2015
R4 <- cbind.data.frame(
    seq(ncell(ldmc15.r4)),
    ldmc15.r4 %>% as.data.frame(na.rm = F),
    ldmc16.r4 %>% as.data.frame(na.rm = F),
    
    n15.r4 %>% as.data.frame(na.rm = F),
    n16.r4 %>% as.data.frame(na.rm = F),
    
    p15.r4 %>% as.data.frame(na.rm = F),
    p16.r4 %>% as.data.frame(na.rm = F),
    
    sla15.r4 %>% as.data.frame(na.rm = F),
    sla16.r4 %>% as.data.frame(na.rm = F)
)
names(R4) <- c('ID','LDMC15', 'LDMC16', 'N15', 'N16', 'P15', 'P16', 'SLA15', 'SLA16')
R4 <- R4 %>% pivot_longer(-ID, names_to = 'TraitYear', values_to = 'TraitValue') %>%
    separate(TraitYear, into = c('Trait', 'Year'), sep = "(?<=\\D)(?=\\d)", convert = T) %>%
    mutate(Year = Year+2000) %>%
    group_by(Trait) %>%
    mutate(Year = as.factor(Year),
           Trait = as.factor(Trait)) %>%
    na.exclude() 

p4 <- ggplot(R4) +
    geom_density(aes(y = TraitValue, fill = Year, color = Year), alpha = 0.3, adjust = 1) +
    facet_wrap(~ Trait, scales = 'free', nrow = 1) +
    scale_fill_manual(values = c('grey30', 'goldenrod4')) +
    scale_color_manual(values = c('grey20', 'goldenrod4')) +
    # scale_y_log10() +
    theme_pubclean(base_size = 30) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          strip.text = element_text(size = 20, face = 'bold'),
          legend.position = 'bottom')+
    labs(y = 'Trait Value'); p4

df <- cbind.data.frame(n = rep(seq(n), each = 4),
                       trait = rep(c('LDMC', 'N','P','SLA'), times = n),
                       delta = rep(0),
                       pdelta = rep(0),
                       p = rep(0))
n <- 500
for (i in seq(n)) {
    
    coef <- summary(lm(TraitValue ~ Year, R4 %>% filter(Trait == 'LDMC') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'LDMC'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'LDMC'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'LDMC'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R4 %>% filter(Trait == 'N') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'N'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'N'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'N'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R4 %>% filter(Trait == 'P') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'P'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'P'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'P'] <- coef[2,4]
    
    coef <- summary(lm(TraitValue ~ Year, R4 %>% filter(Trait == 'SLA') %>% slice_sample(n=n)))$coefficients
    df$delta[df$n == i & df$trait == 'SLA'] <- coef[2,1]
    df$pdelta[df$n == i & df$trait == 'SLA'] <- coef[2,1]/coef[1,1] *100
    df$p[df$n == i & df$trait == 'SLA'] <- coef[2,4]
}

df %>% group_by(trait) %>% summarise(mean(delta), sd(delta), mean(pdelta), sd(pdelta), median(p))

R4 %>% group_by(Trait, Year) %>%
    summarise(CV = sd(TraitValue)/mean(TraitValue)*100) %>%
    group_by(Trait) %>%
    mutate(CV_diff = (CV - lag(CV, default = first(CV)))/lag(CV, default = first(CV))*100) %>%
    filter(CV_diff != 0)



figure(p1, path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/Region1_comp.png', save = T)
figure(p2, path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/Region2_comp.png', save = T)
figure(p3, path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/Region3_comp.png', save = T)
figure(p4, path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/Region4_comp.png', save = T)



# Spectra plots -----------------------------------------------------------

region1 <- (rtr22.r1-rtr15.r1) %>%
    as.data.frame %>% 
    mutate(ID = seq(n())) %>%
    pivot_longer(-ID, names_to = 'wl', values_to = 'Reflectance', names_transform = list(wl = as.integer)) %>%
    group_by(wl) %>%
    filter(Reflectance != 0) %>%
    summarise(mean = mean(Reflectance, na.rm = T), sd = sd(Reflectance, na.rm = T))

region1b <- (rtr22.r1-rtr15.r1) %>%
    as.data.frame %>% 
    mutate(ID = seq(n())) %>%
    pivot_longer(-ID, names_to = 'wl', values_to = 'Reflectance', names_transform = list(wl = as.integer)) %>%
    filter(Reflectance != 0) %>%
    filter(ID %in% sample(length(unique(ID)), 100))

p1 <- ggplot() +
    geom_line(data = region1b, aes(wl, Reflectance, group = ID), alpha = 0.1, color = 'red3') +
    geom_ribbon(data = region1, aes(wl, ymin = mean-sd, ymax = mean+sd), alpha = 0.2)+
    geom_line(data = region1, aes(wl, mean), linewidth = 1) +
    geom_hline(yintercept = 0) +
    theme_pubr(base_size = 25) +
    coord_cartesian(ylim = c(-0.06, 0.06)) +
    scale_y_continuous(breaks = c(-0.05, 0, 0.05)) +
    labs(x = 'Wavelength (nm)', y = bquote(Delta~'Reflectance')); p1

region2 <- (rtr22.r2-rtr15.r2) %>%
    as.data.frame %>% 
    mutate(ID = seq(n())) %>%
    pivot_longer(-ID, names_to = 'wl', values_to = 'Reflectance', names_transform = list(wl = as.integer)) %>%
    group_by(wl) %>%
    filter(Reflectance != 0) %>%
    summarise(mean = mean(Reflectance, na.rm = T), sd = sd(Reflectance, na.rm = T))
region2b <- (rtr22.r2-rtr15.r2) %>%
    as.data.frame %>% 
    mutate(ID = seq(n())) %>%
    pivot_longer(-ID, names_to = 'wl', values_to = 'Reflectance', names_transform = list(wl = as.integer)) %>%
    filter(Reflectance != 0) %>%
    filter(ID %in% sample(length(unique(ID)), 100))
p2 <- ggplot() +
    geom_line(data = region2b, aes(wl, Reflectance, group = ID), alpha = 0.1, color = 'red3') +
    geom_ribbon(data = region2, aes(wl, ymin = mean-sd, ymax = mean+sd), alpha = 0.2)+
    geom_line(data = region2, aes(wl, mean), linewidth = 1) +
    geom_hline(yintercept = 0) +
    theme_pubr(base_size = 25) +
    coord_cartesian(ylim = c(-0.06, 0.06)) +
    scale_y_continuous(breaks = c(-0.05, 0, 0.05)) +
    labs(x = 'Wavelength (nm)', y = bquote(Delta~'Reflectance')); p2

region3 <- (rtr16.r3-rtr15.r3) %>%
    as.data.frame %>% 
    mutate(ID = seq(n())) %>%
    pivot_longer(-ID, names_to = 'wl', values_to = 'Reflectance', names_transform = list(wl = as.integer)) %>%
    group_by(wl) %>%
    filter(Reflectance != 0) %>%
    summarise(mean = mean(Reflectance, na.rm = T), sd = sd(Reflectance, na.rm = T))
region3b <- (rtr16.r3-rtr15.r3) %>%
    as.data.frame %>% 
    mutate(ID = seq(n())) %>%
    pivot_longer(-ID, names_to = 'wl', values_to = 'Reflectance', names_transform = list(wl = as.integer)) %>%
    filter(Reflectance != 0) %>%
    filter(ID %in% sample(length(unique(ID)), 100))
p3 <- ggplot() +
    geom_line(data = region3b, aes(wl, Reflectance, group = ID), alpha = 0.1, color = 'red3') +
    geom_ribbon(data = region3, aes(wl, ymin = mean-sd, ymax = mean+sd), alpha = 0.2)+
    geom_line(data = region3, aes(wl, mean), linewidth = 1) +
    geom_hline(yintercept = 0) +
    theme_pubr(base_size = 25) +
    coord_cartesian(ylim = c(-0.06, 0.06)) +
    scale_y_continuous(breaks = c(-0.05, 0, 0.05)) +
    labs(x = 'Wavelength (nm)', y = bquote(Delta~'Reflectance')); p3

region4 <- (rtr16.r4-rtr15.r4) %>%
    as.data.frame %>% 
    mutate(ID = seq(n())) %>%
    pivot_longer(-ID, names_to = 'wl', values_to = 'Reflectance', names_transform = list(wl = as.integer)) %>%
    group_by(wl) %>%
    filter(Reflectance != 0) %>%
    mutate(mean = mean(Reflectance, na.rm = T), sd = sd(Reflectance, na.rm = T))
region4b <- (rtr16.r4-rtr15.r4) %>%
    as.data.frame %>% 
    mutate(ID = seq(n())) %>%
    pivot_longer(-ID, names_to = 'wl', values_to = 'Reflectance', names_transform = list(wl = as.integer)) %>%
    filter(Reflectance != 0) %>%
    filter(ID %in% sample(length(unique(ID)), 100))
p4 <- ggplot() +
    geom_line(data = region4b, aes(wl, Reflectance, group = ID), alpha = 0.1, color = 'red3') +
    geom_ribbon(data = region4, aes(wl, ymin = mean-sd, ymax = mean+sd), alpha = 0.2)+
    geom_line(data = region4, aes(wl, mean), linewidth = 1) +
    geom_hline(yintercept = 0) +
    theme_pubr(base_size = 25) +
    coord_cartesian(ylim = c(-0.06, 0.06)) +
    scale_y_continuous(breaks = c(-0.05, 0, 0.05)) +
    labs(x = 'Wavelength (nm)', y = bquote(Delta~'Reflectance')); p4


# rtr15df.r1 <- rtr15.r1 %>%
#     as.data.frame %>% 
#     mutate(ID = seq(n())) %>%
#     slice_sample(n = 100) %>%
#     pivot_longer(-ID) %>%
#     mutate(Region = rep('Region1_2015')) %>% 
#     select(-ID)
# 
# rtr15df.r2 <- rtr15.r2 %>%
#     as.data.frame %>% 
#     mutate(ID = seq(n())) %>%
#     slice_sample(n = 100) %>%
#     pivot_longer(-ID)%>%
#     mutate(Region = rep('Region2_2015')) %>% 
#     select(-ID)
# rtr15df.r3 <- rtr15.r3 %>%
#     as.data.frame %>% 
#     mutate(ID = seq(n())) %>%
#     slice_sample(n = 100) %>%
#     pivot_longer(-ID)%>%
#     mutate(Region = rep('Region3_2015')) %>% 
#     select(-ID)
# rtr15df.r4 <- rtr15.r4 %>%
#     as.data.frame %>% 
#     mutate(ID = seq(n())) %>%
#     slice_sample(n = 100) %>%
#     pivot_longer(-ID)%>%
#     mutate(Region = rep('Region4_2015')) %>% 
#     select(-ID)
# 
# rtr16df.r3 <- rtr16.r3 %>%
#     as.data.frame %>% 
#     mutate(ID = seq(n())) %>%
#     slice_sample(n = 100) %>%
#     pivot_longer(-ID)%>%
#     mutate(Region = rep('Region3_2016')) %>% 
#     select(-ID)
# rtr16df.r4 <- rtr16.r4 %>%
#     as.data.frame %>% 
#     mutate(ID = seq(n())) %>%
#     slice_sample(n = 100) %>%
#     pivot_longer(-ID)%>%
#     mutate(Region = rep('Region4_2016')) %>% 
#     select(-ID)
# 
# rtr22df.r1 <- rtr22.r1 %>%
#     as.data.frame %>% 
#     mutate(ID = seq(n())) %>%
#     slice_sample(n = 100) %>%
#     pivot_longer(-ID)%>%
#     mutate(Region = rep('Region1_2022')) %>% 
#     select(-ID)
# 
# rtr22df.r2 <- rtr22.r2 %>%
#     as.data.frame %>% 
#     mutate(ID = seq(n())) %>%
#     slice_sample(n = 100) %>%
#     pivot_longer(-ID)%>%
#     mutate(Region = rep('Region2_2022')) %>% 
#     select(-ID)
# 
# rtr.df <- rbind.data.frame(rtr15df.r1,
#                            rtr15df.r2,
#                            rtr15df.r3,
#                            rtr15df.r4,
#                            rtr16df.r3,
#                            rtr16df.r4,
#                            rtr22df.r1,
#                            rtr22df.r2)
# 
# names(rtr.df) <- c('Wavelength', 'Reflectance', 'RegionYear')
# 
# rtr.df <- rtr.df %>%
#     separate(RegionYear, into = c('Region', 'Year'), sep = "_", convert = T) %>%
#     na.exclude() %>%
#     mutate(Wavelength = as.numeric(Wavelength),
#            Year = as.factor(Year),
#            Region = as.factor(Region))
# 
# R1.df <- rtr.df %>% filter(Region == 'Region1')
# R2.df <- rtr.df %>% filter(Region == 'Region2')
# R3.df <- rtr.df %>% filter(Region == 'Region3')
# R4.df <- rtr.df %>% filter(Region == 'Region4')
# 
# R1.dif <- R1.df %>% group_by(Wavelength, Year) %>%
#     summarise(Avg = mean(Reflectance), sd = sd(Reflectance)) %>%
#     pivot_wider(id_cols = Wavelength, names_from = Year, values_from = c(Avg, sd)) %>%
#     mutate(dif = Avg_2022-Avg_2015,
#            upr = dif+sd_2022,
#            lwr = dif-sd_2022)
# R2.dif <- R2.df %>% group_by(Wavelength, Year) %>%
#     summarise(Avg = mean(Reflectance), sd = sd(Reflectance)) %>%
#     pivot_wider(id_cols = Wavelength, names_from = Year, values_from = c(Avg, sd)) %>%
#     mutate(dif = Avg_2022-Avg_2015,
#            upr = dif+sd_2022,
#            lwr = dif-sd_2022)
# R3.dif <- R3.df %>% group_by(Wavelength, Year) %>%
#     summarise(Avg = mean(Reflectance), sd = sd(Reflectance)) %>%
#     pivot_wider(id_cols = Wavelength, names_from = Year, values_from = c(Avg, sd)) %>%
#     mutate(dif = Avg_2016-Avg_2015,
#            upr = dif+sd_2016,
#            lwr = dif-sd_2016)
# R4.dif <- R4.df %>% group_by(Wavelength, Year) %>%
#     summarise(Avg = mean(Reflectance), sd = sd(Reflectance, na.rm = T)) %>%
#     pivot_wider(id_cols = Wavelength, names_from = Year, values_from = c(Avg, sd)) %>%
#     mutate(dif = Avg_2016-Avg_2015,
#            upr = dif+sd_2016,
#            lwr = dif-sd_2016)


# p1 <- ggplot(region1) +
#     geom_ribbon(aes(Wavelength, ymin = lwr, ymax = upr), alpha = 0.2)+
#     geom_line(aes(Wavelength, dif), linewidth = 1) +
#     geom_hline(yintercept = 0) +
#     theme_pubr(base_size = 25) +
#     labs(y = bquote(Delta~'Reflectance')); p1
# 
# p2 <- ggplot(R2.dif) +
#     geom_ribbon(aes(Wavelength, ymin = lwr, ymax = upr), alpha = 0.2)+
#     geom_line(aes(Wavelength, dif), linewidth = 1) +
#     geom_hline(yintercept = 0) +
#     theme_pubr(base_size = 25) +
#     labs(y = bquote(Delta~'Reflectance')); p2
# 
# p3 <- ggplot(R3.dif) +
#     geom_ribbon(aes(Wavelength, ymin = lwr, ymax = upr), alpha = 0.2)+
#     geom_line(aes(Wavelength, dif), linewidth = 1) +
#     geom_hline(yintercept = 0) +
#     theme_pubr(base_size = 25) +
#     labs(y = bquote(Delta~'Reflectance')); p3
# 
# p4 <- ggplot(R4.dif) +
#     geom_ribbon(aes(Wavelength, ymin = lwr, ymax = upr), alpha = 0.2)+
#     geom_line(aes(Wavelength, dif), linewidth = 1) +
#     geom_hline(yintercept = 0) +
#     theme_pubr(base_size = 25) +
#     labs(y = bquote(Delta~'Reflectance')); p4

# p1 <- ggplot(R1.df) +
#     geom_smooth(aes(Wavelength, Reflectance, color = Year), method = 'loess', span = 0.05, se = F) +
#     scale_color_manual(values = c('grey30', 'red3')) +
#     theme_minimal(base_size = 30)  +
#     theme(legend.position = 'none'); p1
# p2 <- ggplot(R2.df) +
#     geom_smooth(aes(Wavelength, Reflectance, color = Year), method = 'loess', span = 0.05, se = F) +
#     scale_color_manual(values = c('grey30', 'red3')) +
#     theme_minimal(base_size = 30)  +
#     theme(legend.position = 'none'); p2
# p3 <- ggplot(R3.df) +
#     geom_smooth(aes(Wavelength, Reflectance, color = Year), method = 'loess', span = 0.05, se = F) +
#     scale_color_manual(values = c('grey30', 'goldenrod2')) +
#     theme_minimal(base_size = 30)  +
#     theme(legend.position = 'none'); p3
# p4 <- ggplot(R4.df) +
#     geom_smooth(aes(Wavelength, Reflectance, color = Year), method = 'loess', span = 0.05, se = F) +
#     scale_color_manual(values = c('grey30', 'goldenrod2')) +
#     theme_minimal(base_size = 30) +
#     theme(legend.position = 'none'); p4


figure(p1, path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/Region1_comp_spectra.png', height = 3.5, save = T)
figure(p2, path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/Region2_comp_spectra.png', height = 3.5, save = T)
figure(p3, path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/Region3_comp_spectra.png', height = 3.5, save = T)
figure(p4, path.name = '~/Current Projects/SBG/Currey_IS_ClimateExtremes_OSSE/Figures/Region4_comp_spectra.png', height = 3.5, save = T)
