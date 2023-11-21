source('~/R/clean.R')
library(ncdf4)
library(terra)
library(tidyverse)
library(tidyterra)


dp <- '~/Current Projects/SBG/LPJ/Reflectance_Data/version2/'
shp <- vect('~/Current Projects/SBG/LPJ/Misc data/Shapefile/Would_boundaries_noGreenland.shp')

rfl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2015.nc')) %>%
    ncvar_get('DR', start = c(1,1,2,8), count = c(-1,-1,210,1)) %>%
    aperm(c(2,1,3)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) %>%
    as.array()

rdn <- nc_open(file.path(dp, 'lpj-prosail_levelD_TOA-radiance_Version021_m_2015.nc')) %>%
    ncvar_get('Radiance', start = c(2,1,1,8), count = c(210,-1,-1,1)) %>% 
    aperm(c(2,3,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) %>%
    as.array()

rtr <- nc_open(file.path(dp, 'lpj-prosail_levelE_retrieved-HDR_Version021_m_2015.nc')) %>%
    ncvar_get('HDR', start = c(2,1,1,8), count = c(210,-1,-1,1)) %>% 
    aperm(c(2,3,1)) %>%
    rast(ext = c(-180, 180, -90, 90), crs = 'EPSG:4326') %>%
    mask(shp) %>%
    as.array()

lats <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2015.nc')) %>%
    ncvar_get('lat')
lons <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2015.nc')) %>%
    ncvar_get('lon')
wl <- nc_open(file.path(dp, 'lpj-prosail_levelC_DR_Version021_m_2015.nc')) %>%
    ncvar_get('wl', start = c(2), count = c(210))
    
rfl[rfl==0] <- -9999
rdn[rdn==0] <- -9999
rtr[rtr==0] <- -9999

# Create dimensions
lon_dim <- ncdim_def("lon", "degrees_east", lons, longname = 'Longitude')
lat_dim <- ncdim_def("lat", "degrees_north", lats, longname = 'Latitude')
wl_dim <- ncdim_def("wl", paste0("nanometers"), wl, longname = 'Wavelength')



# Reflectance -------------------------------------------------------------

# Define a variable
varInput <- "reflectance"
varname_nc <- "DR"
varUnit <- "unitless"
varLongName <- 'Top-of-Canopy Directional Reflectance'
nc_var <- ncvar_def(varname_nc, varUnit, 
                    list(lat_dim, lon_dim, wl_dim), -9999,
                    longname = varLongName)

nc_name <- file.path(dp, 'Reflectance_8_2015.nc')
nc_out <- nc_create(nc_name, nc_var)

ncvar_put(nc_out, nc_var, rfl)
nc_close(nc_out)


# Radiance ----------------------------------------------------------------

# Define a variable
varInput <- "Radiance_"
varname_nc <- "Radiance"
varUnit <- "w sr-1 m-2"
varLongName <- 'Top-of-Atmosphere Radiance'
nc_var <- ncvar_def(varname_nc, varUnit, 
                    list(lat_dim, lon_dim, wl_dim), -9999,
                    longname = varLongName)

nc_name <- file.path(dp, 'Radiance_8_2015.nc')
nc_out <- nc_create(nc_name, nc_var)

ncvar_put(nc_out, nc_var, rdn)
nc_close(nc_out)


# Retreivals --------------------------------------------------------------


# Define a variable
varInput <- "retrieved_reflectance_"
varname_nc <- "HDR"
varUnit <- "unitless"
varLongName <- 'Top-of-Canopy HDR Reflectance'
nc_var <- ncvar_def(varname_nc, varUnit, 
                    list(lat_dim, lon_dim, wl_dim), -9999,
                    longname = varLongName)

nc_name <- file.path(dp,'Retrievals_8_2015.nc')
nc_out <- nc_create(nc_name, nc_var)

ncvar_put(nc_out, nc_var, rtr)
nc_close(nc_out)

