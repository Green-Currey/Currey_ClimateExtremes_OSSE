global_spectra_comparison = function(
        array1, 
        array2 = NULL, 
        array1_short_name = 'array1',
        array2_short_name = 'array2',
        percent_mask = 0, # in percent
        nc.outpath = 'c:/Users/bryce/OneDrive/Documents/Current Projects/SBG/LPJ/Comparisons', 
        # mode = 'absolute', # absolute or percent
        month # which month to compare. Empty for all months
) {
    
    bn <- function(r) {
        terra::app(r, function(x) x / sqrt(sum(x^2)))
    }
    
    library(terra)
    library(ncdf4)
    library(tidyverse)
    library(neuralnet)
    library(readr)
    
    message('Starting comparison.')
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Check and reorder the dimensions for raster calcs  #
    # dimensions must be lat, lon, wl, time (optional)   #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    
    array1 <- as.array(array1); 
    if (!is.null(array2)) {array2 <- as.array(array2)}
    
    d1 <- dim(array1)
    d2 <- dim(array2)
    
    # array 1
    if (length(d1)==4) {
        array1 <- aperm(array1, match(dim(array1), c(360,720,210,12)))
    } else {
        array1 <- aperm(array1, match(dim(array1), c(360,720,210)))
    }
    
    # array 2
    if (!is.null(array2)) {
        message('Matching arrays.')
        
        if (length(d2)==4) {
            array2 <- aperm(array2, match(dim(array2), c(360,720,210,12)))
        } else {
            array2 <- aperm(array2, match(dim(array2), c(360,720,210)))
        }
        
        
        # if arrays have differing dimensions, subset month to match 
        if (length(d2) != length(d1)) {
            if (length(d2)==4) {
                array1 <- array1
                array2 <- array2[,,,month]
            } else if (length(d1)==4) {
                array1 <- array1[,,,month]
                array2 <- array2
            }
        }
        
        # get new dimensions
        d1 <- dim(array1)
        d2 <- dim(array2)
        
        try(expr = setequal(d1, d2), stop('Incorrect dimensions'))
        
        message('Arrays matched.')
        
    }
    
    
    # possible conditions
    if (length(d1)==4 & !is.numeric(month)) { allmonths <- TRUE; specificmonth <- FALSE }
    if (is.numeric(month)) { allmonths <- FALSE; specificmonth <- TRUE }    
    if (length(d1)==3 & !is.numeric(month)) { allmonths <- FALSE; specificmonth <- FALSE }
    if (length(month) == 2) {month_to_month <- TRUE; allmonths <- FALSE; specificmonth <- FALSE }
    
    
    # NC dimensions and variables
    bands <- seq(410,2500,10)
    
    lon <- ncdim_def( "lon", "degrees_east", seq(-179.75, 179.75, 0.5), create_dimvar = T)
    lat <- ncdim_def( "lat", "degrees_north", seq(89.75, -89.75, -0.5), create_dimvar = T)
    tmonth <- ncdim_def( "Month", "Month", seq(12), create_dimvar = T)
    stats <- ncdim_def('Statistics', 'MaxDif/Wl/Mean/SD/CV/RMSE/%RMSE/MAE/%MAE/%Dif', seq(10), create_dimvar = T)
    dims <- c(720,360,10,12)
    
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    #    if array are four, loop through entire array    #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    if (allmonths) {
        message('Comparing all months.')
        pb <- txtProgressBar(min = 0, max = 12, style = 3)
        for (m in seq(12)) {
            
            array1m <- array1[,,,m]
            array2m <- array2[,,,m]
            
            # ~~~~~~~~~~~~~~~ #
            #  STAT ANALYSES  #
            # ~~~~~~~~~~~~~~~ #
            
            r1 <- rast(array1m, crs = 'EPSG:4326', ext = c(-180,180,-90,90))
            names(r1) <- bands
            remove(array1m)
            
            r2 <- rast(array2m, crs = 'EPSG:4326', ext = c(-180,180,-90,90))
            names(r2) <- bands
            remove(array1m)
            
            # Stats
            dif <- r2-r1
            mean_dif <- mean(dif, na.rm = T)
            sd_dif <- sqrt(sum((dif - mean_dif)^2, na.rm = T)/210)
            
            cv <- sd_dif/abs(mean_dif)
            upr <- global(cv, function(x) quantile(x,0.95,na.rm=T)) %>% as.numeric()
            lwr <- global(cv, function(x) quantile(x,0.05,na.rm=T)) %>% as.numeric()
            cv[cv>upr] <- NA
            cv[cv<lwr] <- NA
            
            pdif <- dif/r1*100
            pdif[is.infinite(pdif)] <- NA; #pdif
            pdif[pdif>50] <- NA; #pdif
            pdif[pdif< -50] <- NA; #pdif
            mean_pdif <- mean(pdif, na.rm = T)
            
            rmse <- sqrt(mean(dif^2, na.rm = T))
            mae <- mean(abs(dif), na.rm = T)
            prmse <- rmse/mean(r1)*100
            pmae <- mae/mean(r1)*100
            
            for (i in as.character(bands)) {
                dif.i <- dif %>% subset(i)
                if (i == '410') {
                    max_dif <- dif.i
                    wl_dif <- max_dif
                    wl_dif[!is.na(wl_dif)] <- as.numeric(i)
                } else {
                    gt <- abs(dif.i) >= abs(max_dif)
                    max_dif[gt] <- dif.i
                    wl_dif[gt] <- as.numeric(i)
                }
            } # end band loop
            
            # max_dif[max_dif<percent_mask & max_dif > percent_mask*-1] <- NA
            wl_dif <- wl_dif %>% mask(max_dif)
            
            # ~~~~~~~~~~~~~~~~~~~~ #
            #  END STATS ANALYSES  #
            # ~~~~~~~~~~~~~~~~~~~~ #
            
            if (m == 1) {stat_array <- array(dim = dims)}
            
            stat_array[,,1,m] <- as.array(t(max_dif))
            stat_array[,,2,m] <- as.array(t(wl_dif))
            stat_array[,,3,m] <- as.array(t(mean_dif))
            stat_array[,,4,m] <- as.array(t(sd_dif))
            stat_array[,,5,m] <- as.array(t(cv))
            stat_array[,,6,m] <- as.array(t(rmse))
            stat_array[,,7,m] <- as.array(t(prmse))
            stat_array[,,8,m] <- as.array(t(mae))
            stat_array[,,9,m] <- as.array(t(pmae))
            stat_array[,,10,m] <- as.array(t(mean_pdif))
            setTxtProgressBar(pb, m)
            
        } # end monthly loop
        
    } else {
        
        
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # If month is specified numerically #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        if (specificmonth) {
            message(paste0('Comparing month ',month,' arrays.'))
            if(length(d1)==4) {
                array1m <- array1[,,,month]
                array2m <- array2[,,,month]
            }else{
                array1m <- array1
                array2m <- array2
            }
        } else if (month_to_month){
            message(paste0('Comparing month ',month[1],' with month ', month[2],' from array 1.'))
            message(paste0('Ignoring array 2.'))
            array1m <- array1[,,,month[1]]
            array2m <- array1[,,,month[2]]
            
            
        } else {
            message('Comparing two arrays, unspecifed month.')
            array1m <- array1
            array2m <- array2
        }
        
        # ~~~~~~~~~~~~~~~ #
        #  STAT ANALYSES  #
        # ~~~~~~~~~~~~~~~ #
        
        r1 <- rast(array1m, crs = 'EPSG:4326', ext = c(-180,180,-90,90))
        r1[r1==0] <- NA
        names(r1) <- bands
        remove(array1m)
        
        r2 <- rast(array2m, crs = 'EPSG:4326', ext = c(-180,180,-90,90))
        r2[r2==0] <- NA
        names(r2) <- bands
        remove(array2m)
        
        # Stats
        dif <- r2-r1
        mean_dif <- mean(dif, na.rm = T)
        sd_dif <- sqrt(sum((dif - mean_dif)^2, na.rm = T)/210)
        
        cv <- sd_dif/abs(mean_dif)
        upr <- global(cv, function(x) quantile(x,0.95,na.rm=T)) %>% as.numeric()
        lwr <- global(cv, function(x) quantile(x,0.05,na.rm=T)) %>% as.numeric()
        cv[cv>upr] <- NA
        cv[cv<lwr] <- NA
        
        pdif <- abs(dif)/r1*100
        pdif[is.infinite(pdif)] <- NA; #pdif
        pdif[pdif>50] <- NA; #pdif
        pdif[pdif< -50] <- NA; #pdif
        mean_pdif <- mean(pdif, na.rm = T)
        
        rmse <- sqrt(mean(dif^2, na.rm = T))
        mae <- mean(abs(dif), na.rm = T)
        prmse <- rmse/mean(r1)*100
        pmae <- mae/mean(r1)*100

        remove(r1); remove(r2); 
        for (i in as.character(bands)) {
            dif.i <- dif %>% subset(i)
            if (i == '410') {
                max_dif <- dif.i
                wl_dif <- max_dif
                wl_dif[!is.na(wl_dif)] <- as.numeric(i)
            } else {
                gt <- abs(dif.i) >= abs(max_dif)
                max_dif[gt] <- dif.i
                wl_dif[gt] <- as.numeric(i)
            }
        } # end band loop
        
        # max_dif[max_dif<percent_mask & max_dif > percent_mask*-1] <- NA
        wl_dif <- wl_dif %>% mask(max_dif)
        
        # ~~~~~~~~~~~~~~~~~~~~ #
        #  END STATS ANALYSES  #
        # ~~~~~~~~~~~~~~~~~~~~ #
        
        stat_array <- array(dim = dims[1:3])
        
        stat_array[,,1] <- as.array(t(max_dif))
        stat_array[,,2] <- as.array(t(wl_dif))
        stat_array[,,3] <- as.array(t(mean_dif))
        stat_array[,,4] <- as.array(t(sd_dif))
        stat_array[,,5] <- as.array(t(cv))
        stat_array[,,6] <- as.array(t(rmse))
        stat_array[,,7] <- as.array(t(prmse))
        stat_array[,,8] <- as.array(t(mae))
        stat_array[,,9] <- as.array(t(pmae))
        stat_array[,,10] <- as.array(t(mean_pdif))
        stat_array[is.na(stat_array)] <- -9999
        
    } # end monthly logical test
    
    remove(array1); remove(array2);
    remove(cv); remove(dif); remove(dif.i); remove(max_dif); remove(mean_dif); remove(sd_dif); remove(wl_dif); remove(gt); remove(rmse);
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    #                    Create NCDF                     #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    
    message('Creating ncdf.')
    var.name <- paste0('Stats')
    
    if (allmonths) {
        nc.name <- paste0("Dif_", array2_short_name, '_', array1_short_name, '.nc')
        
        varNC <- ncvar_def(name = var.name, 
                           units = 'stats',
                           dims = list(lon,lat,stats,tmonth),
                           missval = -9999,
                           prec = "double")
    } else {
        if (specificmonth) {
            nc.name <- paste0('Dif_', array2_short_name, '_', array1_short_name, '_month_', month, '.nc')
            
        } else if (month_to_month) {
            nc.name <- paste0('Dif_months_',month[2], '_and_',month[1], '_', array1_short_name, '.nc')        
            
        } else {
            nc.name <- paste0("Dif_", array2_short_name, '_', array1_short_name, '.nc')        
        }
        
        varNC <- ncvar_def(name = var.name, 
                           units = 'stats',
                           dim = list(lon,lat,stats),
                           missval = -9999,
                           prec = "double")
    }
    
    
    ncnew <- nc_create(filename = file.path(nc.outpath, nc.name),
                       vars = varNC)
    
    ncvar_put(nc = ncnew,
              varid = varNC,
              vals = stat_array)
    
    nc_close(ncnew)
    message('ncdf: ')
    message(nc.outpath)
    message('ouput: ')
    message(nc.name)
    
}

