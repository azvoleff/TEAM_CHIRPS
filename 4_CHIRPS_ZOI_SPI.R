###############################################################################
# Calculates mean precipitation within ZOI of team site from CHIRPS pentad or 
# monthly data.
###############################################################################

source('0_settings.R')

stopifnot(dataset == 'monthly')

library(rgdal)
library(raster)
library(rgeos)
library(lubridate)
library(dplyr)
library(SPEI)

library(doParallel)
library(foreach)
n_cpus <- 8

cl  <- makeCluster(n_cpus)
registerDoParallel(cl)

spi_periods <- c(6, 12, 24, 36)

# Define function to calculate SPI
calc_spi <- function(chirps_mat, spi_period) {
    # Split the chirps_mat into pieces to minimize the number of calls to the
    # spi function
    start_n <- floor(seq(1, ncol(chirps_mat),
                         length.out=min(ncol(chirps_mat), n_cpus + 1)))
    end_n <- start_n[2:length(start_n)]
    start_n <- start_n[1:(length(start_n) - 1)]
    end_n <- end_n - 1
    end_n[length(end_n)] <- end_n[length(end_n)] + 1
    spi_mat <- foreach(start_n=start_n, end_n=end_n,
                       .combine=cbind, .packages=c("SPEI")) %dopar% {
        # Multiply by 1000 and round so results can be stored as INT2S
        round(spi(chirps_mat[, start_n:end_n], spi_period, na.rm=TRUE)$fitted * 1000)
    }
    return(spi_mat)
}

foreach (sitecode=sitecodes) %do% {
    message(sitecode)
    base_name <- file.path(out_folder,
                           paste0(sitecode, '_CHIRPS_', dataset,
                                  '_', start_date, '-', end_date))
    chirps_tif <- paste0(base_name, '.tif')
    stopifnot(file_test('-f', chirps_tif))
    chirps <- brick(chirps_tif)

    chirps <- calc(chirps, function(vals) {
        vals[vals == chirps_NA_value] <- NA
        return(vals)
    })

    chirps_mat <- t(as.matrix(chirps))

    for (spi_period in spi_periods) {
        spi_mat <- calc_spi(chirps_mat, spi_period)
        out_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
        out_rast <- setValues(out_rast, t(spi_mat))
        spi_filename <- paste0(base_name, '_SPI_', spi_period, '.tif')
        out_rast <- writeRaster(out_rast, spi_filename, overwrite=TRUE,
                                datatype="INT2S")
    }
}

stopCluster(cl)
