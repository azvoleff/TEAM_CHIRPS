###############################################################################
# Calculates mean precipitation within ZOI of team site from CHIRPS pentad or 
# monthly data.
###############################################################################

library(rgdal)
library(raster)
library(rgeos)
library(lubridate)
library(dplyr)
library(SPEI)

library(doParallel)
library(foreach)
n_cpus <- 4

registerDoParallel(n_cpus)

overwrite <- TRUE

sites <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')
sitecodes <- sites$sitecode

product <- 'v1p8chirps'
chirps_NA_value <- -9999
dataset <- 'monthly'
date_limits_string <- '198101-201404'
# Note the below code is INCLUSIVE of the start date
chirps_start_date <- as.Date('1981/1/1')
# Note the below code is EXCLUSIVE of the end date
chirps_end_date <- as.Date('2014/5/1')

zoi_folder <- 'H:/Data/TEAM/ZOIs'
in_base_dir <- 'H:/Data/CHIRPS'
out_base_dir <- 'H:/Data/CHIRPS'
in_folder <- file.path(in_base_dir, paste0('global_', dataset))
out_folder <- file.path(out_base_dir, paste0('TEAM_', dataset))
stopifnot(file_test('-d', in_folder))
stopifnot(file_test('-d', out_folder))

yrs <- seq(year(chirps_start_date), year(chirps_end_date))
if (dataset == 'monthly') {
    dates <- seq(chirps_start_date, chirps_end_date, by='months')
    dates <- dates[dates < chirps_end_date]
    num_periods <- 12
} else if (dataset == 'pentad') {
    # Have 72 pentads per year (12 months per year, 6 pentads per month)
    yrs_rep <- rep(yrs, each=12*6)
    days <- c(1, 6, 11, 16, 21, 26)
    mths <- rep(paste(rep(seq(1, 12), each=6), days, sep='/'), length(yrs))
    dates <- as.Date(paste(yrs_rep, mths, sep='/'))
    dates <- dates[dates < chirps_end_date]
    num_periods <- 72
}

# This is the projection of the CHIRPS files, read from the .hdr files 
# accompanying the data
s_srs <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'

foreach(sitecode=iter(sitecodes), .inorder=FALSE,
        .packages=c('dplyr', 'raster', 'rgdal', 'sp', 'rgeos',
                    'SPEI')) %dopar% {
    zoi_file <- dir(zoi_folder, pattern=paste0('^ZOI_', sitecode, '_[0-9]{4}.RData'), full.names=TRUE)
    stopifnot(length(zoi_file) == 1)
    load(zoi_file)
    zoi <- spTransform(zoi, CRS(s_srs))

    chirps_file <- file.path(in_folder,
                          paste0(product, '_', dataset, '_', sitecode, '_', 
                                 date_limits_string, '.dat'))
    stopifnot(file_test('-f', chirps_file))

    chirps <- brick(chirps_file)
    chirps <- calc(chirps, function(vals) {
        vals[vals == chirps_NA_value] <- NA
        return(vals)
    })
    chirps <- mask(chirps, zoi)

    calc_spi <- function(chirps_mat, scale) {
        spi_mat <- spi(ppt_layers_in_cols, scale, na.rm=TRUE)$fitted
        spi_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
        spi_rast <- setValues(spi_rast, t(spi_mat))
        spi_filename <- file.path(out_folder,
                                  paste0(dataset, '_', sitecode, '_SPI_', 
                                         scale, '.tif'))
        writeRaster(spi_rast, spi_filename, overwrite=overwrite)
    }

    # Now calculate the SPI
    if (dataset == 'monthly') {
        ppt_layers_in_cols <- t(as.matrix(chirps))
        calc_spi(ppt_layers_in_cols, 1)
        calc_spi(ppt_layers_in_cols, 6)
        calc_spi(ppt_layers_in_cols, 12)
        calc_spi(ppt_layers_in_cols, 24)
        calc_spi(ppt_layers_in_cols, 36)
    }
}
