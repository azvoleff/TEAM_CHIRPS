###############################################################################
# Calculates mean precipitation within ZOI of team site from CHIRPS pentad or 
# monthly data.
###############################################################################

library(rgdal)
library(raster)
library(rgeos)
library(lubridate)
library(dplyr)

library(doParallel)
library(foreach)
n_cpus <- 4

registerDoParallel(n_cpus)

overwrite <- TRUE

sites <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')
sitecodes <- sites$sitecode

product <- 'v1p8chirps'
chirps_NA_value <- -9999
#dataset <- 'pentad'
#date_limits_string <- '198101-201424'
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

zoi_ppt <- foreach(sitecode=iter(sitecodes), .inorder=FALSE,
                   .packages=c('dplyr', 'raster', 'rgdal', 'sp', 'rgeos'),
                   .combine=rbind) %dopar% {
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

    # Setup a mask of NAs for all areas outside the ZOI
    chirps_mask <- is.na(chirps[[1]])
    chirps_mask[is.na(chirps[[1]])] <- NA
    chirps_mask[!is.na(chirps[[1]])] <- 1

    # Setup a dataframe with the precipitation data so anomalies, etc can be 
    # calculated
    years <- rep(yrs, each=num_periods)[1:nlayers(chirps)]
    years_rep <- rep(years, each=nrow(chirps)*ncol(chirps))
    subyears <- rep(seq(1, num_periods),  length.out=nlayers(chirps))
    subyears_rep <- rep(subyears, each=nrow(chirps)*ncol(chirps))
    pixels_rep <- rep(seq(1:(nrow(chirps)*ncol(chirps))), nlayers(chirps))
    chirps_df <- data.frame(year=years_rep,
                            subyear=subyears_rep, 
                            pixel=pixels_rep,
                            ppt=as.vector(chirps),
                            inzoi=as.vector(!is.na(chirps_mask)))
    chirps_df <- tbl_df(chirps_df)

    # Add 12 month total precip column. Use sides=1 in filter to have it be 12 
    # month total as of that date
    chirps_df <- group_by(chirps_df, pixel) %>%
        arrange(year, subyear) %>%
        mutate(tot_12mth=as.numeric(stats::filter(ppt, rep(1/12, 12), sides=1)*12),
               tot_24mth=as.numeric(stats::filter(ppt, rep(1/24, 24), sides=1)*24))
    chirps_df$tot_12mth[!chirps_df$inzoi] <- NA
    chirps_df$tot_24mth[!chirps_df$inzoi] <- NA

    # First calculate the mean_monthly for each subyear for each pixel
    ppt_mean_monthly <- group_by(chirps_df, pixel, subyear) %>%
        summarize(mean_monthly=mean(ppt, na.rm=TRUE))
    # Use chirps raster as a template
    ppt_mean_monthly_rast <- brick(chirps, values=FALSE, nl=num_periods)
    ppt_mean_monthly_rast <- setValues(ppt_mean_monthly_rast,
                               matrix(ppt_mean_monthly$mean_monthly, 
                                      nrow=nrow(chirps)*ncol(chirps), 
                                      ncol=num_periods, byrow=TRUE))
    # Mask areas outside ZOI
    ppt_mean_monthly_rast <- ppt_mean_monthly_rast * chirps_mask
    writeRaster(ppt_mean_monthly_rast,
                filename=file.path(out_folder, paste0(dataset, "_", sitecode, '_ppt_mean_monthly.tif')), 
                overwrite=overwrite)

    # Calculate the mean annual precipitation for each pixel
    ppt_mean_12mth <- group_by(chirps_df, pixel, year) %>%
        summarize(total_annual=sum(ppt, na.rm=TRUE)) %>%
        group_by(pixel) %>%
        summarize(mean_12mth=mean(total_annual, na.rm=TRUE))
    # Use chirps raster as a template
    ppt_mean_12mth_rast <- brick(chirps, values=FALSE, nl=1)
    ppt_mean_12mth_rast <- setValues(ppt_mean_12mth_rast,
                               matrix(ppt_mean_12mth$mean_12mth, 
                                      nrow=nrow(chirps)*ncol(chirps), 
                                      ncol=1))
    # Mask areas outside ZOI
    ppt_mean_12mth_rast <- ppt_mean_12mth_rast * chirps_mask
    writeRaster(ppt_mean_12mth_rast,
                filename=file.path(out_folder, paste0(dataset, "_", sitecode, '_ppt_mean_12mth.tif')), 
                overwrite=overwrite)

    # Calculate the mean annual precipitation for each pixel
    ppt_mean_24mth <- ppt_mean_12mth
    ppt_mean_24mth$mean_24mth <- ppt_mean_24mth$mean_12mth * 2
    ppt_mean_24mth <- ppt_mean_24mth[names(ppt_mean_24mth) != "mean_12mth"]
    # Use chirps raster as a template
    ppt_mean_24mth_rast <- brick(chirps, values=FALSE, nl=1)
    ppt_mean_24mth_rast <- setValues(ppt_mean_24mth_rast,
                               matrix(ppt_mean_24mth$mean_24mth, 
                                      nrow=nrow(chirps)*ncol(chirps), 
                                      ncol=1))
    # Mask areas outside ZOI
    ppt_mean_24mth_rast <- ppt_mean_24mth_rast * chirps_mask
    writeRaster(ppt_mean_24mth_rast,
                filename=file.path(out_folder, paste0(dataset, "_", sitecode, '_ppt_mean_24mth.tif')), 
                overwrite=overwrite)

    # Calculate 12 and 24 month precipidation anomalies
    chirps_df$anom_12mth <- chirps_df$tot_12mth - ppt_mean_12mth$mean_12mth[match(chirps_df$pixel, ppt_mean_12mth$pixel)]
    chirps_df$anom_24mth <- chirps_df$tot_24mth - ppt_mean_24mth$mean_24mth[match(chirps_df$pixel, ppt_mean_24mth$pixel)]
    #filter(chirps_df, pixel == 1)[1:36,]
    save(chirps_df, file=file.path(out_folder, paste0(dataset, "_", sitecode, '_ppt.RData')))

    anom_12mth_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    anom_12mth_rast <- setValues(anom_12mth_rast,
                                 matrix(chirps_df$anom_12mth, 
                                        nrow=nrow(chirps)*ncol(chirps), 
                                        ncol=nlayers(chirps),
                                        byrow=TRUE))
    writeRaster(anom_12mth_rast,
                filename=file.path(out_folder, paste0(dataset, "_", sitecode, '_ppt_anom_12mth.tif')), 
                overwrite=overwrite)

    anom_24mth_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    anom_24mth_rast <- setValues(anom_24mth_rast,
                                 matrix(chirps_df$anom_24mth, 
                                        nrow=nrow(chirps)*ncol(chirps), 
                                        ncol=nlayers(chirps),
                                        byrow=TRUE))
    writeRaster(anom_24mth_rast,
                filename=file.path(out_folder, paste0(dataset, "_", sitecode, '_ppt_anom_24mth.tif')), 
                overwrite=overwrite)

    data.frame(sitecode=sitecode, date=date_limits_string, 
               ppt_annual_mean=cellStats(ppt_mean_12mth_rast, 'mean'),
               ppt_annual_min=cellStats(ppt_mean_12mth_rast, 'min'),
               ppt_annual_max=cellStats(ppt_mean_12mth_rast, 'max'),
               ppt_annual_sd=cellStats(ppt_mean_12mth_rast, 'sd'))
}

save(zoi_ppt, file=file.path(out_folder, paste0(dataset, '_ZOI_ppt.RData')))
