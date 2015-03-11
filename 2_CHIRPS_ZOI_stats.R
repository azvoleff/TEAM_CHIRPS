###############################################################################
# Calculates mean precipitation within ZOI of team site from CHIRPS pentad or 
# monthly data.
###############################################################################

source("0_settings.R")

library(rgdal)
library(raster)
library(rgeos)
library(lubridate)
library(dplyr)

library(doParallel)
library(foreach)

cl <- makeCluster(4)

registerDoParallel(cl)

zoi_ppt <- foreach(sitecode=iter(sitecodes), .inorder=FALSE,
                   .packages=c('dplyr', 'raster', 'rgdal', 'sp', 'rgeos'),
                   .combine=rbind) %do% {
    print(sitecode)

    zoi_file <- dir(zoi_folder, pattern=paste0('^ZOI_', sitecode, 
                                               '_[0-9]{4}.RData'), 
                    full.names=TRUE)
    stopifnot(length(zoi_file) == 1)
    load(zoi_file)
    zoi <- spTransform(zoi, CRS(s_srs))

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

    # Setup a dataframe with the precipitation data so anomalies, etc can be 
    # calculated
    years_rep <- rep(years, each=nrow(chirps)*ncol(chirps))
    subyears_rep <- rep(subyears, each=nrow(chirps)*ncol(chirps))
    pixel_IDs_rep <- rep(seq(1:(nrow(chirps)*ncol(chirps))), nlayers(chirps))
    chirps_df <- data.frame(year=years_rep,
                            subyear=subyears_rep, 
                            pixel=pixel_IDs_rep,
                            ppt=as.vector(chirps))
    chirps_df <- tbl_df(chirps_df)

    # Add 12 month total precip column. Use sides=1 in filter to have it be 12 
    # month total as of that date
    chirps_df <- group_by(chirps_df, pixel) %>%
        arrange(year, subyear) %>%
        mutate(tot_12mth=as.numeric(stats::filter(ppt, rep(1/12, 12), sides=1)*12))

    # First calculate the mean_monthly for each subyear for each pixel
    ppt_mean_monthly <- group_by(chirps_df, pixel, subyear) %>%
        summarize(mean_monthly=mean(ppt, na.rm=TRUE))
    # Use chirps raster as a template
    ppt_mean_monthly_rast <- brick(chirps, values=FALSE, nl=num_periods)
    ppt_mean_monthly_rast <- setValues(ppt_mean_monthly_rast,
                               matrix(ppt_mean_monthly$mean_monthly, 
                                      nrow=nrow(chirps)*ncol(chirps), 
                                      ncol=num_periods, byrow=TRUE))
    writeRaster(ppt_mean_monthly_rast,
                filename=paste0(base_name, '_ppt_mean_monthly.tif'), 
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
    writeRaster(ppt_mean_12mth_rast,
                filename=paste0(base_name,'_ppt_mean_12mth.tif'), 
                overwrite=overwrite)

    # Calculate 12 month precipitation anomalies
    chirps_df$anom_12mth <- chirps_df$tot_12mth - ppt_mean_12mth$mean_12mth[match(chirps_df$pixel, ppt_mean_12mth$pixel)]
    #filter(chirps_df, pixel == 1)[1:36,]
    save(chirps_df, file=paste0(base_name, '_ppt.RData'))

    anom_12mth_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    anom_12mth_rast <- setValues(anom_12mth_rast,
                                 matrix(chirps_df$anom_12mth, 
                                        nrow=nrow(chirps)*ncol(chirps), 
                                        ncol=nlayers(chirps),
                                        byrow=TRUE))
    writeRaster(anom_12mth_rast,
                filename=paste0(base_name, '_ppt_anom_12mth.tif'), 
                overwrite=overwrite)

    data.frame(sitecode=sitecode, date=paste(start_date, end_date, sep="-"), 
               ppt_annual_mean=cellStats(mask(ppt_mean_12mth_rast, zoi), 'mean'),
               ppt_annual_min=cellStats(mask(ppt_mean_12mth_rast, zoi), 'min'),
               ppt_annual_max=cellStats(mask(ppt_mean_12mth_rast, zoi), 'max'),
               ppt_annual_sd=cellStats(mask(ppt_mean_12mth_rast, zoi), 'sd'))
}
save(zoi_ppt, file=file.path(out_folder,
                             paste0('ALL_ZOIs_CHIRPS_', dataset, '_', 
                                    start_date, '-', end_date, 
                                    '_ZOI_ppt_stats.RData')))

stopCluster(cl)
