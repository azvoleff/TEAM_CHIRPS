###############################################################################
# Calculates mean precipitation within ZOI of team site from CHIRPS pentad or 
# monthly data.
###############################################################################

library(rgdal)
library(raster)
library(rgeos)
library(lubridate)
library(dplyr)
library(Rcpp)
library(inline)

library(doParallel)
library(foreach)
n_cpus <- 4

cl <- makeCluster(n_cpus)
registerDoParallel(cl)

overwrite <- TRUE

sites <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')
sitecodes <- sites$sitecode

product <- 'v1p8chirps'
chirps_NA_value <- -9999
#dataset <- 'pentad'
#date_limits_string <- '198101-201424'
dataset <- 'monthly' # For SPI, use monthly
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

evap <- read.csv("evapotranspirations.csv")

zoi_cwd <- foreach(sitecode=iter(sitecodes), .inorder=FALSE,
                   .packages=c('dplyr', 'raster', 'rgdal', 'sp', 'Rcpp', 
                               'inline'),
                   .combine=rbind) %dopar% {

    # Define a function to calculate cumulative water deficit, given a 
    # precipitation time series and evapotranspiration rate
    src <- '
        Rcpp::NumericVector x = Rcpp::NumericVector(X);
        unsigned sz = x.size();
        Rcpp::NumericVector deficit(sz);
        double evap = as<double>(EVAP);
        if (evap < 0) throw std::range_error("Evapotranspiration rate must be positive");
        deficit[0] = 0;
        for (unsigned n = 1; n < sz; n++) {
            if ((deficit[n - 1] - evap + x[n]) < 0) {
                deficit[n] = deficit[n - 1] - evap + x[n];
            } else {
                deficit[n] = 0;
            }
        }
        deficit[0] = NA_REAL;
        return(deficit);
        '
    cwd <- cxxfunction(signature(X="numeric", EVAP="numeric"), body=src, 
                        plugin="Rcpp")

    # Define a function to calculate cumulative water deficit, given a 
    # precipitation time series and evapotranspiration rate
    #
    # This function calculate a _running_ CWD for the specified period, as 
    # opposed to the above which is cumulative across the whole timeseries
    src <- '
        Rcpp::NumericVector x = Rcpp::NumericVector(X);
        unsigned sz = x.size();
        Rcpp::NumericVector deficit(sz);
        double evap = as<double>(EVAP);
        unsigned period = as<unsigned>(PERIOD);
        deficit[0] = 0;
        if (period > sz) throw std::range_error("Period cannot be greater than length of series");
        if (evap < 0) throw std::range_error("Evapotranspiration rate must be positive");
        for (unsigned n = period - 1; n < sz; n++) {
            Rcpp::NumericVector period_deficit(period);
            period_deficit[0] = 0;
            for (unsigned i = 1; i < period; i++) {
                if ((period_deficit[i - 1] - evap + x[n - period + i]) < 0) {
                    period_deficit[i] = period_deficit[i - 1] - evap + x[n - period + i];
                } else {
                    period_deficit[i] = 0;
                }
            }
            deficit[n - 1] = period_deficit[period - 1];
            //Rf_PrintValue(period_deficit);
        }
        for (unsigned i = 0; i < period - 1; i++) {
            deficit[i] = NA_REAL;
        }
        return(deficit);
        '
    cwd_running <- cxxfunction(signature(X="numeric", EVAP="numeric", 
                                         PERIOD="integer"), body=src, 
                               plugin="Rcpp")

    # Define a function to calculate maximum cumulative water deficit over a 
    # given cwd time series and period
    src <- '
        Rcpp::NumericVector cwd = Rcpp::NumericVector(CWD);
        unsigned sz = cwd.size();
        Rcpp::NumericVector mcwd(sz);
        unsigned period = as<unsigned>(PERIOD);
        if (period > sz) throw std::range_error("Period cannot be greater than length of series");
        for (unsigned n = period - 1; n < sz; n++) {
            double this_mcwd = 0;
            for (unsigned i=0; i < period; i++) {
                if (cwd[n - i] < this_mcwd) this_mcwd = cwd[n - i];
            }
            mcwd[n] = this_mcwd;
        }
        for (unsigned i = 0; i < period - 1; i++) {
            mcwd[i] = NA_REAL;
        }
        return(mcwd);
        '
    mcwd <- cxxfunction(signature(CWD="numeric", PERIOD="integer"), body=src, 
                        plugin="Rcpp")

    zoi_file <- dir(zoi_folder,
                    pattern=paste0('^ZOI_', sitecode, '_[0-9]{4}.RData'), 
                    full.names=TRUE)
    stopifnot(length(zoi_file) == 1)
    load(zoi_file)
    zoi <- spTransform(zoi, CRS(s_srs))

    evapotranspiration <- evap$monthly_rate[evap$sitecode == sitecode]

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

    chirps_df <- group_by(chirps_df, pixel) %>%
        arrange(year, subyear) %>%
        mutate(cwd=cwd(ppt, evapotranspiration),
               cwd_run6=cwd_running(ppt, evapotranspiration, 6),
               cwd_run12=cwd_running(ppt, evapotranspiration, 12))
    chirps_df <- group_by(chirps_df, pixel) %>%
        arrange(year, subyear) %>%
        mutate(mcwd6=mcwd(cwd, 6),
               mcwd12=mcwd(cwd, 12),
               mcwd_run6=mcwd(cwd_run6, 6),
               mcwd_run12=mcwd(cwd_run12, 12))
    chirps_df$cwd[!chirps_df$inzoi] <- NA
    chirps_df$cwd_run6[!chirps_df$inzoi] <- NA
    chirps_df$cwd_run12[!chirps_df$inzoi] <- NA
    chirps_df$mcwd_6[!chirps_df$inzoi] <- NA
    chirps_df$mcwd12[!chirps_df$inzoi] <- NA
    chirps_df$mcwd_run6[!chirps_df$inzoi] <- NA
    chirps_df$mcwd_run12[!chirps_df$inzoi] <- NA
    save(chirps_df, file=file.path(out_folder, paste0(dataset, "_", sitecode, '_ppt.RData')))

    # Use chirps raster as a template
    cwd_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    cwd_rast <- setValues(cwd_rast, matrix(chirps_df$cwd, 
                                           nrow=nrow(chirps)*ncol(chirps), 
                                           ncol=nlayers(chirps),
                                           byrow=TRUE))
    writeRaster(cwd_rast,
                filename=file.path(out_folder, paste0(dataset, "_", sitecode, '_cwd.tif')), 
                overwrite=overwrite)

    cwd_run12_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    cwd_run12_rast <- setValues(cwd_run12_rast, matrix(chirps_df$cwd_run12, 
                                                       nrow=nrow(chirps)*ncol(chirps), 
                                                       ncol=nlayers(chirps), 
                                                       byrow=TRUE))
    writeRaster(cwd_run12_rast,
                filename=file.path(out_folder, paste0(dataset, "_", sitecode, '_cwd_run12.tif')), 
                overwrite=overwrite)

    mcwd12_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    mcwd12_rast <- setValues(mcwd12_rast, matrix(chirps_df$mcwd12, 
                                                 nrow=nrow(chirps)*ncol(chirps), 
                                                 ncol=nlayers(chirps),
                                                 byrow=TRUE))
    writeRaster(mcwd12_rast,
                filename=file.path(out_folder, paste0(dataset, "_", sitecode, '_mcwd12.tif')), 
                overwrite=overwrite)

    mcwd_run12_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    mcwd_run12_rast <- setValues(mcwd_run12_rast, matrix(chirps_df$mcwd_run12, 
                                                         nrow=nrow(chirps)*ncol(chirps), 
                                                         ncol=nlayers(chirps),
                                                         byrow=TRUE))
    writeRaster(mcwd_run12_rast,
                filename=file.path(out_folder, paste0(dataset, "_", sitecode, '_mcwd_run12.tif')), 
                overwrite=overwrite)

    # Calculate maximum cumulative water deficit for each year
    mcwd <- group_by(chirps_df, pixel, year) %>%
        summarize(mcwd=min(cwd))
    # Use chirps raster as a template
    mcwd_rast <- brick(chirps, values=FALSE, nl=num_periods)
    mcwd_rast <- setValues(mcwd_rast, matrix(mcwd$mcwd, 
                                             nrow=nrow(chirps)*ncol(chirps), 
                                             ncol=length(unique(mcwd$year)), 
                                             byrow=TRUE))
    # Mask areas outside ZOI
    mcwd_rast <- mcwd_rast * chirps_mask
    writeRaster(mcwd_rast,
                filename=file.path(out_folder, paste0(dataset, "_", sitecode, '_cwd_annual_max.tif')), 
                overwrite=overwrite)

    data.frame(sitecode=sitecode, date=dates,
               cwd_mean=cellStats(cwd_rast, 'mean'),
               cwd_min=cellStats(cwd_rast, 'min'),
               cwd_max=cellStats(cwd_rast, 'max'),
               cwd_sd=cellStats(cwd_rast, 'sd'),
               mcwd12_mean=cellStats(mcwd12_rast, 'mean'),
               mcwd12_min=cellStats(mcwd12_rast, 'min'),
               mcwd12_max=cellStats(mcwd12_rast, 'max'),
               mcwd12_sd=cellStats(mcwd12_rast, 'sd'),
               cwd_run12_mean=cellStats(cwd_run12_rast, 'mean'),
               cwd_run12_min=cellStats(cwd_run12_rast, 'min'),
               cwd_run12_max=cellStats(cwd_run12_rast, 'max'),
               cwd_run12_sd=cellStats(cwd_run12_rast, 'sd'),
               mcwd_run12_mean=cellStats(mcwd_run12_rast, 'mean'),
               mcwd_run12_min=cellStats(mcwd_run12_rast, 'min'),
               mcwd_run12_max=cellStats(mcwd_run12_rast, 'max'),
               mcwd_run12_sd=cellStats(mcwd_run12_rast, 'sd'))
}
save(zoi_cwd, file=file.path(out_folder, paste0(dataset, '_ZOI_cwd_stats.RData')))

library(ggplot2)
ggplot(zoi_cwd, aes(date, cwd_mean)) + geom_line() + facet_wrap(~sitecode)
ggsave("cwd_mean_by_site.png")

ggplot(zoi_cwd, aes(date, mcwd12_mean)) + geom_line() + facet_wrap(~sitecode)
ggsave("mcwd12_mean_by_site.png")

ggplot(zoi_cwd, aes(date, cwd_run12_mean)) + geom_line() + facet_wrap(~sitecode)
ggsave("cwd_run12_mean_by_site.png")

ggplot(zoi_cwd, aes(date, mcwd_run12_mean)) + geom_line() + facet_wrap(~sitecode)
ggsave("mcwd_run12_mean_by_site.png")

ggplot(filter(zoi_cwd, date > as.Date("1999/12/31")), aes(date, mcwd_run12_mean)) + geom_line() + facet_wrap(~sitecode)

ggplot(filter(zoi_cwd, date > as.Date("1999/12/31")),
       aes(date, mcwd_run12_mean, colour=reorder(sitecode, mcwd_run12_mean, order=TRUE))) +
    geom_line() + scale_colour_discrete("Site")

load(file.path(out_folder, "monthly_ZOI_ppt.RData"))
zoi_ppt$sitecode <- factor(zoi_ppt$sitecode, 
                           levels=zoi_ppt$sitecode[order(zoi_ppt$ppt_annual_mean)])
ggplot(zoi_ppt, aes(sitecode, ppt_annual_mean)) +
    geom_bar(stat="identity")
ggsave("ppt_annual_by_site.png")
