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
library(Rcpp)
library(inline)
library(ggplot2)

library(doParallel)
library(foreach)

n_cpus <- 4
cl <- makeCluster(n_cpus)
registerDoParallel(cl)

evap <- read.csv("evapotranspirations.csv")

zoi_cwd <- foreach(sitecode=sitecodes, .inorder=FALSE,
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
    save(chirps_df, file=paste0(base_name, '_ppt.RData'))

    # Use chirps raster as a template
    cwd_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    cwd_rast <- setValues(cwd_rast, matrix(chirps_df$cwd, 
                                           nrow=nrow(chirps)*ncol(chirps), 
                                           ncol=nlayers(chirps),
                                           byrow=TRUE))
    writeRaster(cwd_rast,
                filename=paste0(base_name, '_cwd.tif'), 
                overwrite=overwrite)

    cwd_run12_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    cwd_run12_rast <- setValues(cwd_run12_rast, matrix(chirps_df$cwd_run12, 
                                                       nrow=nrow(chirps)*ncol(chirps), 
                                                       ncol=nlayers(chirps), 
                                                       byrow=TRUE))
    writeRaster(cwd_run12_rast,
                filename=paste0(base_name, '_cwd_run12.tif'), 
                overwrite=overwrite)

    mcwd12_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    mcwd12_rast <- setValues(mcwd12_rast, matrix(chirps_df$mcwd12, 
                                                 nrow=nrow(chirps)*ncol(chirps), 
                                                 ncol=nlayers(chirps),
                                                 byrow=TRUE))
    writeRaster(mcwd12_rast,
                filename=paste0(base_name, '_mcwd12.tif'), 
                overwrite=overwrite)

    mcwd_run12_rast <- brick(chirps, values=FALSE, nl=nlayers(chirps))
    mcwd_run12_rast <- setValues(mcwd_run12_rast, matrix(chirps_df$mcwd_run12, 
                                                         nrow=nrow(chirps)*ncol(chirps), 
                                                         ncol=nlayers(chirps),
                                                         byrow=TRUE))
    writeRaster(mcwd_run12_rast,
                filename=paste0(base_name, '_mcwd_run12.tif'), 
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
    writeRaster(mcwd_rast,
                filename=paste0(base_name, '_cwd_annual_max.tif'), 
                overwrite=overwrite)

    data.frame(sitecode=sitecode, date=dates,
               cwd_mean=cellStats(mask(cwd_rast, zoi), 'mean'),
               cwd_min=cellStats(mask(cwd_rast, zoi), 'min'),
               cwd_max=cellStats(mask(cwd_rast, zoi), 'max'),
               cwd_sd=cellStats(mask(cwd_rast, zoi), 'sd'),
               mcwd12_mean=cellStats(mask(mcwd12_rast, zoi), 'mean'),
               mcwd12_min=cellStats(mask(mcwd12_rast, zoi), 'min'),
               mcwd12_max=cellStats(mask(mcwd12_rast, zoi), 'max'),
               mcwd12_sd=cellStats(mask(mcwd12_rast, zoi), 'sd'),
               cwd_run12_mean=cellStats(mask(cwd_run12_rast, zoi), 'mean'),
               cwd_run12_min=cellStats(mask(cwd_run12_rast, zoi), 'min'),
               cwd_run12_max=cellStats(mask(cwd_run12_rast, zoi), 'max'),
               cwd_run12_sd=cellStats(mask(cwd_run12_rast, zoi), 'sd'),
               mcwd_run12_mean=cellStats(mask(mcwd_run12_rast, zoi), 'mean'),
               mcwd_run12_min=cellStats(mask(mcwd_run12_rast, zoi), 'min'),
               mcwd_run12_max=cellStats(mask(mcwd_run12_rast, zoi), 'max'),
               mcwd_run12_sd=cellStats(mask(mcwd_run12_rast, zoi), 'sd'))
}
save(zoi_cwd, file=file.path(out_folder,
                             paste0('ALL_ZOIs_CHIRPS_', dataset, '_', 
                                    start_date, '-', end_date, 
                                    '_ZOI_cwd_stats.RData')))
stopCluster(cl)
