###############################################################################
# Extracts pentad precipitation timeseries for TEAM vegetation plots
###############################################################################

library(raster)
library(stringr)
library(reshape2)
library(lubridate)

library(doParallel)
library(foreach)
n_cpus <- 4

registerDoParallel(n_cpus)

overwrite <- TRUE

sites <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')
sitecodes <- sites$sitecode

product <- 'v1p8chirps'
chirps_NA_value <- -9999
dataset <- 'pentad'
date_limits_string <- '198101-201424'
#dataset <- 'monthly' # For SPI, use monthly
#date_limits_string <- '198101-201404'
# Note the below code is INCLUSIVE of the start date
chirps_start_date <- as.Date('1981/1/1')
# Note the below code is EXCLUSIVE of the end date
chirps_end_date <- as.Date('2014/5/1')

#in_base_dir <- 'O:/Data/CHIRPS'
#out_base_dir <- 'O:/Data/CHIRPS'
in_base_dir <- 'H:/Data/CHIRPS'
out_base_dir <- 'H:/Data/CHIRPS'
in_folder <- file.path(in_base_dir, paste0('global_', dataset))
out_folder <- file.path(out_base_dir, paste0('global_', dataset))
stopifnot(file_test('-d', in_folder))
stopifnot(file_test('-d', out_folder))

if (dataset == 'monthly') {
    dates <- seq(chirps_start_date, chirps_end_date, by='months')
    dates <- dates[dates < chirps_end_date]
    num_periods <- 12
} else if (dataset == 'pentad') {
    yrs <- seq(year(chirps_start_date), year(chirps_end_date))
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

load('vg_pts.RData')

vg_plot_ppt <- foreach(sitecode=iter(sitecodes), .inorder=FALSE,
                       .packages=c('stringr', 'raster', 'reshape2', 
                                   'lubridate'),
                       .combine=rbind) %dopar% {
    these_plots <- vg_pts[vg_pts$sitecode == sitecode, ]

    chirps_file <- file.path(out_folder,
                          paste0(product, '_', dataset, '_', sitecode, '_', 
                                 date_limits_string, '.dat'))
    stopifnot(file_test('-f', chirps_file))

    chirps <- brick(chirps_file)
    chirps <- calc(chirps, function(vals) {
        vals[vals == chirps_NA_value] <- NA
        return(vals)
    })

    plot_ppt <- extract(chirps, these_plots, df=TRUE)
    plot_ppt <- plot_ppt[names(plot_ppt) != 'ID']
    plot_ppt <- cbind(sitecode=these_plots$sitecode, 
                      plot_ID=these_plots$Unit_ID, plot_num=these_plots$number, 
                      plot_ppt)
    plot_ppt <- melt(plot_ppt,
                     id.vars=c('sitecode', 'plot_ID', 'plot_num'), 
                     variable.name='date', value.name='ppt')
    plot_ppt$date <- as.numeric(str_extract(plot_ppt$date, '[0-9]*$'))
    if (dataset == 'pentad') {
        plot_ppt$pentad <- plot_ppt$date %% num_periods
    }
    plot_ppt$date <- dates[plot_ppt$date]
    plot_ppt$year <- year(plot_ppt$date)
    plot_ppt$month <- month(plot_ppt$date)

    return(plot_ppt)
}

save(vg_plot_ppt, file=paste0(dataset, '_vg_plot_ppt.RData'))
