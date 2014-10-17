###############################################################################
# Extracts cwd time series for TEAM vegetation plots
###############################################################################

library(raster)
library(stringr)
library(dplyr)
library(reshape2)

library(foreach)

sites <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')
sitecodes <- sites$sitecode

in_base_dir <- 'H:/Data/CHIRPS'
in_folder <- file.path(in_base_dir, 'TEAM_monthly')
stopifnot(file_test('-d', in_folder))

# Note the below code is INCLUSIVE of the start date
chirps_start_date <- as.Date('1981/1/1')
# Note the below code is EXCLUSIVE of the end date
chirps_end_date <- as.Date('2014/5/1')
dates <- seq(chirps_start_date, chirps_end_date, by='months')
dates <- dates[dates < chirps_end_date]
num_periods <- 12

cwd_files <- dir(in_folder, 'cwd.tif$')

load('vg_pts.RData')

get_cwd_data <- function (cwd_type, plots) {
    this_data <- brick(file.path(in_folder, paste0("monthly_", sitecode, "_", 
                                                   cwd_type, ".tif")))
    plot_datas <- extract(this_data, plots, df=TRUE)
    plot_datas <- plot_datas[names(plot_datas) != 'ID']
    plot_datas  <- data.frame(t(plot_datas))
    names(plot_datas) <- plots$Unit_ID
    plot_datas <- cbind(date=dates, plot_datas)
    plot_datas <- melt(plot_datas, id.vars='date', variable.name='plot_ID', 
                       value.name=cwd_type)
    plot_datas <- cbind(sitecode=plots$sitecode, plot_datas)
}

cwds <- foreach(sitecode=sitecodes,
                .packages=c('stringr', 'raster', 'reshape2'),
                .inorder=FALSE, .combine=rbind) %do% {
    these_plots <- vg_pts[vg_pts$sitecode == sitecode, ]

    plot_cwds <- get_cwd_data("cwd", these_plots)
    plot_mcwd12s <- get_cwd_data("mcwd12", these_plots)
    plot_cwd_run12s <- get_cwd_data("cwd_run12", these_plots)
    plot_mcwd_run12s <- get_cwd_data("mcwd_run12", these_plots)

    # Pull out max 12-month running cumulative water deficit
    stopifnot(nrow(plot_cwds) == nrow(plot_mcwd12s))
    stopifnot(nrow(plot_mcwd12s) == nrow(plot_cwd_run12s))
    stopifnot(nrow(plot_cwd_run12s) == nrow(plot_mcwd_run12s))

    plot_cwd_indicators <- cbind(plot_cwds,
                                 mcwd12=plot_mcwd12s$mcwd12, 
                                 cwd_run12=plot_cwd_run12s$cwd_run12, 
                                 mcwd_run12=plot_mcwd_run12s$mcwd_run12)
    return(plot_cwd_indicators)
}
cwds <- arrange(cwds, sitecode, plot_ID, date)
save(cwds, file='vg_plot_cwds.RData')
