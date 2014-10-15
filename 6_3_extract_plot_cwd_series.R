###############################################################################
# Extracts cwd time series for TEAM vegetation plots
###############################################################################

library(raster)
library(stringr)
library(plyr)
library(reshape2)

library(foreach)

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

cwds <- foreach(cwd_file=cwd_files,
                .packages=c('stringr', 'raster', 'reshape2'),
                .inorder=FALSE, .combine=rbind) %do% {
    this_cwd <- brick(file.path(in_folder, cwd_file))
    this_sitecode <- gsub('monthly_', '', str_extract(cwd_file, 
                                                      'monthly_[A-Z]{2,3}'))
    these_plots <- vg_pts[vg_pts$sitecode == this_sitecode, ]
    plot_cwds <- extract(this_cwd, these_plots, df=TRUE)
    names(plot_cwds) <- gsub('monthly_', '', names(plot_cwds))
    plot_cwds <- plot_cwds[names(plot_cwds) != 'ID']
    plot_cwds <- cbind(sitecode=these_plots$sitecode,
                       plot_ID=these_plots$Unit_ID, 
                       plot_num=these_plots$number, plot_cwds)
    plot_cwds <- melt(plot_cwds,
                      id.vars=c('sitecode', 'plot_ID', 'plot_num'), 
                      variable.name='date', value.name='cwd')
    plot_cwds$date <- as.numeric(gsub('[A-Z]{2,3}_cwd.', '', plot_cwds$date))
    plot_cwds$date <- dates[plot_cwds$date]
    return(plot_cwds)
}
save(cwds, file='vg_plot_cwds.RData')
