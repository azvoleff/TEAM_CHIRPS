###############################################################################
# Extracts SPI time series for TEAM vegetation plots
###############################################################################

library(raster)
library(stringr)
library(plyr)
library(reshape2)

library(foreach)
library(doParallel)
registerDoParallel(4)

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

spi_files <- dir(in_folder, 'SPI_[0-9]{1,2}.tif$')

load('vg_pts.RData')

spis <- foreach(spi_file=iter(spi_files),
                .packages=c('stringr', 'raster', 'reshape2'),
                .inorder=FALSE, .combine=rbind) %dopar% {
    this_spi <- brick(file.path(in_folder, spi_file))
    this_sitecode <- gsub('monthly_', '', str_extract(spi_file, 
                                                      'monthly_[A-Z]{2,3}'))
    spi_period <- gsub('SPI_', '', str_extract(spi_file, 'SPI_[0-9]{1,2}'))
    these_plots <- vg_pts[vg_pts$sitecode == this_sitecode, ]
    plot_spis <- extract(this_spi, these_plots, df=TRUE)
    names(plot_spis) <- gsub('monthly_', '', names(plot_spis))
    plot_spis <- plot_spis[names(plot_spis) != 'ID']
    plot_spis <- cbind(sitecode=these_plots$sitecode,
                       plot_ID=these_plots$Unit_ID, 
                       plot_num=these_plots$number,
                       spi_period=spi_period, plot_spis)
    plot_spis <- melt(plot_spis,
                      id.vars=c('sitecode', 'plot_ID', 'plot_num', 
                                'spi_period'), variable.name='date', 
                      value.name='spi')
    plot_spis$date <- as.numeric(gsub('[A-Z]{2,3}_SPI_[0-9]{1,2}.', '', 
                                      plot_spis$date))
    plot_spis$date <- dates[plot_spis$date]
    return(plot_spis)
}
save(spis, file='vg_plot_spis.RData')
