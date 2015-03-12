###############################################################################
# Extracts SPI time series for TEAM vegetation plots
###############################################################################

source('0_settings.R')

stopifnot(dataset == 'monthly')

library(raster)
library(stringr)
library(reshape2)

library(foreach)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

spi_files <- dir(out_folder, 'SPI_[0-9]{1,2}.tif$')

load('vg_pts.RData')

# Only can extract timeseries for sites that have veg plots sampled, so need to 
# exclude spi files from sites without vegetation plots:
sitecodes <- as.character(sitecodes[sitecodes %in% vg_pts$sitecode])
spi_files <- spi_files[str_extract(spi_files, '^[A-Z]{2,3}') %in% sitecodes]

spis <- foreach(spi_file=spi_files,
                .packages=c('stringr', 'raster', 'reshape2'),
                .inorder=FALSE, .combine=rbind) %dopar% {
    this_spi <- brick(file.path(out_folder, spi_file))
    this_sitecode <- str_extract(spi_file, '^[A-Z]{2,3}')
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
    plot_spis$date <- as.numeric(str_extract(plot_spis$date, '[0-9]*$'))
    plot_spis$date <- dates[plot_spis$date]
    return(plot_spis)
}

save(spis, file='vg_plot_spis.RData')
