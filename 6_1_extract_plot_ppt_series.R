###############################################################################
# Extracts precipitation timeseries for TEAM vegetation plots
###############################################################################

source('0_settings.R')

library(stringr)
library(raster)
library(reshape2)
library(lubridate)

library(doParallel)
library(foreach)
n_cpus <- 4
cl <- makeCluster(4)
registerDoParallel(cl)

load('vg_pts.RData')

# Only can extract timeseries for sites that have veg plots sampled:
sitecodes <- as.character(sitecodes[sitecodes %in% vg_pts$sitecode])

vg_plot_ppt <- foreach(sitecode=sitecodes, .inorder=FALSE,
                       .packages=c('stringr', 'raster', 'reshape2', 
                                   'lubridate'),
                       .combine=rbind) %dopar% {
    these_plots <- vg_pts[vg_pts$sitecode == sitecode, ]

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
