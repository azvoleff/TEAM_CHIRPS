###############################################################################
# Extracts cwd time series for TEAM vegetation plots
###############################################################################

source('0_settings.R')

library(raster)
library(stringr)
library(dplyr)
library(reshape2)

library(foreach)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

load('vg_pts.RData')

# Only can extract timeseries for sites that have veg plots sampled, so need to 
# exclude cwd files from sites without vegetation plots:
sitecodes <- as.character(sitecodes[sitecodes %in% vg_pts$sitecode])

get_cwd_data <- function (cwd_type, plots) {
    cwd_file <- file.path(out_folder,
                           paste0(sitecode, '_CHIRPS_', dataset,
                                  '_', start_date, '-', end_date, "_", 
                                  cwd_type, ".tif"))
    this_data <- brick(cwd_file)
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
                .inorder=FALSE, .combine=rbind) %dopar% {
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

stopCluster(cl)
