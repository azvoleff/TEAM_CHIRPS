###############################################################################
# Plots basic stats from CHIRPS pentad precipitation data w/in TEAM ZOIs
###############################################################################

library(ggplot2)
library(lubridate)
library(dplyr)
library(scales) # for date_breaks

in_base_dir <- 'H:/Data/CHIRPS'
in_folder <- file.path(in_base_dir, 'global_pentad')
stopifnot(file_test('-d', in_folder))

width <- 10
height <- 7.5
dpi <- 300
transparent_opts <- theme(legend.position="bottom", 
                          axis.text=element_text(colour='white'), 
                          axis.title.x=element_text(colour='white'), 
                          axis.title.y=element_text(colour='white'), 
                          plot.background=element_rect(fill='transparent', colour=NA))

load(file=file.path(in_folder, 'pentad_ZOI_ppt.RData'))

zoi_ppt <- tbl_df(zoi_ppt)
zoi_ppt$year <- year(zoi_ppt$date)
zoi_ppt$month <- month(zoi_ppt$date)

monthly_totals <- summarize(group_by(zoi_ppt, year, month, sitecode),
                            total=sum(ppt_mean))
monthly_mean_total <- summarize(group_by(monthly_totals, month, sitecode),
                           total=mean(total))
ggplot(monthly_mean_total, aes(month, total)) +
    geom_bar(stat='identity') + facet_wrap(~sitecode) +
    xlab('Month') +
    ylab('Total precipitation (mm)') +
    scale_x_continuous(breaks=seq(1, 12), labels=c('J', 'F', 'M', 'A', 'M', 
                                                   'J', 'J', 'A', 'S', 'O', 
                                                   'N', 'D')) +
    theme_grey(base_size=18) + transparent_opts
ggsave('zoi_ppt_monthly_total.png', width=width, height=height, dpi=dpi,
       bg='transparent')
#ggsave('zoi_ppt_monthly_total_not_transparent.png', width=width, height=height, dpi=dpi)

bif <- ggplot(filter(monthly_mean_total, sitecode == "BIF"), aes(month, total)) +
    geom_bar(stat='identity') + 
    xlab('Month') +
    ylab('Total precipitation (mm)') +
    scale_x_continuous(breaks=seq(1, 12), labels=c('J', 'F', 'M', 'A', 'M', 
                                                   'J', 'J', 'A', 'S', 'O', 
                                                   'N', 'D')) +
    theme_bw(base_size=18)
ggsave('Bwindi_ppt_monthly_total.png', bif, width=width, height=height, dpi=dpi)
ggsave('Bwindi_ppt_monthly_total.eps', bif, width=width, height=height, dpi=dpi)

bbs <- ggplot(filter(monthly_mean_total, sitecode == "BBS"), aes(month, total)) +
    geom_bar(stat='identity') + 
    xlab('Month') +
    ylab('Total precipitation (mm)') +
    scale_x_continuous(breaks=seq(1, 12), labels=c('J', 'F', 'M', 'A', 'M', 
                                                   'J', 'J', 'A', 'S', 'O', 
                                                   'N', 'D')) +
    theme_bw(base_size=18)
ggsave('BukitBarisan_ppt_monthly_total.png', bbs, width=width, height=height, dpi=dpi)
ggsave('BukitBarisan_ppt_monthly_total.eps', bbs, width=width, height=height, dpi=dpi)

annual_total <- summarize(group_by(zoi_ppt, year, sitecode),
                          total=sum(ppt_mean))
ggplot(annual_total[annual_total$year != 2014, ], aes(year, total)) +
    geom_line() + facet_wrap(~sitecode) +
    xlab('Year') +
    ylab('Total precipitation (mm)') +
    theme_grey(base_size=18) + transparent_opts
ggsave('zoi_ppt_annual_total.png', width=width, height=height, dpi=dpi,
       bg='transparent')
#ggsave('zoi_ppt_annual_total_not_transparent.png', width=width, height=height, dpi=dpi)

annual_total <- summarize(group_by(zoi_ppt, year, sitecode),
                          mean_anom=mean(ppt_anom_mean))
ggplot(annual_total[annual_total$year != 2014, ], aes(year, mean_anom)) +
    geom_line() + facet_wrap(~sitecode) +
    xlab('Year') +
    ylab('Mean anomaly (mm)') +
    theme_grey(base_size=18) + transparent_opts
ggsave('zoi_ppt_anom_mean.png', width=width, height=height, dpi=dpi,
       bg='transparent')
#ggsave('zoi_ppt_anom_mean_not_transparent.png', width=width, height=height, dpi=dpi)

ggplot(zoi_ppt, aes(date, ppt_anom_mean)) +
    geom_line() + facet_wrap(~sitecode) +
    xlab('Year') +
    ylab('Mean precipitation anomaly (z-score)') +
    theme_grey(base_size=18) + transparent_opts
ggsave('zoi_ppt_anom_mean.png', width=width, height=height, dpi=dpi,
       bg='transparent')


