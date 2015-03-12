###############################################################################
# Plots basic stats from CHIRPS pentad precipitation data w/in TEAM ZOIs
###############################################################################

source('0_settings.R')

library(ggplot2)
library(lubridate)
library(dplyr)
library(scales) # for date_breaks

#load(file.path(in_folder, 'pentad_ZOI_ppt.RData'))
load(file.path(out_folder, 'ALL_ZOIs_CHIRPS_monthly_198101-201412_ZOI_ppt_stats.RData'))

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


for (this_sitecode in sitecodes) {
    p <- ggplot(filter(monthly_mean_total, sitecode == this_sitecode), aes(month, total)) +
        geom_bar(stat='identity') + 
        xlab('Month') +
        ylab('Total precipitation (mm)') +
        scale_x_continuous(breaks=seq(1, 12), labels=c('J', 'F', 'M', 'A', 'M', 
                                                       'J', 'J', 'A', 'S', 'O', 
                                                       'N', 'D')) +
        theme_bw(base_size=18)
    ggsave(paste0('ppt_monthly_total_', this_sitecode, '.png'), p, width=width, height=height, dpi=dpi)
    #ggsave(paste0('ppt_monthly_total_', sitecode, '.eps'), p, width=width, 
    #height=height, dpi=dpi)
}

annual_total <- summarize(group_by(zoi_ppt, year, sitecode),
                          total=sum(ppt_mean))
ggplot(annual_total, aes(year, total)) +
    geom_line() + facet_wrap(~sitecode) +
    xlab('Year') +
    ylab('Total precipitation (mm)') +
    theme_grey(base_size=18) + transparent_opts
ggsave('zoi_ppt_annual_total.png', width=width, height=height, dpi=dpi,
       bg='transparent')
#ggsave('zoi_ppt_annual_total_not_transparent.png', width=width, height=height, dpi=dpi)

ggplot(zoi_ppt, aes(date, ppt_anom_mean)) +
    geom_line() + facet_wrap(~sitecode) +
    xlab('Year') +
    ylab('Mean anomaly (mm)') +
    theme_grey(base_size=18) + transparent_opts
ggsave('zoi_ppt_anom_mean.png', width=width, height=height, dpi=dpi,
       bg='transparent')
#ggsave('zoi_ppt_anom_mean_not_transparent.png', width=width, height=height, dpi=dpi)
