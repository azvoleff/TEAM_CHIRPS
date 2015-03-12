source("0_settings.R")

library(ggplot2)

load(file.path(out_folder, paste0('ALL_ZOIs_CHIRPS_', dataset, '_', start_date, 
                                  '-', end_date, '_ZOI_cwd_stats.RData')))

ggplot(zoi_cwd, aes(date, cwd_mean)) + geom_line() + facet_wrap(~sitecode)
ggsave("cwd_mean_by_site.png")

ggplot(zoi_cwd, aes(date, mcwd12_mean)) + geom_line() + facet_wrap(~sitecode)
ggsave("mcwd12_mean_by_site.png")

ggplot(zoi_cwd, aes(date, cwd_run12_mean)) + geom_line() + facet_wrap(~sitecode)
ggsave("cwd_run12_mean_by_site.png")

ggplot(zoi_cwd, aes(date, mcwd_run12_mean)) + geom_line() + facet_wrap(~sitecode)
ggsave("mcwd_run12_mean_by_site.png")

ggplot(filter(zoi_cwd, date > as.Date("1999/12/31")), aes(date, mcwd_run12_mean)) + geom_line() + facet_wrap(~sitecode)

ggplot(filter(zoi_cwd, date > as.Date("1999/12/31")),
       aes(date, mcwd_run12_mean, colour=reorder(sitecode, mcwd_run12_mean, order=TRUE))) +
    geom_line() + scale_colour_discrete("Site")

load(file.path(out_folder, paste0('ALL_ZOIs_CHIRPS_', dataset, '_', start_date, 
                                  '-', end_date, '_ZOI_ppt_stats.RData')))
zoi_ppt$sitecode <- factor(zoi_ppt$sitecode, 
                           levels=zoi_ppt$sitecode[order(zoi_ppt$ppt_annual_mean)])
ggplot(zoi_ppt, aes(sitecode, ppt_annual_mean)) +
    geom_bar(stat="identity")
ggsave("ppt_annual_by_site.png")
