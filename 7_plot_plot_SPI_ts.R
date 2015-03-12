source('0_settings.R')

library(ggplot2)
library(dplyr)
library(scales) # for date_breaks

load('vg_plot_spis.RData')

# Add in continent names and countries
spis <- merge(spis, sites, by.x='sitecode', by.y='sitecode')

# Remember spi is scaled by 1000:
spis$spi <- spis$spi / 1000

spis <- spis[spis$date > as.Date('2000/01/01'), ]

spis$sitename_pretty <- ordered(spis$sitename_pretty, levels=sites$sitename_pretty)

spi_means <- summarize(group_by(spis, spi_period, sitecode, date),
                       mean_spi=mean(spi),
                       sitename_pretty=sitename_pretty[1],
                       Continent=continent[1])
spi_means$spi_sign <- ifelse(spi_means$mean_spi >= 0, "wet", "dry")
spi_means$spi_period <- ordered(spi_means$spi_period, level=c(1, 6, 12, 24, 36))

spi_plot_periods <- unique(spis$spi_period)
for (spi_plot_period in spi_plot_periods) {
    scale_breaks <- seq(as.Date('2000/1/1'), as.Date('2015/1/1'), by='5 years')
    ggplot(spi_means[spi_means$spi_period == spi_plot_period, ]) +
        geom_line(aes(date, mean_spi, colour=Continent)) +
        xlab('Year') + scale_x_date(labels=date_format("'%y "), breaks=scale_breaks) +
        ylab('Standardized Precipitation Index (SPI)') +
        facet_wrap(~sitename_pretty) +
        geom_hline(yintercept=seq(-2, 2), alpha=.3) +
        geom_hline(yintercept=0, alpha=.8) +
        theme_grey(base_size=18) + transparent_opts +
        theme(legend.position="bottom") + ylim(c(-3, 3))
    ggsave(paste0('spi_vgplotmean_', spi_plot_period, '_transparent.png'), width=width, 
           height=height, dpi=dpi,
           bg='transparent')

    ggplot(spi_means[spi_means$spi_period == spi_plot_period, ]) +
        geom_line(aes(date, mean_spi, colour=Continent)) +
        xlab('Year') + scale_x_date(labels=date_format("'%y"), breaks=scale_breaks) +
        ylab('Standardized Precipitation Index (SPI)') +
        facet_wrap(~sitename_pretty) +
        geom_hline(yintercept=seq(-2, 2), alpha=.3) +
        geom_hline(yintercept=0, alpha=.8) +
        theme_grey(base_size=18) + theme(legend.position="bottom") +
        ylim(c(-3, 3))
    ggsave(paste0('spi_vgplotmean_', spi_plot_period, '.png'), width=width, 
           height=height, dpi=dpi)
}
