library(ggplot2)
library(dplyr)
library(scales) # for date_breaks

load('vg_plot_spis.RData')

# Add in continent names and countries
sites <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')
spis <- merge(spis, sites, by.x='sitecode', by.y='sitecode')

width <- 10
height <- 7.5
dpi <- 300
transparent_opts <- theme(axis.text=element_text(colour='white'), 
                          axis.title.x=element_text(colour='white'), 
                          axis.title.y=element_text(colour='white'), 
                          legend.background=element_rect(fill='transparent', colour=NA),
                          legend.title=element_text(colour='white'), 
                          legend.text=element_text(colour='white'), 
                          plot.background=element_rect(fill='transparent', colour=NA))

spis <- spis[spis$date > as.Date('2000/01/01'), ]

spis$sitename_pretty <- ordered(spis$sitename_pretty, levels=sites$sitename_pretty)

spi_means <- summarize(group_by(spis, spi_period, sitecode, date),
                       mean_spi=mean(spi),
                       sitename_pretty=sitename_pretty[1],
                       Continent=continent[1])
spi_means$spi_sign <- ifelse(spi_means$mean_spi >= 0, "wet", "dry")
spi_means$spi_period <- ordered(spi_means$spi_period, level=c(1, 6, 12, 24, 36))

spi_plot_periods <- c(1, 6, 12, 24, 36)
for (spi_plot_period in spi_plot_periods) {
    scale_breaks <- c(seq(as.Date('2000/1/1'), as.Date('2014/1/1'),
                          by='5 years'), as.Date('2014/1/1'))
    ggplot(spi_means[spi_means$spi_period == spi_plot_period, ]) +
        geom_line(aes(date, mean_spi, colour=Continent)) +
        xlab('Year') + scale_x_date(labels=date_format("'%y"), breaks=scale_breaks) +
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
