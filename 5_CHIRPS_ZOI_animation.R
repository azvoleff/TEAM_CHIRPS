###############################################################################
# Plots basic stats from CHIRPS pentad precipitation data w/in TEAM ZOIs
###############################################################################

library(ggplot2)
library(lubridate)
library(plyr)
library(scales) # for date_breaks
library(raster)
library(rasterVis)
library(animation)

in_base_dir <- 'H:/Data/CHIRPS'
in_folder <- file.path(in_base_dir, 'global_pentad')
stopifnot(file_test('-d', in_folder))

width <- 10
height <- 7.5
dpi <- 300
size_scale <- 2

# Note the below code is INCLUSIVE of the start date
chirps_start_date <- as.Date('1981/1/1')
# Note the below code is EXCLUSIVE of the end date
chirps_end_date <- as.Date('2014/5/1')
yrs <- seq(year(chirps_start_date), year(chirps_end_date))
# Have 72 pentads per year (12 months per year, 6 pentads per month)
yrs_rep <- rep(yrs, each=12*6)
days <- c(1, 6, 11, 16, 21, 26)
mths <- rep(paste(rep(seq(1, 12), each=6), days, sep='/'), length(yrs))
dates <- as.Date(paste(yrs_rep, mths, sep='/'))

transparent_opts <- theme(axis.text=element_text(colour='white'), 
                          legend.title=element_text(colour='white'), 
                          legend.text=element_text(colour='white'), 
                          plot.title=element_text(colour='white'), 
                          plot.background=element_rect(fill='transparent', colour=NA),
                          legend.background=element_rect(fill='transparent', colour=NA))

nak_anom <- brick(file.path(in_folder, 'pentad_NAK_ZOI_ppt_anom.envi'))

load('C:/Users/azvoleff/Code/TEAM/teamlucc_scripts/AOIs/ZOI_CSA_PA_RData/NAK_ZOI_CSA_PA.RData')
aoi_tr <- spTransform(aois, CRS(proj4string(nak_anom)))
aoi_tr$ID <- row.names(aoi_tr)
aoi_tr$Area <- aoi_tr$label
aoi_tr@data$id <- rownames(aoi_tr@data)
aoi_points <- fortify(aoi_tr, region="id")
aoi_df <- join(aoi_points, aoi_tr@data, by="id")

make_frame <- function(period_num) {
    gplot(nak_anom[[period_num]]) + geom_tile(aes(fill=value)) +
        coord_fixed() + 
        theme_bw(base_size=18) +
        theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
              axis.title.x=element_blank(), axis.title.y=element_blank(),
              panel.background=element_blank(), panel.border=element_blank(),
              panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              plot.background=element_blank(), axis.ticks=element_blank(),
              plot.margin=unit(c(.1, .1, .1, .1), 'cm'),
              legend.key.size=unit(3, "line"),
              legend.key.height=unit(1.5, "line")) +
        geom_path(aes(long, lat, group=ID, linetype=Area, colour=Area), data=aoi_df, 
                  size=1, alpha=.7) +
        geom_path(aes(long, lat, group=ID), data=aoi_df, 
                  size=2, alpha=.2) +
        ggtitle(dates[period_num]) +
        scale_fill_gradientn('Precipitation anomaly', limits=c(-3, 3), 
                             colours=c('red', 'orange', 'white', 'blue', 'purple'),
                             labels=c('-3 sd',
                                      '-2 sd',
                                      '-1 sd',
                                      '0 sd',
                                      '1 sd',
                                      '2 sd',
                                      '3 sd')) +
        transparent_opts
}
make_frame(1)
ggsave('NAK_precip_anom_1981-01-01.png', bg='transparent')
make_frame(1369)
ggsave('NAK_precip_anom_2000-01-01.png', bg='transparent')
make_frame(2397)
ggsave('NAK_precip_anom_2014-04-11.png', bg='transparent')
make_frame(2398)
ggsave('NAK_precip_anom_2014-04-16.png', bg='transparent')
make_frame(2399)
ggsave('NAK_precip_anom_2014-04-21.png', bg='transparent')
make_frame(2400)
ggsave('NAK_precip_anom_2014-04-26.png', bg='transparent')

height <- 4
width <- 4
dpi <- 200
ani.options(ffmpeg="C:/Program Files/ffmpeg/bin/ffmpeg.exe",
            ani.width=width*dpi, ani.height=height*dpi, verbose=TRUE,
            outdir=getwd())
saveVideo({
            #for (period_num in 1369:nlayers(nak_anom)) {
            for (period_num in 1369:1400) {
                make_frame(period_num)
            }
          },
          video.name='NAK_anom_anim.mov',
          # For Mac:
          other.opts='-preset slow -r 25 -pix_fmt yuv420p')
          # For PC:
          #other.opts='-f mp4 -preset slow -r 25')
