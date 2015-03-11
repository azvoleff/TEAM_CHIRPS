###############################################################################
# Crops CHIRPS pentad or monthly precipitation data to cover the spatial extent 
# of the ZOI/CSA/PA boundary of each team site.
###############################################################################

source("0_settings.R")

library(rgdal)
library(raster)
library(stringr)
library(gdalUtils)
library(rgeos)
library(teamlucc)

warp_threads <- 4

# Build a VRT with all dates in a single layer stacked VRT file (this stacks 
# the tifs, but with delayed computation - the actual cropping and stacking 
# computations won't take place until the gdalwarp line below that is run for 
# each aoi)
vrt_file <- extension(rasterTmpFile(), 'vrt')
gdalbuildvrt(file.path(in_folder, tifs), vrt_file, separate=TRUE, 
             overwrite=TRUE)

for (sitecode in sitecodes) {
    timestamp()
    message('Processing ', sitecode, '...')

    load(file.path(prefix, 'TEAM/ZOI_CSA_PAs',
                   paste0(sitecode, '_ZOI_CSA_PA.RData')))
    aoi <- gConvexHull(aois)
    aoi <- spTransform(aoi, CRS(utm_zone(aoi, proj4string=TRUE)))
    aoi <- gBuffer(aoi, width=5000)
    aoi <- spTransform(aoi, CRS(s_srs))
    te <- as.numeric(bbox(aoi))
    # Round extent so that pixels are aligned properly
    te <- round(te * 20) / 20

    base_name <- file.path(out_folder,
                           paste0(sitecode, '_CHIRPS_', dataset,
                                  '_', start_date, '-', end_date))
    chirps_tif <- paste0(base_name, '.tif')

    chirps <- gdalwarp(vrt_file, chirps_tif, s_srs=s_srs, te=te, multi=TRUE, 
                       wo=paste0("NUM_THREADS=", warp_threads), overwrite=TRUE, 
                       output_Raster=TRUE)
}
