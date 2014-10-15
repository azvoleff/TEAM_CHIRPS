###############################################################################
# Crops CHIRPS pentad or monthly precipitation data to cover the spatial extent 
# of the ZOI/CSA/PA boundary of each team site.
###############################################################################

library(rgdal)
library(raster)
library(stringr)
library(gdalUtils)
library(rgeos)
library(teamlucc)

n_cpus <- 4
overwrite <- TRUE

sites <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')
sitecodes <- sites$sitecode

#dataset <- 'pentad'
dataset <- 'monthly'

in_base_dir <- 'D:/CHIRPS_Originals'
out_base_dir <- 'H:/Data/CHIRPS'
in_folder <- file.path(in_base_dir, paste0('global_', dataset))
out_folder <- file.path(out_base_dir, paste0('global_', dataset))
stopifnot(file_test('-d', in_folder))
stopifnot(file_test('-d', out_folder))

bils <- dir(in_folder, pattern='.bil$')

datestrings <- gsub('.bil', '', (str_extract(bils, '[0-9]{6}.bil$')))
years <- as.numeric(str_extract(datestrings, '^[0-9]{4}'))
# The subyears strings are numeric codes referring to either pentads or months, 
# depending on the dataset chosen.
subyears <- as.numeric(str_extract(datestrings, '[0-9]{2}$'))

datestrings <- datestrings[order(years, subyears)]
bils <- bils[order(years, subyears)]

product <- unique(str_extract(bils, '^v[0-9]*p[0-9]*chirps'))
stopifnot(length(product) == 1)

# Build a VRT with all dates in a single layer stacked VRT file (this stacks 
# the bils, but with delayed computation - the actual cropping and stacking 
# computations won't take place until the gdalwarp line below that is run for 
# each aoi)
vrt_file <- file.path(out_folder, paste0(product, '_', dataset, '_stack.vrt'))
gdalbuildvrt(file.path(in_folder, '*.bil'), vrt_file, separate=TRUE, 
             overwrite=overwrite)

# This is the projection of the CHIRPS files, read from the .hdr files 
# accompanying the data
s_srs <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'

for (sitecode in sitecodes) {
    timestamp()
    message('Processing ', sitecode, '...')

    load(file.path('H:/Data/TEAM/ZOI_CSA_PAs',
                   paste0(sitecode, '_ZOI_CSA_PA.RData')))
    aoi <- gConvexHull(aois)
    aoi <- spTransform(aoi, CRS(utm_zone(aoi, proj4string=TRUE)))
    aoi <- gBuffer(aoi, width=5000)
    aoi <- spTransform(aoi, CRS(s_srs))
    te <- as.numeric(bbox(aoi))

    dstfile <- file.path(out_folder,
                          paste0(product, '_', dataset, '_', sitecode, '_', 
                                 datestrings[1], '-', 
                                 datestrings[length(datestrings)], '.dat'))

    # Crop bils for this site
    gdalwarp(vrt_file, dstfile, s_srs=s_srs, t_srs=s_srs, te=te, 
             multi=TRUE, wo=paste0("NUM_THREADS=", n_cpus), 
             overwrite=overwrite)

}
