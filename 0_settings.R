library(ggplot2)
library(stringr)

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

overwrite <- TRUE

sites <- read.csv(file.path(prefix, 'TEAM/Sitecode_Key/sitecode_key.csv'))
sitecodes <- sites$sitecode

chirps_NA_value <- -9999
#dataset <- 'pentad'
dataset <- 'monthly'

zoi_folder <- file.path(prefix, "TEAM", "ZOIs")
in_folder <- file.path(prefix, "CHIRPS-2.0", paste0('global-', dataset))
out_folder <- file.path(prefix, "CHIRPS-2.0", paste0("TEAM-", dataset))
stopifnot(file_test('-d', in_folder))
stopifnot(file_test('-d', out_folder))

tifs <- dir(in_folder, pattern='.tif$')

datestrings <- gsub('.tif', '', (str_extract(tifs, '[0-9]{4}\\.[0-9]{2}.tif$')))
years <- as.numeric(str_extract(datestrings, '^[0-9]{4}'))
# The subyears strings are numeric codes referring to either pentads or months, 
# depending on the dataset chosen.
subyears <- as.numeric(str_extract(datestrings, '[0-9]{2}$'))

datestrings <- datestrings[order(years, subyears)]
tifs <- tifs[order(years, subyears)]

datestrings <- gsub('[.]', '', datestrings)
start_date <- datestrings[1]
end_date <- datestrings[length(datestrings)]

if (dataset == 'monthly') {
    dates <- as.Date(paste0(datestrings, '01'), '%Y%m%d')
    num_periods <- 12
} else if (dataset == 'pentad') {
    # Have 72 pentads per year (12 months per year, 6 pentads per month)
    days <- c(1, 6, 11, 16, 21, 26)
    mths <- rep(paste(rep(seq(1, 12), each=6), days, sep='/'), length(years))
    dates <- as.Date(paste(years, mths, sep='/'))
    num_periods <- 72
}

# This is the projection of the CHIRPS files, read from the .hdr files 
# accompanying the data
s_srs <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'

width <- 10
height <- 7.5
dpi <- 300
transparent_opts <- theme(legend.position="bottom", 
                          axis.text=element_text(colour='white'), 
                          axis.title.x=element_text(colour='white'), 
                          axis.title.y=element_text(colour='white'), 
                          legend.background=element_rect(fill='transparent', colour=NA),
                          legend.title=element_text(colour='white'), 
                          legend.text=element_text(colour='white'), 
                          plot.background=element_rect(fill='transparent', colour=NA))
