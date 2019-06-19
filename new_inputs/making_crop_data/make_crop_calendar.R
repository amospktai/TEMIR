# read in crop calendar maps from 7 separate files (Sacks et al., 2010), interpolate to [144,91] resolution and store them in an R file
library(ncdf4);library(maps); library(fields)
source('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/clm2_data/surfdata_map/making_surfdata/code/get_geo.R')

# interpolate plant and harvest dates from dim=[720, 360] to [144,91]
lon.frame = seq(-180, 177.5, 2.5) # 144
lat.frame = seq(-90, 90, 2) # 91

planting.date = list(lon.frame, lat.frame)
harvest.date = list(lon.frame, lat.frame)

setwd('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/clm2_data/surfdata_map/making_surfdata/Sacks/')

crop.types = c('Rice', 'Rice.2', 'Maize', 'Maize.2', 'Soybeans', 'Wheat', 'Wheat.Winter')
for (crops in crop.types){
  tmp = nc_open(paste0(crops, '.crop.calendar.nc'))
  lon = ncvar_get(tmp, 'longitude') # 720
  lat = ncvar_get(tmp, 'latitude') # 360
  planting = ncvar_get(tmp, 'plant') # dimension is [720, 360]
  harvest = ncvar_get(tmp, 'harvest') # [720, 360]
  nc_close(tmp)
  
  plot.field(planting, lon, lat) # original
  #plot.field(harvest, lon ,lat)
  
  # prepare lon, lat and data for interpolation
  lonlon = replicate(length(lat), lon)
  latlat = t(replicate(length(lon), rev(lat)))
  planting.new = planting[,ncol(planting):1]
  harvest.new   = harvest[,ncol(harvest):1]
  
  planting.array = array(planting.new)
  harvest.array = array(harvest.new)
  
  lon.array = array(lonlon)
  lat.array = array(latlat)
  latlon = cbind(lat.array, lon.array)
  # interpolation, takes a while
  interpol.plant = spavg(planting.array, latlon, lat.frame, lon.frame)
  interpol.harvest = spavg(harvest.array, latlon, lat.frame, lon.frame)
  # plot results of interpolation
  plot.field(interpol.plant, lon.frame, lat.frame) # after interpolation
  plot.field(interpol.harvest, lon.frame, lat.frame)
  
  planting.date[[crops]] = interpol.plant
  harvest.date[[crops]] = interpol.harvest
}
# store plang and harvest dates into an RData
pft_num = 25
crop_planting_date =array(NA, dim=c(length(lon.frame), length(lat.frame), pft_num))
crop_planting_date[,,18] = planting.date$Maize # Maize
crop_planting_date[,,19] = planting.date$Maize.2 # Maize 2
crop_planting_date[,,20] = planting.date$Rice # Rice
crop_planting_date[,,21] = planting.date$Rice.2 # Rice 2
crop_planting_date[,,22] = planting.date$Wheat # spring cereal -- wheat
crop_planting_date[,,23] = planting.date$Wheat.Winter # winter cereal -- wheat winter
crop_planting_date[,,24:25] = planting.date$Soybeans # soybean - soybean

crop_harvest_date =array(NA, dim=c(length(lon.frame), length(lat.frame), pft_num))
crop_harvest_date[,,18] = harvest.date$Maize # Maize
crop_harvest_date[,,19] = harvest.date$Maize.2 # Maize 2
crop_harvest_date[,,20] = harvest.date$Rice # Rice
crop_harvest_date[,,21] = harvest.date$Rice.2 # Rice 2
crop_harvest_date[,,22] = harvest.date$Wheat # spring cereal -- wheat
crop_harvest_date[,,23] = harvest.date$Wheat.Winter # winter cereal -- wheat winter
crop_harvest_date[,,24:25] = harvest.date$Soybeans # soybean - soybean

setwd('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/clm2_data/surfdata_map/making_surfdata/')
save(list=c('lon.frame', 'lat.frame', 'crop_planting_date', 'crop_harvest_date'), file = "crop_calendar.RData")

#for(i in 18:25) plot.field(crop_harvest_date[,,i], lon.frame, lat.frame)
#for(i in 18:25) plot.field(crop_planting_date[,,i], lon.frame, lat.frame)
#for(i in 18:25) plot.field(crop_harvest_date[,,i] - crop_planting_date[,,i], lon.frame, lat.frame, type = 'sign')

