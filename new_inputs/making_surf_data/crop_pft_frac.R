# make crop PFT fraction data using crop_calendar.RData
library(ncdf4);library(maps); library(fields)
source('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/clm2_data/surfdata_map/making_surfdata/code/get_geo.R')

setwd('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/clm2_data/surfdata_map/making_surfdata/')
load('crop_calendar.RData')

lon = seq(-180, 177.5, 2.5) # 144
lat = seq(-90, 90, 2) # 91
PFT_frac = array(0, dim=c(length(lon), length(lat), 25)) # 25 pft types
#for(i in 18:25) plot.field(crop_planting_date[,,i], lon, lat) 
#for(i in 18:25) plot.field(crop_harvest_date[,,i], lon, lat) 

tmp = nc_open('surfdata_1.9x2.5_mp24_simyr2000_c130419_Sadiq.nc', write = T)
pct_pft = ncvar_get(tmp, 'PCT_PFT') # dim ?
lon.frame = ncvar_get(tmp, 'LONGXY') # 0 to 357.5
lat.frame = ncvar_get(tmp, 'LATIXY') # -90 to 90

# create PFT fraction maps and save them in surf_data_regrid_2000.RData
crop_index = c(18:25)
for(k in crop_index){
crop.harvest = crop_harvest_date[,,k]
plot.field(crop.harvest, lon, lat)
crop.frac = array(data = NA, dim=c(length(lon), length(lat)))
  for(i in 1:144)
    for(j in 1:91){
      if(!is.na(crop.harvest[i,j])) crop.frac[i,j] = 0.01 # when there is sig crop in this grid cell, set frac to be 0.01 (to satisfy the criteria to simulate physiology)
      else crop.frac[i,j] = 0
    }
plot.field(crop.frac, lon, lat)
print(sum(crop.frac*100, na.rm=TRUE)) # fraction of grid cells modified

PFT_frac[,,k] = crop.frac
}

# change the output from -180,180 to 0,357.5
PFT_frac.tmp = PFT_frac # PCT_PFT unit is percent
PFT_frac.tmp[1:72,,] = PFT_frac[73:144,,]
PFT_frac.tmp[73:144,,] = PFT_frac[1:72,,]
PFT_frac = PFT_frac.tmp*100 # unit: percent, dim [144,91,25]

# prepare data for interpolation
lon = seq(0, 357.5, 2.5) # 144
lonlon = replicate(length(lat), lon)
latlat = t(replicate(length(lon), lat))
lon.array = array(lonlon)
lat.array = array(latlat)
latlon = cbind(lat.array, lon.array)

# interpolate PFT_frac from lon x lat [144,91] to lon.frame x lat.frame [144,96]
for(k in crop_index){
  pft.frac.array = array(PFT_frac[,,k])
  pct_pft[,,k] = spavg(spdata = pft.frac.array, latlon = latlon, lat.frame = lat.frame[1,], lon.frame = lon.frame[,1])
}
for(k in crop_index) plot.field(pct_pft[,,k], lon.frame[,1], lat.frame[1,], Pacific.centric = T)
pct_pft[,86,18:25] = 0
for(k in crop_index) plot.field(pct_pft[,,k], lon.frame[,1], lat.frame[1,], Pacific.centric = T)

# put outputs in the surfdata
ncvar_put(nc = tmp, varid = 'PCT_PFT', vals = pct_pft)

nc_close(tmp)

