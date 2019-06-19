# 1, Make LAI and SAI maps for all crops using methods in Simpson 2012 (Sect. 5, Table 3 and S7.1)
# and crop calendar data in crop_calendar.RData
# 2, interpolate the resulting maps and store them in surfdata_1.9x2.5_mp24_simyr2000_c130419_Sadiq.nc
# interpolation takes a while, up to 20 minutes

source('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/clm2_data/surfdata_map/making_surfdata/code/get_geo.R')
library(fields); library(ncdf4)

# source the functions needed
source('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/clm2_data/surfdata_map/making_surfdata/code/functions_LAI.R')

# get plant and harvest dates, as start and end of the growing season
setwd('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/clm2_data/surfdata_map/making_surfdata/')
load('crop_calendar.RData') # crop_harvest_date, crop_planting_date, lon.frame, lat.frame

daily_lai = array(0, dim=c(144,91,25,365)) # 144 x 91 = lon x lat, 25 pfts, 365 days
monthly_lai = array(0, dim=c(144,91,25,12)) # 144 x 91 = lon x lat, 25 pfts, 12 months

for(i in 1:144)
  for(j in 1:91)
    for(k in 1:25)
    {
      if(is.na(crop_planting_date[i,j,k])) next
      if(is.na(crop_harvest_date[i,j,k])) next
      
      d_sgs = crop_planting_date[i,j,k] # start of growing season
      d_egs = crop_harvest_date[i,j,k] # end of growing season
      
      if(d_sgs < d_egs) # boreal/summer crops
           lai_365 = lai_summer_crop(d_sgs = d_sgs, d_egs = d_egs, LAI_min = 0, LAI_max = 3.5, l_s = 70, l_e = 20)
      # austral/winter crops
      else lai_365 = lai_winter_crop(d_sgs = d_sgs, d_egs = d_egs, LAI_min = 0, LAI_max = 3.5, l_s = 70, l_e = 20)
      lai_12 = lai_monthly_mean(input = lai_365)
      daily_lai[i,j,k,] = lai_365 # daily LAI
      monthly_lai[i,j,k,] = lai_12 # monthly LAI
    }

# change monthly_lai data from -180,180 to 0,360
lai.tmp = monthly_lai
lai.tmp[1:72,,,] = monthly_lai[73:144,,,]
lai.tmp[73:144,,,] = monthly_lai[1:72,,,]
monthly_lai = lai.tmp

# lon and lat of the above LAI
lon = seq(0, 357.5, 2.5) # 144
lat = seq(-90, 90, 2) # 91
#for(i in 18:25) plot.field(monthly_lai[,,i,8], lon.frame, lat.frame)

# get lon, lat and LAI from surfdata
setwd('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/clm2_data/surfdata_map/making_surfdata/')
tmp = nc_open('surfdata_1.9x2.5_mp24_simyr2000_c130419_Sadiq.nc', write = T)
lon.frame = ncvar_get(tmp, 'LONGXY') # 0 to 357.5
lat.frame = ncvar_get(tmp, 'LATIXY') # -90 to 90
lai.original = ncvar_get(tmp, 'MONTHLY_LAI') #[144,96,25,12]

# prepare data for interpolation, from lon x lat to lon.frame x lat.frame
lonlon = replicate(length(lat), lon)
latlat = t(replicate(length(lon), lat))
lon.array = array(lonlon)
lat.array = array(latlat)
latlon = cbind(lat.array, lon.array)

# interpolation, takes a while
for(i in 1:12)
  for(j in 18:25)
  {
    lai.array = array(monthly_lai[,,j,i])
    lai.original[,,j,i] = spavg(spdata = lai.array, latlon = latlon, lat.frame = lat.frame[1,], lon.frame = lon.frame[,1])
  }

for(k in 18:25) plot.field(lai.original[,,k,8], lon.frame[,1], lat.frame[1,], Pacific.centric = T)
lai.original[,86,18:25,] = 0

# put them in the surface data
ncvar_put(nc = tmp, varid = 'MONTHLY_LAI', vals = lai.original)

# also change the SAI, set it to be 0.5 for pft = 18 to 25 all year long
sai = ncvar_get(tmp, 'MONTHLY_SAI')
sai_new = sai
sai_new[,,18:25,] = 0.5
ncvar_put(nc = tmp, varid = 'MONTHLY_SAI', vals = sai_new)

nc_close(tmp)

# crude way of replacing, replace the processed surf data and daily lai data
# load('surf_data_regrid_2000.RData')
# for(i in 18:25) plot.field(LAI_mon_PFT[,,i,8], lon.frame, lat.frame)
# LAI_mon_PFT = monthly_lai
# save(list = c("LAI_mon_PFT", "PFT_frac", "SAI_mon_PFT", "b_psi_PFT", "exist_regrid",
#        "grid_area", "lat", "lon", "pftname", "pftnum", "psi_sat_PFT",
#        "soil_albedo", "soil_color", "theta_sat_PFT"), 
#      file = 'surf_data_regrid_2000.RData')
# 
# setwd("~/Desktop/Library/2019-05-14/")
# load('daily_LAI_CLM_2000.RData')
# LAI_day_PFT = daily_lai
# save(list = "LAI_day_PFT", "SAI_day_PFT", "day", "lat", "lon", "pftname",
#      "pftnum", file = 'daily_LAI_CLM_2000.RData')
