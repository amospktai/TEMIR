################################################################################
### Module for calculating statistics of default hourly output data in nc files
################################################################################

# Function to find daily statistics (e.g., mean, max, min) from default hourly output data in nc files:

library(ncdf4)
library(fields)
library(maps)
load('~/Desktop/TEMIR/PFT_frac.RData')
load("~/Dropbox/TGABI/R/functions.RData")

f_daily_stat = function(hist_name, start_date, end_date, varid, FUN=mean, hist_data_dir='~/Desktop/TEMIR/Cluster_outputs/') {
  
  # This function requires external functions: make.date.vec (from "tools.R")
  # This function requires R packages: ncdf4
  # There cannot be any missing dates in between "start_date" and "end_date".
  
  # Vector of simulation dates:
  date_vec = make.date.vec(start.date=start_date, end.date=end_date)
  # Number of simulation days:
  n_day = length(date_vec)
  
  # Define new data array:
  filename = paste0(hist_data_dir, hist_name, '/hist_grid_', as.character(start_date), '.nc')
  nc = nc_open(filename)
  lon = ncvar_get(nc, varid='lon')
  lat = ncvar_get(nc, varid='lat')
  pft = ncvar_get(nc, varid='pft')
  hour = ncvar_get(nc, varid='hour')
  nc_close(nc)
  hist_daily = array(NaN, dim=c(length(lon), length(lat), n_day, length(pft)))
  
  # Looping over days:
  for (d in 1:n_day) {
    
    print(paste0('Processing ', varid, ' for date = ', as.character(date_vec[d])), quote=FALSE)
    
    # Extract nc file:
    hist_hourly = array(NaN, dim=c(length(lon), length(lat), length(pft), length(hour)))
    filename = paste0(hist_data_dir, hist_name, '/hist_grid_', as.character(date_vec[d]), '.nc')
    nc = nc_open(filename)
    hist_hourly[,,,] = ncvar_get(nc, varid=varid)
    nc_close(nc)
    
    # Find daily statistics:
    hist_daily[,,d,] = apply(hist_hourly, MARGIN=1:3, FUN=FUN, na.rm=TRUE)
    
  }
  
  return(hist_daily)
  
}

#global_2010_LAI = f_daily_stat('hist_data_2010_global_LAI', 20100101, 20101231, 'A_can')

hist_data_dir = '~/Desktop/TEMIR/Cluster_outputs/'
hist_name = 'hist_data_2010_global_LAI'
start_date = 20100101
end_date = 20101231
date_vec = make.date.vec(start.date=start_date, end.date=end_date)
n_day = length(date_vec)


filename = paste0(hist_data_dir, hist_name, '/hist_grid_', as.character(start_date), '.nc')
nc = nc_open(filename)
lon = ncvar_get(nc, varid='lon')
lat = ncvar_get(nc, varid='lat')
pft = ncvar_get(nc, varid='pft')
hour = ncvar_get(nc, varid='hour')
nc_close(nc)

grid_area = matrix(NaN, nrow=length(lon), ncol=length(lat))
for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    grid_area[i,j] = area.latlon(lat[j]-2/2, lon[i]-2.5/2, lat[j]+2/2, lon[i]+2.5/2)
  }
}
grid_area_flip <- t(grid_area[nrow(grid_area):1,])

global_2010_LAI_total = array(NaN, dim=c(length(lon), length(lat), length(1:n_day), length(pft)))
for (d in 1:n_day) {
  for (ipft in 1:length(pft)) {
    global_2010_LAI_total[, , d, ipft] = global_2010_LAI[, , d, ipft]*grid_area[ , ]
  }
}

global_2010_LAI_GPP_daily = array(NaN, dim=c(length(lon), length(lat), length(1:n_day)))
for (d in 1:n_day) {
  global_2010_LAI_GPP_daily [, , d] = apply(global_2010_LAI_total [, , d, ] * PFT_frac, 1:2, sum, na.rm = TRUE)
}

global_2010_1.1_GPP_mean = apply (global_2010_1.1_GPP_daily, 1:2, mean) / grid_area * 1.03775
filename = '~/Desktop/CUHK/GOPT.2010_GPP/2010_GPP_mean.RData'
save (list = 'global_2010_debug_GPP_mean', file = filename)
load('~/Desktop/CUHK/GOPT.2010_GPP/2010_GPP_mean.RData')
filename = '~/Desktop/CUHK/GOPT.2010_GPP/observed_GPP_plot.RData'
save(list = 'GPP_2010', file = filename)
global_2010_total_obs_Dec = GPP_Dec * grid_area_flip * 1000000 * 31
sum_GPP_Dec = sum(global_2010_total_obs_Dec[,], na.rm = T)
load ('~/Desktop/CUHK/GOPT.2010_GPP/observed_GPP_plot.RData')
GPP_2010_diff_1.1 = global_2010_debug_GPP_mean - global_2010_1.1_GPP_mean 
plot.field(spdata = GPP_2010_diff_1.1, lon.map = lon, lat.map = lat, type = 'sign', Pacific.centric = F)
GPP_2010_diff_1.1_perc = (global_2010_debug_GPP_mean - global_2010_1.1_GPP_mean)/global_2010_1.1_GPP_mean * 100
#GPP_2010_diff_1.1_perc <- ifelse(test = GPP_2010_diff_1.1_perc >= 100, yes = NaN, no = GPP_2010_diff_1.1_perc)
plot.field(spdata = GPP_2010_diff_1.1_perc, lon.map = lon, lat.map = lat, type = 'sign', zlim = c(-100, 100), Pacific.centric = F)
# For GOPT observed data
f_monthly_stat = function(start_date, end_date, varid, observed_data_dir='~/Desktop/CUHK/GOPT.2010_GPP/') {
    
    # Vector of simulation dates:
    date_vec = make.date.vec(start.date=start_date, end.date=end_date)
    # Number of simulation days:
    n_day = length(date_vec)
    
    # Define new data array:
    filename = paste0(observed_data_dir, '/GOPT.gosat.', 201012, '.2.5x2.nc')
    nc = nc_open(filename)
    lon = ncvar_get(nc, varid='lon')
    lat = ncvar_get(nc, varid='lat')
    GPP_Dec = ncvar_get(nc, varid='GPO')
    nc_close(nc)
    hist_daily = array(NaN, dim=c(length(lon), length(lat), 12))
    
sum_GPP_2010 = sum_GPP_Jan + sum_GPP_Feb + sum_GPP_Mar + sum_GPP_Apr + sum_GPP_May + sum_GPP_Jun + sum_GPP_Jul + sum_GPP_Aug + sum_GPP_Sep + sum_GPP_Oct + sum_GPP_Nov + sum_GPP_Dec
GPP_2010 = (GPP_Jan*31 + GPP_Feb*28 + GPP_Mar*31 + GPP_Apr*30 + GPP_May*31 + GPP_Jun*30 + GPP_Jul*31 + GPP_Aug*31 + GPP_Sep*30 + GPP_Oct*31 + GPP_Nov*30 + GPP_Dec*31)/365
GPP_2010 = apply(t(GPP_2010), 2, rev)
GPP_2010 = GPP_2010[c(144:1), , drop = F]
filename = '~/Desktop/CUHK/GOPT.2010_GPP/observed_GPP_plot.RData'
save(list = 'GPP_2010', file = filename)
filename = '~/Desktop/CUHK/GOPT.2010_GPP/2010_GPP.RData'
save(list = 'sum_GPP_2010', file = filename)

    # Looping over months:
    for (d in 1:12) {
      
      print(paste0('Processing ', varid, ' for date = ', as.character(date_vec[d])), quote=FALSE)
      
      # Extract nc file:
      hist_monthly= array(NaN, dim=c(length(lon), length(lat), length(pft))
      filename = paste0('/hist_grid_', as.character(date_vec[d]), '.nc')
      nc = nc_open(filename)
      hist_monthly[,,] = ncvar_get(nc, varid=varid)
      nc_close(nc)
  
    return(hist_daily)
    
  }

################################################################################

# Function to find monthly mean from default hourly output data in nc files, including the option to sum over all PFTs:

f_monthly_mean = function(hist_name, start_date, end_date, varid, PFT_sum=FALSE, PFT_frac=NULL, hist_data_dir='~/TGABI/Tai/TEMIR/hist_data/') {
   
   # This function requires external functions: make.date.vec (from "tools.R")
   # This function requires R packages: ncdf4
   # There cannot be any missing dates in between the first and last days of a month.
   # Therefore, "start_date" should always be the first day of a given month, and "end_date" should be the last day of a given month.
   # If "PFT_sum=TRUE", weighted sum over all PFTs (weighted by "PFT_frac") will be calculated.
   # "PFT_frac": dim1 = lon; dim2 = lat; dim3 = pft
   
   # Vector of simulation months:
   date_vec = make.date.vec(start.date=start_date, end.date=end_date)
   month_vec = unique(floor(date_vec/1e2))*1e2 + 1
   # Number of simulation months:
   n_month = length(month_vec)
   
   # Define new data array:
   filename = paste0(hist_data_dir, hist_name, '/hist_grid_', as.character(start_date), '.nc')
   nc = nc_open(filename)
   lon = ncvar_get(nc, varid='lon')
   lat = ncvar_get(nc, varid='lat')
   pft = ncvar_get(nc, varid='pft')
   hour = ncvar_get(nc, varid='hour')
   if (PFT_sum) hist_monthly = array(NaN, dim=c(length(lon), length(lat), n_month)) else hist_monthly = array(NaN, dim=c(length(lon), length(lat), n_month, length(pft)))
   
   # Looping over months:
   for (m in 1:n_month) {
      
      print(paste0('Processing ', varid, ' for month = ', substr(as.character(month_vec[m]), start=1, stop=6)), quote=FALSE)
      
      # Generate hourly data array for each month:
      date_vec_sub = date_vec[which(floor(date_vec/1e2) == floor(month_vec[m]/1e2))]
      hist_hourly = array(NaN, dim=c(length(lon), length(lat), length(pft), length(date_vec_sub)*length(hour)))
      
      # Looping over days:
      for (d in 1:length(date_vec_sub)) {
         # Extract nc file:
         ind_hr = ((d - 1)*length(hour) + 1):((d - 1)*length(hour) + length(hour))
         filename = paste0(hist_data_dir, hist_name, '/hist_grid_', as.character(date_vec_sub[d]), '.nc')
         nc = nc_open(filename)
         hist_hourly[,,,ind_hr] = ncvar_get(nc, varid=varid)
         nc_close(nc)
      }
        
      # Find monthly mean:
      if (PFT_sum) {
         hist_monthly_PFT = apply(hist_hourly, MARGIN=1:3, FUN=mean, na.rm=TRUE)
         hist_monthly[,,m] = apply(hist_monthly_PFT*PFT_frac, MARGIN=1:2, FUN=sum, na.rm=TRUE)
      } else {
         hist_monthly[,,m,] = apply(hist_hourly, MARGIN=1:3, FUN=mean, na.rm=TRUE)
      }
      
   }
   
   return(hist_monthly)
   
}

# timestamp()
# out2 = f_monthly_mean(hist_name='control_global_2010', start_date=20100601, end_date=20100831, varid='A_can')
# timestamp()
# # It requires ~100 seconds to finish 3 months.


################################################################################

# Function to find monthly mean of any daily statistic (e.g., daily max, min) from default hourly output data in nc files:

f_monthly_mean_stat = function(hist_name, start_date, end_date, varid, FUN=max, hist_data_dir='~/TGABI/Tai/TEMIR/hist_data/') {

   # This function requires external functions: make.date.vec (from "tools.R"), f_daily_stat
   # This function requires R packages: ncdf4
   # There cannot be any missing dates in between the first and last days of a month.
   # Therefore, "start_date" should always be the first day of a given month, and "end_date" should be the last day of a given month.

   # Vector of simulation months:
   date_vec = make.date.vec(start.date=start_date, end.date=end_date)
   month_vec = unique(floor(date_vec/1e2))*1e2 + 1
   # Number of simulation months:
   n_month = length(month_vec)

   # Define new data array:
   filename = paste0(hist_data_dir, hist_name, '/hist_grid_', as.character(start_date), '.nc')
   nc = nc_open(filename)
   lon = ncvar_get(nc, varid='lon')
   lat = ncvar_get(nc, varid='lat')
   pft = ncvar_get(nc, varid='pft')
   hour = ncvar_get(nc, varid='hour')
   hist_monthly = array(NaN, dim=c(length(lon), length(lat), n_month, length(pft)))

   # Looping over months:
   for (m in 1:n_month) {
      # Find daily statistic for each month:
      date_vec_sub = date_vec[which(floor(date_vec/1e2) == floor(month_vec[m]/1e2))]
      hist_daily = f_daily_stat(hist_name=hist_name, start_date=date_vec_sub[1], end_date=tail(date_vec_sub, 1), varid=varid, FUN=FUN, hist_data_dir=hist_data_dir)
      # Find monthly mean of daily statistic:
      hist_monthly[,,m,] = apply(hist_daily, MARGIN=c(1,2,4), FUN=mean, na.rm=TRUE)
   }

   return(hist_monthly)

}

# timestamp()
# out1 = f_monthly_mean(hist_name='control_global_2010', start_date=20100601, end_date=20100831, varid='A_can')
# timestamp()
# # It requires ~390 seconds to finish 3 months.


################################################################################
### End of module
################################################################################
