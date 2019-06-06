################################################################################
### Module for calculating statistics of default hourly output data in nc files
################################################################################

# Function to find daily statistics (e.g., mean, max, min) from default hourly output data in nc files:

f_daily_stat = function(hist_name, start_date, end_date, varid, FUN=mean, hist_data_dir='~/TGABI/Tai/TEMIR/hist_data/') {
   
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
