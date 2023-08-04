################################################################################
### Module for calculating statistics of default hourly output data in nc files
### or in "hist_grid" in the debugging mode.
################################################################################

# Last updated: Apr 2020 (Tai)

################################################################################
### Functions for handling nc data:
################################################################################

library(ncdf4)

# Function to find daily statistics (e.g., mean, max, min) from default hourly output data in nc files or in "hist_grid" in the debugging mode.

f_daily_stat = function(hist_data_dir, start_date, end_date, varid, FUN=mean) {
  
  # This function requires external functions: make.date.vec (from "tools.R")
  # This function requires R packages: ncdf4
  # There cannot be any missing dates in between "start_date" and "end_date".
  
  # Vector of simulation dates:
  date_vec = make.date.vec(start.date=start_date, end.date=end_date)
  # Number of simulation days:
  n_day = length(date_vec)
  
  # Define new data array:
  filename = paste0(hist_data_dir, 'hist_grid_', as.character(start_date), '.nc')
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
    filename = paste0(hist_data_dir, 'hist_grid_', as.character(date_vec[d]), '.nc')
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

f_monthly_mean = function(hist_data_dir, start_date, end_date, varid, PFT_sum=FALSE, PFT_frac=NULL) {
   
   # This function requires external functions: make.date.vec (from "tools.R")
   # This function requires R packages: ncdf4
   # There cannot be any missing dates in between the first and last days of a month.
   # Therefore, "start_date" should always be the first day of a given month, and "end_date" should be the last day of a given month.
   # If "PFT_sum=TRUE", weighted sum over all PFTs (weighted by "PFT_frac") will be calculated.
   # "PFT_frac": dim1 = lon; dim2 = lat; dim3 = pft
   # In most cases, "PFT_frac" here should be the entirety or a subset of the variable "PFT_frac" stored in surface data file "TEMIR_inputs/processed_surf_data/surf_data_regrid_2000.RData".
   
   # Vector of simulation months:
   date_vec = make.date.vec(start.date=start_date, end.date=end_date)
   month_vec = unique(floor(date_vec/1e2))*1e2 + 1
   # Number of simulation months:
   n_month = length(month_vec)
   
   # Define new data array:
   filename = paste0(hist_data_dir, 'hist_grid_', as.character(start_date), '.nc')
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
         filename = paste0(hist_data_dir, 'hist_grid_', as.character(date_vec_sub[d]), '.nc')
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

f_monthly_mean_stat = function(hist_data_dir, start_date, end_date, varid, FUN=max) {

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
   filename = paste0(hist_data_dir, 'hist_grid_', as.character(start_date), '.nc')
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
      hist_daily = f_daily_stat(hist_data_dir=hist_data_dir, start_date=date_vec_sub[1], end_date=tail(date_vec_sub, 1), varid=varid, FUN=FUN)
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
### Functions for analyzing "hist_grid" in debugging mode:
################################################################################

# Function to find daily mean (weighted by PFT fraction over all PFTs) for variables "A_can", "R_can" or "g_can" from "hist_grid":

f_daily_PFT_sum = function(hist_grid, var.name='A_can') {
   # Only works for variables "A_can", "R_can" or "g_can".
   if (var.name == 'A_can') {
      ivar = 1
   } else if (var.name == 'R_can') {
      ivar = 2
   } else if (var.name == 'g_can') {
      ivar = 3
   } else { stop('Variable defined is not valid for this function!') }
   var_daily = array(NaN, dim=c(length(ind_lon), length(ind_lat), n_day_sim))
   for (n_i in 1:length(ind_lon)) {
      for (n_j in 1:length(ind_lat)) {
         PFT_frac_ij = PFT_frac[ind_lon[n_i],ind_lat[n_j],]
         for (d in 1:n_day_sim) {
            # Find daily mean first:
            var_daily_PFT = apply(hist_grid[n_i,n_j,d,,,ivar], MARGIN=1, FUN=mean, na.rm=TRUE)
            # Find weighted sum over all PFTs:
            var_daily[n_i,n_j,d] = sum(PFT_frac_ij*var_daily_PFT, na.rm=TRUE)
         }
      }
   }
   return(var_daily)
}

# Function to find weighted sum over all PFTs for any PFT-level data array:

f_PFT_sum = function(X, PFT_frac=PFT_frac, PFT_dim) {
   # "X" is any multidimensional data array with at least three dimensions for at least longitude, latitude and PFT. The first dimension has to be longitude, and the second dimension has to be latitude. Maximum number of dimensions is six.
   # "PFT_frac" is a three-dimensional array with 1st dim being longitude, 2nd dim being latitude, and 3rd dim being PFT, containing the fractional coverage of different PFTs over the entire grid cell (not land only).
   # "lon_dim", "lat_dim" and "PFT_dim" are the n-th dimensions for longitude, latitude and PFT, respectiviely, of "X".
   dim_X = dim(X)
   if (PFT_dim > length(dim_X) | PFT_dim == 1 | PFT_dim == 2) stop("PFT_dim is not correctly specified!")
   if (length(dim_X) == 3) {
      X_out = apply(X*PFT_frac, MARGIN=1:2, FUN=sum, na.rm=TRUE)
   } else if (length(dim_X) == 4) {
      X_out = array(NaN, dim=dim_X[-PFT_dim])
      extra_dim = (1:4)[-c(1, 2, PFT_dim)]
      for (k in 1:dim_X[extra_dim]) {
         if (PFT_dim == 3) X_sub = X[,,,k] else X_sub = X[,,k,]
         X_out[,,k] = apply(X_sub*PFT_frac, MARGIN=1:2, FUN=sum, na.rm=TRUE)
      }
   } else if (length(dim_X) == 5) {
      X_out = array(NaN, dim=dim_X[-PFT_dim])
      extra_dim = (1:5)[-c(1, 2, PFT_dim)]
      for (k in 1:dim_X[extra_dim[1]]) {
         for (l in 1:dim_X[extra_dim[2]]) {
            if (PFT_dim == 3) X_sub = X[,,,k,l] else if (PFT_dim == 4) X_sub = X[,,k,,l] else X_sub = X[,,k,l,]
            X_out[,,k,l] = apply(X_sub*PFT_frac, MARGIN=1:2, FUN=sum, na.rm=TRUE)
         }
      }
   } else if (length(dim_X) == 6) {
      X_out = array(NaN, dim=dim_X[-PFT_dim])
      extra_dim = (1:6)[-c(1, 2, PFT_dim)]
      for (k in 1:dim_X[extra_dim[1]]) {
         for (l in 1:dim_X[extra_dim[2]]) {
            for (m in 1:dim_X[extra_dim[3]]) {
               if (PFT_dim == 3) X_sub = X[,,,k,l,m] else if (PFT_dim == 4) X_sub = X[,,k,,l,m] else if (PFT_dim == 5) X_sub = X[,,k,l,,m] else X_sub = X[,,k,l,m,]
               X_out[,,k,l,m] = apply(X_sub*PFT_frac, MARGIN=1:2, FUN=sum, na.rm=TRUE)
            }
         }
      }
   } else {
      stop("Input data array has more than six dimensions!")
   }
   return(X_out)
}

################################################################################
### End of module
################################################################################
