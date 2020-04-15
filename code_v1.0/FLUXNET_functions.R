################################################################################
### Module for incorporating FLUXNET site micrometeorology and site surface data
################################################################################

# All data information extracted from http://fluxnet.fluxdata.org/data/fluxnet2015-dataset/fullset-data-product/
# WARNING : UNITS SHOULD BE CHECKED FOR INDIVIDUAL VARIABLES AS THEY MIGHT CHANGE FOR EACH TIME STEP
# See list of FLUXNET variable and corresponding TEMIR inputs in FLUXNET_var_info.xlsx

library(readxl) # For reading excel/csv files

################################################################################
### Functions:
################################################################################

### Function to find time difference between FLUXNET site and model time in UTC ### NOTE : FLUXNET has no daylight saving time therefore use timezone provided in their site information .../FLUXNET/AA-Flx.xlsx
# NOTE : NOT NEEDED in v1.1 BUT LEFT FOR INSURANCE 
if (FALSE) {
   f_FLUXNET_tz_from_UTC = function(date, lon.sim = lon_sim, lat.sim = lat_sim, timezone = FLUXNET_timezone, tz.offline.flag = timezone_offline_flag){
      options(digits = 12)
      model_time = as.POSIXct(strptime(x = toString(date * 10000), format = "%Y%m%d%H%M"), tz = "UTC")
      if (!tz.offline.flag){# Getting timezone ONLINE to shift FLUXNET data appropriately
         # https://developers.google.com/maps/documentation/timezone/
         apiurl = sprintf("https://maps.googleapis.com/maps/api/timezone/%s?location=%s,%s&timestamp=%d&sensor=%s", "xml", lat.sim, lon.sim, as.numeric(model_time), "false")
         library(XML)
         tz = xmlParse(readLines(apiurl))[["string(//time_zone_id)"]]
         data_time = model_time
         attr(data_time, 'tzone') = tz
         if (model_time - data_time != 0) stop("Timezone conversion is not successful!")
      } else {# Shifting timezone from UTC to given FLUXNET timezone
         if (!is.element(timezone, OlsonNames())) stop("Invalid timezone given")
         data_time = model_time
         attr(data_time, 'tzone') = timezone
         if (model_time - data_time != 0) stop("Timezone conversion is not successful!")
      }
      output_time = as.numeric(substr(gsub('-| |:',"",toString(data_time)), start = 1, stop = 12))
      shifted_time = as.POSIXct(strptime(x = toString(output_time), format = "%Y%m%d%H%M"), tz = "UTC")
      time_shift = as.numeric(shifted_time - model_time)
      return(c(output_time, time_shift))
   }
}

################################################################################
### Function to shift model time to site time using UTC_OFFSET given by FLUXNET
f_FLUXNET_UTC_2_local = function(FLUXNET.dir, date, site.id, utc.offset = NA, out.utc.offset = FALSE) {
   # Get current FLUXNET site time difference if not known:
   if (is.na(utc.offset)) {
      # Get current FLUXNET site info:
      site_info = read_excel(path = paste0(FLUXNET.dir,'AA-Flx.xlsx'))
      site_selected = site_info[which(site_info$SITE_ID == site.id),]
      
      # Get current FLUXNET site time difference from UTC:
      if (length(which(site_selected[,"VARIABLE"] == "UTC_OFFSET")) != 0) {
         utc_offset = as.numeric(toString(site_selected[which(site_selected[,"VARIABLE"] == "UTC_OFFSET"),"DATAVALUE"]))
         
      } else if (FALSE && sign(as.numeric(toString(site_selected[which(site_selected[,"VARIABLE"] == "LOCATION_LAT"),"DATAVALUE"]))) == 1) {
         utc_offset = f_FLUXNET_tz_from_UTC(date = as.numeric(paste0(substr(date, 1, 4), '0101')), lon.sim = as.numeric(toString(site_selected[which(site_selected[,"VARIABLE"] == "LOCATION_LONG"),"DATAVALUE"])), lat.sim = as.numeric(toString(site_selected[which(site_selected[,"VARIABLE"] == "LOCATION_LAT"),"DATAVALUE"])), tz.offline.flag = FALSE)[2]
      } else if (FALSE && sign(as.numeric(toString(site_selected[which(site_selected[,"VARIABLE"] == "LOCATION_LAT"),"DATAVALUE"]))) == -1) {
         utc_offset = f_FLUXNET_tz_from_UTC(date = as.numeric(paste0(substr(date, 1, 4), '0701')), lon.sim = as.numeric(toString(site_selected[which(site_selected[,"VARIABLE"] == "LOCATION_LONG"),"DATAVALUE"])), lat.sim = as.numeric(toString(site_selected[which(site_selected[,"VARIABLE"] == "LOCATION_LAT"),"DATAVALUE"])), tz.offline.flag = FALSE)[2]
      }
      
      # Check if current FLUXNET site time was defaulted in UTC:
      utc_offset_comment = toString(site_selected[which(site_selected[,"VARIABLE"] == "UTC_OFFSET_COMMENT"),"DATAVALUE"])
      if (nchar(utc_offset_comment) > 0 & !is.na(utc_offset_comment)) {
         if (substrEND(utc_offset_comment,3) == 'UTC') utc_offset = 0
      }
      
      if (is.na(utc_offset)) stop('FLUXNET site lack time zone information!')
   } else utc_offset = utc.offset
   
   # Use time difference to shift UTC time to site local time:
   output_time = vector(length = length(date))
   for (date_ind in 1:length(date)){
      model_time = as.POSIXct(strptime(x = toString(date[date_ind] * 10000), format = "%Y%m%d%H%M"), tz = "UTC")
      data_time = model_time + utc_offset*60*60
      if (utc_offset == 0) {output_time[date_ind] = as.numeric(substr(gsub('-| |:',"",toString(data_time)), start = 1, stop = 12)) * 10000
      }else output_time[date_ind] = as.numeric(substr(gsub('-| |:',"",toString(data_time)), start = 1, stop = 12))
   }
   if (out.utc.offset) output_time = list('time' = output_time, 'time_diff' = utc_offset)
   return(output_time)
}

################################################################################
### Function to check FLUXNET data availability over the simulation days
f_FLUXNET_date_range = function(FLUXNET.dir, start.or.end.year = c('start', 'end'), whole.date.for.year = c(start_date, end_date), site.id){
   # Get current FLUXNET site:
   site_info = if (file.exists(paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx'))) read_excel(path = paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx')) else read_excel(path = paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx'))
   site_ind = match(FLUXNET_site_id, site_info$SITE_ID)
   if (is.na(site_ind)) stop("Invalid FLUXNET Site ID is given!!!")
   
   # Get current FLUXNET data start and end years:
   if (start.or.end.year == 'start'){
      if (as.numeric(substr(site_info[site_ind,4], start = 8, stop = 11)) <= as.numeric(substr(whole.date.for.year,start = 1, stop = 4))) date_in_range = TRUE else date_in_range = FALSE
   } else if (start.or.end.year == 'end'){
      if (as.numeric(substr(site_info[site_ind,5], start = 8, stop = 11)) >= as.numeric(substr(whole.date.for.year,start = 1, stop = 4))) date_in_range = TRUE
      else date_in_range = FALSE
   } else stop('Unclear input year is provided, Please try again')
   return(date_in_range)
}

################################################################################
### Function to find (lon, lat) from site.id
f_lon_lat_sim_from_FLUXNET = function(FLUXNET.dir, site.id){
   # Get current FLUXNET site:
   site_info = if (file.exists(paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx'))) read_excel(path = paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx')) else read_excel(path = paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx'))
   site_ind = match(FLUXNET_site_id, site_info$SITE_ID)
   
   # Get current FLUXNET site lon lat:
   lon_sim = as.double(site_info[site_ind,7])
   lat_sim = as.double(site_info[site_ind,6])
   remove(site_info, site_ind)
   return(list('lon_sim' = lon_sim, 'lat_sim' = lat_sim))
}

################################################################################
### Function to set appropriate FLUXNET data file directory
f_FLUXNET_file = function(FLUXNET.dir, site.id, hr.part, hourly.data = TRUE){
   # Read most optimal FLUXNET data for current simulation: ( Prefer FLUXNET FULLSET to SUBSET )
   ### NOTE : FLUXNET_nrows has extra_row as daylight saving time complicates data selecting in f_FLUXNET_and_Met_grid
   ###        In summer time, the hour shifted forwards so data is gapfilled with -9999 or duplicated data which needs to be IGNORED
   ###        In winter time, the hour shifted backwards but as f_FLUXNET_and_Met_grid is coded for 24 hours this shift is overlooked
   if (hourly.data) {
      if (length(file.exists(Sys.glob(paste0(FLUXNET.dir, 'FLX_', site.id, '*', '_FULLSET_H', '*','.csv')))) == 1) {
         data_size = 'FULLSET'
      } else data_size = 'SUBSET'
      day_hours = 24/hr.part
      extra_row = 1
      if (hr.part == 0.5){
         FLUXNET_file = Sys.glob(paste0(FLUXNET.dir, 'FLX_', site.id, '*', data_size, '_HH', '*','.csv'))
         row_hour_ratio = 2
         FLUXNET_nrows = (day_hours + extra_row) * row_hour_ratio 
      }else if (length(file.exists(Sys.glob(paste0(FLUXNET.dir, 'FLX_', site.id, '*', data_size, '_HR', '*','.csv')))) == 1){
         FLUXNET_file = Sys.glob(paste0(FLUXNET.dir, 'FLX_', site.id, '*', data_size, '_HR', '*','.csv'))
         row_hour_ratio = 1
         FLUXNET_nrows = (day_hours + extra_row) * row_hour_ratio
      } else {
         FLUXNET_file = Sys.glob(paste0(FLUXNET.dir, 'FLX_', site.id, '*', data_size, '_HH', '*','.csv'))
         row_hour_ratio = 2
         FLUXNET_nrows = (day_hours + extra_row) * row_hour_ratio
      }
   } else {
      if (length(file.exists(Sys.glob(paste0(FLUXNET.dir, 'FLX_', site.id, '*', '_FULLSET_MM', '*','.csv')))) == 1) {data_size = 'FULLSET'
      } else data_size = 'SUBSET'
      FLUXNET_file = Sys.glob(paste0(FLUXNET.dir, 'FLX_', site.id, '*', data_size, '_MM', '*','.csv'))
      FLUXNET_nrows = 1
   }
   return(list(filedir = FLUXNET_file, nrows = FLUXNET_nrows))
}

################################################################################
### Function to count occurrence of logic condition
num.elem = function(condition) {return(length(which(condition)))}

################################################################################
### Function to count get last characters of string
substrEND = function(x, n) {substr(x, nchar(x)-n+1, nchar(x))}

################################################################################
### Function to find row of current date in FLUXNET data
f_FLUXNET_row_skip = function(FLUXNET.dir, current.date, direction, rangeloc = NA, site.id, FLUXNET.nrows, hourly.data = TRUE) {
   # Get current FLUXNET data and date:
   site_info = read_excel(path = paste0(FLUXNET.dir, 'FLUXNET_site_info.xlsx'))
   site_ind = match(site.id, site_info$SITE_ID)
   FLUXNET_data_start_date = as.numeric(substr(site_info[site_ind,"2015_DATA_START"], start = 8, stop = 11))
   FLUXNET_data_end_date = as.numeric(substr(site_info[site_ind, "2015_DATA_END"], start = 8, stop = 11))
   remove(site_info, site_ind)
   
   # Skip to current date for FLUXNET data:
   if (hourly.data){
      if (direction == 'forward') {
         FLUXNET_skip_years = (signif(current.date, digits = 4) - FLUXNET_data_start_date*100000000)/100000000
         if (FLUXNET_skip_years == 0){
            if ((signif(current.date, digits = 4)/100000000)%%4 == 0){
               if (from.yyyymmdd(signif(current.date, digits = 8)/10000) > from.yyyymmdd(signif(current.date, digits = 4)/10000 + 229)) num_leap=1
            } else num_leap = 0 
         } else num_leap = num.elem((FLUXNET_data_start_date : (FLUXNET_data_start_date + FLUXNET_skip_years - 1))%%4 == 0)
         n_day_of_the_year = (from.yyyymmdd(signif(current.date, digits = 8)/10000) - from.yyyymmdd(signif(current.date, digits = 4)/10000 + 101))/24
         if (FLUXNET.nrows < 48) rows_per_day = 24 else rows_per_day = 48
         # Inclusive of the day before or after
         if (is.na(rangeloc)) extra = 0
         else if (rangeloc == 'front') extra = -1
         else if (rangeloc == 'back') extra = 1
         return((FLUXNET_skip_years*365 + num_leap + n_day_of_the_year + extra)*rows_per_day)
      } else if (direction == 'backward') {
         # This section is not complete nor required...
         # FLUXNET_skip_years = (FLUXNET_data_end_date*100000000 - signif(current.date, digits = 4))/100000000
         # if (FLUXNET_skip_years == 0){
         #    if ((signif(current.date, digits = 4)/100000000)%%4 == 0){
         #       if (from.yyyymmdd(signif(current.date, digits = 8)/10000) > from.yyyymmdd(signif(current.date, digits = 4)/10000 + 229)) num_leap=1
         #    } else num_leap = 0 
         # } else num_leap = num.elem(((FLUXNET_data_end_date - FLUXNET_skip_years + 1) : FLUXNET_data_end_date)%%4 == 0)
         # n_day_of_the_year = (from.yyyymmdd(signif(current.date, digits = 4)/10000 + 1231) - from.yyyymmdd(signif(current.date, digits = 8)/10000))/24
         # return((FLUXNET_skip_years*365 + num_leap + n_day_of_the_year)*rows_per_day)
      }
   } else if (!hourly.data){
      if (direction == 'forward') {
         FLUXNET_skip_years = (signif(current.date, digits = 4) - FLUXNET_data_start_date*10000)/10000
         n_months = as.numeric(substr(current.date, start = 5, stop = 6)) - 1
      } else if (direction == 'backward') {
         FLUXNET_skip_years = (FLUXNET_data_end_date*10000 - signif(current.date, digits = 4))/10000
         n_months = as.numeric(substr(current.date, start = 5, stop = 6)) - 1
      }
      return(FLUXNET_skip_years*12 + n_months)
   }
}

################################################################################
### Function to check if FLUXNET data of the site chosen has listed input variables data
f_FLUXNET_variable_check = function(FLUXNET.dir, FLUXNET.header, var.match.df, site.id){
   # Get current FLUXNET data variables:
   data_variable = FLUXNET.header
   err_output = vector()
   counter = 1
   
   # Check FLUXNET variable availability:
   for (irow in seq_along(var.match.df[, "FLUXNET_var_name"])){
      if ((!is.na(var.match.df[irow, "FLUXNET_var_name"])) && (!var.match.df[irow, "FLUXNET_var_name"] %in% data_variable)){
         line1 = paste0('FLUXNET data of site ', site.id, ' does not have variable ', var.match.df[irow ,"FLUXNET_var_name"], '.')
         line2 = paste0('Error resolved by using variable ', var.match.df[irow, "TEMIR_var_name"], ' of ', met_name, ' data as substitute.')
         print((line1), quote = FALSE)
         print(paste0('    ', line2), quote = FALSE)
         err_output[counter] = paste(line1, line2, sep = " ")
         var.match.df[irow, "FLUXNET_var_name"] = NA
         counter = counter + 1
      }
   }
   
   # Print check completion message
   if (length(err_output) > 0){
      print('All other variables entered exist in FLUXNET data of the concerned site.', quote = FALSE)
   }
   return(list(variable_df = var.match.df, global_err = err_output))
}

################################################################################
### Function to check if FLUXNET data of the site chosen has listed input variables data
f_FLUXNET_HH_to_HR = function(FLUXNET.sel.data, method){
   if (method == 'mean') {
      return(unname(tapply(FLUXNET.sel.data, INDEX = rep(1:24, each = 2), FUN = mean)))
   } else if (method == 'sum') {
      return(unname(tapply(FLUXNET.sel.data, INDEX = rep(1:24, each = 2), FUN = sum)))
   } else if (method == 'exclude') {
      return(FLUXNET.sel.data[c(TRUE, FALSE)])
   }
}

################################################################################
### Function to convert relative humidity to specific humidity
f_RH_to_SH = function(data.array, ind.lon, ind.lat, M.w = M_w, M.d = M_da, R = R_uni, T.atm, P.atm) {
   
   # Set lon lat for meteorological fields:
   T_K = T.atm[ind_lon, ind_lat, ]
   P = P.atm[ind_lon, ind_lat, ]
   
   # Calculation loop for each hour:
   for (ihr in 1:dim(data.array)[3]){
      # f_esat saturation vapour pressure (Pa) is found in Farquhar_Ball_Berry.R
      e_sat = f_esat(T_K = T_K[ihr])
      # Saturation volume mixing ratio w_s:
      w_s =  (M.w/M.d) * e_sat/ (P[ihr] - e_sat)
      # Volume mixing ratio w from relative humidity:
      w = data.array[ind_lon, ind_lat, ihr] * w_s / 100
      # Partial vapour pressure (Pa) e:
      e = w * P[ihr] / ((M.w/M.d) + w)
      # Specific humidity calculated from water vapour mass/moist air mass with ideal gas law:
      data.array[ind_lon, ind_lat, ihr] = (M.w/M.d) * e / (P[ihr] - (1 - (M.w/M.d)) * e)
   }
   
   return(data.array)
}

################################################################################
### Function from pressure to sea level pressure in Pa (Barometric law with standard temperature lapse rate) 
# NOTE : NOT NEEDED IN v1.1
if (FALSE) {
   f_PA_to_SLP = function(FLUXNET.dir, data.array, M.da = M_da, g = g_E, R = R_uni, T_K = T2M[ind_lon, ind_lat, ], ind.lon = ind_lon, ind.lat = ind_lat){
      site_info = if (file.exists(paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx'))) read_excel(path = paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx')) else read_excel(path = paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx'))
      site_ind = match(FLUXNET_site_id, site_info$SITE_ID)
      Z_surf = as.double(unname(site_info[site_ind,"LOCATION_ELEV"]))
      rm(site_info, site_ind)
      # L is standard temperature lapse rate (K/m) in ISA
      if (11000 <= Z_surf && Z_surf < 20000) L = 0.0
      else if (Z_surf < 11000) L = -0.0065
      pwr = g*M.da/(R*L)
      for (ihr in 1:dim(data.array)[3]){
         fctr = T_K[ihr]/(T_K[ihr] - Z_surf * L)
         data.array[ind.lon, ind.lat, ihr] =  data.array[ind.lon, ind.lat, ihr] * fctr ^ pwr
      }
      return(data.array)
   }
}

################################################################################
### Function to write FLUXNET data into the site grid of meteorological data
f_FLUXNET_met_grid = function(FLUXNET.data, TEMIR.variable, FLUXNET.variable, FLUXNET.nrows, dt.hr, st.ind, end.ind, ind.lon, ind.lat, input.var.df, hh.to.hr.flag){
   
   # Find variable index:
   index = match(TEMIR.variable, input.var.df[, "TEMIR_var_name"])
   
   # Catch variable error:
   if (is.na(index)) stop('Invalid meteorlogical variable given!!!')
   err_out = NULL
   
   # Get meterological data and field:
   selected_var_data = f_FLUXNET_unit_convert(TEMIR.variable = TEMIR.variable, FLUXNET.var = FLUXNET.variable, FLUXNET.data = FLUXNET.data)
   data_array = get(x = TEMIR.variable)
   
   # Convert FLUXNET data from half-hourly to hourly data accordingly:
   if (hh.to.hr.flag) {
      HR_sum_var = c('P_F')
      HR_mean_var = c()
      if (any(FLUXNET.variable == HR_sum_var)) {
         HH_HR_method = 'sum' 
      } else if (any(FLUXNET.variable == HR_sum_var)) {
         HH_HR_method = 'mean'
      } else HH_HR_method = 'exclude'
      data_select = f_FLUXNET_HH_to_HR(FLUXNET.sel.data = selected_var_data, method = HH_HR_method)
   } else data_select = selected_var_data
   
   # Format FLUXNET data for current simulation day for abnormal time structure:
   if (end.ind - st.ind == 24/dt.hr) {
      # Remove extra hour on current simulation day:
      if (num.elem(FLUXNET.data[, FLUXNET.variable] == -9999) == 1/dt.hr || table(duplicated(data_select))['TRUE'] == 1/dt.hr) {
         data_select = data_select[data_select != -9999]
         if (length(data_select) == (24+1)/dt.hr) data_select = data_select[!duplicated(data_select)]
         data_array[ind_lon, ind_lat,] = data_select
      } else stop('FLUXNET data cannot be converted into grid!')
   } else if (end.ind - st.ind == 22/dt.hr) {
      # Fill in missing hour on current simulation day:
      data_array[ind_lon,ind_lat,] = c(data_select, tail(x = data_select, n = 1/dt.hr))
   } else if (!is.na(match(-9999, FLUXNET.data[, FLUXNET.variable])) || NA %in% FLUXNET.data[, FLUXNET.variable]) {
      # Default to global meteorology as FLUXNET data has unknown time structure for current simulation day:
      line1 = paste0('FLUXNET data ', FLUXNET.variable, ' is incomplete or unreliable for day ', current_date,'!')
      line2 = paste0('Using ', TEMIR.variable, ' of ', met_name, ' data as replacement....')
      print(line1, quote = FALSE)
      print(paste0('    ',line2), quote = FALSE)
      err_out = c('Variable Error' = FLUXNET.variable, 'Variable Replacement' = TEMIR.variable)
      
   } else {
      # NOTE : Conversion to compatible data was done in the function f_FLUXNET_unit_convert
      # should have no problem to directly write into array
      data_array[ind_lon, ind_lat,] = data_select
   }
   
   # Conversion calculations from FLUXNET variables to TEMIR variables
   # Converting relative humidity to specific humidity
   if (TEMIR.variable == 'QV2M') data_array = f_RH_to_SH(data.array = data_array, ind.lon = ind.lon, ind.lat = ind.lat, M.w = M_w, M.d = M_da, R = R_uni, T.atm = T2M, P.atm = ATMP)
   # Converting pressure into sea-level pressure # Not required for FLUXNET with constant pressure profile 
   # if (TEMIR.variable == 'SLP') data_array = f_PA_to_SLP(data.array = data_array)
   
   return(list(data = data_array, err_out = err_out))
}

################################################################################
### Function for variable agreement between FLUXNET and meteorological data
f_FLUXNET_unit_convert = function(TEMIR.variable, FLUXNET.variable, FLUXNET.data){
   # Get particular meteorological data:
   data_return = FLUXNET.data[, FLUXNET.variable]
   # Degree to Kelvin:
   if (any(TEMIR.variable == c('T2M', 'T10M'))) data_return = FLUXNET.data[, FLUXNET.variable] + 273.15
   # kPa to Pa:
   if (TEMIR.variable == 'ATMP') data_return = FLUXNET.data[, FLUXNET.variable] * 1000
   # umol m^-2 s^-1 to W m^-2:
   if (any(TEMIR.variable == c('PARDF', 'PARDR'))) data_return = FLUXNET.data[, FLUXNET.variable] * 0.219
   # mm h-1 to kg m^2 s-1:
   if (TEMIR.variable == 'PRECTOT') data_return = FLUXNET.data[, FLUXNET.variable] / 3600
   # Volumetric soil moisture to soil wetness: done in simulate_IJ (line 365)
   return(data_return)
}

################################################################################
### End of file
################################################################################
