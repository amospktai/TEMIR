################################################################################
### Terrestrial Ecosystem Model in R (TEMIR)
### Execution module for single-site, regional or global simulation
################################################################################

# Load necessary libraries upfront:
library(tools) # For utilities
library(methods) # For basic R utility
library(dotCall64) # For interfacing other languages
library(stringr) # For string manipulations
library(abind) # For array manipulations
library(parallel) # For parallel computing
library(ncdf4) # For reading ncdf4 files
library(filesstrings) # For file and directory manipulations
# Load libraries for post-simulation processing:
library(grid) # For graphics
library(spam) # For matrix manipulations
library(maps) # For map display
library(fields) # For spatial data
timestamp()

################################################################################
### Source input scripts for model configuration:
################################################################################

options(warn = 2)

# Source input scripts:
# *** Please make sure the input scripts (e.g., input_TEMIR_basic_settings.R) are in the same simulation directory as this execution script. ***
source('input_TEMIR_basic_settings.R')

if (file.exists('input_TEMIR_biogeochem_extension.R')) {
   print("Sourcing 'input_TEMIR_biogeochem_extension.R'")
    source('input_TEMIR_biogeochem_extension.R')
} else {
    biogeochem_flag = FALSE
}

print(O3_data_dir)

# Check existence of directory paths:
dir_check = ls(pattern = "_dir$")
NA_dirs = match(optional_dirs, dir_check)[is.na(sapply(X = optional_dirs, FUN = get))]
if (length(NA_dirs) != 0) dir_check = dir_check[-NA_dirs]
for (idir in seq_along(dir_check)) {
   # if (length(Sys.glob(paths = get(dir_check[idir]))) == 0) stop(paste0(dir_check[idir], ' does not exist!'))
   if (!dir.exists(paths = get(dir_check[idir]))) stop(paste0(dir_check[idir], ' does not exist!'))
   # assign(dir_check[idir], paste0(file_path_as_absolute(Sys.glob(get(dir_check[idir]))), '/'))
}
rm(dir_check, idir, optional_dirs, NA_dirs)

# Check whether the current working directory is a simulation directory:
simulation_name.txt = paste0(basename(getwd()), '.txt')
if (!file.exists(simulation_name.txt)) stop('Current working directory is not a simulation directory containing the execution_v*.R script! Please set the working directory to be the simulation directory to run the model.')
rm(simulation_name.txt)

# Set simulation directory:
simulation_dir = paste0(getwd(), '/')
print(paste0('Simulation directory: ', simulation_dir), quote=FALSE)
print('Please make sure this is the simulation that you desire!', quote=FALSE)
cat('\n')

# Current TEMIR execution and input_TEMIR directory path:
# simulation_dir = paste0(file_path_as_absolute('.'), '/')
# if (dirname(simulation_dir) != file_path_as_absolute(TEMIR_dir)) simulation_dir = file_path_as_absolute(dirname(sub('--file=', '', commandArgs(trailingOnly = FALSE)[grep('--file=', commandArgs(trailingOnly = FALSE))])))

# Check TEMIR directory path:
# if (dirname(simulation_dir) != file_path_as_absolute(TEMIR_dir)) {
   # print(paste0('Simulation directory : ', simulation_dir))
   # print(paste0('TEMIR directory : ', TEMIR_dir))
   # stop('Simulation directory is not a subdirectory in the TEMIR directory!')
# }
# setwd(simulation_dir)

# Get simulation configuration:
model_config_vec = setdiff(ls(), ls(pattern = '_dir$'))

# Turn on default TEMIR aerodynamic conductance scheme if Monin_Obukhov_flag=TRUE:
if (Monin_Obukhov_flag) use_TEMIR_ga_flag = TRUE


################################################################################
### Source scripts
################################################################################

# Globally available geophysical constants:
# All are in SI units when applicable.
source(paste0(code_dir, 'geophys_const.R'))

# Useful functions:
source(paste0(code_dir, 'tools.R'))

# Functions to compute aerodynamic conductance, temperature and humidity profiles:
source(paste0(code_dir, 'Monin_Obukhov.R'))

# Functions to compute canopy radiative transfer:
source(paste0(code_dir, 'radiative_transfer.R'))

# Functions to compute photosynthesis:
source(paste0(code_dir, 'Farquhar_Ball_Berry.R'))

# Functions to loop simulation over longitudes and latitudes:
# It includes a function to reshape history data into gridded format.
source(paste0(code_dir, 'simulate_ij.R'))

# Functions to compute dry deposition velocities and required parameters:
source(paste0(code_dir, 'drydep_toolbox.R'))

# Functions to incorporate FLUXNET datasets:
if (!is.na(FLUXNET_dir)) source(paste0(code_dir, 'FLUXNET_functions.R'))

# Functions to incorporate meteorological data:
# Obsolete, but may revive for future purposes (Tai, Jan 2019)
# source(paste0(code_dir, 'met_tools.R'))

# Surface and PFT input parameters (prescribed):
source(paste0(code_dir, 'PFT_surf_data.R'))

# Functions used in biogeochemistry module
if(biogeochem_flag){
   source(paste0(BGC_code_dir, 'biomass_partitioning.R'))
   source(paste0(BGC_code_dir, 'plant_phenology.R'))
   source(paste0(BGC_code_dir, 'plant_physiology.R'))
   source(paste0(BGC_code_dir, 'maintenance_respiration.R'))
   source(paste0(BGC_code_dir, 'Cpool_addition_subtraction.R'))
}

################################################################################
### Define additional dimensional and variable info:
################################################################################

cat('\n')

# Time step (s):
dt = dt_hr*3600

# Vector of simulation dates:
date_vec = make.date.vec(start.date=start_date, end.date=end_date)

# Number of simulation days:
n_day_sim = length(date_vec)

# Global simulation setting:
if (!single_site_flag) {
   # Set FLUXNET_flag as FALSE as single_site_flag is FALSE 
   FLUXNET_site_flag = FALSE
   FLUXNET_flag = FALSE
}

# FLUXNET site setting:
if (FLUXNET_site_flag){
   if (FLUXNET_flag){
      # Check for FLUXNET data availability over the simulation days:
      FLUXNET_time_info = f_FLUXNET_UTC_2_local(FLUXNET.dir = FLUXNET_dir, date = start_date, site.id = FLUXNET_site_id, utc.offset = NA, out.utc.offset = TRUE)
      time_shift = FLUXNET_time_info$time_diff
      shifted_start = FLUXNET_time_info$time ; rm(FLUXNET_time_info)
      shifted_end = f_FLUXNET_UTC_2_local(FLUXNET.dir = FLUXNET_dir, date = to.yyyymmdd(from.yyyymmdd(end_date) + 24), site.id = FLUXNET_site_id, utc.offset = time_shift, out.utc.offset = FALSE)
      day_mod_flag = FALSE
      # Shift start date by one day to see if data is in range
      if (!f_FLUXNET_date_range(FLUXNET.dir = FLUXNET_dir, start.or.end.year = 'start', whole.date.for.year = shifted_start, site.id = FLUXNET_site_id)) {
         temp_date = start_date
         start_date = to.yyyymmdd(from.yyyymmdd(start_date) + 24)
         print(paste0('Time shifted FLUXNET data is not available for start_date!'), quote = FALSE) 
         print(paste0('So delaying start_date ', temp_date,' to the next day ', start_date), quote = FALSE); remove(temp_date)
         day_mod_flag = TRUE
      }
      # Shift end date by one day to see if data is in range
      if (!f_FLUXNET_date_range(FLUXNET.dir = FLUXNET_dir, start.or.end.year = 'end', whole.date.for.year = shifted_end, site.id = FLUXNET_site_id)) {
         temp_date = end_date
         end_date = to.yyyymmdd(from.yyyymmdd(end_date) - 24)
         print(paste0('Time shifted FLUXNET data is not available for end_date!'), quote = FALSE)
         print(paste0('So advancing end_date ', temp_date,' to the previous day ', end_date), quote = FALSE); remove(temp_date)
         day_mod_flag = TRUE
      }
      # Check for FLUXNET data availability over the shifted simulation days:
      if (day_mod_flag) {
         FLUXNET_time_info = f_FLUXNET_UTC_2_local(FLUXNET.dir = FLUXNET_dir, date = start_date, site.id = FLUXNET_site_id, utc.offset = NA, out.utc.offset = TRUE)
         time_shift = FLUXNET_time_info$time_diff
         shifted_start = FLUXNET_time_info$time ; rm(FLUXNET_time_info)
         shifted_end = f_FLUXNET_UTC_2_local(FLUXNET.dir = FLUXNET_dir, date = to.yyyymmdd(from.yyyymmdd(end_date) + 24), site.id = FLUXNET_site_id, utc.offset = time_shift, out.utc.offset = FALSE)
         # Check if FLUXNET data is available for shifted dates
         if (!f_FLUXNET_date_range(FLUXNET.dir = FLUXNET_dir, start.or.end.year = 'start', whole.date.for.year = shifted_start, site.id = FLUXNET_site_id)) stop(paste0('Time shifted FLUXNET data is not available for start_date even after delaying by a day!!! Please check FLUXNET data of site ', FLUXNET_site_id, '....'))
         if (!f_FLUXNET_date_range(FLUXNET.dir = FLUXNET_dir, start.or.end.year = 'end', whole.date.for.year = shifted_end, site.id = FLUXNET_site_id)) stop(paste0('Time shifted FLUXNET data is not available for end_date even after advancing by a day!!! Please check FLUXNET data of site ', FLUXNET_site_id, '....'))
      }
      # FLUXNET data taken as better data, replacing Monin-Obukhov outputs:
      Monin_Obukhov_flag = FALSE
   }
   # Get location (lon, lat) from FLUXNET_site_id:
   FLUXNET_lon_lat = f_lon_lat_sim_from_FLUXNET(FLUXNET.dir = FLUXNET_dir, site.id = FLUXNET_site_id)
   lon_sim = FLUXNET_lon_lat$lon_sim
   lat_sim = FLUXNET_lon_lat$lat_sim
   rm(FLUXNET_lon_lat)
}

# Corresponding indices in global lon/lat grid for the simulated region:
if (single_site_flag) {
   ind_lon = find.lon.lat(lonspec=lon_sim, latspec=lat_sim, lon=lon, lat=lat)[1]
   ind_lat = find.lon.lat(lonspec=lon_sim, latspec=lat_sim, lon=lon, lat=lat)[2]
} else {
   ind_lon = which((lon + dlon/2) > lon_sim[1] & (lon - dlon/2) < lon_sim[2])
   ind_lat = which((lat + dlat/2) > lat_sim[1] & (lat - dlat/2) < lat_sim[2])
}

# Define a list of c(i, j), where i and j are indices of global lon/lat grid to be simulated:
ij = list()
n = 0
for (n_i in 1:length(ind_lon)) {
   i = ind_lon[n_i]
   for (n_j in 1:length(ind_lat)) {
      j = ind_lat[n_j]
      if (single_site_flag) {
         # Runs site simulation unless land fraction = 0
         if (FRLAND[i,j] != 0.0) {
            n = n + 1 ; ij[[n]] = c(i, j)
            if (!FLUXNET_site_flag) {
               print(paste0('Single site simulation using ', met_name, ' data :'), quote = FALSE)
               print(paste0('Longitude = ', lon_sim, ', Latitude = ', lat_sim), quote = FALSE)
               print(paste0('Land cover fraction of site = ', FRLAND[i,j]), quote = FALSE)
            } else { 
               if (FLUXNET_flag){
                  print(paste0('Simulation using FLUXNET data of site ', FLUXNET_site_id, ' (gapfilled with ', met_name, '):'), quote = FALSE)
               } else {
                  print(paste0('Simulation using ', met_name, ' data of FLUXNET site ', FLUXNET_site_id, ' :'), quote = FALSE)
               }
               print(paste0('Longitude = ', lon_sim, ', Latitude = ', lat_sim), quote = FALSE)
               print(paste0('Land cover fraction of site ', FLUXNET_site_id, ' = ', FRLAND[i,j]), quote = FALSE)
            }
         } else stop('Grid cell contains no land!!!')
      } else { 
         if (FRLAND[i,j] < 0.05) { 
            # Too little non-glacial land: skip calculations altogether and do not include in list of c(i, j). This prevents calculations for oceanic and permanently glacial grid cells.
            n = n
         } else {
            # n = (n_i - 1)*length(ind_lat) + n_j
            n = n + 1
            ij[[n]] = c(i, j)
         }
      }   
   }
}

# Print general simulation information:
if (!single_site_flag) {
   if (lon_sim == c(-180, 180) && lat_sim == c(-90, 90)){
      print(paste0('Global simulation using ', met_name, ' data...'), quote = FALSE)
   } else {
      print(paste0('Simulation using ', met_name, ' data : longitude = ', lon_sim[1],' to ', lon_sim[2], ', latitude = ', lat_sim[1], ' to ', lat_sim[2] ,'...'), quote = FALSE)
   }
}

################################################################################

# Define nc data dimensions and variables if nc archiving is turned on:

# Define variable dimensions:
X = ncdim_def(name='lon', units='degrees_east', vals=lon[ind_lon], longname='Longitude')
Y = ncdim_def(name='lat', units='degrees_north', vals=lat[ind_lat], longname='Latitude')
# D = ncdim_def(name='date', units='YYYYMMDD', vals=date_vec, unlim=TRUE, calendar='standard', longname='Date in YYYYMMDD (defined by 00:00-24:00 UTC)')
P = ncdim_def(name='pft', units='unitless', vals=pftnum, longname='Plant function type (0 = bare land)')
H = ncdim_def(name='hour', units='hours', vals=dt_hr*(1:(24/dt_hr)), longname='Hour on YYYYMMDD (starting from 00:00 UTC)')

# Make list of variables selected to archive:
# If each nc file contains data for a day of outputs, "dim=list(X, Y, P, H)" in "ncvar_def". If it contains all days, "dim=list(X, Y, D, P, H)".
var_name = available_outputs_df[na.omit(match(unique(output_variables), available_outputs_df$variable_name)),]
var_list = list()
var_name[] = lapply(var_name, as.character) # for R version 3.1.1 
for (v in 1:nrow(var_name)) {
   var_dim_list = if (var_name[v,'res_level'] == 'PFT') {list(X, Y, P, H)} else if (var_name[v,'res_level'] == 'grid') {list(X, Y, H)} else if (var_name[v,'res_level'] == 'PFT_daily') {list(X, Y, P)}
   var_list[[v]] = ncvar_def(name=var_name[v,'variable_name'], units=var_name[v,'unit'], dim=var_dim_list, longname=var_name[v,'long_name'], prec='float', compression=4)
}
rm(var_dim_list)

cat('\n')

################################################################################
### Simulation:
################################################################################

print(paste0('### Begin simulation for case: ', basename(simulation_dir), ' ###'), quote=FALSE)

# Number of cores to use:
# If running on a cluster instead of a personal computer, overwrite what is set above with these system-dependent settings:
if (cluster_flag) {
   # These settings are dependent on the cluster environment.
   n_core = as.numeric(Sys.getenv('PBS_NUM_PPN'))
   if (is.na(n_core)) n_core = detectCores()
   n_node = as.numeric(Sys.getenv('PBS_NUM_NODES'))
   if (is.na(n_node)) n_node = 1
   if (n_node > 1) stop('CLM-R cannot be run on multiple nodes!')
}
if (multicore_flag) print(paste0('# of cores used: ', as.character(n_core)), quote=FALSE) else print('# of cores used: 1', quote=FALSE)

# Continue run:
if (continue_flag) print('Continuing simulation from previous run...', quote=FALSE)

# Check and reshape data arrays:
if (debug_flag) {
   # Data array holding all days of data will be created:
   hist_grid = array(NaN, dim=c(length(ind_lon), length(ind_lat), n_day_sim, length(pftname), 24/dt_hr, nrow(var_name)))
   # Check error for all days:
   err_hist_ij = list()
   err_msg = NULL
}

# Input (meteorological and vegetation) variables dataframe for simulation from various sources (NA if no corresponding input exists):
# 1st col = TEMIR input variable name; 2nd col = TEMIR input unit; 
# 3rd col = MERRA2 variable name; 4th col = MERRA2 unit conversion needed;
# 5th col = GEOS-FP variable name; 6th col = GEOS-FP unit conversion needed;
# 7th col = FLUXNET variable name; 8th col = FLUXNET conversion needed
input_data_df = `colnames<-`(rbind.data.frame(
   # Cloud fraction (0-1)
   c('CLDTOT', '0-1', 'CLDTOT', FALSE, 'CLDTOT', FALSE, NA, FALSE),
   # Land/water/ice flags (0-1)
   # c('LWI', '0-1', 'LWI', FALSE, 'LWI', FALSE, NA, FALSE),
   # Root zone soil wetness (0-1)
   c('GWETROOT', '0-1', 'GWETROOT', FALSE, 'GWETROOT', FALSE, NA, FALSE),
   # Top soil wetness (0-1)
   c('GWETTOP','0-1', 'GWETTOP', FALSE, 'GWETTOP', FALSE, NA, FALSE),
   # Sea-level pressure (Pa)
   c('SLP', 'Pa', 'SLP', FALSE, 'SLP', TRUE, NA, FALSE),
   # Atmospheric pressure (Pa)
   c('ATMP', 'Pa', NA, FALSE, NA, FALSE, 'PA_F', TRUE), 
   # Temperature at 10 m above displacement height (K)
   c('T10M', 'K', 'T10M', FALSE, 'T10M', FALSE, 'TA_F', FALSE),
   # Temperature at 2 m above displacement height (K)
   c('T2M', 'K', 'T2M', FALSE, 'T2M', FALSE, 'TA_F', FALSE),
   # Sensible heat flux (W m^-2)
   c('HFLUX', 'W m^-2','HFLUX', FALSE,'HFLUX', FALSE, 'H_F_MDS', FALSE),
   # Latent heat flux (W m^-2)
   c('EFLUX', 'W m^-2', 'EFLUX', FALSE, 'EFLUX', FALSE, 'LE_F_MDS', FALSE),
   # Surface downward PAR diffuse flux (W m^-2)
   c('PARDF', 'W m^-2', 'PARDF', FALSE,'PARDF', FALSE, 'PPFD_DIF', TRUE),
   # Surface downward PAR direct beam flux (W m^-2)
   c('PARDR', 'W m^-2', 'PARDR', FALSE, 'PARDR', FALSE, 'PPFD_IN', TRUE),
   # Surface incident shortwave flux (W m^-2)
   c('SWGDN', 'W m^-2', 'SWGDN', FALSE, 'SWGDN', FALSE, 'SW_IN_F', FALSE),
   # Total Precipitation (kg m-2 s-1)
   c('PRECTOT', 'kg m-2 s-1', 'PRECTOT', FALSE, 'PRECTOT', FALSE, 'P_F', TRUE),
   # Snowfall (kg m-2 s-1)
   c('PRECSNO', 'kg m-2 s-1', 'PRECSNO', FALSE, 'PRECSNO', FALSE, NA, FALSE),
   # Snow depth (m)
   c('SNODP', 'm', 'SNODP', FALSE, 'SNODP', FALSE, NA, FALSE),
   # Surface evaporation (kg m^-2 s^-1)
   c('EVAP', 'kg m^-2 s^-1', 'EVAP', FALSE, 'EVAP', FALSE, NA, FALSE),
   # Specific humidity at 2 m above displacement height (kg kg^-1)
   c('QV2M', 'kg kg^-1', 'QV2M', FALSE, 'QV2M', FALSE, 'RH', TRUE),
   # Specific humidity at 10 m above displacement height (kg kg^-1) # Not included in GEOS-Chem met fields
   # c('QV10M', 'kg kg^-1', 'QV10M', FALSE, 'QV10M', FALSE, NA, FALSE),
   # Eastward wind at 10 m above displacement height (m s^-1)
   c('U10M', 'm s^-1', 'U10M', FALSE, 'U10M', FALSE, NA, FALSE),
   # Northward wind at 10 m above displacement height (m s^-1)
   c('V10M','m s^-1', 'V10M', FALSE,'V10M', FALSE, NA, FALSE),
   # Wind speed (m s^-1)
   c('WS','m s^-1', NA, FALSE, NA, FALSE, 'WS_F', FALSE),
   # Friction velocity (m s^-1)
   c('USTAR', 'm s^-1', 'USTAR', FALSE, 'USTAR', FALSE, 'USTAR', FALSE),
   # Roughness length for momentum (m)
   c('Z0M', 'm', 'Z0M', FALSE, 'Z0M', FALSE, NA, FALSE),
   # Dataframe settings
   stringsAsFactors = FALSE),
   c('TEMIR_var_name', 'TEMIR_unit', 'MERRA2_var_name', 'MERRA2_unit_differ', 
     'GEOSFP_var_name', 'GEOSFP_unit_differ', 'FLUXNET_var_name', 'FLUXNET_unit_differ'))

# Set FLUXNET data directory and information if FLUXNET_flag=TRUE:
if (FLUXNET_flag) {
   
   # Get FLUXNET data information:
   FLUXNET_file_settings = f_FLUXNET_file(FLUXNET.dir = FLUXNET_dir, site.id = FLUXNET_site_id, hr.part = dt_hr, hourly.data = TRUE)
   FLUXNET_file = FLUXNET_file_settings$filedir
   FLUXNET_nrows = as.numeric(FLUXNET_file_settings$nrows) ; remove(FLUXNET_file_settings)
   FLUXNET_header = names(read.csv(file = FLUXNET_file, nrows = 1, header = TRUE, stringsAsFactors = FALSE))
   FLUXNET_check = f_FLUXNET_variable_check(FLUXNET.dir = FLUXNET_dir, FLUXNET.header = FLUXNET_header, site.id = FLUXNET_site_id, var.match.df = input_data_df)
   updated_input_df = FLUXNET_check$variable_df
   FLUXNET_global_err = FLUXNET_check$global_err ; remove(FLUXNET_check)
   
   # Get starting row of FLUXNET data for start date:
   data_start = f_FLUXNET_row_skip(FLUXNET.dir = FLUXNET_dir, site.id = FLUXNET_site_id, current.date = shifted_start, direction = 'forward', rangeloc = 'front', FLUXNET.nrows = FLUXNET_nrows, hourly.data = TRUE)
   
   # Get subset of FLUXNET data for relevant dates:
   FLUXNET_selected_data = read.csv(file = FLUXNET_file, skip = data_start, nrows = f_FLUXNET_row_skip(FLUXNET.dir = FLUXNET_dir, site.id = FLUXNET_site_id, current.date = shifted_end, direction = 'forward', rangeloc = 'back', FLUXNET.nrows = FLUXNET_nrows, hourly.data = TRUE) - data_start, header=TRUE)
   colnames(FLUXNET_selected_data) = FLUXNET_header
   FLUXNET_selected_data = FLUXNET_selected_data[,c('TIMESTAMP_START', 'TIMESTAMP_END', unique(na.omit(updated_input_df$FLUXNET_var_name)))]
   
   # Subset half-hourly data into hourly data if required
   if (dt_hr == 1 && substr(basename(FLUXNET_file), 32, 33) == 'HH') {
      print('NOTE : Half-hourly data is subsetted into hourly data BUT precipitation is coverted accordingly')
      HH_to_HR_flag = TRUE
   } else HH_to_HR_flag = FALSE
   
} else updated_input_df = input_data_df

################################################################################

# Start simulation for each day:

for (d in 1:n_day_sim) {
   
   timestamp()
   
   # Current date:
   current_date = to.yyyymmdd(from.yyyymmdd(start_date) + (d - 1)*24)
   # Strings for year, month and day:
   YYYY = substr(x=as.character(current_date), start=1, stop=4)
   MM = substr(x=as.character(current_date), start=5, stop=6)
   DD = substr(x=as.character(current_date), start=7, stop=8)
   print(paste0('Current simulation date = ', YYYY, '/', MM, '/', DD), quote=FALSE)
   
   # Create temporary data directory for each simulation day:
   # if (length(dir(path=paste0(simulation_dir, 'temp_data/'), pattern=paste0('temp_', YYYY, MM, DD))) == 0) system(command=paste0("mkdir '", simulation_dir, 'temp_data/temp_', YYYY, MM, DD, "'"))
   if (!dir.exists(paste0(simulation_dir, 'temp_data/temp_', YYYY, MM, DD))) dir.create(paste0(simulation_dir, 'temp_data/temp_', YYYY, MM, DD))
   
   # Number of days from 00:00 UTC Jan 1 (whole number not including hours):
   leap = is.leap(yyyy = as.numeric(YYYY))
   n_day_whole = date.to.day(yyyymmdd = current_date, leap = leap)
   
   # Time shifting FLUXNET to UTC for current simulation date:
   if (FLUXNET_flag){
      FLUXNET_current = f_FLUXNET_UTC_2_local(FLUXNET.dir = FLUXNET_dir, date = current_date, site.id = FLUXNET_site_id, utc.offset = time_shift, out.utc.offset = FALSE)
      FLUXNET_end = f_FLUXNET_UTC_2_local(FLUXNET.dir = FLUXNET_dir, date = to.yyyymmdd(from.yyyymmdd(current_date) + 24), site.id = FLUXNET_site_id, utc.offset = time_shift, out.utc.offset = FALSE)
      # Time shift obtained daily as daylight saving time is not required as Fluxnet discounts daylight saving time
   }
   
   #############################################################################
   
   # Reload interannually varying inputs on YYYY/01/01 if necessary:
   
   if (paste0(MM, DD) == '0101') {
      
      # Reload LAI data if interannually varying LAI data are used:
      if (LAI_data_flag) {
         if (length(year_vec) > 1) {
            YYYY_LAI = as.character(LAI_used_years[which(year_vec == as.numeric(YYYY))])
            last_YYYY_LAI = as.character(LAI_used_years[which(year_vec == as.numeric(YYYY)) - 1])
            if (YYYY_LAI != last_YYYY_LAI) {
               subfn = paste0('daily_monthly_LAI_data_', YYYY_LAI, '.RData')
               filename = paste0(processed_surf_data_dir, subfn)
               print(paste0('Loading existing daily LAI and SAI data for year ', YYYY_LAI, '...'), quote=FALSE)
               load(filename)
            }
         }
      }
      
      # Reload hourly ozone field:
      # soyFACE temporary 
      if (O3_damage_flag & !O3_fixed_flag) {
         # if (length(year_vec) > 1) {
         #    YYYY_O3 = as.character(O3_used_years[which(year_vec == as.numeric(YYYY))])
         #    last_YYYY_O3 = as.character(O3_used_years[which(year_vec == as.numeric(YYYY)) - 1])
         #    if (YYYY_O3 != last_YYYY_O3) {
         #       filename = paste0(O3_data_dir, O3_subn1, YYYY_O3, O3_subn2)
         #       print(paste0('Loading surface O3 concentrations from ', filename, '...'), quote=FALSE)
         #       nc = nc_open(filename)
         #       lon_O3 = ncvar_get(nc, unname(O3_dim_vec['longitude']))
         #       lat_O3 = ncvar_get(nc, unname(O3_dim_vec['latitude']))
         #       # Surface O3 concentration (ppbv):
         #       O3_hourly = ncvar_get(nc, O3_array_name)
         #       nc_close(nc)
         #       # Regrid to model resolution if input resolution is not consistent:
         #       if (sum(lon != lon_O3) > 0 | sum(lat[2:(length(lat)-1)] != lat_O3[2:(length(lat_O3)-1)]) > 0) {
         #          # Regrid to model resolution:
         #          print('Regridding hourly O3 concentrations for year ', YYYY_O3, '...', quote=FALSE)
         #          O3_hourly = sp.regrid(spdata=O3_hourly, lon.in=lon_O3, lat.in=lat_O3, lon.out=lon, lat.out=lat)
         #       }
         #    }
         # }
         
         # soyFACE temporary 
         if (force_prescribed_o3) {
            O3_hourly = ozone_soyFACE_vec
         }

      }
      
   }
   
   #############################################################################
   
   # Meteorological inputs:
   if (met_name == 'GEOSFP') {
      subdir = 'GEOS_FP/'
      file_ext = 'nc'
   } else if (met_name == 'MERRA2') {
      subdir = 'MERRA2/'
      file_ext = 'nc4'
   } else {
      stop('met_name specified is not available!')
   }
   
   subfn = paste0(subdir, YYYY, '/', MM, '/', met_name, '.', YYYY, MM, DD, '.A1.2x25.', file_ext)
   filename = paste0(met_data_dir, subfn)
   # These met fields are 1-hour average starting from 00:00 UTC of the day.
   
   # Open nc file:
   nc = nc_open(filename)
   
   # Load meteorological data:
   for (imet in 1:nrow(updated_input_df)) {
      
      # Get TEMIR and meteorological variable name:
      TEMIR_variable_name = as.character(updated_input_df$TEMIR_var_name[imet])
      met_variable_name = as.character(updated_input_df[,paste0(met_name, '_var_name')][imet])
      
      # Load particular meteorological data if variable exists in global meteorological field:
      if (!is.na(met_variable_name)) if (is.na(as.logical(met_variable_name))) assign(x = TEMIR_variable_name, value = ncvar_get(nc, met_variable_name)) else next
      
      # Get meteorological field dimension for FLUXNET conformity:
      if (!exists('met_dim') && FLUXNET_flag) met_dim = dim(get(x = TEMIR_variable_name))
      
      # Convert unit if required:
      if (as.logical(updated_input_df[,paste0(met_name, '_unit_differ')][imet])) assign(x = TEMIR_variable_name, value = f_met_unit_convert(met.name = met_name, TEMIR.var = TEMIR_variable_name, met.var = met_variable_name))
      
   }
   
   # Close nc file:
   nc_close(nc)
   
   # Load FLUXNET meteorology if FLUXNET_flag=TRUE:
   if (FLUXNET_flag) {
      
      # Read in FLUXNET data for simulation date in UTC: (Have not implement the data quality check provided by FLUXNET QC flag)
      start_ind = match(FLUXNET_current,FLUXNET_selected_data$TIMESTAMP_START)
      end_ind = match(FLUXNET_end,FLUXNET_selected_data$TIMESTAMP_START) - 1
      # Check if FLUXNET data extraction is successful
      if (start_ind < 0 || is.na(end_ind)) stop(paste('Error in FLUXNET Data input range!!!', start_ind, end_ind))
      
      # Select current FLUXNET data:
      FLUXNET_day_data = FLUXNET_selected_data[start_ind:end_ind,]
      FLUXNET_var_err = NULL
      print('Extracting corresponding FLUXNET data for each variable....', quote = FALSE)
      # NOTE : Calculations of scaling for variable agreement between the FLUXNET and meteorological data is done in the function f_FLUXNET_convert in FLUXNET_functions.R
      
      # Load meteorological data:
      for (imet in 1:nrow(updated_input_df)) {
         
         # Get TEMIR and meteorological variable name:
         TEMIR_variable_name = as.character(updated_input_df$TEMIR_var_name[imet])
         current_FLUXNET_met = as.character(updated_input_df$FLUXNET_var_name[imet])
         
         # Check if FLUXNET meteorlogical data is available for the particular variable
         if (is.na(current_FLUXNET_met)) next
         
         # Set meteorological field array if variable has not been declared:
         if (!exists(TEMIR_variable_name)) assign(x = TEMIR_variable_name, value = array(data = NA, dim = met_dim))
         
         # Load particular meteorological data:
         temporary = f_FLUXNET_met_grid(FLUXNET.data = FLUXNET_day_data, TEMIR.variable = TEMIR_variable_name, FLUXNET.variable = current_FLUXNET_met, FLUXNET.nrows = FLUXNET_nrows, dt.hr = dt_hr, st.ind = start_ind, end.ind = end_ind, input.var.df = updated_input_df, ind.lon = ind_lon, ind.lat = ind_lat, hh.to.hr.flag = HH_to_HR_flag)
         assign(x = updated_input_df$TEMIR_var_name[imet], value = temporary$data)
         
         # Get FLUXNET variable error if exists
         if (!is.null(temporary$err_out)) {
            if (is.null(FLUXNET_var_err)) {
               FLUXNET_var_err = rbind(temporary$err_out)
            } else {
               FLUXNET_var_err = rbind(FLUXNET_var_err, temporary$err_out)
            }
         }
         
      }
      
      # Make FLUXNET error for archive: (not supported by ncdf4 as matrix dimnames are lost when archiving nc)
      if (!is.null(FLUXNET_var_err) && is.na(match('FLUXNET_error', names(var_list)))) {
         err_nchar = ncdim_def(name='error_name_length', units='', vals=1:8, create_dimvar=FALSE)
         err_row = ncdim_def(name='error_number', units='', vals=1:nrow(FLUXNET_var_err),  longname='Number of Variable Error')
         err_col = ncdim_def(name='error_replace', units='', vals=1:2, longname='Error and Replacement Variables')
         FLUXNET_err_nc_def = ncvar_def(name='FLUXNET_data_error', units='', dim=list(err_nchar, err_row, err_col), longname='FLUXNET Variable Error and Replacement', prec='char', compression=4)
         var_list[['FLUXNET_error']] = FLUXNET_err_nc_def
      }
      
   }
   
   #############################################################################
   
   # Soil temperature for biogeochemistry (Pang)
   if (biogeochem_flag){
       # Specify the soil layer depth (unit: m) for soil temperature from the surface layer to the deepest layer, can be a number (single layer) or an array (multiple layers)
       soil_layer_depth = c(0.0988,0.2940,0.6799,1.4425,2.9496)
       
       filename = list.files(path = paste0(soilT_data_dir), pattern = paste0(YYYY,MM,DD), recursive = T)
       filename = paste0(soilT_data_dir, filename)
       
       nc = nc_open(filename)
       # Soil temperature 1st layer (0 - 9.88cm) (K)
       TSOIL1 = ncvar_get(nc,"TSOIL1")
       # Soil temperature 2nd layer (9.88cm - 29.4cm) (K)
       TSOIL2 = ncvar_get(nc,"TSOIL2")
       # Soil temperature 3rd layer (29.4cm - 68.0cm) (K)
       TSOIL3 = ncvar_get(nc,"TSOIL3")
       # Soil temperature 4th layer (68.0cm - 144cm) (K)
       TSOIL4 = ncvar_get(nc,"TSOIL4")
       # Soil temperature 5th layer (144cm - 295cm) (K)
       TSOIL5 = ncvar_get(nc,"TSOIL5")
       nc_close(nc)
   }
   
   # GDD map for determining GDDmat
   if (biogeochem_flag) {
       if (get_GDDmat_method == "CLM4.5") {
           # For simplicity, we decide not to implement the moving average of GDDx like the one in CLM4.5, as the prediction of the change in planting/harvesting date are not very accurate anyway
           # We only use the GDDx map in yr 2000 and calculate the corresponding GDDmat
           filename = paste0(GDDx_map_dir,'MEERA2_year_2000_growing_season_GDDx_map.nc')
           nc = nc_open(filename)
           GDD0_map = ncvar_get(nc,"GDD0")
           GDD8_map = ncvar_get(nc,"GDD8")
           GDD10_map = ncvar_get(nc,"GDD10")
           nc_close(nc)
       } else if (get_GDDmat_method == "Sack") {
           # The data is read in at PFT_surf_data.R
       } else if (get_GDDmat_method == "custom") {
           print("'get_GDDmat_method' is custom, will not read in GDDx map")
       }
   }
   
   #### soyFACE project
   if (biogeochem_flag) {
       if (force_prescribed_LAI) {
           daily_LAI = mean(force_prescribed_LAI_input[((d-1)*24+1):((d-1)*24+24)], na.rm = T)
           if (d != 1) {
               LAI_dayMinus1 = mean(force_prescribed_LAI_input[((d-2)*24+1):((d-2)*24+24)], na.rm = T)
           } else {
               LAI_dayMinus1 = 0
           }
       }
   }
   
   #############################################################################
   environment(f_simulate_ij) = globalenv()
   environment(f_hist_reshape) = globalenv()
   
   # Simulate for each lon/lat:
   print('Simulating for each lon/lat...', quote=FALSE)
   if (multicore_flag) hist_ij = mclapply(ij, FUN=f_simulate_ij, mc.cores=n_core) else hist_ij = lapply(ij, FUN=f_simulate_ij)
   print(paste0('Done on ', Sys.time()), quote=FALSE)
   
   if (debug_flag) {
      
      # "hist_grid" is where all the output data are.
      # Its dimensions: hist_grid = array(NaN, dim=c(length(ind_lon), length(ind_lat), n_day_sim, length(pftname), 24/dt_hr, nrow(var_name)))
      
      output = f_hist_reshape(ij=ij, hist_ij=hist_ij)
      hist_grid[,,d,,,] = output$hist_grid
      err_hist_ij[[d]] = output$err_hist_ij
      err_msg = c(err_msg, output$err_msg)
      
   } else {
      
      print('Reshaping and saving history data into gridded data file...', quote=FALSE)
      output = f_hist_reshape(ij=ij, hist_ij=hist_ij)
      hist_grid = output$hist_grid
      err_hist_ij = output$err_hist_ij
      err_msg = output$err_msg
      
      if (archive_format == 'RData') {
         
         # Save output history data in RData:
         filename = paste0(simulation_dir, 'hist_data/hist_grid_', YYYY, MM, DD, '.RData')
         FLUXNET_error_vec = if (FLUXNET_flag) c('FLUXNET_var_err', 'FLUXNET_global_err') else NULL
         save(list=c('hist_grid', 'err_hist_ij', 'err_msg', FLUXNET_error_vec), file=filename)
         
      } else if (archive_format == 'nc') {
         
         # Get history nc file name:
         filename = paste0(simulation_dir, 'hist_data/hist_grid_', YYYY, MM, DD, '.nc')
         
         print(paste0('nc filename = ', filename))
         
         # Delete previous history nc file:
         if (file.exists(filename)) file.remove(filename)
         
         # Create history nc file:
         nc = nc_create(filename=filename, vars=var_list)
         
         # Put in values of variables:
         for (v in 1:nrow(var_name)) {
            if (var_name[v,4] == 'grid') {
                    ncvar_put(nc=nc, varid=var_list[[v]], vals= hist_grid[,,1,,v])
                } else if (var_name[v,4] == 'PFT') {
                    ncvar_put(nc=nc, varid=var_list[[v]], vals=hist_grid[,,,,v])
                } else if (var_name[v,4] == 'PFT_daily') {
                    ncvar_put(nc=nc, varid=var_list[[v]], vals=hist_grid[,,,24,v])
                }
         }
         
         # Put in dimentions and attributes:
         ncatt_put(nc, varid='lon', attname='axis', attval='X')
         ncatt_put(nc, varid='lat', attname='axis', attval='Y')
         ncatt_put(nc, varid='pft', attname='axis', attval='P')
         ncatt_put(nc, varid='hour', attname='axis', attval='T')
         ncatt_put(nc, varid=0, attname='Title', attval=paste0(basename(simulation_dir), ' ', YYYY, MM, DD))
         ncatt_put(nc, varid=0, attname='Conventions', attval='COARDS')
         ncatt_put(nc, varid=0, attname='History', attval=paste0('Generated on ', Sys.time()))
         
         # Disable at the time being (Pang, Apr 2020)
         # Apparently we cannot wrtie logical values into the ncdf4
         
         # # FLUXNET data warnings:
         # if (single_site_flag && FLUXNET_flag) {
         #    ncatt_put(nc, varid=0, attname='FLUXNET Site Warnings', attval=FLUXNET_global_err)
         #    if (!is.null(FLUXNET_var_err)) {
         #       ncvar_put(nc, varid=FLUXNET_err_nc_def, vals=FLUXNET_var_err)
         #       ncatt_put(nc, varid='error_number', attname='axis', attval='n_error')
         #       ncatt_put(nc, varid='error_replace', attname='axis', attval='error_replace')
         #       rm(FLUXNET_var_err)
         #    }
         # }
         
         # Close nc file:
         nc_close(nc)
         
         # Save error messages:
         if (!is.null(err_msg) || length(err_hist_ij) != 0) {
            filename = paste0(simulation_dir, 'hist_data/hist_err_', YYYY, MM, DD, '.RData')
            save(list=c('err_hist_ij', 'err_msg'), file=filename)
         }
         
      } else stop('Debugging mode is off but data archiving format is not correctly specified.')
      
      print(paste0('Done on ', Sys.time()), quote=FALSE)
      
   }
   
   # Delete temporary data 15 days ago to prevent excessive amount of data:
   if (d > 15) {
      d_prev = 15
      previous_date = to.yyyymmdd(from.yyyymmdd(current_date) - d_prev*24)
      pathname = paste0(simulation_dir, 'temp_data/temp_', as.character(previous_date))
      #system(command = paste0('rm -rf ', pathname))
      dir.remove(pathname)
   }
}

################################################################################

# Save model configuration:

if (!debug_flag) {
   filename = paste0(simulation_dir, 'hist_data/model_config.RData')
   execution_config = c('dt', 'ind_lon', 'ind_lat', 'ij',  'n_day_sim', 'pftname', 'pftnum', 'PFT_frac')
   out_config = c(execution_config, model_config_vec)
   save(list = out_config, file=filename)
}

# End of simulation for all days.
print('### End of Simulation ###', quote=FALSE)
timestamp()

################################################################################
### End of execution
################################################################################
