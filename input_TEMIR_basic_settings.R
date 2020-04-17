################################################################################
### Terrestrial Ecosystem Model in R (TEMIR)
### Input script for single-site, regional or global simulation
################################################################################

# Computing environment:

# Use multiple cores?
multicore_flag = FALSE

# Run in cluster environment (i.e., not on personal computer)?
cluster_flag = TRUE

# Number of cores to use:
# This will be overwritten by system-dependent setting if "cluster_flag=TRUE".
# Use parallel::detectCores() for all available cores
n_core = detectCores()

# Debugging mode?
# If turned on, no output data will be saved externally. Debugging mode should be run only in an R console environment (i.e., not via a UNIX command line), and is convenient for single-site simulation.
debug_flag = FALSE

################################################################################
### Directories:
################################################################################

# Required directories:

# Set TEMIR directory:
TEMIR_dir = '/users/b146986/TEMIR/'
# Set source code directory:
code_dir = paste0(TEMIR_dir, 'code_v1.1/')
# Set meteorological data directory:
met_data_dir = '/lustre/project/TGABI/data/GEOS_2x2.5/'
# Set PFT and surface data directory:
surf_data_dir = paste0(TEMIR_dir, 'TEMIR_inputs/surf_data/')
# Set processed PFT and surface output directory:
processed_surf_data_dir = paste0(TEMIR_dir, 'TEMIR_inputs/processed_surf_data/')

# Optional directories:

# Set FLUXNET directory (set as NA if this is not available):
FLUXNET_dir = '/lustre/project/TGABI/data/FLUXNET/'

# Whatever are specified below, they will be overridden by "O3_data_dir" and "LAI_data_dir" set in the input script. They are preliminarily defined here to avoid problems during checking directory paths (below) when the execution script is re-run in the same console. (Tai, Feb 2019)
# Set user-defined LAI data directory (containing LAI nc files):
LAI_data_dir = NA
# Set ozone data directory:
O3_data_dir = paste0(TEMIR_dir, 'TEMIR_inputs/o3_input/')

optional_dirs = c('FLUXNET_dir', 'LAI_data_dir', 'O3_data_dir')

################################################################################
### Model configuration:
################################################################################

# Choice of meteorological fields:
# It should be either "GEOSFP" or "MERRA2".
# met_name = 'GEOSFP'
met_name = 'MERRA2'

# Model dimensions:
# Time step (hr):
dt_hr = 1
# Global longitude and latitude grid:
# It should always match the resolution of input meteorological fields (which are set in this file and in "PFT_surf_data.R").
dlon = 2.5
lon = seq(-180, 177.5, by=dlon)
dlat = 2.0
lat = seq(-90, 90, by=dlat)

# Model run dates:
# (temporary) read from bash
sim_year = as.numeric(Sys.getenv("bash_year"))
start_date = as.numeric(paste0(sim_year,'0101'))
end_date = as.numeric(paste0(sim_year,'1231')); rm(sim_year)
# start_date = 20030101
# end_date = 20031231

# Continue from previous run?
# If set true, temporary data within 10 days before the start date are needed.
continue_flag = FALSE

# Specify lon/lat for a single site or a given region (including global), or specify site ID if simulation is for a FLUXNET site.

# Single site or not?
single_site_flag = TRUE

# Running FLUXNET site?
FLUXNET_site_flag = TRUE
# Running simulation with FLUXNET met/canopy data?
FLUXNET_flag = TRUE
# Specify FLUXNET site ID (ignored if FLUXNET_site_flag=FALSE):
FLUXNET_site_id = "US-Ne2"
# # Saving FLUXNET meteorological data as ncdf? NOT SUPPORTED YET
# save_FLUXNET_data_flag = TRUE

# Specify location (lon, lat) of local site of interest:
if (single_site_flag && !FLUXNET_site_flag) {
   lon_sim = -88.14
   lat_sim = 40.02
} else {
   # Specify range of lon/lat of a given region of interest:
   # lon_sim = c(30, 40)
   # lat_sim = c(25, 50)
   lon_sim = c(-180,180)
   lat_sim = c(-90,90)
}

################################################################################

### Basic ecosystem model parameters ###

# Ambient CO2 concentration (ppm):
# (temporary) read from bash
CO2_conc = as.numeric(Sys.getenv("bash_co2"))
# CO2_conc = 550

# Fixed photosynthetic parameters:
# Nitrogen extinction coefficient:
K_n = 0.30
# Quantum yield of photosystem II (mol mol^-1):
Phi_PSII = 0.85
# Curvature parameter:
Theta_PSII = 0.70

# Stomatal conductance scheme: Farquhar-Ball-Berry model ('FBB') or Medlyn model ('Medlyn')
gs_scheme = 'FBB'

# Use simplified radiative trasnfer model?
# If not, the default model with two-stream approximation and light exctinction by both leaves and stems will be used to calculate variables "phi_sun", "phi_sha", "LAI_sun",  "LAI_sha" and "K_b", among others.
# If so, these five variables will be replaced by values from a simplified model assuming spherical leaf orientation and extinction by leaves only, to be consistent with Zhang et al. [2003] dry deposition scheme.
simple_radiation_flag = FALSE

# Use Monin-Obukhov theory to infer temperature and humidity in canopy air?
# If not, temperature and humidity at 2 m above displacement height will be used as proxies for temperature and humidity in canopy air.
# NOTE: if FLUXNET_flag=TRUE (so FLUXNET data are used for simulation), Monin_Obukhov_flag=FALSE by default, and what is set here will be ignored.
# NOTE: This method always uses the default TEMIR aerodynamic conductance scheme, never the scheme used in "drydep_toolbox.R"; "use_TEMIR_ga_flag" below will be set TRUE if Monin_Obukhov_flag=TRUE.
Monin_Obukhov_flag = FALSE

################################################################################

### Dry deposition parameters ###

# Which dry deposition scheme to use (either 'Wesely', 'Zhang' or NULL)?
# Set drydep_scheme=NULL if dry deposition computation is to be turned off.
drydep_scheme = 'NULL'

# Use default TEMIR aerodynamic conductance for scalars?
# If so, the default method in "Monin_Obukhov.R" will be used.
# If not, the method used in "drydep_toolbox.R" that is consistent with GEOS-Chem (Wesely) and Zhang et al. [2003] dry deposition schemes will be used.
# Please that this flag will be used for calculating ozone damage as well.
# If Monin_Obukhov=TRUE, this flag to be set TRUE, overriding what is set here.
use_TEMIR_ga_flag = TRUE

# Use default TEMIR stomatal conductance?
# If so, the default method in "Farquhar_Ball_Berry.R" will be used.
# If not, stomatal conductance is calculated using GEOS-Chem (Wesely) or Zhang et al. [2003] methods depending on "drydep_scheme" set above.
use_TEMIR_gs_flag = TRUE

# Use TEMIR water stress to scale stomatal resistance in Zhang et al. [2003]?
# If not, water stress follows Zhang et al. [2003] method.
use_TEMIR_beta_flag = FALSE

# CO2 scaling function for dry deposition velocity calculations:
CO2_scale_flag = TRUE

################################################################################

### POD3 simulation ###

O3_POD = FALSE
if (O3_POD) {
    # primary or secondary growing season used in Sack planting date/ harvesting date dataset
    # only applied for maize and rice
    crop_growing_season = 'primary'
}

################################################################################

### Ozone damage scheme ###

# Turn on ozone damage?
O3_damage_flag = FALSE

# Ozone damage scheme ('Lombardozzi' or 'Sitch'):
O3_damage_scheme = 'Lombardozzi'
# Ozone sensitivity:
# For "Lombardozzi" scheme: "high", "average" or "low" (added 'custom' for soyFACE simultaions, Feb 2020, Pang)
# For "Sitch" scheme: "high" or "low"
O3_sensitivity = 'custom'

# Use default TEMIR aerodynamic conductance for scalars to calculate O3 flux?
# Set use_TEMIR_ga_flag in the dry deposition menu above.

# Fixed ozone concentration?
O3_fixed_flag = FALSE
# If ozone damage scheme is turned off, fixing ozone concentration can save computational time, and no hourly ozone fields need to be provided:
if (!O3_damage_flag) O3_fixed_flag = TRUE

# Set fixed ozone concentration (ppb) if "O3_fixed_flag=TRUE":
if (O3_fixed_flag) O3_conc = 40

if (O3_damage_flag & !O3_fixed_flag) {
   
   # Set ozone data directory that contains the hourly O3 nc files:
   O3_data_dir = paste0(TEMIR_dir, 'TEMIR_inputs/o3_input/')
   
   # Set ozone nc filename:
   # Dimensions of ozone concentration array must be: [longitude, latitude, hours of simulation year]
   # Thus, ozone damage can only be simulated year by year, and ozone nc file has to be specified here manually.
   # E.g., if naming convention is "CTRL_hourlyO3_2010.nc", the following subnames should be specified as:
   # O3_subn1 = 'CTRL_hourlyO3_'
   # O3_subn2 = '.nc'
   # Leaving the year out which is to be replaced by the simulation year.
   O3_subn1 = 'CTRL_hourlyO3_'
   O3_subn2 = '.nc'
   
   # Set ozone array name in the nc file:
   O3_array_name = 'houro3'
   # Set ozone dimensions in a named vector:
   O3_dim_vec = c(longitude = 'lon', latitude = 'lat', time = 'time')
   
   # Please note that if the O3 file for the simulation year is not available, the model would not stop, but would look for an available year that is the closest to the simulation year.
   
}

################################################################################

### Leaf area index ###

# Use interannually varying user-defined LAI data as input?
# If not (LAI_data_flag=FALSE), default LAI used would be satellite phenology (SP) LAI from the Community Land Model (CLM) for year 2000.
# If so (LAI_data_flag=TRUE), the default PFT-level LAI from CLM-SP would be scaled to the user-defined (total, grid-level) LAI so that, when weigthted by PFT fractions and summed to obtain grid-level LAI, it would match the user-defined LAI.
LAI_data_flag = FALSE

if (LAI_data_flag) {
   
   # Set monthly LAI data directory that contains the LAI nc files:
   LAI_data_dir = paste0(TEMIR_dir, 'TEMIR_inputs/LAI_data/MODIS_LAI_201707/For_Olson_2001/')
   
   # LAI nc file naming convention:
   # Dimensions of LAI array must be: [longitude, latitude, months of simulation year]
   # Thus, each nc file is assumed to contain all monthly data in a given year.
   # E.g., if naming convention is "MODIS.LAIv.V5.generic.025x025.2010.nc", the following subnames should be specified as:
   # LAI_subn1 = 'MODIS.LAIv.V5.generic.025x025.'
   # LAI_subn2 = '.nc'
   # Leaving the year out which is to be replaced by the simulation year.
   LAI_subn1 = 'MODIS.LAIv.V5.generic.025x025.'
   LAI_subn2 = '.nc'
   
   # Set LAI array name in the nc file:
   LAI_array_name = 'MODIS'
   # Set LAI dimensions in a named vector:
   LAI_dim_vec = c(longitude = 'lon', latitude = 'lat', time = 'time')
   # Time indices to extract the twelve months of data in case LAI file contains extra data in their time dimension:
   ind_LAI_time = 1:12
   
   # Please note that if the LAI file for the simulation year is not available, the model would not stop, but would look for an available year that is the closest to the simulation year.
   
}

################################################################################
################################################################################

# Output/history data archiving:

# Data format to archive history data:
# No data will be archived if debugging mode is on.
# Only two options: "RData" or "nc"
archive_format = 'nc'

# Variables to archive:
# Variable names must match the variable name given in available_outputs below
# output_variables = c('CUO_can', 'CUO_sun', 'CUO_sha')
output_variables = c()

# Recommended outputs:
if (TRUE) output_variables = c(output_variables, 
                               'A_can', 'R_can', 'g_can', 
                               # 'L_sun', 'L_sha', 'A_nsun', 'A_nsha',
                               'LAI_sun', 'LAI_sha', 'A_nsun', 'A_nsha',
                               'R_dsun', 'R_dsha', 'g_ssun', 'g_ssha', 'g_s',
                               'K_b', 'phi_sun', 'phi_sha',
                               'surf_alb_beam', 'surf_alb_diff',
                               'beta_t', 'T_a', 'q_a'
                               # 'CUO_can', 'CUO_sun', 'CUO_sha'
)

# If new variables are added, both "available_outputs_df" here and "output_assign_df" in "simulate_ij.R" have to be modified.
# All available outputs are shown below:
# Additional outputs shown in the section following are also available when selected scheme is active.
# 1st col = variable name; 2nd col = unit; 3rd col = long name; 4th col = whether variable is resolved at PFT-level ('PFT' if yes, 'grid' if not)
available_outputs_df = `colnames<-`(rbind.data.frame(
   # Canopy-integrated ecophysiological outputs
   c('A_can', 'umol CO2 m^-2 s^-1', 'Canopy photosynthesis', 'PFT'),
   c('R_can', 'umol CO2 m^-2 s^-1', 'Canopy respiration', 'PFT'),
   c('g_can', 'm s^-1', 'Canopy conductance', 'PFT'),
   # Partitioned sun-lit / shaded outputs
   # c('L_sun', 'm^2 m^-2', 'Total sunlit leaf area index', 'PFT'),
   # c('L_sha', 'm^2 m^-2', 'Total shaded leaf area index', 'PFT'),
   c('LAI_sun', 'm^2 m^-2', 'Total sunlit leaf area index', 'PFT'),
   c('LAI_sha', 'm^2 m^-2', 'Total shaded leaf area index', 'PFT'),
   c('A_nsun', 'umol CO2 m^-2 s^-1', 'Sunlit leaf net photosynthesis', 'PFT'),
   c('A_nsha', 'umol CO2 m^-2 s^-1', 'Shaded leaf net photosynthesis', 'PFT'),
   c('R_dsun', 'umol CO2 m^-2 s^-1', 'Sunlit leaf mitochondrial respiration', 'PFT'),
   c('R_dsha', 'umol CO2 m^-2 s^-1', 'Shaded leaf mitochondrial respiration', 'PFT'),
   c('g_ssun', 'm s^-1', 'Sunlit leaf stomatal conductance', 'PFT'),
   c('g_ssha', 'm s^-1', 'Shaded leaf stomatal conductance', 'PFT'),
   c('g_s', 'm s^-1', 'Leaf stomatal conductance weighted by sunlit & shaded fractions', 'PFT'),
   # Canopy radiative transfer variables
   c('K_b', 'unitless', 'Canopy light extinction coefficient', 'PFT'),
   c('phi_sun', 'W m^-2', 'Absorbed photosynthetically active radiation by sunlit leaves', 'PFT'),
   c('phi_sha', 'W m^-2', 'Absorbed photosynthetically active radiation by shaded leaves', 'PFT'),
   c('surf_alb_beam', '0-1', 'Surface albedo for direct beam visible light', 'PFT'),
   c('surf_alb_diff', '0-1', 'Surface albedo for diffuse visible light', 'PFT'),
   # Water stress variable
   c('beta_t', '0-1', 'Soil water stress function', 'PFT'),
   # Canopy micrometeorological variables
   c('T_a', 'K', 'Canopy air temperature', 'grid'),
   c('q_a', 'kg kg^-1', 'Canopy air specific humidity', 'grid'),
   # Dataframe settings
   stringsAsFactors = FALSE), c('variable_name', 'unit', 'long_name', 'res_level'))

# Conditional outputs:
# Additional ozone damage outputs:
if (O3_damage_flag) {
   available_outputs_df = rbind.data.frame(
      available_outputs_df,
      # Canopy ozone uptake
      c('CUO_can', 'mmol m^-2', 'Cumulative O3 uptake', 'PFT'),
      # # Partitioned sun-lit / shaded ozone uptake
      c('CUO_sun', 'mmol m^-2', 'Cumulative O3 uptake for sunlit leaves', 'PFT'),
      c('CUO_sha', 'mmol m^-2', 'Cumulative O3 uptake for shaded leaves', 'PFT')
   )
}

# Additional dry deposition outputs:
if (!is.null(drydep_scheme)) {
   available_outputs_df = rbind.data.frame(
      available_outputs_df,
      # Ozone dry deposition velocity
      c('v_d','m s^-1','Dry deposition velocity','PFT'),
      # Resistances
      c('r_a','s m^-1', 'Aerodynamic resistance','PFT'),
      c('r_b','s m^-1', 'Leaf boundary layer resistance','PFT'),
      c('r_s','s m^-1', 'Stomatal  resistance','PFT'),
      c('r_cut','s m^-1', 'Cuticular resistance','PFT'),
      c('ra_lower_canopy','s m^-1', 'Lower canopy aerodynamic resistance','PFT'),
      c('rc_ground','s m^-1', 'Ground surface resistance','PFT')
   )
   
   # Outputs for Wesely dry deposition scheme:
   if (drydep_scheme == 'Wesely') {
      available_outputs_df = rbind.data.frame(
         available_outputs_df,
         c('rc_lower_canopy','s m^-1', 'Lower canopy surface resistance','PFT'),
         c('ra_ground','s m^-1', 'Ground aerodynamic resistance','PFT')
      )
   }
   
   # Outputs for Zhang dry deposition scheme:
   if (drydep_scheme == 'Zhang') {
      available_outputs_df = rbind.data.frame(
         available_outputs_df,
         c('is_wet','', 'Wet or dry canopy','grid'),
         c('g_par','s m^-1', 'Unstressed canopy stomatal conductance','PFT'),
         c('g_t','0-1', 'Temperature factor','PFT'),
         c('g_w','0-1', 'Water stress factor','PFT'),
         c('g_vpd','0-1', 'Water vapour presssure deficit factor','PFT'),
         c('w_st','0-0.5', 'Stomatal blocking factor','PFT'),
         c('f_snow','0-1', 'Snow cover fraction','PFT')
      )
   }
} 

################################################################################
### End of input specification
################################################################################
