# Temporary
force_prescribed_LAI = FALSE
force_prescribed_o3 = FALSE; O3_damage_flag = FALSE
sim_year = as.numeric(Sys.getenv("bash_year"))
o3_multiplier = 1
simulation_site = 'US_Ne2'   # 'US_Ne2', 'US_Ne3' or 'soyFACE

if (force_prescribed_LAI && simulation_site == 'soyFACE') {
   
    # sim_year = 2003
    
    # load('/users/b146986/TEMIR/TEMIR_inputs/soyFACE_20012008_interpolated_hourly_LAI.RData')
    load('/users/b146986/TEMIR/TEMIR_inputs/soyFACE_20092010_interpolated_hourly_LAI.RData')
    
    ####### variable name
    ###  year 2001 - 2008
    # prescribed_LAI_2001_ctr; _eco2; _eo3; _eco2o3
    # force_prescribed_LAI_input = prescribed_LAI_2008_ctr
    
    # force_prescribed_LAI_input = eval(parse(text = paste0('prescribed_LAI_',sim_year,'_ctr')))
    
    ### year 2009 - 2010
    # prescribed_LAI_2009_amb; _40; _55; _70; 85; 110; 130; 160; 200
    # prescribed_LAI_2010_amb; _55; _70; _85; _110; _130: _150; _170; _190
    
    o3_case_name_index = as.numeric(Sys.getenv("bash_o3_index")) - 3
    if (sim_year == 2009) {
        o3_case_name = c('_amb','_40','_55','_70','_85','_110','_130','_160','_200')
        print(paste0('Year = ', sim_year, ' o3_index = ', o3_case_name_index))
    } else if (sim_year == 2010) {
        o3_case_name = c('_amb','_55','_70','_85','_110','_130','_150','_170','_190')
        print(paste0('Year = ', sim_year, ' o3_index = ', o3_case_name_index))
    }
    force_prescribed_LAI_input = eval(parse(text = paste0('prescribed_LAI_',sim_year,o3_case_name[o3_case_name_index])))
    # force_prescribed_LAI_input = eval(parse(text = paste0('prescribed_LAI_',sim_year,'_amb')))
    print(paste0('Loading LAI vector: ', paste0('prescribed_LAI_',sim_year,'_amb')))
    # print(paste0('Loading LAI vector: ', 'prescribed_LAI_',sim_year,o3_case_name[o3_case_name_index]))
    rm (o3_case_name_index)
}

if (force_prescribed_o3 & O3_damage_flag) {
    
    if (simulation_site == 'soyFACE') {
        
        ### yr 2001 - 2008: O3_hourly = ozone_vec_ppb
        ### yr 2009 - 2010: O3_hourly = ozone_hourly_df_xxxx$amb
        
        if (sim_year >= 2001 & sim_year <= 2008) {
            filename = paste0('/users/b146986/TEMIR/TEMIR_inputs/o3_input/soyFACE/','soyFACE_ozone_',sim_year,'_processed.RData')
            load(filename)
            ozone_soyFACE_vec = ozone_vec_ppb * o3_multiplier
            print(paste0('Loading O3 file ',filename))
            rm(modified_ozone_df, singleSite_ozone_df, ozone_vec_ppb)
            
        } else if (sim_year >= 2009 & sim_year <= 2010) {
            load('/users/b146986/TEMIR/TEMIR_inputs/o3_input/soyFACE/soyFACE_20092010_allRings_ozoneConc.RData')
            print('Loading O3 file /users/b146986/TEMIR/TEMIR_inputs/o3_input/soyFACE/soyFACE_20092010_allRings_ozoneConc.RData')
            
            tmp_df = eval(parse(text = paste0('ozone_hourly_df_',sim_year)))
            
            # 2009: amb 40 55 70 85 110 130 160 200 (index = 4:12)
            # 2010: amb 55 70 85 110 130 150 170 190 (index = 4:12)
            o3_column_index = as.numeric(Sys.getenv("bash_o3_index"))
            ozone_soyFACE_vec = tmp_df[,o3_column_index]
            print(head(ozone_soyFACE_vec))
            rm(tmp_df, ozone_hourly_df_2009, ozone_hourly_df_2010, o3_column_index)
            
        }
    } else if (simulation_site == 'US_Ne2' || simulation_site == 'US_Ne3') {
        
    }


} 

bashIndex = as.numeric(Sys.getenv("bash_index")) + 1

# yr 2001 - 2008
# # yr 2003 is now not omitted
if (simulation_site == 'soyFACE') {
    
    emergence_carbon_to_leaf_input = 0.6
    
    if (sim_year <= 2008) {
        soyFACE_planting_date_vec = c(143,152,157,149,145,145,142,169); soyFACE_planting_date_input = soyFACE_planting_date_vec[bashIndex]
        soyFACE_GDDemer_vec = c(139.62, 68.42, 44.85, 46.42, 49.80, 68.08, 75.70, 89.43); soyFACE_GDDemer_input = soyFACE_GDDemer_vec[bashIndex]
        soyFACE_GDDrepr_vec = c(927.56, 1050.4, 839.88, 710.22, 829.71, 715.25, 856.52, 716.41); soyFACE_GDDrepr_input = soyFACE_GDDrepr_vec[bashIndex]
        soyFACE_GDDmat_vec = c(1564.3, 1680.08, 1448.56, 1403.30, 1833.74, 1571.6, 1719.2, 1267.35); soyFACE_GDDmat_input = soyFACE_GDDmat_vec[bashIndex]    
    }
    
    if (sim_year >= 2009) {
        soyFACE_planting_date_vec = c(160, 147); soyFACE_planting_date_input = soyFACE_planting_date_vec[bashIndex]
        soyFACE_GDDemer_vec = c(51.4, 62.89); soyFACE_GDDemer_input = soyFACE_GDDemer_vec[bashIndex]
        # soyFACE_GDDemer_vec = c(110.4, 62.89); soyFACE_GDDemer_input = soyFACE_GDDemer_vec[bashIndex]
        soyFACE_GDDrepr_vec = c(775.4, 921.4); soyFACE_GDDrepr_input = soyFACE_GDDrepr_vec[bashIndex]
        soyFACE_GDDmat_vec = c(1262.0, 1677.5); soyFACE_GDDmat_input = soyFACE_GDDmat_vec[bashIndex]
    }
} else if (simulation_site == 'US_Ne2') {
    emergence_carbon_to_leaf_input = 0.75
    ### poor variable name here...
    # for year 2002 2004 2006 2008
    soyFACE_planting_date_vec = c(140, 154, 132, 136); soyFACE_planting_date_input = soyFACE_planting_date_vec[bashIndex]
    soyFACE_GDDemer_vec = c(53.36, 77.28, 80.90, 66.46); soyFACE_GDDemer_input = soyFACE_GDDemer_vec[bashIndex]
    # soyFACE_GDDemer_vec = c(110.4, 62.89); soyFACE_GDDemer_input = soyFACE_GDDemer_vec[bashIndex]
    soyFACE_GDDrepr_vec = c(1231.30, 1012.72, 1285.38, 1074); soyFACE_GDDrepr_input = soyFACE_GDDrepr_vec[bashIndex]
    soyFACE_GDDmat_vec = c(1866.61, 1437.94, 1796.36, 1550.70); soyFACE_GDDmat_input = soyFACE_GDDmat_vec[bashIndex]
} else if (simulation_site == 'US_Ne3') {
    emergence_carbon_to_leaf_input = 0.75
    ### poor variable name here...
    soyFACE_planting_date_vec = c(140, 154, 131, 135); soyFACE_planting_date_input = soyFACE_planting_date_vec[bashIndex]
    soyFACE_GDDemer_vec = c(53.36, 104.3, 66.90, 67.98); soyFACE_GDDemer_input = soyFACE_GDDemer_vec[bashIndex]
    # soyFACE_GDDemer_vec = c(110.4, 62.89); soyFACE_GDDemer_input = soyFACE_GDDemer_vec[bashIndex]
    soyFACE_GDDrepr_vec = c(1279.94, 838.77, 1165.55, 987.7); soyFACE_GDDrepr_input = soyFACE_GDDrepr_vec[bashIndex]
    soyFACE_GDDmat_vec = c(1871.04, 1434.67, 1818.86, 1550.57); soyFACE_GDDmat_input = soyFACE_GDDmat_vec[bashIndex]
}



print(paste0('[settings] GDDemer = ', soyFACE_GDDemer_input))
################################################################################
### Terrestrial Ecosystem Model in R (TEMIR)
### Input script for biogeochemistry (with crop module)
################################################################################

# Additional directories for biogeochemistry

# Initial data directory (mandatory)
# Initial data provides the starting condition of LAI and SAI for different PFTs 
# Future developement will involove reading in soil N and other variables
initial_data_dir = paste0(TEMIR_dir, 'TEMIR_inputs/initial_data/')

# Additional R scripts for biogeochemistry (mandatory)
BGC_code_dir = paste0(code_dir,"extension_crop/")

# Soil temperature directory for maintenance respiraion calculation (mandatory)
soilT_data_dir = '/lustre/project/TGABI/data/MERRA2_2x2.5_soilT/'

# Three additional optional directories settings at below
# planting_date_map_dir = ...
# GDDx_map_dir = ...
# Sack_GDDmat_dir = ...

optional_dirs = c(optional_dirs, 'planting_date_map_dir', 'GDDx_map_dir', 'Sack_GDDmat_dir')
################################################################################

### Simulation settings ###

# Select PFT for simulation
# Option can either be 'all_pft' (default), 'all_nat_veg_only', 'all_prog_crop_only', 'custom' 
# 'all_pft' = c(2,25)
# 'all_nat_veg_only' = c(2,17)
# 'all_prog_crop_only' = c(18,25)
# 'custom': enter a number or an array, the number is the index of the 'pftname' from PFT_surf_data.R

BGC_PFT_selection_option = 'custom' 
if (BGC_PFT_selection_option == 'custom') {BGC_pft_selection = c(24)}
if (!any(BGC_PFT_selection_option == c('all_pft','all_nat_veg_only','all_prog_crop_only','custom'))){stop("'BGC_PFT_selection_option' should either be 'all_pft', 'all_nat_veg_only', 'all_prog_crop_only' or 'custom'")}

################################################################################

### General settings for vegetation ###

# All options can either be 'CLM4.5' or 'custom'
# 'CLM4.5' (default): follow the calculation in CLM4.5
# 'custom' write your own scheme in ...
#  root_fraction: maintenance_respiration.R
#  LAI, SAI, canopy height: physiology.R

root_fraction_scheme = 'CLM4.5'
prognostic_LAI_calculation_scheme = 'custom'
prognostic_SAI_calculation_scheme = 'CLM4.5'
prognostic_canopy_height_calculation_scheme = 'CLM4.5'
if (!any(root_fraction_scheme == c('CLM4.5','custom'))){stop("'root_fraction_scheme' should either be 'CLM4.5' or 'custom'")}
if (!any(prognostic_LAI_calculation_scheme == c('CLM4.5','custom'))){stop("'prognostic_LAI_calculation_scheme' should either be 'CLM4.5' or 'custom'")}
if (!any(prognostic_SAI_calculation_scheme == c('CLM4.5','custom'))){stop("'prognostic_SAI_calculation_scheme' should either be 'CLM4.5' or 'custom'")}
if (!any(prognostic_canopy_height_calculation_scheme == c('CLM4.5','custom'))){stop("'prognostic_canopy_height_calculation_scheme' should either be 'CLM4.5' or 'custom'")}

################################################################################

### Natural vegetation settings ###

# Can either be 'CLM4.5' or custom'
# 'CLM4.5' (default): follow the description in CLM4.5 (in developement)
# 'custom' write your own scheme in biomass_partitioning.R
natveg_biomass_partitioning_scheme = 'CLM4.5'
if (!any(natveg_biomass_partitioning_scheme == c('CLM4.5','custom'))){stop("'natveg_biomass_partitioning_scheme' should either be 'CLM4.5' or 'custom'")}

################################################################################

### Crop module settings ###

### Planting date setting
# Can either be 'CLM4.5', 'prescribed-map' or 'prescribed-site'
# 'CLM4.5' (default): planting date is determined by climate (GDD) and 10-day running mean T2m
# 'prescribed-map': read in a .nc file that contains planting date (in julian day), by default we use the Sack planting/harvesting dataset
#   - crop_growing_season can either be 'primary' and 'secondary'
#   - Sack dataset has primary and secondary growing season for maize and rice
# 'prescribed-site': enter a single plant date (in julian day), best for single site / small scale simulations with one planting date
#####prescribed_planting_date_flag 
get_planting_date_option = 'prescribed-site'
if (get_planting_date_option == 'prescribed-map'){
    planting_date_map_dir = paste0(TEMIR_dir,'TEMIR_inputs/Sack_crop_calender/')
    if (exists('crop_growing_season')) {
        warning("'crop_growing_season' for Sack planting dataset is declared in both 'input_TEMIR_basic_settings.R' and 'input_TEMIR_biogeochem_extension.R'\nIt is overwritten by the one decalred in input_TEMIR_biogeochem_extension.R")
    }
    crop_growing_season = 'primary'
} else if (get_planting_date_option == 'prescribed-site') {
    prescribed_planting_date = soyFACE_planting_date_input
} else if (get_planting_date_option == 'CLM4.5'){
   # nothing here
} else {
    stop("'get_planting_date_option' should either be 'CLM4.5', ''prescribed-site' or 'prescribed-map'")
}

### GDD requirement for crop to reach matuiry
# Can either be 'CLM4.5', 'Sack' or 'custom'
# 'CLM4.5': GDDmat are derived from 20-year running mean of GDD0, GDD8 and GDD10 as described in the CLM4.5 technical note. 
# 'Sack': GDDmat are the mean GDD accumulated during a growing season from year 1995 to year 2015, growing season is the period between the planting date and harvesting date from Sack et al. (2010)
# 'custom': enter a single GDD value (in degree C. day), best for single site / small scale simulations with one GDDmat
get_GDDmat_method = 'custom'
if (get_GDDmat_method == 'CLM4.5') {
    GDDx_map_dir = paste0(TEMIR_dir, 'TEMIR_inputs/GDDx_map/')
} else if (get_GDDmat_method == 'Sack') {
    Sack_GDDmat_dir = paste0(TEMIR_dir,'TEMIR_inputs/Sack_GDDmat/')
} else if (get_GDDmat_method == 'custom') {
    prescribed_GDD_mat = soyFACE_GDDmat_input
} else {
    stop("'get_GDDmat_method' should either be 'CLM4.5', Sack' or 'custom'")
}

### GDD requirement for crop to reach reproductive stage
# Can either be 'CLM4.5' (default) or 'custom'
# 'CLM4.5': GDDrepr are derived from GDDmat (i.e., 50 - 70% of GDDmat)
# 'custom': enter a single GDD value (in degree C. day), best for single site / small scale simulations with one GDDrepr
get_GDDrepr_method = 'custom'
if (get_GDDrepr_method == 'CLM4.5'){
    # Nothing
} else if (get_GDDmat_method == 'custom'){
    prescribed_GDD_repr = soyFACE_GDDrepr_input
} else {
    stop("'get_GDDrepr_method' should either be 'CLM4.5' or 'custom'")
}

### GDD requirement for crop to reach vegetative stage
# Can either be 'CLM4.5' or 'custom'
# 'CLM4.5': GDDemer are derived from GDDmat (i.e., 3 - 5% of GDDmat)
# 'custom': enter a single GDD value (in degree C. day), best for single site / small scale simulations with one GDDemer
get_GDDemer_method = 'custom'
if (get_GDDemer_method == 'CLM4.5'){
    # Nothing
} else if (get_GDDemer_method == 'custom'){
    prescribed_GDD_emer = soyFACE_GDDemer_input
} else {
    stop("'get_GDDemer_method' should either be 'CLM4.5' or 'custom'")
}

### GDD accumulation rate calculation
# If (TRUE): write your own GDD accumulation scheme in phenology.R
# If (FALSE): GDD increase rate is based on CLM4.5 crop module (default)
user_defined_crop_GDD_accmulation_flag = FALSE

### Allocation setting option for crops ###

# Can either be 'CLM4.5' (default), 'JULES' or 'custom'
# 'CLM4.5': based on the allocation coefficients in CLM4.5
#  - limit_crop_LAI_flag = T, crop LAI follows the original CLM4.5 in which is constrained by laimx from PFT_surf.data, this flag controls allocation coefficients such that aleaf = 0 when LAI >= laimx
#  - limit_crop_LAI_flag = F (default), this constrain is removed
# 'JULES': based on the allocation coefficients in JULES
# 'custom' implement your biomass partitioning function in biomass_partitioning.R
crop_biomass_partitioning_scheme = 'custom'
if (crop_biomass_partitioning_scheme == 'CLM4.5'){
    limit_crop_LAI_flag = FALSE
}

### Translocation setting option for crops ###
# Can either be 'CLM4.5', 'JULES' (default) or 'custom'
# 'CLM4.5': based on the retranslocation scheme for crops in CLM4.5, this scheme requires nitrogen cycle (in developement). 
# 'JULES': based on the translocation scheme for crops in JULES-crop, this scheme assumes carbon from leaf senescence is transfered to grain during the reproductive stage. Noted that the DVI (phenology) requirement is reduced from 1.5 (in JULES) to 1.0.
# 'custom' write your own translocation requirement, or simply ignore it
crop_translocation_scheme = 'JULES'
if (!any(crop_translocation_scheme == c('CLM4.5', 'JULES', 'custom'))) {stop("'crop_translocation_scheme' should either be 'CLM4.5', JULES' or 'custom'")}

### Leaf senescence scheme for crops ###
# Can either be 'CLM4.5' (default) or 'custom'
# 'CLM4.5': based on the leaf senescence rate described in CLM4.5 Technical Note, leaf_senescence_rate = - (leaf_long * leafC) / 365 * 86400)
# 'custom': write your own leaf senescence scheme in f_crop_phenology() at plant_phenology.R 
crop_leaf_senescence_scheme = 'custom'
if (crop_leaf_senescence_scheme == 'custom') {
    leaf_senescence_DVI = 1.0
} else {
    leaf_senescence_DVI = NA
}

if (!any(crop_leaf_senescence_scheme == c('CLM4.5','custom'))) {stop("'crop_leaf_senescence_scheme' should either be 'CLM4.5' or 'custom'")}

### Starting leafC (i.e., LAI) for crops after emergence ###
# Enter numeric value for leafC (unit: gC m^-2)
# It affects the peak LAI that the crop can attain, potentially influences yield simulation
# Default: 1 gC m^-2 for all crops (from CLM4.5)
emergence_carbon_to_leaf = emergence_carbon_to_leaf_input

################################################################################
# Additional biogeochemistry (crop module) output to RData / ncdf4, modify this if necessary
output_variables = c(output_variables,         # output from input_basic_settings.R
                     'LAI', 'SAI', 'grainC', 'GPP', 'NPP', 'leafC', 'livestemC',
                     'GDDT2m', 'GDDTsoil', 'aleaf', 'astem', 'aroot', 'arepr',
                     # 'crop_live_flag', 'crop_plant_flag', 'leaf_emergence_flag', 'grain_fill_flag', 'harvest_flag',
                     'day_of_planting', 'day_of_harvesting'
                     # 'CUO_grainfill'
                     )

# Additional biogeochemistry (crop module) available output
available_outputs_df = rbind.data.frame(
    available_outputs_df,
    # Physiology variables
    c('LAI', 'm^2 m^-2', 'Total leaf area index', 'PFT_daily'),
    c('SAI', 'm^2 m^-2', 'Total stem area index', 'PFT_daily'),
    c('htop', 'm', 'Canopy height', 'PFT_daily'),
    c('hbot', 'm', 'Canopy bottom', 'PFT_daily'),
    # Plant carbon pools
    c('leafC', 'gC m^-2', 'Leaf carbon stock', 'PFT_daily'),
    c('finerootC', 'gC m^-2', 'Fine root carbon stock', 'PFT_daily'),
    c('livestemC', 'gC m^-2', 'Live stem carbon stock', 'PFT_daily'),
    c('deadstemC', 'gC m^-2', 'Dead stem carbon stock', 'PFT_daily'),
    c('livecoraserootC', 'gC m^-2', 'Live coarse root carbon stock', 'PFT_daily'),
    c('deadcoraserootC', 'gC m^-2', 'Dead coarse root carbon stock', 'PFT_daily'),
    c('grainC', 'gC m^-2', 'Grain carbon stock (crop dry yield)', 'PFT_daily'),
    # GPP, NPP, biomass partitioning fluxes and respirations
    c('GPP', 'gC m^-2 s^-1', 'Daily mean GPP', 'PFT_daily'),
    c('NPP', 'gC m^-2 s^-1', 'Daily mean NPP (GPP - growth resp. - main. resp.)', 'PFT_daily'),
    c('leafC_alloc', 'gC m^-2 s^-1', 'Biomass partitioning rate to leafC', 'PFT_daily'),
    c('finerootC_alloc', 'gC m^-2 s^-1', 'Biomass partitioning rate to finerootC', 'PFT_daily'),
    c('livestemC_alloc', 'gC m^-2 s^-1', 'Biomass partitioning rate to livestemC', 'PFT_daily'),
    c('deadstemC_alloc', 'gC m^-2 s^-1', 'Biomass partitioning rate to deadstemC', 'PFT_daily'),
    c('livecoraserootC_alloc', 'gC m^-2 s^-1', 'Biomass partitioning rate to livecoraserootC', 'PFT_daily'),
    c('deadcoraserootC_alloc', 'gC m^-2 s^-1', 'Biomass partitioning rate to deadcoraserootC', 'PFT_daily'),
    c('grainC_alloc', 'gC m^-2 s^-1', 'Biomass partitioning rate to grainC', 'PFT_daily'),
    c('mr_leaf', 'gC m^-2 s^-1', 'Leaf maintenance respiration', 'PFT'),
    c('mr_fineroot', 'gC m^-2 s^-1', 'Fine root maintenance respiration', 'PFT'),
    c('mr_livestem', 'gC m^-2 s^-1', 'Live stem maintenance respiration', 'PFT'),
    c('mr_livecoraseroot', 'gC m^-2 s^-1', 'Live corase root maintenance respiration', 'PFT'),
    c('mr_grain', 'gC m^-2 s^-1', 'Grain maintenance respiration', 'PFT'),
    c('mr_total', 'gC m^-2 s^-1', 'Total maintenance respiration', 'PFT'),
    # Crop phenology
    c('GDDT2m', 'Degree day', 'Growing degree day of T2m', 'PFT_daily'),
    c('GDDTsoil', 'Degree day', 'Growing degree day of Tsoil', 'PFT_daily'),
    c('GDDmat' , 'Degree day', 'Growing degree day requirement to reach crop maturity', 'PFT_daily'),
    c('GDDDrepr', 'Degree day', 'Growing degree day requirement to reach grain fill (reproductive stage)', 'PFT_daily'),
    c('GDDemer', 'Degree day', 'Growing degree day requirement to reach leaf emergence (vegetative stage)', 'PFT_daily'),
    c('crop_live_flag', '', 'Flag determines whether crop is live', 'PFT_daily'),
    c('crop_plant_flag', '', 'Flag determines whether crop is planted already', 'PFT_daily'),
    c('leaf_emergence_flag', '', 'Flag determines whether crop is in the vegetaive stage', 'PFT_daily'),
    c('grain_fill_flag', '', 'Flag determines whether crop is in the reproductive stage', 'PFT_daily'),
    c('harvest_flag', '', 'Flag determines whether crop is harvested', 'PFT_daily'),
    c('day_of_planting', 'Julian day', 'Crop planting day', 'PFT_daily'),
    c('day_of_grain_filling', 'Julian day', 'Crop flowering day (start of reproduction stage)', 'PFT_daily'),
    c('day_of_harvesting', 'Julian day', 'Crop harvesting day', 'PFT_daily'),
    # Allocation coefficients
    c('aleaf','','Biomass allocation fraction to leaf', 'PFT_daily'),
    c('aleaf_leafem', '', 'aleaf during the vegetative stage', 'PFT_daily'),
    c('astem','','Biomass allocation fraction to livestem', 'PFT_daily'),
    c('astem_leafem', '', 'astem during the vegetative stage', 'PFT_daily'),
    c('aroot','','Biomass allocation fraction to fine root', 'PFT_daily'),
    c('arepr','','Biomass allocation fraction to grain (reproductive)', 'PFT_daily')
)

# Ozone plant dammage, for soyFACE only at the time being
if (O3_damage_flag & O3_damage_scheme == 'Lombardozzi' & O3_sensitivity == 'custom') {
    available_outputs_df = rbind.data.frame(
        available_outputs_df,
        # Canopy ozone uptake
        c('CUO_grainfill', 'mmol m^-2', 'Cumulative O3 uptake when grain fill starts (crops only)', 'PFT_daily')
    )
}

# End of input_TEMIR_biogeochem_extension.R
################################################################################

# Declaring missing flags and values, do not modify this part
biogeochem_flag = TRUE
if (!exists("planting_date_map_dir")){planting_date_map_dir = NA}
if (!exists("prescribed_planting_date")){prescribed_planting_date = NA}
if (!exists("GDDx_map_dir")){GDDx_map_dir = NA}
if (!exists("Sack_GDDmat_dir")){Sack_GDDmat_dir = NA}
if (!exists("prescribed_GDD_mat")){prescribed_GDD_mat = NA}
if (!exists("prescribed_GDD_repr")){prescribed_GDD_mat = NA}
if (!exists("prescribed_GDD_emer")){prescribed_GDD_mat = NA}
if (!exists("limit_crop_LAI_flag")){limit_crop_LAI_flag = FALSE}

if (BGC_PFT_selection_option == 'all_PFT'){
    BGC_pft_selection = 2:length(pftname)
} else if (BGC_PFT_selection_option == 'all_nat_veg_only') {
    BGC_pft_selection = 2:17
} else if (BGC_PFT_selection_option == 'all_prog_crop_only') {
    BGC_pft_selection = 18:25
}

print('Finished setting up BGC simulation')
