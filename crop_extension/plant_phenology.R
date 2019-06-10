f_crop_phenology = function(T_10_d, T_min_10_d, T_soil, T2m,
                            leafC, livestemC, finerootC, grainC, LAI, SAI,
                            GDD_T2m, GDD_Tsoil, GDDmat, GDDrepr, GDDemer,
                            crop_living_flag, crop_planting_flag, leaf_emer_flag, grainfill_flag, har_flag,
                            prescribed_planting_date_readin = NULL, planting_jday, harvest_jday,
                            GDD0_20yr = NULL, GDD8_20yr = NULL, GDD10_20yr = NULL, GDDmat_M = NULL, GDDmat_S = NULL, GDDmat_W = NULL,
                            hybgdd, GDD_baseT, GDD_maxIncrease, max_growing_season_length, leaf_longevity){
  # crop phenology in CLM4.5
  # 1. determining planting date based on climate (GDDx20) and weather, it can be prescribed
  # 2. determining GDD requirement for emergence (GDD_emer) based on GDD_mat, it can be prescribed
  # 3. determining GDD requirement for grain fill (GDD_repr) based on GDD_mat, it can be prescribed
  # 4. determining GDD requirement for reaching physical maturity (a.k.a. harvest date) based on climate (GDDx20), it can be prescribed
  
  if (leap) {day_per_year = 366}
  else {day_per_year = 365}
  
  if (at_NH_flag){
    prescribed_min_plant_day = c(92,92,122)
    prescribed_max_plant_day = c(166,166,166)
    crop_calendar_start_jday = 1
  } else {
    prescribed_min_plant_day = c(275,275,306)
    prescribed_max_plant_day = c(349,349,349)
    crop_calendar_start_jday = 183
  }
  
  if (leap){
    prescribed_min_plant_day = prescribed_min_plant_day + 1
    prescribed_max_plant_day = prescribed_max_plant_day + 1
    if(!at_NH_flag){
      crop_calender_start_jday + 1
    }
  }
  
  current_jday = date.to.day(current_date)
  GDD_update = f_crop_GDD_update(GDD_T_2m_prev = GDD_T2m, GDD_T_soil_prev = GDD_Tsoil, today_T_soil = T_soil, today_T_2m = T2m, GDD_base_T = GDD_baseT, GDD_max_increase = GDD_maxIncrease, crop_live_flag = crop_living_flag)
  GDD_T2m = GDD_update$GDD_T_2m
  GDD_Tsoil = GDD_update$GDD_T_soil
  
  # initializing variables for potential growing season
  if (current_jday == crop_calendar_start_jday) {
    if (!crop_living_flag) {
      crop_planting_flag = FALSE; leaf_emer_flag = FALSE; grain_filling_flag = FALSE; harvesting_flag = FALSE
      planting_jday = NA; harvest_jday = NA
      LAI = 0; SAI = 0; leafC = 0; livestemC; grainC = 0; finerootC = 0; 
    } else {
      warning("[f_crop_phenology] crop_living_flag is TRUE when a potential growing season starts, needs debugging")
    }
  }
  
  # determining planting date
  if (!crop_living_flag && !crop_planting_flag){
    
    if (any(get_planting_date_option == c('prescribed-map','prescribed-site'))){
      if (current_jday == if(is.na(prescribed_planting_date)) prescribed_planting_date_readin else prescribed_planting_date){ # using prescibed planting date
        crop_planting_flag = T; planting_jday = current_jday 
        seed_C_to_leaf_C = emergence_carbon_to_leaf
        print(paste0("[Crop model] Crop is planted on ", current_date))
      }
    } else if (get_planting_date_option == 'CLM4.5') { # using climate and weather condition to predict planting date according to CLM4.5 crop scheme
      # 4 conditions:
      # 1. 10-days running mean daily mean T_2m > T_plant_req  (weather constraint)
      # 2. 10-days running mean daily min. T_2m > T_min_plant_req (weather constraint)
      # 3. GDD8_20yr != 0 dday (climate constraint),this is modified a bit from the original CLM4.5 which is GDD8_20yr >= 50 dday
      # 4. Between 'prescribed_min_plant_day' and 'prescribed_max_plant_day'
      if (T_10_d > T_plant_req && T_min_10_d > T_min_plant_req &&
          GDD8_20yr > 0 &&
          current_jday >= prescribed_min_plant_day[icrop] && current_jday <= prescribed_max_plant_day[icrop]){
          crop_planting_flag = T; planting_jday = current_jday 
        seed_C_to_leaf_C = emergence_carbon_to_leaf
        print(paste0("[Crop model] Crop is planted on ", current_date))
      }
    } # else (!prescribed_planting_date_flag)
    
    # determining GDD_emer, GDD_repr, GDD_mat
    if (crop_living_flag){
      
        if (get_GDDmat_method == "custom"){
            GDD_mat = prescribed_GDD_mat
        } else if (get_GDDmat_method == "Sack") {
            if (ipft == 18 || ipft == 19){GDDmat = GDDmat_M}
            if (ipft == 20 || ipft == 21){GDDmat = GDDmat_W}
            if (ipft == 24 || ipft == 25){GDDmat = GDDmat_S}
        } else if (get_GDDmat_method == "CLM4.5"){
            if (ipft == 18 || ipft == 19){GDD_mat = max(950, min(GDD8_20yr * 0.85, hybgdd))}
            if (ipft == 20 || ipft == 21){GDD_mat = min(GDD0_20yr, hybgdd)}
            if (ipft == 24 || ipft == 25){GDD_mat = min(GDD10_20yr, hybgdd)}
        }
        
        if (get_GDDrepr_method == "custom") {
            GDD_repr = prescribed_GDD_repr
        } else if (get_GDDrepr_method == "CLM4.5") {
            if (ipft == 18 || ipft == 19){
                maize_maturity_rating = max(73, min(135, (GDD_mat + 53.683)/13.882))          # details of this formula can be found in Kucharik (2003), Earth Interactions No.7
                GDD_repr_factor = -0.002 * (maize_maturity_rating - 73) + 0.65
                GDD_repr_factor = min(max(GDD_repr_factor, 0.55), 0.65)
                GDD_repr = GDD_mat * GDD_repr_factor
            }
            if (ipft == 20 || ipft == 21){GDD_repr = 0.6 * GDD_mat}
            if (ipft == 24 || ipft == 25){GDD_repr = 0.7 * GDD_mat}
        }
        
        if (get_GDDemer_method == "custom") {
            GDD_emer = prescribed_GDD_emer
        } else if (get_GDDrepr_method == "CLM4.5") {
            if (ipft == 18 || ipft == 19){GDD_emer = 0.03 * GDD_mat}
            if (ipft == 20 || ipft == 21){GDD_emer = 0.05 * GDD_mat}
            if (ipft == 24 || ipft == 25){GDD_emer = 0.03 * GDD_mat}
        }
    
      # Checking the value of GDD_mat, GDD_repr, GDD_emer
      if (GDD_mat <= GDD_repr || GDD_mat <= GDD_emer || GDD_repr <= GDD_emer){
        warning("The values of GDD_mat or/and GDD_repr or/and GDD_emer have some problems, make sure they follow the order 'GDD_mat > GDD_repr > GDD_emer'")
      }
    } else {
    # crop_living_flag == F, GDD requirments will be determined once the crop is planted
    }
    
  } else {
    # crop is planted already, skipping all preocess calculating GDD requirement
  }
  
  # determining emergence, grain fill and harvesting
  if (crop_living_flag) {
    
    if (current_jday >= planting_jday) { # find the day since planting for determining harvest date
      days_since_planting = current_jday - planting_jday
    } else {
        days_since_planting = current_jday - planting_jday + day_per_year   # for growing season that cross a year (in Southern Hemisphere)
    }
    
    # harvesting conditions
    # 1. GDDT2m > GDDmat   OR
    # 2. growing season is too long that crops are forced to be harvested 
    if (GDD_T_2m >= GDD_mat || days_since_planting > max_growing_season_length) {   
      crop_living_flag = F
      grain_filling_flag = F; leaf_emer_flag = F; harvesting_flag = T
      harvest_jday = current_jday
      
      # at harvest, remove C from all carbon pools
      seed_C_to_leaf_C = 0 
      leaf_C_loss_flux = leafC / 86400
      grain_C_loss_flux = grainC / 86400
      livestem_C_loss_flux = livestemC / 86400
      fineroot_C_loss_flux = finerootC / 86400
      
    } else if (GDD_T_2m >= GDD_repr && GDD_T_2m < GDD_mat) {    # at reproductive stage
      if (!grain_filling_flag){ #first time enter this condition, turn off leaf_emer_flag
        grain_filling_flag = T
        leaf_emer_flag = F
      }
      # leaf senescence for every timestep during grain fill, unit: gC m^-2 s^-1, leaf_longevity = leaf longevity (yr) from PFT_surf_data.R
      leaf_C_loss_flux = 1/(leaf_longevity * 86400 * 365) * leafC  
      seed_C_to_leaf_C_flux = 0; grain_C_loss_flux = 0; livestem_C_loss_flux = 0; fineroot_C_loss_flux = 0
      # CHECK apart from leaf, is there other loss fluxes???
    } else if (GDD_T_soil >= GDD_emer && GDD_T_2m < GDD_repr){   # at vegetative stage
      if(!leaf_emer_flag){ #first time enter this condition, trigger carbon transfer to leaf from seed in this timestep only
        leaf_emer_flag = T
        seed_C_to_leaf_C_flux = seed_C_to_leaf_C / 86400  # this flux has unit gC m^-2 s^-1
        leaf_C_loss_flux = 0; grain_C_loss_flux = 0; livestem_C_loss_flux = 0; fineroot_C_loss_flux = 0
      } else {
        seed_C_to_leaf_C_flux = 0; leaf_C_loss_flux = 0; grain_C_loss_flux = 0; livestem_C_loss_flux = 0; fineroot_C_loss_flux = 0
      }
    } else {
      # between planting and emergence
      # declare variables for returning
      seed_C_to_leaf_C_flux = 0; leaf_C_loss_flux = 0; grain_C_loss_flux = 0; livestem_C_loss_flux = 0; fineroot_C_loss_flux = 0
    } 
      
  } else {
      # crop_living_flag == F, skipping processes to check for emergence, grain fill and harvesting
      seed_C_to_leaf_C_flux = 0; leaf_C_loss_flux = 0; grain_C_loss_flux = 0; livestem_C_loss_flux = 0; fineroot_C_loss_flux = 0
  }
  
  #output
  #1. fluxes: seedC_to_leafC, leafC_loss, grainC_loss, livestemC_loss, finerootC_loss
  #2. phenology flags: crop_living_flag, crop_planting_flag, leaf_emer_flag, grain_filling_flag, harvest_flag
  #3. phenology date:  planting_jday, harvest_jday
  #4. phenology GDD: GDD_T_2m, GDD_T_soil
  #5. C pools and other physiology variables: leaf_C, fineroot_C, grain_C, livestem_C, LAI, SAI
  output = list(seedC_to_leafC_flux = seed_C_to_leaf_C_flux, leafC_loss_flux = leaf_C_loss_flux, grainC_loss_flux = grainC_loss_flux, livestemC_loss_flux = livestem_C_loss_flux, finerootC_loss_flux = fineroot_C_loss_flux,
                croplive_flag = crop_living_flag, cropplant_flag = crop_planting_flag, leafemergence_flag = leaf_emer_flag, grainfill_flag = grain_filling_flag, har_flag = harvesting_flag,
                planting_julianday = planting_jday, harvesting_julianday = harvest_jday,
                GDD_T2m = GDD_T2m, GDD_Tsoil = GDD_Tsoil, GDDmat = GDDmat, GDDrepr = GDDrepr, GDDemer = GDDemer,
                leafC = leafC, grainC = grainC, finerootC = finerootC, livestemC = livestemC, LAI_out = LAI, SAI_out = SAI)
  
  return(output)
  
} # end of f_crop_phenology


f_evergreen_phenology = function(){
  # evergreen phenology in CLM4.5
  # d(leaf_C)/dt = -leaf_long * leaf_C
  
}

f_stress_deciduous_phenology = function(){
  
}

f_seasonal_deciduous_phenology = function(){
  
}



f_onset_growth_fluxes = function(){
  
}

f_litter_fall_fluxes = function(){
  
}

f_wood_turnover_fluxes = function(){
  
}

f_litter_to_soil_fluxes = function(){
  
}

f_crop_GDD_update = function(GDD_T_2m_prev, GDD_T_soil_prev, today_T_soil, today_T_2m, GDD_base_T, GDD_max_increase, crop_live_flag){
  T_soil_degC = today_T_soil - 273.15
  T_2m_degC = today_T_2m - 273.15
  
  if (!user_defined_crop_GDD_accmulation_flag){
    GDD_soil_increase = min(max(T_soil_degC - GDD_base_T,0), GDD_max_increase)
    GDD_T_2m_increase = min(max(T_2m_degC - GDD_base_T,0), GDD_max_increase)
  } else {
    #user-defined GDD
    # GDD_T_2m_increase = ...
    # GDD_soil_increase = ...
  }
  
  # GDD accummulates from planting until harvesting
  if (crop_live_flag){
    GDD_T_soil_output = GDD_T_soil_prev + GDD_soil_increase
    GDD_T_2m_output = GDD_T_2m_prev + GDD_T_2m_increase
  } else {
    GDD_T_2m_output = 0
    GDD_T_soil_output = 0
  }
  
  output = list(GDD_T_2m = GDD_T_2m_output, GDD_T_soil = GDD_T_soil_output)
  
}
