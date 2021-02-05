###############################################################################
### Module for the main simulation function to loop over lon and lat
### Including a function to reshape history data into gridded format
###############################################################################

# Both functions require access to variables in the global environment.

################################################################################
### Revision history
################################################################################

# Aug 2017: v1.0 finalized (Tai)
# Oct 2018: Added new "Medlyn" stomatal conductance scheme in "Farquhar_Ball_Berry.R", and thus here necessary variables (e.g., vpd) are added. (Sun)
# Feb 2019: Modified the way L_sun and L_sha are defined and calculated. We believe that A_can, R_can and g_can should be scaled up by LAI only, not by LAI + SAI. Therefore, L_sun and L_sha here should be sunlit and shaded leaf area index, not plant area index, but we still consider both LAI and SAI when calculating light extinction. Now they are renamed "LAI_sun" and "LAI_sha" and "LAI_sun" has to taken explicitly from canopy radiative transfer model. (Tai)
# Feb 2019: Now aerodynamic conductance for heat, g_ah, is always calculated using either one of the two schemes (default from Monin_Obukhov.R vs. GEOS-Chem method in drydep_toolbox.R), and inputted into f_canopy_photosyn(). (Tai)

###############################################################################

# Conduct hourly simulation for one single day, looping over i (global longitude index) and j (global latitude index):

f_simulate_ij = function(IJ) {
   
   # "IJ = c(i, j)" is a 2-element vector containing the global lon ("i") and lat ("j") index.
   # It requires variables and parameters defined in the global environment.
   # When the list of c(i, j) is defined, oceanic and permanently glacial grid cells are already excluded.
   
   # Redefine new history outputs within function:
   for (ivar in seq_along(var_name$variable_name)){
      current_var = var_name$variable_name[ivar]
      if (var_name$res_level[ivar] == 'PFT') {
         assign(x = paste0(current_var, '_PFT_hist_ijd'), value = array(NA, dim=c(length(pftname), 24/dt_hr)))
      } else if (var_name$res_level[ivar] == 'grid') {
         assign(x = paste0(current_var, '_hist_ijd'), value =  array(NA, dim=24/dt_hr))
      } else if (var_name$res_level[ivar] == 'PFT_daily') {
         assign(x = paste0(current_var, '_PFT_hist_ijd'), value = array(NA, dim = length(pftname)))
      }
   }
   rm(current_var, ivar)
   
   # Redefine new history outputs within function if not existent for next simulation step:
   if (O3_damage_flag) {
      for (out_var in c('CUO_sun_PFT_hist_ijd', 'CUO_sha_PFT_hist_ijd', 'CUO_can_PFT_hist_ijd')) {
         assign(x = out_var, value = array(NA, dim=c(length(pftname), 24/dt_hr)))
      }
      rm(out_var)
   }
   
   # Redefine new history output for temp_data for BGC mode
   if (biogeochem_flag) {
       BGC_temp_var_vec = c('LAI_PFT_hist_ijd', 'SAI_PFT_hist_ijd', 'leafC_PFT_hist_ijd', 'finerootC_PFT_hist_ijd', 'livestemC_PFT_hist_ijd', 'deadstemC_PFT_hist_ijd', 'livecoarserootC_PFT_hist_ijd', 'deadcoarserootC_PFT_hist_ijd', 'grainC_PFT_hist_ijd',
                            'GDDT2m_PFT_hist_ijd', 'GDDTsoil_PFT_hist_ijd', 'GDDmat_PFT_hist_ijd', 'GDDemer_PFT_hist_ijd', 'GDDrepr_PFT_hist_ijd',
                            'crop_live_flag_PFT_hist_ijd', 'crop_plant_flag_PFT_hist_ijd', 'leaf_emergence_flag_PFT_hist_ijd', 'grain_fill_flag_PFT_hist_ijd', 'harvest_flag_PFT_hist_ijd', 'peak_LAI_flag_PFT_hist_ijd',
                            'day_of_planting_PFT_hist_ijd', 'day_of_grain_filling_PFT_hist_ijd', 'day_of_harvesting_PFT_hist_ijd', 'astem_PFT_hist_ijd', 'aleaf_PFT_hist_ijd', 'astem_leafem_PFT_hist_ijd', 'aleaf_leafem_PFT_hist_ijd')
       
       for (out_var in BGC_temp_var_vec) {
           if (!exists(out_var)) {
               # All these variables are calculated daily, no hourly dimension
               assign(x = out_var, value = array(NA, dim = c(length(pftname))))
           }
       }
       
       if (O3_damage_flag & O3_damage_scheme == 'Lombardozzi' & O3_sensitivity == 'custom') {
          for (out_var in c('CUO_grainfill_PFT_hist_ijd')) {
             assign(x = out_var, value = array(NA, dim = c(length(pftname))))
          }
       }
   }
   
   #############################################################################
   
   # Indices for lon/lat:
   i = IJ[1]
   j = IJ[2]
   
   if (nchar(i) > 3) stop('Max i >= 1000!')
   i_str = paste0(paste0(rep(0, 3-nchar(i)), collapse = ''), i)
   if (nchar(j) > 3) stop('Max j >= 1000!')
   j_str = paste0(paste0(rep(0, 3-nchar(j)), collapse = ''), j)
   
   # Longitude (rad):
   lon_rad = lon[i]/180*pi
   # Latitude (rad):
   lat_rad = lat[j]/180*pi
   
   # Geopotential height of surface (m): 
   if (FLUXNET_flag){
      site_info = read_excel(path = paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx'))
      site_ind = match(FLUXNET_site_id, site_info$SITE_ID)
      Z_surf = as.double(unname(site_info[site_ind,"LOCATION_ELEV"]))
      rm(site_info, site_ind)
      # if elevation is not provided by FLUXNET then uses geopotential height of surface (m) as replacement
      if (is.na(Z_surf)) Z_surf = PHIS[i,j]/g_E
   } else {
      Z_surf = PHIS[i,j]/g_E
   }
   
   # Daily mean 2-m air and soil temperature (K):
   T_daily = mean(T2M[i,j,], na.rm=TRUE)
   Tmin_daily = min(T2M[i,j,], na.rm=TRUE)
   if (biogeochem_flag) {T_soil1_daily = mean(TSOIL1[i,j,], na.rm = TRUE)}
   
   # 10-day mean air temperature (K) to calculate temperature acclimation:
   acclimation_flag = ((d > 10) | continue_flag)
   if (acclimation_flag) {
      # Daily mean temperature of the past 10 days:
      T_daily_last10d = rep(NaN, times=10)
      T_dailymin_last10d = rep(NaN, times=10)
      for (d_prev in 10:1) {
         # Date of d_prev days before current date:
         previous_date = to.yyyymmdd(from.yyyymmdd(current_date) - d_prev*24)
         filename = paste0(simulation_dir, 'temp_data/temp_', as.character(previous_date), '/temp_i', i_str, '_j', j_str, '.RData')
         load(filename)
         T_daily_last10d[11-d_prev] = T_daily_ijd
         T_dailymin_last10d[11-d_prev] = T_dailymin_ijd
      }
      # 10-day mean air temperature (K):
      T_10d = mean(T_daily_last10d, na.rm=TRUE)
      Tmin_10d = mean(T_dailymin_last10d, na.rm=TRUE)
   } else {
      # 10-day mean air temperature (K):
      T_10d = T_daily
      Tmin_10d = Tmin_daily
   }
   
   # Soil albedo for PAR for dry and saturated soil:
   # 3rd dim of "soil albedo" = [dry (visible), dry (Near IR), saturated (visible), saturated (Near IR)]
   alpha_soil_dry = soil_albedo[i,j,1]
   alpha_soil_sat = soil_albedo[i,j,3]
   
   # Other model parameters:
   met_cond_flag = TRUE
   colimit_flag = TRUE
   
   #############################################################################
   
   success = FALSE
   
   if (biogeochem_flag) {
       selected_pft = BGC_pft_selection
   } else {
       selected_pft = 2:length(pftname)
   }
   
   
   # PFT-specific parameters:
   for (ipft in selected_pft) {
      if (!biogeochem_flag){
          # Prescribed leaf area index (m^2 m^-2):
          if (leap & as.numeric(MM) > 2) n_PAI = n_day_whole - 1 else n_PAI = n_day_whole
          LAI = LAI_day_PFT[i,j,ipft,n_PAI]
      } else {
          if (d == 1){
              # Read LAI and SAI from initial_data for the first day of simulation
              LAI = LAI_initial_map[i,j,ipft]
              SAI = SAI_initial_map[i,j,ipft]
          } else {
              # Read LAI and SAI from temp_data otherwise
              prev_date = to.yyyymmdd(from.yyyymmdd(current_date) - 24)
              # filename = paste0(sim_dir,'temp_data/temp_',prev_date,'/temp_i',i_str,'_j', j_str, '.RData')
              filename = paste0('temp_data/temp_',prev_date,'/temp_i',i_str,'_j', j_str, '.RData')
              print(paste0("Loading ", filename))
              load(filename)
              
              LAI = tlai_PFT_lasthd_ijd[ipft]
              SAI = tsai_PFT_lasthd_ijd[ipft]
          }
      }

      if ((LAI < 0.01 && !biogeochem_flag) | PFT_frac[i,j,ipft] < 0.01) {
         # Too little vegetation. Skip current PFT calculations.
          print(paste0("Skip calculation for ipft = ", ipft, " - ", pftname[ipft]))
         next
         
      } else {
         if (!biogeochem_flag) {
             
             # Leaf area index range (m^2 m^-2) (shsun)
             LAI_min = min(LAI_day_PFT[i,j,ipft,], na.rm = TRUE)
             LAI_max = max(LAI_day_PFT[i,j,ipft,], na.rm = TRUE)
             
         # Stem area index (m^2 m^-2):
         SAI = SAI_day_PFT[i,j,ipft,n_PAI]
         # Consider LAI only?
         # SAI = 0
         
         # Leaf area index of previous day (m^2 m^-2):
         if (n_PAI == 1) LAI_prev_day = LAI_day_PFT[i,j,ipft,365] else LAI_prev_day = LAI_day_PFT[i,j,ipft,n_PAI-1]
         }
         
         # PFT-level displacement height and roughness length (not used in the current version):
         # Canopy height (m):
         # Needed for dry deposition calculation?
         z_can = htop[ipft]
         # z_can = 5
         
         # # Fractional weight (0-1):
         # V_fw = (1 - exp(-min(c((LAI + SAI), 2), na.rm=TRUE)))/(1 - exp(-2))
         # # PFT-level displacement height (m):
         # z_disp = htop[ipft]*displar[ipft]*V_fw
         # # PFT-level roughness length for momentum (m):
         # z_0m = exp(V_fw*log(htop[ipft]*z0mr[ipft]) + (1 - V_fw)*log(0.01))
         # # 0.01 is z_0m for soil, glacier and wetland. For snow-covered surface it should be 0.0024.
         # # Roughness lengths for heat and water vapor are assumed the same.
         
         #######################################################################
         
         # Photosynthetic and canopy radiative parameters:
         # C3 or C4 plant?
         C3_plant_flag = (c3psn[ipft] == 1)
         # Plant group:
         plant_group = plantgroup[ipft]
         # Evergreen?
         evergreen_flag = (evergreen[ipft] == 1)
         # Characteristic leaf dimension (m):
         d_leaf = dleaf[ipft]
         # Maximum rate of carboxylation at top of canopy at 25 degC (umol CO2 m^-2 s^-1):
         V_cmax25_0 = vcmax25top[ipft]
         # Leaf nitrogen concentration at top of canopy (g N m^-2 leaf area):
         leaf_N_conc = lnctop[ipft]
         # Soil matric potential when stomata are fully closed (mm):
         psi_c = smpsc[ipft]
         # Soil matric potential when stomata are fully opened (mm):
         psi_o = smpso[ipft]
         # Leaf longevity (yr):
         leaflong = leaf_long[ipft]
         # Leaf orientation:
         x_l = xl[ipft]
         # Leaf reflectance for PAR:
         alpha_leaf = rholvis[ipft]
         # Stem reflectance for PAR:
         alpha_stem = rhosvis[ipft]
         # Leaf transmittance for PAR:
         tau_leaf = taulvis[ipft]
         # Stem transmittance for PAR:
         tau_stem = tausvis[ipft]
         
         if (biogeochem_flag) {
             
             # Flags for differentiating the parameters for crops planting season and phenology
             at_NH_flag = ifelse(test = lat[j] >= 0, yes = TRUE, no = FALSE)
             at_tropical_flag = ifelse(test = lat[j] < 30 & lat[j] > -30, yes = TRUE, no = FALSE)

             # Note:
             # C3 unmanaged rainfed/irrigated crops (ipft == 16:17) are belong to both stress deciduous and crop PFT, but the phenology is based on stress deciduous PFT 
             
             # Leaf allocation coefficient parameter for crops
             a_leaf_final = aleaff[ipft]
             a_leaf_alloc_power = allconsl[ipft]
             b_factor = bfact[ipft]
             a_leaf_base = fleafi[ipft]
             # Stem allocation coefficient parameter for crops
             a_stem_final = astemf[ipft]
             a_stem_alloc_power = allconss[ipft]
             GDD_decline_factor = declfact[ipft]
             # Root allocation coefficient parameter for crops
             a_root_initial = arooti[ipft]
             a_root_final = arootf[ipft]
             # Base temperature og GDD accumulation for crops
             # Tropical maize/soybean: GDD_base_T = 10; 'Tropical' wheat: GDD_base_T = 12 - 0.4 * abs(lat)
             GDD_base_T = baset[ipft]; if (at_tropical_flag && any(ipft == c(18,19,24,25))) {GDD_base_T = 10}; if (at_tropical_flag && any(ipft == c(20,21))) {GDD_base_T = 12 - 0.4*abs(lat[j])}
             # Allocation ratio of coarse root : live stem
             allocRatio_croot.stem = croot_stem[ipft]
             # Crop (0 = non-crop; 1 = crop):
             crop_flag = crop[ipft]
             # C:N ratio of dead wood
             deadwd_cn = deadwdcn[ipft]
             # Through canopy (projected area basis) dSLA/dLAI (m^2 gC^-1):
             dsla_dlai = dsladlai[ipft]
             # Allocation flag to storage pools (deciduous) or display pools (evergreen and crops)
             f_cur = fcur[ipft]
             # C:N ratio of fine root during reproductive stage for crops
             froot_cn_final = ffrootcn[ipft]
             # C:N ratio of leaf during reproductive stage for crops
             leaf_cn_final = fleafcn[ipft]
             # Allocation flag to wood
             alloc_livewood = flivewd[ipft]
             # C:N ratio of fine root
             froot_cn = frootcn[ipft]
             # Allocation ratio of fine root : leaf
             allocRatio_froot.leaf = froot_leaf[ipft]
             # C:N ratio of live stem during reproductive stage
             stem_cn_final = fstemcn[ipft]
             # C:N ratio of grain
             grain_cn = graincn[ipft]
             # % of GDD_mat to reach reproductive stage
             # Note that this parameters for temperate soybean is changed from 0.7 in CLM4.5 to 0.5 in CLM5 according to the technical note. Here we adopt the CLM5 paramter
             repr_GDDfraction = grnfill[ipft]; if (at_tropical_flag && any(ipft == c(18,19,24,25))) {repr_GDDfraction = 0.5}
             # Maximum GDD_mat allowed for crops
             GDDmat_max = hybgdd[ipft]; if (at_tropical_flag && any(ipft == c(24,25))) {GDDmat_max = 2100}
             # Prescribed maximum LAI allowed for crops
             LAI_max = laimx[ipft]
             # C:N ratio of leaf
             leaf_cn = leafcn[ipft]
             # % of GDD_mat to reach vegetative stage
             emer_GDDfraction = lfemerg[ipft]
             # C:N ratio of leaf litter
             leafLitter_cn = lflitcn[ipft]
             # C:N ratio of live wood
             liveWood_cn = livewdcn[ipft]
             # Daily min planting temperature requirment (K)
             min_T_planting_req = min_planting_temp[ipft]; if (at_tropical_flag && any(ipft == c(18,19,24,25))) {min_T_planting_req = 283.15}
             # Maximum growing season length allowed
             crop_season_length_max = mxmat[ipft]; if (at_tropical_flag && any(ipft == c(18,19))) {crop_season_length_max = 160}
             # Maximum increase in GDD_T2m allowed
             GDDT2m_change_max = mxtmp[ipft]
             # Daily average planting temperature requirement (K)
             avg_T_planting_req = planting_temp[ipft]; if (at_tropical_flag && any(ipft == c(18,19,24,25))) {avg_T_planting_req = 294.15}
             # CLM rooting distribution parameter a and b and the root fraction
             root_a_par = roota_par[ipft]
             root_b_par = roota_par[ipft]
             get_root_fraction = f_get_root_fraction(soil_depth_array = soil_layer_depth, a_rootfrac = root_a_par, b_rootfrac = root_b_par)
             root_frac_array = get_root_fraction$root_fraction_array
             # Seasonal deciduous PFT:
             season_decid_flag = (season_decid[ipft] == 1)
             # Specific leaf area (SLA) at top of canopy (projected area basis) (m^2 gC^-1):
             sla_top = slatop[ipft]
             # Allocation ratio of live stem : leaf
             allocRatio_stem.leaf = stem_leaf[ipft]
             # Stress deciduous PFT:
             stress_decid_flag = (stress_decid[ipft] == 1)
             # Woody life form (0 = non-woody; 1 = woody):
             woody_flag = (woody[ipft] == 1)
             # Maximum canopy height for crops:
             ztopmax = ztopmx[ipft]; if (at_tropical_flag && any(ipft == c(24,25))) {ztopmax = 1}
             # CLM4.5 planting period for crops
             # Not sure if the phenology in the SH tropical region in CLM5 is also 6 month later than the NH counterpart. Here we assume yes at the time being.
             earliest_planting_doy_NH = earliest_planting_jday_possible_NH[ipft]; if (at_tropical_flag) {if (any(ipft == c(18,19))) {earliest_planting_doy_NH = 79} else if (any(ipft == c(24,25))) {earliest_planting_doy_NH = 105}}
             latest_planting_doy_NH = latest_planting_jday_possible_NH[ipft]; if (at_tropical_flag) {if (any(ipft == c(18,19))) {latest_planting_doy_NH = 105} else if (any(ipft == c(24,25))) {latest_planting_doy_NH = 181}}
             earliest_planting_doy_SH = earliest_planting_jday_possible_SH[ipft]; if (at_tropical_flag) {if (any(ipft == c(18,19))) {earliest_planting_doy_SH = 263} else if (any(ipft == c(24,25))) {earliest_planting_doy_SH = 288}}
             latest_planting_doy_SH = latest_planting_jday_possible_SH[ipft]; if (at_tropical_flag) {if (any(ipft == c(18,19))) {latest_planting_doy_SH = 288} else if (any(ipft == c(24,25))) {latest_planting_doy_SH = 365}}
         }
         
         # [Crop model] read in site and PFT planting date (POD or BGC), harvesting date (POD only) and GDDmat (BGC only)
         if (biogeochem_flag || O3_POD) {
             if (O3_POD || (biogeochem_flag && get_planting_date_option == 'prescribed-map')) {
                 # crop name in the Sack data set
                 # maize (primary growing season), soybean, wheat, winter wheat, rice (primary growing season), rice (secondary growing season), maize (secondary growing season)
                 # The regridded numbers may have decimal, round them to the nearest whole number (Pang, Jun 2019)
                 if (any(ipft == c(18,19))) {
                     if (crop_growing_season == 'primary') {
                         prescribed_planting_date = round(prescribed_planting_date_Sack[i,j,1], digits = 0)
                         print(paste0('[Sack planting date] = ', prescribed_planting_date))
                     } else if (crop_growing_season == 'secondary') {
                         # prescribed_planting_date = round(prescribed_planting_date_Sack[i,j,7], digits = 0)
                     } else {stop("'crop_growing_season' for Sack dataset is not correctly specified")}
                 } else if (any(ipft == c(20,21))) {
                     prescribed_planting_date = round(prescribed_planting_date_Sack[i,j,2], digits = 0)
                 } else if (any(ipft == c(22,23))) {
                     prescribed_planting_date = round(prescribed_planting_date_Sack[i,j,3], digits = 0)
                 } else if (any(ipft == c(24,25))) {
                     prescribed_planting_date = round(prescribed_planting_date_Sack[i,j,4], digits = 0)
                 } else {
                     prescribed_planting_date = NA
                 }
             }
             
             if (O3_POD && !biogeochem_flag) {
                 # read in harvesting date only
                 # ... = prescribed_harvesting_date_Sack[i,j,x]
             }
             
             if (!O3_POD && biogeochem_flag) {
                 # read in GDDmat
                 if (get_GDDmat_method == 'CLM4.5') {
                     GDD0 = GDD0_map[i,j]
                     GDD8 = GDD8_map[i,j]
                     GDD10 = GDD10_map[i,j]
                 } else if (get_GDDmat_method == 'Sack') {
                     GDDmat_maize = GDDmat_maize_map[i,j]
                     GDDmat_springwheat = GDDmat_springwheat_map[i,j]
                     GDDmat_winterwheat = GDDmat_winterwheat_map[i,j]
                     GDDmat_soybean = GDDmat_soybean_map[i,j]
                 }
             }
         }
         
         # Soil parameters:
         # Saturated soil matric potential for bulk root zone (mm):
         psi_sat = psi_sat_PFT[i,j,ipft]
         # Saturated volumetric water content for bulk root zone (fraction):
         theta_sat = theta_sat_PFT[i,j,ipft]
         # Clapp and Homberger parameter for bulk root zone:
         b_psi = b_psi_PFT[i,j,ipft]
         
         #######################################################################
         
         for (h in 1:(24/dt_hr)) {
            # Number of days from 00:00 UTC Jan 1 including hours:
            # "0.5" because hourly met fields are 1-hour average, not instantaneous.
            n_day = n_day_whole - 1 + (h - 0.5)/24
            # Solar declination angle (rad):
            decl = f_decl(n_day=n_day)
            # Cosine of solar zenith angle:
            cos_SZA = f_cosSZA(lat=lat_rad, lon=lon_rad, n_day=n_day)
            
            # Incoming radiation:
            # Incident direct beam PAR (W m^-2):
            PAR_beam = PARDR[i,j,h]
            # Incident diffuse beam PAR (W m^-2):
            PAR_diff = PARDF[i,j,h]
            # Incoming shortwave radiation (W m^-2):
            swr = SWGDN[i,j,h] 
            
            ####################################################################
            
            # Micrometeorological variables:
            
            # Sensible heat flux (W m^-2):
            H_sen = HFLUX[i,j,h]
            # Evapotranspiration flux (kg m^-2 s^-1):
            ET = EVAP[i,j,h]
            # Characteristic velocity scale or friction velocity (m s^-1):
            u_star = USTAR[i,j,h]
            # Grid-level surface roughness for momentum (m):
            Z_0m = Z0M[i,j,h]
            # Grid-level displacement height (m):
            Z_disp = Z_0m*disp_on_z0m[i,j]
            # The zero-plane displacement height where wind speed is extrapolated to zero is Z_disp + Z_0m. In what follows, "displacement height" means the zero-plane displacement height.
            
            # Sea-level pressure (Pa):
            slp = SLP[i,j,h]
            # Atmospheric temperature at 2 m above displacement height (K):
            T_2m = T2M[i,j,h]
            # Atmospheric temperature at 10 m above displacement height (K):
            T_10m = T10M[i,j,h]
            # Specific humidity at 2 m above displacement height (kg kg^-1):
            q_2m = QV2M[i,j,h]
            # # Specific humidity at 10 m above displacement height (kg kg^-1): # Not included in GEOS-Chem met fields
            # q_10m = QV10M[i,j,h]
            # Wind speed at 10 m above displacement height (m s^-1):
            u_10m = if (FLUXNET_flag) WS[i,j,h] else sqrt(U10M[i,j,h]^2 + V10M[i,j,h]^2)
            print(paste0('h = ', h, ' u_10m = ', signif(WS[i,j,h], digits = 3)))
            
            # Variables below are needed for dry deposition (Sun, Oct 2018):
            # Liquid Precipitation (kg m-2 s-1): 
            prec_liq = if (FLUXNET_flag) PRECTOT[i,j,h] else PRECTOT[i,j,h] - PRECSNO[i,j,h] 
            # Snow depth (m):
            if (FLUXNET_flag) d_snow = 0 else d_snow = SNODP[i,j,h]
            # Latent heat (W m^-2):
            L_latent = EFLUX[i,j,h]
            
            # Atmospheric scale height (m):
            Z_scale = R_da*(T_2m + 0.0065*(Z_surf + Z_disp + Z_0m + 2)/2)/g_E
            # Surface pressure (Pa):
            P_surf = if (FLUXNET_flag) ATMP[i,j,h] else slp*exp(-Z_surf/Z_scale)
            # P_surf = slp*(1 - 0.0065*Z_surf/(T_2m + 0.0065*(Z_surf + Z_disp + Z_0m + 2)))^5.257
            # Pressure at zero-plane displacement height (d + z0m) where wind speed is extrapolated to zero (Pa):
            P_disp = if (FLUXNET_flag) ATMP[i,j,h] else slp*exp(-(Z_surf + Z_disp + Z_0m)/Z_scale)
            print(paste0('h = ', h, ' P_disp or P_surf  = ', signif(ATMP[i,j,h], digits = 3)))
            
            
            # At 2 m above displacement height:
            # Define T = theta at the surface.
            # Atmospheric pressure (Pa):
            P_2m = if (FLUXNET_flag) ATMP[i,j,h] else slp*exp(-(Z_surf + Z_disp + Z_0m + 2)/Z_scale)
            # Atmospheric potential temperature (K):
            theta_2m = if (FLUXNET_flag) T_2m + (g_E/c_p)*(Z_disp + Z_0m + 2) else T_2m*(P_surf/P_2m)^(R_da/c_p)
            print(paste0('h = ', h, ' theta_2m = ', signif(theta_2m, digits = 3)))
            
            # Vapor pressure (Pa):
            e_2m = P_2m*q_2m/(0.622 + 0.378*q_2m)
            # Moist air density (kg m^-3):
            rho_2m = (P_2m - 0.378*e_2m)/(R_da*T_2m)
            
            # At 10 m above displacement height:
            # Define T = theta at the surface.
            # Atmospheric pressure (Pa):
            P_10m = if (FLUXNET_flag) ATMP[i,j,h] else slp*exp(-(Z_surf + Z_disp + Z_0m + 10)/Z_scale)
            print(paste0('h = ', h, ' P_10m = ', signif(P_10m, digits = 3)))
            
            # Atmospheric potential temperature (K):
            # Define T = theta at the surface.
            theta_10m = if (FLUXNET_flag) T_10m + (g_E/c_p)*(Z_disp + Z_0m + 10) else T_10m*(P_surf/P_10m)^(R_da/c_p)
            print(paste0('h = ', h, ' theta_10m = ', signif(theta_10m, digits = 3)))
            # If "q_10m" is not provided, needed to scale it from the temperature difference (theta_10m - theta_2m).
            # Specific humidity (kg kg^-1):
            q_10m = if (!exists('QV10M')) q_2m + (theta_10m - theta_2m)*c_p*ET/H_sen
            # Vapor pressure (Pa):
            e_10m = P_10m*q_10m/(0.622 + 0.378*q_10m)
            # Moist air density (kg m^-3):
            rho_10m = (P_10m - 0.378*e_10m)/(R_da*T_10m)
            print(paste0('[IJ] T_10m = ', signif(T_10m, digits = 3), 'P_10m = ', signif(P_10m, digits = 3), 'e_10m = ', signif(e_10m, digits = 3)))
            
            # # Impose that "u_10m" cannot be smaller than "u_star" or 1:
            # if (u_10m < u_star) u_10m = u_star
            # if (u_10m < 1) u_10m = 1
            
            # Use a reference height (above zero-plane displacement height, which already includes roughness length here):
            Z_atm = 10
            # Atmospheric pressure at Z_atm (Pa):
            P_atm = P_10m
            # Atmospheric temperature at Z_atm (K):
            T_atm = T_10m
            # Atmospheric potential temperature at Z_atm (K):
            theta_atm = theta_10m
            # Wind speed at Z_atm (m s^-1):
            u_atm = u_10m
            # Moist air density (kg m^-3):
            rho_atm = rho_10m
            # Conductance for water vapor:
            # Not needed anywhere now (Tai, Feb 2019):
            # g_aw = if (Monin_Obukhov_flag) Monin_Obukhov$g_aw else NULL
            # Conversion from molar to meteorological conductance:
            mol_to_met = 1e-6*R_uni*theta_atm/P_atm
            
            ####################################################################
            
            # Biogeochemistry module initialization (Pang)
            # In the first simulation step, decalre new variables and read initial conditions from initial_map.nc
            # In other time steps, read yesterday simulation result from temp_data/xxx.RData
            
            if (biogeochem_flag){
                if (current_date == start_date && h == 1){
                    
                    ###### Initialization by declaring variables
                    BGC_variable_declaration_names = c("a_stem","a_leaf","a_stem_leafem","a_leaf_leafem","GDDT2m","GDDTsoil","GDD_mat","GDD_repr","GDD_emer",
                                                       "crop_live_flag", "peak_LAI_flag", "crop_plant_flag","leaf_emergence_flag","grain_fill_flag","harvest_flag","day_of_planting","day_of_grain_filling","day_of_harvest",
                                                       "CUO_grainfill")
                    BGC_variable_declaration_value = c(rep(0,9),
                                                       rep(FALSE,6), rep(NA,4))
                    
                    if (length(BGC_variable_declaration_names) != length(BGC_variable_declaration_value)) {stop("[simulate_ij.R] length(BGC_variable_declaration_names) should be the same as length(BGC_variable_declaration_value)")}
                    for (z in seq(BGC_variable_declaration_names)) {assign(x = BGC_variable_declaration_names[z], value = BGC_variable_declaration_value[z])}
                    
                    ###### Turning flags back to logical values
                    BGC_flagvariable_names = c("crop_live_flag","peak_LAI_flag","crop_plant_flag","leaf_emergence_flag","grain_fill_flag","harvest_flag")
                    for (z in seq(BGC_flagvariable_names)) {assign(x = BGC_flagvariable_names[z], value = as.logical(eval(parse(text = BGC_flagvariable_names[z]))))}
                    
                    ###### Initialization by reading from initial_map.nc
                    BGC_variable_readin_names = c("leaf_C","fineroot_C","livestem_C","deadstem_C","livecoarseroot_C","deadcoarseroot_C","grain_C", "LAI", "SAI")
                    BGC_variable_readin_value = c(leafC_initial_map[i,j,ipft], frootC_initial_map[i,j,ipft], livestemC_initial_map[i,j,ipft], deadstemC_initial_map[i,j,ipft], livecrootC_initial_map[i,j,ipft], deadcrootC_initial_map[i,j,ipft], grainC_initial_map[i,j,ipft], LAI_initial_map[i,j,ipft], SAI_initial_map[i,j,ipft])
                    if (length(BGC_variable_readin_names) != length(BGC_variable_readin_value)) {stop("[simulate_ij.R] length(BGC_variable_readin_names) should be the same as length(BGC_variable_readin_value)")}
                    for (z in seq(BGC_variable_readin_names)) {assign(x = BGC_variable_readin_names[z], value = BGC_variable_readin_value[z])}
                    ### Turn NA LAI and SAI values to zero
                    LAI = ifelse(test = is.na(LAI), yes = 0, no = LAI); SAI = ifelse(test = is.na(SAI), yes = 0, no = SAI)
                    
                    rm(z, BGC_variable_declaration_names, BGC_variable_declaration_value, BGC_flagvariable_names, BGC_variable_readin_names, BGC_variable_readin_value)
                    
                } else if (h == 1) {
                    ### LAI_prev for the second day is read from initial_data in BGC mode, otherwise it is read from temp_data from the day before yesterday
                    if (current_date != to.yyyymmdd(from.yyyymmdd(start_date) + 24)){
                        date_before_yesterday = to.yyyymmdd(from.yyyymmdd(current_date) - 48)
                        # filename = paste0(sim_dir,'temp_data/temp_',date_before_yesterday,'/temp_i',i_str,'_j', j_str, '.RData')
                        filename = paste0('temp_data/temp_',date_before_yesterday,'/temp_i',i_str,'_j', j_str, '.RData')
                        load(filename)
                        LAI_prev_day = tlai_PFT_lasthd_ijd[ipft]
                    } else {
                        LAI_prev_day = LAI_initial_map[i,j,ipft]
                    }
                    
                    
                    prev_date = to.yyyymmdd(from.yyyymmdd(current_date) - 24)
                    # filename = paste0(sim_dir,'temp_data/temp_',prev_date,'/temp_i',i_str,'_j', j_str, '.RData')
                    filename = paste0('temp_data/temp_',prev_date,'/temp_i',i_str,'_j', j_str, '.RData')
                    print(paste0("Loading ", filename))
                    load(filename)
                    
                    ### Name in BGC_variable_readin_value should be consistent to the output temp_data 
                    BGC_variable_readin_names = c("a_stem","a_leaf","a_stem_leafem","a_leaf_leafem","GDDT2m","GDDTsoil","GDD_mat","GDD_repr","GDD_emer",
                                                  "crop_live_flag", "peak_LAI_flag", "crop_plant_flag","leaf_emergence_flag","grain_fill_flag","harvest_flag","day_of_planting","day_of_grain_filling","day_of_harvest",
                                                  "leaf_C","fineroot_C","livestem_C","deadstem_C","livecoarseroot_C","deadcoarseroot_C","grain_C", "LAI", "SAI")
                    BGC_variable_readin_value = c(astem_PFT_lasthd_ijd[ipft], aleaf_PFT_lasthd_ijd[ipft], astem_leafem_PFT_lasthd_ijd[ipft], aleaf_leafem_PFT_lasthd_ijd[ipft], GDDT2m_PFT_lasthd_ijd[ipft], GDDTsoil_PFT_lasthd_ijd[ipft], GDDmat_PFT_lasthd_ijd[ipft], GDDrepr_PFT_lasthd_ijd[ipft], GDDemer_PFT_lasthd_ijd[ipft],
                                                  crop_live_flag_PFT_lasthd_ijd[ipft], peak_LAI_flag_lasthd_ijd[ipft], crop_plant_flag_PFT_lasthd_ijd[ipft], leaf_emergence_flag_lasthd_ijd[ipft], grain_fill_flag_lasthd_ijd[ipft], harvest_flag_lasthd_ijd[ipft], planting_jday_lasthd_ijd[ipft], grain_filling_jday_lasthd_ijd[ipft], harvesting_jday_lasthd_ijd[ipft],
                                                  leafC_PFT_lasthd_ijd[ipft], frootC_PFT_lasthd_ijd[ipft], livestemC_PFT_lasthd_ijd[ipft], deadstemC_PFT_lasthd_ijd[ipft], livecrootC_PFT_lasthd_ijd[ipft], deadcrootC_PFT_lasthd_ijd[ipft], grainC_PFT_lasthd_ijd[ipft], tlai_PFT_lasthd_ijd[ipft], tsai_PFT_lasthd_ijd[ipft])
                    
                    if (O3_damage_flag & O3_damage_scheme == 'Lombardozzi' & O3_sensitivity == 'custom') {
                       BGC_variable_readin_names = c(BGC_variable_readin_names, 'CUO_grainfill')
                       BGC_variable_readin_value = c(BGC_variable_readin_value, CUO_grainfill_PFT_lasthd_ijd[ipft])
                    } else {
                       BGC_variable_readin_names = c(BGC_variable_readin_names, 'CUO_grainfill')
                       BGC_variable_readin_value = c(BGC_variable_readin_value, NA)
                    }
                    
                    if (length(BGC_variable_readin_names) != length(BGC_variable_readin_value)) {stop("[simulate_ij.R] length(BGC_variable_readin_names) should be the same as length(BGC_variable_readin_value)")}
                    for (z in seq(BGC_variable_readin_names)) {assign(x = BGC_variable_readin_names[z], value = BGC_variable_readin_value[z])}
                    
                    ### Turn flags back to logical values
                    BGC_flagvariable_names = c("crop_live_flag","peak_LAI_flag","crop_plant_flag","leaf_emergence_flag","grain_fill_flag","harvest_flag")
                    for (z in seq(BGC_flagvariable_names)) {assign(x = BGC_flagvariable_names[z], value = as.logical(eval(parse(text = BGC_flagvariable_names[z]))))}
                    
                    rm(z, BGC_flagvariable_names, BGC_variable_readin_names, BGC_variable_readin_value)
                    
                }
            }
            
            
            #### soyFACE project
            if (biogeochem_flag) {
                if (force_prescribed_LAI) {
                    LAI = daily_LAI
                    LAI_prev_day = LAI_dayMinus1
                    if (h == 1) {
                        print(paste0('starting LAI = ', signif(LAI,3)))
                    }
                }
            }

            ####################################################################
            
            # Calculate aerodynamic conductance and/or temperature and humidity profiles using Monin Obukhov theory:
            # Now g_ah is always computed using either one of the two schemes below. (Tai, Feb 2019)
            if (Monin_Obukhov_flag | use_TEMIR_ga_flag) {
               Monin_Obukhov = f_Monin_Obukhov(Z_0m=Z_0m, Z_atm=Z_atm, 
                                               H_sen=H_sen, ET=ET, 
                                               u_star=u_star, T_2m=T_2m, 
                                               T_atm=T_atm, theta_2m=theta_2m, 
                                               theta_atm=theta_atm, q_2m=q_2m, 
                                               q_atm=NULL, rho_atm=rho_atm, 
                                               u_atm=u_atm, P_disp=P_disp, 
                                               P_surf=P_surf)
               # Obukhov length (m):
               L_Obuk = Monin_Obukhov$L_Obuk
               # Aerodynamic conductance for scalars (m s^-1):
               g_ah = min(c(Monin_Obukhov$g_ah, 1e4), na.rm=TRUE)
            } else {
               # Make use of f_Obuk() and f_aerodyn_cond() in "drydep_toolbox.R" to calculate aerodynamic conductance.
               Monin_Obukhov = NULL
               # Obukhov length (m):
               L_Obuk = f_Obuk(rho=rho_atm, T.s=T_2m, ustar=u_star, H=H_sen)
               # Aerodynamic conductance for scalars (m s^-1):
               g_ah = 1/f_aerodyn_cond(ustar=u_star, z0=Z_0m, obk=L_Obuk, T.k=T_2m, cz=(Z_atm + Z_0m))
            }
            g_ah = max(c(g_ah, 1e-4), na.rm=TRUE)
            
            ####################################################################
            
            # Canopy and soil environment:
            
            # Canopy air temperature (K):
            # Surface temperature (at displacement height) is used as a proxy for canopy air temperature.
            T_a = if (Monin_Obukhov_flag) Monin_Obukhov$T_s else T_2m
            # Vegetation temperature (K):
            # Canopy air temperature is used as a proxy for vegetation temperature.
            T_v = T_a
            
            # Specific humidity in canopy air (kg kg^-1):
            # Surface humidity (at displacement height) is used as a proxy for humidity in canopy air.
            q_a = if (Monin_Obukhov_flag) Monin_Obukhov$q_s else q_2m
            # Vapor pressure in canopy air (Pa):
            e_a = P_atm*q_a/0.622
            
            # Vapor pressure deficit (kPa) (Sun, Oct 2018): 
            vpd = max(0, (f_esat(T_a) - e_a)/1000)
            
            # Ambient (canopy) air CO2 partial pressure (Pa):
            c_a = CO2_conc*1e-6*P_atm
            # Wind speed incident on leaf (equivalent to "u_star") (m s^-1):
            u_leaf = u_star
            
            # Soil wetness for root zone (0-1):
            soil_wetness_root = GWETROOT[i,j,h]
            # Soil wetness for top soil (0-1):
            soil_wetness_top = GWETTOP[i,j,h]
            # Soil volumetric water content for top soil (0-1):
            theta_wtop = theta_sat*soil_wetness_top
            # Cloud fraction (0-1):
            cldtot = CLDTOT[i,j,h]
            # Soil temperature for different layers
            if (biogeochem_flag) {
            T_soil1 = TSOIL1[i,j,h]
            T_soil2 = TSOIL2[i,j,h]
            T_soil3 = TSOIL3[i,j,h]
            T_soil4 = TSOIL4[i,j,h]
            T_soil5 = TSOIL5[i,j,h]
            
            # Stop simulating the grid if Tsoil is missing, especially at the coastal area
            if(any(is.na(c(T_soil1, T_soil2, T_soil3, T_soil4, T_soil5)))) {break}
            T_soil_array_input = array(NA, dim = length(soil_layer_depth))
            
            for (z in 1:length(soil_layer_depth)){
                T_soil_array_input[z] = eval(parse(text = paste0("T_soil",z)))
            }
            }
            
            
            ####################################################################
            
            # Ozone concentration in air
            # if (!O3_fixed_flag) O3_conc = O3_hourly[i,j,(n_day_whole-1)*24+h]
            # soyFACE temporary 
            if (!O3_fixed_flag) O3_conc = O3_hourly[(d-1)*24+h]
            
            # Cumulative ozone update from previous time step (mmol m^-2):
            if (!O3_damage_flag) {
               CUO_prev_sun = 0
               CUO_prev_sha = 0
               CUO_prev_can = 0
            } else {
               if (h == 1 & (d == 1 & !continue_flag)) {
                  CUO_prev_sun = 0
                  CUO_prev_sha = 0
                  CUO_prev_can = 0
               } else if (h == 1 & (d > 1 | continue_flag)) {
                  # Cumulative ozone uptake of the last hour of last day:
                  # Date before current date:
                  previous_date = to.yyyymmdd(from.yyyymmdd(current_date) - 1*24)
                  filename = paste0(simulation_dir, 'temp_data/temp_', as.character(previous_date), '/temp_i', i_str, '_j', j_str, '.RData')
                  load(filename)
                  CUO_prev_sun = CUO_sun_PFT_lasthd_ijd[ipft]
                  CUO_prev_sha = CUO_sha_PFT_lasthd_ijd[ipft]
                  CUO_prev_can = CUO_can_PFT_lasthd_ijd[ipft]
               } else {
                  CUO_prev_sun = CUO_sun_PFT_hist_ijd[ipft,(h-1)]
                  CUO_prev_sha = CUO_sha_PFT_hist_ijd[ipft,(h-1)]
                  CUO_prev_can = CUO_can_PFT_hist_ijd[ipft,(h-1)]
               }
            }
            
            # LAI of previous time step (assume linear interpolation):
            if (!biogeochem_flag) {
            LAI_prev = LAI - (LAI - LAI_prev_day)/(24/dt_hr)
            } else {
                if (current_date != start_date) {
                    # In BGC mode, new LAI is calculated at h == 24, therefore LAI_prev and LAI are taken from temp_data of d = n - 2 and d = n - 1 respectively.
                    # The 'heal' in f_ozone_impact calculated from these two variables will be overestimated (underestimated) if d^2LAI/dt^2 < 0 (> 0)
                    LAI_prev = LAI - (LAI - LAI_prev_day)/(24/dt_hr)
                } else {
                    LAI_prev = 0
                }
            }
            
            ####################################################################
            
            # Find canopy transfer for PAR:
            # Default TEMIR method is used unless otherwise specified.
            if (cos_SZA <= 0) {
               # Canopy light extinction coefficient:
               K_b = 1e6
               # Absorption of direct beam and diffuse radiation by sunlit and shaded leaves (0-1):
               I_beam_sun = 0
               I_diff_sun = 0
               I_beam_sha = 0
               I_diff_sha = 0
               # Absorbed photosynthetically active radiation by sunlit leaves (W m^-2):
               phi_sun = 0
               # Absorbed photosynthetically active radiation by shaded leaves (W m^-2):
               phi_sha = 0
               # Sunlit leaf area index:
               LAI_sun = 0
               # Shaded leaf area index:
               LAI_sha = LAI
               # Surface albedo for visible light:
               surf_alb_beam = 0
               surf_alb_diff = 0
            } else {
               canopy_albedo = f_canopy_albedo(cos_SZA=cos_SZA, x_l=x_l, 
                                               LAI=LAI, SAI=SAI, 
                                               alpha_soil_dry=alpha_soil_dry, 
                                               alpha_soil_sat=alpha_soil_sat, 
                                               theta0=theta_wtop, 
                                               alpha_leaf=alpha_leaf, 
                                               alpha_stem=alpha_stem, 
                                               tau_leaf=tau_leaf, 
                                               tau_stem=tau_stem)
               # Canopy light extinction coefficient:
               K_b = canopy_albedo$K_b
               # Assume spherical leaf orientation if K_b cannot be found:
               if (is.na(K_b)) K_b = min(c(0.5/coz_SZA, 1e6), na.rm=TRUE)
               # Absorption of direct beam and diffuse radiation by sunlit and shaded leaves (0-1):
               # Handle the exceptions for prognostic crop in biogeochem. When LAI and SAI == 0 before planting, these outputs are NaN. (Pang, Jun 2019)
               # I_beam_sun = canopy_albedo$I_beam_sun
               # I_diff_sun = canopy_albedo$I_diff_sun
               # I_beam_sha = canopy_albedo$I_beam_sha
               # I_diff_sha = canopy_albedo$I_diff_sha
               I_beam_sun = ifelse(test = biogeochem_flag && LAI <= 1e-5, yes = 0, no = canopy_albedo$I_beam_sun)
               I_diff_sun = ifelse(test = biogeochem_flag && LAI <= 1e-5, yes = 0, no = canopy_albedo$I_diff_sun)
               I_beam_sha = ifelse(test = biogeochem_flag && LAI <= 1e-5, yes = 0, no = canopy_albedo$I_beam_sha)
               I_diff_sha = ifelse(test = biogeochem_flag && LAI <= 1e-5, yes = 0, no = canopy_albedo$I_diff_sha)
               
               # Find absorbed PAR:
               PAR_absorb = f_PAR_absorb(PAR_beam=PAR_beam, PAR_diff=PAR_diff, 
                                         LAI=LAI, SAI=SAI, K_b=K_b, 
                                         I_beam_sun=I_beam_sun, 
                                         I_diff_sun=I_diff_sun, 
                                         I_beam_sha=I_beam_sha, 
                                         I_diff_sha=I_diff_sha)
               # Absorbed photosynthetically active radiation by sunlit leaves (W m^-2):
               phi_sun = PAR_absorb$phi_sun
               if (is.na(phi_sun) || is.infinite(phi_sun)) phi_sun = 0
               # Absorbed photosynthetically active radiation by shaded leaves (W m^-2):
               phi_sha = PAR_absorb$phi_sha
               if (is.na(phi_sha) || is.infinite(phi_sha)) phi_sha = 0
               # Now calculate LAI_sun and LAI_sha explicitly from plant area index from f_PAR_absorb() (Tai, Feb 2019).
               # Sunlit leaf area index:
               # Handle the exception for prognostic crops. When the crops are not planted, LAI and SAI are zero. (Pang, Jun 2019)
               # LAI_sun = PAR_absorb$PAI_sun*LAI/(LAI + SAI)
               LAI_sun = if (LAI != 0 || SAI != 0) {PAR_absorb$PAI_sun*LAI/(LAI + SAI)} else {0}
               # Shaded leaf area index:
               # LAI_sha = PAR_absorb$PAI_sha*LAI/(LAI + SAI)
               LAI_sha = if (LAI != 0 || SAI != 0) {PAR_absorb$PAI_sha*LAI/(LAI + SAI)} else {0}
               # Surface albedo for visible light:
               # surf_alb_beam = canopy_albedo$I_beam_up
               # surf_alb_diff = canopy_albedo$I_diff_up
               # In theory surf_alb_beam and surf_alb_diff should be NA if LAI and SAI = 0, but it would cause error when saving the array as ncdf4, so we just put a zero here.
               surf_alb_beam = if (LAI != 0 || SAI != 0) {canopy_albedo$I_beam_up} else {NA}
               surf_alb_diff = if (LAI != 0 || SAI != 0) {canopy_albedo$I_diff_up} else {NA}
            }
            
            # Use simplified radiative transfer scheme that is consistent with Zhang et al. (2002) dry deposition mechanisms to override "phi_sun", "phi_sha", "LAI_sun", "LAI_sha" and "K_b" from default model above: (Tai, Feb 2019)
            if (simple_radiation_flag) {
               if (cos_SZA > 0) {
                  # Absorbed photosynthetically active radiation by sunlit vs. shaded leaves (W m^-2):
                  phi_sha = PAR_diff * exp(-0.5 * LAI^0.8) + 0.07 * PAR_beam * (1.1 - 0.1 * LAI) * exp(-cos_SZA)  # PAR_beam corrected to PAR_diff (Tai, Feb 2019)
                  phi_sun = phi_sha + PAR_beam^0.8 * 0.5 / cos_SZA
                  # Canopy light extinction coefficient:
                  K_b = 0.5 / cos_SZA
                  # Sunlit and shaded leaf area index:
                  # Note that this is exactly equivalent to the default method but with K_b = 0.5/cos_SZA for spherical leaf orientation and considering light extinction by leaves only (not by stems). (Tai, Feb 2019)
                  LAI_sun = (1 - exp(-K_b * LAI)) / K_b
                  LAI_sha = LAI - LAI_sun
                  # Use Zhang et al. 2001 to calculate radiative transfer (Wong)
                  # We now use Norman (1982) (c.f. Guether (1995)) to reproduce the result of GEOS-Chem Wesely scheme.
                  if (((LAI < 2.5 | swr < 200) & drydep_scheme == 'Zhang') | (drydep_scheme == 'Wesely')) {
                     # Please note, however, that Wesely scheme as exactly implemented in GEOS-Chem does not need phi_sun and phi_sha, but has its own way of accounting for light, although the original Wesely scheme does use them as input (Sun, Feb 2019).
                     # Absorbed photosynthetically active radiation by sunlit vs. shaded leaves (W m^-2):
                     phi_sha = PAR_diff * exp(-0.5 * LAI^0.7) + 0.07 * PAR_beam * (1.1 - 0.1 * LAI) * exp(-cos_SZA)
                     phi_sun = phi_sha + PAR_beam * 0.5 / cos_SZA
                  }
               } else {
                  # Canopy light extinction coefficient:
                  K_b = 1e6
                  # Absorbed photosynthetically active radiation by sunlit vs. shaded leaves (W m^-2):
                  phi_sun = 0
                  phi_sha = 0
                  # Fraction of sunlit and shaded leaf area
                  LAI_sun = 0
                  LAI_sha = LAI
               }
            }
            
            ####################################################################
            
            # Soil water stress function:
            beta_t = f_water_stress(soil_wetness=soil_wetness_root, 
                                    psi_sat=psi_sat, b_psi=b_psi, 
                                    psi_c=psi_c, psi_o=psi_o)
            
            # beta_t of irrigated crops are set to be one
            beta_t = ifelse(any(ipft == c(19,21,23,25)), yes = 1, no = beta_t)
            
            # Temporary: force al crops to be irrigated
            beta_t = ifelse(ipft >= 18, yes = 1, no = beta_t)
            
            # Find canopy photosynthesis:
            canopy_photosyn = f_canopy_photosyn(c_a=c_a, e_a=e_a, 
                                                phi_sun=phi_sun, 
                                                phi_sha=phi_sha, 
                                                T_v=T_v, P_atm=P_atm, 
                                                theta_atm=theta_atm, 
                                                C3_plant=C3_plant_flag, 
                                                V_cmax25_0=V_cmax25_0, 
                                                Phi_PSII=Phi_PSII, 
                                                Theta_PSII=Theta_PSII, 
                                                beta_t=beta_t, 
                                                acclimation=acclimation_flag, 
                                                T_10d=T_10d, lat=lat_rad, 
                                                decl=decl, 
                                                colimit=colimit_flag, 
                                                LAI=LAI, LAI_sun=LAI_sun, 
                                                K_n=K_n, K_b=K_b, 
                                                O3_damage=O3_damage_flag, 
                                                O3_conc=O3_conc, 
                                                scheme=O3_damage_scheme, 
                                                gs_scheme=gs_scheme, dt=dt, 
                                                LAI_prev=LAI_prev, 
                                                evergreen=evergreen_flag, 
                                                leaf_long=leaflong, 
                                                CUO_prev_sun=CUO_prev_sun, 
                                                CUO_prev_sha=CUO_prev_sha,
                                                CUO_prev_can=CUO_prev_can,
                                                sensitivity=O3_sensitivity, 
                                                plant_group=plant_group, 
                                                # u_atm=u_atm, 
                                                g_ah=g_ah, 
                                                biogeochem=biogeochem_flag, 
                                                leaf_N_conc=leaf_N_conc, 
                                                u_leaf=u_leaf, d_leaf=d_leaf, 
                                                met_cond=met_cond_flag, 
                                                tol=1e-3, g1_med=g1_med)
            
            # Total absorbed PAR (W m^-2):
            PAR_tot = phi_sun * LAI_sun + phi_sha * LAI_sha
            # Canopy-integrated stomatal resistance (s m^-1):
            R_s = min(c(1/(canopy_photosyn$g_s * (LAI_sun + LAI_sha)), 1e7), na.rm=TRUE)
            
            # Calculate dry deposition velocity (now only for ozone):
            if (!is.null(drydep_scheme)) {
               
               if (drydep_scheme == 'Wesely') {
                  
                  # Dry deposition parameter (Anthony)
                  # Base cuticular resistance (s m-1)
                  rcut = RLU_table[WESELY_MAPPING[ipft + 1]]
                  # Ground aerodynamic resistance (s m^-1)
                  ra_g = RAC_table[WESELY_MAPPING[ipft + 1]]
                  # Ground resistance for ozone (s m-1)
                  rg_o = RGO_table[WESELY_MAPPING[ipft + 1]]
                  # Ground resistance for SO2 (s m-1) 
                  rg_s = RGS_table[WESELY_MAPPING[ipft + 1]]
                  # Lower canopy resistance for ozone (s m-1)
                  rcl_o = RCLO_table[WESELY_MAPPING[ipft + 1]]
                  # Lower canopy resistance for SO2 (s m-1)
                  rcl_s = RCLS_table[WESELY_MAPPING[ipft + 1]]
                  # Minimin stomatal resistance resistance for ozone (s m-1)
                  rsmin = RI_table[WESELY_MAPPING[ipft + 1]]
                  # Baldocchi dry deposition polynomial coefficients
                  drycoeff = DRYCOEFF_table
                  
                  v_d = f_drydep_Wesely(rho = rho_atm, T.s = T_2m, H = H_sen, 
                                      z0 = Z_0m, cz = (Z_atm + Z_0m), 
                                      rsmin = rsmin, r_s = R_s, 
                                      use_temir_rs = use_TEMIR_gs_flag, 
                                      r_a = 1/g_ah, 
                                      use_temir_ra = use_TEMIR_ga_flag, 
                                      L = LAI, p = P_atm, mx = 48e-3, 
                                      ustar = u_star, srad = swr, 
                                      par = PAR_tot, hstar = 0.01, f0 = 1, 
                                      h = z_can, SAI = SAI, e_a = e_a, 
                                      cldfrac = cldtot, drycoeff = drycoeff, 
                                      r.cut0 = rcut, rcls = rcl_s, 
                                      rclo = rcl_o, rgs = rg_s, rgo = rg_o,
                                      ra_g = ra_g, cosSZA = cos_SZA, 
                                      co2_scale = CO2_scale_flag)
                  
               }
               
               if (drydep_scheme == 'Zhang') {
                  
                  # Determine whether the canopy is wet
                  # shsun 3/30/2018
                  # 1 mm hr^-1 = 2.78E-7 kg m^-2 s^-1
                  # rain.threshold = 0.2 / 3600 # Convert from mm hr^-1 to kg m^-2 s^-1
                  if (FLUXNET_flag == TRUE) rain.threshold = 0.2 else rain.threshold = 0.2/1000/3600
                  
                  is_wet = (prec_liq > rain.threshold) | (L_latent < 0)   # le: latent heat
                  
                  # Fraction covered by snow:
                  fsnow = min(d_snow / sdmax[mapping_TEMIR[ipft]], 1)
                  # Read in PFT-specific parameters:
                  if(is_wet) r.cut0 = r_cut_wO[mapping_TEMIR[ipft]] else r.cut0 = r_cut_dO[mapping_TEMIR[ipft]]
                  
                  rst.min = rst_min[mapping_TEMIR[ipft]]
                  b.rs = b_rs[mapping_TEMIR[ipft]]
                  t.min.st = t_min_st[mapping_TEMIR[ipft]] 
                  t.max.st = t_max_st[mapping_TEMIR[ipft]]
                  t.opt.st = t_opt_st[mapping_TEMIR[ipft]]
                  b.vpd = b_vpd[mapping_TEMIR[ipft]]
                  psi.min = psi_min[mapping_TEMIR[ipft]]
                  psi.max = psi_max[mapping_TEMIR[ipft]]
                  rgs = rg_s[mapping_TEMIR[ipft]]
                  rgo = rg_o[mapping_TEMIR[ipft]]
                  rac0_min = Rac0_min[mapping_TEMIR[ipft]]
                  rac0_max = Rac0_max[mapping_TEMIR[ipft]]
                  
                  # new PAR calculation (Anthony) 
                  v_d = f_drydep_Zhang(rho = rho_atm, T.s = T_2m, H = H_sen, 
                                     z0 = Z_0m, cz = (Z_atm + Z_0m), r_s = R_s, 
                                     use_temir_rs = use_TEMIR_gs_flag, 
                                     r_a = 1/g_ah, 
                                     use_temir_ra = use_TEMIR_ga_flag, 
                                     use_temir_beta = use_TEMIR_beta_flag, 
                                     rsmin = rst.min, brs = b.rs, 
                                     par.sun = phi_sun, par.sha = phi_sha, 
                                     Lsun = LAI_sun, Lsha = LAI_sha, 
                                     tmin = t.min.st, tmax = t.max.st, 
                                     topt = t.opt.st, vpd = vpd, bvpd = b.vpd, 
                                     e_a = e_a, srad = swr, psi.min = psi.min, 
                                     psi.max = psi.max, use.betat = FALSE, 
                                     betat = beta_t, L = LAI, LAI_min = LAI_min,
                                     LAI_max = LAI_max, rac0_min = rac0_min, 
                                     rac0_max = rac0_max, p = P_atm, 
                                     mx = 48e-3, ustar = u_star, hstar = 0.01, 
                                     f0 = 1, h = z_can, SAI = SAI, 
                                     fsnow = fsnow, is.wet = is_wet, 
                                     r.cut0 = r.cut0, rgs = rgs, rgo = rgo, 
                                     co2_scale = CO2_scale_flag)
                  
               }
               
            
            }
            
            #  Calculate biogeochemistry
            if (biogeochem_flag) {
                
                # Maintenance respiration (hourly)
                # CN ratios depend on phenology for crops,
                if (!crop || (crop && !grain_fill_flag)) {
                    leaf_cn_ratio = leaf_cn
                    fineroot_cn_ratio = froot_cn
                    coarseroot_cn_ratio = liveWood_cn
                    stem_cn_ratio = liveWood_cn
                    grain_cn_ratio = grain_cn
                } else if (crop && grain_fill_flag){
                    leaf_cn_ratio = leaf_cn_final
                    fineroot_cn_ratio = froot_cn_final
                    coarseroot_cn_ratio = liveWood_cn
                    stem_cn_ratio = stem_cn_final
                    grain_cn_ratio = grain_cn
                }
                
                maintenance_respiration_fluxes = f_maintenance_respiration_fluxes(woody.flag = woody_flag, root_frac_cnst_a = root_a_par, root_frac_cnst_b = root_b_par,
                                                                                  soil_depth_array = soil_layer_depth, T_soil_array = T_soil_array_input, T_2M = T_2m,
                                                                                  leaf_N = leaf_C / leaf_cn_ratio, livestem_N = livestem_C / stem_cn_ratio, livecoarseroot_N = livecoarseroot_C / coarseroot_cn_ratio, fineroot_N = fineroot_C / fineroot_cn_ratio, grain_N = grain_C / grain_cn_ratio,
                                                                                  root_fraction_array = root_frac_array)              
                
                single_timestep_mr = maintenance_respiration_fluxes$mr_total
                single_timestep_An = canopy_photosyn$A_can

                if (h == 1){
                    dailyMean_An = 0
                    dailyMean_mr = 0
                    
                }
                
                dailyMean_An = dailyMean_An + single_timestep_An
                dailyMean_mr = dailyMean_mr + single_timestep_mr
                
                # daily mean GPP and MR
                if (h == 24/dt_hr) {
                    dailyMean_An = dailyMean_An / (24/dt_hr)
                    dailyMean_mr = dailyMean_mr / (24/dt_hr)
                    # print(paste0('[simulated_ij.R] dailyMean_An = ', signif(dailyMean_An * 12.011e-6,4), ' gCm-2s-1  dailyMean_mr = ', signif(dailyMean_mr,4),' gCm-2s-1'))
                    
                }
                
                # Phenology, biomass partitioning, plant physiology (daily)
                if (h == 24/dt_hr){
                   
                   ### crop-ozone sensitivity is derived from top canopy leaves, therefore we use CUO_sun rather than CUO_can
                   ### use this value to track CUO at different growing stages
                   CUO_BGC = canopy_photosyn$CUO_sun

                    # Daily mean T2m and Tsoil for GDD calculation
                    
                    
                    # Plant phenology and biomass partitioning
                    if (evergreen_flag) {
                        evergreen_phenology = f_evergreen_phenology()
                        evergreen_allocation_fluxes = f_evergreen_allocation_fluxes()
                    } else if (stress_decid_flag) {
                        stress_deciduous_phenology = f_stress_deciduous_phenology()
                        stress_deciduous_allocation_fluxes = f_stress_deciduous_allocation_fluxes()
                    } else if (season_decid_flag) {
                        seasonal_deciduous_phenology = f_seasonal_deciduous_phenology()
                        seasonal_deciduous_allocation = f_seasonal_deciduous_allocation_fluxes()
                    } else if (crop_flag && !stress_decid_flag){
                        # prescribed_planting_date_readin are required if get_planting_date_option == 'prescribed-map' or 'prescribed-site'
                        # GDDx_20yr are required if get_GDDmat_method == 'CLM4.5'
                        # GDDmat_M/S/W are requried if get_GDDmat_method == 'Sack'
                        
                        crop_phenology = f_crop_phenology(T_10_d = T_10d, T_min_10_d = Tmin_10d, T_soil = T_soil1_daily, T2m = T_daily,
                                                          leafC = leaf_C, livestemC = livestem_C, finerootC = fineroot_C, grainC = grain_C, LAI = LAI, SAI = SAI,
                                                          CUO = CUO_BGC, CUO_grain_filling = CUO_grainfill,
                                                          GDD_T2m = GDDT2m, GDD_Tsoil = GDDTsoil, GDDmat = GDD_mat, GDDrepr = GDD_repr, GDDemer = GDD_emer,
                                                          crop_living_flag = crop_live_flag, crop_planting_flag = crop_plant_flag, leaf_emer_flag = leaf_emergence_flag, grain_filling_flag = grain_fill_flag, harvesting_flag = harvest_flag,
                                                          prescribed_planting_date_readin = if (get_planting_date_option != 'CLM4.5') {prescribed_planting_date} else {NULL}, planting_jday = day_of_planting, grain_filling_jday = day_of_grain_filling, harvest_jday = day_of_harvest,
                                                          GDD0_20yr = ifelse(exists("GDD0"), yes = GDD0, no = NULL), GDD8_20yr = ifelse(exists("GDD8"), yes = GDD8, no = 10), GDD10_20yr = ifelse(exists("GDD10"), yes = GDD10, no = NULL),
                                                          GDDmat_M = ifelse(exists("GDDmat_maize"), yes = GDDmat_maize, no = NULL), GDDmat_S = ifelse(exists("GDDmat_soybean"), yes = GDDmat_soybean, no = NULL), GDDmat_SW = ifelse(exists("GDDmat_springwheat"), yes = GDDmat_wheat, no = NULL), GDDmat_WW = ifelse(exists("GDDmat_winterwheat"), yes = GDDmat_winterwheat, no = NULL),
                                                          hybgdd = GDDmat_max, GDD_baseT = GDD_base_T, GDD_maxIncrease = GDDT2m_change_max, max_growing_season_length = crop_season_length_max, leaf_longevity = leaflong, T_plant_req = avg_T_planting_req, T_min_plant_req = min_T_planting_req, repr_GDDfrac = repr_GDDfraction, emer_GDDfrac = emer_GDDfraction,
                                                          prescribed_min_plant_jday = ifelse(test = at_NH_flag, yes = earliest_planting_doy_NH, no = earliest_planting_doy_SH), prescribed_max_plant_jday = ifelse(test = at_NH_flag, yes = latest_planting_doy_NH, no = latest_planting_doy_SH),
                                                          at_NH_flag = at_NH_flag, ipft = ipft)
                        
                        # Receive the outputs from the function using assign() instead of explictly writing them out if the function has mutliple outputs (>10)
                        # i.e. assign(x = simulateIJ_crop_phenology_variable_name, value = f_crop_phenology[[f_crop_phenology_output_listvariable_name]])
                        simulateIJ_crop_phenology_variable_name = c("crop_emergence_flux","leaf_C_loss_flux","grain_C_loss_flux","livestem_C_loss_flux","fineroot_C_loss_flux",
                                                                    "crop_live_flag","crop_plant_flag","leaf_emergence_flag","grain_fill_flag","harvest_flag",
                                                                    "day_of_planting", "day_of_grain_filling", "day_of_harvest",
                                                                    "GDDT2m","GDDTsoil", "GDD_mat", "GDD_repr", "GDD_emer",
                                                                    "leaf_C","grain_C","fineroot_C","livestem_C","LAI","SAI")
                        f_crop_phenology_output_listvariable_name = c("seedC_to_leafC_flux","leafC_loss_flux","grainC_loss_flux","livestemC_loss_flux", "finerootC_loss_flux",
                                                                      "croplive_flag", "cropplant_flag","leafemergence_flag","grainfill_flag", "har_flag",
                                                                      "planting_julianday", "grain_filling_julianday", "harvesting_julianday",
                                                                      "GDD_T2m","GDD_Tsoil", "GDDmat", "GDDrepr", "GDDemer",
                                                                      "leafC","grainC","finerootC","livestemC","LAI_out","SAI_out")
                        
                        # print(paste0('[simulated.ij] after crop_phenology leafC = ', signif(leaf_C, 5)))
                        
                        if (length(simulateIJ_crop_phenology_variable_name) != length(f_crop_phenology_output_listvariable_name)){stop("[simulate_ij.R f_crop_phenology] length(simulateIJ_crop_phenology_variable_name) != length(f_crop_phenology_output_listvariable_name)")}
                        for (z in seq(simulateIJ_crop_phenology_variable_name)){assign(x = simulateIJ_crop_phenology_variable_name[z], value = crop_phenology[[f_crop_phenology_output_listvariable_name[z]]])}
                        rm(z, simulateIJ_crop_phenology_variable_name, f_crop_phenology_output_listvariable_name)
                        
                        crop_allocation_fluxes = f_crop_allocation_fluxes(A_can_umolm2s1 = dailyMean_An, mr_total = dailyMean_mr, gr_fraction = if (ipft >= 18) {0.25} else {0.3}, 
                                                                          GDDmat = GDD_mat, GDD_T2m = GDDT2m, GDD_Tsoil = GDDTsoil, GDDemer = GDD_emer, GDDrepr = GDD_repr, crop_living_flag = crop_live_flag, peak_lai_flag = peak_LAI_flag, grain_filling_flag = grain_fill_flag,
                                                                          astem_leafem = a_stem_leafem, aleaf_leafem = a_leaf_leafem, aleaf = a_leaf, astem = a_stem,
                                                                          bfact = b_factor, arooti = a_root_initial, arootf = a_root_final, astemf = a_stem_final, declfact = GDD_decline_factor, allconss = a_stem_alloc_power, aleaff = a_leaf_final, allconsl = a_leaf_alloc_power, lfemerg = emer_GDDfraction, fleafi = a_leaf_base)
                        simulateIJ_crop_allocation_variable_name = c("leaf_C_alloc_flux","fineroot_C_alloc_flux","livestem_C_alloc_flux","deadstem_C_alloc_flux","livecoarseroot_C_alloc_flux","deadcoarseroot_C_alloc_flux","grain_C_alloc_flux",
                                                                     "a_stem","a_leaf","a_root","a_repr","a_stem_leafem","a_leaf_leafem")
                        f_crop_allocation_output_listvariable_name = c("leaf_carbon_partitioning_flux_gCm2s1","fineroot_carbon_partitioning_flux_gCm2s1","livestem_carbon_partitioning_flux_gCm2s1","deadstem_carbon_partitioning_flux_gCm2s1","livecoarseroot_carbon_partitioning_flux_gCm2s1","deadcoarseroot_carbon_partitioning_flux_gCm2s1","grain_carbon_partitioning_flux_gCm2s1",
                                                                       "astem","aleaf","aroot","arepr","astem_em","aleaf_em")
                        if (length(simulateIJ_crop_allocation_variable_name) != length(f_crop_allocation_output_listvariable_name)){stop("[simulate_ij.R f_crop_phenology] length(simulateIJ_crop_allocation_variable_name) != length(f_crop_allocation_output_listvariable_name)")}
                        for (z in seq(simulateIJ_crop_allocation_variable_name)) {assign(x = simulateIJ_crop_allocation_variable_name[z], value = crop_allocation_fluxes[[f_crop_allocation_output_listvariable_name[z]]])}
                        rm(z, simulateIJ_crop_allocation_variable_name, f_crop_allocation_output_listvariable_name)
                    }
                    
                    
                    
                    # Plant carbon pool addition and substraction
                    plant_Cpool_addition_subtraction = f_Cpool_addition_subtraction(grain_filling_flag = grain_fill_flag, harvesting_flag = harvest_flag, ipft = ipft,
                                                                                    leaf_C_prev = leaf_C,         leaf_C_biomass_partitioning_flux_gCm2s1 = leaf_C_alloc_flux,                     leaf_C_loss_flux_gCm2s1 = leaf_C_loss_flux,         seedC_to_leafC_flux_gCm2s1 = crop_emergence_flux,
                                                                                    livestem_C_prev = livestem_C, livestem_C_biomass_partitioning_flux_gCm2s1 = livestem_C_alloc_flux,             livestem_C_loss_flux_gCm2s1 = livestem_C_loss_flux,
                                                                                    deadstem_C_prev = 0,          deadstem_C_biomass_partitioning_flux_gCm2s1 = deadstem_C_alloc_flux,             deadstem_C_loss_flux_gCm2s1 = 0,
                                                                                    livecoarseroot_C_prev = 0,    livecoarseroot_C_biomass_partitioning_flux_gCm2s1 = livecoarseroot_C_alloc_flux, livecoarseroot_C_loss_flux_gCm2s1 = 0,
                                                                                    deadcoarseroot_C_prev = 0,    deadcoarseroot_C_biomass_partitioning_flux_gCm2s1 = deadcoarseroot_C_alloc_flux, deadcoarseroot_C_loss_flux_gCm2s1 = 0,
                                                                                    fineroot_C_prev = fineroot_C, fineroot_C_biomass_partitioning_flux_gCm2s1 = fineroot_C_alloc_flux,             fineroot_C_loss_flux_gCm2s1 = fineroot_C_loss_flux,
                                                                                    grain_C_prev = grain_C,       grain_C_biomass_partitioning_flux_gCm2s1 = grain_C_alloc_flux,                   grain_C_loss_flux_gCm2s1 = grain_C_loss_flux)
                    
                    leaf_C = plant_Cpool_addition_subtraction$leaf_C_new
                    # fineroot_C = plant_Cpool_addition_subtraction$fineroot_C_new
                    # livestem_C = plant_Cpool_addition_subtraction$livestem_C_new
                    deadstem_C = plant_Cpool_addition_subtraction$deadstem_C_new
                    # livecoarseroot_C = plant_Cpool_addition_subtraction$livecoarseroot_C_new
                    # deadcoarseroot_C = plant_Cpool_addition_subtraction$deadcoarseroot_C_new
                    # grain_C = plant_Cpool_addition_subtraction$grain_C_new
                    
                    # Plant physiology
                    vegetation_structure = f_vegetation_structure(slatop = sla_top, dsladlai = dsla_dlai, laimx = LAI_max, woody = woody_flag, crop = crop_flag, ipft = ipft, ztopmx = ztopmax,
                                                                  leafC = leaf_C, deadstemC = deadstem_C, tlai = LAI, tsai = SAI, peak_lai_flag = peak_LAI_flag, harvesting_flag = harvest_flag, crop_living_flag = crop_live_flag,
                                                                  CUO = CUO_BGC)
                    
                }
                
            } # End of biogeochemistry simulation
            
            # Outputs (i.e., history data):
            output_assign_df = `colnames<-`(rbind.data.frame(
               # Canopy-integrated variables:
               c('A_can_PFT_hist_ijd', 'canopy_photosyn$A_can'),
               c('R_can_PFT_hist_ijd', 'canopy_photosyn$R_can'),
               c('g_can_PFT_hist_ijd', 'canopy_photosyn$g_can'),
               c('CUO_can_PFT_hist_ijd', 'canopy_photosyn$CUO_can'),
               # Variables for sunlit and shaded leaves separately:
               c('LAI_sun_PFT_hist_ijd', 'LAI_sun'),
               c('LAI_sha_PFT_hist_ijd', 'LAI_sha'),
               c('A_nsun_PFT_hist_ijd', 'canopy_photosyn$A_nsun'),
               c('A_nsha_PFT_hist_ijd', 'canopy_photosyn$A_nsha'),
               c('R_dsun_PFT_hist_ijd', 'canopy_photosyn$R_dsun'),
               c('R_dsha_PFT_hist_ijd', 'canopy_photosyn$R_dsha'),
               c('g_ssun_PFT_hist_ijd', 'canopy_photosyn$g_ssun'),
               c('g_ssha_PFT_hist_ijd', 'canopy_photosyn$g_ssha'),
               c('g_s_PFT_hist_ijd', 'canopy_photosyn$g_s'),
               c('CUO_sun_PFT_hist_ijd', 'canopy_photosyn$CUO_sun'),
               c('CUO_sha_PFT_hist_ijd', 'canopy_photosyn$CUO_sha'),
               # Canopy radiative transfer variables:
               c('K_b_PFT_hist_ijd', 'K_b'),
               c('phi_sun_PFT_hist_ijd', 'phi_sun'),
               c('phi_sha_PFT_hist_ijd', 'phi_sha'),
               c('I_beam_sun_PFT_hist_ijd', 'I_beam_sun'),
               c('I_beam_sha_PFT_hist_ijd', 'I_beam_sha'),
               c('I_diff_sun_PFT_hist_ijd', 'I_diff_sun'),
               c('I_diff_sha_PFT_hist_ijd', 'I_diff_sha'),
               c('surf_alb_beam_PFT_hist_ijd', 'surf_alb_beam'),
               c('surf_alb_diff_PFT_hist_ijd', 'surf_alb_diff'),
               # Water stress variables:
               c('beta_t_PFT_hist_ijd', 'beta_t'),
               # Variables related to dry deposition:
               c('v_d_PFT_hist_ijd', 'v_d$v_d'),
               c('r_a_PFT_hist_ijd', 'v_d$r_a'),
               c('r_b_PFT_hist_ijd', 'v_d$r_b'),
               c('r_s_PFT_hist_ijd', 'v_d$r_s'),
               c('r_cut_PFT_hist_ijd', 'v_d$r_cut'),
               c('ra_lower_canopy_PFT_hist_ijd', 'v_d$ra_dc'),
               c('rc_ground_PFT_hist_ijd', 'v_d$rgx'),
               c('ra_ground_PFT_hist_ijd', 'v_d$ra_g'),
               c('rc_lower_canopy_PFT_hist_ijd', 'v_d$rclx'),
               c('is_wet_hist_ijd', 'is_wet'),
               c('g_par_hist_PFT_ijd', 'v_d$g_par'),
               c('g_t_hist_PFT_ijd', 'v_d$g_t'),
               c('g_w_hist_PFT_ijd', 'v_d$g_w'),
               c('g_vpd_hist_PFT_ijd', 'v_d$g_vpd'),
               c('w_st_hist_PFT_ijd', 'v_d$w_st'),
               c('f_snow_hist_PFT_ijd', 'v_d$fsnow'),
               # Canopy air micrometeorological variables:
               c('T_a_hist_ijd', 'T_a'),
               c('q_a_hist_ijd', 'q_a'),
               # Dataframe settings
               stringsAsFactors = FALSE), c('output_hist_name', 'assign_string'))
            
            # Additional outputs in biogeochemistry mode
            if (biogeochem_flag) {
                if (evergreen_flag) {
                    alloc_prefix = 'evergreen_allocation_fluxes'
                } else if (season_decid_flag) {
                    alloc_prefix = 'seasonal_deciduous_allocation_fluxes'
                } else if (stress_decid_flag) {
                    alloc_prefix = 'stress_deciduous_allocation_fluxes'
                } else if (crop_flag && !stress_decid_flag){
                    alloc_prefix = 'crop_allocation_fluxes'
                }
                
                output_assign_df = rbind.data.frame(output_assign_df,
                                                    # Maintenance respiration
                                                    c('mr_leaf_PFT_hist_ijd','maintenance_respiration_fluxes$mr_leaf'),
                                                    c('mr_fineroot_PFT_hist_ijd','maintenance_respiration_fluxes$mr_fineroot'),
                                                    c('mr_livestem_PFT_hist_ijd','maintenance_respiration_fluxes$mr_livestem'),
                                                    c('mr_livecoarseroot_PFT_hist_ijd','maintenance_respiration_fluxes$mr_coarseroot'),
                                                    c('mr_grain_PFT_hist_ijd','maintenance_respiration_fluxes$mr_grain'),
                                                    c('mr_total_PFT_hist_ijd','maintenance_respiration_fluxes$mr_total'),
                                                    # Physiology
                                                    c('LAI_PFT_hist_ijd','vegetation_structure$tlai'),
                                                    c('SAI_PFT_hist_ijd','vegetation_structure$tsai'),
                                                    c('htop_PFT_hist_ijd','vegetation_structure$canopy_top'),
                                                    c('hbot_PFT_hist_ijd','vegetation_structure$canopy_bottom'),
                                                    # C pools
                                                    c('leafC_PFT_hist_ijd', 'plant_Cpool_addition_subtraction$leaf_C_new'),
                                                    c('finerootC_PFT_hist_ijd', 'plant_Cpool_addition_subtraction$fineroot_C_new'),
                                                    c('livestemC_PFT_hist_ijd', 'plant_Cpool_addition_subtraction$livestem_C_new'),
                                                    c('deadstemC_PFT_hist_ijd', 'plant_Cpool_addition_subtraction$deadstem_C_new'),
                                                    c('livecoarserootC_PFT_hist_ijd', 'plant_Cpool_addition_subtraction$livecoarseroot_C_new'),
                                                    c('deadcoarserootC_PFT_hist_ijd', 'plant_Cpool_addition_subtraction$deadcoarseroot_C_new'),
                                                    c('grainC_PFT_hist_ijd', 'plant_Cpool_addition_subtraction$grain_C_new'),
                                                    # GPP, NPP, biomass partitioning fluxes and respirations
                                                    c('GPP_PFT_hist_ijd', paste0(alloc_prefix,'$daily_mean_GPP')),
                                                    c('NPP_PFT_hist_ijd', paste0(alloc_prefix,'$daily_mean_NPP')),
                                                    c('leafC_alloc_PFT_hist_ijd', paste0(alloc_prefix,'$leaf_carbon_partitioning_flux_gCm2s1')),
                                                    c('finerootC_alloc_PFT_hist_ijd', paste0(alloc_prefix,'$fineroot_carbon_partitioning_flux_gCm2s1')),
                                                    c('livestemC_alloc_PFT_hist_ijd', paste0(alloc_prefix,'$livestem_carbon_partitioning_flux_gCm2s1')),
                                                    c('deadstemC_alloc_PFT_hist_ijd', paste0(alloc_prefix,'$deadstem_carbon_partitioning_flux_gCm2s1')),
                                                    c('livecoarserootC_alloc_PFT_hist_ijd', paste0(alloc_prefix,'$livecoarseroot_carbon_partitioning_flux_gCm2s1')),
                                                    c('deadcoarserootC_alloc_PFT_hist_ijd', paste0(alloc_prefix,'$deadcoarseroot_carbon_partitioning_flux_gCm2s1')),
                                                    c('grainC_alloc_PFT_hist_ijd', paste0(alloc_prefix,'$grain_carbon_partitioning_flux_gCm2s1')),
                                                    # Crop phenology and flags
                                                    c('GDDT2m_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_phenology$GDD_T2m')),
                                                    c('GDDTsoil_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_phenology$GDD_Tsoil')),
                                                    c('GDDmat_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_phenology$GDDmat')),
                                                    c('GDDrepr_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_phenology$GDDrepr')),
                                                    c('GDDemer_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_phenology$GDDemer')),
                                                    c('crop_live_flag_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(FALSE), no = 'crop_phenology$croplive_flag')),
                                                    c('crop_plant_flag_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(FALSE), no = 'crop_phenology$cropplant_flag')),
                                                    c('leaf_emergence_flag_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(FALSE), no = 'crop_phenology$leafemergence_flag')),
                                                    c('grain_fill_flag_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(FALSE), no = 'crop_phenology$grainfill_flag')),
                                                    c('harvest_flag_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(FALSE), no = 'crop_phenology$har_flag')),
                                                    c('peak_LAI_flag_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(FALSE), no = 'vegetation_structure$peak_lai_flag')),
                                                    c('day_of_planting_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_phenology$planting_julianday')),
                                                    c('day_of_grain_filling_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_phenology$grain_filling_julianday')),
                                                    c('day_of_harvesting_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_phenology$harvesting_julianday')),
                                                    # Allocation coefficients
                                                    c('aleaf_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_allocation_fluxes$aleaf')),
                                                    c('aleaf_leafem_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_allocation_fluxes$aleaf_em')),
                                                    c('astem_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_allocation_fluxes$astem')),
                                                    c('astem_leafem_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_allocation_fluxes$astem_em')),
                                                    c('aroot_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_allocation_fluxes$aroot')),
                                                    c('arepr_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_allocation_fluxes$arepr'))
                )
                
                if (O3_damage_flag & O3_damage_scheme == 'Lombardozzi' & O3_sensitivity == 'custom') {
                   output_assign_df = rbind.data.frame(output_assign_df,
                                                       c('CUO_grainfill_PFT_hist_ijd', ifelse(test = ipft < 18, yes = as.character(NA), no = 'crop_phenology$CUO_grainfill')))
                } else {
                   output_assign_df = rbind.data.frame(output_assign_df,
                                                       c('CUO_grainfill_PFT_hist_ijd', 'NA'))
                }
            }
            
            # Assign outputs: 
            for (ivar in seq_along(var_name$variable_name)) {
               output_var =  paste0(var_name$variable_name[ivar], if (any(var_name$res_level[ivar] == c('PFT', 'PFT_daily'))) '_PFT_hist_ijd' else if (var_name$res_level[ivar] == 'grid') '_hist_ijd')
               
               # Assigning PFT_daily variables only at the last timestep, other types of variable are not affected
               if (var_name$res_level[ivar] != 'PFT_daily' || h == 24/dt_hr) {
                   current_assign = output_assign_df[match(output_var, output_assign_df$output_hist_name),]
                   temp_var = get(current_assign$output_hist_name)
                   
                   if (var_name$res_level[ivar] == 'PFT') {
                       temp_var[ipft,h] = eval(parse(text = current_assign$assign_string))
                   } else if (var_name$res_level[ivar] == 'grid') {
                       temp_var[h] = eval(parse(text = current_assign$assign_string))
                   } else if (var_name$res_level[ivar] == 'PFT_daily'){
                       temp_var[ipft] = eval(parse(text = current_assign$assign_string))
                   }    
                   assign(x = current_assign$output_hist_name, value = temp_var)
               }
            }
            
            rm(ivar, output_var, current_assign)
            # Assign missing essential outputs needed for next time step:
            # Cumulative ozone uptake:
            # 
            if (O3_damage_flag) {
               for (out_var in c('CUO_sun_PFT_hist_ijd', 'CUO_sha_PFT_hist_ijd')) {
                  current_assign = output_assign_df[match(out_var, output_assign_df$output_hist_name),]
                  temp_var = get(current_assign$output_hist_name)
                  temp_var[ipft,h] = eval(parse(text = current_assign$assign_string))
                  assign(x = current_assign$output_hist_name, value = temp_var)
               }
               rm(out_var, current_assign)
            }
            
            if (biogeochem_flag && h == 24/dt_hr) {
                BGC_temp_var_vec = c('LAI_PFT_hist_ijd', 'SAI_PFT_hist_ijd', 'leafC_PFT_hist_ijd', 'finerootC_PFT_hist_ijd', 'livestemC_PFT_hist_ijd', 'deadstemC_PFT_hist_ijd', 'livecoarserootC_PFT_hist_ijd', 'deadcoarserootC_PFT_hist_ijd', 'grainC_PFT_hist_ijd',
                                     'GDDT2m_PFT_hist_ijd', 'GDDTsoil_PFT_hist_ijd', 'GDDmat_PFT_hist_ijd', 'GDDemer_PFT_hist_ijd', 'GDDrepr_PFT_hist_ijd',
                                     'crop_live_flag_PFT_hist_ijd', 'crop_plant_flag_PFT_hist_ijd', 'leaf_emergence_flag_PFT_hist_ijd', 'grain_fill_flag_PFT_hist_ijd', 'harvest_flag_PFT_hist_ijd', 'peak_LAI_flag_PFT_hist_ijd',
                                     'day_of_planting_PFT_hist_ijd', 'day_of_grain_filling_PFT_hist_ijd', 'day_of_harvesting_PFT_hist_ijd', 'astem_PFT_hist_ijd', 'aleaf_PFT_hist_ijd', 'astem_leafem_PFT_hist_ijd', 'aleaf_leafem_PFT_hist_ijd')
                for (out_var in BGC_temp_var_vec) {
                    current_assign = output_assign_df[match(out_var, output_assign_df$output_hist_name),]
                    temp_var = get(current_assign$output_hist_name)
                    temp_var[ipft] = eval(parse(text = current_assign$assign_string))
                    assign(x = current_assign$output_hist_name, value = temp_var)
                }
                
                if (O3_damage_flag & O3_damage_scheme == 'Lombardozzi' & O3_sensitivity == 'custom') {
                   BGC_ozone_temp_var_vec = c('CUO_grainfill_PFT_hist_ijd')
                   
                   for (out_var in BGC_ozone_temp_var_vec) {
                      current_assign = output_assign_df[match(out_var, output_assign_df$output_hist_name),]
                      temp_var = get(current_assign$output_hist_name)
                      temp_var[ipft] = eval(parse(text = current_assign$assign_string))
                      assign(x = current_assign$output_hist_name, value = temp_var)
                   }
                   rm(BGC_ozone_temp_var_vec)
                }
                
                rm(out_var, current_assign, BGC_temp_var_vec)
            }
         }
         # End of simulation for all hours.
      }
   }

   # End of simulation for all PFTs for a single lon/lat.
   success = TRUE
   
   # Save temporary data needed for next time step:
   temp_env = new.env()
   temp_env$T_daily_ijd = T_daily
   temp_env$T_dailymin_ijd = Tmin_daily
   if (O3_damage_flag) {
      temp_env$CUO_sun_PFT_lasthd_ijd = CUO_sun_PFT_hist_ijd[,24/dt_hr]
      temp_env$CUO_sha_PFT_lasthd_ijd = CUO_sha_PFT_hist_ijd[,24/dt_hr]
      temp_env$CUO_can_PFT_lasthd_ijd = CUO_can_PFT_hist_ijd[,24/dt_hr]
   }
   
   if (biogeochem_flag) {
       # variables need to be transfered for the next time step in BGC mode
       # Physiology
       temp_env$tlai_PFT_lasthd_ijd = LAI_PFT_hist_ijd
       temp_env$tsai_PFT_lasthd_ijd = SAI_PFT_hist_ijd

       # Cpool
       temp_env$leafC_PFT_lasthd_ijd = leafC_PFT_hist_ijd
       temp_env$frootC_PFT_lasthd_ijd = finerootC_PFT_hist_ijd
       temp_env$livestemC_PFT_lasthd_ijd = livestemC_PFT_hist_ijd
       temp_env$deadstemC_PFT_lasthd_ijd = deadstemC_PFT_hist_ijd
       temp_env$livecrootC_PFT_lasthd_ijd = livecoarserootC_PFT_hist_ijd
       temp_env$deadcrootC_PFT_lasthd_ijd = deadcoarserootC_PFT_hist_ijd
       temp_env$grainC_PFT_lasthd_ijd = grainC_PFT_hist_ijd
       # GDD, phenology flags and day of planting/harvesting
       temp_env$GDDT2m_PFT_lasthd_ijd = GDDT2m_PFT_hist_ijd
       temp_env$GDDTsoil_PFT_lasthd_ijd = GDDTsoil_PFT_hist_ijd
       temp_env$GDDmat_PFT_lasthd_ijd = GDDmat_PFT_hist_ijd
       temp_env$GDDemer_PFT_lasthd_ijd = GDDemer_PFT_hist_ijd
       temp_env$GDDrepr_PFT_lasthd_ijd = GDDrepr_PFT_hist_ijd
       temp_env$crop_live_flag_PFT_lasthd_ijd = crop_live_flag_PFT_hist_ijd
       temp_env$crop_plant_flag_PFT_lasthd_ijd = crop_plant_flag_PFT_hist_ijd
       temp_env$leaf_emergence_flag_lasthd_ijd = leaf_emergence_flag_PFT_hist_ijd
       temp_env$grain_fill_flag_lasthd_ijd = grain_fill_flag_PFT_hist_ijd
       temp_env$harvest_flag_lasthd_ijd = harvest_flag_PFT_hist_ijd
       temp_env$peak_LAI_flag_lasthd_ijd = peak_LAI_flag_PFT_hist_ijd
       temp_env$planting_jday_lasthd_ijd = day_of_planting_PFT_hist_ijd
       temp_env$grain_filling_jday_lasthd_ijd = day_of_grain_filling_PFT_hist_ijd
       temp_env$harvesting_jday_lasthd_ijd = day_of_harvesting_PFT_hist_ijd
       # Allocation coefficients
       temp_env$astem_PFT_lasthd_ijd = astem_PFT_hist_ijd
       temp_env$aleaf_PFT_lasthd_ijd = aleaf_PFT_hist_ijd
       temp_env$astem_leafem_PFT_lasthd_ijd = astem_leafem_PFT_hist_ijd
       temp_env$aleaf_leafem_PFT_lasthd_ijd = aleaf_leafem_PFT_hist_ijd
           
       if (O3_damage_flag & O3_damage_scheme == 'Lombardozzi' & O3_sensitivity == 'custom') {
          temp_env$CUO_grainfill_PFT_lasthd_ijd = CUO_grainfill_PFT_hist_ijd
       }
   }
   filename = paste0(simulation_dir, 'temp_data/temp_', YYYY, MM, DD, '/temp_i', i_str, '_j', j_str, '.RData')
   save(list = ls(temp_env), envir = temp_env, file=filename)
   
   # Save output history data as a list for all c(i, j) for each d:
   # Note that there is an additional output "success" at the end that is not defined in "var_list".
   
   output = `names<-`(x = as.list(rep_len(x = NA, length = length(paste(var_name$variable_name)))), value = paste(var_name$variable_name))
   output_list_name = c()
   for (ivar in seq_along(var_name$variable_name)){
      current_var = var_name$variable_name[ivar]
      output[[var_name$variable_name[ivar]]] = get(paste0(current_var, ifelse(any(var_name$res_level[ivar] == c('PFT', 'PFT_daily')), '_PFT', ''), '_hist_ijd'))
      output_list_name = c(output_list_name, paste0(current_var, ifelse(any(var_name$res_level[ivar] ==  c('PFT', 'PFT_daily')), '_PFT', ''), '_hist'))
   }
   names(output) = output_list_name
   output[["success"]] = success
   rm(current_var, ivar, output_list_name)
   
   # # Save output history data externally for each c(i, j, d):
   # filename = paste0(hist_data_dir, hist_name, '/hist_', YYYY, MM, DD, '/hist_i', i_str, '_j', j_str, '.RData')
   # save(list=c(ls(pattern='_hist_ijd'), 'lon', 'lat', 'i', 'j', 'current_date'), file=filename)
   # output = success
   return(output)
   
}

###############################################################################

# Function to check and reshape a list of hourly history data into lon-lat gridded format for each simulation day:

f_hist_reshape = function(ij, hist_ij, err_check_only=FALSE) {
   
   # "ij" and "hist_ij" are list objects defined in "execution.R" and returned by function "f_simulate_ij", respectively.
   # It requires variables and parameters defined in the global environment.
   # It includes error checking.
   
   err_num = 0
   err_hist_ij = list()
   err_msg = NULL
   
   if (!err_check_only) hist_grid = array(NaN, dim=c(length(ind_lon), length(ind_lat), length(pftname), 24/dt_hr, nrow(var_name)))
   
   # Put in values of variables:
   for (n in 1:length(ij)) {
      
      i = ij[[n]][1]
      j = ij[[n]][2]
      ind_i = i - (ind_lon[1] - 1)
      ind_j = j - (ind_lat[1] - 1)
      if (nchar(i) > 3) stop('Max i >= 1000!')
      i_str = paste0(paste0(rep(0, 3-nchar(i)), collapse = ''), i)
      if (nchar(j) > 3) stop('Max j >= 1000!')
      j_str = paste0(paste0(rep(0, 3-nchar(j)), collapse = ''), j)
      
      if (is.atomic(hist_ij[[n]])) {
         
         # Record error message in "hist_ij":
         err_num = err_num + 1
         err_hist_ij[[err_num]] = hist_ij[[n]]
         
         # Error happened in simulation for this c(i, j, d). Mark down:
         new_msg = paste0('Error on ', YYYY, '/', MM, '/', DD, ': i = ', i_str, '; j = ', j_str, '; n = ', as.character(n))
         # print(new_msg, quote=FALSE)
         err_msg = c(err_msg, new_msg)
         
      } else {
         if (!err_check_only) {
            for (v in 1:nrow(var_name)) {
               if (var_name[v,4] == 'grid') {
                   hist_grid[ind_i,ind_j,1,,v] = hist_ij[[n]][[v]] 
               } else if (var_name[v,4] == 'PFT') {
                   hist_grid[ind_i,ind_j,,,v] = hist_ij[[n]][[v]]
               } else if (var_name[v,4] == 'PFT_daily'){
                   hist_grid[ind_i,ind_j,,24,v] = hist_ij[[n]][[v]]
               }
            }
         }
      }
      
   }
   # Done with all ij.
   # Return output:
   output = list(hist_grid=hist_grid, err_hist_ij=err_hist_ij, err_msg=err_msg)
   return(output)
   
}

###############################################################################

# Function to reshape hourly data array that consists of data for multiple days (applicable only for the debugging mode):

f_hourly_reshape = function(hist_grid, dim_day=3, dim_hr=5) {
   
   # This function reshapes the default hourly data output in the debugging mode by combining the two temporal dimensions (days and hours) into one single dimension of consecutive hours.
   # For this function to work, "dim_day" has to be 3 and "dim_hr" has to be either 4 or 5.
   
   dim_hist = dim(hist_grid)
   dim_hist_new = dim_hist[-dim_day]
   if (dim_day < dim_hr) dim_hr_new = dim_hr - 1 else stop('Dimensions are not correct!')
   dim_hist_new[dim_hr_new] = dim_hist[dim_day]*dim_hist[dim_hr]
   hist_new = array(NaN, dim=dim_hist_new)
   
   # Looping over simulation days:
   for (d in 1:dim_hist[dim_day]) {
      ind_hr = ((d - 1)*dim_hist[dim_hr] + 1):((d - 1)*dim_hist[dim_hr] + dim_hist[dim_hr])
      if (dim_day == 3 & dim_hr_new == 4) {
         hist_new[,,,ind_hr,] = hist_grid[,,d,,,]
      } else if (dim_day == 3 & dim_hr_new == 3) {
         hist_new[,,ind_hr,] = hist_grid[,,d,,]
      } else {
         stop('Dimensions are not correct!')
      }
   }
   
   return(hist_new)
   
}

###############################################################################
### End of module
###############################################################################
