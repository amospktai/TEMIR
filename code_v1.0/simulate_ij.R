###############################################################################
### Module for the main simulation function to loop over lon and lat
### Including a function to reshape history data into gridded format
###############################################################################

# Both functions require access to variables in the global environment.

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
      }
   } ; rm(current_var, ivar)
   
   # Redefine new history outputs within function if not existent for next simulation step:
   if (O3_damage_flag) {
      for (out_var in c('CUO_sun_PFT_hist_ijd', 'CUO_sha_PFT_hist_ijd')) {
         assign(x = out_var, value = array(NA, dim=c(length(pftname), 24/dt_hr)))
      }
      rm(out_var)
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
      site_info = if (file.exists(paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx'))) read_excel(path = paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx')) else read_excel(path = paste0(FLUXNET_dir, 'FLUXNET_site_info.xlsx'))
      site_ind = match(FLUXNET_site_id, site_info$SITE_ID)
      Z_surf = as.double(unname(site_info[site_ind,"LOCATION_ELEV"]))
      rm(site_info, site_ind)
      # if elevation is not provided by FLUXNET then uses geopotential height of surface (m) as replacement
      if (is.na(Z_surf)) Z_surf = PHIS[i,j]/g_E
   } else {
      Z_surf = PHIS[i,j]/g_E
   }
   # Daily mean temperature (K):
   T_daily = mean(T2M[i,j,], na.rm=TRUE)
   
   # 10-day mean air temperature (K) to calculate temperature acclimation:
   acclimation_flag = ((d > 10) | continue_flag)
   if (acclimation_flag) {
      # Daily mean temperature of the past 10 days:
      T_daily_last10d = rep(NaN, times=10)
      for (d_prev in 10:1) {
         # Date of d_prev days before current date:
         previous_date = to.yyyymmdd(from.yyyymmdd(current_date) - d_prev*24)
         filename = paste0(simulation_dir, 'temp_data/temp_', as.character(previous_date), '/temp_i', i_str, '_j', j_str, '.RData')
         load(filename)
         T_daily_last10d[11-d_prev] = T_daily_ijd
      }
      # 10-day mean air temperature (K):
      T_10d = mean(T_daily_last10d, na.rm=TRUE)
   } else {
      # 10-day mean air temperature (K):
      T_10d = T_daily
   }
   # Soil parameters:
   # Saturated soil matric potential for bulk root zone (mm):
   psi_sat = psi_sat_bulk[i,j]
   # Saturated volumetric water content for bulk root zone (fraction):
   theta_sat = theta_sat_bulk[i,j]
   # Clapp and Homberger parameter for bulk root zone:
   b_psi = b_psi_bulk[i,j]
   # Saturated soil matric potential at top soil layer (mm):
   psi_sat_top = psi_sat_bulk_top[i,j]
   # Saturated volumetric water content at top soil layer (fraction):
   theta_sat_top = theta_sat_bulk_top[i,j]
   # Clapp and Homberger parameter at top soil layer (mm):
   b_psi_top = b_psi_bulk_top[i,j]
   # Saturated soil matric potential in bottom soil layer (mm):
   psi_sat_bottom = psi_sat_bulk_bottom[i,j]
   # Saturated volumetric water content in bottom soil layer (fraction):
   theta_sat_bottom = theta_sat_bulk_bottom[i,j]
   # Clapp and Homberger parameter in bottom soil layer:
   b_psi_bottom = b_psi_bulk_bottom[i,j]
   # Soil albedo for PAR for dry and saturated soil:
   # 3rd dim of "soil albedo" = [dry (visible), dry (Near IR), saturated (visible), saturated (Near IR)]
   alpha_soil_dry = soil_albedo[i,j,1]
   alpha_soil_sat = soil_albedo[i,j,3]
   
   # Other model parameters:
   met_cond_flag = TRUE
   colimit_flag = TRUE
   
   #############################################################################
   
   success = FALSE
   
   # PFT-specific parameters:
   for (ipft in (sim_PFT + 1)) { # + 1 to shift counter bare land into the index
      
      # Leaf area index (m^2 m^-2):
      n_PAI = if (leap & as.numeric(MM) > 2) n_day_whole - 1 else n_day_whole
      LAI = LAI_day_PFT[i,j,ipft,n_PAI]
      
      if (LAI < 0.01 | PFT_frac[i,j,ipft] < 0.01) {
         
         # Too little vegetation. Skip current PFT calculations.
         next
         
      } else {
         
         # Leaf area index range (m^2 m^-2) (shsun)
         LAI_min = min(LAI_day_PFT[i,j,ipft,], na.rm = TRUE)
         LAI_max = max(LAI_day_PFT[i,j,ipft,], na.rm = TRUE)
         
         # Stem area index (m^2 m^-2):
         SAI = SAI_day_PFT[i,j,ipft,n_PAI]
         # Consider LAI only?
         # SAI = 0
         # Leaf area index of previous day (m^2 m^-2):
         LAI_prev_day = if (n_PAI == 1) LAI_day_PFT[i,j,ipft,365] else LAI_day_PFT[i,j,ipft,n_PAI-1]
         
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
         # Medlyn model parameter:
         g1_med = g1_med_table[ipft]
         # Fraction of roots in top soil layer:
         root_frac_in_top = fraction_in_top[ipft]
         
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
            
            # Variables below are needed for dry deposition (Sun, Oct 2018):
            # Liquid Precipitation (kg m-2 s-1): 
            prec_liq = if (FLUXNET_flag) PRECTOT[i,j,h] else PRECTOT[i,j,h] - PRECSNO[i,j,h] 
            # Snow depth (m):
            d_snow = if (FLUXNET_flag) 0 else SNODP[i,j,h]
            # Latent heat (W m^-2):
            L_latent = EFLUX[i,j,h]
            # Atmospheric scale height (m):
            Z_scale = R_da*(T_2m + 0.0065*(Z_surf + Z_disp + Z_0m + 2)/2)/g_E
            # Surface pressure (Pa):
            P_surf = if (FLUXNET_flag) ATMP[i,j,h] else slp*exp(-Z_surf/Z_scale)
            # P_surf = slp*(1 - 0.0065*Z_surf/(T_2m + 0.0065*(Z_surf + Z_disp + Z_0m + 2)))^5.257
            # Pressure at zero-plane displacement height (d + z0m) where wind speed is extrapolated to zero (Pa):
            P_disp = if (FLUXNET_flag) ATMP[i,j,h] else slp*exp(-(Z_surf + Z_disp + Z_0m)/Z_scale)
            
            # At 2 m above displacement height:
            # Define T = theta at the surface.
            # Atmospheric pressure (Pa):
            P_2m = if (FLUXNET_flag) ATMP[i,j,h] else slp*exp(-(Z_surf + Z_disp + Z_0m + 2)/Z_scale)
            # Atmospheric potential temperature (K):
            theta_2m = if (FLUXNET_flag) T_2m + (g_E/c_p)*(Z_disp + Z_0m + 2) else T_2m*(P_surf/P_2m)^(R_da/c_p)
            # Vapor pressure (Pa):
            e_2m = P_2m*q_2m/(0.622 + 0.378*q_2m)
            # Moist air density (kg m^-3):
            rho_2m = (P_2m - 0.378*e_2m)/(R_da*T_2m)
            
            # At 10 m above displacement height:
            # Define T = theta at the surface.
            # Atmospheric pressure (Pa):
            P_10m = if (FLUXNET_flag) ATMP[i,j,h] else slp*exp(-(Z_surf + Z_disp + Z_0m + 10)/Z_scale)
            # Atmospheric potential temperature (K):
            # Define T = theta at the surface.
            theta_10m = if (FLUXNET_flag) T_10m + (g_E/c_p)*(Z_disp + Z_0m + 10) else T_10m*(P_surf/P_10m)^(R_da/c_p)
            # If "q_10m" is not provided, needed to scale it from the temperature difference (theta_10m - theta_2m).
            # Specific humidity (kg kg^-1):
            q_10m = if (!exists('q_10m')) { if (H_sen == 0 ) q_2m else { q_2m + (theta_10m - theta_2m)*c_p*ET/H_sen }} else q_10m
            # Vapor pressure (Pa):
            e_10m = P_10m*q_10m/(0.622 + 0.378*q_10m)
            # Moist air density (kg m^-3):
            rho_10m = (P_10m - 0.378*e_10m)/(R_da*T_10m)
            
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
            # g_aw = if (infer_canopy_met_flag) Monin_Obukhov$g_aw else NULL
            # Conversion from molar to meteorological conductance:
            mol_to_met = 1e-6*R_uni*theta_atm/P_atm
            
            # Delete temporary meteorological variables:
            rm(P_10m, T_10m, theta_10m, u_10m, rho_10m)
            
            ####################################################################
            # Calculate aerodynamic conductance and/or temperature and humidity profiles using Monin Obukhov theory:
            # Now g_ah is always computed using either one of the two schemes below. (Tai, Feb 2019)
            if (infer_canopy_met_flag | (ga_scheme == 'CLM4.5')) {
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
            T_a = if (infer_canopy_met_flag) Monin_Obukhov$T_s else T_2m
            # Vegetation temperature (K):
            # Canopy air temperature is used as a proxy for vegetation temperature.
            T_v = T_a
            
            # Specific humidity in canopy air (kg kg^-1):
            # Surface humidity (at displacement height) is used as a proxy for humidity in canopy air.
            q_a = if (infer_canopy_met_flag) Monin_Obukhov$q_s else q_2m
            rm(Monin_Obukhov)
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
            soil_wetness_top = if (FLUXNET_flag) GWETTOP[i,j,h] / theta_sat_top / 100 else GWETTOP[i,j,h]
            # Soil wetness for bottom soil layer (0-1):
            # Since GWETTOP and GWETROOT overlaps, in order to obtain wetness of bottom layer, the wetness of top layer is removed from GWETROOT here. 0.05 and 0.95 are the depths of top and bottom layers, respectively. Volumetric soil moisture is 
            # theta_w = theta_sat*soil_wetness
            # theta_w_top = theta_sat_top*soil_wetness_top 
            # theta_w_bottom = (theta_w - theta_w_top*0.05)/0.95
            soil_wetness_bottom = (theta_sat*soil_wetness_root - theta_sat_top*soil_wetness_top*0.05)/0.95 / theta_sat_bottom 
            # Soil volumetric water content for top soil (0-1):
            # ... (Ma, Sep 2019)
            # theta_wtop = theta_sat*soil_wetness_top
            theta_wtop = theta_sat_top*soil_wetness_top
            # Cloud fraction (0-1):
            cldtot = CLDTOT[i,j,h]
            
            ####################################################################
            
            # Ozone concentration in air
            if (!O3_fixed_flag) O3_conc = O3_hourly[i,j,(n_day_whole-1)*24+h]
            
            # Cumulative ozone update from previous time step (mmol m^-2):
            if (!O3_damage_flag) {
               CUO_prev_sun = 0
               CUO_prev_sha = 0
            } else {
               if (h == 1 & (d == 1 & !continue_flag)) {
                  CUO_prev_sun = 0
                  CUO_prev_sha = 0
               } else if (h == 1 & (d > 1 | continue_flag)) {
                  # Cumulative ozone uptake of the last hour of last day:
                  # Date before current date:
                  previous_date = to.yyyymmdd(from.yyyymmdd(current_date) - 1*24)
                  filename = paste0(simulation_dir, 'temp_data/temp_', as.character(previous_date), '/temp_i', i_str, '_j', j_str, '.RData')
                  load(filename)
                  CUO_prev_sun = CUO_sun_PFT_lasthd_ijd[ipft]
                  CUO_prev_sha = CUO_sha_PFT_lasthd_ijd[ipft]
               } else {
                  CUO_prev_sun = CUO_sun_PFT_hist_ijd[ipft,(h-1)]
                  CUO_prev_sha = CUO_sha_PFT_hist_ijd[ipft,(h-1)]
               }
            }
            
            # LAI of previous time step (assume linear interpolation):
            LAI_prev = LAI - (LAI - LAI_prev_day)/(24/dt_hr)
            
            ####################################################################
            
            # Find canopy transfer for PAR:
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
            } else if (radiative_scheme == 'two-stream') {
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
               I_beam_sun = canopy_albedo$I_beam_sun
               I_diff_sun = canopy_albedo$I_diff_sun
               I_beam_sha = canopy_albedo$I_beam_sha
               I_diff_sha = canopy_albedo$I_diff_sha
               # Find absorbed PAR:
               PAR_absorb = f_PAR_absorb(PAR_beam=PAR_beam, PAR_diff=PAR_diff, 
                                         LAI=LAI, SAI=SAI, K_b=K_b, 
                                         I_beam_sun=I_beam_sun, 
                                         I_diff_sun=I_diff_sun, 
                                         I_beam_sha=I_beam_sha, 
                                         I_diff_sha=I_diff_sha)
               # Absorbed photosynthetically active radiation by sunlit leaves (W m^-2):
               phi_sun = PAR_absorb$phi_sun
               if (is.na(phi_sun)) phi_sun = 0
               # Absorbed photosynthetically active radiation by shaded leaves (W m^-2):
               phi_sha = PAR_absorb$phi_sha
               if (is.na(phi_sha)) phi_sha = 0
               # Now calculate LAI_sun and LAI_sha explicitly from plant area index from f_PAR_absorb() (Tai, Feb 2019).
               # Sunlit leaf area index:
               LAI_sun = PAR_absorb$PAI_sun*LAI/(LAI + SAI)
               # Shaded leaf area index:
               LAI_sha = PAR_absorb$PAI_sha*LAI/(LAI + SAI)
               # Surface albedo for visible light:
               surf_alb_beam = canopy_albedo$I_beam_up
               surf_alb_diff = canopy_albedo$I_diff_up
               rm(canopy_albedo)
            } else if (radiative_scheme == 'Beer') {
               
               # Use simplified radiative transfer scheme that is consistent with Zhang et al. (2002) dry deposition mechanisms to override "phi_sun", "phi_sha", "LAI_sun", "LAI_sha" and "K_b" from default model above: (Tai, Feb 2019)
               
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
            }
            
            ####################################################################
            
            # Soil water stress function:
            beta_t = f_water_stress(soil_wetness=soil_wetness_root,
                                    soil_wetness_bottom=soil_wetness_bottom,
                                    soil_wetness_top=soil_wetness_top,
                                    root_frac_in_top=root_frac_in_top,
                                    theta_w=theta_w, theta_sat=theta_sat,
                                    psi_sat=psi_sat, b_psi=b_psi, 
                                    psi_c=psi_c, psi_o=psi_o,
                                    psi_sat_top=psi_sat_top, b_psi_top=b_psi_top, 
                                    psi_sat_bottom=psi_sat_bottom, b_psi_bottom=b_psi_bottom,
                                    multilayer = soil_layer_scheme)
            
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
            
            # Calculate dry deposition velocity (now ONLY for ozone):
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
                                        use_temir_rs = (gs_scheme_type == 'ecophysiological'), 
                                        r_a = 1/g_ah, 
                                        use_temir_ra = (ga_scheme == 'CLM4.5'), 
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
                  if (FLUXNET_flag) rain.threshold = 0.2 else rain.threshold = 0.2/1000/3600
                  
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
                                       use_temir_rs = (gs_scheme_type == 'ecophysiological'), 
                                       r_a = 1/g_ah, 
                                       use_temir_ra = (ga_scheme == 'CLM4.5'), 
                                       use_temir_beta = (gs_water_stress_scheme == 'CLM4.5'), 
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
            
            # Assign outputs: 
            for (ivar in seq_along(var_name$variable_name)) {
               output_var =  paste0(var_name$variable_name[ivar], if (var_name$res_level[ivar] == 'PFT') '_PFT_hist_ijd' else if (var_name$res_level[ivar] == 'grid') '_hist_ijd')
               current_assign = output_assign_df[match(output_var, output_assign_df$output_hist_name),]
               temp_var = get(current_assign$output_hist_name)
               if (var_name$res_level[ivar] == 'PFT') {
                  temp_var[ipft,h] = eval(parse(text = current_assign$assign_string))
               } else if (var_name$res_level[ivar] == 'grid') {
                  temp_var[h] = eval(parse(text = current_assign$assign_string))
               }
               assign(x = current_assign$output_hist_name, value = temp_var)
            }
            rm(ivar, output_var, current_assign)
            
            # Assign missing essential outputs needed for next time step:
            # Cumulative ozone uptake:
            if (O3_damage_flag) {
               for (out_var in c('CUO_sun_PFT_hist_ijd', 'CUO_sha_PFT_hist_ijd')) {
                  current_assign = output_assign_df[match(out_var, output_assign_df$output_hist_name),]
                  temp_var = get(current_assign$output_hist_name)
                  temp_var[ipft,h] = eval(parse(text = current_assign$assign_string))
                  assign(x = current_assign$output_hist_name, value = temp_var)
               }
               rm(out_var, current_assign)
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
   if (O3_damage_flag) {
      temp_env$CUO_sun_PFT_lasthd_ijd = CUO_sun_PFT_hist_ijd[,24/dt_hr]
      temp_env$CUO_sha_PFT_lasthd_ijd = CUO_sha_PFT_hist_ijd[,24/dt_hr]
   }
   filename = paste0(simulation_dir, 'temp_data/temp_', YYYY, MM, DD, '/temp_i', i_str, '_j', j_str, '.RData')
   save(list = ls(temp_env), envir = temp_env, file=filename)
   
   # Save output history data as a list for all c(i, j) for each d:
   # Note that there is an additional output "success" at the end that is not defined in "var_list".
   
   output = `names<-`(x = as.list(rep_len(x = NA, length = length(paste(var_name$variable_name)))), value = paste(var_name$variable_name))
   output_list_name = c()
   for (ivar in seq_along(var_name$variable_name)){
      current_var = var_name$variable_name[ivar]
      output[[var_name$variable_name[ivar]]] = get(paste0(current_var, ifelse(var_name$res_level[ivar] == 'PFT', '_PFT', ''), '_hist_ijd'))
      output_list_name = c(output_list_name, paste0(current_var, ifelse(var_name$res_level[ivar] == 'PFT', '_PFT', ''), '_hist'))
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
               if (var_name[v,4] != 'PFT') hist_grid[ind_i,ind_j,1,,v] = hist_ij[[n]][[v]] else hist_grid[ind_i,ind_j,,,v] = hist_ij[[n]][[v]]
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
