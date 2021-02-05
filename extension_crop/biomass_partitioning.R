################################################################################
### Module for calculating biomass partitioning flxues
################################################################################

################################################################################
### Functions:
################################################################################

f_evergreen_allocation_fluxes = function(){
    
}

f_seasonal_deciduous_allocation_fluxes = function(){
    
}

f_stress_deciduous_allocation_fluxes = function(){
    
}

f_crop_allocation_fluxes = function(A_can_umolm2s1, mr_total, gr_fraction = 0.3,
                                  GDDmat, GDD_T2m, GDD_Tsoil, GDDemer, GDDrepr, crop_living_flag, peak_lai_flag, grain_filling_flag, 
                                  astem_leafem, aleaf_leafem, aleaf, astem,
                                  bfact, arooti, arootf, astemf, declfact, allconss, aleaff, allconsl, lfemerg, fleafi){
    
    #convert A_can (GPP) from umol m^-2 s^-1 to gC m^-2 s^-1
    A_can_gCm2s1 = A_can_umolm2s1 * 12.011e-6
    
    # growth respiration = 0.3 * (GPP - MR)
    growth_respiration = max(0, (A_can_gCm2s1 - mr_total) * gr_fraction)
    NPP_gCm2s1 = A_can_gCm2s1 - growth_respiration - mr_total
    
    # Ensure allocation fluxes (or NPP) are positive
    NPP_gCm2s1 = max(0, NPP_gCm2s1)
    
    print(paste0('[biomass partitioning] (unit: gCm-2d-1) NPP = ', signif(NPP_gCm2s1 * 86400,4), ' MR_tol = ', signif(mr_total* 86400,4)))
    print(paste0('[biomass partitioning] NPP:GPP ratio = ', signif(NPP_gCm2s1/A_can_gCm2s1, digits = 3)))
    
      if (crop_biomass_partitioning_scheme == 'custom') {
        # derived from soyFACE data, suitable for soybean only
          if (crop_living_flag) {
              
            # Noted that the 'reproductive stage' in CLM4.5 and JULES are different.
            # In CLM4.5, it is called 'grain filling', the GDD requirment is 70% GDDmat (for soybean)
            # in JULES, it is called 'flowering', which is a bit earlier than 'grain filling'
            
              if (GDD_T2m < GDDemer) {
                DVI = -1 + GDD_T2m / GDDemer
              } else if (GDD_T2m >= GDDemer && GDD_T2m < GDDrepr) {
                DVI = 0 + (GDD_T2m - GDDemer) / (GDDrepr - GDDemer)
              } else if (GDD_T2m >= GDDrepr && GDD_T2m < GDDmat) {
                DVI = 1 + (GDD_T2m - GDDrepr) / (GDDmat - GDDrepr)
              }
            
            # caluclation of allocation coefficients is simliar to JULES, except aroot is not calibrated using multinodial logistic regression as soyFACE didn't record root mass
            # aroot in CLM4.5 amd JULES are roughly linear, therefore I suggest using a linear function instead of a logistic function
            # from JULES crop paper, y-int for aroot ~ 0.5, x-int for aroot ~ DVI = 1.4 (for soybean)
            
              f_CO2 = 1; f_O3 = 1
     
              aroot = f_CO2 * f_O3 * 0.5 - (0.5/1.4) * DVI
              aroot = min(1, max(0, aroot))
              
              above_ground_biomass_fraction = 1 - aroot
              
              # leaf_int = 10.773130; leaf_DVI = -11.324911
              # stem_int = 9.689051 ; stem_DVI = -8.874481

              leaf_int = 6.95774; leaf_DVI = -6.381277
              stem_int = 6.338262; stem_DVI = -5.172378

              denominator = 1 + exp(leaf_int + leaf_DVI * DVI) + exp(stem_int + stem_DVI * DVI)
              
              aleaf = exp(leaf_int + leaf_DVI * DVI) / denominator
              astem = exp(stem_int + stem_DVI * DVI) / denominator
              agrain = 1 / denominator
              
              aleaf = min(1, max(0, aleaf * above_ground_biomass_fraction))
              astem = min(1, max(0, astem * above_ground_biomass_fraction))
              agrain = min(1, max(0, agrain * above_ground_biomass_fraction))
              
          } else {
              aroot = 0; aleaf = 0; aroot = 0; agrain = 0
          }

      } else if (crop_biomass_partitioning_scheme == 'JULES'){

          JULES_crop_alloc = f_get_JULES_crop_allocation_coefficients(ipft = ipft, GDD_mat = GDDmat, GDD_T2m = GDD_T2m, GDD_emer = GDDemer, GDD_repr = GDDrepr, croplive_flag = crop_living_flag)
          astem = JULES_crop_alloc$p_stem
          aleaf = JULES_crop_alloc$p_leaf
          aroot = JULES_crop_alloc$p_root
          agrain = JULES_crop_alloc$p_harv
          
        } else if (crop_biomass_partitioning_scheme == 'CLM4.5'){

          CLM4.5_crop_alloc = f_get_CLM_crop_allocation_coefficients(a_stem_leafem = astem_leafem, a_leaf_leafem = aleaf_leafem, a_leaf = aleaf, a_stem = astem,
                                                                     GDDT2m = GDD_T2m, GDDTsoil = GDD_Tsoil, GDD_mat = GDDmat, GDD_repr = GDDrepr, GDD_emer = GDDemer,
                                                                     crop_live_flag = crop_living_flag, peak_LAI_flag = peak_lai_flag, grain_fill_flag = grain_filling_flag,
                                                                     bfact = bfact, arooti = arooti, arootf = arootf, astemf = astemf, declfact = declfact, allconss = allconss, aleaff = aleaff, allconsl = allconsl, lfemerg = lfemerg, fleafi = fleafi)
          astem = CLM4.5_crop_alloc$a_stem
          aleaf = CLM4.5_crop_alloc$a_leaf
          aroot = CLM4.5_crop_alloc$a_root
          agrain = CLM4.5_crop_alloc$a_grain
          astem_leafem = CLM4.5_crop_alloc$a_stem_leafemergence
          aleaf_leafem = CLM4.5_crop_alloc$a_leaf_leafemergence
        }
      
  
    
    # these allocation fluxes increase the C stock of display pools (i.e., C pool that affects physiology directly)
    # e.g., leaf C -> LAI, grain C -> yield
    print(paste0('[biomass partitioning] aleaf = ', signif(aleaf,3),' aroot = ', signif(aroot,3),' astem = ', signif(astem,3), ' agrain = ', signif(agrain,3)))
    leaf_carbon_partitioning_flux = NPP_gCm2s1 * aleaf
    fineroot_carbon_partitioning_flux = NPP_gCm2s1 * aroot
    livestem_carbon_partitioning_flux = NPP_gCm2s1 * astem
    grain_carbon_partitioning_flux = NPP_gCm2s1 * agrain
    
    # # these storage fluxes are only available for deciduous PFTs
    # # carbon is storaged in these pools until next growing season 
    # leaf_storage_carbon_partitioning_flux = NPP_gCm2s1 * a_leaf 
    # fineroot_storage_carbon_partitioning_flux = NPP_gCm2s1 * a_root
    # livestem_storage_carbon_partitioning_flux = NPP_gCm2s1 * a_stem
    # grain_storage_carbon_partitioning_flux = NPP_gCm2s1 * a_grain
    
    output = list(daily_mean_GPP = A_can_gCm2s1,
                  daily_mean_NPP = NPP_gCm2s1,
                  leaf_carbon_partitioning_flux_gCm2s1 = leaf_carbon_partitioning_flux, 
                  fineroot_carbon_partitioning_flux_gCm2s1 = fineroot_carbon_partitioning_flux, 
                  livestem_carbon_partitioning_flux_gCm2s1 = livestem_carbon_partitioning_flux, 
                  deadstem_carbon_partitioning_flux_gCm2s1 = 0, 
                  livecoarseroot_carbon_partitioning_flux_gCm2s1 = 0,
                  deadcoarseroot_carbon_partitioning_flux_gCm2s1 = 0,
                  grain_carbon_partitioning_flux_gCm2s1 = grain_carbon_partitioning_flux,
                  astem = astem,
                  aleaf = aleaf,
                  aroot = aroot,
                  arepr = agrain,
                  astem_em = ifelse(test = exists("astem_leafem"), yes = astem_leafem, no = NA),
                  aleaf_em = ifelse(test = exists("aleaf_leafem"), yes = aleaf_leafem, no = NA)
                  )
    
    return(output)
    
}



f_get_CLM_crop_allocation_coefficients = function(a_stem_leafem, a_leaf_leafem, a_leaf, a_stem,                                            # allocation coefficients
                                              crop_live_flag, peak_LAI_flag, grain_fill_flag,                                           # flags 
                                              GDDT2m, GDDTsoil, GDD_mat, GDD_repr, GDD_emer,                                         # GDD
                                              bfact, arooti, arootf, astemf, declfact, allconss, aleaff, allconsl, lfemerg, fleafi,    # CLM PFT constants
                                              current_date, ipft){                                                                     # model settings
    
    # print(paste0("[f_get_CLM_alloc] GDDT2m = ", signif(GDDT2m, 5), " GDDTsoil = ", signif(GDDTsoil, 5), "  GDD_mat = ", GDD_mat, " GDD_repr = ",GDD_repr, " GDD_emer = ", GDD_emer ))
    
    if (crop_live_flag) {
        
        if (GDDTsoil >= GDD_emer && GDDT2m < GDD_repr) { #vegetative stage
            if (limit_crop_LAI_flag && peak_LAI_flag) {
                a_leaf = 0; astem = 0; aroot = 1; agrain = 0
            } else {
                a_root = arooti - (arooti - arootf)*(GDDT2m/GDD_mat)
                a_leaf = (1-a_root)*(fleafi * (exp(-bfact) - exp(-bfact * GDDT2m/GDD_repr)) / (exp(-bfact)-1))
                a_stem = 1-a_root-a_leaf
                a_grain = 0
            }
          a_stem_leafem = a_stem; a_leaf_leafem = a_leaf
        } else if (GDDT2m >= GDD_repr && GDDT2m < GDD_mat) { #reproductive stage
          a_root = arooti - (arooti - arootf)*(GDDT2m/GDD_mat)
            if (a_leaf_leafem <= aleaff){
              a_leaf = a_leaf_leafem
            } else {
              a_leaf = max(aleaff, min(1, a_leaf * (1-(GDDT2m-GDD_repr)/(GDD_mat*declfact-GDD_repr))^allconsl))
            }
          
            if (a_stem_leafem <= astemf){
              a_stem = a_stem_leafem
            } else {
              a_stem = max(astemf,min(1, a_stem * (1-(GDDT2m-GDD_repr)/(GDD_mat*declfact-GDD_repr))^allconss))
            }
          a_grain = 1-a_root-a_leaf-a_stem
        } else {  # before emergence 
          a_leaf = 0; a_stem = 0; a_root = 0; a_grain = 0    
        }
    } else {      # after harvesting, before planting
        a_leaf = 0; a_stem = 0; a_root = 0; a_grain = 0
    }
    
  output = list(a_leaf = a_leaf, a_stem = a_stem, a_root = a_root, a_grain = a_grain, a_stem_leafemergence = a_stem_leafem, a_leaf_leafemergence = a_leaf_leafem)
  return(output)
}

f_get_JULES_crop_allocation_coefficients = function(ipft, GDD_mat, GDD_T2m, GDD_emer, GDD_repr, croplive_flag){
    
    #PFT specific parameters for allocation coefficients in JULES, Osborne et al. (2015)
    if (ipft == 24 || ipft == 25) { #soybean
        a_root_JULES = 20.0; b_root_JULES = -16.5; a_stem_JULES = 18.5; b_stem_JULES = -14.5; a_leaf_JULES = 19.5; b_leaf_JULES = -15.0
    } else if (ipft == 20 || ipft == 21) { #wheat
        a_root_JULES = 18.5; b_root_JULES = -20.0; a_stem_JULES = 16.0; b_stem_JULES = -15.0; a_leaf_JULES = 18.0; b_leaf_JULES = -18.5
    } else if (ipft == 18 || ipft == 19) { #maize
        a_root_JULES = 13.5; b_root_JULES = -15.5; a_stem_JULES = 12.5; b_stem_JULES = -12.5; a_leaf_JULES = 13.0; b_leaf_JULES = -14.0
    }
    
    #calculating Development Index (DVI) and partitioning coefficients
    if (croplive_flag){
        
        if (GDD_T2m < GDD_emer) {
            DVI_crop = -1 + GDD_T2m / GDD_emer
        } else if (GDD_T2m >= GDD_emer && GDD_T2m < GDD_repr) {
            DVI_crop = 0 + (GDD_T2m - GDD_emer) / (GDD_repr - GDD_emer)
        } else if (GDD_T2m >= GDD_repr && GDD_T2m <= GDD_mat){
            DVI_crop = 1 + (GDD_T2m - GDD_repr) / (GDD_mat - GDD_repr)
        } else {
            DVI_crop = NA
        }
        
        if(DVI_crop >= 0){
            p_root = exp(a_root_JULES + b_root_JULES * DVI_crop) / (exp(a_root_JULES + b_root_JULES * DVI_crop) + exp(a_stem_JULES + b_stem_JULES * DVI_crop) + exp(a_leaf_JULES + b_leaf_JULES * DVI_crop) + 1)
            p_leaf = exp(a_leaf_JULES + b_leaf_JULES * DVI_crop) / (exp(a_root_JULES + b_root_JULES * DVI_crop) + exp(a_stem_JULES + b_stem_JULES * DVI_crop) + exp(a_leaf_JULES + b_leaf_JULES * DVI_crop) + 1)
            p_stem = exp(a_stem_JULES + b_stem_JULES * DVI_crop) / (exp(a_root_JULES + b_root_JULES * DVI_crop) + exp(a_stem_JULES + b_stem_JULES * DVI_crop) + exp(a_leaf_JULES + b_leaf_JULES * DVI_crop) + 1)
            p_harv = 1 / (exp(a_root_JULES + b_root_JULES * DVI_crop) + exp(a_stem_JULES + b_stem_JULES * DVI_crop) + exp(a_leaf_JULES + b_leaf_JULES * DVI_crop) + 1)
        } else {
            p_root = 0; p_leaf = 0; p_stem = 0; p_harv = 0
        }
        
    } else {
        p_root = 0; p_leaf = 0; p_stem = 0; p_harv = 0; DVI_crop = NA
    }
    
    output = list(p_root = p_root, p_leaf = p_leaf, p_stem = p_stem, p_harv = p_harv, DVI_crop = DVI_crop)
    return(output)
}

f_custom_crop_allocation_coefficients = function(){}
