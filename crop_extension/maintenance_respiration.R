################################################################################
### Module for calculating maintenance respirations
################################################################################

################################################################################
### Functions:
################################################################################

f_maintenance_respiration_fluxes = function(woody.flag, soil_depth_array, T_soil_array, T_2M, root_frac_cnst_a, root_frac_cnst_b,
                                            leaf_N, livestem_N, livecoarseroot_N, fineroot_N, grain_N) {
    
    # maintenance respiration rate per unit N at 20 Deg.C
    mr_perN_at20C = 2.525*10^-6  # unit: gC gN-1 s^-1
    # sensitivity (Q10) of respiration to temperature
    q10_mr = 1.5
    
    T_2m_degC = T_2M - 273.15
    T_soil_array_degC = T_soil_array - 273.15
    
    mr_leaf = leaf_N * mr_perN_at20C * q10_mr ^ ((T_2m_degC - 20) / 10)
    mr_livestem = livestem_N * mr_perN_at20C * q10_mr ^ ((T_2m_degC - 20) / 10)
    mr_grain = grain_N * mr_perN_at20C * q10_mr ^ ((T_2m_degC - 20) / 10)
    mr_coraseroot = livecoarseroot_N * mr_perN_at20C * q10_mr ^ ((T_2m_degC - 20) / 10)
    
    # fine root maintenance respiration is the weighted sum of the root fraction from soil layers 
    
    root_fraction = f_get_root_fraction(soil_depth_array = soil_depth_array, a_rootfrac = root_frac_cnst_a, b_rootfrac = root_frac_cnst_b)
    root_fraction_array = root_fraction_array$root_fraction_array
    
    
    mr_fineroot = 0
    for (layer in 1:length(root_fraction_array)){
        mr_fineroot = mr_fineroot + (fineroot_N * mr_perN_at20C * q10_mr ^ ((Tfactor_Tsoil_array[layer] - 20) / 10)) * root_fraction_array[layer]
    }
    
    output = list(mr_leaf = mr_leaf, mr_livestem = mr_livestem, mr_grain = mr_grain, mr_coraseroot = mr_coraseroot, mr_fineroot = mr_fineroot,
                  mr_total = (mr_leaf + mr_livestem + mr_grain + mr_coraseroot + mr_fineroot))
    
    return(output)
}

f_get_root_fraction = function(soil_depth_array, a_rootfrac, b_rootfrac){
    
    if (root_fraction_scheme == 'custom') {
       # user-defined root fraction
      
    } else if (root_fraction_scheme == 'CLM4.5') {
       # calculate root fraction based on the algorithm in CLM4.5
       # cummulative root fraction Y(d) = 1 - ((exp(a*d)) + exp(-b*d)) / 2
       cummulative_root_fraction_array = array(data = NA, dim = length(soil_depth_array))
       root_fraction_array = array(data = NA, dim = length(soil_depth_array))
       
       cummulative_root_fraction_array[1] = 0; root_fraction_array[1] = 0
       
       for (layer in 1:length(cummulative_root_fraction_array)) {
         cummulative_root_fraction_array[layer] = 1 - (exp(a_rootfrac * soil_depth_array[layer]) + exp(-b_rootfrac * soil_depth_array[layer])) / 2
         
         if (layer >= 2){
           root_fraction_array[layer] = cummulative_root_fraction_array[layer] - cummulative_root_fraction_array[layer-1]
         }
       }
    }
    output = list(root_fraction_array = root_fraction_array, cummulative_root_fraction_array = cummulative_root_fraction_array)
}