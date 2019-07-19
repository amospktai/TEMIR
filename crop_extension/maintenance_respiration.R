################################################################################
### Module for calculating maintenance respirations
################################################################################

################################################################################
### Functions:
################################################################################

# Maintenance respiration (gC m^-2 s^-1)
f_maintenance_respiration_fluxes = function(woody.flag, soil_depth_array, T_soil_array, T_2M, root_frac_cnst_a, root_frac_cnst_b,
                                            leaf_N, livestem_N, livecoarseroot_N, fineroot_N, grain_N,
                                            root_fraction_array) {

    # maintenance respiration rate per unit N at 20 Deg.C
    mr_perN_at20C = 2.525*10^-6  # unit: gC gN-1 s^-1
    # sensitivity (Q10) of respiration to temperature
    q10_mr = 1.5

    T_2m_degC = T_2M - 273.15
    T_soil_array_degC = T_soil_array - 273.15

    # print(paste0('[MR] T2m = ', T_2m_degC, " Tsoil = ", T_soil_array_degC[1]))

    mr_leaf = leaf_N * mr_perN_at20C * q10_mr ^ ((T_2m_degC - 20) / 10)
    mr_livestem = livestem_N * mr_perN_at20C * q10_mr ^ ((T_2m_degC - 20) / 10)
    mr_grain = grain_N * mr_perN_at20C * q10_mr ^ ((T_2m_degC - 20) / 10)
    mr_coarseroot = livecoarseroot_N * mr_perN_at20C * q10_mr ^ ((T_2m_degC - 20) / 10)

    mr_fineroot = 0
    for (layer in 1:length(root_fraction_array)){
        mr_fineroot = mr_fineroot + (fineroot_N * mr_perN_at20C * q10_mr ^ ((T_soil_array_degC[layer] - 20) / 10)) * root_fraction_array[layer]
    }

    # print(paste0("MR_leaf = ", signif(mr_leaf,3), " froot = ", signif(mr_fineroot,3), " stem = ", signif(mr_livestem,3), " grain = ", signif(mr_grain,3)))

    output = list(mr_leaf = mr_leaf, mr_livestem = mr_livestem, mr_grain = mr_grain, mr_coarseroot = mr_coarseroot, mr_fineroot = mr_fineroot,
                  mr_total = (mr_leaf + mr_livestem + mr_grain + mr_coarseroot + mr_fineroot))

    return(output)
}

################################################################################

# Root fraction of different PFTs for calculating root maintenance respiration
f_get_root_fraction = function(soil_depth_array, a_rootfrac, b_rootfrac){

    if (root_fraction_scheme == 'custom') {
       # user-defined root fraction
       # write your own scheme here
    } else if (root_fraction_scheme == 'CLM4.5') {
       # calculate root fraction based on the algorithm in CLM4.5
       # cumulative root fraction Y(d) = 1 - ((exp(a*d)) + exp(-b*d)) / 2
       cumulative_root_fraction_array = array(data = NA, dim = length(soil_depth_array))
       root_fraction_array = array(data = NA, dim = length(soil_depth_array))

       cumulative_root_fraction_array[1] = 0; root_fraction_array[1] = 0

       for (layer in 1:length(cumulative_root_fraction_array)) {
         cumulative_root_fraction_array[layer] = 1 - (exp(-a_rootfrac * soil_depth_array[layer]) + exp(-b_rootfrac * soil_depth_array[layer])) / 2

         if (layer >= 2){
           root_fraction_array[layer] = cumulative_root_fraction_array[layer] - cumulative_root_fraction_array[layer-1]
         }
       }
    }
    output = list(root_fraction_array = root_fraction_array, cumulative_root_fraction_array = cumulative_root_fraction_array)
}
