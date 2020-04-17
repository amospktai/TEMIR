f_Cpool_addition_subtraction  = function(second_per_day = 86400, natural_mortality_rate = 0 / (365*86400), ipft, grain_filling_flag, harvesting_flag,
                                         leaf_C_prev, leaf_C_biomass_partitioning_flux_gCm2s1, leaf_C_loss_flux_gCm2s1, seedC_to_leafC_flux_gCm2s1 = NULL,
                                         livestem_C_prev, livestem_C_biomass_partitioning_flux_gCm2s1, livestem_C_loss_flux_gCm2s1,
                                         deadstem_C_prev, deadstem_C_biomass_partitioning_flux_gCm2s1, deadstem_C_loss_flux_gCm2s1,
                                         livecoarseroot_C_prev, livecoarseroot_C_biomass_partitioning_flux_gCm2s1, livecoarseroot_C_loss_flux_gCm2s1,
                                         deadcoarseroot_C_prev, deadcoarseroot_C_biomass_partitioning_flux_gCm2s1, deadcoarseroot_C_loss_flux_gCm2s1,
                                         fineroot_C_prev, fineroot_C_biomass_partitioning_flux_gCm2s1, fineroot_C_loss_flux_gCm2s1,
                                         grain_C_prev, grain_C_biomass_partitioning_flux_gCm2s1, grain_C_loss_flux_gCm2s1){
 
 leaf_C = leaf_C_prev
 livestem_C = livestem_C_prev
 deadstem_C = deadstem_C_prev
 livecoarseroot_C = livecoarseroot_C_prev
 deadcoarseroot_C = deadcoarseroot_C_prev
 fineroot_C = fineroot_C_prev
 grain_C = grain_C_prev 
  
 print(paste0('[Cpool addition before] (unit: gC m-2) leaf_C = ', signif(leaf_C, digits = 3), ' livestem_C = ', signif(livestem_C, digits = 3), ' fineroot_C = ', signif(fineroot_C, digits = 3) , ' grain_C = ', signif(grain_C, digits = 3)))
 
 
  # crop emergence flux
 if (ipft >= 18 && ipft <= 25) {
     if (!is.na(seedC_to_leafC_flux_gCm2s1)){leaf_C = leaf_C + seedC_to_leafC_flux_gCm2s1 * second_per_day}
 }
 
  #biomass partitioning of display pools
  leaf_C = leaf_C + leaf_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  livestem_C = livestem_C + livestem_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  deadstem_C = deadstem_C + deadstem_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  livecoarseroot_C = livecoarseroot_C + livecoarseroot_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  deadcoarseroot_C = deadcoarseroot_C + deadcoarseroot_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  fineroot_C = fineroot_C + fineroot_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  grain_C = grain_C + grain_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  
  #biomass partitioning of storage pools (for deciduous PFTs)
  # leaf_C_storage = leaf_C_storage + leaf_C_storage_biomass_partitioning_flux_gCm2s1 * second_per_day
  # ...
  # ...
  
  #biomass loss due to senescence or harvest (crops)
  leaf_C = leaf_C - leaf_C_loss_flux_gCm2s1 * second_per_day
  livestem_C = livestem_C - livestem_C_loss_flux_gCm2s1 * second_per_day
  deadstem_C = deadstem_C - deadstem_C_loss_flux_gCm2s1 * second_per_day
  livecoarseroot_C = livecoarseroot_C - livecoarseroot_C_loss_flux_gCm2s1 * second_per_day
  deadcoarseroot_C = deadcoarseroot_C - deadcoarseroot_C_loss_flux_gCm2s1 * second_per_day
  fineroot_C = fineroot_C - fineroot_C_loss_flux_gCm2s1 * second_per_day
  grain_C = grain_C - grain_C_loss_flux_gCm2s1 * second_per_day
  
  #biomass translocation (for crops)
  if (ipft >= 18 && ipft <= 25 && crop_translocation_scheme == 'JULES') {
    if (grain_filling_flag && !harvesting_flag){
        # In JULES, 100% of the leaf carbon loss is going to grain after DVI >= 1.5.
        # I have added a translocation efficiency term here as (probably) not all the C loss can be retranslocated
      
        translocation_flux = leaf_C_loss_flux_gCm2s1
        translocation_eff = 1
        grain_C = grain_C + (translocation_flux * translocation_eff * second_per_day)
    }
  }
  
  #natural mortality loss
  #only leaf, fineroot and livestem (display pools) have mortality loss
  leaf_C = leaf_C - natural_mortality_rate * leaf_C * second_per_day
  fineroot_C = fineroot_C - natural_mortality_rate * fineroot_C * second_per_day
  livestem_C = livestem_C - natural_mortality_rate * livestem_C * second_per_day

  #Ensure all the values are positive
  leaf_C = max(0, leaf_C)
  livestem_C = max(0, livestem_C)
  deadstem_C = max(0, deadstem_C)
  livecoarseroot_C = max(0, livecoarseroot_C)
  deadcoarseroot_C = max(0, deadcoarseroot_C)
  fineroot_C = max(0, fineroot_C)
  grain_C = max(0, grain_C)
  
  # print(paste0('[Cpool addition final] leaf_C = ', signif(leaf_C, digits = 3), ' livestem_C = ', signif(livestem_C, digits = 3), ' fineroot_C = ', signif(fineroot_C, digits = 3) , ' grain_C = ', signif(grain_C, digits = 3)))
  
  ##### soyFACE temporary 
  ### forst damage in year 2003, on doy 202. ~60% of leafC was lost
  if (current_date == 20030721) {
    leaf_C = leaf_C * 0.4
  }
  
  output = list(leaf_C_new = leaf_C, fineroot_C_new = fineroot_C, livestem_C_new = livestem_C, grain_C_new = grain_C, deadstem_C_new = deadstem_C, livecoarseroot_C_new = livecoarseroot_C, deadcoarseroot_C_new = deadcoarseroot_C)
  return(output)
}

