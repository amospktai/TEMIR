f_Cpool_addition_subtraction  = function(second_per_day = 86400, natural_mortality_rate = 2 / (365*86400), ipft, grain_filling_flag, harvesting_flag,
                                         leaf_C_prev, leaf_C_biomass_partitioning_flux_gCm2s1, leaf_C_loss_flux_gCm2s1, seedC_to_leafC_flux_gCm2s1 = NULL,
                                         livestem_C_prev, livestem_C_biomass_partitioning_flux_gCm2s1, livestem_C_loss_flux_gCm2s1,
                                         deadstem_C_prev, deadstem_C_biomass_partitioning_flux_gCm2s1, deadstem_C_loss_flux_gCm2s1,
                                         livecoraseroot_C_prev, livecoraseroot_C_biomass_partitioning_flux_gCm2s1, livecoraseroot_C_loss_flux_gCm2s1,
                                         deadcoraseroot_C_prev, deadcoraseroot_C_biomass_partitioning_flux_gCm2s1, deadcoraseroot_C_loss_flux_gCm2s1,
                                         fineroot_C_prev, fineroot_C_biomass_partitioning_flux_gCm2s1, fineroot_C_loss_flux_gCm2s1,
                                         grain_C_prev, grain_C_biomass_partitioning_flux_gCm2s1, grain_C_loss_flux_gCm2s1){
  
 leaf_C = leaf_C_prev
 livestem_C = livestem_C_prev
 deadstem_C = deadstem_C_prev
 livecoraseroot_C = livecoraseroot_C_prev
 deadcoraseroot_C = deadcoraseroot_C_prev
 fineroot_C = fineroot_C_prev
 grain_C = grain_C_prev 
  
  # crop emergence flux
 if (ipft >= 18 && ipft <= 25) {
   leaf_C = leaf_C + seedC_to_leafC_flux_gCm2s1 * second_per_day
 }
 
  #biomass partitioning of display pools
  leaf_C = leaf_C + leaf_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  livestem_C = livestem_C + livestem_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  deadstem_C = deadstem_C + deadstem_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  livecoraseroot_C = livecoraseroot_C + livecoraseroot_C_biomass_partitioning_flux_gCm2s1 * second_per_day
  deadcoraseroot_C = deadcoraseroot_C + deadcoraseroot_C_biomass_partitioning_flux_gCm2s1 * second_per_day
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
  livecoraseroot_C = livecoraseroot_C - livecoraseroot_C_loss_flux_gCm2s1 * second_per_day
  deadcoraseroot_C = deadcoraseroot_C - deadcoraseroot_C_loss_flux_gCm2s1 * second_per_day
  fineroot_C = fineroot_C - fineroot_C_loss_flux_gCm2s1 * second_per_day
  grain_C = grain_C - grain_C_loss_flux_gXm2s1 * second_per_day
  
  #biomass translocation (for crops)
  if (ipft >= 18 && ipft <= 25 && JULES_retranslocation_flag) {
    if (grain_fill_flag && !harvest_flag){
      translocation_flux = leaf_C_loss_flux_gCm2s
      grain_C = grain_C + translocation_flux * second_per_day
    }
  }
  
  #natural mortality loss
  #only leaf, fineroot and livestem (display pools) have mortality loss
  leaf_C = leaf_C - natural_mortality_rate * leaf_C * second_per_day
  fineroot_C = fineroot_C - natural_mortality_rate * fineroot_C * second_per_day
  livestem_C = livestem_C - natural_mortality_rate * livestem_C * second_per_day
  
  output = list(leaf_C_new = leaf_C, fineroot_C_new = fineroot_C, livestem_C_new = livestem_C, grain_C_new = grain_C)
  return(output)
}

