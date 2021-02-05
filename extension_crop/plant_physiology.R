f_vegetation_structure = function(slatop, dsladlai, laimx, woody, crop, ipft, ztopmx,
                                  leafC, deadstemC, tlai = LAI, tsai = SAI, peak_lai_flag, harvesting_flag, crop_living_flag,
                                  CUO){
    # constants used in CLM4.5
    taper = 200            # height/radius of wood
    stocking = 0.1         # no. individual/m2
    dwood = 2.5e5          # wood density (gC/m2)
    dtsmonth = 2592000     # no. of sec. in a 30-days month
    
    # using LAI and SAI at t-1 to calculate SAI later
    tlai_old = tlai
    tsai_old = tsai
    
    print(paste0('[physiology] LAI and SAI = ', signif(tlai_old, digits = 3),'  ', signif(tsai_old, digits = 3)))

    if (!force_prescribed_LAI) {
    # LAI calculation
    if (prognostic_LAI_calculation_scheme == 'CLM4.5'){
      if (dsladlai > 0) { #trees, shrubs
        tlai = slatop * exp(leafC * dsladlai - 1) / dsladlai
      } else { # crops, grasses
        tlai = slatop * leafC
        if (ipft >= 18 && ipft <= 25) {
          if (limit_crop_LAI_flag && tlai > laimx){
              peak_lai_flag == T
          } else {
              peak_lai_flag == F
          }
        }
      }
    } else if (prognostic_LAI_calculation_scheme == 'custom') {
        # soyFACE project:
        # modeling SLA as function of CO2 concentration and O3 uptake for soybean
        # SLA = SLA_375ppm + a * ([CO2] - 375)
        # SLA = SLA_noO3 * f(CUO)
        sla_375 = 0.055695   # from soyFACE data, unit: m^2 gC^-1 yr 2001 - 2008
        # sla_375 = 0.048783   # from soyFACE data, unit: m^2 gC^-1 yr 2009 - 2010
        
        sla_CO2_sensitivity = - 3.1611e-5 # derive from soyFACE data, for soybean only!!!
        # sla_f_O3 = 
        
        sla_co2 = sla_375 + sla_CO2_sensitivity * (CO2_conc - 375)
        
        if (O3_damage_flag & O3_damage_scheme == 'Lombardozzi' & O3_sensitivity == 'custom') {
          if (is.na(CUO)) {CUO = 0}
          sla_o3_factor = 1 + 0.0028 * CUO
          sla_o3_factor = max(0.01, min(1, sla_o3_factor))

        } else {
          sla_o3_factor = 1
        }
        
        
        tlai = sla_co2 * sla_o3_factor * leafC
        
        rm(list = c('sla_375', 'sla_CO2_sensitivity', 'sla_co2', 'sla_o3_factor'))
    }
    }
    
    tlai = max(0,tlai)
    
    # SAI calculation
    if(prognostic_SAI_calculation_scheme == 'CLM4.5'){
      if (ipft == 16 || ipft == 17) { # C3 unmanaged crops
        tsai_alpha = 1 - 1*dt/dtsmonth
        tsai_min = 0.1 * 0.5  # 0.5 is the scale to match MODIS derived value according to CLM
        tsai = max(tsai_alpha * tsai_old + max(tlai_old - tlai,0) ,tsai_min)
      } else if (!(ipft >= 18 && ipft <= 25)){ # other PFTs except prognostic crops
        tsai_alpha = 1 - 0.5*dt/dtsmonth
        tsai_min = 1 * 0.5
        tsai = max(tsai_alpha * tsai_old + max(tlai_old - tlai,0) ,tsai_min)
      } else { # prognostic crops
        if (harvesting_flag && tlai < 1e-4){ # after harvesting, SAI = 0.25
          tsai = 0.25
        } else if (ipft == 18 || ipft == 19){ # maize
          tsai = 0.1 * tlai
        } else {
          tsai = 0.2 * tlai
        }
      }
    } else if (prognostic_SAI_calculation_scheme == 'custom'){ # user-defined SAI calculation
      # tsai = ...
    }
    tsai = max(0,tsai)
    
    # canopy height (htop) and bottom (hbot)
    if (prognostic_canopy_height_calculation_scheme == 'CLM4.5') {
      if (!(ipft >= 18 && ipft <= 25)){
        if (woody) { # trees and shrubs
           if (ipft >= 10 && ipft <= 12) { # shrubs
             taper = 10
           } else {
             taper = 200
           }
           h_top = ((3 * deadstemC * taper^2) / (pi * stocking * dwood))^(1/3)
           h_top = max(h_top,0.01)
           
           h_bottom = max(0,min(3,htop-1))
        } else { # grasses and unmanaged crops
          h_top = max(0.25, tlai * 0.25)
          h_top = max(htop,0.01)
          h_bottom = max(0,min(0.05,h_top-0.2))
        }
      } else { #prognostic crops
        if (crop_living_flag){
          h_top = ztopmx * min(tlai/laimx-1,1)^2
          h_top = max(0.05, h_top)
          h_bottom = 0.02  
        } else {
          h_top = 0
          h_bottom = 0
        }
        
      }
    } else if (prognostic_canopy_height_calculation_scheme == 'custom'){ # user-defined h_top, h_bottom calculation
      # h_top = ....
      # h_bottom = ...
    }
    h_top = max(0, h_top)
    h_bottom = max(0, h_bottom)
    
    if (ipft >= 18 && ipft <= 25) {
      output = list(tlai = tlai, tsai = tsai, canopy_top = h_top, canopy_bottom = h_bottom, peak_lai_flag = peak_lai_flag)
    } else {
      output = list(tlai = tlai, tsai = tsai, canopy_top = h_top, canopy_bottom = h_bottom, peak_lai_flag = F)
    }
    
    # print(paste0('[physiology] LAI = ', signif(tlai, digits = 3)))
    return(output)
}

