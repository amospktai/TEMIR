################################################################################
### Module for computing surface flux parameters based on Monin-Obukhov Similarity Theory
################################################################################

# This module does not aim to calculate the surface fluxes, but to infer surface flux parameters and driving variables from prescribed or observed fluxes.

################################################################################
### Revision history
################################################################################

# Aug 2017: v1.0 finalized (Tai)
# Feb 2019: Revised the way Obukhov length is calculated to be consistent with GEOS-Chem method. (Tai)

################################################################################

# Flux-gradient parameter function for buoyant stability:

f_psi_zeta = function(zeta) {
   # This is the psi(zeta) function for momentum, heat and water vapor to account for deviation from static neutrality when calculating wind, temperature and moisture profiles.
   # zeta = (z - d)/L, where z is height, d is displacement height, and L is the Obukhov length.
   # The function computes "psi_m" for momentum, "psi_h" for heat and "psi_w" for water vapor.
   # A range of c(-100, 2) for "zeta" is imposed.
   if (zeta < 0) {
      if (zeta < -100) zeta = -100
      x = (1 - 16*zeta)^(1/4)
      psi_m = 2*log((1 + x)/2) + log((1 + x^2)/2) - 2*atan(x) + pi/2
      psi_h = 2*log((1 + x^2)/2)
      psi_w = psi_h
   } else {
      if (zeta > 2) zeta = 2
      psi_m = -5*zeta
      psi_h = psi_m
      psi_w = psi_h
   }
   output = list(psi_m=psi_m, psi_h=psi_h, psi_w=psi_w)
   return(output)
}

################################################################################

# Temperature (K) and specific humidity (kg kg^-1) at the zero-plane displacement height:

f_Monin_Obukhov = function(Z_0m, Z_atm=10, H_sen, ET, u_star, T_2m, T_atm, theta_2m, theta_atm, q_2m, q_atm=NULL, rho_atm, u_atm=NULL, P_disp, P_surf) {
   
   # It requires global constants: g_E, R_da, c_p, von_Kar
   # It requires external functions: f_psi_zeta, f_esat
   # "Z_0m" = roughness length for momentum (m)
   # "Z_atm" = reference height (m) above the zero-plane displacement height, which has to be > 2 m.
   # "H_sen" = sensible heat flux (W m^-2); "ET" = evapotranspiration flux (kg s^-1 m^-2)
   # "u_star" = friction velocity or characteristic velocity scale (m s^-1)
   # "T_2m" = temperature at 2 m above zero-plane displacement height (K); "T_atm" = temperature at Z_atm above zero-plane displacement height (K)
   # "q_2m" = specific humidity at 2 m above zero-plane displacement height (kg kg^-1); "q_atm" = specific humidity at Z_atm above zero-plane displacement height (kg kg^-1)
   # "u_atm" = wind speed at Z_atm above zero-plane displacement height (m s^-1)
   # "Monin_Obukhov=TRUE" instructs the function to calculate and return surface flux parameters and variables including surface temperature and specific humidity.
   
   # Surface flux parameters and variables:
   
   # Impose a lower bound on "Z_0m" to ensure conductances do not get too small:
   # Value 0.01 is the default Z_0m for soil, glacier and wetland in CLM4.5.
   if (Z_0m < 0.01) Z_0m = 0.01
   
   # We use values at a reference height of 10 m above displacement height to calculate characteristic scales and aerodynamic conductances.
   
   # Infer Obukhov length (m) from "u_atm" and "u_star":
   # diff_psi_m = log((10 + Z_0m)/Z_0m) - von_Kar*u_atm/u_star
   # f_val = function(zeta) f_psi_zeta(zeta)$psi_m - f_psi_zeta(zeta*Z_0m/(10 + Z_0m))$psi_m - diff_psi_m
   # if (sign(f_val(-100)) == sign(f_val(2))) {
   #    # No solution within likely range of c(-100, 2). Impose value:
   #    if (abs(f_val(-100)) < abs(f_val(2))) zeta = -100 else zeta = 2
   # } else {
   #    zeta = uniroot(f_val, interval=c(-100, 2))$root
   # }
   # L_Obuk = (10 + Z_0m)/zeta
   # Disabled this adhoc way to calculate L_Obuk; use a more conventional approach below. (Tai, Feb 2019)
   
   # Calculate Obukhov length (m) using conventional method consistent with GEOS-Chem and function obk() in "drydep_toolbox.R": (Tai, Feb 2019)
   denom = von_Kar*g_E*H_sen   # denominator
   L_Obuk = if (abs(denom) > 0) -rho_atm*c_p*T_2m*u_star^3/denom else 1e5
   
   # Roughness lengths for heat and water vapor (m):
   # Here we assume roughness lengths for momentum, heat and moisture are the same, following Zeng and Wang (2007).
   Z_0w = Z_0h = Z_0m
   
   # Aerodynamic conductance for momentum (m s^-1):
   if (is.null(u_atm)) {
      g_am = von_Kar*u_star/(log((Z_atm + Z_0m)/Z_0m) - f_psi_zeta((Z_atm + Z_0m)/L_Obuk)$psi_m + f_psi_zeta(Z_0m/L_Obuk)$psi_m)
   } else g_am = u_star^2/u_atm
   
   # Aerodynamic conductance for heat (m s^-1):
   g_ah = von_Kar*u_star/(log((Z_atm + Z_0m)/Z_0h) - f_psi_zeta((Z_atm + Z_0m)/L_Obuk)$psi_h + f_psi_zeta(Z_0h/L_Obuk)$psi_h)
   
   # Aerodynamic conductance for water vapor (m s^-1):
   g_aw = von_Kar*u_star/(log((Z_atm + Z_0m)/Z_0w) - f_psi_zeta((Z_atm + Z_0m)/L_Obuk)$psi_w + f_psi_zeta(Z_0w/L_Obuk)$psi_w)
   
   # Characteristic temperature scale (K):
   theta_star = H_sen/(-rho_atm*c_p*u_star)
   # Characteristic humidity scale (kg kg^-1):
   q_star = ET/(-rho_atm*u_star)
   
   # Calculate temperature (K) at zero-plane displacement height:
   theta_s = theta_atm - H_sen/(-rho_atm*c_p*g_ah)
   # theta_s = theta_atm - theta_star/von_Kar*(log((10 + Z_0m)/Z_0h) - f_psi_zeta((10 + Z_0m)/L_Obuk)$psi_h + f_psi_zeta(Z_0h/L_Obuk)$psi_h)
   # Check that theta_s should not be lower than theta_2m if H_sen > 0, and should not be higher than theta_2m if H_sen < 0:
   if (H_sen > 0 & theta_s < theta_2m) theta_s = theta_2m
   if (H_sen < 0 & theta_s > theta_2m) theta_s = theta_2m
   T_s = theta_s*(P_disp/P_surf)^(R_da/c_p)
   
   # Calculate specific humidity (kg kg^-1) at zero-plane displacement height:
   q_s = q_atm - ET/(-rho_atm*g_aw)
   # q_s = q_atm - q_star/von_Kar*(log((10 + Z_0m)/Z_0w) - f_psi_zeta((10 + Z_0m)/L_Obuk)$psi_w + f_psi_zeta(Z_0w/L_Obuk)$psi_w)
   # Prevent q_s from dropping below zero and being larger than q_sat(T_s):
   q_sat = f_esat(T_s)*0.622/P_disp
   q_s = max(c(0, min(c(q_s, q_sat), na.rm=TRUE)), na.rm=TRUE)
   # Check that q_s should not be lower than q_2m if ET > 0, and should not be higher than q_2m if ET < 0:
   if (ET > 0 & q_s < q_2m) q_s = q_2m
   if (ET < 0 & q_s > q_2m) q_s = q_2m
   
   # Water vapor pressure (Pa) at zero-plane displacement height:
   e_s = P_disp*q_s/(0.622 + 0.378*q_s)
   
   # Output:
   # possible outputs : c(T_s, q_s, e_s, L_Obuk, g_am, g_ah, g_aw, theta_s, theta_star, q_star, q_sat)
   return(list(T_s=T_s, q_s=q_s, e_s=e_s, L_Obuk=L_Obuk, g_am=g_am, g_ah=g_ah, g_aw=g_aw))
   
}

################################################################################
### End of module
################################################################################
