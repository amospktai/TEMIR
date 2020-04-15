###############################################################################
### Script containing relevant constants
###############################################################################

# Physical constants:
sigma_SB = 5.670367e-8  # Stefan-Boltzmann constant (W m^-2 K^-4)
k_B = 1.38065e-23	      # Boltzmann constant (J K^-1)
N_A = 6.02214e23 	      # Avogadro constant (molecules mol^-1)
R_uni = k_B*N_A		   # Universal gas constant (J K^-1 mol^-1)
# R_uni = 8.314468	   # max. sig. fig. = 6
c_l = 2.99792458e8	   # Speed of light (m s^-1)
h_P = 6.62608e-34	      # Planck constant (J s)
h_hat = 1.05457e-34	   # Planck constant (J s)
G_uni = 6.673e-11	      # Universal gravitational constant (N m^2 kg^-2)
# pi = 3.14159265358979323846

# Geophysical constants (consistent with CESM usage):
g_E = 9.80616		      # Gravitational acceleration on Earth (m s^-2)
R_E = 6.37122e6		   # Mean Earth's radius (m)
Omega_E = 7.292115e-5   # Earth's angular velocity (rad s^-1)
d_SE = 1.496e11		   # Mean Sun-Earth distance (m)
M_da = 28.966e-3		   # Molar mass of dry air (kg mol^-1)
M_w = 18.016e-3		   # Molar mass of water (kg mol^-1)
R_da = R_uni/M_da		   # Gas constant for dry air (J K^-1 kg^-1)
# R_da = 287.0423	      # max. sig. fig. = 5
R_wv = R_uni/M_w		   # Gas constant for water vapor (J K^-1 kg^-1)
# R_wv = 461.5046	      # max. sig. fig. = 5
c_p = 1004.64		      # Specific heat of dry air at constant pressure (J K^-1 kg^-1)
c_lw = 4188			      # Specific heat of liquid water (J K^-1 kg^-1)
c_ice = 2117.27		   # Specific heat of ice (J K^-1 kg^-1)
lambda_vap = 2.501e6	   # Latent heat of vaporization (J kg^-1)
atm = 101325		      # Standard sea level pressure (Pa)
M_atm = 5.1352e18	      # Total mass of atmosphere (dry air only) (kg)
von_Kar = 0.4           # von Karman constant (unitless)

###############################################################################
### End of script
###############################################################################
