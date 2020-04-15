###############################################################################
### Module for canopy radiative transfer to calculate important inputs for the photosynthesis-stomatal conductance model
################################################################################

################################################################################
### Revision history
################################################################################

# Feb 2019: Now we clearly differentiate between PAI_sun/sha vs. LAI_sun/sha. See related updates in "Farquhar_Ball_Berry.R". (Tai)

################################################################################
### Functions:
################################################################################

# Solar declination angle:
f_decl = function(n_day, epsilon=23.44/180*pi, ecc=0.0167) {
   # This function calculates the solar declination angle (rad) as a function of Earth's obliquity ("epsilon" = 23.44 deg currently), eccentricity ("ecc" = 0.0167 currently), and the number of days ("n_day") since 00:00 UTC as Jan 1 begins.
   ang1 = 360/365.24*(n_day - 2)                                  # deg
   ang2 = 360/365.24*(n_day + 10) + 360/pi*ecc*sin(ang1/180*pi)   # deg
   decl = asin(sin(-epsilon)*cos(ang2/180*pi))                    # rad
   return(decl)
}

################################################################################

# Solar zenith angle:
f_cosSZA = function(lat, lon, n_day) {
   # This function uses external function "f_decl" to find the solar declination angle.
   # It calculates the cosine of the solar zenith angle as a function of latitude (rad, positive in NH), longitude (rad, positive east of Greenwich meridian), and solar declination angle (rad).
   # Hour angle (rad):
   hr_ang = 2*pi*n_day + lon
   # Solar declination angle:
   decl = f_decl(n_day=n_day)
   # Cosine of solar zenith angle:
   cos_SZA = sin(lat)*sin(decl) - cos(lat)*cos(decl)*cos(hr_ang)
   return(cos_SZA)
}

################################################################################

# Solar radiation absorption and albedo of canopy for vegetated land:

f_canopy_albedo = function(cos_SZA, x_l, LAI, SAI=0, alpha_soil_dry, alpha_soil_sat, theta0, alpha_leaf, alpha_stem=0, tau_leaf, tau_stem=0) {
   
   # This function calculates the canopy light extinction coefficient "K_b", and fractional absorption of direct beam and diffuse solar radiation by sunlit and shaded leaves, i.e., "I_beam_sun", "I_beam_sha", "I_diff_sun", "I_diff_sha".
   # Inputs: "cos_SZA" = cosine of solar zenith angle; "x_l" = leaf orientation; "LAI" = leaf area index; "SAI" = stem area index; "alpha_soil_dry" = albedo for dry soil; "alpha_soil_sat" = albedo for saturated soil; "theta0" = volumetric water content of top soil
   # Input parameters also include: reflectances ("alpha...") and transmittances ("tau...") for leaf ("...leaf") and stem ("...stem") of different plant functional types (PFTs)
   # "x_l" = +1 for horizontal leaves, 0 for random leaves, and -1 for vertical leaves. It can be treated as a constant for a given PFT.
   
   ### We assume a snow-free land surface. The capacity to include snow cover is under development. ###
   
   mu = cos_SZA
   chi_l = max(c(-0.4, min(c(x_l, 0.6), na.rm=TRUE)), na.rm=TRUE)
   phi1 = 0.5 - 0.633*chi_l - 0.33*chi_l^2
   phi2 = 0.877*(1 - 2*phi1)
   # Relative projected area of leaves and stems in the direction of SZA:
   G_mu = phi1 + phi2*mu
   # Canopy light extinction coefficient:
   K_b = G_mu/mu
   # Average inverse diffuse optical depth per unit leaf and stem area:
   mu_bar = (1/phi2)*(1 - phi1/phi2*log((phi1 + phi2)/phi1))
   
   # Weighted combination of leaf and stem reflectances:
   alpha = alpha_leaf*LAI/(LAI + SAI) + alpha_stem*SAI/(LAI + SAI)
   # Weighted combination of leaf and stem transmittances:
   tau = tau_leaf*LAI/(LAI + SAI) + tau_stem*SAI/(LAI + SAI)
   # Scattering coefficient:
   omega = alpha + tau
   # Cosine of mean leaf inclination angle relative to horizontal plane:
   cos_thetabar = 0.5*(1 + x_l)
   # Upscatter parameter for diffuse radiation:
   omega_beta = 0.5*(alpha + tau + (alpha - tau)*cos_thetabar^2)
   # Upscatter parameter for direct beam radiation:
   omega_beta0 = (1 + mu_bar*K_b)/(mu_bar*K_b)*(omega/2)*G_mu/(mu*phi2 + G_mu)*(1 - mu*phi1/(mu*phi2 + G_mu)*log((mu*phi1 + mu*phi2 + G_mu)/(mu*phi1)))
   
   # Ground albedo:
   Delta = max(c(0, (0.11 - 0.40*theta0)), na.rm=TRUE)
   alpha_g = min(c((alpha_soil_sat + Delta), alpha_soil_dry), na.rm=TRUE)
   
   # Two-stream approximation parameters:
   b = 1 - omega + omega_beta
   cc = omega_beta
   d = omega_beta0*mu_bar*K_b
   f = omega*mu_bar*K_b - omega_beta0*mu_bar*K_b
   h = sqrt(b^2 - cc^2)/mu_bar
   sigma = (mu_bar*K_b)^2 + cc^2 - b^2
   u1 = b - cc/alpha_g
   u2 = b - cc*alpha_g
   u3 = f + cc*alpha_g
   s1 = exp(-min(c(h*(LAI + SAI), 40), na.rm=TRUE))
   s2 = exp(-min(c(K_b*(LAI + SAI), 40), na.rm=TRUE))
   p1 = b + mu_bar*h
   p2 = b - mu_bar*h
   p3 = b + mu_bar*K_b
   p4 = b - mu_bar*K_b
   d1 = p1*(u1 - mu_bar*h)/s1 - p2*(u1 + mu_bar*h)*s1
   d2 = (u2 + mu_bar*h)/s1 - (u2 - mu_bar*h)*s1
   h1 = -d*p4 - cc*f
   h2 = (1/d1)*((d - p3*h1/sigma)*(u1 - mu_bar*h)/s1 - p2*(d - cc - (h1/sigma)*(u1 + mu_bar*K_b))*s2)
   h3 = (-1/d1)*((d - p3*h1/sigma)*(u1 + mu_bar*h)*s1 - p1*(d - cc - (h1/sigma)*(u1 + mu_bar*K_b))*s2)
   h4 = -f*p3 - cc*d
   h5 = (-1/d2)*((h4/sigma)*(u2 + mu_bar*h)/s1 + (u3 - (h4/sigma)*(u2 - mu_bar*K_b))*s2)
   h6 = (1/d2)*((h4/sigma)*(u2 - mu_bar*h)*s1 + (u3 - (h4/sigma)*(u2 - mu_bar*K_b))*s2)
   h7 = (cc/d1)*(u1 - mu_bar*h)/s1
   h8 = (-cc/d1)*(u1 + mu_bar*h)*s1
   h9 = (1/d2)*(u2 + mu_bar*h)/s1
   h10 = (-1/d2)*(u2 - mu_bar*h)*s1
   
   # Upward diffuse fluxes per unit incident direct beam and diffuse flux (i.e., surface albedo) (0-1):
   I_beam_up = h1/sigma + h2 + h3
   I_diff_up = h7 + h8
   
   # Downward diffuse fluxes per unit incident direct beam and diffuse radiation (0-1):
   I_beam_down = (h4/sigma)*exp(-K_b*(LAI + SAI)) + h5*s1 + h6/s1
   I_diff_down = h9*s1 + h10/s1
   
   # Direct beam and diffuse fluxes absorbed by vegetation per unit incident flux (0-1):
   I_beam = 1 - I_beam_up - (1 - alpha_g)*I_beam_down - (1 - alpha_g)*exp(-K_b*(LAI + SAI))
   I_diff = 1 - I_diff_up - (1 - alpha_g)*I_diff_down
   
   # Absorption of direct beam radiation by sunlit and shaded leaves (0-1):
   a1 = (h1/sigma)*(1 - s2^2)/(2*K_b) + h2*(1 - s2*s1)/(K_b + h) + h3*(1 - s2/s1)/(K_b - h)
   a2 = (h4/sigma)*(1 - s2^2)/(2*K_b) + h5*(1 - s2*s1)/(K_b + h) + h6*(1 - s2/s1)/(K_b - h)
   I_beam_sun = (1 - omega)*(1 - s2 + (a1 + a2)/mu_bar)
   I_beam_sha = I_beam - I_beam_sun
   
   # Absorption of diffuse radiation by sunlit and shaded leaves (0-1):
   a1 = h7*(1 - s2*s1)/(K_b + h) + h8*(1 - s2/s1)/(K_b - h)
   a2 = h9*(1 - s2*s1)/(K_b + h) + h10*(1 - s2/s1)/(K_b - h)
   I_diff_sun = ((1 - omega)/mu_bar)*(a1 + a2)
   I_diff_sha = I_diff - I_diff_sun
   
   output = list(K_b=K_b, alpha_g=alpha_g, I_beam_up=I_beam_up, I_beam_down=I_beam_down, I_diff_up=I_diff_up, I_diff_down=I_diff_down, I_beam_sun=I_beam_sun, I_beam_sha=I_beam_sha, I_diff_sun=I_diff_sun, I_diff_sha=I_diff_sha)
   return(output)
   
}

################################################################################

# Absorbed photosynthetically active (visible) radiation averaged over sunlit and shaded canopy (per unit plant area):

f_PAR_absorb = function(PAR_beam, PAR_diff, LAI, SAI, K_b, I_beam_sun, I_diff_sun, I_beam_sha, I_diff_sha) {
   
   # This function calculates the absorbed PAR averaged over sunlit and shaded canopy (per unit plant area): "phi_sun", "phi_sha" (W m^-2)
   # Inputs: "PAR_beam" = incident direct beam visible radiation (PAR) (W m^-2); "PAR_diff" = incident diffuse visible radiation (PAR) (W m^-2); "LAI" = leaf area index; "SAI" = stem area index; "K_b" = canopy light extinction coefficient; "I_beam_sun", "I_beam_sha", "I_diff_sun", "I_diff_sha" = fractional absorption of direct beam and diffuse PAR by sunlit and shaded leaves (0-1)
   
   # Sunlit and shaded plant area index:
   # Now we clearly differentiate between PAI_sun/sha vs. LAI_sun/sha (Tai, Feb 2019):
   # L_sun = (1 - exp(-K_b*(LAI + SAI)))/K_b
   # L_sha = (LAI + SAI) - L_sun
   PAI_sun = (1 - exp(-K_b*(LAI + SAI)))/K_b
   PAI_sha = (LAI + SAI) - PAI_sun
   
   # Absorbed PAR averaged over sunlit and shaded canopy (W m^-2):
   phi_sun = (I_beam_sun*PAR_beam + I_diff_sun*PAR_diff)/PAI_sun
   phi_sha = (I_beam_sha*PAR_beam + I_diff_sha*PAR_diff)/PAI_sha
   
   output = list(PAI_sun=PAI_sun, PAI_sha=PAI_sha, phi_sun=phi_sun, phi_sha=phi_sha)
   return(output)
   
}

################################################################################
