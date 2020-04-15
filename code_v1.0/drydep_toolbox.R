# Dry deposition simulator (with reference to Wesely scheme in GEOS-Chem)
# Coupled to TEMIR. R_a and g_s directly taken from TEMIR

# 9/7/2018 incorporate Anthony's 6_7_2018 version

################################################################################

# Monin-Obhukov Length (m)
# OBK = (- Air density * Cp * T(surface air) * Ustar^3)/(Kappa * g * Sensible Heat flux)
obk = function(rho, T.s, ustar, H) {
   # Uses global constants: "g_E = 9.80616" Gravitational acceleration on Earth (m s^-2); "von_Kar= 0.4" von Karman constant (unitless)"; "c_p = 1004.64" Specific heat of dry air at constant pressure (J K^-1 kg^-1);
   # Air density "rho" (molec m^-3), surface temperature "T.s" (K), friction velocity "ustar" (m s^-1), Sensible heat flux "H"
   den = von_Kar * g_E * H   # denominator
   OBK = if (abs(den) > 0) {- rho* c_p * T.s * ustar^3 / von_Kar / g_E / H} else 1E5
   return(OBK)
}

f_Obuk = obk

################################################################################

# Aerodynamic resistance (follows GEOS-Chem)
r.a = function(ustar, z0, obk, T.k, cz) {
   
   # Uses global constants:"von_Kar= 0.4" von Karman constant (unitless)"
   # Friction velocity "ustar" (m s^-1), roughness height "z0" (m), Monin-Obhukov Length "obk" (m), surface temperature "T.k" (K), reference height "cz" (m)
   # Calculate the kinematic viscosity xnu (m^2 s^-1) of air as a function of temperature
   xnu = 0.151 * (T.k/273.15) ^ 1.77 * 10^-4
   ckustar = von_Kar * ustar
   reyno = ustar * z0 / xnu # Reynolds number
   corr1 = cz/obk
   z0obk = z0/obk
   lrgera = FALSE # lrgera = T -> stable atmosphere; a high aerodynamic resistance (Ra = 1E4 m/s) is imposed, else Ra is calculated
   
   if (reyno < 0.1) {
      r_a = 1E4 # aerodynamically rough or smooth surface
   } else {
      if (corr1 < 0){
         # unstable conditions
         c.1 = (1 - 15 * corr1)^0.5
         c.2 = (1 - 15 * z0obk)^0.5
         c.3 = abs((c.1 - 1) / (c.1 + 1))
         c.4 = abs((c.2 - 1) / (c.2 + 1))
         r_a = 1 * (1 / ckustar) * log(c.3 / c.4)
      } else if ((corr1 >= 0) & (corr1 <= 1)) {
         r_a = (1/ckustar) * (log(corr1/z0obk) + 5 * (corr1 - z0obk))
      } else r_a = (1/ckustar) * (5 * log(corr1/z0obk) + (corr1 - z0obk))
   }
   
   return(min(r_a, 1E4))
   
}

f_aerodyn_cond = r.a

################################################################################

# Molecular diffusivity (m^2 s^-1)
diffg = function(T.k, p, mx){
   
   # Use global constants: "N_A = 6.02214e23" Avogadro constant (molecules mol^-1); "R_uni = 8.31446" universal gas constant (J K^-1 mol^-1); "M_da = 28.966e-3" Molar mass of dry air (kg mol^-1)
   # Surface temperature "T.k" (K), Pressure "p" (Pa), molecular weight "mx" (kg)
   
   # Air density "rho" (molec m^-3)
   rho = p * N_A / R_uni / T.k
   # "d" is the collision diameter for gas X in air
   d = 2.7 * 10^-10
   # Calculate the mean free path for gas X in air
   z = mx / M_da
   frpath = 1 / (pi * sqrt(1 + z) * rho * d^2)
   # Calculate average speed of gas X
   speed = sqrt(8 * R_uni * T.k / (pi * mx))
   # Calculate diffusion coefficient of gas X in air
   diff.g = (3 * pi / 32) * (1 + z) * frpath * speed
   return(diff.g)
   
}

################################################################################

# Boundary layer resistance
r.b = function(T.k, p, mx, ustar){
   
   # Uses external function "diff.g"
   # Uses global constants: "von_Kar= 0.4" von Karman constant (unitless)"
   # Temperature "T.k" (K), Pressure "p" (Pa), molecular weight "mx" (kg), friction velocity "ustar"  (m/s)
   
   # dair is the thermal diffusivity of air; value of 0.2e-4 m^2 s^-1
   dair = 0.2e-4
   diff.g = diffg(T.k, p, mx)
   r_b = 2 / (von_Kar * ustar) * (dair / diff.g)^0.667
   
   return(r_b)
}

################################################################################

# Henry's law factor for surface resistance
hstar.frac = function(hstar, f0){
   # Henry's law constant "hstar", reactivity factor for oxidation of biological substances "f0"
   return(hstar * 10^-5 +f0)
}

################################################################################

# Stomatal resistance  (follows GEOS-Chem v10.01, shsun)
r.s.wesely.cld = function(rsmin, L, cosSZA, CLDFRAC, drydep.coef, T.c, par){
   
   ## drydep_toolbox_mod.F90 & drydep_mod.F
   ## functions: BIOFIT, SUNPARAM_R4;SUNPARAM_R8
   ## BIOFIT_R4(drydep.coef, L, cosSZA, CLDFRAC, NPOLY) RESULT(BIO_FIT)
   
   # cosSZA Cosine, solar zenith angle
   # CLDFRAC cloud fraction (unitless)
   # drydep.coef (NPOLY) Baldocchi polynomial drydep coefficients
   # L leaf area index
   # GFACI, resultant light correction
   # par solar radiation in W m-2
   
   ### SUNPARAM
   nn = 3  ## number of variables (LAI, SUNCOS, CLDFRC)
   nd = c(55,20,11) ## scaling factor for each variable
   x0 = c(11,1,1) ## maximum of each variable
   X = c(L, cosSZA, CLDFRAC)
   for (i in 1:nn) {
      X[i] = min(X[i], x0[i])
      ## XLOW = minimum for each variable
      
      xlow = if (i != 3) x0[i]/nd[i] else 0
      
      X[i] = max(X[i], xlow)/x0[i]
   }
   
   ### BIOFIT_R4
   kk = 4
   term = c(1, L, cosSZA, CLDFRAC)
   term[2] = X[1] ## call SUNPARAM(term2)
   k = 0
   realterm = seq_along(drydep.coef)
   for (k3 in 1:kk) {
      for (k2 in k3:kk) {
         for (k1 in k2:kk) {
            k = k + 1
            realterm[k] = term[k1]*term[k2]*term[k3]
         }
      }
   }
   
   bio_fit = 0
   for (k in seq_along(drydep.coef)) {
      bio_fit = bio_fit + drydep.coef[k]*realterm[k]
   }
   print(paste0('Are bio_fit and drydep.coef*realterm identical? ', identical(bio_fit, drydep.coef * realterm)), quote=FALSE)
   if (bio_fit <= 0.1) bio_fit = 0.1
   
   ### drydep_mod.F
   if (rsmin < 9999) {
      
      gfact = if (T.c > 0 & T.c < 40) 400/T.c/(40-T.c) else 100
      gfaci = if (par > 0 & L > 0) 1/bio_fit else 100
      rsmin = rsmin*gfaci*gfact
   }
   
   return(rsmin)
}


################################################################################

# Stomatal resistance (Zhang et al. 2003)
r.s.zhang = function(rsmin, brs, par.sun, par.sha, Lsun, Lsha, T.c, tmin, tmax, 
                     topt, vpd, bvpd, srad, psi.min, psi.max, use.betat, betat) {
   
   # Minimum stomatal resistance "rsmin" (s m^-1), empirical light response constant for stomatal resistance "brs" (W m^-2)
   # Absorbed PAR averaged over sunlit and shaded canopy (per unit plant area): "phi_sun", "phi_sha" (W m^-2), sunlit and shaded plant area index: "Lsun", "Lsha" (m^2 m^-2)
   # surface temperature "T.c" (degC), minimum temperature "tmin"(degC), maximum temperature "tmax" (degC), optimum temperature "topt" (degC)
   # Vapour pressure deficit "vpd" (kPa^-1), water-vapour-pressure-deficit constant "bvpd" (kPa^-1)
   # "srad" = solar radiation (W m^-2), parameters specify leaf-water-potential dependency "psi.min" "psi.max" (MPa)
   # "use.betat = TRUE" specify to use TEMIR calculated water stress, water stress "betat", "is.wet = TRUE" canopy is wet
   
   # Unstressed canopy stomatal conductance G.st
   g.par = Lsun / (rsmin * (1 + brs / par.sun)) + Lsha / (rsmin * (1 + brs / par.sha))
   
   ## Temperature factor
   g.t = ((tmax - T.c) / (tmax - topt)) ^ ((tmax - topt) / (topt - tmin))
   g.t = g.t * (T.c - tmin) / (topt - tmin)
   g.t = max(g.t, 0.05)
   
   ## Water-vapour-pressure deficit factor
   g.d = max(0.05, 1 - vpd * bvpd)
   
   # Water stress (leaf water potential) factor
   if (use.betat) {
      g.w = betat
   } else {
      psi = -0.72 - 0.0013 * srad
      g.w = (psi - psi.min) / (psi.max - psi.min)
   }
   g.w = max(0.05, min(1, g.w))
   
   # Stomatal resistance
   r_s = 1 / (g.par * g.t * g.d * g.w)
   
   return(list(r_s = r_s, g_par = g.par, g_t = g.t, g_w = g.w, g_vpd = g.d))
   
}


################################################################################

# Cuticular resistance
r.cut.wesely = function(rcut0, T.c, L, rh, ustar){
   # Base cuticular resistance "r.cut0" (s m^-1), surface temperature "T.c" (degC), Leaf area index "L" (m^2 m^-2)
   RT = 1000 * exp(-T.c - 4)
   r_cut = if (T.c < -1) rcut0 * min(2, exp(0.2 * (-1 - T.c))) else rcut0 / L + RT
   return(r_cut)
}

################################################################################

# Canopy cuticular resistance (Zhang et al. 2003, eq.9a, eq.9b)
r.cut.zhang = function(rcut0, T.c, L, is.wet, rh, ustar){
   
   # Base cuticular resistance "r.cut0" (s m^-1), surface temperature "T.c" (degC), Leaf area index "L" (m^2 m^-2), "is.wet = TRUE" canopy is wet
   # Relative humidity "rh", friction velocity "ustar"  (m/s)
   r.cut = if (is.wet) rcut0 * ustar^-1 * L^-0.5 else r.cut = rcut0 * exp(-0.03 * rh) * L^-0.25 * ustar^-1
   r.cut = if (T.c < -1) r.cut * min(2, exp(0.2 * (-1 - T.c))) else r.cut
   return(r.cut)
}

################################################################################

# Lower canopy aerodynamic resistance (GC, wesely 1989)
r.dc = function(srad) {
   # A gas-phase transfer affected by buoyant convection in canopies
   # determined by the effects of mixing forced by buoyant convection when sunlight heats the ground or lower canopy and by penetration
   # of wind into canopies on the sides of hills.
   # RADIAT, solar radiation W/m2  RAD0 = RADIAT(IJLOOP)
   return(100 * (1 + 1000/(srad + 10)))
}

################################################################################

# In-canopy aerodynamic resistance (Zhang et al. eq.7 & eq.7a)
r.ac = function(rac0_min, rac0_max, LAI, LAI_min, LAI_max, ustar) {
   
   # leaf area index "LAI" (m^2 m^-2), friction velocity "ustar"  (m/s)
   rac0 = rac0_min + ((LAI - LAI_min)/(LAI_max - LAI_min))*(rac0_max - rac0_min)
   rac = (rac0*LAI^0.25)/(ustar^2)
   return(rac)
}

################################################################################

# Lower canopy surface resistance
rclx = function(hstar, f0, rcls, rclo){
   
   # Henry's law constant "hstar", reactivity factor for oxidation of biological substances "f0"
   # lower canopy resistance for SO2 "rcls" (s m^-1), lower canopy resistance for ozone "rclo" (s m^-1)
   return((hstar / (1E5 * rcls) + f0 / rclo) ^ -1)
}

################################################################################

# Ground surface resistance (GC)
rgx = function(hstar, f0, rgs, rgo, T.c){
   
   # Henry's law constant "hstar", reactivity factor for oxidation of biological substances "f0"; for O3, f0 = 1
   # ground resistance for SO2 "rgs" (s m^-1), ground resistance for ozone "rgo" (s m^-1)
   return((hstar / (1E5 * rgs) + f0 / rgo) ^ -1)
}

################################################################################

# Use wesely resistance model to compute dry deposition velocity
f_drydep_Wesely = function(rho = 0, T.s, H = 0, z0 = 0, cz = 0, rsmin = 0, r_s=NULL, use_temir_rs=FALSE, r_a=NULL, use_temir_ra=FALSE, drycoeff, cldfrac, cosSZA, L, p, mx, ustar, hstar, f0, h, SAI, srad, par, r.cut0, rcls, rclo, rgs, rgo, ra_g, e_a, co2_scale = FALSE) {
   
   # Uses external functions "r.a", r.b", "r.s.wesely", "obk", "r.cut", "hstar.frac", "ra.dc", "rclx", "rgx"
   # Aerodynamic resistance "r_a" (s m^-1), "use_temir_ra = TRUE" specifies to use TEMIR calculated ra, air density "rho" (molec m^-3)
   # Surface temperature "T.s" (K), sensible heat flux "H" (W m^-2), roughness height "z0" (m), altitude "cz" (m)
   # stomatal resistance "r_s" (s m^-1),"use_temir_rs" specified whether to use stomatal resistance calculated from TEMIR, minimum stomatal resistance "rsmin" (s m^-1)
   # Photosynthetic active radiation "phi_sun" "phi_sha"(W m^-2) (not used now)
   # Leaf area index "L" (m^2 m^-2), temperature "T.k" (K), Pressure "p" (Pa), molecular weight "mx" (kg), friction velocity "ustar"  (m/s)
   # Henry's law constant "hstar", reactivity factor for oxidation of biological substances "f0"
   # Canopy height "h" (m), stem area index "SAI" (m^2 m^-2)
   # lower canopy resistance for O3 "rclo" (s m^-1), ground surface resistance "rgs" (s m^-1)
   # Base cuticular resistance "r.cut0" (s m^-1), lower canopy resistance for SO2 "rcls" (s m^-1), lower canopy resistance for ozone "rclo" (s m^-1), ground resistance for SO2 "rgs" (s m^-1),
   # Ground resistance for ozone "rgo" (s m^-1), ground aerodynamic resistance "ra_g" (s m^-1)
   
   T.k = T.s
   T.c = T.k - 273.15
   
   rh = max(0, min(1, e_a/f_esat(T.k))) * 100
   
   # Raw stomatal resistance
   if (!use_temir_rs) {
      if ((L > 0) & (par > 0)){
         r_s = r.s.wesely.cld(rsmin = rsmin, cosSZA = cosSZA, T.c = T.c, L = L, CLDFRAC = cldfrac, par = par, drydep.coef = drycoeff)
         if (co2_scale & L > 0) r_s = r_s * (CO2_conc - 40)*440/(CO2_conc + 80)/320*CO2_conc/360 ### apply co2 effect
         r_s = min(r_s, 1E7)
      } else r_s = 1E7
   } else {
      if (is.null(r_s)) stop('r_s has not been provided as input!') else r_s = min(r_s, 1E7)
   }
   
   r_s = r_s * diffg(T.k, p, 18 * 10^-3) / diffg(T.k, p, mx)
   
   # Raw cuticular resistance
   if (L > 0) {
      r_cut = r.cut.wesely(r.cut0, T.c, L, rh, ustar) / hstar.frac(hstar, f0)
      r_cut = min(r_cut, 1E7)
   } else {
      r_cut = 1E7 # if LAI is too low, shut down stomatal and cuticular uptake
   }
   
   ra_dc = r.dc(srad) # in-canopy aerodynamic resistance
   rclx = rclx(hstar, f0, rcls, rclo) # lower canopy resistance (leaves, twig, bark, etc.)
   rgx = rgx(hstar, f0, rgs, rgo, T.c) # ground surface resistance (soil, leaf litter, etc.)
   
   # aerodynamic resistance
   if (!use_temir_ra) {
      obk = obk(rho, T.s, ustar, H)
      r_a = r.a(ustar, z0, obk, T.k, cz)
   } else {
      if (is.null(r_a)) stop('r_a has not been provided as input!')
      # Take r_a as provided.
   }
   
   # quasi-laminar resistance
   r_b = r.b(T.k, p, mx, ustar)
   
   ## Shsun (follows GC)
   g_c = 1/r_s + 1/r_cut + 1/(ra_dc + rclx) + 1/(rgx + ra_g)
   v_d = 1/(r_a + r_b + 1/g_c) #Calculate the dry deposition velocity
   out = list(r_a = r_a, r_b = r_b, r_s = r_s, r_cut = r_cut, ra_dc = ra_dc, rclx = rclx, ra_g = ra_g, rgx = rgx, v_d = v_d)
   return(out)
   
}

################################################################################

# Use Zhang et al. (2003) model to comput ozone dry deposition velocity
f_drydep_Zhang = function(rho = 0, T.s, H = 0, z0 = 0, cz =0, rsmin = 0, r_s=NULL, use_temir_rs=FALSE, r_a=NULL, use_temir_ra=FALSE, use_temir_beta, brs, par.sun, par.sha, Lsun, Lsha, tmin, tmax, topt, vpd, bvpd, e_a, srad, psi.min, psi.max, use.betat, betat, is.wet, L, LAI_min, LAI_max, rac0_min, rac0_max, p, mx, ustar, hstar, f0, h, SAI, r.cut0, rgs, rgo, fsnow, co2_scale = FALSE){
   
   # Use external function "r.s.zhang", "r.a", "r.b", "diffg", "r.cut.zhang", "ra.dc", "rgx"
   # Aerodynamic resisatance "r_a", "use_temir_ra = TRUE" specifies to use TEMIR calculated ra, air density "rho" (molec m^-3),
   # Surface temperature "T.s" (K), sensible heat flux "H" (W m^-2), roughness height "z0" (m), altitude "cz" (m)
   # Stomatal resistance "r_s"(s m^-1), "use_temir_rs = TRUE" specifies to use TEMIR calculated rs, minimum stomatal resistance "rsmin" (s m^-1)
   # Empirical light response conefficient "brs" (W m^-2),
   # Absorbed PAR averaged over sunlit and shaded canopy (per unit plant area): "phi_sun", "phi_sha" (W m^-2), sunlit and shaded plant area index: "Lsun", "Lsha" (m^2 m^-2)
   # Minimum temperature "tmin"(degC), maximum temperature "tmax" (degC), optimum temperature "topt" (degC)
   # Vapour pressure deficit "vpd" (kPa^-1), water vapour pressure deficit constant "bvpd" (kPa^-1), canopy air vapor pressure "e_a" (Pa)
   # solar radiation "srad" (W m^-2), parameters specify leaf-water-potential dependency "psi.min" "psi.max" (MPa)
   # "use.betat = TRUE" specify to use TEMIR calculated water stress, water stress "betat", "is.wet = TRUE" canopy is wet
   # Leaf area index "L" (m^2 m^-2), temperature "T.k" (K), Pressure "p" (Pa), molecular weight "mx" (kg), friction velocity "ustar"  (m/s)
   # Henry's law constant "hstar", reactivity factor for oxidation of biological substances "f0"
   # Canopy height "h" (m), stem area index "SAI" (m^2 m^-2)
   # Base cuticular resistance "r.cut0" (s m^-1), ground resistance for SO2 "rgs" (s m^-1), ground resistance for ozone "rgo" (s m^-1), snow cover fraction (0.0-1.0) "fsnow"
   
   T.k = T.s
   T.c = T.k - 273.15 # Kelvin to Celsius temperature
   
   rh = max(0, min(1, e_a/f_esat(T.k))) * 100
   
   ## leaf stomatal resistance
   if (!use_temir_beta) betat = 1 
   tmprs = r.s.zhang(rsmin, brs, par.sun, par.sha, Lsun, Lsha, T.c, tmin, tmax, topt, vpd, bvpd, srad, psi.min, psi.max, use.betat = use_temir_beta, betat)
   if (!use_temir_rs) {
      r_s = tmprs$r_s
      if (co2_scale & L > 0) r_s = r_s * (CO2_conc - 40)*440/(CO2_conc + 80)/320*CO2_conc/360 ### apply co2 effect
   } else {
      if (is.null(r_s)) stop('r_s has not been provided as input!')
      # Leave r_s as it is.
   }
   
   g_par = tmprs$g_par ## Zhang et al. conductance reducing factors
   g_w = tmprs$g_w
   g_t = tmprs$g_t
   g_vpd = tmprs$g_vpd
   
   r_s = r_s * diffg(T.k, p, 18E-3) / diffg(T.k, p, mx) # molecular diffusivities for water vapour and the pollutant gas
   r_s = min(r_s, 1E7)
   
   ## Cuticle resistance
   r_cut = r.cut.zhang(r.cut0, T.c, L, is.wet, rh, ustar) / hstar.frac(hstar, f0)
   
   ## In-canopy aerodynamic resistance
   ra_dc = r.ac(rac0_min, rac0_max, L, LAI_min, LAI_max, ustar)
   
   ## Ground resistance for O3
   rgx = rgo
   
   # Scale ground and cuticular resistance according to snow cover
   if (fsnow > 0) {
      rgs_snow = min(500, max(100, 70*(2-T.c)))
      r_snow = rgx(hstar, f0, rgs_snow, 2000) #?
      r_cut = 1 / ((1 - fsnow) / r_cut + fsnow / r_snow)
      rgx = 1 / ((1 - min(2*fsnow, 1)) / rgx + min(2*fsnow, 1) / r_snow)
   }
   
   # Prevent overestimation of r_cut for sparsely vegetated/snow-buried places
   LSAI = L + SAI
   r_cut = if (LSAI < 0.5) 1E7 else r_cut
   
   ## Aerodynamic resistance
   if (!use_temir_ra) {
      obk = obk(rho, T.s, ustar, H)
      r_a = r.a(ustar, z0, obk, T.k, cz)
   } else {
      if (is.null(r_a)) stop('r_a has not been provided as input!')
      # Take r_a as provided.
   }
   
   ## Sub-layer quasi-laminar resistance
   r_b = r.b(T.k, p, mx, ustar)
   
   ## Stomatal blocking factor
   if (is.wet) { ## for wet canopies
      if (srad < 200) {
         w.st = 0
      } else if (srad < 600) {
         w.st = (srad - 200)/800
      } else {
         w.st = 0.5
      }
   } else { ## for dry canopies
      w.st = 0
   }
   
   ### canopy conductance
   r_c = 1/((1-w.st)/r_s + 1/(ra_dc + rgx) + 1/r_cut)
   
   ### dry deposition velocity
   v_d = 1/(r_a + r_b + r_c)
   
   out = list(r_a = r_a, r_b = r_b, r_s = r_s, g_par = g_par, g_vpd = g_vpd, w_st = w.st, g_t = g_t, g_w = g_w, r_cut = r_cut, fsnow = fsnow, ra_dc = ra_dc, rgx = rgx, v_d = v_d)
   
   return(out)
}

###############################################################################
### End of module
###############################################################################
