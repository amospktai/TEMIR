################################################################################
### Module for the Farquhar-Ball-Berry model to calculate leaf and canopy photosynthesis and conductance
################################################################################

################################################################################
### Revision history
################################################################################

# Aug 2017: v1.0 finalized (Tai)
# Feb 2018: Definition of "A_can" is modified to be consistent with literature. It now means GROSS canopy photosynthesis, i.e., GPP. (Tai)
# Oct 2018: Added new "Medlyn" stomatal conductance scheme. (Sun)
# Feb 2019: Modified the way L_sun and L_sha are defined and calculated. We believe that A_can, R_can and g_can should be scaled up by LAI only, not by LAI + SAI. Therefore, L_sun and L_sha here should be sunlit and shaded leaf area index, not plant area index, but we still consider both LAI and SAI when calculating light extinction. Now they are renamed "LAI_sun" and "LAI_sha" and "LAI_sun" has to taken explicitly from canopy radiative transfer model. (Tai)
# Feb 2019: Now modified such that "g_ah" (aerodynamic conductance for heat, umol m^-2 s^-1 or m s^-1) is used instead of "g_am". (Tai, Feb 2019)
# Apr 2020: Added "f_Kn_Kb" to scale "R_d25" in f_leaf_photosyn(). Not sure if "R_d25" needs to be adjusted for daylength in this case. (Tai, Apr 2020)

################################################################################
### Functions:
################################################################################

# Roots of quadratic equation:
quadroot = function(a, b, c) {
    # This function finds the roots of quadratic equation a*x^2 + b*x + c = 0.
    D = b^2 - 4*a*c
    if (is.na(D)) {
        x1 = NaN
        x2 = NaN
    } else {
        if (a == 0) {
            # If a is exactly zero, simply solve for linear equation b*x + c = 0
            x1 = -c/b
            x2 = NaN
        } else {
            if (D >= 0) {
                # Roots are real.
                x1 = (-b+sqrt(D))/(2*a)
                x2 = (-b-sqrt(D))/(2*a)
            }
            else {
                # Roots are complex.
                x1 = complex(real=(-b/2/a), imaginary=(sqrt(-D)/2/a))
                x2 = complex(real=(-b/2/a), imaginary=(-sqrt(-D)/2/a))
            }
        }
    }
    return(c(x1, x2))
}

################################################################################

# Saturation water vapor pressure (Pa):
f_esat = function(T_K, derivative=FALSE) {
    # This function calculates the saturation water vapor pressure (Pa) for a given temperature "T_K" (K) or a vector/matrix of "T_K" (K) based on Lowe and Ficke (1974).
    T_C = T_K - 273.15			# Convert K to degC.
    a0 = 6.107799961
    a1 = 4.436518521e-1
    a2 = 1.428945805e-2
    a3 = 2.650648471e-4
    a4 = 3.031240396e-6
    a5 = 2.034080948e-8
    a6 = 6.136820929e-11
    if (derivative) {
        desat_dT = (a1 + 2*a2*T_C + 3*a3*T_C^2 + 4*a4*T_C^3 + 5*a5*T_C^4 + 6*a6*T_C^5)*100
        return(desat_dT)
    } else {
        esat = (a0 + a1*T_C + a2*T_C^2 + a3*T_C^3 + a4*T_C^4 + a5*T_C^5 + a6*T_C^6)*100
        return(esat)
    }
}

################################################################################

# Leaf boundary-layer conductance (umol m^-2 s^-1 or m s^-1):
f_boundary_cond = function(u_leaf, d_leaf=0.04, met_cond=FALSE, P_atm=101325, theta_atm=NULL) {
    
    # Uses global constants: "R_uni=8.31446" universal gas constant (J K^-1 mol^-1)
    # Wind speed incident on leaf "u_leaf" (m s^-1) must be specified. "u_leaf" is equivalent to the friction velocity in this formulation.
    # Characteristic leaf dimension is typically "d_leaf=0.04" (m).
    # If "met_cond=TRUE", conductance will be given in the common meteorological unit of m s^-1.
    # "P_atm" = atmospheric pressure at reference height (Pa); "theta_atm" =  atmospheric potential temperature at reference height (K).
    
    # Turbulent transfer coefficient (m s^-0.5):
    C_v = 0.01
    # Leaf boundary-layer conductance (m s^-1):
    g_b = C_v*(u_leaf/d_leaf)^0.5
    
    # Convert to stomatal conductance (umol m^-2 s^-1 ) if necessary:
    if (met_cond) {
        # Do nothing.
    } else {
        if (is.null(theta_atm)) stop('theta_atm needs to be specified.')
        mol_to_met = 1e-6*R_uni*theta_atm/P_atm
        g_b = g_b/mol_to_met
    }
    
    return(g_b)
}

################################################################################

# Leaf stomatal conductance (umol m^-2 s^-1 or m s^-1):
f_stomatal_cond = function(A_n, c_s=NULL, e_s=NULL, c_a=NULL, e_a=NULL, g_b=NULL, T_v=298.15, P_atm=101325, beta_t=1, C3_plant=TRUE, gs_scheme='FBB', met_cond=FALSE, theta_atm=NULL, g1_med=NULL) {
    
    # Uses external functions "f_esat" and "quadroot".
    # Uses global constants: "R_uni=8.31446" universal gas constant (J K^-1 mol^-1)
    # Net photosynthesis rate "A_n" (umol CO2 m^-2 s^-1) must be specified.
    # There are two ways to calculate stomatal conductance "g_s":
    # 1. If leaf surface vapor pressure "e_s" (Pa) and leaf surface CO2 partial pressure "c_s" (Pa) are given, it calculates "g_s" directly.
    # 2. If canopy air vapor pressure "e_a" (Pa), atmospheric CO2 partial pressure "c_a" (Pa) and leaf boundary-layer conductance "g_b" (umol m^-2 s^-1 or m s^-1) are given, it calculates "g_s" via solving a quadratic equation assuming A_n = (g_b/1.4)*(c_a - c_s)/P_atm.
    # "T_v" = vegetation temperature (K); "P_atm" = atmospheric pressure at reference height (Pa); "theta_atm" = atmospheric potential temperature at reference height (K); "beta_t" = soil water stress function (0 to 1)
    # "C3_plant=TRUE" specifies C3 photosynthesis; C4 otherwise.
    # If "met_cond=TRUE", conductance will be given in the common meteorological unit of m s^-1, and if "g_b" is specified it has to be in m s^-1 too.
    
    # "gs_scheme" sets the stomatal conductance scheme to be used. It can either be the conventional Ball-Berry model ("FBB"), or the newer Medlyn model ("Medlyn") (see below).
    
    # shsun 18/10/2017
    # Medlyn model to calculate stomatal conductance "g_s":
    # "e_vpd" (kPa) = water vapor pressure deficit (VPD) at leaf surface;
    # 1. If "e_s" (Pa) and "c_s" (Pa) are given, it calculates "g_s" directly;
    # 2. If "e_a" (Pa) adn "c_a" (Pa) and "g_b" (umol m^-2 s^-1 or m s^-1) are given, it calculates "g_s" via solving a quadratic equation using:
    #   e_s = (e_a/g_s+e_i/g_b)/(1/g_b+1/g_s)
    #   e_vpd = e_i-e_s
    #   g_s = 1.6*(1+m/e_vpd^0.5)*A_n*P_atm/c_s + b
    # Medlyn et al. 2011
    
    # Saturation vapor pressure at leaf temperature (Pa):
    e_sat = f_esat(T_v)
    
    if (met_cond) {
        if (is.null(theta_atm)) stop('theta_atm needs to be specified.')
        mol_to_met = 1e-6*R_uni*theta_atm/P_atm
    } else {
        # Do nothing.
    }
    
    if (gs_scheme == 'Medlyn') {
        m = g1_med*beta_t
        b = 100
        if (is.na(g1_med)) {  # if no parameter available, use previous C3/C4 param
            if (C3_plant) {
                m = 3.37*beta_t    # kPa^0.5
                b = 100            # umol m^-2 s^-1
            }else {
                m = 1.1*beta_t      # kPa^0.5
                b = 100             # umol m^-2 s^-1
            }
        }
        
        if (A_n <= 0) {
            g_s = b
        } else {
            if (!is.null(c_s) & !is.null(e_s)) {
               
                # put constraints on RH/vpd in MED mode
                e_vpd = max(e_sat - e_s, 50)*0.001
                g_s = 1.6*(1 + m/(e_vpd^0.5))*A_n/(c_s/P_atm) + b
                
            } else {
                
                if (is.null(c_a) | is.null(e_a) | is.null(g_b)) stop('All of c_a, e_a and g_b need to be specified.')
                if (met_cond) g_b = g_b/mol_to_met
                tmp_c_s = max(c(1e-6, (c_a - (1.4/g_b)*P_atm*A_n)), na.rm=TRUE)
                c_s = tmp_c_s/P_atm
                e_vpd = max(e_sat - e_a, 50)*0.001
                
                term = 1.6*A_n/c_s
                tmp_a = 1
                tmp_b = -(2*b+2*term+term^2*m^2/(g_b*e_vpd))
                tmp_c = b^2+2*b*term+term^2-m^2*term^2/e_vpd
                judge = tmp_b^2-4*tmp_a*tmp_c
                
                if (judge > 0) {
                    gs_roots = quadroot(a=tmp_a, b=tmp_b, c=tmp_c)
                    g_s = max(gs_roots, na.rm = TRUE)
                } else {
                    ## no solution for optimal conductance then switch to fbb scheme
                    if (C3_plant) {
                        m = 9                # unitless
                        b = 10000*beta_t     # umol m^-2 s^-1
                    } else {
                        m = 4                # unitless
                        b = 40000*beta_t     # umol m^-2 s^-1
                    }
                    c_s = max(c(1e-6, (c_a - (1.4/g_b)*P_atm*A_n)), na.rm=TRUE)
                    gs_roots = quadroot(a=c_s, b=(c_s*(g_b - b) - m*A_n*P_atm), c=(-g_b*(c_s*b + m*A_n*P_atm*e_a/e_sat)))
                    g_s = max(gs_roots, na.rm=TRUE)
                }
                
            }
        }
    } else if (gs_scheme == 'FBB') {
        if (C3_plant) {
            m = 9                # unitless
            b = 10000*beta_t     # umol m^-2 s^-1
        } else {
            m = 4                # unitless
            b = 40000*beta_t     # umol m^-2 s^-1
        }
        if (A_n <= 0) {
            g_s = b
        } else {
            # Two methods to calculate stomatal conductance (umol m^-2 s^-1):
            if (!is.null(c_s) & !is.null(e_s)) {
                g_s = m*A_n*(e_s/e_sat)/(c_s/P_atm) + b
            } else {
                if (is.null(c_a) | is.null(e_a) | is.null(g_b)) stop('All of c_a, e_a and g_b need to be specified.')
                if (met_cond) g_b = g_b/mol_to_met
                c_s = max(c(1e-6, (c_a - (1.4/g_b)*P_atm*A_n)), na.rm=TRUE)
                gs_roots = quadroot(a=c_s, b=(c_s*(g_b - b) - m*A_n*P_atm), c=(-g_b*(c_s*b + m*A_n*P_atm*e_a/e_sat)))
                g_s = max(gs_roots, na.rm=TRUE)
            }
        }
    } else {
       stop('Stomatal conductance scheme is not correctly specified!')
    }
    
    # Convert to stomatal conductance (m s^-1) if necessary:
    if (met_cond) g_s = g_s*mol_to_met
    return(g_s)
    
}

################################################################################

# Temperature scaling function 1 for photosynthetic parameters:
# T_v = vegetation temperature (K)
# DH_a (J mol^-1)
f_Tv = function(T_v, DH_a) exp(DH_a/(298.15*R_uni)*(1 - 298.15/T_v))

################################################################################

# Temperature scaling function 2 for photosynthetic parameters:
# T_v = vegetation temperature (K)
# DH_d (J mol^-1)
# DS (J mol^-1 K^-1)
f_H_Tv = function(T_v, DH_d, DS) (1 + exp((298.15*DS - DH_d)/(298.15*R_uni)))/(1 + exp((DS*T_v - DH_d)/(R_uni*T_v)))

################################################################################

# Leaf net photosynthesis rate (umol CO2 m^-2 s^-1):
f_leaf_photosyn = function(c_i, phi, T_v=298.15, P_atm=101325, C3_plant=TRUE, V_cmax25_0, Phi_PSII=0.85, Theta_PSII=0.70, beta_t=1, acclimation=FALSE, T_10d=NULL, lat=0, decl=0, colimit=TRUE, canopy_avg=FALSE, sunlit=NULL, LAI=NULL, LAI_sun=NULL, K_n=0.30, K_b=0.50, biogeochem=FALSE, leaf_N_conc=NULL) {
    
    # Uses external functions "f_Tv", "f_H_Tv" and "quadroot".
    # Uses global constants: "R_uni=8.31446" universal gas constant (J K^-1 mol^-1)
    # Intercellular CO2 partial pressure "c_i" (Pa), absorbed photosynthetically active radiation "phi" (W m^-2), and maximum rate of carboxylation at 25 degC (at top of canopy) "V_cmax25_0" (umol CO2 m^-2 s^-1) should always be provided.
    # "T_v" = vegetation temperature (K); "P_atm" = atmospheric pressure at reference height (Pa); "Phi_PSII=0.85" = quantum yield of photosystem II; "Theta_PSII=0.70" = curvature parameter; "beta_t" = soil water stress function (0 to 1)
    # "C3_plant=TRUE" specifies C3 photosynthesis; C4 otherwise.
    # If "colimit=TRUE", colimitation among A_c, A_j and A_p is considered.
    # If "acclimation=TRUE", temperature acclimation is considered, and "J_max" and "V_cmax" would vary with plant growth temperature ("T_10d" (K) = 10-day mean air temperature), and "lat" and "decl" are the latitude (rad) and solar declination angle (rad).
    # If "canopy_avg=TRUE", photosynthetic parameters are scaled according to canopy characteristics and the canopy-averaged values are calculated. Leaf and stem area indices, "LAI" and "SAI", should be provided. Rates for sunlit vs. shaded leaves are calculated separately, and "sunlit=TRUE/FALSE" has to be explicitly defined. Sunlit vs. shaded "phi" should be provided correspondingly.
   # *** Now modified such that if "canopy_avg=TRUE", "SAI" need not be provided, but the sunlit leaf area index, "LAI_sun", has to explicitly provided. (Tai, Feb 2019)
    # Canopy scaling uses the nitrogen and light extinction coefficient, "K_n" and "K_b", respectively.
    # If "biogeochem=TRUE", leaf mitochondrial respiration at 25 degC "R_d25" is calculated using leaf nitrogen concentration "leaf_N_conc" (g N m^-2 leaf area).
    
    # Canopy scaling:
    if (canopy_avg) {
        # if (is.null(LAI) | is.null(SAI) | is.null(sunlit)) stop('LAI, SAI or sunlit parameters have to be provided to calculate canopy averages.')
        if (is.null(LAI) | is.null(LAI_sun) | is.null(sunlit)) stop('LAI or sunlit parameters have to be provided to calculate canopy averages.')
        # Sunlit and shaded plant area index:
        # Make sure K_b > 0 always. At night, set K_b = 1e6.
        # We believe that A_can, R_can and g_can should be scaled up by LAI only, not by LAI + SAI. Therefore, L_sun and L_sha here should be sunlit and shaded leaf area index, not plant area index, but we still consider both LAI and SAI when calculating light extinction. Now they are renamed "LAI_sun" and "LAI_sha" and "LAI_sun" has to taken explicitly from canopy radiative transfer model. (Tai, Feb 2019)
        # L_sun = (1 - exp(-K_b*(LAI + SAI)))/K_b
        # L_sha = (LAI + SAI) - L_sun
        # Shaded leaf area index:
        LAI_sha = LAI - LAI_sun
        # Light and nitrogen extinction scaling functions:
        f_sun = (1 - exp(-(K_n + K_b)*LAI))/(K_n + K_b)
        if (sunlit) f_Kn_Kb = f_sun/LAI_sun else f_Kn_Kb = ((1 - exp(-K_n*LAI))/K_n - f_sun)/LAI_sha
        if (f_Kn_Kb == 'Inf') f_Kn_Kb = 0
    } else {
        f_Kn_Kb = 1
    }
    
    # Daylength (s):
    arg = -sin(lat)*sin(decl)/(cos(lat)*cos(decl))
    arg = min(c(max(c(-1, arg), na.rm=TRUE), 1), na.rm=TRUE)
    DYL = 2*13750.9871*acos(arg)
    # Maximum daylength (s):
    if (lat >= 0) decl_max = 0.409571 else decl_max = -0.409571
    arg = -sin(lat)*sin(decl_max)/(cos(lat)*cos(decl_max))
    arg = min(c(max(c(-1, arg), na.rm=TRUE), 1), na.rm=TRUE)
    DYL_max = 2*13750.9871*acos(arg)
    # Daylength adjustment function:
    f_DYL = min(c(max(c(0.01, DYL^2/DYL_max^2), na.rm=TRUE), 1), na.rm=TRUE)
    
    # Maximum rate of carboxylation at 25 degC adjusted for daylength (umol CO2 m^-2 s^-1):
    V_cmax25 = V_cmax25_0*f_Kn_Kb*f_DYL
    
    # Calculate leaf mitochondrial respiration both day and night:
    # Leaf mitochondrial respiration at 25 degC (umol CO2 m^-2 s^-1):
    if (biogeochem) {
        if (is.null(leaf_N_conc)) stop('leaf_N_conc has to be specified.')
        R_d25 = 0.2577*leaf_N_conc*f_Kn_Kb
        # Added "f_Kn_Kb" to scale "R_d25" above. Not sure if it needs to be adjusted for daylength. (Tai, Apr 2020)
    } else {
        if (C3_plant) R_d25 = 0.015*V_cmax25 else R_d25 = 0.025*V_cmax25
    }
    # Leaf mitochondrial respiration rate adjusted for soil water stress for C3 or C4 plants (umol CO2 m^-2 s^-1):
    if (C3_plant) {
        R_d = R_d25*f_Tv(T_v=T_v, DH_a=46390)*f_H_Tv(T_v=T_v, DH_d=150650, DS=490)*beta_t
    } else {
        R_d = R_d25*(2^((T_v - 298.15)/10))/(1 + exp(1.3*(T_v - 328.15)))*beta_t
    }
    R_d = max(c(R_d, 0), na.rm=TRUE)
    
    # Calculate leaf photosynthesis only for day time:
    if (phi <= 0) {
        # It is at night. Skip calculations.
        A_c = 0
        A_j = 0
        A_p = 0
        A = 0
    } else {
        
        # Michaelis-Menton constant for carboxylation at 25 degC (Pa):
        K_c25 = 404.9e-6*P_atm
        # Michaelis-Menton constant for oxygenation at 25 degC (Pa):
        K_o25 = 278.4e-3*P_atm
        # CO2 compensation point at 25 degC (Pa):
        Gamma_star25 = 42.75e-6*P_atm
        # Intercellular oxygen concentration (Pa):
        o_i = 0.209*P_atm
        
        # Product utilization rate at 25 degC (umol CO2 m^-2 s^-1):
        T_p25 = 0.167*V_cmax25
        # Initial slope of C4 CO2 response curve at 25 degC (Pa/Pa):
        k_p25 = 20000*V_cmax25
        
        # Temperature acclimation:
        if (acclimation) {
            if (is.null(T_10d)) stop('T_10d needs to be provided to calculate acclimation.')
            # Maximum potential electron transport rate at 25 degC(umol CO2 m^-2 s^-1):
            J_max25 = V_cmax25*(2.59 - 0.035*min(c(max(c((T_10d - 273.15), 11), na.rm=TRUE), 35), na.rm=TRUE))
            # Other temperature dependence parameters: DH_a (J mol^-1), DH_d (J mol^-1), DS (J mol^-1 K^-1)
            DH_a_Vcmax = 72000
            DH_a_Jmax = 50000
            DH_d_Vcmax = 200000
            DH_d_Jmax = 200000
            DS_Vcmax = 668.39 - 1.07*min(c(max(c((T_10d - 273.15), 11), na.rm=TRUE), 35), na.rm=TRUE)
            DS_Jmax = 659.70 - 0.75*min(c(max(c((T_10d - 273.15), 11), na.rm=TRUE), 35), na.rm=TRUE)
        } else {
            # Maximum potential electron transport rate at 25 degC(umol CO2 m^-2 s^-1):
            J_max25 = 1.97*V_cmax25
            # Other temperature dependence parameters: DH_a (J mol^-1), DH_d (J mol^-1), DS (J mol^-1 K^-1)
            DH_a_Vcmax = 65330
            DH_a_Jmax = 43540
            DH_d_Vcmax = 149250
            DH_d_Jmax = 152040
            DS_Vcmax = 485
            DS_Jmax = 495
        }
        
        if (C3_plant) {
            
            # Maximum rate of carboxylation adjusted for soil water stress (umol CO2 m^-2 s^-1):
            V_cmax = V_cmax25*f_Tv(T_v=T_v, DH_a=DH_a_Vcmax)*f_H_Tv(T_v=T_v, DH_d=DH_d_Vcmax, DS=DS_Vcmax)*beta_t
            # CO2 compensation point (Pa):
            Gamma_star = Gamma_star25*f_Tv(T_v=T_v, DH_a=37830)
            # Michaelis-Menton constant for carboxylation (Pa):
            K_c = K_c25*f_Tv(T_v=T_v, DH_a=79430)
            # Michaelis-Menton constant for oxygenation (Pa):
            K_o = K_o25*f_Tv(T_v=T_v, DH_a=36380)
            # Rubisco-limited photosynthesis rate (umol CO2 m^-2 s^-1):
            A_c = V_cmax*(c_i - Gamma_star)/(c_i + K_c*(1 + o_i/K_o))
            
            # Light utilized in electron transport by photosystem II (umol m^-2 s^-1):
            # I_PSII = 0.5*Phi_PSII*alpha_l*I_l
            I_PSII = 0.5*Phi_PSII*(4.6*phi)
            # Electron transport rate (umol m^-2 s^-1):
            # Maximum potential electron transport rate (for C3 plants only) (umol m^-2 s^-1):
            J_max = J_max25*f_Tv(T_v=T_v, DH_a=DH_a_Jmax)*f_H_Tv(T_v=T_v, DH_d=DH_d_Jmax, DS=DS_Jmax)
            J_roots = quadroot(a=Theta_PSII, b=-(I_PSII + J_max), c=I_PSII*J_max)
            J = min(J_roots, na.rm=TRUE)
            # CO2 compensation point (Pa):
            Gamma_star = Gamma_star25*f_Tv(T_v=T_v, DH_a=37830)
            # RuBP-limited photosynthesis rate (umol CO2 m^-2 s^-1):
            A_j = J*(c_i - Gamma_star)/(4*c_i + 8*Gamma_star)
            
            # Triose phosphate utilization rate (for C3 plants only) (umol m^-2 s^-1):
            T_p = T_p25*f_Tv(T_v=T_v, DH_a=DH_a_Vcmax)*f_H_Tv(T_v=T_v, DH_d=DH_d_Vcmax, DS=DS_Vcmax)
            # Product-limited photosynthesis rate (umol CO2 m^-2 s^-1):
            A_p = 3*T_p
            
        } else {
            
            # Maximum rate of carboxylation adjusted for soil water stress (umol CO2 m^-2 s^-1):
            V_cmax = V_cmax25*(2^((T_v - 298.15)/10))/((1 + exp(0.3*(T_v - 313.15)))*(1 + exp(0.2*(288.15 - T_v))))*beta_t
            # Rubisco-limited photosynthesis rate (umol CO2 m^-2 s^-1):
            A_c = V_cmax
            
            # RuBP-limited photosynthesis rate (umol CO2 m^-2 s^-1):
            A_j = 0.05*(4.6*phi)
            
            # Initial slope of C4 CO2 response curve (Pa/Pa):
            k_p = k_p25*2^((T_v - 298.15)/10)
            # Product-limited photosynthesis rate (umol CO2 m^-2 s^-1):
            A_p = k_p*c_i/P_atm
            
        }
        
        A_c = max(c(A_c, 0), na.rm=TRUE)
        A_j = max(c(A_j, 0), na.rm=TRUE)
        A_p = max(c(A_p, 0), na.rm=TRUE)
        
        # Gross photosynthesis rate (umol CO2 m^-2 s^-1):
        if (colimit) {
            if (C3_plant) Theta_cj = 0.98 else Theta_cj = 0.80
            Theta_ip = 0.95
            Ai_roots = quadroot(a=Theta_cj, b=-(A_c + A_j), c=A_c*A_j)
            A_i = min(Ai_roots, na.rm=TRUE)
            A_roots = quadroot(a=Theta_ip, b=-(A_i + A_p), c=A_i*A_p)
            A = min(A_roots, na.rm=TRUE)
        } else {
            A = min(c(A_c, A_j, A_p), na.rm=TRUE)
        }
        
    }
    
    # Net photosynthesis rate (umol CO2 m^-2 s^-1):
    A_n = A - R_d
    
    output = list(A_c=A_c, A_j=A_j, A_p=A_p, A=A, R_d=R_d, A_n=A_n)
    return(output)
    
}

################################################################################

# Ozone impact factors to adjust photosynthesis and/or stomatal conductance (0 to 1):
f_ozone_impact = function(O3_conc, g_s, g_b, g_ah, P_atm=101325, theta_atm=298.15, met_cond=TRUE, scheme='Lombardozzi', dt=3600, LAI=0.5, LAI_prev=0.5, evergreen=FALSE, leaf_long=1e3, CUO_prev=0, sensitivity='high', plant_group='broadleaf') {
    
    # Uses global constants: "R_uni=8.31446" universal gas constant (J K^-1 mol^-1)
    # Atmospheric ozone concentration "O3_conc" (ppb), stomatal "g_s", leaf boundary-layer "g_b" and aerodynamic "g_am" conductances (for momentum) (umol m^-2 s^-1 or m s^-1) have to be always specified.
    # *** Now modified such that "g_ah" (aerodynamic conductance for heat, umol m^-2 s^-1 or m s^-1) is used instead. (Tai, Feb 2019)
    # If "met_cond=TRUE", conductances are specified in m s^-1 instead of umol m^-2 s^-1.
    # "P_atm" = atmospheric pressure at reference height (Pa); "theta_atm" =  atmospheric potential temperature at reference height (K).
    # "scheme" should either be set to 'Lombardozzi' (Lombardozzi et al., 2015) or 'Sitch' (Sitch et al., 2007).
    # If "scheme='Lombardozzi'", "dt" = model time step (s); "LAI" = leaf area index at current time step; "LAI_prev" = LAI at previous time step; "leaf_long" = leaf longevity (yr); and "CUO_prev" = cumulative uptake of ozone until previous time step (mmol m^-2).
    # If "evergreen=TRUE", plant type considered is evergreen.
    # "sensitivity" should either be set to "high", "average" or "low" for the Lombardozzi scheme, or just "high" or "low" for the Sitch scheme.
    # "plant_group" should be one of "broadleaf", "needleleaf", "C3_grass", "C4_grass" or "shrub".
    
    if (met_cond) {
        # Do nothing.
    } else {
        # Convert from umol m^-2 s^-1 to m s^-1:
        mol_to_met = 1e-6*R_uni*theta_atm/P_atm
        g_s = g_s*mol_to_met
        g_b = g_b*mol_to_met
        # g_am = g_am*mol_to_met
        g_ah = g_ah*mol_to_met
    }
    
    # Convert from ppb to nmol m^-3:
    O3_conc_nmol = O3_conc*P_atm/(theta_atm*R_uni)
    # O3:H2O resistance ratio defined by Sitch et al. (2007):
    k_O3 = 1.67
    # Instantaneous ozone flux (nmol m^-2 s^-1):
    O3_flux = O3_conc_nmol/(k_O3/g_s + 1/g_b + 1/g_ah)
    # This follows Lombardozzi et al. [2015], but it seems more appropriate to use g_ah or g_aw instead of g_am to calculate ozone flux.
    # Now modified such that "g_ah" is used instead. (Tai, Feb 2019)
    
    if (scheme == 'Lombardozzi') {
        
        # Add a threshold of uptake (nmol m^-2 s^-1):
        O3_flux_th = 0.8
        O3_flux_crit = max(c(0, (O3_flux - O3_flux_th)), na.rm=TRUE)
        # Model time step in hr:
        dt_hr = dt/3600
        # Ozone flux per time step (mmol m^-2):
        O3_flux_dt = O3_flux_crit*1e-6*dt
        
        if (LAI > 0.5) {
            # Minimize O3 damage to new leaves:
            if (LAI > LAI_prev) heal = max(c(0, (LAI - LAI_prev)/LAI*O3_flux_dt), na.rm=TRUE) else heal = 0
            # Leaf turnover in hr^-1 by PFT:
            if (evergreen) turnover = 1/(leaf_long*365*24) else turnover = 0
            # O3 uptake decay based on leaf lifetime for evergreen plants:
            decay = CUO_prev*turnover*dt_hr
            # New CUO (mmol m^-2):
            CUO = max(c(0, (CUO_prev + O3_flux_dt - decay - heal)), na.rm=TRUE)
        } else {
            CUO = 0
        }
        
        # PFT-specific ozone damage parameters:
        if (sensitivity == 'high') {
            needle_An_int = 1.083
            needle_An_slope = -0.038
            broad_An_int = 0.8502
            broad_An_slope = 0
            grass_An_int = 0.8564
            grass_An_slope = -0.0007
            needle_gs_int = 0.8874
            needle_gs_slope = -0.0144
            broad_gs_int = 0.8900
            broad_gs_slope = 0
            grass_gs_int = 0.7074
            grass_gs_slope = 0
        } else if (sensitivity == 'low') {
            needle_An_int = 0.8595
            needle_An_slope = 0
            broad_An_int = 0.9798
            broad_An_slope = 0
            grass_An_int = 0.7159
            grass_An_slope = 0
            needle_gs_int = 0.7574
            needle_gs_slope = 0.0067
            broad_gs_int = 0.9425
            broad_gs_slope = 0
            grass_gs_int = 0.4621
            grass_gs_slope = 0.0229
        } else {
            needle_An_int = 0.8390
            needle_An_slope = 0
            broad_An_int = 0.8752
            broad_An_slope = 0
            grass_An_int = 0.8021
            grass_An_slope = -0.0009
            needle_gs_int = 0.7823
            needle_gs_slope = 0.0048
            broad_gs_int = 0.9125
            broad_gs_slope = 0
            grass_gs_int = 0.7511
            grass_gs_slope = 0
        }
        
        # Ozone impact factors:
        if (CUO == 0) {
            O3_coef_An = 1
            O3_coef_gs = 1
        } else {
            if (plant_group == 'C3_grass' | plant_group == 'C4_grass' ) {
                O3_coef_An = max(c(0, min(c((grass_An_int + grass_An_slope*CUO), 1), na.rm=TRUE)), na.rm=TRUE)
                O3_coef_gs = max(c(0, min(c((grass_gs_int + grass_gs_slope*CUO), 1), na.rm=TRUE)), na.rm=TRUE)
            } else if (plant_group == 'broadleaf' | plant_group == 'shrub') {
                O3_coef_An = max(c(0, min(c((broad_An_int + broad_An_slope*CUO), 1), na.rm=TRUE)), na.rm=TRUE)
                O3_coef_gs = max(c(0, min(c((broad_gs_int + broad_gs_slope*CUO), 1), na.rm=TRUE)), na.rm=TRUE)
            } else if (plant_group == 'needleleaf') {
                O3_coef_An = max(c(0, min(c((needle_An_int + needle_An_slope*CUO), 1), na.rm=TRUE)), na.rm=TRUE)
                O3_coef_gs = max(c(0, min(c((needle_gs_int + needle_gs_slope*CUO), 1), na.rm=TRUE)), na.rm=TRUE)
            } else {
                stop('plant_group is not correctly specified!')
            }
        }
        
    } else if (scheme == 'Sitch') {
        
        # Add a threshold of uptake (nmol m^-2 s^-1):
        if (plant_group == 'C3_grass' | plant_group == 'C4_grass') O3_flux_th = 5.0 else O3_flux_th = 1.6
        O3_flux_crit = max(c(0, (O3_flux - O3_flux_th)), na.rm=TRUE)
        
        # Ozone sensitivity parameter:
        a = 0
        if (sensitivity == 'high') {
            if (plant_group == 'broadleaf') a = 0.15
            if (plant_group == 'needleleaf') a = 0.075
            if (plant_group == 'C3_grass') a = 1.40
            if (plant_group == 'C4_grass') a = 0.735
            if (plant_group == 'shrub') a = 0.10
            if (a == 0) stop('plant_group is not correctly specified!')
        } else if (sensitivity == 'low') {
            if (plant_group == 'broadleaf') a = 0.04
            if (plant_group == 'needleleaf') a = 0.02
            if (plant_group == 'C3_grass') a = 0.25
            if (plant_group == 'C4_grass') a = 0.13
            if (plant_group == 'shrub') a = 0.03
            if (a == 0) stop('plant_group is not correctly specified!')
        } else {
            stop('sensitivity is not correctly specified!')
        }
        
        # Ozone impact factors:
        O3_coef_An = max(c(0, (1 - a*O3_flux_crit)), na.rm=TRUE)
        O3_coef_gs = 1
        CUO = 0
        
    } else {
        stop('scheme is not correctly specified!')
    }
    
    output = list(O3_coef_An=O3_coef_An, O3_coef_gs=O3_coef_gs, CUO=CUO)
    return(output)
    
}

################################################################################

# Iterative function for updating intercellular CO2 partial pressure "c_i" (for any situation except when "Sitch" ozone damage scheme is used:
f_ci = function() {
    # Used only internally in "f_Farquhar_Ball_Berry", and all arguments for functions "f_leaf_photosyn", and "f_stomatal_cond" have to be specified therein.
    # "c_i" would be updated at every step.
    # Net photosynthesis rate (umol m^-2 s^-1):
    leaf_photosyn = f_leaf_photosyn(c_i=c_i, phi=phi, T_v=T_v, P_atm=P_atm, C3_plant=C3_plant, V_cmax25_0=V_cmax25_0, Phi_PSII=Phi_PSII, Theta_PSII=Theta_PSII, beta_t=beta_t, acclimation=acclimation, T_10d=T_10d, lat=lat, decl=decl, colimit=colimit, canopy_avg=canopy_avg, sunlit=sunlit, LAI=LAI, LAI_sun=LAI_sun, K_n=K_n, K_b=K_b, biogeochem=biogeochem, leaf_N_conc=leaf_N_conc)
    A_n = leaf_photosyn$A_n
    R_d = leaf_photosyn$R_d
    # Leaf stomatal conductance (umol m^-2 s^-1):
    g_s = f_stomatal_cond(A_n=A_n, c_a=c_a, e_a=e_a, g_b=g_b, T_v=T_v, P_atm=P_atm, beta_t=beta_t, C3_plant=C3_plant, met_cond=FALSE, g1_med=g1_med)
    # Revised "c_i":
    ci_new = c_a - (1.4/g_b + 1.6/g_s)*P_atm*A_n
    if (leaf_photosyn$A <= 0) fval = 0 else fval = c_i - ci_new
    output = list(fval=fval, c_i=ci_new, g_s=g_s, A_n=A_n, R_d=R_d)
    return(output)
}

################################################################################

# Iterative function for updating intercellular CO2 partial pressure "c_i" and the associated ozone damage factor for the "Sitch" scheme:
# Would only be used if "O3_damage=TRUE" and "scheme='Sitch'".
f_ci2 = function() {
    # Used only internally in "f_Farquhar_Ball_Berry", and all arguments for functions "f_leaf_photosyn", "f_stomatal_cond" and "f_ozone_impact" have to be specified therein.
    # Both "c_i" and "O3_coef_An" would be updated at every step.
    # Net photosynthesis rate (umol m^-2 s^-1):
    leaf_photosyn = f_leaf_photosyn(c_i=c_i, phi=phi, T_v=T_v, P_atm=P_atm, C3_plant=C3_plant, V_cmax25_0=V_cmax25_0, Phi_PSII=Phi_PSII, Theta_PSII=Theta_PSII, beta_t=beta_t, acclimation=acclimation, T_10d=T_10d, lat=lat, decl=decl, colimit=colimit, canopy_avg=canopy_avg, sunlit=sunlit, LAI=LAI, LAI_sun=LAI_sun, K_n=K_n, K_b=K_b, biogeochem=biogeochem, leaf_N_conc=leaf_N_conc)
    A_n = leaf_photosyn$A_n
    R_d = leaf_photosyn$R_d
    # With ozone damage ("O3_coef_An" = 1 always if "O3_damage" is not turned on or "Sitch" scheme is not used):
    # Assuming that ozone damage only affects A_n when A_n > 0.
    if (A_n > 0) A_n = A_n*max(c(0, min(c(O3_coef_An, 1), na.rm=TRUE)), na.rm=TRUE)
    # R_d = R_d*max(c(0, min(c(O3_coef_An, 1), na.rm=TRUE)), na.rm=TRUE)
    # Leaf stomatal conductance (umol m^-2 s^-1):
    g_s = f_stomatal_cond(A_n=A_n, c_a=c_a, e_a=e_a, g_b=g_b, T_v=T_v, P_atm=P_atm, beta_t=beta_t, C3_plant=C3_plant, met_cond=FALSE, g1_med=g1_med)
    # Ozone impact factor on photosynthesis:
    ozone_impact = f_ozone_impact(O3_conc=O3_conc, g_s=g_s, g_b=g_b, g_ah=g_ah, P_atm=P_atm, theta_atm=theta_atm, met_cond=FALSE, scheme='Sitch', sensitivity=sensitivity, plant_group=plant_group)
    O3_coef_An_new = ozone_impact$O3_coef_An
    # Revised "c_i":
    ci_new = c_a - (1.4/g_b + 1.6/g_s)*P_atm*A_n
    if (leaf_photosyn$A <= 0) fval = 0 else fval = c_i - ci_new
    if (leaf_photosyn$A <= 0) gval = 0 else gval = O3_coef_An - O3_coef_An_new
    output = list(fval=fval, gval=gval, c_i=ci_new, O3_coef_An=O3_coef_An_new, g_s=g_s, A_n=A_n, R_d=R_d)
    return(output)
}


################################################################################

# Solving for leaf photosynthesis (umol CO2 m^-2 s^-1), stomatal conductance and boundary-layer conductance (umol m^-2 s^-1 or m s^-1) with the Farquhar-Ball-Berry model:
f_Farquhar_Ball_Berry = function(c_a, e_a, phi, T_v=298.15, P_atm=101325, theta_atm=298.15, C3_plant=TRUE, V_cmax25_0, Phi_PSII=0.85, Theta_PSII=0.70, beta_t=1, acclimation=FALSE, T_10d=NULL, lat=0, decl=0, colimit=TRUE, canopy_avg=FALSE, sunlit=NULL, LAI=NULL, LAI_sun=NULL, K_n=0.30, K_b=0.50, O3_damage=FALSE, O3_conc=NULL, scheme='Lombardozzi', gs_scheme = 'FBB', dt=3600, LAI_prev=NULL, evergreen=FALSE, leaf_long=1e3, CUO_prev=0, sensitivity='high', plant_group='broadleaf', g_ah=NULL, biogeochem=FALSE, leaf_N_conc=NULL, u_leaf, d_leaf=0.04, met_cond=FALSE, tol=1e-3, g1_med=NULL) {
    
    # Uses external functions "f_boundary_cond", "f_leaf_photosyn", "f_stomatal_cond", "f_ozone_impact", and "f_ci".
    # Uses global constants: "R_uni = 8.31446" universal gas constant (J K^-1 mol^-1)
    # Atmospheric CO2 partial pressure "c_a" (Pa), vapor pressure in canopy air "e_a" (Pa), absorbed photosynthetically active radiation "phi" (W m^-2), maximum rate of carboxylation at 25 degC (at top of canopy) "V_cmax25_0" (umol CO2 m^-2 s^-1), "sunlit=TRUE/FALSE", leaf area index "LAI" (m^2 m^-2), and wind speed incident on leaf "u_leaf" (m s^-1) must be specified.
    # "T_v" = vegetation temperature (K); "P_atm" = atmospheric pressure at reference height (Pa); "theta_atm" = atmospheric potential temperature at reference height (K); "Phi_PSII=0.85" = quantum yield of photosystem II; "Theta_PSII=0.70" = curvature parameter; "beta_t" = soil water stress function (0 to 1)
    # "C3_plant=TRUE" specifies C3 photosynthesis; C4 otherwise.
    # If "colimit=TRUE", colimitation among A_c, A_j and A_p is considered.
    # If "acclimation=TRUE", temperature acclimation is considered, and "J_max" and "V_cmax" would vary with plant growth temperature ("T_10d" (K) = 10-day mean air temperature), and "lat" and "decl" are the latitude (rad) and solar declination angle (rad).
    # If "canopy_avg=TRUE", photosynthetic parameters are scaled according to canopy characteristics and the canopy-averaged values are calculated. Leaf and stem area indices, "LAI" and "SAI", should be provided. Rates for sunlit vs. shaded leaves are calculated separately, and "sunlit=TRUE/FALSE" has to be explicitly defined. Sunlit vs. shaded "phi" should be provided correspondingly.
    # *** Now modified such that if "canopy_avg=TRUE", "SAI" need not be provided, but the sunlit leaf area index, "LAI_sun", has to explicitly provided. (Tai, Feb 2019)
    # Canopy scaling uses the nitrogen and light extinction coefficient, "K_n" and "K_b", respectively.
    # If "O3_damage=TRUE", ozone damage on photosynthesis and stomatal conductance is considered. When it is turned on, "O3_conc" = atmospheric ozone concentration (ppb); "scheme" should either be set to 'Lombardozzi' (Lombardozzi et al., 2015) or 'Sitch' (Sitch et al., 2007).
    # If "scheme='Lombardozzi'", "dt" = model time step (s); "LAI" = leaf area index at current time step; "LAI_prev" = LAI at previous time step; "leaf_long" = leaf longevity (yr); and "CUO_prev" = cumulative uptake of ozone until previous time step (mmol m^-2).
    # If "evergreen=TRUE", plant type considered is evergreen.
    # "sensitivity" should either be set to "high", "average" or "low" for the Lombardozzi scheme, or just "high" or "low" for the Sitch scheme.
    # "plant_group" should either be set to "broadleaf", "needleleaf" or "grass" for the Lombardozzi scheme, or "broadleaf", "needleleaf", "C3_grass", "C4_grass" or "shrub" for the Sitch scheme.
    # "u_atm" = wind speed at reference height adjusted for stability (m)
    # *** Now modified such that "u_atm" is not needed anymore because we directly read in "g_ah" (aerodynamic conductance for heat, umol m^-2 s^-1 or m s^-1) as input. (Tai, Feb 2019)
    # If "biogeochem=TRUE", leaf mitochondrial respiration at 25 degC "R_d25" is calculated using leaf nitrogen concentration "leaf_N_conc" (g N m^-2 leaf area).
    # Characteristic leaf dimension is typically "d_leaf=0.04" (m).
    # If "met_cond=TRUE", conductance will be given in the common meteorological unit of m s^-1.
    # "tol" specifies the absolute difference smaller than which iteration can stop.
    
    # Leaf boundary-layer conductance (umol m^-2 s^-1):
    g_b = f_boundary_cond(u_leaf=u_leaf, d_leaf=d_leaf, met_cond=FALSE, P_atm=P_atm, theta_atm=theta_atm)
    
    if (O3_damage) {
        
        # # Aerodynamic conductance (for momentum) (m s^-1):
        # # Assume that "u_leaf" is equivalent to friction velocity "u_star":
        # g_am = u_leaf^2/u_atm
        # # Aerodynamic conductance (for momentum) (umol m^-2 s^-1):
        # mol_to_met = 1e-6*R_uni*theta_atm/P_atm
        # g_am = g_am/mol_to_met
        
        # Now modified such that "g_ah" (aerodynamic conductance for heat) is used instead, which is an input now. (Tai, Feb 2019)
        mol_to_met = 1e-6*R_uni*theta_atm/P_atm
        # Aerodynamic conductance (for heat) (umol m^-2 s^-1) (conversion is needed because in f_ci() and f_ci2() conductances are all in molar unit):
        if (is.null(g_ah)) stop('g_ah is not provided!')
        if (met_cond) g_ah = g_ah/mol_to_met
        
    }
    
    if (O3_damage & scheme == 'Sitch') {
        
        environment(f_ci2) = environment()
        # Iterations:
        # Use bivariate secant updating method...
        # Initial intercellular CO2 partial pressure guess 0 (Pa):
        if (C3_plant) c_i0 = 0.7*c_a else c_i0 = 0.5*c_a
        O3_coef_An0 = 1
        c_i = c_i0
        O3_coef_An = O3_coef_An0
        Fval = f_ci2()
        fval0 = Fval$fval
        gval0 = Fval$gval
        n_itr = 0
        if (fval0 == 0 & gval0 == 0) {
            # A_n <= 0: g_s = b*beta_t and R_d as calculated.
        } else {
            # Start iterations:
            B0 = cbind(c(1, 0), c(0, 1))
            repeat {
                x0 = c(c_i0, O3_coef_An0)
                Fval0 = c(fval0, gval0)
                s0 = solve(B0, -Fval0)
                x1 = x0 + s0
                c_i1 = x1[1]
                O3_coef_An1 = x1[2]
                c_i = c_i1
                O3_coef_An = O3_coef_An1
                Fval = f_ci2()
                fval1 = Fval$fval
                gval1 = Fval$gval
                if (fval1 == 0 & gval1 == 0) break
                Fval1 = c(fval1, gval1)
                y0 = Fval1 - Fval0
                B1 = B0 + (y0 - B0%*%s0)%*%s0/as.vector(s0%*%s0)
                if (abs(c_i1 - c_i0) < tol & abs(O3_coef_An1 - O3_coef_An0) < tol*0.1) break
                if (n_itr > 100) {
                    print('Too many iterations!')
                    # Take initial value instead...
                    if (C3_plant) c_i = 0.7*c_a else c_i = 0.5*c_a
                    O3_coef_An = 1
                    Fval = f_ci2()
                    break
                }
                c_i0 = c_i1
                O3_coef_An0 = O3_coef_An1
                fval0 = fval1
                gval0 = gval1
                B0 = B1
                n_itr = n_itr + 1
            }
        }
        # Final values for "A_n", "R_d" and "g_s":
        A_n = Fval$A_n
        R_d = Fval$R_d
        g_s = Fval$g_s
        
    } else {
        
        environment(f_ci) = environment()
        # Iterations:
        # Use univariate secant method...
        # Initial intercellular CO2 partial pressure guess 0 (Pa):
        if (C3_plant) c_i = 0.7*c_a else c_i = 0.5*c_a
        fval = f_ci()
        c_i0 = c_i
        fval0 = fval$fval
        n_itr = 0
        if (fval0 == 0) {
            # A_n <= 0: g_s = b*beta_t and R_d as calculated.
        } else {
            # Initial intercellular CO2 partial pressure guess 1 (Pa):
            if (C3_plant) c_i = 0.7*c_a*0.99 else c_i = 0.5*c_a*0.99
            fval = f_ci()
            c_i1 = c_i
            fval1 = fval$fval
            if (!is.na(fval1) & fval1 == 0) {
                # A_n <= 0: g_s = b*beta_t and R_d as calculated.
            } else {
                # Start iterations:
                # # Debug only:
                # ci_vec = c_i
                repeat {
                    n_itr = n_itr + 1
                    # Step 0:
                    c_i = c_i1 - fval1*(c_i1 - c_i0)/(fval1 - fval0)
                    fval = f_ci()
                    c_i0 = c_i
                    fval0 = fval$fval
                    if (fval0 == 0) break
                    # Step 1:
                    c_i = c_i0 - fval0*(c_i0 - c_i1)/(fval0 - fval1)
                    fval = f_ci()
                    c_i1 = c_i
                    fval1 = fval$fval
                    if (fval1 == 0) break
                    if (abs(c_i1 - c_i0) < tol) break
                    if (n_itr > 100) {
                        print('Too many iterations!')
                        # Take initial value instead...
                        if (C3_plant) c_i = 0.7*c_a else c_i = 0.5*c_a
                        fval = f_ci()
                        break
                    }
                    # # Debug only:
                    # ci_vec = c(ci_vec, c_i)
                }
            }
        }
        # Final values for "A_n" and "g_s":
        A_n = fval$A_n
        R_d = fval$R_d
        g_s = fval$g_s
        
    }
    
    # CO2 partial pressure at leaf surface (Pa):
    c_s = max(c(1e-6, (c_a - (1.4/g_b)*P_atm*A_n)), na.rm=TRUE)
    
    if (O3_damage & scheme == 'Lombardozzi') {
        ozone_impact = f_ozone_impact(O3_conc=O3_conc, g_s=g_s, g_b=g_b, g_ah=g_ah, P_atm=P_atm, theta_atm=theta_atm, met_cond=FALSE, scheme='Lombardozzi', dt=dt, LAI=LAI, LAI_prev=LAI_prev, evergreen=evergreen, leaf_long=leaf_long, CUO_prev=CUO_prev, sensitivity=sensitivity, plant_group=plant_group)
        O3_coef_An = ozone_impact$O3_coef_An
        O3_coef_gs = ozone_impact$O3_coef_gs
        CUO = ozone_impact$CUO
        A_n = A_n*O3_coef_An
        g_s = g_s*O3_coef_gs
    } else if (O3_damage & scheme == 'Sitch') {
        # A_n and g_s taken as they are.
        CUO = 0
        # O3_coef_An take as it is.
        O3_coef_gs = 1
        # But run f_ozone_impact() again to obtain O3 flux:
        # Will write soon... (Tai, Feb 2019)
    } else {
        CUO = 0
        O3_coef_An = 1
        O3_coef_gs = 1
    }
    
    if (met_cond) {
        mol_to_met = 1e-6*R_uni*theta_atm/P_atm
        g_b = g_b*mol_to_met
        g_s = g_s*mol_to_met
    } else {
        # Do nothing.
    }
    
    output = list(c_i=c_i, c_s=c_s, g_b=g_b, g_s=g_s, A_n=A_n, R_d=R_d, CUO=CUO, O3_coef_An=O3_coef_An, O3_coef_gs=O3_coef_gs, n_itr=n_itr)
    return(output)
    
}


################################################################################

# Find canopy-integrated photosynthesis (umol CO2 m^-2 s^-1) and conductance (umol m^-2 s^-1 or m s^-1):
f_canopy_photosyn = function(c_a, e_a, phi_sun, phi_sha, T_v=298.15, P_atm=101325, theta_atm=298.15, C3_plant=TRUE, V_cmax25_0, Phi_PSII=0.85, Theta_PSII=0.70, beta_t=1, acclimation=TRUE, T_10d=298.15, lat=0, decl=0, colimit=TRUE, LAI, LAI_sun, K_n=0.30, K_b=0.50, O3_damage=FALSE, O3_conc=NULL, scheme='Lombardozzi', gs_scheme = 'FBB', dt=3600, LAI_prev=NULL, evergreen=FALSE, leaf_long=1e3, CUO_prev_sun=0, CUO_prev_sha=0, sensitivity='high', plant_group='broadleaf', g_ah=NULL, biogeochem=FALSE, leaf_N_conc=NULL, u_leaf, d_leaf=0.04, met_cond=FALSE, tol=1e-3,g1_med=NULL) {
    
    # Uses external function "f_Farquhar_Ball_Berry".
    # Uses global constants: "R_uni = 8.31446" universal gas constant (J K^-1 mol^-1)
    # Atmospheric CO2 partial pressure "c_a" (Pa), vapor pressure in canopy air "e_a" (Pa), absorbed photosynthetically active radiation by sunlit and shaded leaves, "phi_sun" and "phi_sha" (W m^-2), maximum rate of carboxylation at 25 degC (at top of canopy) "V_cmax25_0" (umol CO2 m^-2 s^-1), leaf area index "LAI" (m^2 m^-2), and wind speed incident on leaf "u_leaf" (m s^-1) must be specified.
    # "T_v" = vegetation temperature (K); "P_atm" = atmospheric pressure at reference height (Pa); "theta_atm" = atmospheric potential temperature at reference height (K); "Phi_PSII=0.85" = quantum yield of photosystem II; "Theta_PSII=0.70" = curvature parameter; "beta_t" = soil water stress function (0 to 1)
    # "C3_plant=TRUE" specifies C3 photosynthesis; C4 otherwise.
    # If "colimit=TRUE", colimitation among A_c, A_j and A_p is considered.
    # If "acclimation=TRUE", temperature acclimation is considered, and "J_max" and "V_cmax" would vary with plant growth temperature ("T_10d" (K) = 10-day mean air temperature), and "lat" and "decl" are the latitude (rad) and solar declination angle (rad).
    # Canopy scaling uses the nitrogen and light extinction coefficient, "K_n" and "K_b", respectively.
    # If "biogeochem=TRUE", leaf mitochondrial respiration at 25 degC "R_d25" is calculated using leaf nitrogen concentration "leaf_N_conc" (g N m^-2 leaf area).
    # Characteristic leaf dimension is typically "d_leaf=0.04" (m).
    # If "met_cond=TRUE", conductance will be given in the common meteorological unit of m s^-1.
   # *** Now modified such that "u_atm" is not needed anymore because we directly read in "g_ah" (aerodynamic conductance for heat, umol m^-2 s^-1 or m s^-1) as input. (Tai, Feb 2019)
    
    # Sunlit leaves:
    FBB_sun = f_Farquhar_Ball_Berry(c_a=c_a, e_a=e_a, phi=phi_sun, T_v=T_v, 
                                    P_atm=P_atm, theta_atm=theta_atm, 
                                    C3_plant=C3_plant, V_cmax25_0=V_cmax25_0, 
                                    Phi_PSII=Phi_PSII, Theta_PSII=Theta_PSII, 
                                    beta_t=beta_t, acclimation=acclimation, 
                                    T_10d=T_10d, lat=lat, decl=decl, 
                                    colimit=colimit, canopy_avg=TRUE, 
                                    sunlit=TRUE, LAI=LAI, LAI_sun=LAI_sun, 
                                    K_n=K_n, K_b=K_b, O3_damage=O3_damage, 
                                    O3_conc=O3_conc, scheme=scheme, 
                                    gs_scheme=gs_scheme, dt=dt, 
                                    LAI_prev=LAI_prev, evergreen=evergreen, 
                                    leaf_long=leaf_long, CUO_prev=CUO_prev_sun, 
                                    sensitivity=sensitivity, 
                                    plant_group=plant_group, 
                                    g_ah=g_ah, biogeochem=biogeochem, 
                                    leaf_N_conc=leaf_N_conc, u_leaf=u_leaf, 
                                    d_leaf=d_leaf, met_cond=met_cond, tol=tol,g1_med=g1_med)
    A_nsun = FBB_sun$A_n
    R_dsun = FBB_sun$R_d
    g_ssun = FBB_sun$g_s
    g_b = FBB_sun$g_b
    CUO_sun = FBB_sun$CUO
    
    # Shaded leaves:
    FBB_sha = f_Farquhar_Ball_Berry(c_a=c_a, e_a=e_a, phi=phi_sha, T_v=T_v, 
                                    P_atm=P_atm, theta_atm=theta_atm, 
                                    C3_plant=C3_plant, V_cmax25_0=V_cmax25_0, 
                                    Phi_PSII=Phi_PSII, Theta_PSII=Theta_PSII, 
                                    beta_t=beta_t, acclimation=acclimation, 
                                    T_10d=T_10d, lat=lat, decl=decl, 
                                    colimit=colimit, canopy_avg=TRUE, 
                                    sunlit=FALSE, LAI=LAI, LAI_sun=LAI_sun, 
                                    K_n=K_n, K_b=K_b, O3_damage=O3_damage, 
                                    O3_conc=O3_conc, scheme=scheme, 
                                    gs_scheme=gs_scheme, dt=dt, 
                                    LAI_prev=LAI_prev, evergreen=evergreen, 
                                    leaf_long=leaf_long, CUO_prev=CUO_prev_sha, 
                                    sensitivity=sensitivity, 
                                    plant_group=plant_group, 
                                    g_ah=g_ah, biogeochem=biogeochem, 
                                    leaf_N_conc=leaf_N_conc, u_leaf=u_leaf, 
                                    d_leaf=d_leaf, met_cond=met_cond, tol=tol,g1_med=g1_med)
    A_nsha = FBB_sha$A_n
    R_dsha = FBB_sha$R_d
    g_ssha = FBB_sha$g_s
    CUO_sha = FBB_sha$CUO
    
    # Sunlit and shaded plant area index:
    # Make sure K_b > 0 always. At night, set K_b = 1e6.
    # We believe that A_can, R_can and g_can should be scaled up by LAI only, not by LAI + SAI. Therefore, L_sun and L_sha here should be sunlit and shaded leaf area index, not plant area index, but we still consider both LAI and SAI when calculating light extinction. Now they are renamed "LAI_sun" and "LAI_sha" and "LAI_sun" has to taken explicitly from canopy radiative transfer model. (Tai, Feb 2019)
    # L_sun = (1 - exp(-K_b*(LAI + SAI)))/K_b
    # L_sha = (LAI + SAI) - L_sun
    # Shaded leaf area index:
    LAI_sha = LAI - LAI_sun
    
    # Canopy photosynthesis (umol CO2 m^-2 s^-1) and conductance (umol m^-2 s^-1 or m s^-1):
    # "A_can" used to mean NET canopy photosynthesis, but to be consistent with literature, we decided to modify its definition to directly mean GROSS canopy photosynthesis, i.e., GPP, so that GPP = A_can instead of GPP = A_can + R_can as previously done. (Tai, Feb 2018)
    # A_can = A_nsun*L_sun + A_nsha*L_sha
    A_can = (A_nsun + R_dsun)*LAI_sun + (A_nsha + R_dsha)*LAI_sha
    R_can = R_dsun*LAI_sun + R_dsha*LAI_sha
    g_can = LAI_sun/(1/g_b + 1/g_ssun) + LAI_sha/(1/g_b + 1/g_ssha)
    
    # The way to calculate CUO below is not right... We will have to correct it . (Tai, Feb 2019)
    CUO_can = CUO_sun*LAI_sun + CUO_sha*LAI_sha
    
    # Stomatal conductance weighted by sunlit and shaded fraction (umol m^-2 s^-1 or m s^-1):
    g_s = (g_ssun*LAI_sun + g_ssha*LAI_sha)/(LAI_sun + LAI_sha)
    
    output = list(A_can=A_can, R_can=R_can, g_can=g_can, CUO_can=CUO_can, A_nsun=A_nsun, A_nsha=A_nsha, R_dsun=R_dsun, R_dsha=R_dsha, g_ssun=g_ssun, g_ssha=g_ssha, g_s=g_s, g_b=g_b, CUO_sun=CUO_sun, CUO_sha=CUO_sha)
    return(output)
    
}

################################################################################

# Soil water stress function (0 to 1):
# This is to be applied to V_cmax, R_d and minimum stomatal conductance (b) when calculating photosynthesis.
f_water_stress = function(soil_wetness=NULL, soil_wetness_top=NULL, soil_wetness_bottom=NULL, theta_w=NULL, theta_sat=NULL, psi_sat, b_psi, psi_sat_top, b_psi_top, psi_sat_bottom, b_psi_bottom, psi_c, psi_o, root_frac_in_top, multilayer) {
  # Soil wetness (i.e. % saturation): soil_wetness
    # Volumetric soil water content (fraction): theta_w
    # Saturated volumetric water content (fraction): theta_sat
    # Either soil_wetness or (theta_w, theta_sat) have to be specified.
    ### Saturated soil matric potential (mm): psi_sat = -478 ### (wrong line; pointed out by Lam, Jan 2019)
    # Saturated soil matric potential (mm): psi_sat
    # Clapp and Homberger parameter: b_psi
    # Soil matric potential when stomata are fully closed (mm): psi_c
    # Soil matric potential when stomata are fully opened (mm): psi_o
    # We assume an ice-free soil. The capacity to consider ice is under development.
    if (multilayer == 'multilayer') {
        # A multilayer soil column model to calculate soil physics is under development and not available at the moment.
        stop('Multilayer model is not available in the current version.')
    } else if (multilayer == 'two-layer') {
        if (is.null(soil_wetness_bottom)) soil_wetness_bottom = theta_w/theta_sat
        if (is.null(soil_wetness_top)) soil_wetness_top = theta_w/theta_sat
        soil_wetness_bottom = max(c(0.01, soil_wetness_bottom), na.rm=TRUE)
        soil_wetness_top = max(c(0.01, soil_wetness_top), na.rm=TRUE)
        # We use a bulk parameterization for soil aggregate root zone.
        psi = max(c(psi_c, psi_sat_bottom*soil_wetness_bottom^-b_psi_bottom), na.rm=TRUE)
        # Soil matric potential of top soil layer
        psi_top = max(c(psi_c, psi_sat_top*soil_wetness_top^-b_psi_top), na.rm=TRUE)
        # Plant wilting factor (bulk):
        wilt_factor = min(c(max(c(0, (psi_c - psi)/(psi_c - psi_o)), na.rm=TRUE), 1), na.rm=TRUE)
        # Plant wilting factor (top soil layer):
        wilt_factor_top = min(c(max(c(0, (psi_c - psi_top)/(psi_c - psi_o)), na.rm=TRUE), 1), na.rm=TRUE)
        # Multiply by a scaling factor for match multilayer model results. It is set to be one for now.
        #gamma = 1
        #beta_t = wilt_factor*gamma
        beta_t = wilt_factor_top * root_frac_in_top + wilt_factor * (1 - root_frac_in_top)
        return(beta_t)
    } else if (multilayer == 'bulk') {
        if (is.null(soil_wetness)) soil_wetness = theta_w/theta_sat
        soil_wetness = max(c(0.01, soil_wetness), na.rm=TRUE)
        # We use a bulk parameterization for soil aggregate root zone.
        psi = max(c(psi_c, psi_sat*soil_wetness^-b_psi), na.rm=TRUE)
        # Plant wilting factor (bulk):
        wilt_factor = min(c(max(c(0, (psi_c - psi)/(psi_c - psi_o)), na.rm=TRUE), 1), na.rm=TRUE)
        # Multiply by a scaling factor for match multilayer model results. It is set to be one for now.
        gamma = 1
        beta_t = wilt_factor*gamma
        return(beta_t)
    }
}
###############################################################################
