################################################################################
### Input parameters for PFT and land surface from NCDF files:
################################################################################

# This module requires functions from "tools.R".
# It requires data directories, "met_name", "lon", "lat", "dlon" and "dlat" defined and R packages loaded in "execution.R".

################################################################################
### Revision history
################################################################################

# Aug 2017: v1.0 finalized (Tai)
# Oct 2017: Bug found in "psi_sat_PFT = sp.regrid(spdata=zero.to.NaN(PSI_SAT_BULK), ...)". Since PSI_SAT_BULK is always negative, it is always turned into NaN and then zero. Now fixed to "psi_sat_PFT = -sp.regrid(spdata=zero.to.NaN(-PSI_SAT_BULK), ...)". (Wong, Tai)
# Feb 2018: "PCT_PFT" was previously incorrectly interpreted as % PFT of grid cell. Now it is corrected - it should % PFT on land, so "PFT_frac" before regridding should be "PCT_PFT/100*LANDFRAC_PFT". (Tai)
# Feb 2018: A new variable for areas of grid cells, "grid_area", was added, and now stored in "surf_data_regrid_2000.RData" too. (Tai)
# Feb 2019: Changed code for reading newer versions (ie vBNU) of GEOS-Chem MODIS LAI (Yung)
# Feb 2019: Completely revamped the input methods for reading in interannually varying LAI data and O3 data.

################################################################################
### PFT-specific parameters from NCDF file:
################################################################################

print('Loading PFT-specific parameters...', quote=FALSE)

filename = paste0(surf_data_dir, 'pftdata/pft-physiology.c130503.nc')

# PFT parameters:
nc = nc_open(filename)
# PFT description:
pftname = ncvar_get(nc, 'pftname')
# PFT number:
pftnum = ncvar_get(nc, 'pftnum')
# Photosynthetic pathway (1 = C3; 0 = C4):
c3psn = ncvar_get(nc, 'c3psn')
# Crop (0 = non-crop; 1 = crop):
crop = ncvar_get(nc, 'crop')
# Ratio of displacement height to canopy top height (0-1):
displar = ncvar_get(nc, 'displar')
# Characteristic leaf dimension (m):
dleaf = ncvar_get(nc, 'dleaf')
# Through canopy (projected area basis) dSLA/dLAI (m^2 gC^-1):
dsladlai = ncvar_get(nc, 'dsladlai')
# Evergreen leaf habit (1 = evergreen; 0 = not):
evergreen = ncvar_get(nc, 'evergreen')
# Fraction of leaf N in Rubisco enzyme (0-1):
flnr = ncvar_get(nc, 'flnr')
# Irrigated (0 = rainfed; 1 = irrigrated):
irrigated = ncvar_get(nc, 'irrigated')
# Leaf longevity (yr):
leaf_long = ncvar_get(nc, 'leaf_long')
# Leaf C:N (gC gN^-1):
leafcn = ncvar_get(nc, 'leafcn'); leafcn[18:19] = 10; leafcn[20:23] = 15
# Leaf reflectance for near-IR (0-1):
rholnir = ncvar_get(nc, 'rholnir')
# Leaf reflectance for visible (0-1):
rholvis = ncvar_get(nc, 'rholvis')
# Stem reflectance for near-IR (0-1):
rhosnir = ncvar_get(nc, 'rhosnir')
# Stem reflectance for visible (0-1):
rhosvis = ncvar_get(nc, 'rhosvis')
# CLM rooting distribution parameter a (m^-1):
roota_par = ncvar_get(nc, 'roota_par')
# CLM rooting distribution parameter b (m^-1):
rootb_par = ncvar_get(nc, 'rootb_par')
# Rooting beta parameter, for C and N vertical discretization:
rootprof_beta = ncvar_get(nc, 'rootprof_beta')
# Seasonal deciduous (0 = not; 1 = seasonal-deciduous):
season_decid = ncvar_get(nc, 'season_decid')
# Specific leaf area (SLA) at top of canopy (projected area basis) (m^2 gC^-1):
slatop = ncvar_get(nc, 'slatop')
# Soil water potential at full stomatal closure (mm):
smpsc = ncvar_get(nc, 'smpsc')
# Soil water potential at full stomatal opening (mm):
smpso = ncvar_get(nc, 'smpso')
# Stress deciduous (0 = not; 1 = stress-deciduous):
stress_decid = ncvar_get(nc, 'stress_decid')
# Leaf transmittance for near-IR (0-1):
taulnir = ncvar_get(nc, 'taulnir')
# Leaf transmittance for visible (0-1):
taulvis = ncvar_get(nc, 'taulvis')
# Stem transmittance for near-IR (0-1):
tausnir = ncvar_get(nc, 'tausnir')
# Stem transmittance for visible (0-1):
tausvis = ncvar_get(nc, 'tausvis')
# Woody life form (0 = non-woody; 1 = woody):
woody = ncvar_get(nc, 'woody')
# Leaf/stem orientation index:
xl = ncvar_get(nc, 'xl')
# Ratio of momentum roughness length to canopy top height (0-1):
z0mr = ncvar_get(nc, 'z0mr')

if (biogeochem_flag){
    # Leaf allocation coefficient parameter for crops
    aleaff = ncvar_get(nc, 'aleaff')
    allconsl = ncvar_get(nc, 'allconsl')
    bfact = ncvar_get(nc, 'bfact')
    fleafi = ncvar_get(nc, 'fleafi')
    # Stem allocation coefficient parameter for crops
    astemf = ncvar_get(nc, 'astemf'); astemf[18:19] = 0
    allconss = ncvar_get(nc, 'allconss')
    declfact = ncvar_get(nc, 'declfact')
    # Root allocation coefficient parameter for crops
    arooti = ncvar_get(nc, 'arooti')
    arootf = ncvar_get(nc, 'arootf'); arootf[20:23] = 0      # The value should be 0 instead of NA
    # Base temperature og GDD accumulation for crops
    baset = ncvar_get(nc, 'baset')
    # Allocation ratio of coarse root : live stem 
    croot_stem = ncvar_get(nc, 'croot_stem')
    # C:N ratio of dead wood
    deadwdcn = ncvar_get(nc, 'deadwdcn')
    # Allocation flag to storage pools (deciduous) or display pools (evergreen and crops)
    fcur = ncvar_get(nc, 'fcur')
    # C:N ratio of fine root during reproductive stage for crops
    ffrootcn = ncvar_get(nc, 'ffrootcn'); ffrootcn[c(18:19,24:25)] = 42
    # C:N ratio of leaf during reproductive stage for crops
    fleafcn = ncvar_get(nc, 'fleafcn')
    # Allocation flag to wood
    flivewd = ncvar_get(nc, 'flivewd')
    # C:N ratio of fine root
    frootcn = ncvar_get(nc, 'frootcn'); frootcn[20:23] = 30
    # Allocation ratio of fine root : leaf
    froot_leaf = ncvar_get(nc, 'froot_leaf')
    # C:N ratio of live stem during reproductive stage
    fstemcn = ncvar_get(nc, 'fstemcn')
    # C:N ratio of grain
    graincn = ncvar_get(nc, 'graincn'); graincn[1:17] = 999; graincn[20:23] = 40; graincn[24:25] = 60  # Prevent dividing by NA for natural vegetation
    # % of GDD_mat to reach reproductive stage
    grnfill = ncvar_get(nc, 'grnfill')
    # Maximum GDD_mat allowed for crops
    hybgdd = ncvar_get(nc, 'hybgdd')
    # Prescribed maximum LAI allowed for crops (can be removed)
    laimx = ncvar_get(nc, 'laimx')
    if(!limit_crop_LAI_flag){laimx[18:25] = 999}
    # Leaf longevity (same as above)
    leaf_long[18:19] = 1/6                     #maize leaf longetivity can be reduced to ~60 days (Wolfe et al. (1988), Agronomy Journal)
    leaf_long[24:25] = 1/15                    #soybean leaf longetivity maybe be reduced to 20-30 days (Miyaji (1984), New Phytologist)
    # % of GDD_mat to reach vegetative stage
    lfemerg = ncvar_get(nc,'lfemerg')
    # C:N ratio of leaf litter
    lflitcn = ncvar_get(nc, 'lflitcn')
    # C:N ratio live wood
    livewdcn = ncvar_get(nc, 'livewdcn')
    # Daily min planting temperature requirment (K)
    min_planting_temp = ncvar_get(nc, 'min_planting_temp')
    # Maximum growing season length allowed
    mxmat = ncvar_get(nc, 'mxmat')
    # Maximum increase in GDD_T2m allowed
    mxtmp = ncvar_get(nc, 'mxtmp')
    # Daily average planting temperature requirement (K)
    planting_temp = ncvar_get(nc, 'planting_temp'); planting_temp[22:23] = 0
    # Allocation ratio of live stem : leaf
    stem_leaf = ncvar_get(nc, 'stem_leaf')
    # Maximum canopy height for crops (m)
    ztopmx = ncvar_get(nc, 'ztopmx')
    # Earliest/lastest possible planting day in the Northern Hemisphere (julian day)
    earliest_planting_jday_possible_NH = c(rep(NA,17), rep(91,4), rep(244,2), rep(121,2))
    latest_planting_jday_possible_NH = c(rep(NA,17), rep(166,4), rep(334,2), rep(166,2))
    # Earliest/lastest possible planting day in the Southern Hemisphere (julian day)
    earliest_planting_jday_possible_SH = c(rep(NA,17), rep(274,4), rep(121,2), rep(305,2))
    latest_planting_jday_possible_SH = c(rep(NA,17), rep(349,4), rep(151,2), rep(349,2))
}

nc_close(nc)

# Canopy top height (m):
htop = c(0, 17, 17, 14, 35, 35, 18, 20, 20, rep(0.5, times=16))
# Canopy bottom height (m):
hbot = c(0, 8.5, 8.5, 7, 1, 1, 10, 11.5, 11.5, 0.1, 0.1, 0.1, rep(0.01, times=13))
# Plant group (for use in ozone-vegetation scheme):
plantgroup = c('none', 'needleleaf', 'needleleaf', 'needleleaf', 'broadleaf', 'broadleaf', 'broadleaf', 'broadleaf', 'broadleaf', 'shrub', 'shrub', 'shrub', 'C3_grass', 'C3_grass', 'C4_grass', 'C3_grass', 'C3_grass', 'C4_grass', 'C4_grass', 'C3_grass', 'C3_grass', 'C3_grass', 'C3_grass', 'C3_grass', 'C3_grass')

# Mass ratio of total Rubisco molecular mass to nitrogen in Rubisco (g Rubisco gN^-1):
fnr = 7.16
# Specific activity of Rubisco (umol CO2 g^-1 Rubisco s^-1):
ar25 = 60
# Leaf nitrogen concentration at top of canopy (g N m^-2 leaf area):
lnctop = 1/(leafcn*slatop)
# Maximum rate of carboxylation at top of canopy at 25 degC (umol CO2 m^-2 s^-1):
vcmax25top = lnctop*flnr*fnr*ar25
vcmax25top[1] = 0

# Parameters for Medlyn model
g1_med_table = c(NA, 2.35, 2.35, 2.35, 4.12, 4.12, 4.45, 4.45, 4.45, 4.7, 4.7, 4.7, 2.22, 5.25, 1.62, NA, NA, 1.79, 1.79, NA, NA, NA, NA, 5.79, 5.79)

################################################################################
### Surface data from CLM default ncdf file:
################################################################################

print('Loading surface data...', quote=FALSE)

filename = paste0(surf_data_dir, 'surfdata_map/surfdata_1.9x2.5_mp24_simyr2000_c130419_CLM5crop_coverage.nc')

nc = nc_open(filename)
# Dimensions:
LONGXY = ncvar_get(nc, 'LONGXY')
LON_CLM = LONGXY[,1]
LATIXY = ncvar_get(nc, 'LATIXY')
LAT_CLM = LATIXY[1,]
nlevsoi = 10
# # Indices of natural PFT:
# natpft = ncvar_get(nc, 'natpft')
# # Indices of crop PFT:
# cft = ncvar_get(nc, 'cft')
# Grid cell area (km^2):
AREA = ncvar_get(nc, 'AREA')
# Monthly leaf area index (m^2 m^-2):
MONTHLY_LAI = ncvar_get(nc, 'MONTHLY_LAI')
# Monthly stem area index (m^2 m^-2):
MONTHLY_SAI = ncvar_get(nc, 'MONTHLY_SAI')
# Monthly canopy top height (m) (same every month):
MONTHLY_HEIGHT_TOP = ncvar_get(nc, 'MONTHLY_HEIGHT_TOP')
# Monthly canopy bottom height (m) (same every month):
MONTHLY_HEIGHT_BOT = ncvar_get(nc, 'MONTHLY_HEIGHT_BOT')
# Organic matter density at soil levels (kg m^-3):
ORGANIC = ncvar_get(nc, 'ORGANIC')
# Percent clay (%):
PCT_CLAY = ncvar_get(nc, 'PCT_CLAY')
# Percent sand (%):
PCT_SAND = ncvar_get(nc, 'PCT_SAND')
# Maximum numbers of soil colors:
mxsoil_color = ncvar_get(nc, 'mxsoil_color')
# Soil color:
SOIL_COLOR = ncvar_get(nc, 'SOIL_COLOR')

# Percent plant functional type on land (% of land area in grid cell):
# It was previously incorrectly interpreted as % PFT of grid cell. Now it is corrected. (Tai, Feb 2018)
PCT_PFT = ncvar_get(nc, 'PCT_PFT')
# Land fraction from PFT dataset (0-1):
LANDFRAC_PFT = ncvar_get(nc, 'LANDFRAC_PFT')
# Land mask from PFT dataset, indicative of real/fake points:
PFTDATA_MASK = ncvar_get(nc, 'PFTDATA_MASK')
nc_close(nc)

# The following are variables from file "surfdata_1.9x2.5_mp24_simyr2000_c141219.nc".
# # Precent crop functional type on crop land unit (% of land unit):
# PCT_CFT = ncvar_get(nc, 'PCT_CFT')
# # Total percent crop land unit (%):
# PCT_CROP = ncvar_get(nc, 'PCT_CROP')
# # Percent plant functional type on natural vegetation land unit (% of land unit):
# PCT_NAT_PFT = ncvar_get(nc, 'PCT_NAT_PFT')
# # Total percent natural vegetation land unit (%):
# PCT_NATVEG = ncvar_get(nc, 'PCT_NATVEG')

# Flip longitude if it is in the range of 0 to 360:
if (length(which(LON_CLM > 180)) > 0) {
   print('Flipping longitudes of input data...', quote=FALSE)
   AREA = flip.lon(spdata=AREA, lon=LON_CLM)$spdata
   MONTHLY_LAI = flip.lon(spdata=MONTHLY_LAI, lon=LON_CLM)$spdata
   MONTHLY_SAI = flip.lon(spdata=MONTHLY_SAI, lon=LON_CLM)$spdata
   MONTHLY_HEIGHT_TOP = flip.lon(spdata=MONTHLY_HEIGHT_TOP, lon=LON_CLM)$spdata
   MONTHLY_HEIGHT_BOT = flip.lon(spdata=MONTHLY_HEIGHT_BOT, lon=LON_CLM)$spdata
   ORGANIC = flip.lon(spdata=ORGANIC, lon=LON_CLM)$spdata
   PCT_CLAY = flip.lon(spdata=PCT_CLAY, lon=LON_CLM)$spdata
   PCT_SAND = flip.lon(spdata=PCT_SAND, lon=LON_CLM)$spdata
   SOIL_COLOR = flip.lon(spdata=SOIL_COLOR, lon=LON_CLM)$spdata
   PCT_PFT = flip.lon(spdata=PCT_PFT, lon=LON_CLM)$spdata
   LANDFRAC_PFT = flip.lon(spdata=LANDFRAC_PFT, lon=LON_CLM)$spdata
   PFTDATA_MASK = flip.lon(spdata=PFTDATA_MASK, lon=LON_CLM)$spdata
   LON_CLM = flip.lon(spdata=AREA, lon=LON_CLM)$lon
} else {
   # Do nothing.
}

################################################################################

# Derive soil hydrological parameters for calculating soil water stress and soil albedo:

# Soil organic matter fraction (0-1):
om_frac = ORGANIC/max(ORGANIC, na.rm=TRUE)
# Porosity of mineral soil:
theta_satmin = 0.489 - 0.00126*PCT_SAND
# Porosity of organic matter (0-1):
theta_satom = 0.9
# Saturated volumetric water content (0-1):
theta_sat_xyz = (1 - om_frac)*theta_satmin + om_frac*theta_satom
# Clapp and Homberger parameter for mineral soil:
b_psimin = 2.91 + 0.159*PCT_CLAY
# Clapp and Homberger parameter for organic matter:
b_psiom = 2.7
# Clapp and Homberger parameter:
b_psi_xyz = (1 - om_frac)*b_psimin + om_frac*b_psiom
# Saturated mineral soil matric potential (mm):
psi_satmin = -10.0*10^(1.88 - 0.0131*PCT_SAND)
# Saturated organic matter matric potential (mm):
psi_satom = -10.3
# Saturated soil matric potential (mm):
psi_sat_xyz = (1 - om_frac)*psi_satmin + om_frac*psi_satom

# Depth at layer interface for top 10 layers plus the soil surface (m):
z_hi = c(0, 0.0175, 0.0451, 0.0906, 0.1655, 0.2891, 0.4929, 0.8289, 1.3828, 2.2961, 3.8019)
# Root fraction (0-1) (dim1 = PFT; dim2 = soil level):
root_frac = matrix(0, nrow=length(pftname), ncol=nlevsoi)
for (i in 1:nlevsoi) {
   if (i < nlevsoi) {
      root_frac[,i] = 0.5*(exp(-roota_par*z_hi[i]) + exp(-rootb_par*z_hi[i]) - exp(-roota_par*z_hi[i+1]) - exp(-rootb_par*z_hi[i+1]))
   } else {
      root_frac[,i] = 0.5*(exp(-roota_par*z_hi[i]) + exp(-rootb_par*z_hi[i]))
   }
}

# Calculate bulk soil parameters by PFT for single soil column:
subfn = 'bulk_soil_PFT_2000.RData'
filename = paste0(processed_surf_data_dir, subfn)

if (length(dir(path=processed_surf_data_dir, pattern=subfn)) == 0) {
   # Root zone bulk saturated volumetric water content (0-1), Clapp and Homberger parameter, saturated soil matric potential (mm), percent of sand, clay and organic matter fraction:
   THETA_SAT_BULK = array(NaN, dim=c(length(LON_CLM), length(LAT_CLM), length(pftname)))
   B_PSI_BULK = array(NaN, dim=c(length(LON_CLM), length(LAT_CLM), length(pftname)))
   PSI_SAT_BULK = array(NaN, dim=c(length(LON_CLM), length(LAT_CLM), length(pftname)))
   PCT_SAND_BULK = array(NaN, dim=c(length(LON_CLM), length(LAT_CLM), length(pftname)))
   PCT_CLAY_BULK = array(NaN, dim=c(length(LON_CLM), length(LAT_CLM), length(pftname)))
   OM_FRAC_BULK = array(NaN, dim=c(length(LON_CLM), length(LAT_CLM), length(pftname)))
   for (i in 1:length(LON_CLM)) {
      if (i/10 == floor(i/10)) print(paste0('Calculating bulk soil parameters for lon = ', as.character(LON_CLM[i]), '...'), quote=FALSE)
      for (j in 1:length(LAT_CLM)) {
         for (ipft in 1:length(pftname)) {
            THETA_SAT_BULK[i,j,ipft] = sum(theta_sat_xyz[i,j,]*root_frac[ipft,], na.rm=TRUE)
            B_PSI_BULK[i,j,ipft] = sum(b_psi_xyz[i,j,]*root_frac[ipft,], na.rm=TRUE)
            PSI_SAT_BULK[i,j,ipft] = sum(psi_sat_xyz[i,j,]*root_frac[ipft,], na.rm=TRUE)
            PCT_SAND_BULK[i,j,ipft] = sum(PCT_SAND[i,j,]*root_frac[ipft,], na.rm=TRUE)
            PCT_CLAY_BULK[i,j,ipft] = sum(PCT_CLAY[i,j,]*root_frac[ipft,], na.rm=TRUE)
            OM_FRAC_BULK[i,j,ipft] = sum(om_frac[i,j,]*root_frac[ipft,], na.rm=TRUE)
         }
      }
   }
   exist_bulk_soil = TRUE
   save(list=c('THETA_SAT_BULK', 'B_PSI_BULK', 'PSI_SAT_BULK', 'PCT_SAND_BULK', 'PCT_CLAY_BULK', 'OM_FRAC_BULK', 'LON_CLM', 'LAT_CLM', 'pftname', 'pftnum', 'exist_bulk_soil'), file=filename)
} else {
   if (length(ls(pattern='exist_bulk_soil')) == 0) {
      print('Loading bulk soil parameters...', quote=FALSE)
      load(filename)
   }
}

# Saturated volumetric water content (0-1) of top soil (z = 0):
THETA_SAT_TOP = theta_sat_xyz[,,1]

# Soil albedo parameters (which varies with soil color class):
soil_col_albedo = `row.names<-.data.frame`(`colnames<-`(rbind.data.frame(
   c(0.60, 0.40, 0.60, 0.40),
   c(0.36, 0.61, 0.25, 0.50),
   c(0.34, 0.57, 0.23, 0.46),
   c(0.32, 0.53, 0.21, 0.42),
   c(0.31, 0.51, 0.20, 0.40),
   c(0.30, 0.49, 0.19, 0.38),
   c(0.29, 0.48, 0.18, 0.36),
   c(0.28, 0.45, 0.17, 0.34),
   c(0.27, 0.43, 0.16, 0.32),
   c(0.26, 0.41, 0.15, 0.30),
   c(0.25, 0.39, 0.14, 0.28),
   c(0.24, 0.37, 0.13, 0.26),
   c(0.23, 0.35, 0.12, 0.24),
   c(0.22, 0.33, 0.11, 0.22),
   c(0.20, 0.31, 0.10, 0.20),
   c(0.18, 0.29, 0.09, 0.18),
   c(0.16, 0.27, 0.08, 0.16),
   c(0.14, 0.25, 0.07, 0.14),
   c(0.12, 0.23, 0.06, 0.12),
   c(0.10, 0.21, 0.05, 0.10),
   c(0.08, 0.16, 0.04, 0.08),
   stringsAsFactors = FALSE), 
   c('dry/vis', 'dry/nir', 'saturated/vis', 'saturated/nir')), paste0('soil_col_', 0:20))

# Soil albedo (which varies with soil color class):
# 3rd dim = [dry/vis, dry/nir, saturated/vis, saturated/nir]
SOIL_ALBEDO = array(0, dim=c(length(LON_CLM), length(LAT_CLM), 4))
for (i in 1:length(LON_CLM)) {
   for (j in 1:length(LAT_CLM)) {
      # Grid cells that are not land:
      if (PFTDATA_MASK[i,j] == 0) {
         alpha_soil = rep(0, times=4)
      } else {
      soil_col = SOIL_COLOR[i,j]
      alpha_soil = soil_col_albedo[paste0('soil_col_', as.character(soil_col)), ]
      }
      SOIL_ALBEDO[i,j,] = as.numeric(alpha_soil)
   }
}

# # Compare with another method of calculation:
# var_to_plot = (1 - om_frac_xypft)*(2.91 + 0.159*PCT_CLAY_xypft) + om_frac_xypft*b_psiom
# plot.field(var_to_plot[,,5], lon_CLM, lat_CLM, Pacific.centric=TRUE)
# var_to_plot = (1 - om_frac_xypft)*(-10.0*10^(1.88 - 0.0131*PCT_SAND_xypft)) + om_frac_xypft*psi_satom
# plot.field(var_to_plot[,,5], lon_CLM, lat_CLM, Pacific.centric=TRUE)
# # They are very similar anyway...

################################################################################
### Surface invariant fields from GEOS-FP or MERRA2 default ncdf file:
################################################################################

# GEOS-FP or MERRA2 should have the same resolution as current model resolution.
if (met_name == 'GEOSFP') {
   filename = paste0(met_data_dir, 'GEOS_FP/2011/01/', met_name, '.20110101.CN.2x25.nc')
} else if (met_name == 'MERRA2') {
   filename = paste0(met_data_dir, 'MERRA2/2015/01/', met_name, '.20150101.CN.2x25.nc4')
} else {
   stop('met_name specified is not available!')
}

print(paste0('Loading invariant ', met_name, ' surface fields from ', met_data_dir, '...'), quote=FALSE)

nc = nc_open(filename)
# Dimensions:
lon_met = ncvar_get(nc, 'lon')
lat_met = ncvar_get(nc, 'lat')
# Fraction of lake in grid cell (0-1):
FRLAKE = ncvar_get(nc, 'FRLAKE')
# Fraction of land in grid cell (0-1):
FRLAND = ncvar_get(nc, 'FRLAND')
# Fraction of land ice in grid cell (0-1):
FRLANDIC = ncvar_get(nc, 'FRLANDIC')
# Fraction of ocean in grid cell (0-1):
FROCEAN = ncvar_get(nc, 'FROCEAN')
# Surface geopotential (m^2 s^-2):
PHIS = ncvar_get(nc, 'PHIS')
nc_close(nc)

# Check if resolution of meteorological fields (GEOS-FP or MERRA2) matches that of TEMIR model:
if (!identical(as.numeric(lon_met),as.numeric(lon)) | !identical(head(as.numeric(lat_met),-1)[-1], head(as.numeric(lat),-1)[-1])) stop('Model resolution does not match that of input meteorological data!')

# Compute areas of grid cells [km^2]:
# Newly added. (Tai, Feb 2018)
grid_area = matrix(NaN, nrow=length(lon), ncol=length(lat))
for (i in 1:length(lon)) {
   for (j in 1:length(lat)) {
      grid_area[i,j] = area.latlon(lat[j]-dlat/2, lon[i]-dlon/2, lat[j]+dlat/2, lon[i]+dlon/2)
   }
}

print('Model resolution:', quote=FALSE)
print(paste0('Longitudes from ', as.character(head(lon, 1)), ' to ', as.character(tail(lon, 1)), ' by ', as.character(dlon), ' deg'), quote=FALSE)
print(paste0('Latitudes from ', as.character(head(lat, 1)), ' to ', as.character(tail(lat, 1)), ' by ', as.character(dlat), ' deg'), quote=FALSE)

################################################################################
### Regrid surface data:
################################################################################

# Regrid surface data to match with met field resolution:

subfn = 'surf_data_regrid_2000.RData'
filename = paste0(processed_surf_data_dir, subfn)

if (!file.exists(filename)) {
   
   # Monthly LAI and SAI by PFT (m^2 m^-2):
   LAI_mon_PFT = array(NaN, dim=c(length(lon), length(lat), length(pftname), 12))
   SAI_mon_PFT = array(NaN, dim=c(length(lon), length(lat), length(pftname), 12))
   for (m in 1:12) {
      print(paste0('Regridding LAI and SAI for month = ', as.character(m), '...'), quote=FALSE)
      LAI_mon_PFT[,,,m] = sp.regrid(spdata=zero.to.NaN(MONTHLY_LAI[,,,m]), lon.in=LON_CLM, lat.in=LAT_CLM, lon.out=lon, lat.out=lat)
      SAI_mon_PFT[,,,m] = sp.regrid(spdata=zero.to.NaN(MONTHLY_SAI[,,,m]), lon.in=LON_CLM, lat.in=LAT_CLM, lon.out=lon, lat.out=lat)
   }
   LAI_mon_PFT = NaN.to.zero(LAI_mon_PFT)
   SAI_mon_PFT = NaN.to.zero(SAI_mon_PFT)
   
   print('Regridding PFT fraction...', quote=FALSE)
   # PFT fraction by PFT (0-1):
   # "PCT_PFT" was previously incorrectly interpreted as % PFT of grid cell. Now it is corrected - it should % PFT on land, so "PFT_frac" before regridding should be "PCT_PFT/100*LANDFRAC_PFT". (Tai, Feb 2018)
   # PFT_frac = sp.regrid(spdata=PCT_PFT/100, lon.in=LON_CLM, lat.in=LAT_CLM, lon.out=lon, lat.out=lat)
   PCT_PFT_GRID = PCT_PFT*array(LANDFRAC_PFT, dim=c(dim(LANDFRAC_PFT), dim(PCT_PFT)[3]))
   PFT_frac = sp.regrid(spdata=PCT_PFT_GRID/100, lon.in=LON_CLM, lat.in=LAT_CLM, lon.out=lon, lat.out=lat)
   
   print('Regridding theta_sat, theta_sat0, b_psi, psi_sat, soil albedo...', quote=FALSE)
   # Saturated volumetric water content by PFT (0-1):
   theta_sat_PFT = sp.regrid(spdata=zero.to.NaN(THETA_SAT_BULK), lon.in=LON_CLM, lat.in=LAT_CLM, lon.out=lon, lat.out=lat)
   theta_sat_PFT = NaN.to.zero(theta_sat_PFT)
   # Saturated volumetric water content of top soil (0-1):
   theta_sat0 = sp.regrid(spdata=zero.to.NaN(THETA_SAT_TOP), lon.in=LON_CLM, lat.in=LAT_CLM, lon.out=lon, lat.out=lat)
   theta_sat0 = NaN.to.zero(theta_sat0)
   # Clapp and Homberger parameter by PFT:
   b_psi_PFT = sp.regrid(spdata=zero.to.NaN(B_PSI_BULK), lon.in=LON_CLM, lat.in=LAT_CLM, lon.out=lon, lat.out=lat)
   b_psi_PFT = NaN.to.zero(b_psi_PFT)
   # Saturated soil matric potential by PFT (mm):
   # It is always negative, so can't use the "zero.to.Nan" function in the same way... (fixed by Tai on 25 Oct 2017)
   psi_sat_PFT = -sp.regrid(spdata=zero.to.NaN(-PSI_SAT_BULK), lon.in=LON_CLM, lat.in=LAT_CLM, lon.out=lon, lat.out=lat)
   psi_sat_PFT = NaN.to.zero(psi_sat_PFT)
   # Soil albedo:
   # 3rd dim = [dry/vis, dry/nir, saturated/vis, saturated/nir]
   soil_albedo = sp.regrid(spdata=zero.to.NaN(SOIL_ALBEDO), lon.in=LON_CLM, lat.in=LAT_CLM, lon.out=lon, lat.out=lat)
   soil_albedo = NaN.to.zero(soil_albedo)
   
   print('Regridding soil color...', quote=FALSE)
   # Soil color (1-20 classes):
   soil_color = sp.regrid(spdata=SOIL_COLOR, lon.in=LON_CLM, lat.in=LAT_CLM, lon.out=lon, lat.out=lat, method='mode')
   
   exist_regrid = TRUE
   
   save(list=c('LAI_mon_PFT', 'SAI_mon_PFT', 'PFT_frac', 'theta_sat_PFT', 'b_psi_PFT', 'psi_sat_PFT', 'soil_color', 'soil_albedo', 'grid_area', 'lon', 'lat', 'pftname', 'pftnum', 'exist_regrid'), file=filename)
   
} else {
   if (length(ls(pattern='exist_regrid')) == 0) {
      print('Loading regridded surface data...', quote=FALSE)
      load(filename)
   }
}

cat('\n')

################################################################################
### Regrid planting and harvest date data from Sack et al. (2010): (Pang and Sadiq, Jun 2019)
################################################################################

# Regrid the planting and harvesting date data if the simulation is BGC (with crop)
# The resolution of the original data is 0.5 x 0.5

if ((biogeochem_flag && get_planting_date_option == 'prescribed-map') || O3_POD) {
    
    if (dlon < 0.5 || dlat < 0.5){
        # Is this necessary???? (Pang, Jun 2019)
        warning('The default planting and harvesting data has a resolution of 0.5x0.5. Simulation resoltuion is too high that there may be some problems in regridding.')
    }
    
    if (dlon == 0.5 && dlat == 0.5){
        # Simulation resolution is the same as the data, can be read in directly from .nc
        filename = 'Sack_crop_calendar.nc'
        nc = nc_open(paste0(planting_date_map_dir, filename))
        prescribed_planting_date_Sack = ncvar_get(nc, 'planting')
        prescribed_harvesting_date_Sack = ncvar_get(nc, 'harvest')
        nc_close(nc)
    } else {
        # Like the prescribed LAI and soil data, read in the RData contains the regridded data if it exists. Otherwise, regrid the nc data to a suitable resolution
        filename = paste0(planting_date_map_dir,'crop_planting_harvesting_Sack_',dlat,'x',dlon,'.RData')
        
        if (!file.exists(filename)){
            # Read in the nc file and regrid the data to the simulation resolution
            nc_filename = 'Sack_crop_calendar.nc'
            nc = nc_open(paste0(planting_date_map_dir, nc_filename))
            tmp_planting_date_Sack_nc = ncvar_get(nc, 'planting')
            tmp_harvesting_date_Sack_nc = ncvar_get(nc, 'harvest')
            
            lat_Sack = ncvar_get(nc, 'lat')
            lon_Sack = ncvar_get(nc, 'lon')
            nc_close(nc)
            
            lat_diff_sack = lat_Sack[1] - lat_Sack[2]
            lat_diff = lat[1] - lat[2]
            
            lon_diff_sack = lon_Sack[1] - lon_Sack[2]
            lon_diff = lon[1] - lon[2]
            
            if (lat_diff_sack * lat_diff < 0) {
                print('Fliping the order of lat_Sack for sp.regrid()')
                lat_Sack = rev(lat_Sack)
                tmp_planting_date_Sack_nc[,(1:length(lat_Sack)),] = tmp_planting_date_Sack_nc[,rev((1:length(lat_Sack))),]
                tmp_harvesting_date_Sack_nc[,(1:length(lat_Sack)),] = tmp_harvesting_date_Sack_nc[,rev((1:length(lat_Sack))),]
                
            }
            
            if (lon_diff_sack * lon_diff < 0) {
                print('Fliping the order of lon_Sack for sp.regrid()')
                lon_Sack = rev(lon_Sack)
                tmp_planting_date_Sack_nc[(1:length(lon_Sack)),,] = tmp_planting_date_Sack_nc[rev((1:length(lon_Sack))),,]
                tmp_harvesting_date_Sack_nc[(1:length(lon_Sack)),,] = tmp_harvesting_date_Sack_nc[rev((1:length(lon_Sack))),,]
            }
            
            crop_name_vec = c('maize (primary growing season)', 'wheat', 'winter wheat', 'soybean', 'rice (primary growing season)', 'rice (secondary growing season)', 'maize (secondary growing season)')
            prescribed_planting_date_Sack = array(data = NA, dim = c(length(lon), length(lat), length(crop_name_vec)))
            prescribed_harvesting_date_Sack = array(data = NA, dim = c(length(lon), length(lat), length(crop_name_vec)))
            for (z in seq(crop_name_vec)){
                tmp_plant = sp.regrid(spdata = tmp_planting_date_Sack_nc[,,z], lon.in = lon_Sack, lat.in = lat_Sack, lon.out = lon, lat.out = lat, method = 'mode')
                tmp_har = sp.regrid(spdata = tmp_harvesting_date_Sack_nc[,,z], lon.in = lon_Sack, lat.in = lat_Sack, lon.out = lon, lat.out = lat, method = 'mode')
                prescribed_planting_date_Sack[,,z] = tmp_plant
                prescribed_harvesting_date_Sack[,,z] = tmp_har
                print(paste0('Finished regridding planting and harvesting date for ', crop_name_vec[z]))
            }
            
            save(list = c('lon', 'lat', 'prescribed_planting_date_Sack', 'prescribed_harvesting_date_Sack', 'crop_name_vec'), file = filename)
            rm(tmp_plant , tmp_har, crop_name_vec, lon_Sack, lat_Sack, lon_diff, lat_diff, z)
        } else {
            print(paste0('Crop calendar with a suitable resolution is already existed'), quote = FALSE)
            print(paste0('Loading ', filename, ' ...'), quote = FALSE)
            load(filename)
        }
    }
    
}

################################################################################
### Regrid GDDmat derived from Sack et al. (2010) crop calendar 
################################################################################

if (biogeochem_flag) {
    if (get_GDDmat_method == 'Sack') {
        if (dlon == 360/540 && dlat == 0.5) {
            # The simulation resolution is the same as the input
            filename = 'Sack_crop_calendar_GDDmat_0.667x0.5.nc'
            nc = nc_open(paste0(Sack_GDDmat_dir, filename))
            GDDmat_maize_map = ncvar_get(nc, 'Sack_crop_calendar_GDDmat')[,,1]
            GDDmat_soybean_map = ncvar_get(nc, 'Sack_crop_calendar_GDDmat')[,,2]
            GDDmat_springwheat_map = ncvar_get(nc, 'Sack_crop_calendar_GDDmat')[,,3]
            GDDmat_winterwheat_map = ncvar_get(nc, 'Sack_crop_calendar_GDDmat')[,,4]
            nc_close(nc)
        } else if (dlon < 360/540 || dlat < 0.5) {
            warning('The GDDmat map has a resolution of 0.667x0.5. Simulation resoltuion is too high that there may be some problems in regridding.')
        } else {
            filename = paste0('regridded_crop_GDDmat_Sack_',dlon,'x',dlat,'.RData')
            if (!file.exists(paste0(Sack_GDDmat_dir,filename))) {
                print('Regridding GDDmat data derived from Sack et al. (2010)')
                # Regrid the nc input
                nc_filename = 'Sack_crop_calendar_GDDmat_0.667x0.5.nc'
                nc = nc_open(paste0(Sack_GDDmat_dir, nc_filename))
                tmp_GDDmat_mean = ncvar_get(nc, 'GDDmat_mean') # dim() = 540lon x 361lat x 4crops
                tmp_GDDmat_sd = ncvar_get(nc, 'GDDmat_sd') # this variable "may be" useful for adjusting GDDmat in a particular region (e.g., reduce the matuirty requirement by say: GDDmat = GDDmat_mean - GDDmat_sd)
                lat_Sack = ncvar_get(nc, 'lat')
                lon_Sack = ncvar_get(nc, 'lon')
                nc_close(nc)
                
                lat_diff_sack = lat_Sack[1] - lat_Sack[2]
                lat_diff = lat[1] - lat[2]
                
                lon_diff_sack = lon_Sack[1] - lon_Sack[2]
                lon_diff = lon[1] - lon[2]
                
                if (lat_diff_sack * lat_diff < 0) {
                    print('Fliping the order of lat_Sack for sp.regrid()')
                    lat_Sack = rev(lat_Sack)
                    tmp_GDDmat_mean[,(1:length(lat_Sack)),] = tmp_GDDmat_mean[,rev((1:length(lat_Sack))),]
                    tmp_GDDmat_sd[,(1:length(lat_Sack)),] = tmp_GDDmat_sd[,rev((1:length(lat_Sack))),]
                }
                
                if (lon_diff_sack * lon_diff < 0) {
                    print('Fliping the order of lon_Sack for sp.regrid()')
                    lon_Sack = rev(lon_Sack)
                    tmp_GDDmat_mean[(1:length(lon_Sack)),,] = tmp_GDDmat_mean[rev((1:length(lon_Sack))),,]
                    tmp_GDDmat_sd[(1:length(lon_Sack)),,] = tmp_GDDmat_sd[rev((1:length(lon_Sack))),,]
                }
                
                crop_name_vec = c('maize', 'soybean', 'springwheat', 'winterwheat')
                Sack_derived_GDDmat_mean = array(data = NA, dim = c(length(lon), length(lat), length(crop_name_vec)))
                Sack_derived_GDDmat_sd = array(data = NA, dim = c(length(lon), length(lat), length(crop_name_vec)))
                
                for (z in seq(crop_name_vec)) {
                    tmp_mean = sp.regrid(spdata = tmp_GDDmat_mean[,,z], lon.in = lon_Sack, lat.in = lat_Sack, lon.out = lon, lat.out = lat, method = 'mode')
                    tmp_sd = sp.regrid(spdata = tmp_GDDmat_sd[,,z], lon.in = lon_Sack, lat.in = lat_Sack, lon.out = lon, lat.out = lat, method = 'mode')
                    Sack_derived_GDDmat_mean[,,z] = tmp_mean
                    Sack_derived_GDDmat_sd[,,z] = tmp_sd
                    assign(x = paste0('GDDmat_',crop_name_vec[z],'_map'), value = tmp_mean)
                    print(paste0('Finished regridding derived GDDmat for ', crop_name_vec[z]))
                }
                
                save(list = c('lon', 'lat', 'Sack_derived_GDDmat_mean', 'Sack_derived_GDDmat_sd', 'crop_name_vec'), file = paste0(Sack_GDDmat_dir,filename))
                
                rm(z, tmp_mean, tmp_sd, lat_Sack, lon_Sack, lat_diff_sack, lon_diff_sack, tmp_GDDmat_mean, tmp_GDDmat_sd, nc_filename)
            } else {
                print(paste0('GDDmat map with suitable resolution is already existed'), quote = FALSE)
                print(paste0('Loading ', filename, ' ...'), quote = FALSE)
                load(paste0(Sack_GDDmat_dir,filename))
                for (z in seq(crop_name_vec)) {
                    assign(x = paste0('GDDmat_',crop_name_vec[z],'_map'), value = Sack_derived_GDDmat_mean[,,z])
                }
            }
        }

        
    }
}


################################################################################
### Rescale and interpolate PFT-level LAI and SAI data: (Yung & Tai, Feb 2019)
################################################################################

# Rescale default CLM PFT-level monthly LAI (but not SAI) according to user-defined total LAI (that can be interannually varying):

if (LAI_data_flag) {
   
   # Vector of simulation years:
   year_vec = as.numeric(unique(substr(make.date.vec(start.date=start_date, end.date=end_date), 1, 4)))
   
   # Find LAI data available:
   # Please make sure among the available years for LAI data, there are no missing years in between.
   LAI_nc_files = list.files(LAI_data_dir, pattern = paste0('^', LAI_subn1, '.*', LAI_subn2, '$'))
   LAI_avail_years = sort(as.numeric(str_remove(str_remove(LAI_nc_files, pattern = paste0('^', LAI_subn1)), pattern = paste0(LAI_subn2, '$'))))
   print(paste0('Available years for user-defined LAI data: ', paste(LAI_avail_years, collapse = ', ')), quote = FALSE)
   
   if (length(LAI_avail_years) == 0) exist_LAI_data = FALSE else {
      
      # Match simulation years to LAI data available:
      LAI_used_years = year_vec
      # Checking if LAI data are available for the simulation year:
      if (min(year_vec) < min(LAI_avail_years) | max(year_vec) > max(LAI_avail_years)) {
         ind_start_yrs = which(year_vec < min(LAI_avail_years))
         LAI_used_years[ind_start_yrs] = LAI_avail_years[1]
         ind_end_yrs = which(year_vec > max(LAI_avail_years))
         LAI_used_years[ind_end_yrs] = tail(LAI_avail_years, 1)
         # Print note:
         if (length(c(ind_start_yrs, ind_end_yrs)) > 1) out_string = 'years' else out_string = 'year'
         print(paste0('Simulation ', out_string, ' ', paste(year_vec[c(ind_start_yrs, ind_end_yrs)], collapse = ', '), ' does not have the corresponding LAI data so closest available year(s) are used.'), quote = FALSE)
      }
      
      # Rescale LAI_day_PFT with user-defined LAI for years required:
      LAI_required_years = unique(LAI_used_years)
      if (length(LAI_required_years) > 1) out_string = 'years' else out_string = 'year'
      print(paste0('Daily SAI and LAI data for ', out_string, ' ', paste(LAI_required_years, collapse = ', '), ' are needed for this simulation.'), quote=FALSE)
      
      for (iyear in seq_along(LAI_required_years)) {
         
         output_subfn = paste0('daily_monthly_LAI_data_', as.character(LAI_required_years[iyear]), '.RData')
         output_filename = paste0(processed_surf_data_dir, output_subfn)
         
         # Check if scaled LAI is already available:
         if (!file.exists(output_filename)) {
            
            # Get corresponding LAI nc file:
            filename = paste0(LAI_data_dir, LAI_subn1, as.character(LAI_required_years[iyear]), LAI_subn2)
            print(paste0('User-defined LAI data from ', filename, ' are used.'), quote=FALSE)
            nc = nc_open(filename)
            # Extract variables from user-defined LAI nc file:
            LAI_data_mon = ncvar_get(nc, LAI_array_name)[,,ind_LAI_time]
            lon_LAI = ncvar_get(nc, unname(LAI_dim_vec['longitude']))
            lat_LAI = ncvar_get(nc, unname(LAI_dim_vec['latitude']))
            nc_close(nc)
            
            timestamp()
            
            # Regrid LAI to simulation lon/lat:
            print(paste0('Regridding user-defined LAI for year ', LAI_required_years[iyear], '...'), quote=FALSE)
            LAI_data_mon = NaN.to.zero(sp.regrid(spdata=zero.to.NaN(LAI_data_mon), lon.in=lon_LAI, lat.in=lat_LAI, lon.out=lon, lat.out=lat))
            
            # Rescale CLM PFT-level monthly LAI according to user-defined LAI:
            # Define scaling factor:
            Z = LAI_data_mon/zero.to.NaN(apply(LAI_mon_PFT*array(data = PFT_frac, dim = c(dim(PFT_frac), 12)), MARGIN = c(1,2,4), FUN = sum))
            Z[which(is.na(Z))] = 1
            Z[which(Z > 10)] = 10   # Max it out at a factor of 10
            LAI_scaling_array = Z
            # Rescale:
            LAI_data_mon_PFT = array(data = 0, dim = dim(LAI_mon_PFT))
            # SAI_data_mon_PFT = array(data = 0, dim = dim(SAI_mon_PFT))
            for (ipft in 1:length(pftnum)) {
               LAI_data_mon_PFT[,,ipft,] = LAI_mon_PFT[,,ipft,]*LAI_scaling_array
               # Rescale SAI as well (not turned on):
               # SAI_data_mon_PFT[,,ipft,] = NaN.to.zero(SAI_mon_PFT[,,ipft,]*LAI_scaling_array)
            }
            
            # Interpolate PFT-level monthly LAI and SAI to obtain daily LAI and SAI:
            day = 1:365
            day_midmon = c(-16, 15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349, 380)
            LAI_day_PFT = array(NaN, dim=c(length(lon), length(lat), length(pftname), length(day)))
            SAI_day_PFT = array(NaN, dim=c(length(lon), length(lat), length(pftname), length(day)))
            print(paste0('Interpolating daily LAI and SAI values for year ', LAI_required_years[iyear], '...'), quote=FALSE)
            # Looping over lon/lat/PFT:
            for (i in 1:length(lon)) {
               for (j in 1:length(lat)) {
                  for (ipft in 1:length(pftname)) {
                     # LAI:
                     subdata = c(LAI_data_mon_PFT[i,j,ipft,12], LAI_data_mon_PFT[i,j,ipft,], LAI_data_mon_PFT[i,j,ipft,1])
                     interpol = approx(x=day_midmon, y=subdata, xout=day, method='linear')
                     LAI_day_PFT[i,j,ipft,] = interpol$y
                     # SAI:
                     subdata = c(SAI_mon_PFT[i,j,ipft,12], SAI_mon_PFT[i,j,ipft,], SAI_mon_PFT[i,j,ipft,1])
                     # subdata = c(SAI_data_mon_PFT[i,j,ipft,12], SAI_data_mon_PFT[i,j,ipft,], SAI_data_mon_PFT[i,j,ipft,1])
                     interpol = approx(x=day_midmon, y=subdata, xout=day, method='linear')
                     SAI_day_PFT[i,j,ipft,] = interpol$y
                  }
               }
            }
            LAI_day_PFT = NaN.to.zero(LAI_day_PFT)
            SAI_day_PFT = NaN.to.zero(SAI_day_PFT)
            
            timestamp()
            
            # Save variables and parameters:
            save(list=c('LAI_data_mon', 'LAI_scaling_array', 'LAI_day_PFT', 'SAI_day_PFT', 'lon', 'lat', 'pftname', 'pftnum', 'day'), file=output_filename)
            
         }
         
      }
      
      # Load the interpolated daily LAI and SAI values for the first simulation year:
      output_subfn = paste0('daily_monthly_LAI_data_', as.character(LAI_used_years[1]), '.RData')
      output_filename = paste0(processed_surf_data_dir, output_subfn)
      print(paste0('Loading existing daily LAI and SAI data for year ', as.character(LAI_used_years[1]), '...'), quote=FALSE)
      load(output_filename)
      
      exist_LAI_data = TRUE
      
   }
   
} else exist_LAI_data = FALSE

################################################################################

# Linearly interpolate CLM PFT-level monthly LAI and SAI to obtain daily LAI and SAI if no user-defined LAI data are available:

if (!exist_LAI_data && !biogeochem_flag) {
   
   print("No user-defined LAI data are available so default CLM 2000 LAI and SAI data are used.", quote=FALSE)
   
   subfn = 'daily_LAI_CLM_2000.RData'
   filename = paste0(processed_surf_data_dir, subfn)
   
   if (!file.exists(filename)) {
      
      day = 1:365
      day_midmon = c(-16, 15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349, 380)
      LAI_day_PFT = array(NaN, dim=c(length(lon), length(lat), length(pftname), length(day)))
      SAI_day_PFT = array(NaN, dim=c(length(lon), length(lat), length(pftname), length(day)))
      print(paste0('Interpolating daily LAI and SAI values...'), quote=FALSE)
      for (i in 1:length(lon)) {
         for (j in 1:length(lat)) {
            for (ipft in 1:length(pftname)) {
               # LAI:
               subdata = c(LAI_mon_PFT[i,j,ipft,12], LAI_mon_PFT[i,j,ipft,], LAI_mon_PFT[i,j,ipft,1])
               interpol = approx(x=day_midmon, y=subdata, xout=day, method='linear')
               LAI_day_PFT[i,j,ipft,] = interpol$y
               # SAI:
               subdata = c(SAI_mon_PFT[i,j,ipft,12], SAI_mon_PFT[i,j,ipft,], SAI_mon_PFT[i,j,ipft,1])
               interpol = approx(x=day_midmon, y=subdata, xout=day, method='linear')
               SAI_day_PFT[i,j,ipft,] = interpol$y
            }
         }
      }
      LAI_day_PFT = NaN.to.zero(LAI_day_PFT)
      SAI_day_PFT = NaN.to.zero(SAI_day_PFT)
      
      save(list=c('LAI_day_PFT', 'SAI_day_PFT', 'lon', 'lat', 'pftname', 'pftnum', 'day'), file=filename)
      
   } else {
      if (!exists('LAI_day_PFT') | !exists('SAI_day_PFT')) {
         print('Loading daily LAI and SAI data...', quote=FALSE)
         load(filename)
      }
   }
}

################################################################################
### Load hourly surface ozone concentration field:
################################################################################

# Moved to this file from "execution_vX.Y.R" (Tai, Feb 2019):

# if (O3_damage_flag & !O3_fixed_flag) {
#    
#    # Vector of simulation years:
#    # This is also set above if LAI_data_flag=TRUE.
#    if (!exists('year_vec')) year_vec = as.numeric(unique(substr(make.date.vec(start.date=start_date, end.date=end_date), 1, 4)))
#    
#    # Set surface ozone field path:
#    # Please make sure among the available years for O3 data, there are no missing years in between.
#    O3_nc_files = list.files(O3_data_dir, pattern = paste0('^', O3_subn1, '.*', O3_subn2, '$'))
#    O3_avail_years = sort(as.numeric(str_remove(str_remove(O3_nc_files, pattern = paste0('^', O3_subn1)), pattern = paste0(O3_subn2, '$'))))
#    print(paste0('Available years for O3 data: ', paste(O3_avail_years, collapse = ', ')), quote = FALSE)
#    
#    if (length(O3_avail_years) == 0) stop('Surface O3 data are not found!') else {
#       # Match simulation years to O3 data available:
#       O3_used_years = year_vec
#       # Checking if O3 data are available for the simulation year:
#       if (min(year_vec) < min(O3_avail_years) | max(year_vec) > max(O3_avail_years)) {
#          ind_start_yrs = which(year_vec < min(O3_avail_years))
#          O3_used_years[ind_start_yrs] = O3_avail_years[1]
#          ind_end_yrs = which(year_vec > max(O3_avail_years))
#          O3_used_years[ind_end_yrs] = tail(O3_avail_years, 1)
#          # Print note:
#          if (length(c(ind_start_yrs, ind_end_yrs)) > 1) out_string = 'years' else out_string = 'year'
#          print(paste0('Simulation ', out_string, ' ', paste(year_vec[c(ind_start_yrs, ind_end_yrs)], collapse = ', '), ' does not have the corresponding O3 data so closest available year(s) are used.'), quote = FALSE)
#       }
#    }
#    
#    # Get O3 nc file for the first simulation year:
#    # O3 fields for subsequent simulation years will be read in at the beginning of each year.
#    
#    # Load hourly ozone field:
#    filename = paste0(O3_data_dir, O3_subn1, as.character(O3_used_years[1]), O3_subn2)
#    print(paste0('Loading surface O3 concentrations from ', filename, '...'), quote=FALSE)
#    nc = nc_open(filename)
#    lon_O3 = ncvar_get(nc, unname(O3_dim_vec['longitude']))
#    lat_O3 = ncvar_get(nc, unname(O3_dim_vec['latitude']))
#    # Surface O3 concentration (ppbv):
#    O3_hourly = ncvar_get(nc, O3_array_name)
#    nc_close(nc)
#    
#    # Regrid to model resolution if input resolution is not consistent:
#    if (sum(lon != lon_O3) > 0 | sum(lat[2:(length(lat)-1)] != lat_O3[2:(length(lat_O3)-1)]) > 0) {
#       # Regrid to model resolution:
#       print('Regridding hourly O3 concentrations for year ', as.character(O3_used_years[1]), '...', quote=FALSE)
#       O3_hourly = sp.regrid(spdata=O3_hourly, lon.in=lon_O3, lat.in=lat_O3, lon.out=lon, lat.out=lat)
#    }
#    
# }

################################################################################
### Wesely dry deposition module constants:
################################################################################

if (!is.null(drydep_scheme)) {
   
   if (drydep_scheme == 'Wesely') {
      
      print('Loading Wesely dry deposition parameters...', quote=FALSE)
      
      filename = paste0(TEMIR_dir, 'TEMIR_inputs/Wesely_const/Olson_and_Biome_Drydep_Inputs_PFT_24PFT.nc')
      
      # PFT parameters:
      nc = nc_open(filename)
      # Mapping from CLM landtype to Wesely dry deposition landtype
      WESELY_MAPPING = ncvar_get(nc, 'IDEPB')
      # RCLO
      RCLO_table = ncvar_get(nc, 'IRCLO')
      # RCLS
      RCLS_table = ncvar_get(nc, 'IRCLS')
      # RGSO
      RGO_table = ncvar_get(nc, 'IRGSO')
      # RGS
      RGS_table = ncvar_get(nc, 'IRGSS')
      # RI
      RI_table = ncvar_get(nc, 'IRI')
      # RLU
      RLU_table = ncvar_get(nc, 'IRLU')
      # RAC
      RAC_table = ncvar_get(nc, 'IRAC')
      # DRYCOEFF
      DRYCOEFF_table = ncvar_get(nc, 'DRYCOEFF')
      
      nc_close(nc)
      
   }
   
}

################################################################################
### Zhang dry deposition module constants:
################################################################################

if (!is.null(drydep_scheme)) {
   
   if (drydep_scheme == 'Zhang') {
      
      print('Loading Zhang dry deposition parameters...', quote=FALSE)
      
      filename = paste0(TEMIR_dir, 'TEMIR_inputs/Wesely_const/Zhang_drydep_param_table_ssh.nc')
      
      # PFT parameters:
      nc = nc_open(filename)
      # Mapping from CLM landtype to Wesely dry deposition landtype
      mapping_TEMIR = ncvar_get(nc, 'mapping_TEMIR')
      # b_r for g(PAR)
      b_rs = ncvar_get(nc, 'b_rs')
      # b_vpd for g(VPD)
      b_vpd = ncvar_get(nc, 'bvpd')
      # psi_max for g(psi)
      psi_max = ncvar_get(nc, 'psi_c1')
      # psi_max for g(psi)
      psi_min = ncvar_get(nc, 'psi_c2')
      # r_cut_dO
      r_cut_dO = ncvar_get(nc, 'Rcutd0_O3')
      # r_cut_WO
      r_cut_wO = ncvar_get(nc, 'Rcutw0_O3')
      # rg_o
      rg_o = ncvar_get(nc, 'Rg_O3')
      # rg_s for g(PAR)  ### NOTE: based on temperature for some LUC
      rg_s = ncvar_get(nc, 'Rgd_SO2')
      # rst_min
      rst_min = ncvar_get(nc, 'rstmin')
      # maximum snow depth that the whole canopy is coverred by snow
      sdmax = ncvar_get(nc, 'Sdmax')
      # maximum temperature for g(T)
      t_max_st = ncvar_get(nc, 'Tmax')
      # minimum temperature for g(T)
      t_min_st = ncvar_get(nc, 'Tmin')
      # optimum temperature for g(T)
      t_opt_st = ncvar_get(nc, 'Topt')
      # in-canopy maximum aerodynamic resistance
      Rac0_min = ncvar_get(nc, 'Rac0_min')
      # in-canopy minimum aerodynamic resistance
      Rac0_max = ncvar_get(nc, 'Rac0_max')
      
      nc_close(nc)
      
   }
   
}

################################################################################
### Ratio between grid-level displacement height and surface roughness:
################################################################################

disp_on_z0m = array(NaN, dim=c(length(lon), length(lat)))
for (i in 1:length(lon)) {
   for (j in 1:length(lat)) {
      disp_on_z0m[i,j] = sum(displar * htop * PFT_frac[i,j,], na.rm=TRUE) / sum(z0mr * htop * PFT_frac[i,j,], na.rm=TRUE)
   }
}
disp_on_z0m[which(is.na(disp_on_z0m))] = 0

################################################################################
### Biogeochemistry starting condition
################################################################################

if (biogeochem_flag) {
  
  subfn = 'processed_clmi.ICRUCLM45BGCCROPmp24.0241-01-01.1.9x2.5_g1v6_simyr2000_c130515.nc'
  filename = paste0(initial_data_dir, subfn)
  nc = nc_open(filename)
  LAI_initial_map = ncvar_get(nc, 'tlai'); LAI_initial_map = ifelse(test = is.na(LAI_initial_map), yes = 0, no = LAI_initial_map)
  SAI_initial_map = ncvar_get(nc, 'tsai'); SAI_initial_map = ifelse(test = is.na(SAI_initial_map), yes = 0, no = SAI_initial_map)

  leafC_initial_map = ncvar_get(nc, 'leafC'); leafC_initial_map = ifelse(test = is.na(leafC_initial_map), yes = 0, no = leafC_initial_map)
  frootC_initial_map = ncvar_get(nc, 'frootC'); frootC_initial_map = ifelse(test = is.na(frootC_initial_map), yes = 0, no = frootC_initial_map)
  livestemC_initial_map = ncvar_get(nc, 'livestemC'); livestemC_initial_map = ifelse(test = is.na(livestemC_initial_map), yes = 0, no = livestemC_initial_map)
  deadstemC_initial_map = ncvar_get(nc, 'deadstemC'); deadstemC_initial_map = ifelse(test = is.na(deadstemC_initial_map), yes = 0, no = deadstemC_initial_map)
  livecrootC_initial_map = ncvar_get(nc, 'livecrootC'); livecrootC_initial_map = ifelse(test = is.na(livecrootC_initial_map), yes = 0, no = livecrootC_initial_map)
  deadcrootC_initial_map = ncvar_get(nc, 'deadcrootC'); deadcrootC_initial_map = ifelse(test = is.na(deadcrootC_initial_map), yes = 0, no = deadcrootC_initial_map)
  grainC_initial_map = ncvar_get(nc, 'grainC'); grainC_initial_map = ifelse(test = is.na(grainC_initial_map), yes = 0, no = grainC_initial_map)

  # storage C pools x7 for deciduous
  # leafC_storage_initial_map = ncvar_get(nc, 'leafC_storage')
  # ...

  nc_close(nc)
  print(paste0("Finish loading initial data from ",subfn))
}
################################################################################
### End of module
################################################################################
