################################################################################
### Terrestrial Ecosystem Model in R (TEMIR) v1.1
### Initialization script for single-site, regional or global simulation
################################################################################

# TEMIR version:
TEMIR_version = '1.1'

# Set TEMIR directory:
TEMIR_dir = '/users/b146986/TEMIR/'

################################################################################
### TEMIR simulation setup:
################################################################################

# Case name read from bash script (temporary)
TEMIR_caseName = as.character(Sys.getenv("bash_caseName"))

# Create a name for this simulation:
simulation_name = TEMIR_caseName
# simulation_name = 'soyFACE_set2_forceSP_yr2003_eco2'

# Simulation types:

# Simulating biogeochemistry? (default: FALSE)
# If TRUE, leaf area index and stem area index are now prognostic instead of prescribed.
biogeochem_flag = TRUE

################################################################################
### Check TEMIR model availability:
################################################################################

# Check if TEMIR directory exists:
if (!dir.exists(paths = TEMIR_dir)) stop('TEMIR directory does not exist!')

# Check if TEMIR version exists:
if (!dir.exists(paths = paste0(TEMIR_dir, '/code_v', TEMIR_version))) stop('TEMIR version does not exist!')

################################################################################
### Initialize TEMIR:
################################################################################

# Set simulation directory name:
sim_dir = paste0(TEMIR_dir,'TEMIR_run/' ,simulation_name, '/')

# Check if simulation directory already exists:
# if (length(Sys.glob(paths = sim_dir)) == 1) stop(paste0('Simulation directory "', sim_dir,'" already exists!'))
if (dir.exists(paths = sim_dir)) stop(paste0('Simulation directory "', sim_dir,'" already exists!'))
dir.create(path = sim_dir)

# Create output directory:
# system(command=paste0("mkdir -p ", sim_dir, "hist_data"))
dir.create(path = paste0(sim_dir, 'hist_data'), recursive = TRUE)

# Create temporary directory:
# system(command=paste0("mkdir ", sim_dir, "temp_data"))
dir.create(path = paste0(sim_dir, 'temp_data'))

# Set simulation directory as working directory:
setwd(sim_dir)

# Copy execution script:
# system(command=paste0("cp ", TEMIR_dir, 'code_v', TEMIR_version, '/execution_v', TEMIR_version, ".R ", sim_dir))
file.copy(from = paste0(TEMIR_dir, 'code_v', TEMIR_version, '/execution_v', TEMIR_version, '.R'), to = sim_dir)

# Copy basic settings and ecophysiology input script:
file.copy(from = paste0(TEMIR_dir, 'code_v', TEMIR_version, '/input_TEMIR_basic_settings.R'), to = sim_dir)

# Copy biogeochemistry input script:
if (biogeochem_flag) file.copy(from = paste0(TEMIR_dir, 'code_v', TEMIR_version, '/input_TEMIR_biogeochem_extension.R'), to = sim_dir)

# Copy simulation name text file:
file.copy(from = paste0(TEMIR_dir, 'code_v', TEMIR_version, '/simulation_name.txt'), to = paste0(sim_dir, simulation_name, '.txt'))

# Copy script that contains functions to analyze outputs:
file.copy(from = paste0(TEMIR_dir, 'code_v', TEMIR_version, '/find_hist_stat.R'), to = sim_dir)


################################################################################
### End of initialization
################################################################################
