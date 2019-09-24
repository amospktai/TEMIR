################################################################################
### Terrestrial Ecosystem Model in R (TEMIR) v1.0
### Initialization script for single-site, regional or global simulation
################################################################################

# TEMIR version:
TEMIR_version = '1.0'

# Set TEMIR directory:
TEMIR_dir = '~/Dropbox/TGABI/TEMIR/'

# Set simulation parent directory:
# Users can specify their customized parent directory for their simulations here.
sim_parent_dir = TEMIR_dir

################################################################################
### TEMIR simulation setup:
################################################################################

# Create a name for this simulation:
simulation_name = 'test_v1.0'

################################################################################
### Check TEMIR model availability:
################################################################################

# Set simulation parent directory:
sim_dir = paste0(sim_parent_dir, simulation_name, '/')

# Check if TEMIR directory exists:
if (!dir.exists(paths = TEMIR_dir)) stop('TEMIR directory does not exist!')

# Check if TEMIR version exists:
if (!dir.exists(paths = paste0(TEMIR_dir, 'code_v', TEMIR_version))) stop('TEMIR version does not exist!')

################################################################################
### Initialize TEMIR:
################################################################################

# Check if simulation directory already exists:
if (dir.exists(paths = sim_dir)) stop(paste0('Simulation directory "', sim_dir,'" already exists!'))

# Create output directory:
dir.create(path = paste0(sim_dir, 'hist_data'), recursive = TRUE)

# Create temporary directory:
dir.create(path = paste0(sim_dir, 'temp_data'))

# Set simulation directory as working directory:
setwd(sim_dir)

# Copy execution script:
file.copy(from = paste0(TEMIR_dir, 'code_v', TEMIR_version, '/execution_v', TEMIR_version, '.R'), to = sim_dir)

# Copy baseline input script:
file.copy(from = paste0(TEMIR_dir, 'code_v', TEMIR_version, '/input_TEMIR.R'), to = sim_dir)

# Copy script that contains functions to analyze outputs:
file.copy(from = paste0(TEMIR_dir, 'code_v', TEMIR_version, '/find_hist_stat.R'), to = sim_dir)

################################################################################
### End of initialization
################################################################################