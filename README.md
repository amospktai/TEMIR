# Terrestrial Ecosystem Model in R (TEMIR)

## Introduction
The **Terrestrial Ecosystem Model in R (TEMIR)** is a ecosystem model developed by the [Tai Group for Atmosphere-Biosphere Interactions](http://www.cuhk.edu.hk/sci/essc/tgabi/index.html). It computes ecophysiological processes and responses of local and global terrestrial ecosystems (e.g., canopy radiative transfer, photosynthesis, conductance, dry deposition of air pollutants) driven by prescribed meteorological and land surface input data. The primary purpose is to evaluate changes in ecosystem functions in response to changes in the terrestrial and atmospheric environment (e.g., CO<sub>2</sub>, ozone, temperature, humidity, soil moisture, plant type distribution), and to evaluate how changes in ecosystem functions influence atmospheric chemistry and climate.

The default meteorological inputs for TEMIR are from NASA/GAMO [MERRA2](http://wiki.seas.harvard.edu/geos-chem/index.php/List_of_MERRA-2_met_fields) and [GEOS-FP](http://wiki.seas.harvard.edu/geos-chem/index.php/List_of_GEOS-FP_met_fields) reanalysis products, so that ecosystem simulations are entirely consistent and can be asynchronously coupled with atmospheric chemistry simulations by the [GEOS-Chem](http://wiki.seas.harvard.edu/geos-chem/index.php/Main_Page) chemical transport model. Land surface (e.g., plant function type, soil type) inputs are mostly from the [Community Land Model (CLM)](http://www.cesm.ucar.edu/models/clm/) version 4.5, with 24 plant function types (PFTs) including several individual midlatitude crops. The model is, however, highly customizable, and users can liberally replace the default inputs with their own data. The model can be run parallelly on multiple cores, and computational speed is scalable almost linearly to the number of cores used.

## Version and history (Last update: Jun, 2019)
- July 2019, v1.1a with ... etc.
- June 2019, v1.1 with dry drposition

## Future development
### Biogenic emission
### Soil biogeochehmistry, N cycle

## Getting start
The following components are required to run the model.
- Model code
- Required input data
- MERRA2 meteorological fields
- GEOS-Chem-simulated hourly ozone concentrations (optional; to be put into "/TEMIR_input/O3_data/"; see [Fu and Tai, 2015](https://www.atmos-chem-phys.net/15/10093/2015/acp-15-10093-2015.html))
- Monthly leaf area index (LAI) used by GEOS-Chem (optional; to be put into "/TEMIR_inputs/LAI_data/"; see [here](http://wiki.seas.harvard.edu/geos-chem/index.php/Leaf_area_indices_in_GEOS-Chem))

## Publications
- Tai and Yung (2019) GMD ...
- Pang and Tai (2019) GMD ...
- Sun et al. (2019) ...
