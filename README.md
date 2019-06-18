# Terrestrial Ecosystem Model in R (TEMIR)

## Introduction
The **Terrestrial Ecosystem Model in R (TEMIR)** is a ecosystem model developed by the [Tai Group for Atmosphere-Biosphere Interactions](http://www.cuhk.edu.hk/sci/essc/tgabi/index.html). It computes ecophysiological processes and responses of local and global terrestrial ecosystems (e.g., canopy radiative transfer, photosynthesis, conductance, dry deposition of air pollutants) driven by prescribed meteorological and land surface input data. The primary purpose is to evaluate changes in ecosystem functions in response to changes in the terrestrial and atmospheric environment (e.g., CO<sub>2</sub>, ozone, temperature, humidity, soil moisture, plant type distribution), and to evaluate how changes in ecosystem functions influence atmospheric chemistry and climate.

The default meteorological inputs for TEMIR are from NASA/GAMO [MERRA2](http://wiki.seas.harvard.edu/geos-chem/index.php/List_of_MERRA-2_met_fields) and [GEOS-FP](http://wiki.seas.harvard.edu/geos-chem/index.php/List_of_GEOS-FP_met_fields) reanalysis products, so that ecosystem simulations are entirely consistent and can be asynchronously coupled with atmospheric chemistry simulations by the [GEOS-Chem](http://wiki.seas.harvard.edu/geos-chem/index.php/Main_Page) chemical transport model. Land surface (e.g., plant function type, soil type) inputs are mostly from the [Community Land Model (CLM)](http://www.cesm.ucar.edu/models/clm/) version 4.5, with 24 plant function types (PFTs) including several individual midlatitude crops. The model is, however, highly customizable, and users can liberally replace the default inputs with their own data. The model can be run parallelly on multiple cores, and computational speed is scalable almost linearly to the number of cores used.

## Version and history (Last update: June, 2019)
- Jul 2019,	 v1.1a	; soil water stress representation with double layers ....
- 18 Jun 2019,	 v1.1	; add an ozone dry deposition module and a crop module (in branch *biogeochem* only) (First release on GitHub!)

## Future development
### Biogenic emission
### Soil biogeochehmistry, N cycle ( in branch *biogeochem*)
### Natural vegetation (in branch *biogeochem*)

## Getting start
The following components are required to run the model. You can follow *TEMIR_manual_v1.1* after cloning the repository.
#### R
The latest version of R can be downloaded [here](https://www.r-project.org/), the minimum requirement is 3.3.1.
#### Model code
To download the code, simply use `git clone` at your local git repository
```
git init
git clone https://github.com/amospktai/TEMIR.git
```
#### Required input data
The required input can be downloaded [here](https://www.dropbox.com/sh/qbcpcrc9f71ks2x/AACuI1iAiBamp4vP7lw2eubRa?dl=0). It includes:
- PFT input from '/clm2_data/pftdata/pft-physiology.c130503.nc'
- Land surface data from '/clm2_data/surfdata_map/surfdata_1.9x2.5_mp24_simyr2000_c130419.nc'
- [FLUXNET](https://fluxnet.fluxdata.org/) input data from '/FLUXNET/'
- Meteorological inputs for two days and a MERRA2 constant field input from '/met_data/', more data is available below
- Input for the dry deposition module from '/Wesely_const/'
- Processed land surface data input in '/processed_surf_data' (It saves time from regirdding the land surface data if you run with the default 1.9x2.5 resolution)
- Two empty directories 'LAI_data' and 'O3_data', put your LAI and ozone field in these directories if you want to run with custom LAI input and ozone-plant damage respectively.
#### MERRA2 and GEOS-FP meteorological fields (very huge!)
Meteorological inputs can be downloaded from [here](https://www.dropbox.com/sh/dqqs4g8lhdu2z27/AAB9WkkzPYSElD2ouur-yz3ba?dl=0). OR
You can download the latest inputs by following the tutorial in the [GEOS-Chem website](http://wiki.seas.harvard.edu/geos-chem/index.php/Downloading_GEOS-Chem_data_directories). Noted that we only use the A1 field as the input.
#### GEOS-Chem-simulated hourly ozone concentrations (optional; to be put into "/TEMIR_input/O3_data/"; see [Fu and Tai, 2015](https://www.atmos-chem-phys.net/15/10093/2015/acp-15-10093-2015.html))
#### Monthly leaf area index (LAI) used by GEOS-Chem (optional; to be put into "/TEMIR_inputs/LAI_data/"; see [here](http://wiki.seas.harvard.edu/geos-chem/index.php/Leaf_area_indices_in_GEOS-Chem))


## Publications
- Tai and Yung (2019) GMD ...
- Pang and Tai (2019) GMD ...
- Sun et al. (2019) ...

## Other references
- [*Ecological Climatology: Concepts and Applications (3<sup>rd</sup> Ed)*](http://www.cambridge.org/catalogue/catalogue.asp?isbn=9781107043770) by Gordon Bonan
- [Technical note of CLM4.5](https://www.dropbox.com/s/tlno1pbh4axdn6f/CLM45_Tech_Note.pdf?dl=0)
- [Dry deposition module description](http://wiki.seas.harvard.edu/geos-chem/index.php/Dry_deposition)
