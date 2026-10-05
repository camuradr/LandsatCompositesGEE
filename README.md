# Resilience Sensing of Forests with Landsat

A methodological framework specifically designed for the computation of early warning indicators.
It allows to generate Landsat composites over large extents of forest ecosystems, while reducing sensor discontinuities, 
enhancing spatial consistency and preserving ecologically meaningful temporal  dynamics to enable reproducible, robust and 
comparable assessments of forest resilience across biomes.

![](https://github.com/camuradr/LandsatCompositesGEE/blob/main/Data/Figures/GraphicalAbstract.svg)

##

## Landsat composite generator in GEE (Google Earth Egine)

The GEE workflow integrates robust processing and a consistent sensor cross-calibration to improve the reliability of the Landsat imagery used for resilience analysis. 

![General criteria to address practical and multispectral limitations](https://github.com/camuradr/LandsatCompositesGEE/blob/main/Data/Figures/GEEworkflow.svg)

##

The link below corresponds to the scrip available in GEE for generating high-quality, cross-calibrated Landsat composites.It can also be found
in the Scripts folder as a .js file. It allows for the exploration and generation of Landsat composites that cover large extents of terrestrial biomes. 

https://code.earthengine.google.com/c22d0250e5741dfe7f1351351f52aa12

Input paramaters can be modified according to the preferences, which include the criteria outlined in the framework. 

## Computing STEWS (RStudio)

The rasters exported from GEE can be processed with RStudio to compute early warning signals using available packages.
Several vegetation indices can be derived from the exported surface reflectance bands inside R, which is done within the script.
It also allows for the computation of the spatial early warning signals (SEWS) through a parallel processing plan.

![Multiscale grid approach to tackle the scale limitations ](https://github.com/camuradr/LandsatCompositesGEE/blob/main/Data/Figures/HDGGworkflow.png)

##

The code to manage the raster outputs from GEE exports to generate the composite dataset is available inside the Scripts folder.
The instructions, along with examples, for generating the hexagonal discrete global grids are inside the HDGGs folder.