# Resilience Sensing of Forests with Landsat

A methodological framework specifically designed for the computation of early warning indicators.
It allows to generate Landsat composites over large extents of forest ecosystems, while reducing sensor discontinuities, 
enhancing spatial consistency and preserving ecologically meaningful temporal  dynamics to enable reproducible, robust and 
comparable assessments of forest resilience across biomes.

## Landsat composite generator in GEE (Google Earth Egine)

The available code in GEE for generating high-quality, cross-calibrated Landsat composites.
It allows for ... that cover large extents of terrestrial biomes

https://code.earthengine.google.com/c22d0250e5741dfe7f1351351f52aa12

![General criteria to address practical and multispectral limitations](https://raw.githubusercontent.com/camuradr/LandsatCompositesGEE/refs/heads/main/Data/Figures/GEEworkflow_v2.svg?token=GHSAT0AAAAAAEJLF7OX5W5Y6QOOKYMWC5YU2WDR55Q)

## Computing STEWS (RStudio)

Code to manage the raster outputs from GEE exports to generate the composite dataset.
It includes the code for computing several vegetation indices when surface reflectance bands are exported.
It also allows for the computation of the spatial early warning signals (SEWS) through a parallel processing plan.

![Multiscale grid approach to tackle the scale limitations ](https://raw.githubusercontent.com/camuradr/LandsatCompositesGEE/refs/heads/main/Data/Figures/HDGGworkflow_v2.svg?token=GHSAT0AAAAAAEJLF7OXLQH65WLX4BXCXXMQ2WDR7OA)