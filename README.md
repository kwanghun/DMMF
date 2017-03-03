## DMMF: Daily Based Morgan--Morgan--Finney (DMMF) Soil Erosion Model

#### Package description
Implements the daily based Morgan--Morgan--Finney (DMMF) soil erosion model for estimating surface runoff and sediments budgets from a field or a catchment on a daily basis.
#### Details
This package is the implementation of the daily based Morgan--Morgan--Finney (DMMF) soil erosion model to estimate surface runoff and sediment budget of a field and a catchment on a daily basis.
The DMMF model is one of the variant of the widely used Morgan--Morgan--Finney soil erosion model and largely based on the modified MMF model with several modifications (see details in Choi et al. (2017)). 
This R implementation of the DMMF model is suitable for estimating surface runoff and sediment budgets of fields or catchments represented by raster. The package provides the `DMMF` function that estimates surface runoff and soil erosion using DMMF model, `DMMF_Simple` function which is the simpler version of `DMMF`, and `SinkFill` function that generates sink free map for the DMMF model based on the sink fill algorithm from Wang and Liu (2006) and the SAGA-GIS module of "Module Fill Sinks (Wang & Liu)" by Volker Wichmann (2007).
#### Package maintainer
Kwanghun Choi <kwanghun.choi@yahoo.com>
