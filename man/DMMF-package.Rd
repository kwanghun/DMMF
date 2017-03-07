\name{DMMF-package}
\alias{DMMF-package}
\docType{package}
\title{
Daily Based Morgan--Morgan--Finney (DMMF) Soil Erosion Model
}
\description{
Implements the daily based Morgan--Morgan--Finney (DMMF) soil erosion model for estimating surface runoff and sediments budgets from a field or a catchment on a daily basis.
}
\details{
This package is the implementation of the daily based Morgan--Morgan--Finney (DMMF) soil erosion model to estimate surface runoff and sediment budget of a field and a catchment on a daily basis.

The DMMF model is one of the variant of the widely used Morgan--Morgan--Finney soil erosion model and largely based on the modified MMF model with several modifications (see details in Choi et al. (2017)). 

This \R implementation of the DMMF model is suitable for estimating surface runoff and sediment budgets of fields or catchments represented by raster.

The package provides the \code{\link{DMMF}} function that estimates surface runoff and soil erosion using DMMF model, \code{\link{DMMF_Simple}} function which is the simpler version of the \code{\link{DMMF}} and \code{\link{SinkFill}} function that generates sink free map for the DMMF model based on the sink fill algorithm from Wang and Liu (2006) and the SAGA-GIS module of ``Module Fill Sinks (Wang \& Liu)'' by Volker Wichmann (2007).
}
\author{
    Kwanghun Choi, Sebastian Arnhold, Bjoern Reineking

Maintainer: Kwanghun Choi (\email{kwanghun.choi@yahoo.com}) 
}
\references{
    Choi, K., Arnhold, S., Huwe, B., Reineking, B. (2017). 
    Daily based Morgan--Morgan--Finney (DMMF) model: A spatially distributed conceptual soil erosion model to simulate complex soil surface configurations. \emph{Submitted to Water in Feb, 2017}. 

    Shin, S. S., Park, S. D., and Choi, B. K. (2016). 
    Universal power law for relationship between rainfall kinetic energy and rainfall intensity. \emph{Advances in Meteorology}, Article ID 2494681, 11 pages.

    Arnhold, S., Ruidisch, M., Bartsch, S., Shope, C., Huwe, B. (2013). 
    Simulation of runoff patterns and soil erosion on mountainous farmland with and without plastic-covered ridge-furrow cultivation in South Korea. \emph{Transactions of the ASABE}, 56(2):667--679.

    Morgan, R. P. C. and Duzant, J. H. (2008). 
    Modified MMF (Morgan--Morgan--Finney) model for evaluating effects of crops and vegetation cover on soil erosion. \emph{Earth Surface Processes and Landforms}, 32(1):90--106.

    Seibert, J., McGlynn, B. L. (2007). 
    A new triangular multiple flow direction algorithm for computing upslope areas from gridded digital elevation models. \emph{Water Resources Research}, 43(4):W04501. W04501.

    Volker Wichmann (2007). 
    Module Fill Sinks (Wang \& Liu). SAGA-GIS Module Library Documentation (v2.1.3) [ cited 2015. 08. 15 ], Available from: \url{http://www.saga-gis.org/saga_tool_doc/2.1.3/ta_preprocessor_4.html}

    Wang, L. and Liu, H. (2006). 
    An efficient method for identifying and filling surface depressions in digital elevation models for hydrologic analysis and modelling. \emph{International Journal of Geographical Information Science}, 20(2):193--213.

    Morgan, R. P. C. (2001) 
    A simple approach to soil loss prediction: a revised Morgan--Morgan--Finney model. \emph{Catena}, 44(4):305--322.

    Morgan, R. P. C., Morgan, D. D. V., and Finney, H. J. (1984) 
    A predictive model for the assessment of soil erosion risk. \emph{Journal of Agricultural Engineering Research}, 30(0):245--253.

}
\keyword{Soil erosion, surface runoff, Hydrology, MMF, DMMF, Hillslope}
\seealso{
    \code{\link{DMMF}} to run the DMMF model, 
    \code{\link{DMMF_Simple}} to run the DMMF model for simple one element during a day,
    \code{\link{SinkFill}} to generate sink free map for the DMMF model
}
