\name{Potato.Concave}
\alias{Potato.Concave}
\title{Example data from the concave potato field from Arnhold et al. (2013)}
\concept{datasets}
\description{
    Example input data for running the \code{\link{DMMF}} function. This data is based on the field data from the concave potato field described in Arnhold et al. (2013).with optimized parameters described in Choi et al. (2017).
}   
\usage{
Potato.Concave
}
\value{
    \code{Potato.Concave} contains input data for the \code{\link{DMMF}} including map type and numeric type variables. Static variables are fixed and do not change with time (or period). Dynamic variables change with time (or period).\cr
    \code{s.map}: Static input maps for the \code{\link{DMMF}} model.\cr
    \itemize{
        \item \code{DEM.original}: RasterLayer object of original digital elevation model (DEM) of the concave potato field before filling sinks.
        \item \code{DEM}: RasterLayer object of DEM after applying sinkfill algorithm.
        \item \code{IMP}: RasterLayer object of the proportion of area covered with impervious surface cover.
        \item \code{n_s}: RasterLayer object of the Manning's roughness coefficient of the soil surface.
    }
    \code{d.map}: Dynamic input maps for the \code{\link{DMMF}} model.\cr
    \itemize{
        \item \code{theta_init}: RasterStack object of the initial soil water content of entire soil profile for each rainfall period.
    }
    \code{s.var}: Static input numeric variables for the \code{\link{DMMF}} model.\cr
    \itemize{
        \item \code{P_c}: Numeric object of proportion of clay particles in surface area.
        \item \code{P_z}: Numeric object of proportion of silt particles in surface area.
        \item \code{P_s}: Numeric object of proportion of sand particles in surface area.
        \item \code{theta_sat}: Numeric object of saturated soil water content of entire soil profile. 
        \item \code{theta_fc}: Numeric object of soil water content at field capacity of entire soil profile. 
        \item \code{SD}: Numeric object of soil depth.
        \item \code{K}: Numeric object of saturated lateral hydraulic conductivity which is obtained from model calibration described in Choi et al. (2017).
        \item \code{P_I}: Numeric object of proportion of permanent interception area of rainfall.  
        \item \code{GC}: Numeric object of proportion of area with pervious vegetated ground cover.
        \item \code{D}: Numeric object of average diameter of individual plant elements at the surface.  
        \item \code{NV}: Numeric object of number of individual plant elements per unit area. 
        \item \code{d_a}: Numeric object of typical flow depth of surface runoff.
        \item \code{DK_c}: Numeric of detachability of clay particles by rainfall.
        \item \code{DK_z}: Numeric of detachability of silt particles by rainfall.
        \item \code{DK_s}: Numeric of detachability of sand particles by rainfall.
        \item \code{DR_c}: Numeric of detachability of clay particles by surface runoff.
        \item \code{DR_z}: Numeric of detachability of silt particles by surface runoff.
        \item \code{DR_s}: Numeric of detachability of sand particles by surface runoff.
    }
    \code{d.var}: Dynamic input numeric variables for the \code{\link{DMMF}} model.\cr
    \itemize{
        \item \code{R}: Numeric vector object of daily rainfall per unit area.
        \item \code{RI}: Numeric vector object of mean hourly rainfall intensity during a day.
        \item \code{ET}: Numeric vector object of evapotranspiration per unit area.
        \item \code{CC}: Numeric vector object of proportion of area with canopy cover.
        \item \code{PH}: Numeric vector or object of average height of vegetation or crop cover where leaf drainage start to fall
        \item \code{Breaking}: Numeric vector of starting points of each simulation period.
        \item \code{Init_point}: Numeric vector of starting points of each rainfall event.
    }
}
\author{
    Kwanghun Choi and Sebastian Arnhold
}
\references{
    Choi, K., Arnhold, S., Huwe, B., Reineking, B. (2017) 
    Daily based Morgan--Morgan--Finney (DMMF) model: A spatially distributed conceptual soil erosion model to simulate complex soil surface configurations. \emph{Submitted to Water in Feb, 2017}. 

    Arnhold, S., Ruidisch, M., Bartsch, S., Shope, C., Huwe, B. (2013). 
    Simulation of runoff patterns and soil erosion on mountainous farmland with and without plastic-covered ridge-furrow cultivation in South Korea. \emph{Transactions of the ASABE}, 56(2):667--679.
}
\seealso{
    \code{\link{DMMF}} for details of variables.
}
