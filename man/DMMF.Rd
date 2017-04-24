\name{DMMF}
\alias{DMMF}
\title{Daily based Morgan--Morgan--Finney (DMMF) soil erosion model}
\concept{soil erosion model}
\description{
    Estimating surface runoff and sediment budget using the Daily based Morgan--Morgan--Finney soil erosion model using the algorithm from Choi et al. (2017)
}
\usage{
DMMF(DEM, R, RI, ET, P_c, P_z, P_s, theta_init, theta_sat, theta_fc, 
     SD, K, P_I, n_s, CC, GC, IMP, PH, D, NV, d_a = 0.005, 
     DK_c = 0.1, DK_z = 0.5, DK_s = 0.3, 
     DR_c = 1.0, DR_z = 1.6, DR_s = 1.5, 
     Breaking, Init_point, Sinks, R_Type = 0, slpMode = 2, ALL = TRUE)
}
\arguments{
    \item{DEM}{RasterLayer object of digital elevation model (DEM)}
    \item{R}{Numeric vector or RasterBrick object of daily rainfall per unit area (unit: mm)}
    \item{RI}{Numeric vector or RasterBrick object of mean hourly rainfall intensity during a day (unit: mm/h)}
   \item{ET}{Numeric vector or RasterBrick object of evapotranspiration per unit area (unit: mm)}
   \item{P_c}{Numeric vector or RasterLayer object of proportion of clay particles in surface soil (0-1)}
   \item{P_z}{Numeric vector or RasterLayer object of proportion of silt particles in surface soil (0-1)}
   \item{P_s}{Numeric vector or RasterLayer object of proportion of sand particles in surface soil (0-1)}
   \item{theta_init}{Numeric vector or RasterBrick object of initial soil water content of entire soil profile per unit area (unit: vol/vol)}
   \item{theta_sat}{Numeric vector or RasterLayer object of saturated soil water content of entire soil profile per unit area (unit: vol/vol)}
   \item{theta_fc}{Numeric vector or RasterLayer object of soil water content at field capacity of entire soil profile per unit area (unit: vol/vol)}
   \item{SD}{Numeric vector or RasterLayer object of soil depth of entire soil profile (unit: m)}
   \item{K}{Numeric vector or RasterLayer object of saturated lateral hydraulic conductivity of entire soil profile (unit: m/d)}
   \item{P_I}{Numeric vector or Raster* object of proportion of permanent interception area of rainfall (0-1)}
   \item{n_s}{Numeric vector or Raster* object of Manning's roughness coefficient of the soil surface (unit: \eqn{\mathrm{s}/{\mathrm{m}^{1/3}}}{s/m^(1/3)})}
   \item{CC}{Numeric vector or Raster* object of proportion of area with canopy cover (0-1)}
   \item{GC}{Numeric vector or Raster* object of proportion of area with pervious vegetated ground cover (0-1)}
   \item{IMP}{Numeric vector or Raster* object of proportion of area with impervious ground cover (0-1)}
   \item{PH}{Numeric vector or Raster* object of average height of vegetation or crop cover where leaf drainage start to fall (unit: m)}
   \item{D}{Numeric vector or Raster* object of average diameter of individual plant elements at the surface (unit: m)}
   \item{NV}{Numeric vector or Raster* object of number of individual plant elements per unit area (unit: \eqn{\mathrm{number}/\mathrm{m^2}}{number/m^2})}
   \item{d_a}{Numeric vector or Raster* object of typical flow depth of surface runoff (unit: m)}
   \item{DK_c}{Numeric vector or RasterLayer object of detachability of clay particles by rainfall (unit: g/J)}
   \item{DK_z}{Numeric vector or RasterLayer object of detachability of silt particles by rainfall (unit: g/J)}
   \item{DK_s}{Numeric vector or RasterLayer object of detachability of sand particles by rainfall (unit: g/J)}
   \item{DR_c}{Numeric vector or RasterLayer object of detachability of clay particles by surface runoff (unit: g/mm)}
   \item{DR_z}{Numeric vector or RasterLayer object of detachability of silt particles by surface runoff (unit: g/mm)}
   \item{DR_s}{Numeric vector or RasterLayer object of detachability of sand particles by surface runoff (unit: g/mm)}
   \item{Breaking}{Numeric vector of starting points of each simulation period (0: during the period points, 1: breaking points), If \code{breaking} is missing, by default the first day is set to 1 and the other days are set to 0.}
   \item{Init_point}{Numeric vector of starting points of each rainfall event (0: during the period points, 1: breaking points), If \code{Init_point} is missing, \code{Init_point} is set to \code{Breaking} by default.}
   \item{Sinks}{RasterLayer object of sinks (e.g., streams, reservoirs, and lakes) of sediments and surface runoff (optional).}
   \item{R_Type}{Integer object of each rainfall type for estimating kinetic energy of direct throughfall (DT) (0-8) (see details)}
   \item{slpMode}{Integer object of slope calculation algorithm (2: second order, 3: third order)}
   \item{ALL}{Logical object of selection for the entire output results}
}
\details{
       \code{R_Type} can be chosen among options appropriate for each regions according to Morgan (2001) and Morgan (2005). Default option adopts recent universal relationships between kinetic energy density and mean hourly rainfall intensity (\code{RI}) from Shin et al. (2016). Default option is recommended when users have measured \code{RI}. Detailed descriptions are described below.
       \tabular{cll}{
        Option\tab Region             \tab Kinetic energy density of throughfall [\eqn{\mathrm{J/m^{2}/mm}}{J/m^2/mm}] \cr
        0     \tab Universal (Default)\tab \eqn{10.3 \cdot \mathrm{RI}^{2/9}}{10.3 * RI^(2/9)} \cr 
        1     \tab North America      \tab \eqn{11.87 + 8.73 \cdot \log_{10}(\mathrm{RI})}{11.87 + 8.73 * log10(RI)} \cr
        2     \tab North-west Europe  \tab \eqn{8.95 + 8.44 \cdot \log_{10}(\mathrm{RI})}{8.95 + 8.44 * log10(RI)} \cr
        3     \tab Mediterranean      \tab \eqn{9.81 + 11.25 \cdot \log_{10}(\mathrm{RI})}{9.81 + 11.25 * log10(RI)} \cr
        4     \tab West Mediterranean \tab \eqn{35.9 \cdot (1.0 - 0.56 \cdot \exp(-0.034 \cdot \mathrm{RI}))}{35.9 * (1.0 - 0.56 * exp(-0.034 * RI))} \cr
        5     \tab Tropics            \tab \eqn{29.8 - ( 127.5 / \mathrm{RI} )}{29.8 - (127.5 / RI)} \cr
        6     \tab East Asia          \tab \eqn{9.81 + 10.60 \cdot \log_{10}( \mathrm{RI} )}{9.81 + 10.60 * log10(RI)} \cr
        7     \tab Temperate Southern hemisphere \tab \eqn{29.0 \cdot ( 1.0 - 0.6 \cdot \exp( -0.04 \cdot \mathrm{RI} ) )}{29.0 * (1.0 - 0.6 * exp(-0.04 * RI))} \cr
        8    \tab Universal          \tab  \eqn{28.3 \cdot ( 1.0 - 0.52 \cdot \exp( -0.042 \cdot \mathrm{RI} ) )}{28.3 * (1.0 - 0.52 * exp(-0.042 * RI))}
       }
       Multi-flow direction algorithm (MD\eqn{\infty}) from Seibert and McGlynn (2007) is built in the \code{\link{DMMF}} model. 
}
\value{
    The output of the function \code{DMMF} is a list of RasterLayer or RasterBrick objects containing the following elements:\cr
    \itemize{
        \item \code{A}: Surface area size of each element considering slope (unit: \eqn{\mathrm{m^2}}{m^2})
        \item \code{Rf}: The amount of effective rainfall reaching each element (unit: \eqn{\mathrm{mm}}{mm})
        \item \code{SW_c}: Surface water infiltration capacity of each element (unit: \eqn{\mathrm{mm}}{mm}
        \item \code{theta_r}: Remaining soil water content of each element (unit: vol/vol)
        \item \code{TC}: Transport capacity of the runoff of each element (unit: \eqn{\mathrm{kg/m^2}}{kg/m^2})
        \item \code{Q_in}: Volume of surface runoff flowing into each element (unit: L)
        \item \code{Q_out}: Volume of surface runoff flowing from each element (unit: L)
        \item \code{IF_in}: Volume of subsurface water flowing into each element (unit: L)
        \item \code{IF_out}: Volume of subsurface water flowing from each element (unit: L)
        \item \code{SS_c}: Area density of clay delivered to surface runoff in unit area of each element (unit: \eqn{\mathrm{kg/m^2}}{kg/m^2})  
        \item \code{SS_z}: Area density of silt delivered to surface runoff in unit area of each element (unit: \eqn{\mathrm{kg/m^2}}{kg/m^2})  
        \item \code{SS_s}: Area density of sand delivered to surface runoff of each element (unit: \eqn{\mathrm{kg/m^2}}{kg/m^2})  
        \item \code{G_c}: Area density of available clay for transport by surface runoff of each element (unit: \eqn{\mathrm{kg/m^2}}{kg/m^2})
        \item \code{G_z}: Area density of available silt for transport by surface runoff of each element (unit: \eqn{\mathrm{kg/m^2}}{kg/m^2})
        \item \code{G_s}: Area density of available sand for transport by surface runoff of each element (unit: \eqn{\mathrm{kg/m^2}}{kg/m^2})
        \item \code{SL_c_in}: Mass of clay inputs into each element (unit: kg)
        \item \code{SL_z_in}: Mass of silt inputs into each element (unit: kg)
        \item \code{SL_s_in}: Mass of sand inputs into each element (unit: kg) 
        \item \code{SL_in}: Mass of sum of total sediments inputs into each element (unit: kg)
        \item \code{SL_c_out}: Mass of clay outputs from each element (unit: kg)
        \item \code{SL_z_out}: Mass of silt outputs from each element (unit: kg)
        \item \code{SL_s_out}: Mass of sand outputs from each element (unit: kg) 
        \item \code{SL_out}: Mass of sum of total sediments outputs from each element (unit: kg)
    }
}
\references{
    Choi, K., Arnhold, S., Huwe, B., Reineking, B. (2017). 
    Daily based Morgan--Morgan--Finney (DMMF) model: A spatially distributed conceptual soil erosion model to simulate complex soil surface configurations. \emph{Water}, 9(4), 278. 

    Shin, S. S., Park, S. D., and Choi, B. K. (2016). 
    Universal power law for relationship between rainfall kinetic energy and rainfall intensity. \emph{Advances in Meteorology}, Article ID 2494681, 11 pages.

    Seibert, J., McGlynn, B. L. (2007). 
    A new triangular multiple flow direction algorithm for computing upslope areas from gridded digital elevation models. \emph{Water Resources Research}, 43(4):W04501.

    Morgan, R. P. C. (2005). 
    Soil erosion and conservation. Blackwell Publishing, Malden, MA, 3rd ed. ISBN 1-4051-1781-8

    Morgan, R. P. C. (2001)
    A simple approach to soil loss prediction: a revised Morgan--Morgan--Finney model. \emph{Catena}, 44(4):305--322.
}
\author{
    Kwanghun Choi and Bjoern Reineking
}
\examples{
\dontrun{
## Load example data for test running DMMF model
data(Potato.Convex)
attach(Potato.Convex)
## Run DMMF function using Potato.Convex data
Result <- DMMF(DEM = s.map$DEM, R = d.var$R, RI = d.var$RI, ET = d.var$ET, 
               P_c = s.var$P_c, P_z = s.var$P_z, P_s = s.var$P_s, 
               theta_init = d.map$theta_init, theta_sat = s.var$theta_sat, 
               theta_fc = s.var$theta_fc,
               SD = s.var$SD, K = s.var$K, P_I = s.var$P_I, n_s = s.map$n_s, 
               CC = d.var$CC, GC = s.var$GC, IMP = s.map$IMP, PH = d.var$PH, 
               D = s.var$D, NV = s.var$NV, d_a = s.var$d_a,
               DK_c = s.var$DK_c, DK_z = s.var$DK_z, DK_s = s.var$DK_s, 
               DR_c = s.var$DR_c, DR_z = s.var$DR_z, DR_s = s.var$DR_s, 
               Breaking = d.var$Breaking, Init_point = d.var$Init_point, 
               R_Type = 0, slpMode = 2, ALL = TRUE)
## Check results
## 1-a. Total volume of runoff generated from each element
plot(Result$Q_out)
## 1-b. Total volume of runoff per unit surface area generated from each element
plot(Result$Q_out/Result$A)

## 2-a. Total mass of eroded soil from each element
plot(Result$SL_out)
## 2-a. Total mass of eroded soil per unit surface area from each element
plot(Result$SL_out/Result$A)
}
}
\seealso{
    \code{\link{DMMF_Simple}} for one element during a day.
    \code{\link{Potato.Convex}} for the description of the data.
}
