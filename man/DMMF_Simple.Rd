\name{DMMF_Simple}
\alias{DMMF_Simple}
\title{Simplified Daily based Morgan--Morgan--Finney (DMMF) soil erosion model for one element during a day}
\concept{soil erosion model}
\description{
    This is the simplified version of \code{\link{DMMF}} for simple one element (cell) during a day
}
\usage{
DMMF_Simple( W, L = W/cos(S), S, R, RI, ET, 
            P_c, P_z, P_s, theta_init, theta_sat, theta_fc, 
            SD, K, P_I, n_s, CC, GC, IMP, PH, D, NV, d_a = 0.005, 
            DK_c = 0.1, DK_z = 0.5, DK_s = 0.3, DR_c = 1.0, DR_z = 1.6, DR_s = 1.5, 
            Q_in = 0, IF_in = 0, SL_c_in = 0, SL_z_in = 0, SL_s_in = 0, R_type = 0)
}
\arguments{
    \item{W}{Width of an element (cell) (unit: m)}
    \item{L}{Length of an element (cell) (unit: m). If \code{L} is not specified, the function assumes \code{L} as \code{W}/cos(\code{S}) (i.e., The function assumes a cell as a subset of DEM with slope angle of \code{S}).}
    \item{S}{Slope angle of an element (cell) (unit: rad)}
    \item{R}{Daily rainfall per unit area (unit: mm)}
    \item{RI}{Mean hourly rainfall intensity during a day (unit: mm/h)}
    \item{ET}{Evapotranspiration per unit area (unit: mm)}
    \item{P_c}{Proportion of clay particles in surface soil (0-1)}
    \item{P_z}{Proportion of silt particles in surface soil (0-1)}
    \item{P_s}{Proportion of sand particles in surface soil (0-1)}
    \item{theta_init}{Initial soil water content of entire soil profile per unit area (unit: vol/vol)}
    \item{theta_sat}{Saturated soil water content of entire soil profile per unit area (unit: vol/vol)}
    \item{theta_fc}{Soil water content at field capacity of entire soil profile per unit area (unit: vol/vol)}
    \item{SD}{Soil depth of entire soil profile (unit: m)}
    \item{K}{Saturated lateral hydraulic conductivity of entire soil profile (unit: m/d)}
    \item{P_I}{Proportion of permanent interception area of rainfall (0-1)}
    \item{n_s}{Manning's roughness coefficient of the soil surface (unit: \eqn{\mathrm{s/m^{1/3}}}{s/m^(1/3)})}
    \item{CC}{Proportion of area with canopy cover (0-1)}
    \item{GC}{Proportion of area with pervious vegetated ground cover (0-1)}
    \item{IMP}{Proportion of area with impervious ground cover (0-1)}
    \item{PH}{Average height of vegetation or crop cover where leaf drainage start to fall (unit: m)}
    \item{D}{Average diameter of individual plant elements at the surface (unit: m)}
    \item{NV}{A number of individual plant elements per unit area (unit: \eqn{\mathrm{number/m^2}}{number/m^2})}
    \item{d_a}{Typical flow depth of surface runoff (unit: m)}
    \item{DK_c}{Detachability of clay particles by rainfall (unit: g/J)}
    \item{DK_z}{Detachability of silt particles by rainfall (unit: g/J)}
    \item{DK_s}{Detachability of sand particles by rainfall (unit: g/J)}
    \item{DR_c}{Detachability of clay particles by surface runoff (unit: g/mm)}
    \item{DR_z}{Detachability of silt particles by surface runoff (unit: g/mm)}
    \item{DR_s}{Detachability of sand particles by surface runoff (unit: g/mm)}
    \item{Q_in}{The volume of surface runoff from outside of an element (cell) (unit: L)}
    \item{IF_in}{The volume of subsurface interflow from outside of an element (cell) (unit: L)}
    \item{SL_c_in}{The mass of clay particles from outside of an element (cell) (unit: kg)}
    \item{SL_z_in}{The mass of silt particles from outside of an element (cell) (unit: kg)}
    \item{SL_s_in}{The mass of sand particles from outside of an element (cell) (unit: kg)}
    \item{R_type}{Integer object of each rainfall type for estimating kinetic energy of direct throughfall (DT) (0-8) (see details)}
}
\details{
    Detailed information about options of \code{R_type} can be found in details of \code{\link{DMMF}}.\cr
    This function is suitable for projecting the DMMF model for a field represented in one element during a day. 
}
\value{
    The output of the function \code{DMMF_Simple} is a data frame that contains following numeric elements:\cr
    \itemize{
        \item \code{Q_out}: Volume of surface runoff flowing from the element (unit: L)
        \item \code{IF_out}: Volume of subsurface water flowing from the element (unit: L)
        \item \code{theta_r}: Remaining soil water content of the element (unit: vol/vol) 
        \item \code{SL_c_out}: Mass of clay outputs from the element (unit: kg) 
        \item \code{SL_z_out}: Mass of silt outputs from the element (unit: kg)
        \item \code{SL_s_out}: Mass of sand outputs from the element (unit: kg)
        \item \code{A}: Surface area of the element (unit: \eqn{\mathrm{m^2}}{m^2})
    }
}
\references{
    Choi, K., Arnhold, S., Huwe, B., Reineking, B. (2017). 
    Daily based Morgan--Morgan--Finney (DMMF) model: A spatially distributed conceptual soil erosion model to simulate complex soil surface configurations. \emph{Water}, 9(4), 278. 

    Shin, S. S., Park, S. D., and Choi, B. K. (2016). 
    Universal power law for relationship between rainfall kinetic energy and rainfall intensity. \emph{Advances in Meteorology}, Article ID 2494681, 11 pages.

    Arnhold, S., Ruidisch, M., Bartsch, S., Shope, C., Huwe, B. (2013). 
    Simulation of runoff patterns and soil erosion on mountainous farmland with and without plastic-covered ridge-furrow cultivation in South Korea. \emph{Transactions of the ASABE}, 56(2):667--679.

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
## Load example data for test running DMMF_Simple function
data(Potato.Convex)
attach(Potato.Convex)
## Make toy dataset using summary statistics of field datasets of "Potato.Convex" data.
## The width and slope of the field are set according to the slope length of the convex field 
## from Arnhold et al. (2013)
W = 25
L = 25 
## The slope of the field is set as the average slope angle of the field from Arnhold et al. (2013)
S <- pi/180 * 9
## The proportion of impervious areas are estimated as the mean value of the IMP map.
IMP <- cellStats(s.map$IMP, stat='mean', na.rm=TRUE)
## "n_s" estimated using the guide value of RFR using paraplough (10cm/m) 
## and conversion equation from RFR to Manning's n from Morgan and Duzant (2008).
n_s <- 0.171  
## Using the ridge height of the field from Arnhold et al. (2013) 
## as the hydrological radius (flow depth) of the field.
d_a <- 0.15
## We use dynamic variables of the 12th day that has enough rainfall to produce soil erosion.
d.var_12<- d.var[12,]

## Run DMMF_Simple function
Output <- DMMF_Simple( W = W, L = L, S = S, R = d.var_12$R, RI = d.var_12$RI, ET = d.var_12$ET, 
                      P_c = s.var$P_c, P_z = s.var$P_z, P_s = s.var$P_s, 
                      theta_init = s.var$theta_fc, theta_sat = s.var$theta_sat, 
                      theta_fc = s.var$theta_fc, SD = s.var$SD, K = s.var$K, P_I = s.var$P_I, 
                      n_s = n_s, CC = d.var_12$CC, GC = s.var$GC, IMP = IMP, PH = d.var_12$PH, 
                      D = s.var$D, NV = s.var$NV, d_a = d_a, 
                      DK_c = s.var$DK_c, DK_z = s.var$DK_z, DK_s = s.var$DK_s, 
                      DR_c = s.var$DR_c, DR_z = s.var$DR_z, DR_s = s.var$DR_s, 
                      Q_in = 0, IF_in = 0, SL_c_in = 0, SL_z_in = 0, SL_s_in = 0, R_type = 0)

## Print outputs.
Output
## Calculate per unit surface area
## Runoff per surface area of the field
Output$Q_out / Output$A 
## Interflow per surface area of the field
Output$IF_out / Output$A 
## Soil loss per surface area of the field
(Output$SL_c_out + Output$SL_z_out + Output$SL_s_out) / Output$A 
}
}
\seealso{
    \code{\link{DMMF}} for fields and catchments with more than one element.
}
