\name{MapChecker}
\alias{MapChecker}
\title{Map checker for raster layers}
\description{
    Simple boundary, sink and stand checker for DEMs (Digital Elevation Models).
}
\usage{
    MapChecker( DEM )
}
\arguments{
    \item{DEM}{RasterLayer object of digital elevation model (DEM) which needs sinkfill for DMMF model}
}
\details{
    Simple boundary, sink and stand checker for DEMs (Digital Elevation Models).
    The code defines a cell as a
    \itemize{
        \item \code{boundary}: when there are at least one NA values and at least one cell with valid value in adjacent cells.
        \item \code{sink}: when all adjacent cells of it are higher than the target cell. 
        \item \code{stand}: when all adjacent cells of it are higher than or equal to the target cell. 
    }
}
\value{
    The output of the function \code{MapChecker} is a RasterLayer object of the \code{boundary}, \code{sink}, and \code{stand}.
}
\author{
    Kwanghun Choi
}
\examples{
\dontrun{
## Load example data for test running B_Checker function
data(Potato.Convex)
attach(Potato.Convex)
## Run B_Checker function with original DEM which has internal sinks.
DEM <- s.map$DEM_original
MapCheck <- MapChecker( DEM )
## Check maps
par(mfrow=c(2,2))
plot(DEM)
plot(MapCheck$boundary)
plot(MapCheck$sink)
plot(MapCheck$stand)
}
}
\seealso{
    \code{\link{Potato.Convex}} for the description of the data.
}
