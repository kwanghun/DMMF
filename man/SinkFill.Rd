\name{SinkFill}
\alias{SinkFill}
\title{Sinkfill algorithm based on Wang and Liu (2006)}
\description{
    Implementation of sinkfill algorithm based on Wang and Liu (2006).
}
\usage{
    SinkFill( DEM, Boundary, min_angle = 0.00001)
}
\arguments{
    \item{DEM}{RasterLayer object of digital elevation model (DEM) which needs sinkfill for DMMF model}
    \item{Boundary}{RasterLayer object for setting the locations of local sinks in the DEM, if \code{Boundary} is missing, every boundary cells of the map are set as potential sinks of the map.}
    \item{min_angle}{Numeric object of angle to make slope between cells (unit: rad)}
}
\details{
    The function is the implementation of sink-fill algorithm from Wang and Liu (2006), which is already implemented as a SAGA-GIS module (see Wichmann (2007)).\cr
  The function \code{SinkFill} is modified from the original algorithm with including the option to set the initial boundary manually (e.g., global sink of the field).  
}
\value{
    If \code{Boundary} is missing, function automatically identifies internal and external boundary of the input DEM which is in contact with cells having NaN value. 
    The output of the function \code{SinkFill} is a list of RasterLayer objects of the sinkfilled DEM (\code{$nosink}) and partitions (\code{$partition}) of the area for the \code{\link{DMMF}} model.
}
\references{
    Volker Wichmann (2007). 
    Module Fill Sinks (Wang & Liu). SAGA-GIS Module Library Documentation (v2.1.3) [ cited 2015. 08. 15 ], Available from: \url{http://www.saga-gis.org/saga_tool_doc/2.1.3/ta_preprocessor_4.html}

    Wang, L. and Liu, H. (2006). 
    An efficient method for identifying and filling surface depressions in digital elevation models for hydrologic analysis and modelling. \emph{International Journal of Geographical Information Science}, 20(2):193--213.

}
\author{
    Kwanghun Choi
}
\examples{
\dontrun{
## Load example data for test running SinkFill function
data(Potato.Convex)
attach(Potato.Convex)
## Run SinkFill function with original DEM which has internal sinks.
DEM_with_sink <- s.map$DEM_original
DEM_sinkfilled <- SinkFill( DEM_with_sink )
DEM_nosink <- DEM_sinkfilled$nosink
partitions <- DEM_sinkfilled$partition
## Check maps
par(mfrow=c(2,2))
plot(DEM_with_sink)
plot(DEM_nosink)
plot(DEM_with_sink - DEM_nosink)
plot(partitions)
}
}
\seealso{
    \code{\link{Potato.Convex}} for the description of the data.
}
