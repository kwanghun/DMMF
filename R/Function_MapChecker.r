MapChecker <- 
    function ( DEM )
    {
        if ( missing( DEM ) || !is( DEM, "RasterLayer" ) )
        { 
            stop("Invalid type of DEM") 
        }else if ( xres( DEM ) == yres( DEM ) ) 
        {
            # Extract basic infomation of a DEM
            res     <- xres( DEM )
            xmn     <- extent(DEM)[1]
            xmx     <- extent(DEM)[2]
            ymn     <- extent(DEM)[3]
            ymx     <- extent(DEM)[4]
            crs     <- proj4string(DEM)
            # Convert format of a DEM from raster to matrix for FORTRAN subroutine.  
            DEM[is.na(DEM)] <- -999999 
            DEM_m   <- as.matrix( DEM )
        }else{
            stop("Resolution of DEM is not compatible")
        }

        # Run mapchecking algorithm from FORTRAN
        Outputs <- .Fortran(C_mapchecker, DEM = DEM_m, nr = nrow(DEM), nc = ncol(DEM), boundary = DEM_m, sink = DEM_m, stand = DEM_m, NAOK = TRUE)

        # Results
        # 1. call from the outputs
        Boundary_m <- Outputs$boundary
        Sink_m <- Outputs$sink 
        Stand_m <- Outputs$stand 
        # 2. rasterize the outputs
        Boundary <- raster(Boundary_m, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
        Sink <- raster(Sink_m, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
        Stand <- raster(Stand_m, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
        # 3. Make the list of the outputs
        Result <- list(Boundary, Sink, Stand)
        names(Result) <- c("boundary", "sink", "stand")
        # 4. return Result
        return(Result)
    }
