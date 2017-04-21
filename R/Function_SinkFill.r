SinkFill <- 
    function( DEM, Boundary, min_angle = 0.00001)
    {
        # Input variable check
        # 1. DEM check
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

        # 2. Boundary check
        if ( missing( Boundary ) )
        { 
            Boundary_list <- .Fortran(C_checkboundary, DEM = DEM_m, nr = nrow(DEM), nc=ncol(DEM), boundary = DEM_m, NAOK = T)
            Boundary_m <- Boundary_list$boundary 
            Boundary <- raster(Boundary_m, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
        }
        Boundary[is.na(Boundary)] <- -999999 
        Boundary_m   <- as.matrix( Boundary )

        # Run SinkFill module
        SinkFill_result <- .Fortran( C_sinkfill, DEM = DEM_m, nr = nrow( DEM ), 
                                    nc = ncol( DEM ), res = res, boundary = Boundary_m,
                                    min_angle = min_angle, DEM_nosink = DEM_m, partition = DEM_m,
                                    NAOK = T )

        # Output
        DEM_nosink <- raster(SinkFill_result$DEM_nosink, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
        partition <- raster(SinkFill_result$partition, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
        result <- list(DEM_nosink, partition)
        names(result) <- c("DEM_nosink", "partition")
        return( result )
    }



