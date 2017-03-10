SinkFill <- 
    function( DEM, Boundary, min_angle = 0.00001)
    {
        # Input variable check
        # 1. DEM check
        if ( missing( DEM ) || !is( DEM, "RasterLayer" ) )
        { 
            stop("Invalid type of DEM") 
        }else{
            if ( xres( DEM ) == yres( DEM ) )
            {
                DEM_m   <- as.matrix( DEM )
                res     <- xres( DEM )
            }else{
                stop("Resolution of DEM is not compatible")
            }
        }

        # 2. Boundary check
        if ( missing( Boundary ) )
        { 
            Boundary_list <- .Fortran(C_checkboundary, DEM = DEM_m, nr = nrow(DEM), nc=ncol(DEM), boundary = DEM_m*0, NAOK=T)
            Boundary_m <- Boundary_list$boundary 
        }else{
            Boundary_m <- as.matrix( Boundary ) 
        }
        
        # 3. Result map
        DEM_nosink_m <- 0 * DEM_m

        # Run SinkFill module
        SinkFill_result <- .Fortran( C_sinkfill, DEM = DEM_m, nr = nrow( DEM ), 
                                    nc = ncol( DEM ), res = res, boundary = Boundary_m,
                                    min_angle = min_angle, DEM_nosink = DEM_nosink_m,
                                    NAOK = T )
        
        # Output
        DEM_nosink <- DEM; DEM_nosink[] <- SinkFill_result$DEM_nosink
        return( DEM_nosink )
    }



