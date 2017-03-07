subroutine sinkfill(DEM, nr, nc, res, boundary, min_angle, DEM_nosink)
    ! This is the Fortran code to produce sink free DEM using an algorithm
    ! proposed by Wang & Liu (2006), and code is created by Kwanghun Choi (2015). 
    ! References
    ! [1] Wang, L. and Liu, H. (2006). An efficient method for identifying and
    ! filling surface depressions in digital elevation models for hydrologic
    ! analysis and modelling. International Journal of Geographical Information
    ! Science, 20(2):193–213.
    ! [2] Volker Wichmann (2007) Module Fill Sinks (Wang & Liu). SAGA-GIS Module 
    ! Library Documentation (v2.1.3) [ cited 2015. 08. 15 ], Available from:
    ! http://www.saga-gis.org/saga_module_doc/2.1.3/ta_preprocessor_4.html
    implicit none
    ! The "nr", and the "nc" are the number of row and column of the DEM.
    integer, intent(in) :: nr, nc
    ! The "res" is the reolution of the DEM.
    ! The "min_angle" is the minimum angle between the cell and the adjacent
    ! cells.
    double precision, intent(in) :: res, min_angle
    ! The "DEM" is the original DEM before sink filling.
    ! The "DEM_t" is the temporary DEM for sink filling algorithm.
    ! The "boundary" is the boundary of the DEM"
    ! The "DEM_nosink" is the sink filled DEM.
    double precision, dimension( nr, nc ) :: DEM, DEM_t, boundary, DEM_nosink
    ! Dummy variables and the minimum and maximum row and column number.
    integer :: i, j, k , r, c, mnr, mxr, mnc, mxc
    ! Variables of locations for maximum value and minumum value.
    integer, dimension(2) :: min_loc
    ! Parameters for the NaN value.
    double precision, parameter ::&
        NaN = transfer(z'7ff8000000000000', 1.0d0)

    DEM_t = DEM

    do i = 1, size(DEM)
        min_loc = minloc( DEM_t,&
            ( ( boundary .eq. 1.0d0 ) .and. (.not. isnan( DEM_t ) ) ) )
        r = min_loc(1)
        c = min_loc(2)
        if ( isnan( DEM_t( r, c ) ) ) then
            exit
        else
            mnr = max0( 1, r - 1 )
            mxr = min0( nr, r + 1 )
            mnc = max0( 1, c - 1 )
            mxc = min0( nc, c + 1 )
            forall( j= mnr:mxr, k= mnc:mxc, .not.isnan( DEM_t( j, k ) ) ) 
                DEM_t( j, k ) = dmax1( DEM_t( j, k ), DEM_t( r, c ) +&
                    res * dsqrt( dble( (j - r)**2.0d0 + (k - c)**2.0d0 ) ) *&
                    dtan( min_angle ) )
                boundary( j, k ) = 1.0d0
            end forall
            DEM_nosink(r, c) = DEM_t(r, c)
            ! order(r, c) = dble(i)
            DEM_t(r, c) = NaN
        end if
    end do
end subroutine sinkfill

