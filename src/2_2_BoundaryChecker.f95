subroutine checkboundary(DEM, nr, nc, boundary)
    ! This is the code to identify boundaries of the DEM, created by Kwanghun
    ! Choi (2015).
    implicit none
    ! The "nr", and the "nc" are the number of row and column of the DEM.
    integer, intent( in ) :: nr, nc
    ! The "DEM" is the original DEM before sink filling.
    double precision, dimension( nr, nc ), intent( in ) :: DEM
    ! The "buffer" is the DEM with buffer to identify the boundary of the DEM.
    double precision, dimension(0:(nr+1),0:(nc+1)) :: buffer
    ! The "boundary" is the boundary of the DEM"
    double precision, dimension( nr, nc ) :: boundary
    ! The "boundary_t" is the temporary bondary for sink filling algorithm.
    logical, dimension( nr, nc ) :: boundary_t
    ! Dummy variables and the minimum and maximum row and column number.
    integer :: i, j
    ! Parameters for the NaN value.
    double precision, parameter ::&
        NaN = transfer(z'7ff8000000000000', 1.0d0)

    buffer = NaN; buffer(1:nr, 1:nc) = DEM
    boundary = NaN; boundary_t = .FALSE.
    forall( i = 1:nr, j = 1:nc )
        boundary_t(i,j) = any( isnan( buffer( i-1:i+1, j-1:j+1 ) ) ) .and. &
            ( any( .not. isnan( buffer( i-1:i+1, j-1:j+1 ) ) ) ) .and. &
            ( .not. isnan( buffer( i, j ) ) )
    end forall
    
    where( boundary_t ) boundary = 1.0d0

end subroutine checkboundary



