function slope ( m_block, res, order )
    implicit none
    ! The "m_block" represents matrix block for flow direction algorithm. 
    double precision, dimension( -1:1, -1:1 ), intent( in ) :: m_block
    double precision, intent( in ) :: res
    integer, intent( in ) :: order
    double precision :: dx=0.0d0, dy=0.0d0, slope
    integer :: i, dxcnt, dycnt
    if ( order == 2 ) then
    ! This is the 2nd order finite difference method (Fleming and Hoffer, 1979;
    ! Zhou and Liu, 2004), modified to calculate the slope even the cell at
    ! edges and surrounded by adjacent cells of NaN values.
        dx = 0.0d0; dy = 0.0d0 ! Initialize the dummy variables as null.
        if ( m_block( 1, 0 ) == m_block( 1, 0 ) ) then
        ! NaN value check for the cell below the center cell.
            if ( m_block( -1, 0 ) == m_block( -1, 0 ) ) then
            ! NaN value check for the cell above the center cell.
                dx = ( m_block( 1, 0 ) - m_block( -1, 0 ) ) * 0.5d0
            else
            ! If the above cell from the center is NaN, only use the center cell
            ! and the below cell from the center
                dx = m_block( 1, 0 ) - m_block( 0, 0 ) 
            end if
        else
            if ( m_block( -1, 0 ) == m_block( -1, 0 ) ) then
            ! If the below cell from the center is NaN, only use the center cell
            ! and the above cell from the center
                dx = ( m_block( 0, 0 ) - m_block( -1, 0 ) )
            end if
        end if
        if ( m_block( 0, 1 ) == m_block( 0, 1 ) ) then
        ! NaN value check for the right cell from the center cell.
            if ( m_block( 0, -1 ) == m_block( 0, -1 ) ) then
            ! NaN value check for the left cell from the center cell.
                dy = ( m_block( 0, 1 ) - m_block( 0, -1 ) ) * 0.5d0
            else
            ! If the left cell from the center is NaN, only use the center cell
            ! and the right cell from the center
                dy = m_block( 0, 1 ) - m_block( 0, 0 ) 
            end if
        else
            if ( m_block( 0, -1 ) == m_block( 0, -1 ) ) then
            ! If the right cell from the center is NaN, only use the center cell
            ! and the left cell from the center
                dy = ( m_block( 0, 0 ) - m_block( 0, -1 ) )
            end if
        end if
    else if ( order == 3 ) then
    ! This is the 3rd order finite difference method (Sharpnack and Akin, 1969;
    ! Zhou and Liu, 2004), modified to calculate the slope even the cell at
    ! edges and surrounded by adjacent cells of NaN values.
    ! The procedure is similar to 2nd order finite difference method which is
    ! described above, but it calculates row average and column average using
    ! upper to lower for column ( 3 ), right to left for row ( 3 ) 
        dx = 0.0d0; dy = 0.0d0; dxcnt = 3; dycnt = 3
        ! Initialize the dummy variables.
        do i = -1, 1
            if ( m_block( 1, i ) == m_block( 1, i ) ) then
                if ( m_block( -1, i ) == m_block( -1, i ) ) then
                    dx = dx + ( m_block( 1, i ) - m_block( -1, i ) ) * 0.5d0
                else
                    dx = dx + m_block( 1, i ) - m_block( 0, i ) 
                end if
            else
                if ( m_block( -1, i ) == m_block( -1, i ) ) then
                    dx = dx + ( m_block( 0, i ) - m_block( -1, i ) )
                else 
                    dxcnt = dxcnt - 1
                end if
            end if
            if ( m_block( i, 1 ) == m_block( i, 1 ) ) then
                if ( m_block( i, -1 ) == m_block( i, -1 ) ) then
                    dy = dy + ( m_block( i, 1 ) - m_block( i, -1 ) ) * 0.5d0
                else
                    dy = dy + m_block( i, 1 ) - m_block( i, 0 ) 
                end if
            else
                if ( m_block( i, -1 ) == m_block( i, -1 ) ) then
                    dy = dy + ( m_block( i, 0 ) - m_block( i, -1 ) )
                else
                    dycnt = dycnt - 1
                end if
            end if
        end do
        if ( dxcnt .ne. 0 ) then
            dx = dx / dxcnt
        else
        ! If there is no valud cell in the left and right side of the center,
        ! then set dx as null
            dx = 0.0d0
        end if
        if ( dycnt .ne. 0 ) then
            dy = dy / dycnt
        else
        ! If there is no valud cell in the up and down side of the center,
        ! then set dy as null
            dy = 0.0d0
        end if
    end if
    ! Slope can be calculated from the equation described below.
    slope = datan( dsqrt( dx * dx + dy * dy ) / res )  
    ! The case coded above is for quadratic grid cell of same x, y resolution,
    ! if we want to consider different x,y resolution, divide resolution of its
    ! direction to dx and dy.
end function slope
