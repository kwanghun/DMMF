function MDInf ( m_block )
    ! The purpose of this function is to calculate flow direction using the
    ! triangular multiple flow direction algorithm by Seibert and McGlynn
    ! (2007) ( MDInf )
    implicit none
    ! The "m_block" represents matrix block for flow direction algorithm. 
    double precision, dimension( -1:1 , -1:1 ), intent( in ) :: m_block
    ! The "m_rel_h" is the relative height matrix to the center. 
    ! The "MDInf" is the result of the MDInf flow direction algorithm.
    double precision, dimension( -1:1, -1:1 ) :: m_rel_h, MDInf
    ! The "m_index" is for the indexing.
    integer, dimension( -1:1, -1:1 ) :: m_index
    ! V1 and V2 are the coordinate of the adjacent cell from the center.
    ! N is the normal vector of the plane which V1 and V2 belong to.
    double precision, dimension (3) :: V1, V2, N
    integer :: i, x1, y1, x2, y2
    double precision :: a1, a2, a3, slp
    double precision, parameter :: pi = 3.141592653589793239d0
    ! This is the interface for the cross product function. 
    interface
        function crossproduct ( a, b )
            double precision, dimension( 3 ), intent( in ) :: a, b
            double precision, dimension( 3 ) :: crossproduct
        end function crossproduct
    end interface
    ! Initialize the MDInf matrix as 0.0 and m_index as 0.
    MDInf = 0.0d0
    m_index = 0
    ! Calculate m_rel_h which is the relative height to the center.
    m_rel_h = m_block - m_block( 0, 0 )
    ! Because the range of atan2 is (-pi, pi], we calculate pair of triangular
    ! facet from -pi.
    do i = -4, 3, 1
        ! The angle of the first element (in the polar coordinate)
        a1 = dble(i) * pi / 4.0d0
        ! The angle of the second element (in the polar coordinate)
        a2 = dble(i + 1) * pi / 4.0d0
        ! The x, y, and z values of the first element (in the cartesian)
        x1 = nint( dcos( a1 ) ); y1 = nint( dsin( a1 ) )
        V1( 1 ) = dble(x1); V1( 2 ) = dble(y1); V1( 3 ) = m_rel_h( x1, y1 )
        ! The x, y, and z values of the second element (in the cartesian)
        x2 = nint( dcos( a2 ) ); y2 = nint( dsin( a2 ) )
        V2( 1 ) = dble(x2); V2( 2 ) = dble(y2); V2( 3 ) = m_rel_h( x2, y2 )
        ! The case that both values are not NaN
        if ( V1(3) == V1(3) .and. V2(3) == V2(3) ) then 
            ! The case that at least one of the values is less than 0 
            if ( ( V1(3) .lt. 0.0d0 ) .or. ( V2(3) .lt. 0.0d0 ) ) then   
                ! N is the cross product of a pair of two vectors.
                N = crossproduct( V1, V2 )
                ! a3 is the angle of N projeted on the x-y plane.
                a3 = atan2 ( N(2), N(1) )
                ! If a3 is in the range of [a1, a2],
                if ( ( a3 .ge. a1 ) .and. ( a3 .le. a2 ) ) then
                    ! slope is calculated from Seibert and McGlynn (2007)
                    slp = dabs( dsqrt( N(1)*N(1) + N(2)*N(2) ) / N(3) )
                    ! the ratio of inflow to the adjacent cell.
                    ! = slope * ( (pi/4) - (a3-a1) ) / (pi / 4)
                    ! = 4 * slope * (a2 - a3) / pi # <a2-a1 = pi/4>
                    MDInf( x1, y1 ) = MDInf( x1, y1 ) + &
                        4.0d0 * slp * ( a2 - a3 ) / pi  
                    MDInf( x2, y2 ) = MDInf( x2, y2 ) + &
                        4.0d0 * slp * ( a3 - a1 ) / pi 
                ! If a3 is not in the range of [a1, a2],
                else
                    ! 1. Decide which of the element has smaller value
                    ! 2. Add 1 to the index matrix of the element.
                    ! 3. If index reaches 2, calculate slope of the element
                    ! 4. and calculate the ratio of the flow.
                    if ( V1(3) .lt. V2(3) ) then
                        m_index( x1, y1 ) = m_index( x1, y1 ) + 1
                        if ( m_index ( x1, y1 ) == 2 ) then
                            slp = dabs( dsqrt( V1(1)*V1(1) + V1(2)*V1(2) ) / V1(3) )
                            MDInf( x1, y1 ) = MDInf( x1, y1 ) + slp 
                        end if
                    else
                        m_index( x2, y2 ) = m_index( x2, y2 ) + 1
                        if ( m_index ( x2, y2 ) == 2 ) then
                            slp = dabs( dsqrt( V2(1)*V2(1) + V2(2)*V2(2) ) / V2(3) )
                            MDInf( x2, y2 ) = MDInf( x2, y2 ) + slp 
                        end if
                    end if
                end if
            else
                continue
            end if
        else 
            ! 1. Decide which of the element has value except NaN
            ! 2. Add 1 to the index matrix of the element.
            ! 3. If index reaches 2, calculate slope of the element
            ! 4. and calculate the ratio of the flow.
            if ( ( V1(3) == V1(3) ) .and. ( V1(3) .lt. 0.0d0 ) ) then
                m_index( x1, y1 ) = m_index( x1, y1 ) + 1
                if ( m_index ( x1, y1 ) == 2 ) then
                    slp = dabs( dsqrt( V1(1)*V1(1) + V1(2)*V1(2) ) / V1(3) )
                    MDInf( x1, y1 ) = MDInf( x1, y1 ) + slp 
                end if
            else if ( ( V2(3) == V2(3) ) .and. ( V2(3) .lt. 0.0d0 ) ) then
                m_index( x2, y2 ) = m_index( x2, y2 ) + 1
                if ( m_index ( x2, y2 ) == 2 ) then
                    slp = dabs( dsqrt( V2(1)*V2(1) + V2(2)*V2(2) ) / V2(3) )
                    MDInf( x2, y2 ) = MDInf( x2, y2 ) + slp 
                end if
            end if
        end if
    end do
    ! If there is at least one flow from the center ( it means that sum of the
    ! matrix is greater than 0), divide the matrix by the sum of it for
    ! normalization, if not just return the matrix of value 0.
    if ( sum( MDInf ) .gt. 0.0d0 ) then
        MDInf = MDInf / sum( MDInf )
    else
        MDInf = 0.0d0
    end if
end function MDInf
function crossproduct ( a, b )
    ! The purpose of this function is to calculate crossproduct of two vectors.
    implicit none
    double precision, dimension( 3 ), intent( in ) :: a, b
    double precision, dimension( 3 ) :: crossproduct
    crossproduct( 1 ) = a( 2 ) * b( 3 ) - a( 3 ) * b( 2 )
    crossproduct( 2 ) = a( 3 ) * b( 1 ) - a( 1 ) * b( 3 )
    crossproduct( 3 ) = a( 1 ) * b( 2 ) - a( 2 ) * b( 1 )
end function crossproduct
