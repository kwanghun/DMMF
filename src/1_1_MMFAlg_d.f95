subroutine DMMF( DEM, nr, nc, res, option, days, R, RI, R_Type, ET,&
        P_c, P_z, P_s, theta_init, theta_sat, theta_fc, SD, K,& 
        P_I, n_s, d_a, CC, GC, IMP, PH, D, NV,&
        DK_c, DK_z, DK_s, DR_c, DR_z, DR_s,&
        Breaking, N_out, vc, Init_point,&
        A, Rf_r, SW_c_r, theta_r_r, TC_r,&
        Q_in_r, Q_out_r, IF_in_r, IF_out_r,& 
        SS_c_r, SS_z_r, SS_s_r, G_c_r, G_z_r, G_s_r,&
        SL_c_in_r, SL_z_in_r, SL_s_in_r, SL_in_r,&
        SL_c_out_r, SL_z_out_r, SL_s_out_r, SL_out_r )
! This is the program to run corrected Morgan-Morgan-Finney model
! suggested by Choi et al. (2015)
    implicit none
    ! Interface of the functions.
    interface
        function KE( RI, Climate, Rf, CC, PH )
            implicit none
            double precision, intent( in ) :: RI, Rf, CC, PH
            integer, intent ( in ) :: Climate
            double precision :: KE
        end function KE
        function MDInf ( m_block )
            implicit none
            double precision, dimension( -1:1 , -1:1 ), intent( in ) :: m_block
            double precision, dimension( -1:1,  -1:1 ) :: MDInf
        end function MDInf
        function slope ( m_block, res, order )
            implicit none
            double precision, dimension( -1:1, -1:1 ), intent( in ) :: m_block
            double precision, intent( in ) :: res
            integer, intent( in ) :: order
            double precision :: slope
        end function slope
    end interface
    ! These are the input variables.
    ! nr: the number of rows of the input rasters
    ! nc: the number of columns of the input rasters    
    ! option: option for the slope algorithm ( see slope function above )
    ! R_Type: the type of a climate region ( integer )
    integer, intent( in ) :: nr, nc, option, days, R_Type, N_out, vc
    ! res: the resolution of the input rasters.
    double precision, intent( in ) :: res
    ! DEM: Digital elevation model ( raster )
    ! R: the annual rainfall per unit square m ( mm/m^2 ).
    ! days: the number of rain days of the year.
    ! RI: typical intensity of erosive rain ( mm/h )
    ! P_x: proportion of the clay(c), the silt(z), and the sand(s) in the soil.
    ! theta_sat: Volumetric soil water content when soil is saturated (v/v) 
    ! theta_fc: Volumetric soil water content at field capacity ( v/v )
    ! SD: Depth of the soil ( m )
    ! K: the saturated lateral permeability of the soil ( m/day )
    ! ST: Proportion of rock fragments on the soil surface
    ! P_I: Permanent interception expressed as the proportion of rainfall
    ! RFR: Surface roughness ( cm/m )
    ! n_s: manning coefficient for the surface conditions assuming that there is
    ! no vegetations.
    ! R_Depth: flow or rill depth ( m )
    ! ET: daily mean evapotranspiration of the element in the period.
    ! E_t_E_0: Ratio of actual to potential evapotranspiration
    ! CC: Canopy cover of the soil surface protected by the vegetation canopy
    ! GC: Ground cover of the soil surface protected by vegetation cover on the ground
    ! PH: Plant height ( m )    
    ! D: Average diameter (m) of the individual plant elements (stems, leaves)
    ! at the ground surface
    ! NV: Number of plant elements per unit area ( number/m^2 ) at the ground
    ! surface
    ! d_a: rill depth for the real surface condition
    integer, dimension( days ) :: Breaking, Init_point
    double precision, dimension( nr, nc ) :: DEM
    double precision, dimension( nr, nc, days ) :: R, RI, ET
    double precision, dimension( nr, nc ) :: P_c, P_z, P_s,& 
        theta_sat, theta_fc, SD, K
    double precision, dimension( nr, nc, days ) :: P_I, n_s, d_a
    double precision, dimension( nr, nc, days ) :: CC, GC, IMP, PH, D, NV 
    double precision, dimension( nr, nc ) :: DK_c, DK_z, DK_s, DR_c, DR_z, DR_s
    ! These are the outputs
    ! Q: the volume of runoff generated on the element per unit surface area ( mm )
    ! Q_in: the total volume of runoff flows in to the element ( mm * m^2 )
    ! Q_in_u: Runoff influx from upper adjacent elements per unit area (mm)
    ! Q_out: The volume of runoff flows out from the element ( mm * m^2 )
    ! SW_c: the infiltrable soil water storage capacity of the element (mm)
    ! IF_in_u: the volume of interflow inputs on unit surface area ( mm )
    ! IF_in: the total volume of interflow flow in to the element ( mm * m^2 )
    ! IF_out: The volume of interflow flows out from the element ( mm * m^2 )
    ! TC: the transport capacity of the runoff ( kg/m^2 )
    ! SL: the weight of soil loss from the unit surface area ( kg/m^2 )
    ! SLt: the total weight of soil loss from the element ( kg )
    ! x_in: the amount of the sediment inflow from higher adjacent elements as a
    ! type of clay( c ), silt( z ), and sand( s ) particles. ( kg/m^2 )
    ! SS_x: Suspended sediments delivered to runoff ( kg/m^2 )
    double precision :: Q_in_u, Q, IF_in_u, IF_out_u,& 
        SL_c, SL_z, SL_s, SL_c_in_u, SL_z_in_u, SL_s_in_u
    double precision, dimension( nr, nc ) :: SW_c, TC
    double precision, dimension( 0:( nr + 1 ), 0:( nc + 1 ) ) :: &
        Q_in, IF_in, SL_c_in, SL_z_in, SL_s_in
    double precision, dimension( nr, nc ) :: Q_out, IF_out, SS_c, SS_z, SS_s,&
        G_c, G_z, G_s, SL_c_out, SL_z_out, SL_s_out, SWr_d
    ! S: slope of the element.
    ! L, W, A: length( L ), width( W ) and the area size( A ) of the element.
    ! Rf: effective rainfall of the element.
    ! R_sum: The total daily input of the surface water (mm) 
    ! R_0: the amount of rainfall per rain day.
    ! SW_sat: Maximum soil moisture storage capacity of the element (mm)
    ! SW_fc: the quantity of soil water at field capacity in the element (mm)
    ! SW_in: Soil water existing or flowing in prior to runoff process (mm)
    ! SW_ex: Soil water exceed maximum soil water storage capacity of the element (mm)
    ! SW: Average daily reserved water in the vadose zone (mm)
    double precision, dimension( nr, nc ) :: R_0, RI_0, Rf, ET_0, SW_init
    double precision, dimension( nr, nc ) :: A
    double precision :: S, W, R_sum,&
        SW_sat, SW_fc, SW_in, SW, SW_t
    !, SW_ex
    ! F_t: Detachment capacity of soil particles by raindrop impact ( kg/m^2 )
    ! F_x: Detachment of clay( c ), silt( z ) and sand( s ) particles by raindrop 
    ! impact ( kg/m^2 )
    ! DK_x: Detachbility of clay( c ), silt( z ) and sand( s ) particles by
    ! raindrop impact ( g/J )
    ! H_t: Detachment capacity of soil particles by surface runoff ( kg/m^2 )
    ! H_x: Detachment of clay( c ), silt( z ) and sand( s ) particles by surface
    ! runoff ( kg/m^2 )
    ! DR_x: Detachbility of clay( c ), silt( z ) and sand( s ) particles by
    ! surface runoff ( g/mm )
    !double precision, parameter :: DK_c = 0.1d0, DK_z = 0.5d0, DK_s = 0.3d0
    !double precision, parameter :: DR_c = 1.0d0, DR_z = 1.6d0, DR_s = 1.5d0
    double precision :: F_t, F_c, F_z, F_s
    double precision :: H_t, H_c, H_z, H_s
    ! G_x:  the delivery of detached particles to runoff ( kg/m^2 )
    ! G_x1: the delivery of detached particles to runoff after deposition 
    ! TC_x: the transport capacity of surface runoff ( kg/m^2 )
    ! DEP_x: the percentage of detached particles that is deposited ( % )       
    ! c: clay, z: silt, s: sand
    double precision :: TC_c, TC_z, TC_s, DEP_c, DEP_z, DEP_s
    
    ! Factors for surface runoff velocity and relevant factor.
    ! v_b: Reference velocity of the surface runoff of the smooth bare soil 
    ! v: generalized velocity considering all surface conditions by Petryk's eq.
    ! C: surface factor silmilar to RMMF-C factor.
    ! d_b: rill depth for the reference condition    
    ! n_t: Modified n' for Petryk's equation considering vegetation effects.
    double precision :: v_b, n_t, v, C
    double precision, parameter :: d_b = 0.005d0, n_b = 0.015d0 
    ! Factors for particle fall number and fall velocities.  
    ! Nf_x: Particle fall number of the clay(c), silt(z) and sand(s)
    ! vs_x: Fall velociy of clay(c), silt(z) and sand(s)
    double precision :: Nf_c, Nf_z, Nf_s, vs_c, vs_z, vs_s
    ! EPA: The ratiao of Erosion Protected Area
    double precision :: EPA
    double precision, parameter :: rho_s = 2650.0d0, rho = 1000.0d0,&
        eta = 0.0015d0, di_c = 0.000002d0, di_z = 0.00006d0, di_s = 0.0002d0 
    ! These are the dummy variables for calculation.
    ! mask: logical matrix for maxloc function to apply algorithm from highest
    ! altitude to the lowest.
    logical, dimension( nr, nc ) :: mask
    ! b_DEM: DEM matrix with buffers.
    double precision, dimension( 0:( nr + 1 ), 0:( nc + 1 ) ) :: b_DEM
    ! m_block: sliced matrix of 3 X 3 from DEM matrix.
    ! m_weight: result of flow direction algorithm using m_block (3 X 3 matrix).
    double precision, dimension( 3, 3 ) :: m_block, m_weight
    ! max_loc: the current center location of DEM matrix.
    integer, dimension(2) :: max_loc
    ! row: current calculation row number.
    ! col: current calculation col number.
    integer :: row, col
    ! i: dummy variable for the loop.
    integer :: i, j, t, init_counter
    ! Parameters
    double precision, parameter :: NaN = transfer(z'7ff8000000000000', 1.0d0)
    double precision, parameter :: pi = 3.141592653589793239d0
    double precision, parameter :: g = 9.80665d0
    ! Input of initial soil water content of whole soil profile
    double precision, dimension( nr, nc, N_out ) :: theta_init
    ! declaration and initialization of result matrix
    double precision, dimension( nr, nc, N_out ) :: Rf_r,& 
        SW_c_r, theta_r_r, TC_r,&
        Q_in_r, Q_out_r, IF_in_r, IF_out_r,& 
        SS_c_r, SS_z_r, SS_s_r, G_c_r, G_z_r, G_s_r,&
        SL_c_in_r, SL_z_in_r, SL_s_in_r, SL_in_r,&
        SL_c_out_r, SL_z_out_r, SL_s_out_r, SL_out_r
    integer, dimension( vc ) :: row_s, col_s
    double precision, dimension( vc ) :: slp, L
    double precision, dimension( 3, 3, vc ) :: weight_m

mask = .true.
! Initialize b_DEM using DEM and make buffers as NaN.
b_DEM = NaN; b_DEM( 1:nr, 1:nc ) = DEM
! Initialize the area and the width of elements (cells)
! Area is varying according to the slope of an elements.
! Width is fixed as a resolution of DEM.
A = NaN; W = res
do i = 1, vc
    max_loc = maxloc( DEM, mask = mask )
    row = max_loc( 1 ); col = max_loc( 2 )
    row_s( i ) = row; col_s( i ) = col
    mask( row, col ) = .false.
    ! m_block is the 3 X 3 sliced matrix around the element. 
    m_block = b_DEM( ( row - 1 ):( row + 1 ), ( col - 1 ):( col + 1 ) )
    ! m_weight is the 3 X 3 weight matrix for runoff and sediment
    ! distribution from MDInf flow direction algorithm.
    weight_m( 1:3, 1:3, i ) = MDInf( m_block )
    ! S is the slope of the element.
    slp( i ) = slope( m_block, res, option )
    L( i ) = W / dcos( slp( i ) )
    A( row, col ) = W * L( i )
end do

! Initialize the result matrix
SW_c_r = 0.0d0; theta_r_r = 0.0d0; TC_r = 0.0d0
Q_in_r = 0.0d0; Q_out_r = 0.0d0
IF_in_r = 0.0d0; IF_out_r = 0.0d0
SS_c_r = 0.0d0; SS_z_r = 0.0d0; SS_s_r = 0.0d0
G_c_r = 0.0d0; G_z_r = 0.0d0; G_s_r = 0.0d0
SL_c_in_r = 0.0d0; SL_z_in_r = 0.0d0; SL_s_in_r = 0.0d0; SL_in_r = 0.0d0
SL_c_out_r = 0.0d0; SL_z_out_r = 0.0d0; SL_s_out_r = 0.0d0; SL_out_r = 0.0d0

t = 1
init_counter = 1
SW_init = 1000.0d0 * theta_init( 1:nr, 1:nc, init_counter ) * SD

do j = 1, days
    R_0 = R( 1:nr, 1:nc, j )
    RI_0 = RI( 1:nr, 1:nc, j )
    ET_0 = ET( 1:nr, 1:nc, j )

    ! Matrix initialization
    Rf = NaN
    SW_c = NaN
    TC = NaN
    SS_c = NaN; SS_z = NaN; SS_s = NaN
    G_c = NaN; G_z = NaN; G_s = NaN
    Q_in = 0.0d0; IF_in = 0.0d0
    Q_out = NaN; IF_out = NaN 
    SL_c_in = 0.0d0; SL_z_in = 0.0d0; SL_s_in = 0.0d0
    SL_c_out = NaN; SL_z_out = NaN; SL_s_out = NaN
! all per unit quantity of Q and SL is the volume/mass of quantity per unit area.
! Therefore declare them as numeric and converting them using A(row, col).
! Q = NaN; Q_in_u = NaN; IF_in_u = NaN; IF_out_u = NaN 
! SL_c = NaN; SL_z = NaN; SL_s = NaN 
! SL_c_in_u = NaN; SL_z_in_u = NaN; SL_s_in_u = NaN

    ! Do loop from the hightest to the lowest elements.
    do i = 1, vc
        ! Declare row and col number
        row = row_s( i ); col = col_s( i )
        ! Call weight matrix of row, col
        m_weight = weight_m( 1:3, 1:3, i )
        ! S is the slope of the element.
        S = slp( i )
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Begin water quantity calculation module
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! R_0 is the rainfall of rain day.
            ! Rf is modified from Original version ( / cos( S ) -> * cos( S ) ).
            Rf( row, col ) = R_0( row, col ) * ( 1.0d0 - P_I( row, col, j ) ) * dcos( S )
            ! the amount of the runoff and the interflow on unit surface area 
            ! from other elements.
            Q_in_u = Q_in( row, col ) / A( row, col )
            IF_in_u = IF_in( row, col ) / A( row, col )
            ! ET is the average daily evapotranspiration from MOD16ET
            ! Modified
            ! To modify model for shorter term, it is better to use MOD16ET.
            ! SW_c: the infiltrable soil water storage capacity of the element ( mm )
            SW_sat = 1000.0d0 * theta_sat( row, col ) * SD( row, col )
            SW_fc = 1000.0d0 * theta_fc( row, col ) * SD( row, col ) 
            SW_in = SW_init( row, col ) + IF_in_u
            !SW_in =  dmax1( SW_sat - SW_fc - IF_in_u, 0.0d0 )
            !SW_ex = dmax1( 0.0d0, SW_in - SW_sat )
            SW_c( row, col ) = ( 1.0d0 - IMP( row, col, j ) ) *& 
                ( SW_sat - SW_in )
            ! Runoff per unit area from the element 
            R_sum = Rf( row, col ) + Q_in_u
            ! Modified as erasing empirical term related to the slope
            ! length.
            ! Also applied impermeable area as reduced soil moisture storage
            ! capacity in terms of area.
            if ( R_sum .gt. SW_c( row, col ) ) then
                Q = R_sum - SW_c( row, col )
            else
                Q = 0.0d0
            end if
            Q_out( row, col ) = Q * A( row, col )
            ! SW: overall soil water of whole soil profile ( mm )
            ! The overall amount of interflow from the element 
            SW = dmax1( 0.0d0, SW_in + R_sum -& 
                Q - ET_0( row, col ) )
            ! transferrable soil water
            SW_t = dmax1( 0.0d0, SW - SW_fc )

            if ( SW_t .gt. 0.0d0 ) then
                IF_out( row, col ) = dmin1( K( row, col ) * dsin( S ) *& 
                    SW_t * W, SW_t * A( row, col ) )
            else
                IF_out( row, col ) = 0.0d0
            end if
            ! Interflow from an element per unit surface area ( mm )
            IF_out_u = IF_out( row, col ) / A( row, col )
            ! Total amount of captured water in the soil ( mm )
            SWr_d( row, col ) = dmax1( 0.0d0, SW - IF_out_u )
            ! Distribute the water quantities from the element to lower adjacent
            ! elements.
            Q_in( ( row - 1 ) : ( row + 1 ), ( col - 1 ): ( col + 1 ) ) = &
                Q_in( ( row - 1 ) : ( row + 1 ), ( col - 1 ): ( col + 1 ) ) + &
                Q_out( row, col ) * m_weight
            IF_in( ( row - 1 ) : ( row + 1 ), ( col - 1 ): ( col + 1 ) ) = &
                IF_in( ( row - 1 ) : ( row + 1 ), ( col - 1 ): ( col + 1 ) ) + &
                IF_out( row, col ) * m_weight

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! End water quantity calculation module
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Begin sediment loss calculation module
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! The ratiao of Erosion Protected Area ( EPA )
            EPA = IMP( row, col, j ) + ( 1 - IMP( row, col, j ) ) * GC( row, col, j )
            ! Detachment of soil particles by raindrop impact ( kg/m^2 )
            F_t = 0.001d0 * dmax1( 0.0d0, ( 1.0d0 - EPA ) ) *& 
                KE( RI_0( row, col ), R_Type, Rf( row, col ),&
                CC( row, col, j ), PH( row, col, j ) )
            F_c = DK_c( row, col ) * P_c( row, col ) * F_t
            F_z = DK_z( row, col ) * P_z( row, col ) * F_t
            F_s = DK_s( row, col ) * P_s( row, col ) * F_t

            ! Detachment of soil particles by surface runoff ( kg/m^2 )
            H_t = 0.001d0 * ( Q**1.5d0 ) * ( dsin( S )**0.3d0 ) *&
                dmax1( 0.0d0, ( 1.0d0 - EPA ) )  
            H_c  = DR_c( row, col ) * P_c( row, col ) * H_t
            H_z  = DR_z( row, col ) * P_z( row, col ) * H_t
            H_s  = DR_s( row, col ) * P_s( row, col ) * H_t

            ! Influx of detached soil particles from adjacent elements ( kg/m^2 ).
            SL_c_in_u = SL_c_in( row, col ) / A( row, col )
            SL_z_in_u = SL_z_in( row, col ) / A( row, col )
            SL_s_in_u = SL_s_in( row, col ) / A( row, col )

            ! Suspended sediments which are delivered to runoff (SS) ( kg/m^2 )
            SS_c( row, col ) = F_c + H_c + SL_c_in_u
            SS_z( row, col ) = F_z + H_z + SL_z_in_u
            SS_s( row, col ) = F_s + H_s + SL_s_in_u

            ! Surface runoff velocities ( m/s )
            ! Reference velocities on the standard bare soil. (Reference)
            v_b = ( d_b ** 0.67d0 ) * dsqrt( dtan( S ) ) / n_b
            ! n' for flow velocity (n_t), n_s is the Manning's coefficients for
            ! surface condition.
            n_t = dsqrt( ( n_s( row, col, j ) * n_s( row, col, j ) ) +&
                ( 0.5d0 * D( row, col, j ) * NV( row, col, j ) *&
                ( d_a( row, col, j ) ** ( 4.0d0 / 3.0d0 ) ) ) / g )  
            ! Flow velocity
            v = ( ( d_a( row, col, j ) ** ( 2.0d0 / 3.0d0 ) ) *&
                dsqrt( dtan( S ) ) ) / n_t
            
            ! Particle fall velocities (vs_x) and particle fall number (Nf_x)
            if ( v .gt. 0.0d0 ) then 
                vs_c = di_c * di_c * ( rho_s - rho ) * g / ( 18.0d0 * eta )
                Nf_c = vs_c * L( i ) / ( v * d_a( row, col, j ) )
                vs_z = di_z * di_z * ( rho_s - rho ) * g / ( 18.0d0 * eta )
                Nf_z = vs_z * L( i ) / ( v * d_a( row, col, j ) )
                vs_s = di_s * di_s * ( rho_s - rho ) * g / ( 18.0d0 * eta )
                Nf_s = vs_s * L( i ) / ( v * d_a( row, col, j ) )

                ! Percentage of detached particles that is deposited ( DEP_x ) (%)
                DEP_c = dmin1( ( 0.441d0 * ( Nf_c**0.29d0 ) ), 1.0d0 )
                DEP_z = dmin1( ( 0.441d0 * ( Nf_z**0.29d0 ) ), 1.0d0 )
                DEP_s = dmin1( ( 0.441d0 * ( Nf_s**0.29d0 ) ), 1.0d0 )
            else
                DEP_c = 1.0d0; DEP_z = 1.0d0; DEP_s = 1.0d0
            end if

            ! Delivery of detached particles to runoff ( kg/m^2 )
            G_c( row, col ) = SS_c( row, col ) * ( 1.0d0 - DEP_c ) 
            G_z( row, col ) = SS_z( row, col ) * ( 1.0d0 - DEP_z ) 
            G_s( row, col ) = SS_s( row, col ) * ( 1.0d0 - DEP_s ) 

            ! RMMF-C factor with velocities 
            ! the how the surface runoff velocity is faster or slower than that
            ! of standard condition.
            if ( v_b == 0.0d0 ) then
                C = 0.0d0
            else
                C = ( v / v_b )  
            end if

            ! Transport capacity of the runoff ( kg/m^2 )
            TC( row, col ) = 0.001d0 * C * Q * Q *&
                dsin( S )
            TC_c = P_c( row, col ) * TC( row, col )
            TC_z = P_z( row, col ) * TC( row, col )
            TC_s = P_s( row, col ) * TC( row, col )

            ! Sediment balance
            if ( TC_c .ge. G_c( row, col ) ) then
                SL_c = G_c( row, col )
            else
                SL_c = TC_c
            end if
            if ( TC_z .ge. G_z( row, col ) ) then
                SL_z = G_z( row, col )
            else
                SL_z = TC_z
            end if
            if ( TC_s .ge. G_s( row, col ) ) then
                SL_s = G_s( row, col )
            else
                SL_s = TC_s
            end if
            SL_c_out( row, col ) = SL_c * A( row, col )
            SL_z_out( row, col ) = SL_z * A( row, col )
            SL_s_out( row, col ) = SL_s * A( row, col )
            ! Distribute the soil particles from the element to lower adjacent
            ! elements.
            SL_c_in( ( row - 1 ) : ( row + 1 ), ( col - 1 ): ( col + 1 ) ) = &
                SL_c_in( ( row - 1 ) : ( row + 1 ), ( col - 1 ): ( col + 1 ) ) +&
                SL_c_out( row, col ) * m_weight
            SL_z_in( ( row - 1 ) : ( row + 1 ), ( col - 1 ): ( col + 1 ) ) = &
                SL_z_in( ( row - 1 ) : ( row + 1 ), ( col - 1 ): ( col + 1 ) ) +&
                SL_z_out( row, col ) * m_weight
            SL_s_in( ( row - 1 ) : ( row + 1 ), ( col - 1 ): ( col + 1 ) ) = &
                SL_s_in( ( row - 1 ) : ( row + 1 ), ( col - 1 ): ( col + 1 ) ) +&
                SL_s_out( row, col ) * m_weight
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! End sediment loss calculation module
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end do
    ! Outputs of a period
    ! 1. Total effective rainfall
    Rf_r( 1:nr, 1:nc, t ) = Rf_r( 1:nr, 1:nc, t ) + Rf( 1:nr, 1:nc )
    ! 2. Soil water infiltration capacity
    SW_c_r( 1:nr, 1:nc, t ) = SW_c_r(1:nr, 1:nc, t ) + SW_c( 1:nr, 1:nc )
    ! 3. Remained soil water content
    theta_r_r( 1:nr, 1:nc, t ) = SWr_d / SD / 1000.0d0  
    ! 4. Accumulated transport capacity of surface runoff
    TC_r( 1:nr, 1:nc, t ) = TC_r( 1:nr, 1:nc, t ) + TC( 1:nr, 1:nc )
    ! 5. Total inflow of surface runoff
    Q_in_r( 1:nr, 1:nc, t ) = Q_in_r( 1:nr, 1:nc, t ) + Q_in( 1:nr, 1:nc ) 
    ! 6. Total outflow of surface runoff
    Q_out_r( 1:nr, 1:nc, t ) = Q_out_r( 1:nr, 1:nc, t ) + Q_out( 1:nr, 1:nc )
    ! 7. Total inflow of interflow
    IF_in_r( 1:nr, 1:nc, t ) = IF_in_r( 1:nr, 1:nc, t ) + IF_in( 1:nr, 1:nc ) 
    ! 8. Total outflow of interflow
    IF_out_r( 1:nr, 1:nc, t ) = IF_out_r( 1:nr, 1:nc, t ) + IF_out( 1:nr, 1:nc )
    ! 9. Total suspended sediments of each soil particle size class.
    SS_c_r( 1:nr, 1:nc, t ) = SS_c_r( 1:nr, 1:nc, t ) + SS_c( 1:nr, 1:nc )
    SS_z_r( 1:nr, 1:nc, t ) = SS_z_r( 1:nr, 1:nc, t ) + SS_z( 1:nr, 1:nc )
    SS_s_r( 1:nr, 1:nc, t ) = SS_s_r( 1:nr, 1:nc, t ) + SS_s( 1:nr, 1:nc )
    ! 10. Total available sediments for transport 
    G_c_r( 1:nr, 1:nc, t ) = G_c_r( 1:nr, 1:nc, t ) + G_c( 1:nr, 1:nc )
    G_z_r( 1:nr, 1:nc, t ) = G_z_r( 1:nr, 1:nc, t ) + G_z( 1:nr, 1:nc )
    G_s_r( 1:nr, 1:nc, t ) = G_s_r( 1:nr, 1:nc, t ) + G_s( 1:nr, 1:nc )
    ! 11. Total sediment inputs
    SL_c_in_r( 1:nr, 1:nc, t ) = SL_c_in_r( 1:nr, 1:nc, t ) +&
        SL_c_in( 1:nr, 1:nc )
    SL_z_in_r( 1:nr, 1:nc, t ) = SL_z_in_r( 1:nr, 1:nc, t ) +& 
        SL_z_in( 1:nr, 1:nc )
    SL_s_in_r( 1:nr, 1:nc, t ) = SL_s_in_r( 1:nr, 1:nc, t ) +& 
        SL_s_in( 1:nr, 1:nc )
    ! 12. Total sediment outputs
    SL_c_out_r( 1:nr, 1:nc, t ) = SL_c_out_r( 1:nr, 1:nc, t ) +&
        SL_c_out( 1:nr, 1:nc )
    SL_z_out_r( 1:nr, 1:nc, t ) = SL_z_out_r( 1:nr, 1:nc, t ) +& 
        SL_z_out( 1:nr, 1:nc )
    SL_s_out_r( 1:nr, 1:nc, t ) = SL_s_out_r( 1:nr, 1:nc, t ) +& 
        SL_s_out( 1:nr, 1:nc )

    if ( Breaking( j ) .gt. 0.0d0 ) then
        t = t + 1
    end if

    if ( Init_point( j+1 ) .gt. 0.0d0 ) then
        init_counter = init_counter + 1
        SW_init = 1000.0d0 * theta_init( 1:nr, 1:nc, init_counter ) * SD
    else
        ! Daily update of remaining soil water 
        ! Remaining soil water is the last value of the event
        ! (Remaining soil water will be initial soil water of the next day.)
        SW_init = SWr_d
    end if
end do
! 13. Total sediment inputs and outputs of all sediment particle size classes
SL_in_r = SL_c_in_r + SL_z_in_r + SL_s_in_r
SL_out_r = SL_c_out_r + SL_z_out_r + SL_s_out_r
end subroutine DMMF
