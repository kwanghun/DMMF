function KE( RI, Climate, Rf, CC, PH )
    ! This function is the calculator of the kinetic energy of the direct
    ! throughfall of a selected climate region with rainfall intensity. 
    implicit none
    ! The "RI" and the "Climate" are input parameters.
    ! The "RI" means the "Rainfall intensity" which has a type of real (8).
    ! The "Climate" means the type of climate region which has a type of integer.
    ! 0: Default, 1: North America, 2: North-west Europe, 3: Mediterranean,
    ! 4: West Mediterranean, 5: Tropic, 6: East Asia, 
    ! 7: Temperate Southern Hemisphere
    ! The "Rf" is the effective rainfall on the element.
    ! The "CC" is the canopy cover ratio of the element (0 - 1).
    ! The "PH" is the average plant height of the element.
    double precision, intent( in ) :: RI, Rf, CC, PH
    integer, intent ( in ) :: Climate
    ! DT is the "Direct rainfall", the ratio of effective rainfall on bare soil.
    ! KE_DT is the "Kinetic energy of the direct throughfall".
    ! KE_LD is the "Kinetic energy of the leaf drainage from plant".
    ! KE is the "Total kinetic energy" which is the sum of KE_DT and KE_LD. 
    double precision :: DT, LD, KE_DT, KE_LD, KE
    ! Direct rainfall is the amount of rainfall on the surface without canopy
    ! cover.
    LD = Rf * CC
    DT = Rf * ( 1.0d0 - CC )    
    ! KE_DT is the kinetic energy of the direct rainfall.
    select case ( Climate )
    case ( 1 )   ! 1: North America    
        KE_DT = DT * ( 11.87d0 + 8.73d0 * log10( RI ) )
    case ( 2 )   ! 2: North-west Europe    
        KE_DT = DT * ( 8.95d0 + 8.44d0 * log10( RI ) )
    case ( 3 )    ! 3: Mediterranean 
        KE_DT = DT * ( 9.81d0 + 11.25d0 * log10( RI ) )
    case ( 4 )    ! 4: West Mediterranean
        KE_DT = DT * ( 35.9d0 * ( 1.0d0 - 0.56d0 * exp( -0.034d0 * RI ) ) )
    case ( 5 )    ! 5: Tropic 
        KE_DT = DT * ( 29.8d0 - ( 127.5d0 / RI ) )
    case ( 6 )    ! 6: East Asia
        KE_DT = DT * ( 9.81d0 + 10.60d0 * log10( RI ) )
    case ( 7 )    ! 7: Temperate Southern Hemisphere
        KE_DT = DT * ( 29.0d0 * ( 1.0d0 - 0.6d0 * exp( -0.04d0 * RI ) ) )
    case ( 8 )    ! 8: Universial equation by Van Dijk et al. (2002)
        KE_DT = DT * ( 28.3d0 * ( 1.0d0 - 0.52d0 * exp( -0.042d0 * RI ) ) )
    case default   ! Default: Universial power law eq. by Shin et al. (2015) 
        KE_DT = DT * ( 10.3d0 * ( RI ** (2.0d0/9.0d0) ) )
    end select
    if ( PH .lt. 0.15d0) then
        KE_LD = 0.0d0
    else
        KE_LD = LD * ( ( 15.8d0 * sqrt( PH ) ) - 5.87d0 ) ! Modified from OV.
    end if
    KE = KE_DT + KE_LD
end function KE

