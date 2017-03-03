DMMF_Simple <- 
function( W, L = W/cos(S), S, R, RI, ET, P_c, P_z, P_s, theta_init, theta_sat, theta_fc, SD, K, P_I, n_s, CC, GC, IMP, PH, D, NV, d_a = 0.005, DK_c = 0.1, DK_z = 0.5, DK_s = 0.3, DR_c = 1.0, DR_z = 1.6, DR_s = 1.5, Q_in = 0, IF_in = 0, SL_c_in = 0, SL_z_in = 0, SL_s_in = 0, R_type = 0)
{
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate surface area of an element
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    A = W * L
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Begin water quantity calculation module
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Because we only consider daily mean rainfall for sensitivity analysis, we don't need to specify number of raindays in certain period.
    R_eff = R * ( 1.0 - P_I ) * cos( S )
    # Soil moisture parameters
    SW_sat = 1000.0 * theta_sat * SD
    SW_fc = 1000.0 * theta_fc * SD
    SW_init = 1000.0 * theta_init * SD
    SW_c = ( 1.0 - IMP ) * ( SW_sat - SW_init - IF_in / A)
    # Total surface water before runoff process.
    R_sum = R_eff + Q_in / A
    if ( R_sum >= SW_c ) 
    {
        Q = R_sum - SW_c
    }else{
        Q = 0.0
    }
    Q_out = Q * A
    # Soil water budget and interflow 
    SW = max( 0.0, SW_init + IF_in / A + R_sum - Q - ET )
    SW_t <- max(SW - SW_fc, 0.0)
    IF_out = min( K * sin( S ) * SW_t * W, SW_t * A )

    # Remaining soil water content of entire soil profile
    theta_r = ( SW - IF_out / A ) / SD / 1000
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# End water quantity calculation module
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Begin sediment loss calculation module
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Detachment of soil particles by raindrop impact ( kg/m^2 )
######################################################################################
    if ( R_type == 0 )
    { 
        KE_dt = R_eff * ( 10.3 * RI^(2/9) )
    }else if( R_type == 1 )
    { 
        KE_dt = R_eff * ( 11.87 + 8.73 * log10 ( RI ) ) 
    }else if( R_type == 2 )
    {
        KE_dt = R_eff * ( 8.95 + 8.44 * log10( RI ) )
    }else if( R_type == 3 )
    {
        KE_dt = R_eff * ( 9.81 + 11.25 * log10 ( RI ) )
    }else if( R_type == 4 )
    {
        KE_dt = R_eff * ( 35.9 * (1.0 - 0.56 * exp( -0.034 * RI ) ) )
    }else if( R_type == 5 )
    {
         KE_dt = R_eff * ( 29.8 - ( 127.5 / RI ) )
    }else if( R_type == 6 )
    {
        KE_dt = R_eff * ( 9.81 + 10.60 * log10 ( RI ) )
    }else if( R_type == 7 )
    {
        KE_dt = R_eff * ( 29.0 * ( 1.0 - 0.6 * exp( -0.04 * RI ) ) )
    }else if( R_type == 8 )
    {
     KE_dt = R_eff * ( 28.3 * ( 1.0 - 0.52 * exp( -0.042 * RI ) ) )
    }
######################################################################################
    KE_ld = R_eff * max( 0.0, 15.8 * sqrt(PH) - 5.87 )
    KE = ( 1 - CC ) * KE_dt + CC * KE_ld 
    F_t = 0.001 * max( 0.0, 1.0 - ( GC * ( 1 - IMP ) + IMP ) ) * KE
    F_c = DK_c * P_c * F_t
    F_z = DK_z * P_z * F_t
    F_s = DK_s * P_s * F_t

    # Detachment of soil particles by surface runoff ( kg/m^2 )
    H = 0.001 * Q^1.5 * ( sin( S ) ) ^ 0.3 * max( 0.0, 1.0 - ( GC * ( 1 - IMP ) + IMP ) )  
    H_c  = DR_c * P_c * H
    H_z  = DR_z * P_z * H
    H_s  = DR_s * P_s * H

    # Suspended sediments which are delivered to runoff (SS) ( kg/m^2 )
    SS_c = F_c + H_c + SL_c_in / A
    SS_z = F_z + H_z + SL_z_in / A
    SS_s = F_s + H_s + SL_s_in / A

    # Surface runoff velocities ( m/s )
    # Reference velocities on the standard bare soil. (Reference)
    v_b = ( 0.005^0.67 ) * sqrt( tan( S ) ) / 0.015     
    # n' for flow velocity (n_t), n_s is the Manning's coefficients for
    # surface condition.
    g = 9.80665
    n_t = sqrt( ( n_s * n_s ) + ( 0.5 * D * NV * ( d_a^(4.0/3.0) ) ) / g )  
    # Flow velocity
    v = ( ( d_a^( 2.0 / 3.0 ) ) * sqrt( tan( S ) ) ) / n_t

    # Particle fall velocities (vs_x) and particle fall number (Nf_x)
    if ( v  > 0.0 ) 
    {

        rho_s = 2650.0
        rho = 1000.0
        eta = 0.0015
        di_c = 0.000002
        di_z = 0.00006 
        di_s = 0.0002 
        vs_c = di_c * di_c * ( rho_s - rho ) * g / ( 18.0 * eta )
        Nf_c = vs_c * L / ( v * d_a )
        vs_z = di_z * di_z * ( rho_s - rho ) * g / ( 18.0 * eta )
        Nf_z = vs_z * L / ( v * d_a )
        vs_s = di_s * di_s * ( rho_s - rho ) * g / ( 18.0 * eta )
        Nf_s = vs_s * L / ( v * d_a )

        # Percentage of detached particles that is deposited ( DEP_x ) (%)
        DEP_c = min( ( 0.441 * ( Nf_c^0.29 ) ), 1.0 )
        DEP_z = min( ( 0.441 * ( Nf_z^0.29 ) ), 1.0 )
        DEP_s = min( ( 0.441 * ( Nf_s^0.29 ) ), 1.0 )
    }else{
        DEP_c = 1.0
        DEP_z = 1.0
        DEP_s = 1.0
    }

    # Delivery of detached particles to runoff ( kg/m^2 )
    G_c = SS_c * ( 1.0 - DEP_c ) 
    G_z = SS_z * ( 1.0 - DEP_z ) 
    G_s = SS_s * ( 1.0 - DEP_s ) 

    # RMMF-C factor with velocities 
    # the how the surface runoff velocity is faster or slower than that
    # of standard condition.
    if ( v_b == 0.0 )
    {    
        C = 0.0
    }else{
        C = ( v / v_b )  
    }

    # Transport capacity of the runoff ( kg/m^2 )
    TC = 0.001 * C * Q * Q * sin( S )
    TC_c = P_c * TC
    TC_z = P_z * TC
    TC_s = P_s * TC

    # Sediment balance
    if ( TC_c >= G_c )
    {
        SL_c = G_c
    }else{
        SL_c = TC_c
    }
    if ( TC_z >= G_z )
    {
        SL_z = G_z
    }else{
        SL_z = TC_z
    }
    if ( TC_s >= G_s )
    {
        SL_s = G_s
    }else{
        SL_s = TC_s
    }

    SL = SL_c + SL_z + SL_s
# Total weight of sediment loss for each soil particle size classes
    SL_c_out = SL_c * A
    SL_z_out = SL_z * A
    SL_s_out = SL_s * A
# Outputs
    Output <- data.frame("Q_out" = Q_out, "IF_out" = IF_out, "theta_r" = theta_r, "SL_c_out" = SL_c_out, "SL_z_out" = SL_z_out, "SL_s_out" = SL_s_out, "A" = A)
    return(Output)
}

