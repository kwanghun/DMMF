#ifndef R_DMMF_H
#define R_DMMF_H

/* definitions not involving SEXPs, plus _() */
#include <R_ext/RS.h>
#include <Rinternals.h>

void 
F77_SUB(dmmf)(double *dem, int *nr, int *nc, double *res, int *option, int *days, double *r, double *ri, int *r_type, double *et, double *p_c, double *p_z, double *p_s, double *theta_init, double *theta_sat, double *theta_fc, double *sd, double *k, double *p_i, double *n_s, double *d_a, double *cc, double *gc, double *imp, double *ph, double *d, double *nv, double *dk_c, double *dk_z, double *dk_s, double *dr_c, double *dr_z, double *dr_s, int *breaking, int *n_out, int *vc, int *init_point, double *a, double *rf_r, double *sw_c_r, double *theta_r_r, double *tc_r, double *q_in_r, double *q_out_r, double *if_in_r, double *if_out_r, double *ss_c_r, double *ss_z_r, double *ss_s_r, double *g_c_r, double *g_z_r, double *g_s_r, double *sl_c_in_r, double *sl_z_in_r, double *sl_s_in_r, double *sl_in_r, double *sl_c_out_r, double *sl_z_out_r, double *sl_s_out_r, double *sl_out_r);

void 
F77_SUB(ke)(double *ri, int *climate, double *rf, double *cc, double *ph);

void 
F77_SUB(mdinf)(double *m_block);

void 
F77_SUB(slope)(double *m_block, double *res, int *order);

void 
F77_SUB(sinkfill)(double *dem, int *nr, int *nc, double *res, double *boundary, double *min_angle, double *dem_nosink);

void 
F77_SUB(checkboundary)(double *dem, int *nr, int *nc, double *boundary);

#endif
