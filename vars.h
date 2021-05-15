#ifndef BINARY_STARS_VARS_H
#define BINARY_STARS_VARS_H

#include "config.h"

#define NPH_MAX 50


void init_test_values();

/*input files*/
void cooling_read(char *fname, double *elambda, double *lambda, double *nergrid, int *nlambda);
void planar_read(char *fname, double *e0min, double *e0bin, int *ne0, double *emin, double *ebin, int *ne, double *ear, double sp_t[NE0_MAX][NE_MAX]);
void opacity_read(char *fname, double *eabs, double *absorp, int *nabs);
/* functions for printing results in output files */

void three_d_results(/*FILE *output_file_2*/);

double coeff;


/* Iterations' limits and regulation*/
int kb1, kv1, km1;
int kb2, kv2, km2;
int n1, n2;
int flag_eta,flag_abs,k;


/* main Parameters */
double dl,ds,ds1;
double xxc,yyc,zzc;
double d1,d2,dd1,dd2,dc;
double cosb,sinb;
double v1n0,v2n0;
double rho10,rho20;
double theta, dtheta;
double coldens11,coldens12,coldens21,coldens22;
double r1,r1x,r1y,r1z,r2,r2x,r2y,r2z,nx,ny,nz,nx1,ny1,nz1,lx1,lx,lxint;
double coseta,coseta1;
double psi, psimax;
double ccphi, ssphi;
double t;
double yy,zz;
double ps1,ps2;
double dum1,dum2,dummy,dummy1,dummy2, dumdum1,dumdum2,dumdum;
double tau1,tau2;
double e1,e2;
double photons[NE_MAX];
double Ierror;


/* Parameters for winds' geometry */
int N;  // number of knots
int nc;// number of knots
double xc[N_MAX], yc[N_MAX], dxdy_c[N_MAX], d2xdy2_c[N_MAX], d2xcdy[N_MAX], x22c[N_MAX];
double cphi[N_MAX], sphi[N_MAX], v1t[N_MAX], v1n[N_MAX], v2t[N_MAX], v2n[N_MAX];
double d1l[N_MAX], d2l[N_MAX];
double sig1[N_MAX], sig2[N_MAX];
double vel1,vel2;
double d2xdy20;
double dy;

/* Stellar Parameters */
int verbosity;
int nb_ph;
int nb_el;
double ph[NPH_MAX];
double incl;
double phase;
double e;
double omega;
double d0;
double vt;
double vorb_0;
double vorb;
double orb_per;

double p; //impact parameter
double z1,z2;
double d1s, d2s;
double lsx, lsy, lsz;
double x10;
double d;
double rstar2;
double rstar1;
int nstar;


/* wind Parameters */
double mdot1_min, mdot1_max, dmdot1; // Minimal, maximal Mdot for 1st component, grid step.
double vinf1_min, vinf1_max, dvinf1; // Minimal, maximal V_infty for 1st component, grid step.
double beta1_min, beta1_max, dbeta1; // Minimal, maximal beta in beta-low for 1st component, grid step.
double mdot1,beta1,vinf1;

double mdot2_min, mdot2_max, dmdot2; // Minimal, maximal Mdot for 2nd component, grid step.
double vinf2_min, vinf2_max, dvinf2; // Minimal, maximal V_infty for 2nd component, grid step.
double beta2_min, beta2_max, dbeta2; // Minimal, maximal beta in beta-low for 2nd component, grid step.
double mdot2,beta2,vinf2;

double eta0;




/* Parameters of the spectra */
double e0min1, e0bin1, emin1, ebin1, e0min2, e0bin2, emin2, ebin2;
static double sp_t1[NE0_MAX][NE_MAX], sp_t2[NE0_MAX][NE_MAX];
int ne01, ne1, ne02, ne2;
double e0max1, e0max2;
double ear1[NE_MAX];
double eee1,eee2;
double cross_abs1[NE_MAX], cross_abs2[NE_MAX],spectrum[NE_MAX], spectr1[NE_MAX], spectr2[NE_MAX];
double spectrum_int[NE_MAX];
double ekin1, ekin2, ekin11,ekin12;

/* opacity parameters */
double eabs1[NE_MAX],absorp1[NE_MAX], eabs2[NE_MAX],absorp2[NE_MAX];
int nabs1, nabs2;

/* cooling data file */
int nlambda1, nlambda2;
double elambda1[N_MAX], lambda1[N_MAX], nergrid1[N_MAX];
double elambda2[N_MAX], lambda2[N_MAX], nergrid2[N_MAX];

/* atomic weight */
double angr[noel], abund[noel], at[noel];
double sum_nz1, sum_nz2, mu_av1, mu_av2, sum_ions;
double nhr1, nhr2, mu1, mu2;

/* 3d model coordinates*/
double x_skewed[NE0_MAX], y_skewed[NE0_MAX];
double x_3d[NE0_MAX][NE0_MAX], y_3d[NE0_MAX][NE0_MAX], z_3d[NE0_MAX][NE0_MAX];
double theta_3d[NE0_MAX];
double x_3d_1[NE0_MAX][NE0_MAX], y_3d_1[NE0_MAX][NE0_MAX], z_3d_1[NE0_MAX][NE0_MAX];
double x_3d_2[NE0_MAX][NE0_MAX], y_3d_2[NE0_MAX][NE0_MAX], z_3d_2[NE0_MAX][NE0_MAX];
double r_3d[NE0_MAX][NE0_MAX],r_2[NE0_MAX], r_3d_2[NE0_MAX][NE0_MAX], r_3d_3[NE0_MAX][NE0_MAX];
/* 3d model parameters*/
double psi_deg;

/* 3d model velocities*/
double v1t_3d[NE0_MAX][NE0_MAX], v2t_3d[NE0_MAX][NE0_MAX], v1n_3d[NE0_MAX][NE0_MAX]
           ,v2n_3d[NE0_MAX][NE0_MAX], v1_x, v1_y, v1_z, v2_x, v2_y, v2_z;

/* 3d model densities */
double sig1_3d[NE0_MAX][NE0_MAX], sig2_3d[NE0_MAX][NE0_MAX];

double r_3d_top[NE0_MAX],dis_to_x_stag_top[NE0_MAX], v1n_3d_top[NE0_MAX],dxdy_3d_top[NE0_MAX],v1t_3d_top[NE0_MAX], yc_top[NE0_MAX], sig1_3d_top[NE0_MAX];
double  y_3d_top[NE0_MAX],x_3d_top[NE0_MAX],z_3d_top[NE0_MAX];

int kkk, stag_index;
double dis_to_x_stag[NE0_MAX], d2xdy2_3d[NE0_MAX];


double r_3d_bot[NE0_MAX],dis_to_x_stag_bot[NE0_MAX], v1n_3d_bot[NE0_MAX],dxdy_3d_bot[NE0_MAX],v1t_3d_bot[NE0_MAX], yc_bot[NE0_MAX], sig1_3d_bot[NE0_MAX];
double  y_3d_bot[NE0_MAX],x_3d_bot[NE0_MAX],z_3d_bot[NE0_MAX];


double r_3d_top_2[NE0_MAX],dis_to_x_stag_top_2[NE0_MAX], v2n_3d_top[NE0_MAX],dxdy_3d_top_2[NE0_MAX],v2t_3d_top[NE0_MAX], yc_top_2[NE0_MAX], sig2_3d_top[NE0_MAX];
double  y_3d_top_2[NE0_MAX],x_3d_top_2[NE0_MAX],z_3d_top_2[NE0_MAX];

int kkk_2, stag_index_2;
double dis_to_x_stag_2[NE0_MAX], d2xdy2_3d_2[NE0_MAX];


double r_3d_bot_2[NE0_MAX],dis_to_x_stag_bot_2[NE0_MAX], v2n_3d_bot[NE0_MAX],dxdy_3d_bot_2[NE0_MAX],v2t_3d_bot[NE0_MAX], yc_bot_2[NE0_MAX], sig2_3d_bot[NE0_MAX];
double  y_3d_bot_2[NE0_MAX],x_3d_bot_2[NE0_MAX],z_3d_bot_2[NE0_MAX];

/* 3d model widths of cooling layers */
double d1l_3d[NE0_MAX][NE0_MAX], d2l_3d[NE0_MAX][NE0_MAX];

/*3d normal vectors*/
double nor1_x[NE0_MAX][NE0_MAX], nor1_y[NE0_MAX][NE0_MAX], nor1_z[NE0_MAX][NE0_MAX];
double nor2_x[NE0_MAX][NE0_MAX], nor2_y[NE0_MAX][NE0_MAX], nor2_z[NE0_MAX][NE0_MAX];

/*directional cosi and sin*/
double cos_x1, cos_y1, cos_z1;
double cos_x2, cos_y2, cos_z2;

/*aux params*/
double sig1_interp[NE0_MAX], sig2_interp[NE0_MAX], y_interp[NE0_MAX];
double x_aux[NE0_MAX], y_aux[NE0_MAX],z_aux[NE0_MAX];

double cone_x[NE0_MAX], cone_y[NE0_MAX];


#endif //BINARY_STARS_VARS_H
