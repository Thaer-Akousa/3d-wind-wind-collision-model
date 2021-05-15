#ifndef BINARY_STARS_LIB_H
#define BINARY_STARS_LIB_H

/*calculates the  x-coordinate of the apex of the shock
*/
int shock_distance();


/* calculates the velocity of the wind of the first star
 */
double v1(double x);


/* Calculates the velocity of the wind of the second star
 */
double v2(double x);


/* The momentum equilibrium of the two winds
*/
double dm(double x);


/* Derivatives for the contact surface
*/
double dxdy(double x, double y);
double d2xdy2(int i);


/* Derivatives for the contact surface (equations 24 and 25 in notes)
*/
double dsigdy(double sigma, double y, int i);

double dsigdy1(double y, double sigma);
double dsigdy2(double y, double sigma);


/* An aux function for the both winds
*/
double dvdy(double x, double y, int i);


/* Derivative dL/dE for both winds (inverted equation 15) */
double dldE(double E, int i);

double dldE1(double E, double l);
double dldE2(double E, double l);


/* Cooling functions of both of the cooling layers consecutively
*/
double Q1(double E);
double Q2(double E);


/* Density along the line of sight
*/
//double density_z(double z, int i);
double density_z(double z);


/* Function that calculates the mean atomic weight of the medium
*/
void atomic_weight(double abundance[], double temperature, double* energy_grid, double* ren_grid, int size_lambda, double* mean_en, double* mean_n, double* sum_nz, double* rnh );


/* Calculation of the true anomaly vt and the distance between the components
*/
void anomaly(double phase, double* true_anomaly, double* separation);


/* Computes the shape of the contact surface, the normal and tangential
components of velocities of both winds in every point on the grid, the width of the two cooling layers
*/
void shock(double x0, double dy);


/* Direct integration of the left-hand side of equation (20)
*/
double lhs_20(double y, int nstar);

double make_3d();

double shift();

double skew();

double sub_x10(double x);

double add_x10(double x);

double calc_velocities();

double calc_sigma1();
double calc_sigma2();



double get_r();

double lhs_20_3d_1(double y);
double lhs_20_3d_2(double y);
double lhs_20_3d_1_xy(double y);
double lhs_20_3d_1_xy_bot(double y);
double lhs_20_3d_2_xy(double y);
double lhs_20_3d_2_xy_bot(double y);
double calc_width_1();
double calc_width_2();
double calc_psi();


#endif //BINARY_STARS_LIB_H


