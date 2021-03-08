#ifndef BINARY_STARS_NUMERICAL_H
#define BINARY_STARS_NUMERICAL_H


/* The Numerical Approximation 'numerical' module provides functions for:
 * - ODE solving
 * - Spline Approximation
 * - Polynomial Interpolation
 * - Linear Approximation
 * These functions can be wrapper for functions from the gsl module
 * or implemented directly, depending on our need.
 */

double interpolate(double x[], double y[], int N, double x0);


void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);


double splint(double xa[], double ya[], double y2a[], int n, double x);


int odeint(double *y, double *t, const double tf, double (*ode_f)(double t, double y),
           double *h, const double epsabs, const double epsrel);
//int ode_function(double t, const double Y[], double f[], void *params);

double normalize_angle(double theta);

double quad(double (*f)(double x), const double a, const double b,
            const double epsabs, const double epsrel, int size_limit, double* abserr);

float qsimp(double (*func)(double), double a, double b);


float trapzd(double (*func)(double), double a, double b, int n);


#endif //BINARY_STARS_NUMERICAL_H
