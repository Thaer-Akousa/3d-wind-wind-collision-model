#ifndef BINARY_STARS_CONFIG_INTERN_H
#define BINARY_STARS_CONFIG_INTERN_H

#include <math.h>
#include <string.h>
#include <float.h>
#include <gsl/gsl_matrix.h>      // GSL routines for matrices
#include <gsl/gsl_integration.h> // GSL routines for integration
#include <gsl/gsl_odeiv2.h>      // GSL routines for ODE.
#include <gsl/gsl_interp.h>      // GSL routines for interpolation.
#include <gsl/gsl_spline.h>      // GSL routines for spline interpolation.
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#define noel 15
#define N_MAX 3000
#define MAX_ITER 1000
#define LMAX 100
#define TINY 1.0e-4
#define EPS 1.0e-5
#define NE0_MAX 1000
#define NE_MAX 12000
#define RMAX 100

#define T1 = 15000 * 1.38065812 / (1.60217725 * 1.0E-7)  ////! in keV [ kT/1.6e-9 ]

#define PI M_PI
// !!!!! #define PI_2 2*PI
// !!!!! #define PI_4 4*PI
#define PI_2 (2*PI)
#define PI_4 (4*PI)
#define sq(x)       ( (x) * (x) )
#define cube(x)     ( (x) * (x) * (x) )
#define dis2(x,y)   ( ( sq(x) ) + ( sq(y) ) )
#define dis(x,y)    (sqrt(dis2((x), (y))))

#define MIN(x,y) ((x >= y) ? y:x)
#define MAX(x,y) ((x >= y) ? x:y)

#define loop(i, n)  for(i=0; i<n; i++)
#define for_i(n)    for (int i=0; i<n; i++)
#define for_j(n)    for (int j=0; j<n; j++)


/* constants */

#define mp (1.6726231e-24)
#define me (9.109389754e-28)
#define rsol (6.96e10)
#define JMAX 20
#define FUNC(x) ((*func)(x))
#endif //BINARY_STARS_CONFIG_INTERN_H
