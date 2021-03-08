#include "config.h"
#include "numerical.h"

double interpolate(double x[], double y[], int N, double x0){
    N--;
    if (x0 <= x[0]) {
        return y[0] + (y[1] - y[0]) * (x0 - x[0])/(x[1] - x[0]);
    }
    if (x0 >= x[N-1]) {
        return y[N-1] + (y[N] - y[N-1]) * (x0 - x[N-1])/(x[N] - x[N-1]);
    }
    int i;
    for (i=0; i <= N; i++) {
        if (!(x0 > x[i])) break;
    }
    i--;
    return y[i] + (y[i+1] - y[i]) * (x0 - x[i]) / (x[i+1] - x[i]);
}

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]) {
    int i,k;
    double p,qn,sig,un;
    double u[n-1];

    if (yp1 > 0.99e30) {
        y2[1] = u[0] = 0.0;
    } else {
        y2[1] = -0.5;
        u[0] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) {
        sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i-1]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i-1]=(6.0*u[i-1]/(x[i+1]-x[i-1])-sig*u[i-2])/p;
    }
    if (ypn > 0.99e30) {
        qn=un=0.0;
    } else {
        qn=0.5;
        un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n] = (un-qn*u[n-2])/(qn*y2[n-1]+1.0);
    for (k=n-1;k>=1;k--) y2[k] = y2[k] * y2[k+1] + u[k-1];
}

double splint(double xa[], double ya[], double y2a[], int n, double x) {
    int klo, khi, k;
    double h, b, a, y;
    klo = 1;
    khi = n;
    while (khi - klo > 1) {
        k = (khi+klo) >> 1;
        if (xa[k] > x) khi=k;
        else klo=k;
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0) puts("Bad xa input to routine splint");
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;

    return a * ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

int odeint(double *y, double *t, const double tf, double (*ode_f)(double t, double y),
           double *h, const double epsabs, const double epsrel) {

    int status;

    int ode_function(double t, const double Y[], double f[], void *params) {
        (void)(t); // avoid unused parameter warning
        f[0] = ode_f(t, Y[0]);
        return GSL_SUCCESS;
    }


    const gsl_odeiv2_step_type *method = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(method, 1);
    gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(epsabs, epsrel);
    gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(1);
    gsl_odeiv2_system system = {ode_function, NULL, 1, NULL};

    double Y[1];
    Y[0] = *y;

// !!!!!!!
    while( *t < tf ) {

      status = gsl_odeiv2_evolve_apply(evolve, control, step, &system, t, tf, h, Y);
      if (status != GSL_SUCCESS) {
        printf( "Error in odeint. Exiting.\n" );
        exit( 0 );
      }

    }
// !!!!!!
    
    *y = Y[0];

    gsl_odeiv2_evolve_free(evolve);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);

    return status;
}

double normalize_angle(double theta) {
    if (theta > PI_2) return normalize_angle(theta - PI_2);
    if (theta < 0) return normalize_angle(theta + PI_2);
    return theta;
}

double quad(double (*f)(double x), const double a, const double b,
            const double epsabs, const double epsrel, int size_limit, double* abserr) {

    double function(double x, void *params){
        return f(x);
    }

    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(size_limit);
    double result;
    gsl_function F;
    F.function = &function;
    F.params = NULL;

    gsl_integration_qags(&F, a, b, epsabs, epsrel, size_limit, ws, &result, abserr);
    gsl_integration_workspace_free(ws);

    return result;
}


float qsimp(double (*func)(double), double a, double b)
{
	float trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	int j;
	float s,st,ost=0.0,os=0.0;

	for (j=1;j<=JMAX;j++) {
		st=trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (j > 5)
			if (fabs(s-os) < EPS*fabs(os) ||
				(s == 0.0 && os == 0.0)) return s;
		os=s;
		ost=st;
	}
	//nrerror("Too many steps in routine qsimp");
	return 0.0;
}

float trapzd(double (*func)(double), double a, double b, int n)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}