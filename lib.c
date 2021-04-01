#include "config.h"
#include "lib.h"
#include "vars.h"
#include "numerical.h"

double v1(double x) {
    return vinf1 * pow((1.0 - rstar1/x), beta1);
}

double v2(double x) {
    return vinf2 * pow((1.0 - rstar2/x), beta2);
}

double dm(double x) {
    return (mdot1 * v1(x) / sq(x)) - (mdot2 * v2(d-x) / sq(d-x));
}

int shock_distance() {
    double x1,x2;
    x1 = 0.5 * d;
    x2 = x1;
    double tmp;

    if (dm(x2) > 0.0) {
        do {

            tmp = d - rstar2 - x2;
            x2 += tmp / 10.0;
        } while ( dm(x2) > 0.0 && tmp / rstar2 > TINY );
    } else {
        do {
            tmp = x2 - rstar1;
            x2 -= tmp / 10;
        } while ( dm(x2) < 0 && tmp / rstar1 > TINY );

    }


    if ((dm(x1) > 0) && (dm(x2) > 0)) {
        eta0 = (mdot1 * vinf1) / (mdot2 * vinf2);
        printf("Eta is too large: %e - contact surface apex lies on 2nd star\n", eta0);
        flag_eta = 2;
        return flag_eta;
    }

    if ((dm(x1) < 0) && (dm(x2) < 0)) {
        eta0 = (mdot1 * vinf1) / (mdot2 * vinf2);
        printf("Eta is too small %e - contact surface apex lies on 1st star\n", eta0);
        flag_eta = 1;
        return flag_eta;
    }

    // Make sure x1 < x2. dm at x < x_sol is positive, at x > x_sol negative.
    if (x1 > x2) {
        x10 = x2;
        x2 = x1;
        x1 = x10;
    }

    do {
        x10 = (x1 + x2) / 2.0;
        if (dm(x10) > 0) x1 = x10;
        else x2 = x10;
    } while (fabs((x1 - x2) / x2) >= EPS);
    x10 = (x1 + x2) / 2.0;

    return 0;
}


double dxdy(double y, double x) {

/*
 1st derivative of the contact surface.
*/

    double eta;

    if (fabs(eta0 - 1.0) < DBL_MIN) return 0.0;

    d1s = dis2(x, y);
    d2s = dis2(d-x,y);


    vel1 = v1(sqrt(d1s));
    vel2 = v2(sqrt(d2s));


    eta = (mdot1*vel1)/(mdot2*vel2);
    if (y < 0.01*d) {
        return d2xdy20 * y;
    }

    return (x - d1s*d/(d1s+sqrt(eta)*d2s)) / y;
}

double d2xdy2(int i) {

/*
 2nd derivative of the contact surface at grid knots.
*/

    double x, y, r, dxdy_y, d1s, d2s, eta, tmp1, tmp2;

    if (fabs(eta0 - 1.0) < DBL_MIN) return 0.0;

    x = xc[i];
    y = yc[i];
    dxdy_y = (y <= 0.01*d) ? dxdy_c[0] : dxdy_c[i] / y; // Equiv. to above.
    d1s = sq(x) + sq(y);
    d2s = sq(d-x) + sq(y);
    vel1 = v1(sqrt(d1s));
    vel2 = v2(sqrt(d2s));
    eta = (mdot1 * vel1) / (mdot2 * vel2);

    tmp1 = sq(d1s + sqrt(eta)*d2s);
// This is the term when v1, v2 are constant.
    r = -2.0 * d * sqrt(eta) * (d2s - d1s + (d2s*x +d1s*(d-x)) * dxdy_y);
    r /= tmp1;
// This is the term nonequal to zero when v(r) are variable.
    tmp1 = d * d1s * d2s * sqrt(eta) / 2.0 / tmp1;
    tmp2 = dvdy(x,y,1)*(1 + x*dxdy_y) - dvdy(x,y,2)*(1 - (d-x)*dxdy_y);
    return r + tmp1*tmp2;
}

double dsigdy(double y, double sig, int i) {

// Derivative dsig/dy for the contact surface.

    double x;
    double dxdy;

    double di;
    double vel;
    double cd;
    double sd;
    double ctheta;
    double stheta;
    double phi;
    double cphi,sphi;
    double delta;
    double rho;
    double dum1,dum2;
    double d2xdy2;
    double dcthetady;

/*
 NOTES: interpolate below could (and possibly should) be replaced by
 spline-interpolation. Leave for the future.

 Important note from my fortran code:

 In principle, spline work as this:

 call "spline" function which computes the coefficients of spline
 interpolation.

 E.g. (in our case of contact surface yc is argument, xc is function):

 spline( yc, xc, nc, 1.e30, 1.e30, spline_coeff )

 Here 1.e30 are parameters which tell spline to use "free" ends of spline
 function.

 To use the coefficients to actually compute the interpolated
 function at an arbitrary point yy, call

 splint( yc, xc, spline_coeff, nc, yy, xx )

 However, in my code I do NOT  derive these coefficients with "spline" as
 the asymptotics I use near y=0 makes a discontinuity in the second
 derivative which results in a weird spline near y=0. Instead, I just use
 the second derivative computed myself, which avoids these problems.
*/

//    x = splint(yc, xc, d2xdy2_c, N, y);
    x = interpolate(yc, xc, nc, y);
//    dxdy = splint(yc, dxdy_c, d2xcdy, N, y);
    dxdy = interpolate(yc, dxdy_c, nc, y);

    if( i == 1 ) {

      di = sqrt( sq(x)+sq(y) );
      vel = v1(di);

      sphi = 1.0/sqrt( sq(dxdy)+1.0 );

      if( fabs(dxdy) < DBL_MIN ){
          cphi = 0.0;
      } else {
          cphi = dxdy*sphi;
      }

      cd = x/di;
      sd = y/di;

      ctheta = cd*cphi + sd*sphi;
      stheta = cd*sphi - sd*cphi;

      phi = asin(sphi);
      delta = asin(sd);

      rho = mdot1/( PI_4*sq(di)*vel );

// Right-hand side of the equation.
      if( y <= 0.01*d ) {
        dum1 = x10*d2xdy20;
        dum1 = 1.0+dum1/(1.0+dum1);
        dum1 = -rho*y/6.0/x10*dum1;
//        printf("i= %2d x=%15.9e y=%15.9e dum1= %15.9e\n", i, x, y, dum1 );
        return dum1;
      }

// Linear interpolation of the second derivative.
      d2xdy2 = interpolate(yc, d2xdy2_c, nc, y);

      dcthetady = sin(delta-phi) * ( (y*dxdy-x)/sq(di) - d2xdy2/(1.0+sq(dxdy)) );

      dum2 = dvdy(x,y,i)*( y+x*dxdy );

      dum1 = rho*stheta/sphi/ctheta - sig*( dcthetady/ctheta+dum2+1.0/y );

      return dum1;

    }

    if( i == 2 ) {

      di = sqrt( sq(d-x)+sq(y) );
      vel = v2(di);

      sphi = 1.0/sqrt( sq(dxdy)+1.0 );

      if( fabs(dxdy) < DBL_MIN ) {
        cphi = 0.0;
      } else {
        cphi = dxdy*sphi;
      }

      cd = (d-x)/di;
      sd = y/di;

      ctheta = -cd*cphi + sd*sphi;
      stheta = cd*sphi + sd*cphi;

      phi = asin(sphi);
      delta = asin(sd);

      rho = mdot2/( PI_4*sq(di)*vel);

      if (y <= 0.01*d) {
        dum1 = (x10-d)*d2xdy20;
        dum1 = 1.0+dum1/( 1.0+dum1 );
        dum1 = -rho*y/6.0/(d-x10)*dum1;
        return dum1;
      }

// Linear interpolation of the second derivative. Gives some (not big) errors
// in the resulting sigma2 at large y.
      d2xdy2 = interpolate(yc, d2xdy2_c, nc, y);

      dcthetady = sin(delta+phi)*( (y*dxdy-(x-d))/sq(di) - d2xdy2/(1.0+sq(dxdy)) );

      dum2 = dvdy(x,y,i) * (y+(x-d)*dxdy);

      dum1 = rho*stheta/sphi/ctheta - sig*(dcthetady/ctheta+dum2+1.0/y);

      return dum1;

    }

    printf( "dsigdy: i must be eq. to 1 or 2. Actual value %d. Exiting.\n", i );
    exit(0);

    return 0.0;
}


double dsigdy1(double y, double sig) {
    return dsigdy(y, sig, 1);
}

double dsigdy2(double y, double sig) {
    return dsigdy(y, sig, 2);
}



double dvdy(double x, double y, int i) {

/*
 First derivative of velocity by y dv/dy.
*/

    double xi;
    double dis;
    double beta;
    double rstar;
    double dummy;

    xi = (i==1) ? x:d-x; // x1 = x; x2 = x -d
    dis = sqrt((sq(xi)) + (sq(y)));
    beta = (i==1) ? beta1:beta2;
    rstar = (i==1) ? rstar1:rstar2;

    if (fabs(beta) < DBL_MIN) return 0.0;

    dummy = beta/(1-rstar/dis)*rstar/cube(dis);

    return dummy;
}

double dldE(double E, int i) {
    double sum_total, sum_xz, QE, mu_av, sum_nz, ner, mu;
    double tmp;
    if (i==1) {
        if (E < exp(elambda1[0])) ner = nergrid1[0];
        else ner = interpolate(elambda1, nergrid1, nlambda1, log(E));
    } else {
        if (E < exp(elambda2[0])) ner = nergrid2[0];
        else ner = interpolate(elambda2, nergrid2, nlambda2, log(E));
    }

    sum_nz = (i==1) ? sum_nz1:sum_nz2;
    mu_av = (i==1) ? mu_av1:mu_av2;
    sum_ions = mu_av * sum_nz;
    sum_total = sum_ions + ner * me / mp;
    mu = sum_total / (sum_nz + ner);

    sum_xz = (i==1) ? nhr1:nhr2;
    sum_xz *= ner;

    QE = (i==1) ? Q1(E):Q2(E);
    tmp = sq(E) / (QE * coeff * cube(mu) * sum_xz * nhr1);
    return sq(E) / (QE * coeff * cube(mu) * sum_xz * nhr1);
}

double dldE1(double E, double l) {
    return dldE(E, 1);


}

double dldE2(double E, double l) {
    return dldE(E, 2);
}

double Q1(double E) {
    if (E < exp(elambda1[0])) return exp(lambda1[0]);
    return exp(interpolate(elambda1, lambda1, nlambda1, log(E)));
}

double Q2(double E) {
    if (E < exp(elambda2[0])) return exp(lambda2[0]);
    return exp(interpolate(elambda2, lambda2, nlambda2, log(E)));
}

/*
double density_z(double z, int i) {
    double r = dis(z, p);
    double mdot, vel;
    if (i==1) {
        mdot = mdot1;
        vel = v1(r);
    } else {
        mdot = mdot2;
        vel = v2(r);
    }
    return mdot/ (PI_4 * sq(r) * vel);
}
*/

double density_z(double z) {
    double r = dis(z, p);
    double mdot, vel;
    if (nstar==1) {
        mdot = mdot1;
        vel = v1(r);
    } else {
        mdot = mdot2;
        vel = v2(r);
    }
    return mdot/ (PI_4 * sq(r) * vel);
}
// WHEN CALLING: atomic_weight(double abund[], double t, double elambda_array1, double nergrid_array1, double nlambda1, double &mu1, double &mu_av1, double &sum_nz1, double &nhr1)
// OR: atomic_weight(double abund[], double t, double elambda_array1, double nergrid_array1, double nlambda1, double &mu1, double &mu_av1, double &sum_nz1, double &nhr1)
//atomic_weight(double abund[], double t, double *elambda, double *nergrid, double nb_lambda, double *mu1, double *mu_av1, double *sum_nz1, double *nhr1) maybe this will work

// Calling:
// atomic_weight( abund, t, elambda, nergrid, nlambda, &mu, &mu_av, &sum_nz, &nhr );


void atomic_weight(double *abund, double t, double *elambda, double *nergrid,
                   int nlambda, double *mu, double *mu_av, double *sum_nz, double *nhr) {
    double tmp, ner;

    *sum_nz = 0.0;
    *mu = 0.0;
    *mu_av = 0.0;

    for_i(noel) {
        tmp = pow(10.0, angr[i] - 12.0) * abund[i];
        *sum_nz += tmp;
        *mu_av += tmp * at[i];
    }
    *nhr = abund[0] / (*sum_nz);
    *mu_av /= (*sum_nz);

    if (t < exp(elambda[0])) ner = nergrid[0];
    else ner = interpolate(elambda, nergrid, nlambda, log(t));

    sum_ions = (*mu_av) * (*sum_nz);
    *mu = (sum_ions + ner * me / mp) / ((*sum_nz) + ner);

    return;

}

void anomaly(double ps, double *vt, double *d12) {
    double m, ea1, ea2, d21, ea, b;

    ea1 = m = ps * PI_2;
    do {
        ea2 = m + e * sin(ea1);
        d21 = ea2 - ea1;
        ea1 = ea2;
    } while (fabs(d21) >= EPS);
    ea = ea1;
    b = sin(ea*0.5) / cos(ea*0.5) * sqrt((1+e)/(1-e));
    *vt = 2 * atan(b);
    *d12 = (1 - sq(e)) / (1 + e * cos(*vt));
}

void shock(double x10, double dy) {
    double ymax;
    int i;
    double eps1, eps, h1, dd1, dd2, tmp;

    double y1, y2, x2;
    double vel1;
    double vel2;
    double a,b,c;
    double dummy1, dummy2;
    int status1,status2,status3,status4;

    short fast=0;

    vel1 = v1(x10);
    vel2 = v2(d - x10);
    eta0 = mdot1 * vel1 / (mdot2 * vel2);
    if (fabs(eta0 - 1) <= 1.0e-5) eta0 = 1.0;


/*

 Vars:
 xc[],yc[] - coordinates of the contact surface at grid knots.
 dxdy_c[]  - first derivative of the contact surface dx/dy at grid knots.
 d2xdy2_c  - second derivative of the contact surface at grid knots.

*/

// Second derivative of the contact surface at y=0.
    if( fabs(eta0-1.0) < DBL_MIN ) {
      d2xdy20 = 0.0;
    } else {
      a = d*sqrt(eta0)/2.0/sq( 1.0+sqrt(eta0) );
      b = dvdy(x10, 0.0, 1);
      c = dvdy(x10, 0.0, 2);

      dummy1 = 2.0*(eta0-1.0)/d/sqrt(eta0) + a*(b-c);
      dummy2 = 3.0 - a*(b*x10+c*(d-x10));
      d2xdy20 = dummy1/dummy2;
    }

    xc[0] = x10;
    yc[0] = 0.0;
    d2xdy2_c[0] = d2xdy20;
    ymax = yc[0] + dy * (nc - 1);
    eps = 1.0e-6;
    y2 = yc[0];
    x2 = xc[0];
    for (i=1; i<nc; i++) {
        y1 = y2;
        y2 = y1 + dy;
        if (fabs(eta0 - 1.0) < DBL_MIN) {
            xc[i] = x10;
            yc[i] = y2;
            dxdy_c[i] = 0.0;
            d2xdy2_c[i] = 0.0;
        } else {
            h1 = (y2 - y1) / 20;
            status1 = odeint(&x2, &y1, y2, dxdy, &h1, eps, 0.0);

            xc[i] = x2;
            yc[i] = y2;
            dxdy_c[i] = dxdy(yc[i], xc[i]);
        }
        d2xdy2_c[i] = d2xdy2(i);
    }


    spline(yc, dxdy_c, nc, d2xdy2_c[0], d2xdy2_c[nc-1], d2xcdy);
    spline(yc, d2xdy2_c, nc, 2.0e30, 2.0e30, x22c);
    //spline(yc, xc, N, 0.0, 2.0e30, xx2c);

   //for_i(nc) printf("%e\t\t%e\t\t%e\t\t%e\n", yc[i], xc[i], dxdy_c[i], d2xdy2_c[i]);

    double rho1, rho2;
//first wind
    rho1 = mdot1 / (PI_4 * sq(x10) * v1(x10));    //wind denisty
    sig1[0] = rho1 * x10 / (2 * (1 + x10 * d2xdy20));   //surface density
//second wind
    rho2 = mdot2 / (PI_4 * sq(d - x10) * v2(d - x10));
    sig2[0] = rho2 * (d - x10) / (2 * (1 - (d - x10) * d2xdy20));

    eps = 1.0e-5;
    y2 = yc[0];
    dd1 = sig1[0];
    dd2 = sig2[0];
    h1 = (yc[1] - yc[0]) / 20;

    for (int i=1; i<nc; i++) {
        y1 = y2;
        y2 = y1 + dy;
        tmp = y1;
        h1 = (yc[1] - yc[0]) / 20;
        status2 = odeint(&dd1, &y1, y2, dsigdy1, &h1, eps, 0.0);
        sig1[i] = dd1;
        h1 = (yc[1] - yc[0]) / 20;
        status2 = odeint(&dd2, &tmp, y2, dsigdy2, &h1, eps, 0.0);
        sig2[i] = dd2;
    }

    double d1, d2, cosb, sinb, rho10, e0, ll1, ll2, e_start, dl2, dl1;
    int iter;
    double xxc, yyc, v1n0, dl, rho20, v2n0;

// Compute the widths of the two shocks.

    for( i=0; i<nc; i++) {
        if (fabs(dxdy_c[i]) <DBL_MIN){     //if the first derivative = 0
            cphi[i] = 0;
        } else {
            cphi[i] = 1 / sqrt(1 + 1/sq(dxdy_c[i]));
        }

        sphi[i] = sqrt(1 - sq(cphi[i]));

        d1 = dis(xc[i], yc[i]);  // Distance from the point on c.s. to center of 1st star.
        d1s = sq(d1);
        vel1 = v1(d1); // Velocity of 1st wind at point of c.s.
        cosb = xc[i] / d1;
        sinb = yc[i] / d1;
        v1t[i] = vel1 * (cosb * cphi[i] + sinb * sphi[i]);
        //Tangential and normal velocities.
        v1n[i] = vel1 * (cosb * sphi[i] - sinb * cphi[i]);

        d2 = dis((d - xc[i]), yc[i]);  // Distance from the point on c.s. to center of 2nd star.
        d2s = sq(d2);
        vel2 = v2(d2);                // Same as above, for 2nd wind.
        cosb = (d - xc[i]) / d2;
        sinb = yc[i] / d2;
        v2t[i] = vel2 * (- cosb * cphi[i] + sinb * sphi[i]);
        v2n[i] = vel2 * (cosb * sphi[i] + sinb * cphi[i]);

        rho10 = mdot1 / (PI_4 * d1s * vel1); // Density of 1st wind at point of c.s.
        e0 = 1.957441607 * mu1 * sq(v1n[i]); // Energy of normal component of the wind.

        coeff = 8.3245e3/d0/sq(mu_av1)*rho10*cube(v1n[i]);

// Compute first estimate of shock width based of energy e0 estimated by v1n0 at the c.s.
// As the shock front is at some distance from the c.s., v1n is different so I have to refine
// this estimate - see below.
        ll1 = 0.0;
        h1 = e0 / 300.0;
        eps = 1.0e-5;
        eps1 = 1.0e-5;
        e_start = 0.0;
        status3 = odeint(&ll1, &e_start, e0, dldE1, &h1, eps, 0.0);
        dl2 = ll1;


/*
 Make iterations to refine the estimate of the shock width:

 1. At the current width dl2, compute velocity, energy e0, end coeff at the
    shock front position;
 2. With these, solve diff eq. to find a new width.
 3. Set the new width as average of the old and new ones.
 4. Repeat 1-3 until convergency.

*/

        if (fast) { // If fast, just take this first estimate as an answer.
            d1l[i] = dl2;
        } else {    // Refine the estimate.
            dl2 = MIN(dl2, 0.9 * (x10 - rstar1));

            eps = 1.0e-5;
            iter = 0;
            dl1 = 0.0;
            while ((fabs((dl2-dl1)/dl2) > eps) && (iter < MAX_ITER)) {
                iter++;
                dl2 = (dl1 + dl2) / 2.0;
                xxc = xc[i] - dl2 * sphi[i];
                yyc = yc[i] + dl2 * cphi[i];
                d1 = dis(xxc, yyc);
                cosb = xxc/d1;
                sinb = yyc/d1;
                vel1 = v1(d1);
                v1n0 = vel1 * (cosb * sphi[i] - sinb * cphi[i]);
                rho10 = mdot1 / (PI_4 * sq(d1) * vel1);
                e0 = 1.957441607 * mu1 * sq(v1n0);
                coeff = 8.3245e3 / (d0 * sq(mu_av1)) * rho10 * cube(v1n0);
                e_start = 0.0;
                h1 = e0/300.0;
                ll1 = 0.0;

                status3 = odeint(&ll1, &e_start, e0, dldE1, &h1, eps, 0.0);

                if (ll1 > dl2) {
                    dl1 = dl2;
                    dl2 = ll1;
                }
            }


            if ((iter == MAX_ITER) && (fabs(dl2 - dl1)/dl1) >= eps) {
                puts("Too many iterations in 1st shock width!");
                exit(EXIT_FAILURE);
            }

            dl = (dl1 + dl2) / 2.0;  // Make final refimenent.
            d1l[i] = dl;

        }

// Width of the second shock. I skip comments as they are the same as for the 1st shock.

        rho20 = mdot2 / (PI_4 * d2s * vel2);
        e0 = 1.957441607 * mu2 * sq(v2n[i]);
        coeff = 8.3245e3 / (d0 * sq(mu_av2)) * rho20 * cube(v2n[i]);

        eps = 1.0e-5;
        ll2 = 0.0;
        h1 = e0/300.0;
        e_start = 0.0;
        status4 = odeint(&ll2, &e_start, e0, dldE2, &h1, eps, 0.0);
        dl2 = ll2;
//        printf("2 dl2= %e\n",dl2);

        if (fast) {
            d2l[i] = dl2;
        } else {
            dl2 = MIN(dl2, 0.9 * (d-x10-rstar2));
            eps = 1.0e-5;
            iter = 0;
            dl1 = 0.0;
//            while ((fabs((dl2-dl1)/dl2) > eps) && (iter <= MAX_ITER)) {
            while ((fabs((dl2-dl1)/dl2) > eps1) && (iter < MAX_ITER)) {
                iter++;
                dl2 = (dl1 + dl2) / 2.0;
                xxc = d - (xc[i]+dl2*sphi[i]);
                yyc = yc[i]-dl2*cphi[i];
                d2 = dis(xxc, yyc);
                cosb = xxc/d2;
                sinb = yyc/d2;
                vel2 = v2(d2);
                v2n0 = vel2*( cosb*sphi[i]+sinb*cphi[i] );
                rho20 = mdot2/PI_4/sq(d2)/vel2;
                e0 = 1.957441607*mu2*sq(v2n0);
                coeff = 8.3245e3/d0/sq(mu_av2)*rho20*cube(v2n0);
                h1 = e0/300.0;
                ll2 = 0.0;
                e_start = 0.0;
//                printf( "iter2= %d dl1= %le dl2= %le\n", iter, dl1, dl2 );
                status4 = odeint(&ll2, &e_start, e0, dldE2, &h1, eps, 0.0);

                if (ll2 > dl2) {
                    dl1 = dl2;
                    dl2 = ll2;
                }
            }

            if ((iter == MAX_ITER) && (fabs(dl2 - dl1)/dl1) >= eps1) {
                puts("Too many iterations in 2nd shock width!");
                exit(EXIT_FAILURE);
            }

            dl = (dl1 + dl2) / 2.0;
            d2l[i] = dl;

        }  // End refining the width of 2nd shock.

    } // End of loop on knots of grid on the c.s.
    return;
}

double lhs_20(double y, int nstar){
    /* Use spline interpolation to get x at given y */
    double x = splint(yc, xc, d2xdy2_c, N, y);

    /* Use spline interpolation to get dx/dy at given y*/
    double dxdy;

    double cphi, sphi, mdot, r, cosb, sinb, vel, vn, rho;

    dxdy = splint(yc, dxdy_c, d2xcdy, N, y);

    if( fabs(dxdy_c[0]) < DBL_MIN ) cphi = 0.0;
    else cphi = dxdy_c[0]/sqrt(sq(dxdy_c[0])+1.0);
    sphi = 1.0/sqrt(sq(dxdy_c[0])+1.0);
    if (nstar == 1) {
        mdot = mdot1;
        r = dis(x,y);
        cosb = x/r;
        sinb = y/r;
        vel = v1(r);
        vn = cosb*sphi-sinb*cphi;
        vn = vel*vn;
    } else {
        mdot = mdot2;
        r = dis(x,y);
        cosb = (d-x)/r;
        sinb = y/r;
        vel = v2(r);
        vn = cosb*sphi+sinb*cphi;
        vn = vel*vn;
    }

    rho = mdot/(PI_4*sq(r)*vel);
    return rho*vn*PI_2*y/sphi;
}

double make_3d(){
    int i,j;
    double delta_theta = PI_2/(nc-1);
    //alpha[0] = 0.0;

    for(int k = 0; k < nc; k++){
        theta_3d[k] = k * delta_theta;
        //printf("%f\n",cos(angle[k]));
    }
    //printf("i\tj\t x\t y\t z\t cos(a)\n");

    for(i = 0; i < nc; i++){
        for(j = 0; j< nc; j++){
            x_3d_2[i][j] = (xc[i]-x10);
            y_3d_2[i][j] = yc[i]* cos(theta_3d[j]);
            z_3d_2[i][j] = yc[i]* sin(theta_3d[j]);
            /*
            if(z_3d_2[i][j]<= 1.e-14){
                z_3d_2[i][j] = 0.0;
            }
            */
        //printf("%i %i %13.8le %13.8le %13.8le, %13.8le\n",i, j, x_3d[i][j], y_3d[i][j],z_3d[i][j], cos(angle[j]));
        }
    }
    return 0;

}

double shift(){
    int i,j;
    for(i = 0; i < nc; i++){
        for(j = 0; j< nc; j++){
            x_3d_1[i][j] = x_3d_2[i][j]*cos(psi) - y_3d_2[i][j]*sin(psi);
        
            y_3d_1[i][j] = x_3d_2[i][j]*sin(psi) + y_3d_2[i][j]*cos(psi);
           
            
            z_3d_1[i][j] = z_3d_2[i][j];
        }

        
        

            
    }
    //printf("vt = %e \t cos(vt)= %e", (*at),cos(*at));
    return 0;
}


double skew(){
    int i,j;
    double psi;
    psi = psi_deg * PI_2 / 360;
    for(i = 0; i < nc; i++){
        for(j = 0; j< nc; j++){
            x_3d[i][j] = ((x_3d_2[i][j])*cos(psi) - y_3d_2[i][j]*sin(psi)) +x10;
        
            y_3d[i][j] = (x_3d_2[i][j])*sin(psi) + y_3d_2[i][j]*cos(psi);
           
            
            z_3d[i][j] = z_3d_2[i][j];
        }
            
    }
    //printf("vt = %e \t cos(vt)= %e", (*at),cos(*at));
    return 0;
}

double add_x10(double x){
    int i,j;
    for(i = 0; i < nc; i++){
        for(j = 0; j< nc; j++){
           
            x_3d[i][j] = x_3d[i][j] + x;


        }
    }
    return 0;
}

double sub_x10(double x){
    int i,j;
    for(i = 0; i < nc; i++){
        for(j = 0; j< nc; j++){
           
            x_3d[i][j] = x_3d[i][j] - x;
    


        }
    }
    return 0;
}
/*  
double calc_velocities(){
    int i,j;

    double vel1, vel2;
    double sinb_3d, cosb_3d;
    double d1, d2, d1s, d2s;


    
    for(i = 0; i < nc; i++){
        for(j = 0; j< nc; j++){

            d1 = sqrt(sq(x_3d_s[i][j]) + sq(y_3d_s[i][j]) + sq(z_3d_s[i][j]));  // Distance from the point on c.s. to center of 1st star.
            d1s = sq(d1);
            //printf("%e\n", d1);
            vel1 = v1(d1); // Velocity of 1st wind at point of c.s.
            cosb_3d = x_3d_s[i][j] / d1;
            sinb_3d = y_3d_s[i][j] / d1;
            //printf("%e\n", cosb_3d);
            v1t_3d[i][j] = vel1 * ( cosb_3d * cphi[i] + sinb_3d * sphi[i] );
            //Tangential and normal velocities.
            v1n_3d[i][j] = vel1 * ( cosb_3d * sphi[i] - sinb_3d * cphi[i] );
        }
    }
     for(i = 0; i < nc; i++){
        for(j = 0; j< nc; j++){
            d2 = sqrt(sq(d-x_3d_s[i][j]) + sq(y_3d_s[i][j]) + sq(z_3d_s[i][j]));  // Distance from the point on c.s. to center of 2st star.
            d2s = sq(d2);
            vel2 = v2(d2); // Velocity of 2st wind at point of c.s.
            cosb_3d = (d-x_3d_s[i][j]) / d1;
            sinb_3d = y_3d_s[i][j] / d1;
            v2t_3d[i][j] = vel2 * (- cosb_3d * cphi[i] + sinb_3d * sphi[i] );
            //Tangential and normal velocities.
            v2n_3d[i][j] = vel2 * ( cosb_3d * sphi[i] + sinb_3d * cphi[i] );

        }
     }
     return 0;

}
*/
double calc_velocities(){
    int i,j;

    double vel1, vel2;
    double cos_x1, cos_y1, cos_z1, cos_x2, cos_y2, cos_z2;
    double d1, d2, d1s, d2s, d1_xyz;
    double v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z;
    double psi;
    psi = psi_deg * (PI_2/360);
    double x_star2, y_star2;
    
    FILE *test_file;
    test_file = fopen("test.dat", "w");
    fprintf(test_file,"x\t\ty\t\tz\t\tn_x\t\tn_y\t\tn_z\t\tv1_x\t\tv1_y\t\tv1_z\t\tcos_x1\t\tcos_y1\t\tcos_z1\n");
/**********explanation************
           these are the coordinates in x,y,z; unshifted and unskewed frame of refrence.

           
                  | -cos(psi)sin(phi) - sin(psi)cos(phi)cos(theta) |   unit normal vector for
           n_xyz =|  -sin(psi)sin(phi) + cos(psi)cos(phi)cos(theta) |  ~the transition from 
                  |              cos(phi)sin(theta)                |   x2y2z2 frame of refrence to xyz

                    | |v1|*x/r | 
           v1_xyz = | |v1|*y/r |    ~projection of the velocity vector on each axis
                    | |v1|*z/r |
             
           v1n = (v1_xyz . n_xyz)   ~normal component of velocity as a dot product

           v1t = sqrt( |v1|^2 - v1n^2 )   ~tangential component of velocity

           did not resort to matrix-vector multiplication funcitons and dot product function.
           resorted to saving each coordinate system in seperate arrays.
           i did this because the amount of funcitons that i would have needed for all the
           procedures without using "malloc" dynamic memory is large and the difficulty and slowing 
           effects far exceed the memory waste.
*************************************/ 
    
    for(i = 0; i < nc; i++){
        for(j = 0; j< nc; j++){

            d1 = sqrt(sq(x_3d[i][j]) + sq(y_3d[i][j])+ sq(z_3d[i][j]) );  // Distance from the point on c.s. to center of 1st star.
            d1s = sq(d1);
            
            //printf("%e\n", d1);
            vel1 = v1(d1); // Velocity of 1st wind at point of c.s.
            
            
            
           
            
            cos_x1 = (x_3d[i][j] / d1);
            cos_y1 = (y_3d[i][j] / d1);
            cos_z1 = (z_3d[i][j] / d1);

            v1_x = (vel1 * cos_x1);
            v1_y = (vel1 * cos_y1);
            v1_z = (vel1 * cos_z1);
            /*
            n1_x = (-cos(psi)*sin(PI/4) - sin(psi)*cos(PI/4)*cos(theta_3d[j]));
            n1_y = (-sin(psi)*sin(PI/4) + cos(psi)*cos(PI/4)*cos(theta_3d[j]));
            n1_z = (cos(PI/4) * sin(theta_3d[j]));
            */
            n1_x = (-cos(psi)*sphi[i] - sin(psi)*cphi[i]*cos(theta_3d[j]));
            n1_y = (-sin(psi)*sphi[i] + cos(psi)*cphi[i]*cos(theta_3d[j]));
            n1_z = (cphi[i] * sin(theta_3d[j]));
            
            
            
            fprintf(test_file,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",x_3d[i][j],y_3d[i][j],z_3d[i][j], 
                        n1_x,n1_y,n1_z,v1_x,v1_y,v1_z,cos_x1,cos_y1,cos_z1);
            //printf("%e %e %e \n", v_x1,v_y1,v_z1);
            
/*
             
                  | -cos(psi)sin(phi) - sin(psi)cos(phi)cos(theta) |   unit normal vector for
           n_xyz =|  -sin(psi)sin(phi) + cos(psi)cos(phi)cos(theta) |  ~the transition from 
                  |              cos(phi)sin(theta)                |   x2y2z2 frame of refrence to xyz

                    | |v1|*x/r | 
           v1_xyz = | |v1|*y/r |    ~projection of the velocity vector on each axis
                    | |v1|*z/r |
             
           v1n = -(v1_xyz . n_xyz)   ~normal component of velocity as a dot product

           */
            v1n_3d[i][j] =  -( (v1_x * n1_x) + (v1_y * n1_y) + (v1_z * n1_z) ) ;
            v1t_3d[i][j] = sqrt(sq(fabs(vel1)) - sq(v1n_3d[i][j]));
            //Tangential and normal velocities.
            
    
            
            //printf("%e\n", cosb_3d);
            
        }
    }
    y_star2 =  d*sin(vt);
    x_star2 =  d*cos(vt);
     for(i = 0; i < nc; i++){
        for(j = 0; j< nc; j++){


            

            d2 = sqrt(sq(x_3d[i][j]-d) + sq(y_3d[i][j])+ sq(z_3d[i][j]) );  // Distance from the point on c.s. to center of 1st star.
            d2s = sq(d1);
            //printf("%e\n", d1);
            vel2 = v2(d2); // Velocity of 1st wind at point of c.s.
            
            
            cos_x2 = ((x_3d[i][j]-d) / d2);
            cos_y2 = (y_3d[i][j]/ d2);
            cos_z2 = (z_3d[i][j] / d2);

            v2_x = (vel2 * cos_x2);
            v2_y = (vel2 * cos_y2);
            v2_z = (vel2 * cos_z2);

            n2_x = (-cos(psi)*sphi[i] - sin(psi)*cphi[i]*cos(theta_3d[j]));
            n2_y = (-sin(psi)*sphi[i] + cos(psi)*cphi[i]*cos(theta_3d[j]));
            n2_z = (cphi[i] * sin(theta_3d[j]));

   /*           
                  | -cos(psi)sin(phi) + sin(psi)cos(phi)cos(theta) |   unit normal vector for
           n_xyz =|  sin(psi)sin(phi) + cos(psi)cos(phi)cos(theta) |  ~the transition from 
                  |              cos(phi)sin(theta)                |   x2y2z2 frame of refrence to xyz

                    | |v2|*x/r | 
           v2_xyz = | |v2|*y/r |    ~projection of the velocity vector on each axis
                    | |v2|*z/r |
             
           v2n = (v1_xyz . n_xyz)   ~normal component of velocity as a dot product
*/
            
            v2n_3d[i][j] =  ( (v2_x * n2_x) + (v2_y * n2_y) + (v2_z * n2_z) ) ;
            v2t_3d[i][j] = sqrt(sq(fabs(vel2)) - sq(v2n_3d[i][j]));
            //Tangential and normal velocities.
            
    
            
            //printf("%e\n", cosb_3d);
            
        }
    }
     return 0;

}

double calc_sigma1(){
    double eps, h1, dd1, dd2, tmp;
    double rho1, rho2;
    double r_xy, m;
    int i,j,k,cc,bb,loop;
    int status3d;
    double r1, r2;
    eps = 1.0e-6;
    double dummy1, dummy2;
    double a,b,c;
    double vel1;
    double vel2;
    double psi;
    double y1,y2;
    int index;
    double mn = v1t_3d[0][nc/2];
    double dum_sig1[NE0_MAX], dum_sig2[NE0_MAX];
   double velocity_3d[NE0_MAX];
    
    psi = psi_deg * PI_2 / 360;
    
    vel1 = v1(x10);
    vel2 = v2(d - x10);
    eta0 = mdot1 * vel1 / (mdot2 * vel2);
    bb = 0;
    for(int cc=0; cc<nc; cc++){
        
        if(v1t_3d[cc][nc/2]<mn){
            mn = v1t_3d[cc][nc/2];
            stag_index = cc;
        }
        
    }
   

    m = y_3d[stag_index][nc/2]/ x_3d[stag_index][nc/2];
   
    /* this is where I concatenate arrays*/
    index = 0;
    for(loop = stag_index; loop > 0; loop--) {
        dis_to_x_stag_top[index] = dis_to_x_stag[loop];
        yc_top[index] = y_3d_2[loop][(nc-1)/2];
        v1t_3d_top[index] = v1t_3d[loop][(nc-1)/2];
        x_3d_top[index] = x_3d[loop][(nc-1)/2];
        y_3d_top[index] = y_3d[loop][(nc-1)/2];
        z_3d_top[index] = z_3d[loop][(nc-1)/2];
        r_3d_top[index] = r_3d[loop][(nc-1)/2];
        v1n_3d_top[index] = v1n_3d[loop][(nc-1)/2];
        dxdy_3d_top[index] = dxdy_c[loop];
        d2xdy2_3d[index] = d2xdy2_c[loop];
       //velocity_3d[index] = v1()
        
        index++;
    }
     
      
   for(loop = 0; loop < nc; loop++) {
        dis_to_x_stag_top[index] = dis_to_x_stag[loop];
        yc_top[index] = y_3d[loop][0];
        v1t_3d_top[index] = v1t_3d[loop][0];
        x_3d_top[index] = x_3d[loop][0];
        y_3d_top[index] = y_3d[loop][0];
        z_3d_top[index] = z_3d[loop][0];
        r_3d_top[index] = r_3d[loop][0];
        v1n_3d_top[index] = v1n_3d[loop][0];
        dxdy_3d_top[index] = dxdy_c[loop];
        d2xdy2_3d[index] = d2xdy2_c[loop];
        index++;
   }
   rho1 = mdot1 / (PI_4 * sq(r_3d_top[0]) * v1(r_3d_top[0]));
   //printf("v1t_stag = %e, index = %i\n",v1t_3d[stag_index][nc/2], stag_index );
    sig1_3d_top[0]  = rho1 * r_3d_top[0]/ (2 * (1 +  r_3d_top[0]* d2xdy2_c[stag_index]));
   //printf("sig1_3d_stag = %e, %e\n",r_3d[stag_index][nc/2],d2xdy2_c[stag_index] );
   
   printf("stag_index_1 = %i slope = %e sig1_stag = %e d2xdy2 = %e\n",stag_index, m, sig1_3d_top[0], d2xdy2_c[stag_index]);

    

    if (fabs(eta0 - 1) <= 1.0e-5) eta0 = 1.0;
    if( fabs(eta0-1.0) < DBL_MIN ) {
      d2xdy20 = 0.0;
    } else {
      a = d*sqrt(eta0)/2.0/sq( 1.0+sqrt(eta0) );
      b = dvdy(x10, 0.0, 1);
      c = dvdy(x10, 0.0, 2);

      dummy1 = 2.0*(eta0-1.0)/d/sqrt(eta0) + a*(b-c);
      dummy2 = 3.0 - a*(b*x10+c*(d-x10));
      d2xdy20 = dummy1/dummy2;
    }

    for(int i = 0; i<index; i++){
        dis_to_x_stag_top[i] = fabs(-m* x_3d_top[i] + y_3d_top[i] )/sqrt(sq(m)+1);
        //printf("%e\n",yc_top[i]);
    }
  
    
    for(k = 1; k<index; k++){
        dum_sig1[0] = dummy1;
        kkk = k;
        y1 = dis_to_x_stag_top[k-1];
        y2 = dis_to_x_stag_top[k];
        //printf("y1 = %e, y2 = %e\n", y1,y2);
        dum_sig1[k] = quad(&lhs_20_3d_1_xy, y1, y2, 0.0, 0.000001, 1000,&Ierror);
        //dum_sig1[k] = qsimp(&lhs_20_3d_1_xy, y1, y2);  
        
        for(i = 1; i<k+1; i++){
            sig1_3d_top[k] += dum_sig1[i]; 
            
        }
        
        sig1_3d_top[k] = sig1_3d_top[k]/(v1t_3d_top[k]*dis_to_x_stag_top[k]);
        
    }
    //printf("%i\n", index);
/*   
    for(k = 0; k<index; k++){
        printf("sig3d= %e v1t = %e dis_to_stag = %e z = %e\n",
        sig1_3d_top[k],v1t_3d_top[k] ,dis_to_x_stag_top[k],z_3d_top[k]);
    }
*/  

/////////////////////////////////////////////////////////

    index = 0;
    for(loop = stag_index; loop < nc; loop++) {
        dis_to_x_stag_bot[index] = dis_to_x_stag[loop];
        yc_bot[index] = y_3d[loop][(nc-1)/2];
        v1t_3d_bot[index] = v1t_3d[loop][(nc-1)/2];
        x_3d_bot[index] = x_3d[loop][(nc-1)/2];
        y_3d_bot[index] = y_3d[loop][(nc-1)/2];
        z_3d_bot[index] = z_3d[loop][(nc-1)/2];
        r_3d_bot[index] = r_3d[loop][(nc-1)/2];
        v1n_3d_bot[index] = v1n_3d[loop][(nc-1)/2];
        dxdy_3d_bot[index] = dxdy_c[loop];
        d2xdy2_3d[index] = d2xdy2_c[loop];
       //velocity_3d[index] = v1()
        
        index++;
    }
     
      
   rho1 = mdot1 / (PI_4 * sq(r_3d_bot[0]) * v1(r_3d_bot[0]));
   //printf("v1t_stag = %e, index = %i\n",v1t_3d[stag_index][nc/2], stag_index );
    sig1_3d_bot[0]  = rho1 * r_3d_bot[0]/ (2 * (1 +  r_3d_bot[0]* d2xdy2_c[stag_index]));
   ///printf("sig1_3d_stag = %e, %e\n",r_3d[stag_index][nc/2],d2xdy2_c[stag_index] );
   

    

    if (fabs(eta0 - 1) <= 1.0e-5) eta0 = 1.0;
    if( fabs(eta0-1.0) < DBL_MIN ) {
      d2xdy20 = 0.0;
    } else {
      a = d*sqrt(eta0)/2.0/sq( 1.0+sqrt(eta0) );
      b = dvdy(x10, 0.0, 1);
      c = dvdy(x10, 0.0, 2);

      dummy1 = 2.0*(eta0-1.0)/d/sqrt(eta0) + a*(b-c);
      dummy2 = 3.0 - a*(b*x10+c*(d-x10));
      d2xdy20 = dummy1/dummy2;
    }

    for(int i = 0; i<index; i++){
        dis_to_x_stag_bot[i] = fabs(-m* x_3d_bot[i] + y_3d_bot[i] )/sqrt(sq(m)+1);
        //printf("%e\n",yc_bot[i]);
    }
 
    for(k = 1; k<index; k++){
        dum_sig1[0] = dummy1;
        kkk = k;
        y1 = fabs(dis_to_x_stag_bot[k-1]);
        y2 = fabs(dis_to_x_stag_bot[k]);
        //printf("y1 = %e, y2 = %e\n", y1,y2);
        dum_sig1[k] = quad(&lhs_20_3d_1_xy_bot, y1, y2, 0.0, 0.000001, 1000,&Ierror);
        //dum_sig1[k] = qsimp(&lhs_20_3d_1_xy, y1, y2);  
        
        for(i = 1; i<k+1; i++){
            sig1_3d_bot[k] += dum_sig1[i]; 
            
        }
        
        sig1_3d_bot[k] = sig1_3d_bot[k]/(v1t_3d_bot[k]*dis_to_x_stag_bot[k]);
        
    }
    
    //printf("///////////////\n");
/*
    for(k = 0; k<index; k++){
            printf("sig3d= %e sig2d = %e delta(sig) = %e v1t = %e dxdy3d = %e dxdy = %e dis_to_stag = %e z = %e\n",
            sig1_3d_bot[k],sig1[k],fabs(sig1_3d_bot[k]-sig1[k])/sig1_3d_bot[k],v1t_3d_bot[k],dxdy_3d_bot[k],dxdy_c[k] ,dis_to_x_stag_bot[k],z_3d_bot[k]);
        }
*/    
    for(i = 0; i<nc; i++){
        sig1_3d[i][0] = sig1_3d_top[i+stag_index];
    }
    for(i = 0; i < stag_index; i++){
        sig1_3d[i][(nc-1)/2] = sig1_3d_top[stag_index-i];
    }
    for(i = 0; i < nc; i++){
        sig1_3d[i+stag_index][(nc-1)/2] = sig1_3d_bot[i];
    } 

    


    for(i = 0; i<nc; i++){
        for(j = 0; j<nc; j++){
            //sig1(i,j) = [sig1(n,j)-sig1(0,j)]/2*[cos(0)-cos(theta(i)] + sig1(0,j)
            sig1_3d[j][i] = (sig1_3d[j][(nc-1)/2] - sig1_3d[j][0])/2*(cos(theta_3d[0])-cos(theta_3d[i])) +sig1_3d[j][0]; 
        }
    }
/*
    for(k = 0; k<nc; k++){
            printf("%e  %e \n" ,sig1_3d[k][0],sig1_3d[k][(nc-1)/2]);
        }

*/
        return 0;

}


double calc_sigma2(){
    double eps, h1, dd1, dd2, tmp;
    double rho1, rho2;
    double r_xy, m;
    int i,j,k,cc,bb,loop;
    int status3d;
    double r1, r2;
    eps = 1.0e-6;
    double dummy1, dummy2;
    double a,b,c;
    double vel1;
    double vel2;
    double psi;
    double y1,y2;
    int index;
    double mn = v2t_3d[0][0];
    double dum_sig1[NE0_MAX], dum_sig2[NE0_MAX];
    //double velocity_3d[NE0_MAX];
    
    psi = psi_deg * PI_2 / 360;
    
    vel1 = v1(x10);
    vel2 = v2(d - x10);
    eta0 = mdot1 * vel1 / (mdot2 * vel2);
    bb = 0;
    cc = 0;
    for(int cc=0; cc<nc; cc++){
        
        if(v2t_3d[cc][0]<mn){
            mn = v2t_3d[cc][0];
            stag_index_2 = cc;
        }
        
    }
    //stag_index_2 = stag_index_2 ;
   
    m = y_3d[stag_index_2][0]/ (x_3d[stag_index_2][0]-d);
    
    /* this is where I concatenate arrays*/
    index = 0;
    for(loop = stag_index_2; loop > 0; loop--) {
        dis_to_x_stag_bot_2[index] = dis_to_x_stag_2[loop];
        yc_bot_2[index] = y_3d_2[loop][0];
        v2t_3d_bot[index] = v2t_3d[loop][0];
        x_3d_bot_2[index] = x_3d[loop][0];
        y_3d_bot_2[index] = y_3d[loop][0];
        z_3d_bot_2[index] = z_3d[loop][0];
        r_3d_bot_2[index] = r_3d_2[loop][0];
        v2n_3d_bot[index] = v2n_3d[loop][0];
        dxdy_3d_bot_2[index] = dxdy_c[loop];
        d2xdy2_3d_2[index] = d2xdy2_c[loop];
       //velocity_3d[index] = v1()
        
        index++;
    }
     
      
   for(loop = 0; loop < nc; loop++) {
        dis_to_x_stag_bot_2[index] = dis_to_x_stag_2[loop];
        yc_bot_2[index] = y_3d_2[loop][(nc-1)/2];
        v2t_3d_bot[index] = v2t_3d[loop][(nc-1)/2];
        x_3d_bot_2[index] = x_3d[loop][(nc-1)/2];
        y_3d_bot_2[index] = y_3d[loop][(nc-1)/2];
        z_3d_bot_2[index] = z_3d[loop][(nc-1)/2];
        r_3d_bot_2[index] = r_3d_2[loop][(nc-1)/2];
        v2n_3d_bot[index] = v2n_3d[loop][(nc-1)/2];
        dxdy_3d_bot_2[index] = dxdy_c[loop];
        d2xdy2_3d_2[index] = d2xdy2_c[loop];
        index++;
   }
   rho2 = mdot2 / (PI_4 * sq(r_3d_bot_2[0]) * v2(r_3d_bot_2[0]));
   //printf("v1t_stag = %e, index = %i\n",v1t_3d[stag_index][nc/2], stag_index );
    sig2_3d_bot[0]  = rho2 * r_3d_bot_2[0]/ (2 * (1-r_3d_bot_2[0]* d2xdy2_c[stag_index]));
   //printf("sig1_3d_stag = %e, %e\n",r_3d_2[stag_index][0],d2xdy2_c[stag_index_2] );
   printf("stag_index_2 = %i slope = %e sig2_stag = %e d2xdy2 = %e\n",stag_index_2, m, sig2_3d_bot[0], d2xdy2_c[stag_index]);
   
    if (fabs(eta0 - 1) <= 1.0e-5) eta0 = 1.0;
    if( fabs(eta0-1.0) < DBL_MIN ) {
      d2xdy20 = 0.0;
    } else {
      a = d*sqrt(eta0)/2.0/sq( 1.0+sqrt(eta0) );
      b = dvdy(x10, 0.0, 1);
      c = dvdy(x10, 0.0, 2);

      dummy1 = 2.0*(eta0-1.0)/d/sqrt(eta0) + a*(b-c);
      dummy2 = 3.0 - a*(b*x10+c*(d-x10));
      d2xdy20 = dummy1/dummy2;
    }
    

    for(int i = 0; i<index; i++){
        dis_to_x_stag_bot_2[i] = fabs(-m* (x_3d_bot_2[i]-d) + y_3d_bot_2[i] )/sqrt(sq(m)+1);
        //printf("%e\n",yc_top[i]);
    }
    
    
    for(k = 1; k<index; k++){
        //dum_sig2[0] = dummy1;
        kkk_2 = k;
        y1 = dis_to_x_stag_bot_2[k-1];
        y2 = dis_to_x_stag_bot_2[k];
        //printf("y1 = %e, y2 = %e\n", y1,y2);
        dum_sig2[k] = quad(&lhs_20_3d_2_xy_bot, y1, y2, 0.0, 0.000001, 1000,&Ierror);
        //dum_sig1[k] = qsimp(&lhs_20_3d_2_xy_bot, y1, y2);  
        
        for(i = 1; i<k+1; i++){
            sig2_3d_bot[k] += dum_sig2[i]; 
            
        }
        
        sig2_3d_bot[k] = sig2_3d_bot[k]/(v2t_3d_bot[k]*dis_to_x_stag_bot_2[k]);
        
    }
    /*
    for(k = stag_index_2+1; k<index; k++){
        //dum_sig2[0] = dummy1;
        kkk_2 = k;
        y1 = dis_to_x_stag_bot_2[k-1];
        y2 = dis_to_x_stag_bot_2[k];
        //printf("y1 = %e, y2 = %e\n", y1,y2);
        dum_sig2[k] = quad(&lhs_20_3d_2_xy_bot, y1, y2, 0.0, 0.000001, 1000,&Ierror);
        //dum_sig1[k] = qsimp(&lhs_20_3d_2_xy_bot, y1, y2);  
        
        for(i = 1; i<k+1; i++){
            sig2_3d_bot[k] += dum_sig2[i]; 
            
        }
        
        sig2_3d_bot[k] = sig2_3d_bot[k]/(v2t_3d_bot[k]*dis_to_x_stag_bot_2[k]);
        
    }
    */
   
    //printf("%i\n", index);
/*  
    for(k = 0; k<index; k++){
        printf("sig3d= %e v1t = %e dis_to_stag = %e z = %e\n",
        sig2_3d_bot[k],v2t_3d_bot[k] ,d2xdy2_3d_2[k],z_3d_bot_2[k]);
    }
  
*/
/////////////////////////////////////////////////////////

    index = 0;
    for(loop = stag_index_2; loop < nc; loop++) {
        dis_to_x_stag_top_2[index] = dis_to_x_stag_2[loop];
        yc_top_2[index] = y_3d_2[loop][0];
        v2t_3d_top[index] = v2t_3d[loop][0];
        x_3d_top_2[index] = x_3d[loop][0];
        y_3d_top_2[index] = y_3d[loop][0];
        z_3d_top_2[index] = z_3d[loop][0];
        r_3d_top_2[index] = r_3d_2[loop][0];
        v2n_3d_top[index] = v2n_3d[loop][0];
        dxdy_3d_top_2[index] = dxdy_c[loop];
        d2xdy2_3d_2[index] = d2xdy2_c[loop];
       //velocity_3d[index] = v1()
        
        index++;
    }

     
      
    rho2 = mdot2 / (PI_4 * sq(r_3d_top_2[0]) * v2(r_3d_top_2[0]));
   //printf("v1t_stag = %e, index = %i\n",v1t_3d[stag_index][nc/2], stag_index );
    sig2_3d_top[0]  = rho2 * r_3d_top_2[0]/ (2 * (1-r_3d_top_2[0]* d2xdy2_c[stag_index]));
    ///printf("sig1_3d_stag = %e %e %e %e\n",r_3d_top_2[0], r_3d_top[0], d2xdy2_c[stag_index_2], d2xdy2_3d[0]);
   



    for(int i = 0; i<index; i++){
        dis_to_x_stag_top_2[i] = fabs(-m*(x_3d_top_2[i]-d) + y_3d_top_2[i] )/sqrt(sq(m)+1);
        //printf("%e\n",yc_top[i]);

    }
  
    
    for(k = 1; k<index; k++){
        //dum_sig2[0] = sig2_3d_top[0];
        kkk_2 = k;
        y1 = fabs(dis_to_x_stag_top_2[k-1]);
        y2 = fabs(dis_to_x_stag_top_2[k]);
        //printf("y1 = %e, y2 = %e\n", y1,y2);
        dum_sig2[k] = quad(&lhs_20_3d_2_xy, y1, y2, 0.0, 0.000001, 1000,&Ierror);
        //dum_sig1[k] = qsimp(&lhs_20_3d_2_xy, y1, y2);  
        
        for(i = 1; i<k+1; i++){
            sig2_3d_top[k] += dum_sig2[i]; 
            
        }
        
        sig2_3d_top[k] = sig2_3d_top[k]/(v2t_3d_top[k]*dis_to_x_stag_top_2[k]);
        
    }
    
///    printf("///////////////\n");
/*
    for(k = 0; k<index; k++){
            printf("sig3d= %e  v1t = %e  dis_to_stag = %e z = %e\n",
            sig2_3d_top[k],v2t_3d_top[k] ,dis_to_x_stag_top_2[k],z_3d_top_2[k]);
        }
*/    
    for(i = 0; i<nc; i++){
        sig2_3d[i][(nc-1)/2] = sig2_3d_bot[i+stag_index_2];
    }
    for(i = 0; i < stag_index_2; i++){
        sig2_3d[i][0] = sig2_3d_bot[stag_index_2-i];
    }
    for(i = 0; i < nc; i++){
        sig2_3d[i+stag_index_2][0] = sig2_3d_top[i];
    } 
/*
    for(j = 1; j<nc; j++){
            //sig1(i,j) = [sig1(n,j)-sig1(0,j)]/2*[cos(0)-cos(theta(i)] + sig1(0,j)
            sig2_3d[j][0] = (sig2_3d[j][0] - sig2_3d[j-1][0])/2*(y_3d[j][0]-y_3d[j-1][0]) + y_3d[j][0];
            sig2_3d[j][(nc-1)/2] = (sig2_3d[j][(nc-1)/2] - sig2_3d[j-1][(nc-1)/2])/2*(y_3d[j][(nc-1)/2]-y_3d[j-1][(nc-1)/2]) + y_3d[j][(nc-1)/2]; 
        }
    
*/
 //   sig2_3d[stag_index_2][0] = (sig2_3d[stag_index_2][0]+sig2_3d[stag_index_2-1][0])/2;
    
    for(i = 0; i<nc; i++){
        for(j = 0; j<nc; j++){
            //sig1_3d[j][i] = (sig1_3d[j][(nc-1)/2] - sig1_3d[j][0])/2*(cos(theta_3d[0])-cos(theta_3d[i])) +sig1_3d[j][0]; 
            sig2_3d[j][i] = (sig2_3d[j][(nc-1)/2] - sig2_3d[j][0])/2*(cos(theta_3d[0])-cos(theta_3d[i])) +sig2_3d[j][0]; 
        }
    }
    

    
    //printf("%e %e %e %e\n", sig2_3d_top[0] ,r_3d_bot_2[0], r_3d_top_2[0],  dxdy_3d_bot_2[0]);
    
/*    for(k = 0; k<nc; k++){
            printf("%e\n" ,sig2_3d[k][0] - sig2_3d[k][(nc-1)/2]);
        }
*/
    return 0;

}




double get_r(){
    
    

    for(int i = 0; i < nc; i++){
        
        
        for(int j = 0; j< nc; j++){
            
            r_3d[i][j] = sqrt( sq(x_3d[i][j]) + sq(y_3d[i][j]) + sq(z_3d[i][j]) );
            r_3d_2[i][j] = sqrt( sq(x_3d[i][j]-d) + sq(y_3d[i][j]) + sq(z_3d[i][j]) );
            r_3d_3[i][j] = sqrt( sq(x_3d_2[i][j]) + sq(y_3d_2[i][j]) + sq(z_3d_2[i][j]) );
            
            //r_3d[i][j] = y_3d[i][j];

        }
    }

    return 0;


}



double lhs_20_3d_1_xy(double y){
    
    int i,j,index,loop;
    

    double x, dxdy; 

    double cphi, sphi, r, vel, vn, vnn, rho;
  
   

   
    if( fabs(dxdy_3d_top[kkk]) < DBL_MIN ) cphi = 0.0;
    else cphi = dxdy_3d_top[kkk]/sqrt(sq(dxdy_3d_top[kkk])+1.0);
    sphi = 1.0/sqrt(sq(dxdy_3d_top[kkk])+1.0);

    ///yc_top[kkk] = (yc_top[kkk]-yc_top[kkk-1])/(dis_to_x_stag_top[kkk]-dis_to_x_stag_top[kkk-1])*(y-dis_to_x_stag_top[kkk-1]) + yc_top[kkk-1];
    r = (r_3d_top[kkk]-r_3d_top[kkk-1])/(dis_to_x_stag_top[kkk]-dis_to_x_stag_top[kkk-1])*(y-dis_to_x_stag_top[kkk-1]) + r_3d_top[kkk-1];
    vel = v1(r);
    vnn =(v1n_3d_top[kkk]-v1n_3d_top[kkk-1])/(dis_to_x_stag_top[kkk]-dis_to_x_stag_top[kkk-1])*(y-dis_to_x_stag_top[kkk-1]) + v1n_3d_top[kkk-1];

    vn = vel * vnn;


    rho = mdot1/(PI_4*sq(r)*vel);

    return rho*vnn*y/sphi; 
   
    
}

double lhs_20_3d_1_xy_bot(double y){
    
    int i,j,index,loop;
    

    double x, dxdy; 

    double cphi, sphi, r, vel, vn, vnn, rho;
  
   

   
    if( fabs(dxdy_3d_bot[kkk]) < DBL_MIN ) cphi = 0.0;
    else cphi = dxdy_3d_bot[kkk]/sqrt(sq(dxdy_3d_bot[kkk])+1.0);
    sphi = 1.0/sqrt(sq(dxdy_3d_bot[kkk])+1.0);

    ///yc_bot[kkk] = (yc_bot[kkk]-yc_bot[kkk-1])/(dis_to_x_stag_bot[kkk]-dis_to_x_stag_bot[kkk-1])*(y-dis_to_x_stag_bot[kkk-1]) + yc_bot[kkk-1];
    r = (r_3d_bot[kkk]-r_3d_bot[kkk-1])/(dis_to_x_stag_bot[kkk]-dis_to_x_stag_bot[kkk-1])*(y-dis_to_x_stag_bot[kkk-1]) + r_3d_bot[kkk-1];
    vel = v1(r);
    vnn =(v1n_3d_bot[kkk]-v1n_3d_bot[kkk-1])/(dis_to_x_stag_bot[kkk]-dis_to_x_stag_bot[kkk-1])*(y-dis_to_x_stag_bot[kkk-1]) + v1n_3d_bot[kkk-1];

    vn = vel * vnn;


    rho = mdot1/(PI_4*sq(r)*vel);

    return rho*vnn*y/sphi; 
   
    
}


double lhs_20_3d_2_xy(double y){
    
    int i,j,index,loop;
    

    double x, dxdy; 

    double cphi, sphi, r, vel, vn, vnn, rho;
  
   

   
    if( fabs(dxdy_3d_top_2[kkk_2]) < DBL_MIN ) cphi = 0.0;
    else cphi = dxdy_3d_top_2[kkk_2]/sqrt(sq(dxdy_3d_top_2[kkk_2])+1.0);
    sphi = 1.0/sqrt(sq(dxdy_3d_top_2[kkk_2])+1.0);

    ///yc_top[kkk] = (yc_top[kkk]-yc_top[kkk-1])/(dis_to_x_stag_top[kkk]-dis_to_x_stag_top[kkk-1])*(y-dis_to_x_stag_top[kkk-1]) + yc_top[kkk-1];
    r = (r_3d_top_2[kkk_2]-r_3d_top_2[kkk_2-1])/(dis_to_x_stag_top_2[kkk_2]-dis_to_x_stag_top_2[kkk_2-1])*(y-dis_to_x_stag_top_2[kkk_2-1]) + r_3d_top_2[kkk_2-1];
    vel = v2(r);
    vnn =(v2n_3d_top[kkk_2]-v2n_3d_top[kkk_2-1])/(dis_to_x_stag_top_2[kkk_2]-dis_to_x_stag_top_2[kkk_2-1])*(y-dis_to_x_stag_top_2[kkk_2-1]) + v2n_3d_top[kkk_2-1];
    vn = vel * vnn;


    rho = mdot2/(PI_4*sq(r)*vel);

    return rho*vnn*y/sphi; 
   
    
}

double lhs_20_3d_2_xy_bot(double y){
    
    int i,j,index,loop;
    

    double x, dxdy; 

    double cphi, sphi, r, vel, vn, vnn, rho;
  
   

   
    if( fabs(dxdy_3d_bot_2[kkk_2]) < DBL_MIN ) cphi = 0.0;
    else cphi = dxdy_3d_bot_2[kkk_2]/sqrt(sq(dxdy_3d_bot_2[kkk_2])+1.0);
    sphi = 1.0/sqrt(sq(dxdy_3d_bot_2[kkk_2])+1.0);

    ///yc_bot[kkk] = (yc_bot[kkk]-yc_bot[kkk-1])/(dis_to_x_stag_bot[kkk]-dis_to_x_stag_bot[kkk-1])*(y-dis_to_x_stag_bot[kkk-1]) + yc_bot[kkk-1];
    r = (r_3d_bot_2[kkk_2]-r_3d_bot_2[kkk_2-1])/(dis_to_x_stag_bot_2[kkk_2]-dis_to_x_stag_bot_2[kkk_2-1])*(y-dis_to_x_stag_bot_2[kkk_2-1]) + r_3d_bot_2[kkk_2-1];
    vel = v2(r);
    vnn =(v2n_3d_bot[kkk_2]-v2n_3d_bot[kkk_2-1])/(dis_to_x_stag_bot_2[kkk_2]-dis_to_x_stag_bot_2[kkk_2-1])*(y-dis_to_x_stag_bot_2[kkk_2-1]) + v2n_3d_bot[kkk_2-1];
    vn = vel * vnn;


    rho = mdot2/(PI_4*sq(r)*vel);

    return rho*vnn*y/sphi; 
   
    
}

double calc_width_1(){
    double d1, d2, rho10, e0, ll1, ll2, e_start, dl2, dl1, h1;
    int iter, status3d;
    double xxc, yyc,zzc, v1n0, dl, rho20, v2n0;
    double  eps = 1.0e-5;
    double vel1, vel2;
    double cos_x1, cos_y1, cos_z1, cos_x2, cos_y2, cos_z2;
    double d1s, d2s, d1_xyz;
    double v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z;
    double psi;
    short fast=0;
    psi = psi_deg * (PI_2/360);
    
    for(int j = 0; j<nc; j++){
        for(int i = 0; i<nc; i++){

            vel1 = v1(r_3d[i][j]);
            

        
            rho10 = mdot1 / (PI_4 * sq(r_3d[i][j]) * vel1); // Density of 1st wind at point of c.s.
            e0 = 1.957441607 * mu1 * sq(v1n_3d[i][j]); // Energy of normal component of the wind.

            coeff = 8.3245e3/d0/sq(mu_av1)*rho10*cube(v1n_3d[i][j]);

        // Compute first estimate of shock width based of energy e0 estimated by v1n0 at the c.s.
        // As the shock front is at some distance from the c.s., v1n is different so I have to refine
        // this estimate - see below.

            ll1 = 0.0;
            h1 = e0 / 300.0;
            e_start = 0.0;
            status3d = odeint(&ll1, &e_start, e0, dldE1, &h1, eps, 0.0);
            dl2 = ll1;


        /*
        Make iterations to refine the estimate of the shock width:

        1. At the current width dl2, compute velocity, energy e0, end coeff at the
        shock front position;
        2. With these, solve diff eq. to find a new width.
        3. Set the new width as average of the old and new ones.
        4. Repeat 1-3 until convergency.

        */

            if (fast) { // If fast, just take this first estimate as an answer.
                d1l_3d[i][j] = dl2;
            } else {    // Refine the estimate.
                dl2 = MIN(dl2, 0.9 * (x10 - rstar1));
                iter = 0;
                dl1 = 0.0;
                while ((fabs((dl2-dl1)/dl2) > eps) && (iter < MAX_ITER)) {
                    iter++;
                    dl2 = (dl1 + dl2) / 2.0;
                    
                

                    //printf("cphi = %e sphi = %e\n", cphi[i],sphi[i]);


                    xxc = (x_3d[i][j] - dl2 * sphi[i]);
                    yyc = (y_3d[i][j] + dl2 * cphi[i]);
                    zzc = (z_3d[i][j]);

                    d1 = sqrt(sq(xxc) + sq(yyc) + sq(zzc));

                    vel1 = v1(d1);

                    cos_x1 = xxc / d1;
                    cos_y1 = yyc / d1;
                    cos_z1 = zzc / d1;

                    v1_x = (vel1 * cos_x1);
                    v1_y = (vel1 * cos_y1);
                    v1_z = (vel1 * cos_z1);
                    /*
                    n1_x = (-cos(psi)*sin(PI/4) - sin(psi)*cos(PI/4)*cos(theta_3d[j]));
                    n1_y = (-sin(psi)*sin(PI/4) + cos(psi)*cos(PI/4)*cos(theta_3d[j]));
                    n1_z = (cos(PI/4) * sin(theta_3d[j]));
                    */
                    n1_x = (-cos(psi)*sphi[i] - sin(psi)*cphi[i]*cos(theta_3d[j]));
                    n1_y = (-sin(psi)*sphi[i] + cos(psi)*cphi[i]*cos(theta_3d[j]));
                    n1_z = (cphi[i] * sin(theta_3d[j]));
                    nor1_x[i][j] = n1_x;
                    nor1_y[i][j] = n1_y;
                    nor1_z[i][j] = n1_z;
                

                    
                    v1n0 = -( (v1_x * n1_x) + (v1_y * n1_y) + (v1_z * n1_z) );

                    rho10 = mdot1 / (PI_4 * sq(d1) * vel1);

                    e0 = 1.957441607 * mu1 * sq(v1n0);
                    coeff = 8.3245e3 / (d0 * sq(mu_av1)) * rho10 * cube(v1n0);
                    e_start = 0.0;
                    h1 = e0/300.0;
                    ll1 = 0.0;
                    //puts("before first ODEINT");
                    status3d = odeint(&ll1, &e_start, e0, dldE1, &h1, eps, 0.0);

                    if (ll1 > dl2) {
                        dl1 = dl2;
                        dl2 = ll1;
                    }
                }


                if ((iter == MAX_ITER) && (fabs(dl2 - dl1)/dl1) >= eps) {
                    puts("Too many iterations in 1st shock width!");
                    exit(EXIT_FAILURE);
                }

                dl = (dl1 + dl2) / 2.0;  // Make final refimenent.
                //printf("dl = %e\n", dl);
                d1l_3d[i][j] = dl;

            }
        }
    }
    puts("3d first cooling layer widths done!");
    return 0;
}

double calc_width_2(){
    double d1, d2, e0, ll1, ll2, e_start, dl2, dl1, h1;
    int iter, status3d;
    double xxc, yyc,zzc, dl, rho20, v2n0;
    double eps = 1.0e-5;
    double vel1, vel2;
    double cos_x2, cos_y2, cos_z2;
    double d1s, d2s, d1_xyz;
    double v2_x, v2_y, v2_z, n2_x, n2_y, n2_z;
    double psi;
    short fast=0;
    psi = psi_deg * (PI_2/360);
    
    for(int j = 0; j<nc; j++){
        for(int i = 0; i<nc; i++){

            vel2 = v2(r_3d_2[i][j]);
            

        
            rho20 = mdot2 / (PI_4 * sq(r_3d_2[i][j]) * vel2); // Density of 1st wind at point of c.s.
            e0 = 1.957441607 * mu2 * sq(v2n_3d[i][j]); // Energy of normal component of the wind.

            coeff = 8.3245e3/d0/sq(mu_av2)*rho20*cube(v2n_3d[i][j]);

        // Compute first estimate of shock width based of energy e0 estimated by v1n0 at the c.s.
        // As the shock front is at some distance from the c.s., v1n is different so I have to refine
        // this estimate - see below.

            ll1 = 0.0;
            h1 = e0 / 300.0;
            e_start = 0.0;
            status3d = odeint(&ll1, &e_start, e0, dldE2, &h1, eps, 0.0);
            dl2 = ll1;


        /*
        Make iterations to refine the estimate of the shock width:

        1. At the current width dl2, compute velocity, energy e0, end coeff at the
        shock front position;
        2. With these, solve diff eq. to find a new width.
        3. Set the new width as average of the old and new ones.
        4. Repeat 1-3 until convergency.

        */

            if (fast) { // If fast, just take this first estimate as an answer.
                d2l_3d[i][j] = dl2;
            } else {    // Refine the estimate.
                dl2 = MIN(dl2, 0.9 * (d-x10 - rstar2));
                iter = 0;
                dl1 = 0.0;
                while ((fabs((dl2-dl1)/dl2) > eps) && (iter < MAX_ITER)) {
                    iter++;
                    dl2 = (dl1 + dl2) / 2.0;
                    
    
                    

                    
                
                  

                    xxc = ((x_3d[i][j] -d) + dl2 * sphi[i]);
                    yyc = (y_3d[i][j] - dl2 * cphi[i]);
                    zzc = (z_3d[i][j]);

                    d2 = sqrt(sq(xxc) + sq(yyc) + sq(zzc));

                    vel2 = v2(d2);

                    cos_x2 = xxc / d2;
                    cos_y2 = yyc / d2;
                    cos_z2 = zzc / d2;

                    v2_x = (vel2 * cos_x2);
                    v2_y = (vel2 * cos_y2);
                    v2_z = (vel2 * cos_z2);
                    /*
                    n1_x = (-cos(psi)*sin(PI/4) - sin(psi)*cos(PI/4)*cos(theta_3d[j]));
                    n1_y = (-sin(psi)*sin(PI/4) + cos(psi)*cos(PI/4)*cos(theta_3d[j]));
                    n1_z = (cos(PI/4) * sin(theta_3d[j]));
                    */
                    n2_x = (-cos(psi)*sphi[i] - sin(psi)*cphi[i]*cos(theta_3d[j]));
                    n2_y = (-sin(psi)*sphi[i] + cos(psi)*cphi[i]*cos(theta_3d[j]));
                    n2_z = (cphi[i] * sin(theta_3d[j]));

                    nor2_x[i][j] = n2_x;
                    nor2_y[i][j] = n2_y;
                    nor2_z[i][j] = n2_z;
                

                

                    
                    v2n0 = ( (v2_x * n2_x) + (v2_y * n2_y) + (v2_z * n2_z) ) ;
                    rho20 = mdot2 / (PI_4 * sq(d2) * vel2);

                    e0 = 1.957441607 * mu2 * sq(v2n0);
                    coeff = 8.3245e3 / (d0 * sq(mu_av2)) * rho20 * cube(v2n0);
                    e_start = 0.0;
                    h1 = e0/300.0;
                    ll1 = 0.0;
                    //puts("before first ODEINT");
                    status3d = odeint(&ll1, &e_start, e0, dldE2, &h1, eps, 0.0);

                    if (ll1 > dl2) {
                        dl1 = dl2;
                        dl2 = ll1;
                    }
                }


                if ((iter == MAX_ITER) && (fabs(dl2 - dl1)/dl1) >= eps) {
                    puts("Too many iterations in 2st shock width!");
                    exit(EXIT_FAILURE);
                }

                dl = (dl1 + dl2) / 2.0;  // Make final refimenent.
                //printf("dl = %e\n", dl);
                d2l_3d[i][j] = dl;

            }
        }
    }
    puts("3d second cooling layer widths done!");
    return 0;
}

double calc_cross_obsorb(){
    ear1[ne1+1] = 40.000;
    for( k=0; k< ne1; k++){
        eee1 = (ear1[k]+ear1[k+1])/2.0;
        if(eee1 <= eabs1[0]){
            cross_abs1[k]= absorp1[0];    
        }else if(eee1 >= eabs1[nabs1-1]){    
            cross_abs1[k]= absorp1[nabs1-1];
        }else{
            cross_abs1[k]= interpolate(eabs1, absorp1, nabs1, eee1);    
        }


        if(eee1 <= eabs2[0]){
            cross_abs2[k]= absorp2[0];
        }else if(eee1 >= eabs2[nabs2-1]){
            cross_abs2[k]= absorp2[nabs2-1];
        }else{
            cross_abs2[k]= interpolate(eabs2, absorp2, nabs2, eee1);        
        }
 
    }
    e0max1 = exp(e0min1 + ne01*e0bin1);
    e0max2 = exp(e0min2 + ne02*e0bin2);
}

/*
calc_spectra(){
    
    for( jj=0 ; jj<1; jj++) {
        phase = 0.3;

        anomaly(phase, &vt, &d);

        phi_angle = vt + PI*0.5 + omega;

        phi_angle = normalize_angle(phi_angle);

        printf("d= %e  vt= %e phi_angle= %e\n",d, vt, phi_angle*180/PI);

        printf("%e %e %0.9e\n", vt, d, phi_angle);

        lsz = cos(incl);
        lsx = sin(incl) * cos(phi_angle);
        lsy = sin(incl) * sin(phi_angle);

        flag_eta = 0;
        exit_status = shock_distance();
        if (exit_status) exit(exit_status);

        eta0 = (mdot1 * v1(x10)) / (mdot2 * v2(d - x10));

        shock(x10, dy);

        anomaly(phase, &vt, &d);

        
        //calc_phi_3d(psi_deg);
        init_coordinates_in_x2y2z2 = make_3d();
        init_coordinates_in_x1y1z1 = shift();
        init_coordinates_in_xyz = skew();
        get_r3d = get_r();

        calculate_velocities = calc_velocities();
        calculate_sigma1 = calc_sigma1();
        calculate_sigma2 = calc_sigma2();
        calculate_widths_1 = calc_width_1();
        calculate_widths_2 = calc_width_2();

        for(i = 0; i<ne1-1; i++){
            spectrum[i]= 0.0;
            spectrum_int[i]= 0.0;
        }
        ekin1= 0.0;
        ekin2= 0.0;
//start of cycle for the C.S.
        for(i = 0; i<nc-1; i++){ //cycle of points on a line
            for( j = 0; j<nc; j++ ){// cycle of lines

                ekin11 = 0.0;
                ekin12 = 0.0;
                ds1 = 0.0;
                //first shock
                dl= d1l_3d[i][j];
                xxc = (x_3d[i][j] - dl2 * sphi[i]);
                yyc = (y_3d[i][j] + dl2 * cphi[i]);
                zzc = (z_3d[i][j]);

                d1 = sqrt(sq(xxc) + sq(yyc) + sq(zzc));

                vel1 = v1(d1);

                cos_x1 = xxc / d1;
                cos_y1 = yyc / d1;
                cos_z1 = zzc / d1;

                v1_x = (vel1 * cos_x1);
                v1_y = (vel1 * cos_y1);
                v1_z = (vel1 * cos_z1);

                v1n0 = -( (v1_x * n1_x[i][j]) + (v1_y * n1_y[i][j]) + (v1_z * n1_z[i][j]) );
                rho10 = mdot1 / (PI_4 * sq(d1) * vel1);
                eee1 = 1.957441607*mu1*sq(v1n0);
                printf("p0: %e \n",phase);
                printf("P1: %e %e %e \n",d,d1,vel1);
                printf("P2: %e %e %e \n",d1l[i],xxc,yyc);
                printf("P3: %e %e %e \n",v1n0,rho10,mu1);

                if( eee1 > e0max1 ){
                    printf("First cooling layer: Too high energy");
                    printf("beyound the planar grid\n E= %e  %e",eee1, e0max1);
                }
                //second shock
                dl = d2l_3d[i][j];

                xxc = ((x_3d[i][j] -d) + dl2 * sphi[i]);
                yyc = (y_3d[i][j] - dl2 * cphi[i]);
                zzc = (z_3d[i][j]);

                d2 = sqrt(sq(xxc) + sq(yyc) + sq(zzc));

                vel2 = v2(d2);

                cos_x2 = xxc / d2;
                cos_y2 = yyc / d2;
                cos_z2 = zzc / d2;

                v2_x = (vel2 * cos_x2);
                v2_y = (vel2 * cos_y2);
                v2_z = (vel2 * cos_z2);

                v2n0 = ( (v2_x * n2_x[i][j]) + (v2_y * n2_y[i][j]) + (v2_z * n2_z[i][j]) ) ;
                rho20 = mdot2 / (PI_4 * sq(d2) * vel2);
                
                eee2 = 1.957441607*mu2*sq(v2n0);
                if( eee1 > e0max1 ){
                    printf("second cooling layer: Too high energy");
                    printf("beyound the planar grid\n E= %e  %e",eee2, e0max2);
                }

                n1 = trunc((log(eee1)-e0min1)/e0bin1)+1;
                n2 = trunc((log(eee2)-e0min2)/e0bin2)+1;

printf("P4: %d %e %d %e %e %e %e %e\n", n1, eee1, n2, eee2, rho10, rho20, e0bin1, e0min1);


                for( k=0; k<ne1; k++){
                    if(n1 < 1){
                        spectr1[k] = 0.0;
                    }else{
                        spectr1[k]= sp_t1[n1-1][k]+(sp_t1[n1][k]-sp_t1[n1-1][k])/e0bin1*(log(eee1)-e0min1 - (n1-1)*e0bin1);                        
                    }

                    if(n2 < 1){
                        spectr2[k]= 0.0;
                        //puts("imhere");
                    }else{
                        spectr2[k]= sp_t2[n2-1][k]+(sp_t2[n2][k] -sp_t2[n2-1][k])/e0bin2  *(log(eee2)-e0min2 - (n2-1)*e0bin2);
                    }

                    spectr1[k] = spectr1[k]*rho10/sq(d0)*0.507686e10;
                    spectr2[k] = spectr2[k]*rho20/sq(d0)*0.507686e10;

                }

                if(flag_abs == 0){
                    //puts("flag_abs =0");
                    coldens11 = 0.0;
                    coldens12 = 0.0;
                    coldens21 = 0.0;
                    coldens22 = 0.0;
                    puts("flag abs = 0,");
                    goto skip;
                }

                r1 = sqrt(sq(x_3d[i][j])+sq(yc_3d[i][j])+sq(z_3d[i][j]));
                r1z = z_3d[i][j]/r1;
                r1y = y_3d[i][j]/r1;
                r1x = x_3d[i][j]/r1;

                coseta = lsx*n1_x[i][j] + lsy*n1_y[i][j] + lsz*n1_z[i][j];

                if(fabs(coseta)<=0.01){
                    puts("fabs(coseta)<0.01, tangential ray full absorbtion");
                    coldens11 = 9.9e30;
                    coldens21 = 9.9e30;
                    coldens22 = 9.9e30;
                    coldens12 = 9.9e30;
                    goto skip;
                }
//Check if the line of sight goes through the 1st wind only.
                if(coseta >= 0.0){

                    //  puts("coseta >= 0.0,line of site goes through 1st wind only");

                    ccphi = lsx*r1x + lsy*r1y + lsz*r1z;
                    if(fabs(ccphi) <= 1.0){
                        ssphi = sqrt(1.0 - sq(ccphi));
                    }else{
                        ssphi = 0.0;
                    }

                    p = r1*ssphi;
                    //printf("%e\n",p);

                    if(p <= rstar1){ //absorption by the body of star 1
                        puts("absorbtion by the body of star 1, skipped");
                        coldens11 = 9.9e30;
                        coldens12 = 9.9e30;
                        coldens21 = 9.9e30;
                        coldens22 = 9.9e30;
                        goto skip;
                    }
                    z1 = r1*ccphi;
                    z2 = sqrt(sq(RMAX)-sq(p));

                    nstar = 1;
                    coldens12 = 0.0;
                    //coldens11 = quad(&density_z,z1,z2,0.0,1.0e-2,1000,&Ierror);
                    //coldens11 = qsimp(&density_z,z1,z2);
                    coldens11 = qage(&density_z,z1,z2,0.0,1.0e-2,1000,1,&Ierror);
                    coldens22 = 0.0;
                    coldens21 = coldens11 + sig1_3d[i][j]/fabs(coseta);
                }else{ //line of sight goes into the second wind, trace the ray

                    r2 = sqrt(sq(d-x_3d[i][j])+sq(yc_3d[i][j])+sq(z_3d[i][j]));// Find p, z1 first.
                    r2z = z_3d[i][j]/r2;
                    r2y = y_3d[i][j]/r2;
                    r2x = -(d-x_3d[i][j])/r2;
//Angle between the radius-vector from the 2nd star to the point on the c.s.
//and the line of sight phi.
                    ccphi = lsx*r2x+lsy*r2y+lsz*r2z;

                    if(fabs(ccphi)<= 1.0){
                        ssphi = sqrt(1.0-sq(ccphi));
                    }else{
                        ssphi = 0.0;
                    }
                    p = r2*ssphi;
                    if(p <= rstar2){  //absorption by the body of star 2
                        coldens11 = 9.9e30;
                        coldens12 = 9.9e30;
                        coldens21 = 9.9e30;
                        coldens22 = 9.9e30;
                        goto skip;
                    }
                    z1 = r2*ccphi;
                    d2 = sqrt(sq(d-xc[nc][j])+sq(y_3d[nc][j])+sq(z_3d[nc][j]));
                    psimax = acos((d-x_3d[nc][j])/d2);
                    psi = acos((d-x_3d[i][j])/r2);
                    z2 = z1;
                    t = 0.0;
                    while(psi <= psimax){
                        t = t + d/100.0;
                        xxc = r1x*r1 + t*lsx;
                        yy = r1y*r1 + t*lsy;
                        zz = r1z*r1 + t*lsz;
                        yyc = sqrt(sq(yy)+sq(zz));
                        d2 = sqrt(sq(d-xxc)+sq(yyc));
                        psi = acos((d-xxc)/d2);

                        ps2 = 0.0;
                        dd2 = d-x_3d[0][j];
                        for( k=1; k<nc; k++){
                            dd1 = dd2;
                            ps1 = ps2;
                            dd2 = sqrt(sq(d-x_3d[k][j])+sq(y_3d[k][j])+sq(z_3d[k][j]));
                            ps2 = acos((d-x_3d[k][j])/dd2);
                            if(ps2 >= psi) break;
                        }

                        dc = dd1 + (dd2-dd1)/(ps2-ps1) * (psi-ps1);
                        if(d2 <= dc) break;
                    }
                    if(psi > psimax){ //No intersection, we are in the 2nd wind.
                        nstar = 2;
                        z2 = sqrt(sq(RMAX)-sq(p));
                        coldens21 = 0.0;
                        
                        //coldens22 = quad(&density_z,z1,z2,0.0,1.0e-2,1000,&Ierror);
                        //coldens22 = qsimp(&density_z,z1,z2);
                        coldens22 = qage(&density_z,z1,z2,0.0,1.0e-2,1000,1,&Ierror);
                        coldens11 = 0.0;
                        coldens12 = coldens22 + sig2_3d[i][j]/fabs(coseta);
                    }else{ //There is intersection.

//puts("psi is smaller than psimax: intersected with the second wind");
                        nstar = 2;
                        r2 = sqrt(sq(d-xxc)+sq(yy)+sq(zz));
                        r2x = -(d-xxc)/d2;
                        r2y = yy/d2;
                        r2z = zz/d2;
                        ccphi = lsx*r2x + lsy*r2y +lsz*r2z;
                        p = d2*sqrt(1.0-sq(ccphi));
                        z2 = d2*ccphi;
                        //dummy2 = quad(&density_z,z1,z2,0.0,1.0e-2,1000,&Ierror);
                        //dummy2 = qsimp(&density_z,z1,z2);
                        dummy2 = qage(&density_z,z1,z2,0.0,1.0e-2,1000,1,&Ierror);
                        
//1st wind. Find p, z of the intersection point for the 1st wind.
                        nstar = 1;
                        r1 = sqrt(sq(xxc)+sq(yy)+sq(zz));
                        r1x = xxc/r1;
                        r1y = yy/r1;
                        r1z = zz/r1;
                        ccphi = lsx*r1x + lsy*r1y + lsz*r1z;
                        if(fabs(ccphi) <= 1.0){
                            ssphi = sqrt(1.0-sq(ccphi));
                        }else{
                            ssphi = 0.0;
                        }
                        p = r1*ssphi;
                        if(p <= rstar1){
                            coldens11 = 9.9e30;
                            coldens12 = 9.9e30;
                            coldens21 = 9.9e30;
                            coldens22 = 9.9e30;
                            goto skip;
                        }
                        z1 = r1*ccphi;
                        z2 = sqrt(sq(RMAX)-sq(p));
                        //dummy1 = quad(&density_z,z1,z2,0.0,1.0e-2,1000,&Ierror);
                        //dummy1 = qsimp(&density_z,z1,z2);
                        dummy1 = qage(&density_z,z1,z2,0.0,1.0e-2,1000,1,&Ierror);
                        
//Column density of the cooling layers in the intersection point.
//c Linear interpolation for now.
                        for(int qqq = 0; qqq < nc; qqq++){
                            y_interp[qqq] = y_3d[qqq][j];
                            sig1_interp[qqq] = sig1_3d[qqq][j];
                            sig2_interp[qqq] = sig2_3d[qqq][j];
                        }
                        dum1 = interpolate(y_interp, sig1_interp, nc, yyc);
                        dum2 = interpolate(y_interp, sig2_interp, nc, yyc);
                        
                        nx1 = -sphi[k];
                        ny1 = cphi[k]*yy/yyc;
                        coseta1 = lsx*nx1 + lsy*ny1 + lsz*nz1;
//absorbtion my the material of the 1st wind
                        if(fabs(coseta1)< 0.01){
                            coldens11 = 9.9e30;
                            coldens21 = 9.9e30;
                        }else{
                            coldens11 = dum1/fabs(coseta1) + dummy1;
                            coldens21 = coldens11;
                        }
                        //Absorption by the material of the 2nd wind.

//1st intersection
                        coldens12 = sig2_3d[i][j]/fabs(coseta);
//2nd intersection.
                        if(fabs(coseta1) < 0.01){
                            coldens12 = 9.9e30;
                            coldens22 = 9.9e30;
                        }else{
                            coldens12 = coldens12 + dum2/fabs(coseta1);
                            coldens22 = dum2/fabs(coseta1);
                        }
                        coldens12 = coldens12 + dummy2;
                        coldens22 = coldens22 + dummy2;
                        //printf("%e  %e\n",coldens11,coldens22);
                        //printf("%e  %e  %e  %e\n",coldens11,coldens12,coldens21 ,coldens22);
                    }

                    //printf("%e  %e  %e  %e\n",coldens11,coldens12,coldens21 ,coldens22);
                }

                coldens11 = coldens11*541.6787246/mu_av1/d0*nhr1;
                coldens12 = coldens12*541.6787246/mu_av2/d0*nhr2;
                coldens21 = coldens21*541.6787246/mu_av1/d0*nhr1;
                coldens22 = coldens22*541.6787246/mu_av2/d0*nhr2;

                if(coldens11 < 0.0){
                    printf("coldens11 = %e", coldens11);
                }else if(coldens21 < 0.0){
                    printf("coldens21 = %e", coldens21);
                }else if(coldens12 < 0){
                    printf("coldens12 = %e", coldens12);
                }else if(coldens22 < 0){
                    printf("coldens22 = %e", coldens22);
                }
//continue;//this is where the program jumps if no wind absorbtion is accounted for
skip:
                if(i == 0){
                    ds = delta_theta*sq(y_3d[1][j])/4;
                    //puts("i==1");

                }else{

                    //puts("i else");
                    dy = (y_3d[i][j]+y_3d[i+1][j])/2.0 - (y_3d[i][j]+y_3d[i-1][j])/2.0;

                    dummy = dy/sphi[i];
                    ds = delta_theta*y_3d[i][j]*dummy;
                }

                lx1 = 0.0;
                    for (k = 1; k<ne1+1 ; k++){
                        //eee1 = (ear1[k]+ear1[k+1])/2.0;
                        //puts("after skipping");
                        if(k == ne1){
                          cross_abs1[k] = cross_abs1[k-1];
                          cross_abs2[k] = cross_abs2[k-1];
                        }
                        tau1 = cross_abs1[k]*coldens11+cross_abs2[k]*coldens12;
                        dummy1 = spectr1[k-1]*exp(-tau1);
                        tau2 = cross_abs1[k]*coldens21+cross_abs2[k]*coldens22;
                        dummy2 = spectr2[k-1]*exp(-tau2);
                        //printf("%e  %e\n",cross_abs1[k],cross_abs2[k]);
                        //printf("%e  %e\n", dummy1, dummy2);
                        //printf("%.70e\n",exp(-tau1));
                        dummy1 = dummy1*ds*sq(d0*rsol);
                        dummy2 = dummy2*ds*sq(d0*rsol);
                        //printf( "%d %d %e %e %e %e %e\n", j, k, eee1, xc[i], yc[i], p, coseta);
                        //printf("%e  %e  %e\n",ds,d0,rsol);

                        spectrum[k-1] = spectrum[k-1] + dummy1 + dummy2;
                        //if (spectrum[k] < 0){
                          //fprintf(test_file,"%e  %e  %e  %e %e %e %e %e \n", dummy1, dummy2, tau1, tau2, spectr1[k], spectr2[k], ds, dy);
//                        spectrum[k] = spectrum[k-1]+ dummy1 + dummy2;
//                         if(spectrum[k-1] < 0){
//                           spectrum[k] = spectrum[k+1] + dummy1 + dummy2;
//                         }
//                         if(spectrum[k] < 0){
//                           spectrum[k] = spectrum[k+1];
//                          }
//
                        //}


                        //printf("d0= %e ds= %e\n", d0, ds );
                        spectrum_int[k] = spectrum_int[k]+ spectr1[k-1]*ds*sq(d0*rsol)+ spectr2[k-1]*ds*sq(d0*rsol);
                        lx1 = lx1+(dummy1+dummy2)*(ear1[k+1]-ear1[k]);
                        //printf("lx1= %e\n", lx1);
                        //printf("%e %e %e %e\n", dummy1, ds, ds*sq(d0*rsol), dummy);
                        //printf("%d %d %e %e  %e  %e  %e\n", j, k, spectrum[k], dummy1, dummy2, ds, ear1[k]);
                    }
                    e1 = 6.31e5*rho10*cube(v1n0)/2.0*ds;
                    e2 = 6.31e5*rho20*cube(v2n0)/2.0*ds;
                    ekin1 = ekin1 + e1;
                    ekin2 = ekin2 + e2;

                    ekin11 = ekin11 + e1;
                    ekin12 = ekin12 + e2;
                    ds1 = ds1 + ds;
                
                
            }//end of cycle for lines

        }//end of cycle for points
        fclose(test_file);

            lx = 0.0;
            lx1 = 0.0;
            lxint = 0.0;
            for(k = 0; k<ne1 ; k++){
                //puts("final loop");
                eee1 = (ear1[k] + ear1[k+1])/2.0 *1.6e-9;//photon energy in erg
                photons[k] = spectrum[k] * (ear1[k+1]-ear1[k])/eee1;
                lx = lx + spectrum[k]*(ear1[k+1]-ear1[k]);
                lxint = lxint + spectrum_int[k]*(ear1[k+1]-ear1[k]);
                //printf("%e  %e  %e  %e\n",lx,lxint,spectrum[k],ear1[k]);
                if(ear1[k] >= 0.5 && ear1[k] < 10.0){
                    lx1 = lx1 + spectrum[k]*(ear1[k+1]-ear1[k]);
                    //puts("imhere");
                    //printf("%e\n",lx1,spectrum[k]);
                }
                //printf("%e  %e  %e\n",lx,lxint,ekin1+ekin2);
                //printf("%e\n",spectrum[k]);
            }



        
    }

}


*/