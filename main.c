#include "config.h"
#include "vars.h"
#include "lib.h"
#include "numerical.h"

void run() {
    flag_abs = 1;
    int k1, k2, k3, k4, k5, k6;
    int exit_status, init_coordinates_in_x2y2z2, 
        init_coordinates_in_x1y1z1, init_coordinates_in_xyz, 
        calculate_velocities,calculate_sigma1,calculate_sigma2, 
        get_r3d, get_dxdy_angle, calculate_widths_1, calculate_widths_2;
    double t1;
    double phi_angle;
    int i,jj,j;

    

    t1= 15000;
    t1 = 15000*1.38065812*1.0e-7/1.60217725;
    k = k1 = k2 = k3 = k4 = k5 = k6 = 0;

    nc = 101;
    dy = 0.04;

    rstar1 /= d0;
    rstar2 /= d0;
    incl *= PI / 180;
    omega *= PI / 180;

    dvinf1 = (log(vinf1_max) - log(vinf1_min)) / kv1; // step in log space
    //dbeta1 = (beta1_max - beta1_min) / (kb1-1);
    dmdot1 = (log(mdot1_max) - log(mdot1_min)) / km1;

    dvinf2 = (log(vinf2_max) - log(vinf2_min)) / kv2; // step in log space
    //dbeta2 = (beta2_max - beta2_min) / (kb2-1);
    dmdot2 = (log(mdot2_max) - log(mdot2_min)) / km2;

    cooling_read("./cooling_solar_ascii.dat", elambda1, lambda1, nergrid1, &nlambda1);
    cooling_read("./cooling_solar_ascii.dat", elambda2, lambda2, nergrid2, &nlambda2);
    planar_read("./planar_spectr_solar_ascii.dat", &e0min1, &e0bin1, &ne01, &emin1, &ebin1, &ne1, ear1, sp_t1);
    planar_read("./planar_spectr_solar_ascii.dat", &e0min2, &e0bin2, &ne02, &emin2, &ebin2, &ne2, ear1, sp_t2);
    ear1[ne1+1] = 40.000;
    puts("file opened4");
    opacity_read("./opacity_warm_solar_ascii.dat", eabs1,absorp1, &nabs1);
    puts("file opened5");
    opacity_read("./opacity_warm_solar_ascii.dat", eabs2,absorp2, &nabs2);
    puts("file opened6");




    e0max1 = exp(e0min1 + ne01*e0bin1);
    e0max2 = exp(e0min2 + ne02*e0bin2);
    printf("e0max1= %e  e0max2= %e\n", e0max1, e0max2);
    atomic_weight( abund, t1, elambda1, nergrid1, nlambda1, &mu1, &mu_av1, &sum_nz1, &nhr1 );
    atomic_weight( abund, t1, elambda2, nergrid2, nlambda2, &mu2, &mu_av2, &sum_nz2, &nhr2 );
/*
    loop(k1, kb1) { loop(k2, kv1) { loop(k3, km1) { loop(k4, kb2) { loop(k5, kv2) { loop(k6, km2) {


        double coefv1 = dvinf1 * k2;
        double coefm1 = dmdot1 * k3;
        double coefv2 = dvinf2 * k5;
        double coefm2 = dmdot2 * k6;

        beta1 = beta1_min;// + dbeta1 * k1;
        vinf1 = vinf1_min * exp(coefv1);
        mdot1 = mdot1_min * exp(coefm1);

        beta2 = beta2_min;// + dbeta2 * k4;
        vinf2 = vinf2_min * exp(coefv2);
        mdot2 = mdot2_min * exp(coefm2);
*/
        vinf1 /= 1000;
        vinf2 /= 1000;
        rstar1 /= d0;
        rstar2 /= d0;
        incl *= PI / 180;
        omega *= PI / 180;
        
        


//        puts("main 1");
        for( jj=0 ; jj<1; jj++) {
            phase = 0.3;

            anomaly(phase, &vt, &d);

            phi_angle = vt + PI*0.5 + omega;

            phi_angle = normalize_angle(phi_angle);

            printf("d= %e  vt= %e phi_angle= %e\n",d, vt, phi_angle*180/PI);

            printf("%e %e %0.9e\n", vt, d, phi_angle);


            flag_eta = 0;
            exit_status = shock_distance();
            if (exit_status) exit(exit_status);

            eta0 = (mdot1 * v1(x10)) / (mdot2 * v2(d - x10));
            /*
            cone_x[0] = 0;
            cone_y[0] = 0;
            for(int i = 1; i < nc; i++){
                cone_x[i] = cone_x[i-1] + dy;
                cone_y[i] = cone_y[i-1] + dy;
                printf("%e\t%e\n",cone_x[i], cone_y[i]);
            }

            */
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
            
            //sub_x10(x10);
            
            



          }
}

/* main program */
int main(int argc, char *argv[]) {
    // parse arguments (fast, output mode, etc)
    init_test_values();

    run(); // this function should take the arguments and parameters

    output_results_c();
    three_d_results();
    return EXIT_SUCCESS;
}
