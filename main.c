#include "config.h"
#include "vars.h"
#include "lib.h"
#include "numerical.h"

void run(double ph, double v_orb) {
    FILE *output_file_2;
    output_file_2 = fopen("model_4_full.dat", "a");

    flag_abs = 1;
    int k1, k2, k3, k4, k5, k6;
    int exit_status, init_coordinates_in_x2y2z2, 
        init_coordinates_in_x1y1z1, init_coordinates_in_xyz, 
        calculate_velocities,calculate_sigma1,calculate_sigma2, 
        get_r3d, get_dxdy_angle, calculate_widths_1, calculate_widths_2,
        calc_cross_obs, calculate_psi, calculate_vorb;
    double t1;
    double phi_angle, vt;
    int i,jj,j;
    

    //for ph in 0.0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.42 0.44 0.46 0.48 0.5 0.52 0.54 0.56 0.58 0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98 1.0
    //0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0


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

    double delta_theta = PI_2/100;

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
    //FILE *test_file;
    //test_file = fopen("colddens_3d.dat", "w");


    for(k=0; k< ne1; k++){

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
        
        lsz = cos(incl);
        


//        puts("main 1");
    // for( jj=0 ; jj<1; jj++) {
        //periastron at phase = 0
        phase = ph;
       
        anomaly(phase, &vt, &d);

        phi_angle = -vt + (PI*0.5) + omega;
        // while(phi_angle >= PI_2){
        //     phi_angle = phi_angle - PI_2;
        // }
        phi_angle = normalize_angle(phi_angle);

        // printf("d= %e  vt= %e phi_angle= %e\n",d, vt, phi_angle*180/PI);

        // printf("%e %e %0.9e %e\n", vt, d, phase, vorb_0);
        //exit(0);
       
        // lsx = sin(incl) * cos(phi_angle+(psi_deg*PI/180));
        // lsy = sin(incl) * sin(phi_angle-(psi_deg*PI/180));
        lsx = sin(incl) * cos(phi_angle);
        lsy = sin(incl) * sin(phi_angle);
        printf("lsx= %e lsy= %e\n", lsx , lsy);
        flag_eta = 0;
        exit_status = shock_distance();
        if (exit_status) exit(exit_status);

        eta0 = (mdot1 * v1(x10)) / (mdot2 * v2(d - x10));

        shock(x10, dy);

        // vorb_0 = v_orb;
        orb_per = 10;
        //calculate_vorb = calc_vorb();
        // vorb_0 = vorb;
        //vorb = vorb/10000000;
        calculate_psi = calc_psi();
        psi_deg = -psi_deg;
        printf("seperation= %e\n", d);
        //exit(0);
        //calc_phi_3d(psi_deg);
        init_coordinates_in_x2y2z2 = make_3d();
        init_coordinates_in_x1y1z1 = shift();
        init_coordinates_in_xyz = skew();
        get_r3d = get_r();
        //printf("x(nc)= %e\n", x_3d[nc-1][(nc-1)/2]);
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
        int aux_j = 0;
        for(i = 0; i<nc-1; i++){
        //for(i = 1; i<2; i++){
            for(j = 0; j<nc; j++ ){
             //cycle of points on a line
       
            
            // 3d aux cycle


                ekin11 = 0.0;
                ekin12 = 0.0;
                ds1 = 0.0;
                //first shock
                // dl= d1l_3d[i][aux_j];
                // xxc = (x_3d[i][aux_j] - dl * sphi[i]);
                // yyc = (y_3d[i][aux_j] + dl * cphi[i]);
                // zzc = (z_3d[i][aux_j]);

                dl= d1l_3d[i][j];
                xxc = (x_3d[i][j] - dl * sphi[i]);
                yyc = (y_3d[i][j] + dl * cphi[i]);
                zzc = (z_3d[i][j]);

                d1 = sqrt(sq(xxc) + sq(yyc) + sq(zzc));

                vel1 = v1(d1);

                cos_x1 = xxc / d1;
                cos_y1 = yyc / d1;
                cos_z1 = zzc / d1;

                v1_x = (vel1 * cos_x1);
                v1_y = (vel1 * cos_y1);
                v1_z = (vel1 * cos_z1);

                // v1n0 = -( (v1_x * nor1_x[i][aux_j]) + (v1_y * nor1_y[i][aux_j]) + (v1_z * nor1_z[i][aux_j]) );

                v1n0 = -( (v1_x * nor1_x[i][j]) + (v1_y * nor1_y[i][j]) + (v1_z * nor1_z[i][j]) );

                rho10 = mdot1 / (PI_4 * sq(d1) * vel1);
                eee1 = 1.957441607*mu1*sq(v1n0);
                // printf("p0: %e \n",phase);
                // printf("P1: %e %e %e \n",d,d1,vel1);
                // printf("P2: %e %e %e \n",d1l[i],xxc,yyc);
                // printf("P3: %e %e %e \n",v1n0,rho10,mu1);

                if( eee1 > e0max1 ){
                    printf("First cooling layer: Too high energy");
                    printf("beyound the planar grid\n E= %e  %e",eee1, e0max1);
                }
                //second shock
                // dl = d2l_3d[i][aux_j];
                // xxc = ((x_3d[i][aux_j]-d) + dl * sphi[i]);
                // yyc = (y_3d[i][aux_j] - dl * cphi[i]);
                // zzc = (z_3d[i][aux_j]);

                dl = d2l_3d[i][j];

                xxc = ((x_3d[i][j]-d) + dl*sphi[i]);
                yyc = (y_3d[i][j] - dl*cphi[i]);
                zzc = (z_3d[i][j]);

                d2 = sqrt(sq(xxc) + sq(yyc) + sq(zzc));

                vel2 = v2(d2);

                cos_x2 = xxc / d2;
                cos_y2 = yyc / d2;
                cos_z2 = zzc / d2;

                v2_x = (vel2 * cos_x2);
                v2_y = (vel2 * cos_y2);
                v2_z = (vel2 * cos_z2);
                

                // v2n0 = ( (v2_x * nor2_x[i][aux_j]) + (v2_y * nor2_y[i][aux_j]) + (v2_z * nor2_z[i][aux_j]) ) ;
                v2n0 = ( (v2_x * nor2_x[i][j]) + (v2_y * nor2_y[i][j]) + (v2_z * nor2_z[i][j]) ) ;
            
                rho20 = mdot2 / (PI_4 * sq(d2) * vel2);
                
                eee2 = 1.957441607*mu2*sq(v2n0);
                if( eee1 > e0max1 ){
                    printf("second cooling layer: Too high energy");
                    printf("beyound the planar grid\n E= %e  %e",eee2, e0max2);
                }

                n1 = trunc((log(eee1)-e0min1)/e0bin1)+1;
                n2 = trunc((log(eee2)-e0min2)/e0bin2)+1;

// printf("P4: %d %e %d %e %e %e %e %e\n", n1, eee1, n2, eee2, rho10, rho20, e0bin1, e0min1);


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
            //aux_j++;
             // cycle for lines
                if(flag_abs == 0){
                    //puts("flag_abs =0");
                    coldens11 = 0.0;
                    coldens12 = 0.0;
                    coldens21 = 0.0;
                    coldens22 = 0.0;
                    puts("flag abs = 0,");
                    goto skip;
                }

                r1 = sqrt(sq(x_3d[i][j])+sq(y_3d[i][j])+sq(z_3d[i][j]));
                r1z = z_3d[i][j]/r1;
                r1y = y_3d[i][j]/r1;
                r1x = x_3d[i][j]/r1;


                coseta = lsx*nor1_x[i][j] + lsy*nor1_y[i][j] + lsz*nor1_z[i][j];
                //fprintf(test_file,"%e\n",coseta);

                if(fabs(coseta)<=0.01){
                    // puts("fabs(coseta)<0.01, tangential ray full absorbtion");
                    coldens11 = 9.9e30;
                    coldens21 = 9.9e30;
                    coldens22 = 9.9e30;
                    coldens12 = 9.9e30;
                    goto skip;
                }
                //puts("didnt skip coseta");
//Check if the line of sight goes through the 1st wind only.
                if(coseta >= 0.0){

                    // puts("coseta >= 0.0,line of site goes through 1st wind only");

                    ccphi = lsx*r1x + lsy*r1y + lsz*r1z;
                    if(fabs(ccphi) <= 1.0){
                        ssphi = sqrt(1.0 - sq(ccphi));
                    }else{
                        ssphi = 0.0;
                    }

                    p = r1*ssphi;
                    //printf("%e\n",p);

                    if(p <= rstar1){ //absorption by the body of star 1
                        // puts("absorbtion by the body of star 1, skipped");
                        coldens11 = 9.9e30;
                        coldens12 = 9.9e30;
                        coldens21 = 9.9e30;
                        coldens22 = 9.9e30;
                        goto skip;
                    }
                    //puts("didnt skip psi");
                    z1 = r1*ccphi;
                    z2 = sqrt(sq(RMAX)-sq(p));

                    nstar = 1;
                    coldens12 = 0.0;
                    //coldens11 = quad(&density_z,z1,z2,0.0,1.0e-2,1000,&Ierror);
                    //coldens11 = qsimp(&density_z,z1,z2);
                    coldens11 = qage(&density_z,z1,z2,0.0,1.0e-2,1000,1,&Ierror);
                    coldens22 = 0.0;
                    coldens21 = coldens11 + sig1_3d[i][j]/fabs(coseta);
                    // fprintf(test_file,"z1= %e, z2= %e\n", z1, z2);
                }else{ //line of sight goes into the second wind, trace the ray
                    // puts("line of site goes through 2nd wind too");
                    r2 = sqrt(sq(d-x_3d[i][j])+sq(y_3d[i][j])+sq(z_3d[i][j]));// Find p, z1 first.
                    r2x = -(d-x_3d[i][j])/r2;
                    r2z = z_3d[i][j]/r2;
                    r2y = y_3d[i][j]/r2;
                    


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
                    // z1 = r2*ccphi;
                    // d2 = sqrt(sq(d-xc[nc-1])+sq(yc[nc-1]));
                    // //printf("d2= %e\n", d2);
                    // psimax = acos((d-xc[nc-1])/d2);
                    // psi = acos((d-xc[i])/r2);
                    // z2 = z1;
                    // t = 0.0;
                    //try here to swap j with (nc-1)/2
                    z1 = r2*ccphi;
                    d2 = sqrt(sq(d-x_3d[nc-1][j])+sq(y_3d[nc-1][j])+sq(z_3d[nc-1][j]));
                    psimax = acos((d-x_3d[nc-1][j])/d2);
                    psi = acos((d-x_3d[i][j])/r2);
                    z2 = z1;
                    t = 0.0;
                    
                    int flag_intersection;
                    while(psi <= psimax){    
                        //printf("while psi <= psimax i= %i j= %i psi= %e psimax= %e\n", i, j, psi, psimax);
                        t = t + d/100.0;
                        xxc = r1x*r1 + t*lsx;
                        yy = r1y*r1 + t*lsy;
                        zz = r1z*r1 + t*lsz;
                        yyc = sqrt(sq(yy)+sq(zz));
                        d2 = sqrt(sq(d-xxc)+sq(yyc));
                        if(d2>sqrt(sq(d-x_3d[nc-1][j])+sq(y_3d[nc-1][j])+sq(z_3d[nc-1][j]))){
                            flag_intersection = 0;
                            break;
                        }
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
                        if(d2 >= dc){
                            flag_intersection = 1;
                            break;
                        } 
                    
                    }
                    
                    if(psi > psimax){
                        flag_intersection = 0;
                    }
                    if(flag_intersection == 0){ //No intersection, we are in the 2nd wind.
                        // puts("psi> psimax");
                        nstar = 2;
                        z2 = sqrt(sq(RMAX)-sq(p));
                        coldens21 = 0.0;
                        //fprintf(test_file,"z1= %e, z2= %e\n", z1, z2);
                        //coldens22 = quad(&density_z,z1,z2,0.0,1.0e-2,1000,&Ierror);
                        //coldens22 = qsimp(&density_z,z1,z2);
                        coldens22 = qage(&density_z,z1,z2,0.0,1.0e-2,1000,1,&Ierror);
                        coldens11 = 0.0;
                        coldens12 = coldens22 + sig2_3d[i][j]/fabs(coseta);
                    }else{ //There is intersection.
                        // puts("psi< psimax");
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
                        //fprintf(test_file,"z1= %e, z2= %e\n", z1, z2);
//Column density of the cooling layers in the intersection point.
//c Linear interpolation for now.
                        for(int qqq = 0; qqq < nc; qqq++){
                            y_interp[qqq] = sqrt(sq(y_3d[qqq][j])+sq(z_3d[qqq][j]));
                            sig1_interp[qqq] = sig1_3d[qqq][j];
                            sig2_interp[qqq] = sig2_3d[qqq][j];
                        }
                        dum1 = interpolate(y_interp, sig1_interp, nc, yyc);
                        dum2 = interpolate(y_interp, sig2_interp, nc, yyc);
                     //not sure if this is right   
              
                        nx1 = nor1_x[k][j];
                        ny1 = nor1_y[k][j]*yy/yyc;
                        nz1 = nor1_z[k][j]*zz/yyc;
                        coseta1 = lsx*nor1_x[k][j] + lsy*nor1_y[k][j]+ lsz*nor1_z[k][j];
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
                    printf("coldens11 = %e j= %i i= %i\n", coldens11, j, i);
                }else if(coldens21 < 0.0){
                    printf("coldens21 = %e j= %i i= %i\n", coldens21, j, i);
                }else if(coldens12 < 0){
                    printf("coldens12 = %e j= %i i= %i\n", coldens12, j, i);
                }else if(coldens22 < 0){
                    printf("coldens22 = %e j= %i i= %i\n", coldens22, j, i);
                }
//continue;//this is where the program jumps if no wind absorbtion is accounted for
skip:
                if(i == 0){
                    ds = delta_theta*sq(y_3d_2[1][0])/4;
                    //puts("i==1");

                }else{

                    //puts("i else");
                    dy = (y_3d_2[i][0]+y_3d_2[i+1][0])/2.0 - (y_3d_2[i][0]+y_3d_2[i-1][0])/2.0;

                    dummy = dy/sphi[i];
                    ds = delta_theta*y_3d_2[i][0]*dummy;
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
                        dummy1 = dummy1 *ds*sq(d0*rsol);
                        dummy2 = dummy2 *ds*sq(d0*rsol);

                        //printf( "%d %d %e %e %e %e %e\n", j, k, eee1, xc[i], yc[i], p, coseta);
                        //printf("%e  %e  %e\n",ds,d0,rsol);

                        spectrum[k-1] = spectrum[k-1] + dummy1  + dummy2;
                        //if (spectrum[k] < 0){
                          //fprintf(test_file,"%e  %e  %e  %e %e %e %e %e \n", dummy1, dummy2, tau1, tau2, spectr1[k], spectr2[k], ds, dy);
/*                        spectrum[k] = spectrum[k-1]+ dummy1 + dummy2;
                          if(spectrum[k-1] < 0){
                            spectrum[k] = spectrum[k+1] + dummy1 + dummy2;
                          }
                          if(spectrum[k] < 0){
                            spectrum[k] = spectrum[k+1];
                          }
*/

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
                
                    //fprintf(test_file,"%0.7e %0.7e %0.7e %0.7e\n",tau1,tau2,dumdum1 ,dumdum2);
              
            }//end of cycle for lines
            //fprintf(test_file,"%0.7e %0.7e %0.7e %0.7e\n",coldens11,coldens12,coldens21 ,coldens22);
            // printf("eee1= %e eee2= %e\n", eee1, eee2);
        }//end of cycle for points


        //fclose(test_file);

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

            three_d_results(output_file_2);

        
    // }//cycle for phase
    fclose(output_file_2);
}

/* main program */
int main(int argc, char *argv[]) {
    if (argc < 3) {
        return -1;
    }
    double ph = strtod(argv[1], NULL);
    double v_orb = strtod(argv[2], NULL);
    // parse arguments (fast, output mode, etc)
    init_test_values();
    // clear file
    // for(int jj = 0; jj<20; jj++){
    //     // open file (O_APPEND)
    //     // clear memory
    //     // for loop with setting zeros
    //     run(file); // this function should take the arguments and parameters
    //     // close file
    // }
    run(ph, v_orb); // this function should take the arguments and parameters

    return EXIT_SUCCESS;
}
