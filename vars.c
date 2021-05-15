#include "vars.h"




void cooling_read(char *fname, double *elambda, double *lambda, double *nergrid, int *nlambda){
    char str[200];
    int count;


    FILE *fp;
    if( !(fp = fopen(fname, "r")) ) {
        printf("Error opening file 1\n");
        exit(1);

    }
    //puts("imhere");
    fscanf( fp, "%*s %d\n", nlambda );
    fgets(str, LMAX, fp); // Skip one line.

    count =0 ;

    while (fscanf(fp, "%lf %lf %lf\n", &elambda[count], &lambda[count], &nergrid[count]) != EOF) {
        count++;
    }

    fclose( fp );

    if( *nlambda != count ) {
        printf( "Number of lines in input file %d not equal to nlambda= %d. Exiting.\n", count, *nlambda);
        exit( 0 );
    }

}

void planar_read(char *fname, double *e0min, double *e0bin, int *ne0, double *emin, double *ebin, int *ne, double ear[NE_MAX], double sp_t[NE0_MAX][NE_MAX]){


    FILE *fp;
    if( !(fp = fopen(fname, "r")) ) {
        printf("Error opening file 2\n");
        exit(1);
    }
    puts("imhere");
    fscanf( fp, "%*s %lf %*s %lf %*s %d %*s %lf %*s %lf %*s %d\n", e0min, e0bin, ne0, emin, ebin, ne);
    *ne0 = *ne0+1;

    for_i((*ne)+1){
        ear[i] = exp((*emin) + (i-1)*(*ebin));
    }

    //int k = 0;
    double tmp;
    for_i(*ne0) {
      if( exp((*e0min)+(i-1)*(*e0bin)) > 1.6 && k == 1 ){
        k = 3;
      }

        for_j(*ne) {
            if (fscanf(fp, "%lf\n", &tmp) != EOF) {
                sp_t[i][j] = tmp;
            }
            else {
                printf("Error: Unexpected EOF\n");
                break; break;
            }
        }
    }

    fclose( fp );

}



void opacity_read(char *fname, double *eabs, double *absorp, int *nabs){
    char str[200];
    int count1;


    FILE *fp;
    if( !(fp = fopen(fname, "r")) ) {
        printf("Error opening file 3\n");
        exit(1);

    }
    //puts("imhere");
    fgets(str, LMAX, fp); // Skip one line.
    fgets(str, LMAX, fp);

    count1 =0 ;
/*
    while ( fgets.(str..) {
        if(str[0] == '#') continue;
        sscanf( str, ....)
    }
*/
    while (fscanf(fp, "%lf %lf\n", &eabs[count1], &absorp[count1]) != EOF) {
        count1++;
    }

    fclose( fp );
    *nabs = count1;
}



void init_test_values() {
    N = 401;
//    ph[]= { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
   
    


    angr[0] = 12.00;
    angr[1] = 10.99;
    angr[2] = 8.56;
    angr[3] = 8.05;
    angr[4] = 8.93;
// !!!!!    angr[5] = 6.33;
    angr[5] = 8.09;
    angr[6] = 6.33;
    angr[7] = 7.58;
    angr[8] = 6.47;
    angr[9] = 7.55;
    angr[10] = 7.21;
    angr[11] = 6.56;
    angr[12] = 6.36;
    angr[13] = 7.67;
    angr[14] = 6.25;

    at[0] = 1.00794 ;
    at[1] = 4.002602;
    at[2] = 12.0107;
    at[3] = 14.0067;
    at[4] = 15.9994;
    at[5] = 20.1797;
    at[6] = 22.989770;
    at[7] = 24.3050;
    at[8] = 26.981538;
    at[9] = 28.0855;
    at[10] = 32.065;
    at[11] = 39.948;
    at[12] = 40.078;
    at[13] = 55.845;
    at[14] = 58.6934;

    abund[0] = 1;
    abund[1] = 1;
    abund[2] = 1;
    abund[3] = 1;
    abund[4] = 1;
    abund[5] = 1;
    abund[6] = 1;
    abund[7] = 1;
    abund[8] = 1;
    abund[9] = 1;
    abund[10] = 1;
    abund[11] = 1;
    abund[12] = 1;
    abund[13] = 1;
    abund[14] = 1;

    






    nb_ph = 10;

    d0 = 60.000;
    incl = 90.0;
    omega = 270.00;
    e = 0.3;

    mdot1 = 1;
    vinf1 = 2000.0;
    beta1 = 2.2;
    rstar1 = 10;

    vinf1_min = 2000.0;
    vinf1_max = 2000.0;
    kv1 = 1;
    beta1_min = 1.1;
    beta1_max = 1.0;
    kb1 = 1;
    mdot1_min = 2.0;
    mdot1_max = 2.0;
    km1 = 1;

    mdot2 = 0.5;
    vinf2 = 2000.0;
    beta2 = 1.3;
    rstar2 = 10;
   
    vinf2_min = 2000.0;
    vinf2_max = 2000.0;
    kv2 = 1;
    beta2_min = 1.0;
    beta2_max = 1.0;
    kb2 = 1;
    mdot2_min = 1.0;
    mdot2_max = 1.0;
    km2 = 1;
}



void three_d_results(FILE *output_file_2){

   
    fprintf(output_file_2,"x10= %lf eta0= %lf\n", x10, eta0);
    fprintf(output_file_2,"incl= %lf omega= %lf d0= %lf phase= %lf eccentricity= %lf Skew= %lf\n", incl, omega, d0, phase, e, psi_deg);
    fprintf(output_file_2,"mdot1= %lf beta1= %lf vinf1= %lf\n", mdot1, beta1, vinf1);
    fprintf(output_file_2,"mdot2= %lf beta2= %lf vinf2= %lf\n", mdot2, beta2, vinf2);
    fprintf(output_file_2,"mu1= %lf mu_av1= %lf sum_nz1= %lf nhr1= %lf\n", mu1, mu_av1, sum_nz1, nhr1);
    fprintf(output_file_2,"mu2= %lf mu_av2= %lf sum_nz2= %lf nhr2= %lf\n", mu2, mu_av2, sum_nz2, nhr2);
    fprintf(output_file_2,"contact surface:\n");
    fprintf(output_file_2,"x\t\ty\t\t z\t\tv1n\t\tv1t\t\tv2n\t        v2t\t\tsigma1\t\tsigma2\t\tr\t\ttheta [degrees]\td1l\t\td2l\n");
    //fprintf(output_file_2,"x_3d\t\ty_3d\t\tz_3d\n");
    for_j(nc){
        for_i(nc){
            fprintf(output_file_2,"%13.8le\t%13.8le\t%13.8le\t%13.8le\t%13.8le\t%13.8le\t%13.8le\t%13.8le\t%13.8le\t%13.8le\t%13.8le\t%13.8le\t%13.8le\n" 
                    ,x_3d[i][j],y_3d[i][j],z_3d[i][j],v1n_3d[i][j], v1t_3d[i][j],v2n_3d[i][j], v2t_3d[i][j]
                    , sig1_3d[i][j], sig2_3d[i][j], r_3d[i][j], theta_3d[j],d1l_3d[i][j], d2l_3d[i][j]);
            //fprintf(output_file_2,"%13.8le\t%13.8le\t%13.8le\n", x_3d[i][j],y_3d[i][j],z_3d[i][j]);
        }
         //fprintf(output_file_2,"%13.8le\t%13.8le\t%13.8le\n", x_skewed[i],y_skewed[i]);
    }
    fprintf(output_file_2, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

    fprintf(output_file_2, "Lxint= %lf [10^30 erg/s]   Lx= %lf [10^30 erg/s]   Lx[0.5-10.0]= %lf [10^30 erg/s]   Lx/Lkin= %lf [10^30 erg/s]  Lkin= %lf[10^30 erg/s]\n" , lxint, lx, lx1, lx/(ekin1+ekin2), ekin1+ekin2);


    fprintf(output_file_2,"(keV)\t     (keV)\t  (10^30 erg/keV/s)\t  (10^30 ph/bin/s)\t (10^30 erg/keV/s)\n");

    for(int i = 1; i < ne2+1; i++){
        fprintf(output_file_2,"%.4le   %.4le    %13.8le\t  %13.8le\t %13.8le\n",
                ear1[i], ear1[i+1], spectrum[i-1], photons[i-1], spectrum_int[i]);
    }
   
    puts("done writing 3D results");
}
