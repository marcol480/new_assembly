#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>


#define N_PASSI_CAMPIONAMENTO_PER_CICLO 5000
#define N_STEP  N_PASSI_CAMPIONAMENTO_PER_CICLO*2*60 //1min. (10000 points is a second)
#define dt  1./((double)N_PASSI_CAMPIONAMENTO_PER_CICLO) //1/5000 sec

#define N_CELL 100
#define N_BARCODE N_CELL
#define VEC_COMP 2

#define R_CELL 10.
#define R_BARCODE 20.
#define L_BOX 3000

#define VISCH2O 0.001 //#define VISC (6.23)*VISCH2O
#define VISC 5*VISCH2O

#define TEMP 300
#define kBT 0.0000138*TEMP


#define OSEEN  1/(4*M_PI*VISC)
#define STOKES_CELL  1/(6*M_PI*VISC*R_CELL)
#define STOKES_BARCODE  1/(6*M_PI*VISC*R_BARCODE)

#define LAMBDA_CELL kBT/STOKES_CELL //fluctuation strength
#define LAMBDA_BARCODE kBT/STOKES_BARCODE //fluctuation strength

//0.004*0.19 //fluctuation strength

#define DENSTIY_WATER 1.
#define DENSTIY_CELL 1.05
#define DENSTIY_BARCODE 1.18

#define GRAVITY 9.81
#define V_CELL  (0.001*GRAVITY*4*M_PI*R_CELL*R_CELL)*(DENSTIY_CELL-DENSTIY_WATER)/(18*VISC)
#define V_BARCODE   (0.001*GRAVITY*4*M_PI*R_BARCODE*R_BARCODE)*(DENSTIY_BARCODE-DENSTIY_WATER)/(18*VISC)

#define strobo  200*5*1754691
#define RPM 10
//#define freq RPM*(M_PI/30)
#define freq 0.2

#define SEED 999999999

#define PRINT_CONFIG 0

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
// NEW DATA STRUCTURE
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

struct posforce{

    double x;
    double y;
    double fx;
    double fy;

};

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
// DECLAR FUNCTIONS
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

double compute_pos(double, double );

double compute_dist(struct posforce , struct posforce );

void randomize_pos(struct posforce *, double , int );

void randomize_pos_b(struct posforce *, double , int );

void update_cell_pos(struct posforce *, int *, double);

int control_cell_barcode_dist(int , int *, int *, int *, struct posforce *, struct posforce *);

int control_barcode_dist(int , int *, struct posforce * );


int main(){ // ---  MAIN

    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
    // DECLAR VAR
    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
    FILE  *xt, *xt_b, *out ;

    int a, i, j, k, l, amax;
    int counter;
    int active_cells[N_CELL], active_barcodes[N_BARCODE], count_attachments[N_BARCODE];
    int count_cell, count_barcode, count_mix;
    int count_one, count_double, count_triple, count_quadruple;

    double deltax, deltay;
    double dist_cell[N_CELL][N_CELL], dist_barcode[N_BARCODE][N_BARCODE], dist_mix[N_CELL][N_BARCODE];

    struct posforce PosForce_Cell[N_CELL], PosForce_Barcode[N_BARCODE];

    char  file_namext[350], file_namext_b[350], file_name[350];

    double sedim_cell, sedim_barcode;

    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
    // SETTINGS
    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

    printf("dt %lf \n", dt);

    printf("STOKES %lf\t OSEEN  %lf\n", STOKES_CELL, OSEEN );

    printf("Vol Fract: %lf\n", M_PI*(N_CELL*(R_CELL*R_CELL)+N_BARCODE*(R_BARCODE*R_BARCODE))/(L_BOX*L_BOX)  );

    /*  random  gen  */
    gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs0);

    /* sid setting */
    gsl_rng_set(r,SEED);

    if(PRINT_CONFIG==1) amax = 1;
    else amax=20;

    for(a=0;a<amax; a++){

        printf("run %d \n", a);

    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
    // DATA FILES
    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

        if(PRINT_CONFIG){

            sprintf(file_namext,"numerics/freq%0.1lf/config_cells_%dcells_%dbarcodes_noise-var%0.3lf_sep.dat", freq, N_CELL, N_BARCODE, LAMBDA_CELL);
            xt=fopen(file_namext, "w");

            sprintf(file_namext_b,"numerics/freq%0.1lf/config_barcodes_%dcells_%dbarcodes_noise-var%0.3lf_sep.dat", freq, N_CELL, N_BARCODE, LAMBDA_CELL);
            xt_b=fopen(file_namext_b, "w");
        }

        sprintf(file_name,"numerics/freq%0.1lf/events_%dcells_%dbarcodes_noise-var%0.3lf_run%d_sin_sep.dat", freq, N_CELL, N_BARCODE, LAMBDA_CELL, a);
        out = fopen(file_name, "w");



    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
    // INIZIALIZATION
    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

        for (k=0; k< N_CELL; k++) {

                PosForce_Cell[k].x = 0.0;
                PosForce_Cell[k].y = 0.0 ;

            }

        for (k=0; k< N_BARCODE; k++) {

            PosForce_Barcode[k].x = 0.0 ;
            PosForce_Barcode[k].y = 30.0 ;

            }

    randomize_pos( &PosForce_Cell[0], 1+a, N_CELL);
    randomize_pos_b( &PosForce_Barcode[0], 3+a, N_BARCODE);

    counter = 0;
    count_barcode=0;
    count_cell = 0;
    count_mix =0;

    for(l=0; l< N_CELL; l++) active_cells[l] = 1;
    for(l=0; l< N_BARCODE; l++){

        active_barcodes[l] = 1;
        count_attachments[l] = 0;

        }

    for (i=0;i<=N_STEP; i++){//-- being temporal evolution

        sedim_cell = -V_CELL*cos(freq*(double)i*dt);
        sedim_barcode = -V_BARCODE*cos(freq*(double)i*dt);

        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
        // CELL DYN
        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

        //update_cell_pos( PosForce_Cell, active_cells,  dt);

        for(l=0; l< N_CELL; l++){

            if( active_cells[l] == 1 ){

                deltax = compute_pos( sqrt(LAMBDA_CELL/(dt))*gsl_ran_gaussian(r,1.), 0.0 )*dt ;
                deltay = compute_pos( sqrt(LAMBDA_CELL/(dt))*gsl_ran_gaussian(r,1.), sedim_cell  )*dt ;

                PosForce_Cell[l].x += deltax ;
                PosForce_Cell[l].y += deltay ;

            }
        }

        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
        // BARCODE DYN
        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
        for(l=0; l< N_BARCODE; l++){

            if( active_barcodes[l] == 1 ){

                deltax = compute_pos( sqrt(LAMBDA_BARCODE/(dt))*gsl_ran_gaussian(r,1.), 0.0  )*dt ;
                deltay = compute_pos( sqrt(LAMBDA_BARCODE/(dt))*gsl_ran_gaussian(r,1.), sedim_barcode  )*dt ;

                PosForce_Barcode[l].x += deltax ;
                PosForce_Barcode[l].y += deltay ;

            }
        }


        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
        // CONTROL CELL-BARCODE DIST
        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

        count_mix = control_cell_barcode_dist(count_mix, &active_cells[0], &active_barcodes[0], &count_attachments[0], &PosForce_Cell[0], &PosForce_Barcode[0]);

        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
        // CONTROL MULTIPLE ATTACHMENTS
        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
        count_one = 0 ;
        count_double = 0;
        count_triple = 0;
        count_quadruple = 0 ;
        for(l=0; l< N_BARCODE; l++){

            if(count_attachments[l] == 1 )count_one++;

            if(count_attachments[l] == 2 )count_double++;

            if(count_attachments[l] == 3 )count_triple++;

            if(count_attachments[l] ==4 )count_quadruple++;

        }

        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
        // PRINT ON FILES
        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

        //if(i>N_STEP -10*40000)
        if(i%strobo==0){

            fprintf(out,"%lf\t %lf\t %lf\t %lf\t %lf\n", (double)i*dt, (100.*count_one)/((double)N_CELL), (100.*count_double)/((double)N_CELL), (100.*count_triple)/((double)N_CELL),  (100.*count_quadruple)/((double)N_CELL));


        if(PRINT_CONFIG){
            // cell config file print
            fprintf(xt, "%lf\t ", dt*((double)i) );
            for(l=0;l<N_CELL; l++) fprintf(xt, "%lf\t %lf\t", PosForce_Cell[l].x, PosForce_Cell[l].y );
            fprintf(xt, "\n");

            // barcode config file print
            fprintf(xt_b, "%lf\t ", dt*((double)i) );
            for(l=0;l<N_BARCODE; l++) fprintf(xt_b, "%lf\t %lf\t", PosForce_Barcode[l].x, PosForce_Barcode[l].y );
            fprintf(xt_b, "\n");

            }
        }

    } // -- end temporal evolution

    fclose(out);

    }// -- end realizations

    gsl_rng_free(r);

return(0);

} // --- END MAIN



/* - $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ - */
//
//                   --- FUNCTIONS --
//
/* - $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ - */


/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
// Randomize positions
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */


void randomize_pos(struct posforce *X, double fraction, int N){
    int k,l;
    double seed;

    gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs0);
    gsl_rng_set(r,12345*((int)fraction));


    for (k=0; k< N; k++) {

            X[k].x =  L_BOX*gsl_rng_uniform(r) ;
            X[k].y = 0.5*L_BOX*gsl_rng_uniform(r) ;

        }
    }

void randomize_pos_b(struct posforce *X, double fraction, int N){
    int k,l;
    double seed;

    gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs0);
    gsl_rng_set(r,12345*((int)fraction));


    for (k=0; k< N; k++) {

        X[k].x =  L_BOX*gsl_rng_uniform(r) ;
        X[k].y = 0.5*L_BOX*(1.+gsl_rng_uniform(r)) ;

    }
}



/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
// Compute distance
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

double compute_dist(struct posforce X1, struct posforce X2){

    double dist ;

    dist  = (X1.x - X2.x)*(X1.x - X2.x)  ;
    dist += (X1.y - X2.y)*(X1.y - X2.y)  ;
    dist = sqrt(dist);

    return dist;
}



/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
// Integrate pos
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

double compute_pos(double v , double coeff ){

    double pos;

    pos  = v + coeff;

    if(pos*dt > 0.3) printf("position = %g troppo!\n", pos), exit(0);
    if(isnan(pos) == 1) printf("position cell NANNA!"), exit(0);

    return pos;
}

double compute_pos_barcode(double v , double coeff){

    double pos;

    pos  = v + coeff;

    if(pos*dt > 0.3) printf("position = %g troppo!\n", pos), exit(0);
    if(isnan(pos) == 1) printf("position barcode NANNA!"), exit(0);

    return pos;
}


/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
// Check separations
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

// -- CELL
int control_cell_barcode_dist(int count, int  *active_cells, int *active_barcodes, int *count_attachments, struct posforce *X, struct posforce *Y){

    double dist_mix[N_CELL][N_BARCODE];
    int l, k;

    for(l=0; l< N_CELL; l++){
        for(k=0; k< N_BARCODE; k++){

            if( active_cells[l]==1 && active_barcodes[k]==1 ){

                dist_mix[l][k] = compute_dist(X[l], Y[k]);

                if(dist_mix[l][k] <  0.5*(R_CELL+R_BARCODE) ){


                    active_cells[l] = 0;
                    //active_barcodes[k] = 0;
                    count_attachments[k] += 1;
                    if(count_attachments[k]> 1); //printf("multiple attachment, N= %d, barcode %d, cell %d \n", count_attachments[k], k, l );
                    else count ++;
                }
            }
        }
    }

    return count;
}
