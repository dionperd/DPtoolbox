/*RKF45 with additive gaussian noise general purpose integrator*/


#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

typedef int (* function) (double t, const double y[], double dydt[], void * params);
/*We give NULL for Jacobian-we don't need it for RKF45*/
/*otherwise: 
 *typedef int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params);*/

/*Inputs
*function F: the function to be integrated of type function defined above
*int D: the dimensionality of F
*void *params: a structure containing all parameters of F
*double x[]: the state variables, a linearly indexed array of size D*iters*size(double)
*double dxdt[]: the rates of change of the state variables, size(x)
*double x0[]: the initial condition of the state variables of size D*size(double)
*double t[]: a vector of time points (mononically increasing) where x should be calculated, array of size iters*size(double)
*int iters: number of iterations
*double s_noise: noise amplitude, standard deviation of gaussian distribution with mean 0
*double abstol: margin of absolute error
*double reltol: margin of relative error
*/   

typedef struct _par_ {
	int D;
	/*other parameters*/
} par;


void integrator(function F, int D, void *params, double x[], double dxdt[], double x0[], double t[], int iters, double s_noise, double abstol, double reltol) {

	/*Temporary variables*/
    double *y = (double*) malloc( sizeof(double)*D ); /*state variable at time t-1 (input) and then at time t(output)*/
    double *dydt_in = (double*) malloc( sizeof(double)*D ); /*rate of change at time point t-1*/
    double *dydt_out= (double*) malloc( sizeof(double)*D ); /*rate of change at time point t*/
    double *yerr = (double*) malloc( sizeof(double)*D );/*error*/
    double t0,tf,tc,dt=(t[1]-t[0])/2,noise;/*initial time point, final time point, current time point, current time step, noise*/
    int j,ii;/*State variable and iteration indexes*/
    int status;/*integrator success flag*/
    
    /*Prepare noise generator*/
    const gsl_rng_type *Q;
    gsl_rng *r;
    gsl_rng_env_setup();
    Q = gsl_rng_default;
    r = gsl_rng_alloc(Q);

    /*Create a stepping function*/
    const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
    gsl_odeiv_step *s  = gsl_odeiv_step_alloc(T, D);
    /*Create an adaptive control function*/
    gsl_odeiv_control *c  = gsl_odeiv_control_y_new(abstol, reltol);
    /*Create the system to be integrated (with NULL jacobian)*/
    gsl_odeiv_system sys = {F, NULL, D, params};
    
    /*Initialize*/
    /*Calculate dx/dt for x0*/
    tc=t[0];
    /* initialise dydt_in from system parameters */
    /*GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);*/
    GSL_ODEIV_FN_EVAL(&sys, tc, x0, dydt_in);
    for (j=0;j<D;j++) {
        y[j] = x[j] = x0[j];
    }
            
    /*Integration*/
    for (ii=1; ii<iters; ii++) {
        
        /*Call the integrator*/
        /*int gsl_odeiv_step_apply(gsl_odeiv_step * s, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt)*/        
        t0=t[ii-1];
        tf=t[ii];
        tc=t0;
        
        while (tc<tf) {
            
            /*Constraint time step h such as that tc+h<=tf*/
            if (tc+dt>tf) dt=tf-tc;
            
            /*Advance a h time step*/
            status=gsl_odeiv_step_apply(s, tc, dt, y, yerr, dydt_in, dydt_out, &sys);
            if (status != GSL_SUCCESS) break;
            
            /*Modify time sep*/
            gsl_odeiv_control_hadjust(c,s,y,yerr,dydt_in,&dt);

            /*Increase current time*/
            tc += dt;

            /*Add noise*/
            for (j=0;j<D;j++) {
                noise=gsl_ran_gaussian_ziggurat(r, s_noise);
                y[j] += dt*noise;
                dydt_in[j]+=noise;
            }
        
        }
        
        /*Unpack and store result for this time point adding some gaussian noise with sigma s_noise*/
       if (status != GSL_SUCCESS) break;
        
        for (j=0;j<D;j++) { 
            
            x[ii*D+j] = y[j];
            dxdt[(ii-1)*D+j] = dydt_in[j];
            
            /*Prepare next step*/
            dydt_in[j] = dydt_out[j];
            
        }
 
    }
    /*Get dxdt for the last time point*/
    for (j=0;j<D;j++) dxdt[(iters-1)*D+j] = dydt_out[j];
        
    gsl_odeiv_control_free(c);
    printf("c freed\n");
    gsl_odeiv_step_free(s);
    printf("s freed\n");
    gsl_rng_free(r);
    printf("rng freed\n");
    free(yerr);
    printf("yerr freed\n");
    free(dydt_out);
    printf("dydt_out freed\n");
    free(dydt_in);
    printf("dydt_in freed\n");
    free(y);
    printf("y freed\n");

}