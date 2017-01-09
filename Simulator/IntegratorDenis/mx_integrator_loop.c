/*Runge Kutta or Euler integrator
%The MATLAB function: [x,dxdt]=integrator(System,param,x0,t,ip)
%Inputs:
%-System:  system to be integrated, 'lorenz, 'rossler' or 'linear' (default)
%-param: structure with the parameters of System that may also change in time
%-x0: state initial condition
%-t: time vector, 
%    i.e. time points where the system is going to be evaluated 
%    (for time step "fix": t(i+1)-t(i)=dt) 
%-ip: integraton parameters structure:
%   -ip.mode: "Euler" or "RK_45". Default: "RK_45" (Dormand-Prince method for adaptive 
%       step)
%   -ip.step: "var", "fix". 
%       Default: "fix" ("Euler" integration can be only with "fix" time
%       step)
%   -ip.min_step: minimum step size. Default: 10^(-6)
%   -ip.max_step: maximum step size. Default: 1
%   -ip.s_noise: standard deviation of gaussian noise. Default s_noise=0.
%                (NOTE!: for more accurate results use Runge-Kutta methods
%                for stochastic differential equations!)
%   -ip.err_tol: maximum absolute error allowed per calculation (only for
%                "var" step. Default: 10^(-6).
%
%Outputs:
%-x: solution system's state vector
%-dxdt: function's F evaluations as used for the integration 
*/

#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>


//Dimensionality of output time series as global vatiables
int N; //Number of time points
int D;//Dimensionality of the integrated system


/*Functions' types' definitions:*/
typedef void (* function) (double t, const double *x, double *dxdt,  void * params);


/*A structure type definition (parameters of the integration):*/
typedef struct _intgr_par_ {
    int mode;/*1 for "RK_45", 0 for "Euler". Default: "RK_45"*/
    int step; /*1 for "fix" time step, "0" for "var". Default: 0*/
    double min_step; /*minimum step size. Default: 10^(-6)*/
    double max_step; /*maximum step size. Default: 1*/
    double s_noise; /*noise amplitude (standard deviation of gaussian distribution with mean 0)*/
    double err_tol; /*maximum absolute error per iteration allowed   */
} intgr_par;



/*A structure type definition for linear system (parameters of the integrated function):*/
typedef struct _paramLinear_ {
    double tau;
    double* x0;
} paramLinear;


/*A structure type definition for excitator system (parameters of the integrated function):*/
typedef struct _paramExc_ {
    double tau;
    double mu;
    double a;
    double b;
    double c;
    double *I;
} paramExc;


/*A structure type definition for Lorenz system(parameters of the integrated function):*/
/*%For chaos: a=10, b=28, c=8/3
% p.tau=1;
% p.sigma=10;
% p.rho=28;
% p.beta=8/3;
 */
typedef struct _paramLorenz_ {
    double tau;
    double sigma;
    double rho;
    double beta;
    
} paramLorenz;

/*A structure type definition for Rossler system (parameters of the integrated function):*/
/*%For chaos: a=0.1, b=0.1, c=14
%Rossler studied it also for a=0.2, b=0.2, c=5.7
 */
typedef struct _paramRossler_ {
    double tau;
    double a;
    double b;
    double c;
} paramRossler;



/* The function to be integrated, a D-dimensional linear point attractor*/ 
/*dx/dt= -(1/tau)*x;*/
void linear( double t, const double *x, double *dxdt, void  *param) 
{
    
  //double *dxdt = malloc(sizeof(x));
  
  paramLinear *p = (paramLinear*) param;
  

  //Loop index
  int j;
  
  for (j=0; j<D; j++) {
     dxdt[j] = -1/(p->tau)*(x[j]-(p->x0)[j]);
  }
  
  //return dxdt;
  
}


/* The function to be integrated, a 2-dimensional excitator system*/ 
/*dxdt(1) = (x(2) + x(1) - x(1)^3)/p.tau;
  dxdt(2) = -p.mu*(p.a*x(2) + p.b*x(1) + p.c)/p.tau;*/
void excitator( double t, const double *x, double *dxdt, void  *param) 
{
    
  //double *dxdt = malloc(sizeof(x));
  
  paramExc *p = (paramExc*) param;

  dxdt[0] =  ( x[1] + x[0] - pow(x[0],3) )/p->tau;
  dxdt[1] = -(p->mu)*( (p->a)*x[1] + (p->b)*x[0] + p->c + (p->I[0]) )/p->tau;
  
  //return dxdt;
  
}

/* The function to be integrated, a Lorenz system*/ 
/*f(1) = p.sigma*( x(2)-x(1) )/p.tau;
f(2) = ( x(1)*(p.rho-x(3)) - x(2) )/p.tau;
f(3) = ( x(1)*x(2) - p.beta*x(3) )/p.tau;*/
void lorenz( double t, const double *x,  double *dxdt, void  *param) 
{
  //double *dxdt = malloc(sizeof(x));
   mexPrintf("311\n");
  paramLorenz *p = (paramLorenz*) param;
   mexPrintf("312\n");
  dxdt[0] = (p->sigma)*( x[1]-x[0] )/(p->tau);
   mexPrintf("313\n");
  dxdt[1] = ( x[0]*( (p->rho)-x[2] ) - x[1] )/(p->tau); 
   mexPrintf("313\n");
  dxdt[2] = ( x[0]*x[1] - (p->beta)*x[2] )/(p->tau);
   mexPrintf("315\n");
  
  //return dxdt;
  
}

/* The function to be integrated, a Rossler system*/ 
/*f(1) = ( -x(2)-x(3) )/p.tau;
f(2) = ( x(1) + p.a*x(2) )/p.tau;
f(3) = ( p.b + x(3)*(x(1)-p.c) )/p.tau;*/
void rossler( double t, const double *x, double *dxdt, void  *param) 
{
    
  //double *dxdt = malloc(sizeof(x));
  
  paramRossler *p = (paramRossler*) param;
  
  dxdt[0] = ( -x[1]-x[2] )/(p->tau);
  dxdt[1] = ( x[0] + (p->a)*x[1] )/(p->tau); 
  dxdt[2] = ( (p->b) + x[2]*(x[0]-(p->c)) )/(p->tau);
  
  //return dxdt;
  
}


//Linear system's function and parameters' structure construction
void linear_constr(function f, paramLinear* p, const mxArray* mxStruct) {
   
     int i;           
              
     
     f = &linear;
    
    
//     paramLinear *p=mxCalloc(1,sizeof(paramLinear));
//     mexPrintf("2\n");
          
    p->tau = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 0));
    mexPrintf("tau = %f\n", p->tau);
    
    p->x0 = (double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 1));
    
    //Temporary pointer to mxArray
//    int D;    
//    mxArray *tmp = mxGetFieldByNumber(mxStruct, 0, 1);
    
//     const mwSize *Dim=mxGetDimensions(tmp);
//     if (Dim[0]>Dim[1]) {
//         D=(int)Dim[0]; 
//         mexPrintf("D=M=%d\n",D);}
//     else {
//         D=(int)Dim[1]; ;
//         mexPrintf("D=N=%d\n",D);}
    
//     p->x0 = (double*)mxGetPr(tmp);

    for (i=0;i<D;i++) {
        printf("i = %d\n", i);
        printf("x0[i] = %f\n", p->x0[i]);
    }
    
    
    //return p;
    
}   

//Excitator system's function and parameters' structure construction
void exc_constr(function f, paramExc* p, const mxArray* mxStruct) {
   
     int i;           
              
     
     f = &excitator;
    
    
//     paramExc *p=mxCalloc(1,sizeof(paramExc));
//     mexPrintf("2\n");
          
    p->tau = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 0));
    mexPrintf("tau = %f\n", p->tau);
    
    p->mu = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 1));
    mexPrintf("mu = %f\n", p->mu);
    
    p->a = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 2));
    mexPrintf("a = %f\n", p->a);
    
    p->b = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 3));
    mexPrintf("b = %f\n", p->b);
    
    p->c = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 4));
    mexPrintf("c = %f\n", p->c);
    
    p->I = (double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 5));
    mexPrintf("I[t] = %f\n", p->I[0]);
    
    //Temporary pointer to mxArray
//    int D;    
//    mxArray *tmp = mxGetFieldByNumber(mxStruct, 0, 5);
    
//     const mwSize *Dim=mxGetDimensions(tmp);
//     if (Dim[0]>Dim[1]) {
//         D=(int)Dim[0]; 
//         mexPrintf("D=M=%d\n",D);}
//     else {
//         D=(int)Dim[1]; ;
//         mexPrintf("D=N=%d\n",D);}
    
//     p->x0 = (double*)mxGetPr(tmp);

//     for (i=0;i<N;i++) {
//         printf("i = %d\n", i);
//         printf("I[i] = %f\n", p->I[i]);
//     }
    
    
    //return p;
    
}  

//Lorenz system's function and parameters' structure construction    
void lorenz_constr(function f, paramLorenz* p, const mxArray* mxStruct) {
    

    f = &lorenz;
// 
//     paramLorenz *p=mxCalloc(1,sizeof(paramLorenz));
//     mexPrintf("2\n");

    
    p->tau = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 0));
    mexPrintf("tau = %f\n", p->tau);
    
    p->sigma = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 1));
    mexPrintf("sigma = %f\n", p->sigma);
    
    p->rho = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 2));
    mexPrintf("rho = %f\n", p->rho);
    
    p->beta = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 3));
    mexPrintf("beta = %f\n", p->beta);
      
    //return p;
    
}


//Rossler system's function and parameters' structure construction
void rossler_constr(function f, paramRossler* p, const mxArray* mxStruct) {
    
    
    f = &rossler;
    
//     paramRossler *p=mxCalloc(1,sizeof(paramRossler));;
//     mexPrintf("2\n");
    
    p->tau = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 0));
    mexPrintf("tau = %f\n", p->tau);
     
    p->a = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 1));
    mexPrintf("a = %f\n", p->a);
    
    p->b = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 2));
    mexPrintf("b = %f\n", p->b);
    
    p->c = *(double*)mxGetPr(mxGetFieldByNumber(mxStruct, 0, 3));
    mexPrintf("c = %f\n", p->c);
    
    
    //return p;
    
}

/*Inputs:
*function f: the function to be integrated of type function defined above
*int D: the dimensionality of F
*void *params: a structure containing all parameters of F
*double x[]: the state variables, a linearly indexed array of size D*iters*size(double)
*double dxdt[]: the rates of change of the state variables, size(x)
*double x0[]: the initial condition of the state variables of size D*size(double)
*double t[]: a vector of time points (mononically increasing) where x should be calculated, array of size iters*size(double)
*int N: number of time points in t[], where F has to be integrated
*par *ip: structure of integration parameters
*/      

/*This is the main function of the integration program that calls the gsl Runge-Kutta integrator:*/
void integrator(function f, void *p, double *x, double *dxdt, double x0[], double t[], intgr_par ip) {
    
     mexPrintf("0\n");
     
    /*State variable and iteration indexes*/
    int j, ii=-1;
    
     
     
    /*(initial) time step, noise*/
    double dt=t[1]-t[0], noise=0;

     
    /*Prepare noise generator*/
    const gsl_rng_type *Q;
    gsl_rng *r;
    gsl_rng_env_setup();
    Q = gsl_rng_default;
    r = gsl_rng_alloc(Q);
    
   
    
    /*Integration*/
    
    /*if variable step, use Dormand-Prince method (RK(4,5))...*/
    if (ip.step==1) {
        
        /*Temporary variables*/
        double *xc = (double*) malloc( sizeof(double)*D ); /*state variable at current time point*/
        double error=ip.err_tol;/*error initialization*/
        double tc=t[0], tr=t[N-1]-t[0];/*current time point, remaining time*/
        int CALCULATION_DONE=0;/* Flag for the successfull termination of an iteration's calculation*/
        
        /*Dormand-Prince method related temporary variables*/
        double *xk = (double*) malloc( sizeof(double)*D );
        double *k1 = (double*) malloc( sizeof(double)*D );
        double *k2 = (double*) malloc( sizeof(double)*D );
        double *k3 = (double*) malloc( sizeof(double)*D );
        double *k4 = (double*) malloc( sizeof(double)*D );
        double *k5 = (double*) malloc( sizeof(double)*D );
        double *k6 = (double*) malloc( sizeof(double)*D );
        double *k7 = (double*) malloc( sizeof(double)*D );
        double *dxdt5 = (double*) malloc( sizeof(double)*D );
        double *dxdt4 = (double*) malloc( sizeof(double)*D );
        
        /*Butcher's tableau for Dormand-Prince method*/
        double a[7][7] = {
            {0,           0,          0,           0,        0,          0    },
            {1/5,         0,          0,           0,        0,          0    },
            {3/40,        9/40,       0,           0,        0,          0    },
            {44/45,      -56/15,      32/9,        0,        0,          0    },
            {19372/6561, -25360/2187, 64448/6561, -212/729,  0,          0    },
            {9017/3168,  -355/33,     46732/5247,  49/176,  -5103/18656, 0    },
            {35/384,	  0,          500/1113,	   125/192,	-2187/6784,  11/84}
        };
        //%b5 = [a(7,:) 0] = [35/384 0 500/1113 125/192 -2187/6784 11/84 0]
        //%b4:
        double b[7] = {5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40};
        double c[7] = {0, 1/5, 3/10, 4/5, 8/9, 1, 1};
        
        
        /*Initialization*/
        for (j=0;j<D;j++) {
            //x[j] = xc[j] = x0[j];
            x[j] = xc[j] = x0[j];
        }
        
        
        /*Condition for the end of the integration*/
        while (tc<t[N-1])
            
            /*Increase iterations*/
            ii++;
        
        /*Set current time*/
        tc=t[ii];
        
        while (tc<t[ii+1]) {
            
            /*Remaining time*/
            tr = t[ii+1]-tc;
            
            CALCULATION_DONE=0;
            
            /*Calculation of next step*/
            while (CALCULATION_DONE==0) {
                
                /*Update time step*/
                dt = dt*pow(ip.err_tol/error, 0.2);
                
                /*Constraint time step in the allowed maximum limit*/
                if (dt>ip.max_step) {
                    dt=ip.max_step;}
                
                /*Constraint time step dt such as that tc+dt<=t[ii+1]*/
                if (dt>tr) {
                    dt=tr;}
                
                /*Constraint time step in the allowed minimum limit*/
                if (dt<ip.min_step) {
                    dt=ip.min_step;}
                
                
                /*Dormand-Prince method (Runge-Kutta 4-5):*/
                f(tc, xc, k1, p);
                
                for (j=0; j<D; j++) {
                    xk[j] = xc[j]+dt* k1[j]*a[2][1];
                }
                f(tc+c[2]*dt, xk, k2, p);
                
                for (j=0; j<D; j++) {
                    xk[j] = xc[j]+dt*( k1[j]*a[3][1] + k2[j]*a[3][2] );
                }
                f(tc+c[3]*dt, xk, k3, p);
                
                for (j=0; j<D; j++) {
                    xk[j] = xc[j]+dt*( k1[j]*a[4][1] + k2[j]*a[4][2] + k3[j]*a[4][3] );
                }
                f(tc+c[4]*dt, xk, k4, p);
                
                for (j=0; j<D; j++) {
                    xk[j] = xc[j]+dt*( k1[j]*a[5][1] + k2[j]*a[5][2] + k3[j]*a[5][3] + k4[j]*a[5][4] );
                }
                f(tc+c[5]*dt, xk, k5, p);
                
                for (j=0; j<D; j++) {
                    xk[j] = xc[j]+dt*( k1[j]*a[6][1] + k2[j]*a[6][2] + k3[j]*a[6][3] + k4[j]*a[6][4] + k5[j]*a[6][5]);
                }
                f(tc+c[6]*dt, xk, k6, p);
                
                /*k7 = f(tc(ii)+c(7)*dt, x(ii,:)+dt*( k1*a(7,1)+k2*a(7,2)+k3*a(7,3)+k4*a(7,4)+k5*a(7,5)+k6*a(7,6) )
                 * , param)
                 *%but also a(7,2)=0... and b5=[a(7,:) 0]...
                 */
                
                for (j=0; j<D; j++) {
                    dxdt5[j] = k1[j]*a[7][1] +k3[j]*a[7][3]+k4[j]*a[7][4]+k5[j]*a[7][5]+k6[j]*a[7][6];
                    xk[j] = xc[j] + dt*dxdt5[j];
                }
                f(tc+c[7]*dt, xk, k7, p);
                
                /*dxdt4 = b(1)*k1 + b(2)*k2 + b(3)*k3 + b(4)*k4 +b(5)*k5 +b(6)*k6 +b(7)*k7;
                 * but b(2)=0;*/
                for (j=0; j<D; j++) {
                    dxdt4[j] = b[1]*k1[j]  + b[3]*k3[j] + b[4]*k4[j] +b[5]*k5[j] +b[6]*k6[j] +b[7]*k7[j];
                }
                
                /*Absolute (euclidean distance) error*/
                error = 0;
                for (j=0; j<D; j++) {
                    error += pow(dxdt5[j]-dxdt4[j], 2);}
                
                error = dt*sqrt(error);
                
                if (error < ip.err_tol) {
                    CALCULATION_DONE = 1;}
                
            }
            
            /*Increase current time*/
            tc += dt;
            
            /*Add noise*/
            for (j=0;j<D;j++) {
                noise=gsl_ran_gaussian_ziggurat(r, ip.s_noise);
                xc[j] += dt*( dxdt4[j] + noise );
                
            }
            
        }
        /*Unpack and store result for this time point*/
        for (j=0;j<D;j++) {
            //x[(ii+1)*D+j] = xc[j];
            //dxdt[ii*D+j] = ( xc[(ii+1)*D+j]-x[ii*D+j] ) / ( t[ii+1]-t[ii] );
            x[ii*D+j] = xc[j];
            dxdt[ii*D+j] = ( x[(ii+1)*D+j]-x[ii*D+j] ) / ( t[ii+1]-t[ii] );
        }
        
        /*Free dynamically allocated memory*/
        free(dxdt4);
        free(dxdt5);
        free(k7);
        free(k6);
        free(k5);
        free(k4);
        free(k3);
        free(k2);
        free(k1);
        free(xk);
        free(xc);
        
    }
    
    
    
    
    /*...if fixed step...*/
    else {
        
        mexPrintf("1\n");
        
        /*and if Euler method...*/
        if (ip.mode==0) {
                        
            /*Initialization*/
            for (j=0;j<D;j++) {
                //x[j] = x0[j];
                x[j] = x0[j];
            }
            
            
            /*...integration...*/
            for (ii=0; ii<N; ii++) {
                
                /*Evaluate function*/
                f(t[ii], &(x[ii*D]), &(dxdt[ii*D]), p);
                
                for (j=0; j<D; j++) {
                    
                    /*Add noise */
                    noise=gsl_ran_gaussian_ziggurat(r, ip.s_noise);
                    //dxdt[ii*D+j]+= noise;
                    dxdt[ii*D+j]+= noise;
                    
                    /*Calculate next state*/
                    //x[(ii+1)*D+j] = x[ii*D+j] +  dxdt[ii*D+j];
                    x[(ii+1)*D+j] = x[ii*D+j] +  dt*dxdt[ii*D+j];
                    
                }
                
            }
            
            
            
        }
        /*and if Runge Kutta 4 method...*/
        else {
            mexPrintf("2\n");
        
            /*Temporary variables*/
            double *xc = (double*) malloc( sizeof(double)*D ); /*state variable at time t (input) and then at time t+1 (output)*/
            double *k1 = (double*) malloc( sizeof(double)*D );
            double *k2 = (double*) malloc( sizeof(double)*D );
            double *k3 = (double*) malloc( sizeof(double)*D );
            double *k4 = (double*) malloc( sizeof(double)*D );
            
            /*Butcher's tableau for Runge-Kutta 4rth order method*/
            double a[4][4] = {
                {0,   0,   0},
                {1/2, 0,   0},
                {0,   1/2, 0},
                {0,   0,   1},
            };
            double b[4] = {1/6, 1/3, 1/3, 1/6};
            double c[4] = {0, 1/2, 1/2, 1};
            
            
            /*Initialization*/
            for (j=0;j<D;j++) {
                //x[j] = xc[j] = x0[j];
                x[j] = xc[j] = x0[j];
            }
            
             mexPrintf("3\n");
             
            /*...integration...*/
            for (ii=0; ii<N; ii++) {
                
                 mexPrintf("31\n");
                 
                /*Evaluate function*/
                f(t[ii], &(x[ii*D]), k1, p);
                mexPrintf("32\n");
                for (j=0; j<D; j++) {
                    xc[j] = x[ii*D+j]+dt*k1[j]*a[2][1];
                }
                mexPrintf("33\n");
                f(t[ii]+c[2]*dt, xc, k2, p);
                mexPrintf("34\n");
                for (j=0; j<D; j++) {
                    xc[j] = x[ii*D+j]+dt*( k1[j]*a[3][1] + k2[j]*a[3][2] );
                }
                mexPrintf("35\n");
                f(t[ii]+c[3]*dt, xc, k3, p);
                for (j=0; j<D; j++) {
                    xc[j] = x[ii*D+j]+dt*( k1[j]*a[4][1] + k2[j]*a[4][2] + k3[j]*a[4][3] );
                }
                mexPrintf("36\n");
                f(t[ii]+c[4]*dt, xc, k4, p);
                mexPrintf("37\n");
                
                for (j=0; j<D; j++) {
                    mexPrintf("38\n");
                    /*Calculate vector field with noise*/
                    noise=gsl_ran_gaussian_ziggurat(r, ip.s_noise);
                    mexPrintf("39\n");
                    //dxdt[ii*D+j] = b[1]*k1[j] + b[2]*k2[j] + b[3]*k3[j] + b[4]*k4[j] + noise;
                    dxdt[ii*D+j] = b[1]*k1[j] + b[2]*k2[j] + b[3]*k3[j] + b[4]*k4[j] + noise;

                    /*Calculate next state*/
                    x[(ii+1)*D+j] = x[ii*D+j] +  dt*dxdt[ii*D+j];
                }
                
            }
             
            mexPrintf("4\n");
            
            free(k4);
            free(k3);
            free(k2);
            free(k1);
            free(xc);
        }
        
    }
    
    mexPrintf("5\n");
    
    /*Free dynamically allocated memory*/
    gsl_rng_free(r);
    /*printf("rng freed\n");*/
    
    mexPrintf("6\n");
}


    
    
/* The matlab command line that calls this mex file:*/ 
/*[x dxdt] = integrator_loop(int32(System),param,x0,t,ip);
 *prhs: 
 *-int System: 1 for 'lorenz', 2 for 'rossler', 0 for 'linear'. Default: 0.
 %-param: structure with the parameters of System that may also change in time
 %-x0: state initial condition
 %-t: time vector, 
%    i.e. time points where the system is going to be evaluated 
%    (for time step "fix": t(i+1)-t(i)=dt) 
 %-ip: integraton parameters structure:
%   -ip.mode: 0 for "Euler" or 1 for "RK_45". Default: 1 (Dormand-Prince method for adaptive 
%       step).
%   -ip.step: 1 for "var", 0 for "fix". Deafult: 0.
%   -ip.min_step: minimum step size. Default: 10^(-6)
%   -ip.max_step: maximum step size. Default: 1
%   -ip.s_noise: standard deviation of gaussian noise. Default s_noise=0.
%                (NOTE!: for more accurate results use Runge-Kutta methods
%                for stochastic differential equations!)
%   -ip.err_tol: maximum absolute error allowed per calculation (only for
%                "var" step. Default: 10^(-6).
*plhs: x, dxdt
 *
* The gateway function (this connects mex file to MATLAB) */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
 {
    



    /* Retrieve the input data */ 
    
    /* General form: type variablename = mxGetPr(prhs(index))*/
    
    //1.Identity of system to be integrated
    int System  = *((int*)mxGetData(prhs[0]));
    mexPrintf("System = %d\n", System);    /* Print in command line to make sure that everything works OK */ 
    
        

    /*Generate and initialize the system's parameters structure*/
    
    //Define a void pointer (for the parameters' structure) 
    //and a pointer to a function to be integrated
        //Initial definitions
//     void *p;
    function f;//=&lorenz;
    void *p;
    
    switch (System) {
    
        case 1:
            
        {//Construct pointers to function and parameters of Lorenz system
            
            paramExc *s=mxCalloc(1,sizeof(paramExc));
            p=s;
            exc_constr(f, p, prhs[1]);
        }
        
        break;
        
        case 2:
            
        {//Construct pointers to function and parameters of Lorenz system
            
            paramLorenz *s=mxCalloc(1,sizeof(paramLorenz));
            p=s;
            lorenz_constr(f, p, prhs[1]);
        }
        
        break;
        
        case 3:
            
        {//Construct pointers to function and parameters of rossler system
            paramRossler *s=mxCalloc(1,sizeof(paramRossler));
            p=s;
            rossler_constr(f, p, prhs[1]);}
        
        break;
        
        default:
            
        {//Construct pointers to function and parameters of linear system
            paramLinear *s=mxCalloc(1,sizeof(paramLinear));
            p=s;
            linear_constr(f, p, prhs[1]);
        }
        
    }
    
    
    //2.Parameters' structure
    
   

    //3.Initial condition vector 
    double *x0 = (double*)mxGetPr(prhs[2]);
    mexPrintf("got x0\n");
    //System's dimensionality
    D=mxGetM(prhs[2]);
    mexPrintf("D = %d\n", D);
    
    //4.Time vector 
    double *t = (double*)mxGetPr(prhs[3]);
    mexPrintf("got t\n");
    //System's dimensionality
    N=mxGetM(prhs[3]);
    mexPrintf("N = %d\n", N);
    
    
    //5. Integration parameters' structure
    
    intgr_par ip;// = (intgr_par)mxCalloc(1,sizeof(intgr_par)); 
    
    //Get 'mode'
    ip.mode = *( (int*)( mxGetData( mxGetFieldByNumber(prhs[4], 0, 0) ) ) );
    mexPrintf("mode = %d\n", ip.mode);
    
    //Get 'step'
    ip.step = *( (int*)( mxGetData( mxGetFieldByNumber(prhs[4], 0, 1) ) ) );
    mexPrintf("step = %d\n", ip.step);
    
    //Get 'min_step'
    ip.min_step = *( (double*)( mxGetData( mxGetFieldByNumber(prhs[4], 0, 2) ) ) );
    mexPrintf("min_step = %f\n", ip.min_step);
    
    //Get 'max_step'
    ip.max_step = *( (double*)( mxGetData( mxGetFieldByNumber(prhs[4], 0, 3) ) ) );
    mexPrintf("max_step = %f\n", ip.max_step);
    
    //Get 's_noise'
    ip.s_noise = *( (double*)( mxGetData( mxGetFieldByNumber(prhs[4], 0, 4) ) ) );
    mexPrintf("s_noise = %f\n", ip.s_noise);
    
    //Get 'err_tol'
    ip.err_tol = *( (double*)( mxGetData( mxGetFieldByNumber(prhs[4], 0, 5) ) ) );
    mexPrintf("err_tol = %f\n", ip.err_tol);
    

        
    /* Prepare/define the output data */
        
    /* Create a pointer to the output data */
    //x: state vector time series=solution of integration
    plhs[0] = mxCreateDoubleMatrix(N, D, mxREAL);/* This creates MATLAB variable/matrix variables in order to return/output them back to MATLAB */ 
    double *x=mxGetPr(plhs[0]);
    printf("made x\n");

    //dxdt: vector field time series=solution of integration
    plhs[1] = mxCreateDoubleMatrix(N-1, D, mxREAL);
    double *dxdt=mxGetPr(plhs[1]);
    printf("made dxdt\n");

    
       
    /* Now the show begins in C. Call the function integrator() */ 
    
    /* void integrator(function f, int D, void *p, double** x, double** dxdt, double x0[], double t[], int N, intgr_par ip) */
    integrator(f, p, x, dxdt, x0, t, ip);

    
    //Free allocated memory for parameters' structure
    mxFree(p);
                             
}

