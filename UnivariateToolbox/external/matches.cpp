#include <math.h>
#include "mex.h"


void filter2(double* y, double* y2, int N);
int runs(double* y, double r, int N, int* M, int *R);
double getm1(double* y, int n);
double getm2(double* y, int n);

void mexFunction (int nlhs, mxArray **plhs, int nrhs,
		  const mxArray **prhs)
{
	int mrows, ncols, N, *M, *R, *Mptr,*Rptr, i;
	double *Nptr, *rptr, r, *y2, *y, stdev;
	int dims[2];
	if(nrhs != 3)
		mexErrMsgTxt("three inputs required.\n");
	if(nlhs != 2)
		mexErrMsgTxt("two outputs required.\n");
	Nptr = mxGetPr(prhs[1]);
	rptr = mxGetPr(prhs[2]);
	r = *rptr;
	N = (int)(*Nptr);
	y = mxGetPr(prhs[0]);
	M = (int*)mxCalloc(N+2, sizeof(int));
	R = (int*)mxCalloc(N+2, sizeof(int));
	N = runs(y, r, N, M, R);
	dims[0] = N;
	dims[1] = 1;
	plhs[0] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
	Mptr = (int*)mxGetPr(plhs[0]);
	Rptr = (int*)mxGetPr(plhs[1]);
	M = M+1;
	R = R+1;
	for(i = 0; i < N; i++)
	{
		Mptr[i] =  M[i];
		Rptr[i] = R[i];
	}
}


int runs(double* y, double r, int N, int* M, int *R)
{
	for(int i = 0;i < N; i++)
	{	R[i] = 0; M[i] = 0;}
	
	for(int j = 1; j < N-1; j++)
	{	
		int k = 0;
		for(int i = 0; i < N-j; i++)
		{
			if(y[i+j]-y[i] < r && y[i+j]-y[i]>-r)
				k++;
			else if(k != 0)
			{ R[k] = R[k]+1; k = 0; }
		}
		if(k != 0)
			R[k]=R[k]+1;
	}
	M[0] = 0;
	for(i = 1; i < N; i++)
	{
		M[i] = R[i];
		for(j = i+1; j < N; j++)
			M[i] = M[i]+(j+1-i)*R[j];
	}
	i = 1;
	while(M[i++] != 0);
	return i;
}
