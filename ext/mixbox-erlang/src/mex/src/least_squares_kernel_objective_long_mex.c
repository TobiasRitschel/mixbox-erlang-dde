/* Evaluate a mixed Erlang kernel */
#include "mex.h"

/* Use long doubles for computing the kernel (this is the /only/ difference to least_squares_kernel_objective_mex) */
#define LONG_KERNEL_COMPUTATIONS

/* Include source code (a bit of a hack to reuse code in multiple MEX files with compiling library files) */
#include "least_squares_kernel_objective_faster.c"

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Number of kernels */
    const int nz = mxGetM(prhs[2]);

    /* Number of time points */
    const int N = mxGetN(prhs[1]);

    /* Highest polynomial order minus one */
    const int M = mxGetNumberOfElements(prhs[0]) - 2;
    
    /* Compute gradient? */
    int ComputeGradient = (nlhs > 1);

    /* Declare inputs and outputs */
    double *theta, *c, *a, *t;
    double *phi, *grad;

    /* Measurements of true kernel */
    double *alphameas;
    
    /* Check for proper number of arguments */
    if (nrhs != 3 && nrhs != 4) {
        mexErrMsgIdAndTxt( "MATLAB:evaluate_kernel:invalidNumInputs",
                "Three or four input arguments required.");
    } else if(nlhs > 2) {
        mexErrMsgIdAndTxt( "MATLAB:evaluate_kernel:maxlhs",
                "Too many output arguments.");
    } /* End if number of inputs and outputs are incorrect */
    
    /* Assign pointers to the inputs */
    theta     = mxGetPr(prhs[0]);
    t         = mxGetPr(prhs[1]);
    alphameas = mxGetPr(prhs[2]);
    
    /* Access individual elements of theta */
    c = theta;
    a = theta + M+1;
    
    /* Allocate memory */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);

    /* Assign pointer */
    phi = mxGetPr(plhs[0]);

    if(ComputeGradient){
        /* Allocate memory */
        plhs[1] = mxCreateDoubleMatrix(M+2, 1, mxREAL);
        
        /* Assign pointer */
        grad = mxGetPr(plhs[1]);
    } /* End if ComputeGradient */
    
    /* Do the actual computations in a subroutine */
    least_squares_kernel_objective_faster(t, c, a, alphameas, nz, M, N, ComputeGradient,
        phi, grad);
    
    /* Return nothing */
    return;
}