/* Evaluate a mixed Erlang kernel */
#include "mex.h"

/* Include source code (a bit of a hack to reuse code in multiple MEX files with compiling library files) */
/*#include "evaluate_kernel_fast.c"*/
#include "least_squares_kernel_objective_hessian_faster.c"

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Number of kernels */
    const int nz = mxGetM(prhs[3]);

    /* Number of time points */
    const int N = mxGetN(prhs[2]);

    /* Highest polynomial order minus one */
    const int M = mxGetNumberOfElements(prhs[0]) - 2;
    
    /* Declare inputs and outputs */
    double *theta, *c, *a, *t;
    double *hess;

    /* Measurements of true kernel */
    double *alphameas;
    
    /* Check for proper number of arguments */
    if (nrhs != 4 && nrhs != 5) {
        mexErrMsgIdAndTxt( "MATLAB:evaluate_kernel:invalidNumInputs",
                "Four or five input arguments required.");
    } else if(nlhs > 1) {
        mexErrMsgIdAndTxt( "MATLAB:evaluate_kernel:maxlhs",
                "Too many output arguments.");
    } /* End if number of inputs and outputs are incorrect */
    
    /* Assign pointers to the inputs */
    theta     = mxGetPr(prhs[0]);
    t         = mxGetPr(prhs[2]);
    alphameas = mxGetPr(prhs[3]);
    
    /* Access individual elements of theta */
    c = theta;
    a = theta + M+1;

    /* Allocate memory */
    plhs[0] = mxCreateDoubleMatrix(M+2, M+2, mxREAL);

    /* Assign pointer */
    hess = mxGetPr(plhs[0]);

    /* Do the actual computations in a subroutine */
    least_squares_kernel_objective_hessian_faster(t, c, a, alphameas, nz, M, N,
    //least_squares_kernel_objective_hessian_fast(t, c, a, alphameas, nz, M, N,
    //least_squares_kernel_objective_hessian(t, c, a, alphameas, nz, M, N,
        hess);
    
    /* Return nothing */
    return;
}