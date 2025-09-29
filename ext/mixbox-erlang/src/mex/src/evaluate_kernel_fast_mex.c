/* Evaluate a mixed Erlang kernel */
#include "mex.h"

/* Include source code (a bit of a hack to reuse code in multiple MEX files with compiling library files) */
#include "evaluate_kernel_fast.c"

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Number of kernels */
    const int nz = mxGetM(prhs[2]);

    /* Number of time points */
    const int N = mxGetN(prhs[0]);

    /* Highest polynomial order minus one */
    const int M = mxGetN(prhs[1]) - 1;
    
    /* Dimensions of outputs */
    const mwSize d[4] = {nz, N, M+2, M+2};

    /* Compute gradient and Hessian? */
    int ComputeGradient = (nlhs > 1);
    int ComputeHessian  = (nlhs > 2);

    /* Declare inputs and outputs */
    double *t, *c, *a;
    double *alpha, *dalpha, *d2alpha;
    
    /* Check for proper number of arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt( "MATLAB:evaluate_kernel:invalidNumInputs",
                "Three input arguments required.");
    } else if(nlhs > 3) {
        mexErrMsgIdAndTxt( "MATLAB:evaluate_kernel:maxlhs",
                "Too many output arguments.");
    } /* End if number of inputs and outputs are incorrect */
    
    /* Assign pointers to the inputs */
    t = mxGetPr(prhs[0]);
    c = mxGetPr(prhs[1]);
    a = mxGetPr(prhs[2]);
    
    /* Allocate memory */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nz, (mwSize)N, mxREAL);

    /* Assign pointer */
    alpha = mxGetPr(plhs[0]);

    if(ComputeGradient){
        /* Allocate memory */
        plhs[1] = mxCreateNumericArray(3, d, mxDOUBLE_CLASS, mxREAL);
        
        /* Assign pointer */
        dalpha = mxGetPr(plhs[1]);

        if(ComputeHessian){
            /* Allocate memory */
            plhs[2] = mxCreateNumericArray(4, d, mxDOUBLE_CLASS, mxREAL);

            /* Assign pointer */
            d2alpha = mxGetPr(plhs[2]);
        } /* End if ComputeHessian */
    } /* End if ComputeGradient */
    
    /* Do the actual computations in a subroutine */
    evaluate_kernel_fast(t, c, a, nz, M, N, ComputeGradient, ComputeHessian,
        alpha, dalpha, d2alpha);
    
    /* Return nothing */
    return;
}