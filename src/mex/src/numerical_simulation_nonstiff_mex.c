/* Evaluate a mixed Erlang kernel */
#include "mex.h"

/* Convenience functions for matrix elements */
#define X(           i, k) (X           [(i) + nx*(k)])
#define Z(           i, k) (Z           [(i) + nz*(k)])
#define R(           i, k) (R           [(i) + nz*(k)])
#define alphatkp1msj(i, k) (alphatkp1msj[(i) + nz*(k)])

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Number of states */
    const int nx = mxGetM(prhs[2]);

    /* Number of time points */
    const int N = mxGetNumberOfElements(prhs[1]) - 1;

    /* Auxiliary variables */
    int i, j, k;

    /* Declare Matlab function inputs and outputs */
    mxArray *fout[1], *finp[5], *hout[1], *hinp[3], *alphaout[1], *alphainp[2];

    /* Check for proper number of arguments */
    if(nrhs != 7) {
        mexErrMsgIdAndTxt( "MATLAB:evaluate_kernel:invalidNumInputs",
            "Seven input arguments required.");
    } else if(nlhs > 4) {
        mexErrMsgIdAndTxt( "MATLAB:evaluate_kernel:maxlhs",
            "Too many output arguments.");
    } /* End if number of inputs and outputs are incorrect */
    
    /* Rename inputs */
    mxArray *const x0_data = prhs[2];

    /* Function inputs that remain unchanged throughout the function */
    finp    [0] = mxDuplicateArray(prhs[0]); /* f     */
    finp    [4] = mxDuplicateArray(prhs[3]); /* p     */
    hinp    [0] = mxDuplicateArray(prhs[5]); /* h     */
    hinp    [2] = mxDuplicateArray(prhs[3]); /* p     */
    alphainp[0] = mxDuplicateArray(prhs[4]); /* alpha */

    /* Assign pointers to the inputs */
    double *const tspan = (double *)mxGetPr(prhs[1]);
    double *const p     = (double *)mxGetPr(prhs[3]);
    int           Nh    = (int)    *mxGetPr(prhs[6]);

    /* Time step size */
    const double dt = tspan[1] - tspan[0];

    /* Evaluate delayed quantity */
    hinp[1] = x0_data;
    mexCallMATLAB(1, hout, 3, hinp, "feval");

    /* Pointers to initial state and initial delayed state */
    const double *const x0 = (const double *const)mxGetPr(x0_data);
    const double *const r0 = (const double *const)mxGetPr(hout[0]);

    /* Number of kernels */
    const int nz = mxGetM(hout[0]);

    /* Inputs to right-hand side function */
    finp[1] = mxCreateDoubleMatrix((mwSize)  1, (mwSize) 1, mxREAL); /* tk */
    finp[2] = mxCreateDoubleMatrix((mwSize) nx, (mwSize) 1, mxREAL); /* xk */
    finp[3] = mxCreateDoubleMatrix((mwSize) nz, (mwSize) 1, mxREAL); /* zk */
    
    hinp[1] = mxCreateDoubleMatrix((mwSize) nx, (mwSize) 1, mxREAL); /* xk+1 */

    double *const tk_ptr = (double *const)mxGetPr(finp[1]);
    double *const xk_ptr = (double *const)mxGetPr(finp[2]);
    double *const zk_ptr = (double *const)mxGetPr(finp[3]);

    double *const xkp1_ptr = (double *const)mxGetPr(hinp[1]);

    /* Times for evaluating the kernel */
    mxArray *const t_alpha_data = mxCreateDoubleMatrix((mwSize) 1, (mwSize) Nh, mxREAL);
    double  *const t_alpha      = mxGetPr(t_alpha_data);
    for(k = 0; k < Nh; k++){ t_alpha[k] = (k+1)*dt; }

    /* Evaluate kernel */
    alphainp[1] = t_alpha_data;
    mexCallMATLAB(1, alphaout, 2, alphainp, "feval");

    /* Pointer to kernel evaluations */
    const double *const alphatkp1msj = (const double *const)mxGetPr(alphaout[0]);

    /* Allocate memory */
    plhs[0] = mxCreateDoubleMatrix((mwSize)  1, (mwSize) (Nh+N), mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize) nx, (mwSize) (Nh+N), mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize) nz, (mwSize) (Nh+N), mxREAL);
    plhs[3] = mxCreateDoubleMatrix((mwSize) nz, (mwSize) (Nh+N), mxREAL);

    /* Assign pointers */
    double *const T = (double *const)mxGetPr(plhs[0]);
    double *const X = (double *const)mxGetPr(plhs[1]);
    double *const Z = (double *const)mxGetPr(plhs[2]);
    double *const R = (double *const)mxGetPr(plhs[3]);
    
    /* Auxiliary variables */
    double *tk, *tkp1;
    double *xk, *zk, *rk, *fk, *hk, *xkp1, *zkp1, *rkp1;

    /* Number of columns in initial state function */
    int x0_cols = mxGetN(x0_data);

    /* Throw an error if an entire initial state function is passed */
    if(x0_cols > 1){
        mexErrMsgIdAndTxt("numerical_simulation_nonstiff_mex:InitialState", "x0 cannot have more than 1 column.");
    } /* End if x0_cols > 1 */

    for(k = 0; k < Nh; k++){
        /* Time */
        T[k] = (k-Nh+1)*dt;

        for(i = 0; i < nx; i++){
            /* Initial state */
            X(i, k) = x0[i];
        } /* End for i */

        for(i = 0; i < nz; i++){
            /* Initial delayed state quantity */
            Z(i, k) = r0[i];
            R(i, k) = r0[i];
        } /* End for i */
    } /* End for k */

    /* Initialize */
    tk = T + (Nh-1);
    xk = X + (Nh-1)*nx;
    zk = Z + (Nh-1)*nz;
    rk = R + (Nh-1)*nz;

    for(k = Nh-1; k < Nh+N-1; k++){
        /* Assign pointers */
        tkp1 = tk +  1;
        xkp1 = xk + nx;
        zkp1 = zk + nz;
        rkp1 = rk + nz;

        /* Time */
        tkp1[0] = tk[0] + dt;
        
        /* Delayed states */
        for(j = 0; j < Nh; j++){
            for(i = 0; i < nz; i++){
                zkp1[i] += alphatkp1msj(i, j)*R(i, k - j)*dt;
            } /* End for i */
        } /* End for j */

        /* Evaluate right-hand side function (TOBK: Can we avoid all of this copying?) */
                                 tk_ptr[0] = tk[0];
        for(i = 0; i < nx; i++){ xk_ptr[i] = xk[i]; }
        for(i = 0; i < nz; i++){ zk_ptr[i] = zk[i]; }
        mexCallMATLAB(1, fout, 5, finp, "feval");

        /* Increment states */
        fk = (double *)mxGetPr(fout[0]);
        for(i = 0; i < nx; i++){
            xkp1[i] = xk[i] + fk[i]*dt;
        } /* End for i */

        /* Evaluate delay function (TOBK: Can we avoid all of this copying?) */
        for(i = 0; i < nx; i++){ xkp1_ptr[i] = xkp1[i]; }
        mexCallMATLAB(1, hout, 3, hinp, "feval");

        /* Store delayed quantity (TOBK: Can we avoid all of this copying?) */
        hk = (double *)mxGetPr(hout[0]);
        for(i = 0; i < nz; i++){
            rkp1[i] = hk[i];
        } /* End for i */

        /* Deallocate memory from mexCallMATLAB inputs */
        mxDestroyArray(fout[0]);
        mxDestroyArray(hout[0]);

        /* Update pointers */
        tk = tkp1;
        xk = xkp1;
        zk = zkp1;
        rk = rkp1;
    } /* End for k */

    /* Deallocate memory from mexCallMATLAB inputs */
    mxDestroyArray(finp    [0]);
    mxDestroyArray(finp    [4]);
    mxDestroyArray(hinp    [0]);
    mxDestroyArray(hinp    [2]);
    mxDestroyArray(alphainp[0]);

    /* Return nothing */
    return;
}