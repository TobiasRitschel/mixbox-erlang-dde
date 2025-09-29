/* Headers */
#include <math.h>
#include <string.h>

/* Min. and max. functions */
#if !defined(MAX)
#define	MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

/* Convenience functions for matrix elements */
#define alphameas(i, n   )  (alphameas[(i) + nz*(n)               ])
#define alphahat( i, n   )  ( alphahat[(i) + nz*(n)               ])
#define dalphahat(i, n, j)  (dalphahat[(i) + nz*(n) + nz*(N-1)*(j)])

static void least_squares_kernel_objective(
    /* Inputs */
    const double *const t,
    const double *const c,
    const double *const a,
    const double *const alphameas, 
    const int           nz,
    const int           M,
    const int           N,
    const int           ComputeGradient,
    
    /* Output */
    double *const phi,
    double *const grad
    ){
    /* Indices */
    int i, n, m;

    /* Declare variables */
    double e, dt, de;

    /* Do not compute Hessian in this routine */
    const int ComputeHessian = false;

    /* Allocate memory */
    double *const  alphahat = calloc(nz*(N-1),       sizeof(double));
    double *const dalphahat = calloc(nz*(N-1)*(M+2), sizeof(double));

    /* Evaluate kernel */
    evaluate_kernel_faster(t, c, a, nz, M, N-1, ComputeGradient, ComputeHessian,
        alphahat, dalphahat, NULL);

    /* Initialize */
    *phi = 0.0;

    for(i = 0; i < nz; i++){
        for(n = 0; n < N-1; n++){
            /* Compute error */
            e = alphameas(i, n) - alphahat(i, n);

            /* Time step size */
            dt = t[n+1] - t[n];

            /* Erlang probability density function */
            *phi += e*e*dt;

            if(ComputeGradient){
                /* Loop over Erlang order */
                for(m = 0; m < M+2; m++){
                    /* Derivative of error */
                    de = -dalphahat(i, n, m);
                    
                    /* Add to derivative */
                    grad[m] += e*de*dt;
                } /* End for m */
            } /* End if ComputeGradient */
        } /* End for n */
    } /* End for i */

    /* Scale objective function */
    *phi *= 0.5;
    
    /* Free memory */
    free( alphahat);
    free(dalphahat);

    /* Return nothing */
    return;
}