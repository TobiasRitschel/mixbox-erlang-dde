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
#define alphameas( i, n      )  (alphameas [(i) + nz*(n)                                    ])
#define alphahat(  i, n      )  ( alphahat [(i) + nz*(n)                                    ])
#define dalphahat( i, n, j   )  (dalphahat [(i) + nz*(n) + nz*(N-1)*(j)                     ])
#define d2alphahat(i, n, j, k)  (d2alphahat[(i) + nz*(n) + nz*(N-1)*(j) + nz*(N-1)*(M+2)*(k)])
#define hess(      m, l      )  (hess      [(m) + (M+2)*(l)                                 ])

static void least_squares_kernel_objective_hessian_fast(
    /* Inputs */
    const double *const t,
    const double *const c,
    const double *const a,
    const double *const alphameas, 
    const int           nz,
    const int           M,
    const int           N,
    
    /* Output */
    double *const hess
    ){
    /* Indices */
    int i, n, m, l;

    /* Declare variables */
    double e, dt, dem, del, d2e;

    /* Do not compute Hessian in this routine */
    const int ComputeGradient = true;
    const int ComputeHessian  = true;

    /* Allocate memory */
    double *const   alphahat = calloc(nz*(N-1),             sizeof(double));
    double *const  dalphahat = calloc(nz*(N-1)*(M+2),       sizeof(double));
    double *const d2alphahat = calloc(nz*(N-1)*(M+2)*(M+2), sizeof(double));

    /* Evaluate kernel */
    evaluate_kernel_fast(t, c, a, nz, M, N-1, ComputeGradient, ComputeHessian,
        alphahat, dalphahat, d2alphahat);

    /* Loop over parameters */
    for(l = 0; l < M+2; l++){
        for(m = 0; m < M+2; m++){
            for(n = 0; n < N-1; n++){
                for(i = 0; i < nz; i++){
                    /* Compute error */
                    e = alphameas(i, n) - alphahat(i, n);

                    /* Time step size */
                    dt = t[n+1] - t[n];

                    /* Derivatives of error */
                    dem = -dalphahat (i, n, m   );
                    del = -dalphahat (i, n, l   );
                    d2e = -d2alphahat(i, n, m, l);

                    /* Add to derivative */
                    hess(m, l) += (dem*del + e*d2e)*dt;
                } /* End for i */
            } /* End for n */
        } /* End for m */
    } /* End for i */

    /* Free memory */
    free(  alphahat);
    free( dalphahat);
    free(d2alphahat);    

    /* Return nothing */
    return;
}