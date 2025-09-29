/* Headers */
#include <math.h>

/* Min. and max. functions */
#if !defined(MAX)
#define	MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

/* Convenience functions for matrix elements */
#define b(      i, m      )  (      b[(i) + nz*(m)                            ])
#define c(      i, m      )  (      c[(i) + nz*(m)                            ])
#define alpha(  i, n      )  (  alpha[(i) + nz*(n)                            ])
#define dalpha( i, n, j   )  ( dalpha[(i) + nz*(n) + nz*N*(j)                 ])
#define d2alpha(i, n, j, k)  (d2alpha[(i) + nz*(n) + nz*N*(j) + nz*N*(M+2)*(k)])

static void evaluate_kernel_fast(
    /* Inputs */
    const double *const t,
    const double *const c,
    const double *const a,
    const int           nz,
    const int           M,
    const int           N,
    const int           ComputeGradient,
    const int           ComputeHessian,
    
    /* Output */
    double *const alpha,
    double *const dalpha,
    double *const d2alpha
    ){
    /* Indices */
    int i, m, n, idx;

    /* Declare variables */
    double alpham, dalpham_a, d2alpham_a;

    /* Allocate memory */
    double *const b = (double *)calloc(nz*(M+1), sizeof(double));

    /* Initialize */
    for(i = 0; i < nz; i++){
        /* Initialize normalization factor */
        b(i, 0) = a[i];

        for(m = 1; m <= M; m++){
            /* Normalization factor */
            b(i, m) = a[i]*b(i, m-1)/m;
        } /* End for m */
    } /* End for i */

    for(n = 0; n < N; n++){
        for(m = 0; m <= M; m++){
            for(i = 0; i < nz; i++){
                /* Erlang probability density function */
                alpham = b(i, m)*pow(t[n], m)*exp(-a[i]*t[n]);

                /* Add to kernel */
                alpha(i, n) += c(i, m)*alpham;
            } /* End for i */
        } /* End for m */
    } /* End for n */

    if(ComputeGradient){
        for(m = 0; m <= M; m++){
            for(n = 0; n < N; n++){
                for(i = 0; i < nz; i++){
                    /* Erlang probability density function (TOBK: REPEATED COMPUTATION) */
                    alpham = b(i, m)*pow(t[n], m)*exp(-a[i]*t[n]);
                    
                    /* Derivative of density wrt. a */
                    dalpham_a = ((m+1)/a[i] - t[n])*alpham;
                    
                    /* Derivative of density wrt. cm */
                    dalpha(i, n, m) = alpham;

                    /* Derivative of kernel wrt. a */
                    dalpha(i, n, M+1) += c(i, m)*dalpham_a;
                } /* End for i */
            } /* End for n */
        } /* End for m */

        if(ComputeHessian){
            for(m = 0; m <= M; m++){
                for(n = 0; n < N; n++){
                    for(i = 0; i < nz; i++){
                        /* Erlang probability density function (TOBK: REPEATED COMPUTATION) */
                        alpham = b(i, m)*pow(t[n], m)*exp(-a[i]*t[n]);
                        
                        /* Derivative of density wrt. a (TOBK: REPEATED COMPUTATION) */
                        dalpham_a = ((m+1)/a[i] - t[n])*alpham;

                        /* Derivative of density wrt. a */
                        d2alpham_a = -(m+1)/(a[i]*a[i])*alpham + ((m+1)/a[i] - t[n])*dalpham_a;

                        /* Mixed derivatives */
                        d2alpha(i, n, m, M+1) = d2alpha(i, n, M+1, m) = dalpham_a;

                        /* Derivative of kernel wrt. a */
                        d2alpha(i, n, M+1, M+1) += c(i, m)*d2alpham_a;
                    } /* End for i */
                } /* End for n */
            } /* End for m */
        } /* End if ComputeHessian */
    } /* End if ComputeGradient */

    /* Free memory */
    free(b);
    
    /* Return nothing */
    return;
}