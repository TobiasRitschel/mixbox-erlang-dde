/* Headers */
#include <math.h>

#ifndef LONG_SUBKERNEL_COMPUTATIONS
    /* Exponential function for double computations */
    #define exponential_function(x) (exp(x))

    /* Data type for subkernel computations */
    typedef double real;
#else
    /* Exponential function for double computations */
    #define exponential_function(x) (expl(x))

    /* Data type for subkernel computations */
    typedef long double real;
#endif

/* Min. and max. functions */
#if !defined(MAX)
#define	MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

/* Convenience functions for matrix elements */
#define b(           i, m      ) (           b[(i) + nz*(m)                            ])
#define c(           i, m      ) (           c[(i) + nz*(m)                            ])
#define alpha(       i, n      ) (       alpha[(i) + nz*(n)                            ])
#define dalpha(      i, n, j   ) (      dalpha[(i) + nz*(n) + nz*N*(j)                 ])
#define d2alpha(     i, n, j, k) (     d2alpha[(i) + nz*(n) + nz*N*(j) + nz*N*(M+2)*(k)])
#define alpham_store(i, n, m   ) (alpham_store[(i) + nz*(n) + nz*N*(m)                 ])

static void evaluate_kernel_faster(
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
    real alpham, dalpham_a, d2alpham_a;

    /* Allocate memory */
    real *const alpham_store = (real *)calloc(N*nz*(M+1), sizeof(real));

    for(n = 0; n < N; n++){
        for(i = 0; i < nz; i++){
            /* Initialize subkernels */
            alpham =  a[i]*exponential_function(-a[i]*t[n]);

            /* Initialize kernel */
            alpha(i, n) = c(i, 0)*alpham;

            /* Store subkernel */
            alpham_store(i, n, 0) = alpham;

            for(m = 1; m <= M; m++){
                /* Erlang probability density function */
                alpham *= a[i]*t[n]/m;

                /* Add to kernel */
                alpha(i, n) += c(i, m)*alpham;

                /* Store subkernel */
                alpham_store(i, n, m) = alpham;
            } /* End for m */
        } /* End for i */
    } /* End for n */

    if(ComputeGradient){
        for(m = 0; m <= M; m++){
            for(n = 0; n < N; n++){
                for(i = 0; i < nz; i++){
                    /* Erlang probability density function */
                    alpham = alpham_store(i, n, m);
                    
                    /* Derivative of density wrt. a */
                    dalpham_a = ((m+1)/a[i] - t[n])*alpham;
                    
                    /* Derivative of density wrt. cm */
                    dalpha(i, n, m) = alpham;

                    /* Derivative of kernel wrt. a */
                    dalpha(i, n, M+1) += c(i, m)*dalpham_a;

                    if(ComputeHessian){
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
    free(alpham_store);
    
    /* Return nothing */
    return;
}

static void evaluate_kernel_faster_without_derivatives(
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
    double *const alpha
    ){
    /* Indices */
    int i, m, n, idx;

    /* Declare variables */
    real alpham, dalpham_a, d2alpham_a;

    for(n = 0; n < N; n++){
        for(i = 0; i < nz; i++){
            /* Initialize subkernels */
            alpham =  a[i]*exponential_function(-a[i]*t[n]);

            /* Initialize kernel */
            alpha(i, n) = c(i, 0)*alpham;

            for(m = 1; m <= M; m++){
                /* Erlang probability density function */
                alpham *= a[i]*t[n]/m;

                /* Add to kernel */
                alpha(i, n) += c(i, m)*alpham;
            } /* End for m */
        } /* End for i */
    } /* End for n */

    /* Return nothing */
    return;
}