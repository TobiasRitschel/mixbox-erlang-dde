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
#define c(      i, m      )  (      c[(i) + nz*(m)                            ])
#define alpha(  i, n      )  (  alpha[(i) + nz*(n)                            ])
#define dalpha( i, n, j   )  ( dalpha[(i) + nz*(n) + nz*N*(j)                 ])
#define d2alpha(i, n, j, k)  (d2alpha[(i) + nz*(n) + nz*N*(j) + nz*N*(M+2)*(k)])

static void evaluate_kernel(
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
    double bm, alpham, dalpham_a, d2alpham_a;

    for(i = 0; i < nz; i++){
        /* Reset normalization factor */
        bm = 1.0;

        /* Loop over Erlang order */
        for(m = 0; m <= M; m++){
            /* Normalization factor */
            bm *= a[i]/MAX(1, m);

            for(n = 0; n < N; n++){
                /* Erlang probability density function */
                alpham = bm*pow(t[n], m)*exp(-a[i]*t[n]);

                /* Add to kernel */
                alpha(i, n) += c(i, m)*alpham;

                if(ComputeGradient){
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
                    } /* End if ComputeHessian */
                } /* End if ComputeGradient */
            } /* End for n */
        } /* End for m */
    } /* End for i */
    
    /* Return nothing */
    return;
}