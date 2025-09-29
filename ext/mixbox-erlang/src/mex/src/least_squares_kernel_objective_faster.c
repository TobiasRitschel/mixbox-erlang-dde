/* Headers */
#include <math.h>
#include <string.h>

#ifndef LONG_KERNEL_COMPUTATIONS
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

/* Convenience functions for matrix elements */
#define c(        i, m) (        c[(i) + nz*(m)])
#define alphameas(i, n) (alphameas[(i) + nz*(n)])

static void least_squares_kernel_objective_faster(
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
    real alpha, dalpha;
    real alpham, dalpham;
    double e, de, dt;

    /* Allocate memory */
    real *const alpham_store = (real *) calloc(M+1, sizeof(real));

    /* Loop over parameters */
    for(n = 0; n < N-1; n++){
        for(i = 0; i < nz; i++){
            /* Initialize subkernels */
            alpham = a[i]*exponential_function(-a[i]*t[n]);

            /* Initialize kernel */
            alpha = c(i, 0)*alpham;

            if(ComputeGradient){
                /* Derivative of subkernel */
                dalpham = (1.0/a[i] - t[n])*alpham;
                
                /* Derivative of kernel */
                dalpha = c(i, 0)*dalpham;

                /* Store values */
                alpham_store[0] = alpham;
            } /* End if ComputeGradient */

            for(m = 1; m <= M; m++){
                /* Erlang probability density function */
                alpham *= a[i]*t[n]/m;

                /* Add to approximate kernel */
                alpha += c(i, m)*alpham;

                if(ComputeGradient){
                    /* Derivative of density wrt. rate parameter */
                    dalpham = ((m+1)/a[i] - t[n])*alpham;
            
                    /* Add to derivative approximate kernel */
                    dalpha += c(i, m)*dalpham;

                    /* Store values */
                    alpham_store[m] = alpham;
                } /* End if ComputeGradient */
            } /* End for m */

            /* Time step size */
            dt = t[n+1] - t[n];

            /* Compute error */
            e = alphameas(i, n) - alpha;

            /* Add to objective function */
            *phi += e*e*dt;

            if(ComputeGradient){
                for(m = 0; m <= M; m++){
                    /* Derivative of error wrt. rate parameter */
                    de = -alpham_store[m];
                    
                    /* Gradient wrt. coefficients */
                    grad[m] += e*de*dt;
                } /* End for m */
                
                /* Derivative of error wrt. rate parameter */
                de = -dalpha;

                /* Gradient wrt. rate parameter */
                grad[M+1] += e*de*dt;
            } /* End if ComputeGradient */
        } /* End for i */
    } /* End for n */

    /* Scale objective function */
    *phi *= 0.5;

    /* Free memory */
    free(alpham_store);

    /* Return nothing */
    return;
}