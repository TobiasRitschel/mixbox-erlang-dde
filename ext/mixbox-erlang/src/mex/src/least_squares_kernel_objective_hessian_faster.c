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
#define c(        i, m)  (   c     [(i) +  nz  *(m)])
#define hess(     m, l)  (hess     [(m) + (M+2)*(l)])
#define alphameas(i, n)  (alphameas[(i) +  nz  *(n)])

static void least_squares_kernel_objective_hessian_faster(
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
    real alpha, dalpha, d2alpha;
    real alpham, dalpham, d2alpham;
    double e, dt;

    /* Do not compute Hessian in this routine */
    const int ComputeGradient = true;
    const int ComputeHessian  = true;

    /* Allocate memory */
    real *const  alpham_store = (real *) calloc(M+1, sizeof(real));
    real *const dalpham_store = (real *) calloc(M+1, sizeof(real));

    /* Loop over parameters */
    for(n = 0; n < N-1; n++){
        for(i = 0; i < nz; i++){
            /* Initialize subkernels */
            alpham   =  a[i]*exponential_function(-a[i]*t[n]);
            dalpham  =                           (1.0/a[i] - t[n])* alpham;
            d2alpham = -1.0/(a[i]*a[i])*alpham + (1.0/a[i] - t[n])*dalpham;

            /* Initialize kernel */
            alpha   = c(i, 0)*  alpham;
            dalpha  = c(i, 0)* dalpham;
            d2alpha = c(i, 0)*d2alpham;

            /* Store values */
            alpham_store [0] =  alpham;
            dalpham_store[0] = dalpham;

            for(m = 1; m <= M; m++){
                /* Erlang probability density function */
                alpham *= a[i]*t[n]/m;

                /* Derivative of density wrt. rate parameter */
                dalpham = ((m+1)/a[i] - t[n])*alpham;
        
                /* Derivative of density wrt. rate parameter */
                d2alpham = -(m+1)/(a[i]*a[i])*alpham + ((m+1)/a[i] - t[n])*dalpham;

                /* Add to approximate kernel */
                alpha += c(i, m)*alpham;

                /* Add to derivative approximate kernel */
                dalpha += c(i, m)*dalpham;
                
                /* Add to derivative approximate kernel */
                d2alpha += c(i, m)*d2alpham;

                /* Store values */
                alpham_store [m] =  alpham;
                dalpham_store[m] = dalpham;
            } /* End for m */

            /* Time step size */
            dt = t[n+1] - t[n];

            /* Compute error */
            e = alphameas(i, n) - alpha;

            /* Add to derivative */
            hess(M+1, M+1) += (dalpha*dalpha - e*d2alpha)*dt;

            for(l = 0; l <= M; l++){
                for(m = l; m <= M; m++){
                    /* Mixed derivative wrt. coefficients */
                    hess(m, l) += alpham_store[m]*alpham_store[l]*dt;
                } /* End for m */

                /* Mixed derivative wrt. coefficients and rate parameter */
                hess(l, M+1) += (dalpha*alpham_store[l] - e*dalpham_store[l])*dt;
            } /* End for l */
        } /* End for i */
    } /* End for n */

    /* Copy entries*/
    for(m = 0; m <= M; m++){
        for(l = m+1; l <= M; l++){
            /* Mixed derivative wrt. coefficients */
            hess(m, l) = hess(l, m);
        } /* End for m */

        /* Mixed derivative wrt. coefficients and rate parameter */
        hess(M+1, m) = hess(m, M+1);
    } /* End for l */

    /* Free memory */
    free( alpham_store);
    free(dalpham_store);

    /* Return nothing */
    return;
}