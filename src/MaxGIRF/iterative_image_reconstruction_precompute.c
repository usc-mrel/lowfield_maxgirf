/*=======================================================================*/
/*                                                                       */
/*  SOURCE_FILE:    ITERATIVE_IMAGE_RECONSTRUCTION_PRECOMPUTE.C          */
/*                                                                       */
/*  Copyright 2020: Nam Gyun Lee                                         */
/*  Author: Nam Gyun Lee                                                 */
/*  USC (University of Southern California)                              */
/*  namgyunl@usc.edu, ggang56@gmail.com (preferred)                      */
/*=======================================================================*/
/*=======================================================================*/
/*  I N C L U D E S                                                      */
/*=======================================================================*/
/*-----------------------------------------------------------------------*/
/*  Common                                                               */
/*-----------------------------------------------------------------------*/
#define _USE_MATH_DEFINES // To use MATH Constants
#include <time.h>
#include <math.h>
#include <xmmintrin.h> // For _mm_malloc
#include <immintrin.h> // For Intel intrinsics
#include <string.h>    // For memset

/*-----------------------------------------------------------------------*/
/*  MATLAB                                                               */
/*-----------------------------------------------------------------------*/
/* Macro to define the correct function prototype */
#if defined (_WIN32) || defined(__APPLE__)
#define FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID 1
#endif

#include "mex.h"
#include "blas.h"

/*=======================================================================*/
/*  G L O B A L   R E F E R E N C E S                                    */
/*=======================================================================*/
/*=======================================================================*/
/*  G L O B A L   D E F I N I T I O N S                                  */
/*=======================================================================*/
/*=======================================================================*/
/*  L O C A L   S Y M B O L   D E F I N I T I O N S                      */
/*=======================================================================*/
/*=======================================================================*/
/*  L O C A L   D A T A   D E F I N I T I O N S                          */
/*=======================================================================*/
/*=======================================================================*/
/*  L O C A L   F U N C T I O N   P R O T O T Y P E S                    */
/*=======================================================================*/
/*=======================================================================*/
/*  F U N C T I O N   P R O T O T Y P E S                                */
/*=======================================================================*/
void iterative_image_reconstruction_precompute(mxComplexSingle *d,
                                               mxComplexSingle *E,
                                               mxComplexSingle *csm,
                                               float beta,
                                               size_t maxiter,
                                               size_t Nk,
                                               size_t Nc,
                                               size_t N,
                                               size_t Ni,
                                               mxComplexSingle *m_tilde,
                                               double *delta);

/* complex, Hadamard product (basic) */
void chad(const char *trans, const size_t *n, const float *alpha, const float *x, const float *y, const float *beta, float *z);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const size_t *d_size;  /* pointer to the dimensions array of d     */
    size_t Nk;             /* number of k-space samples per spiral arm */
    size_t Ni;             /* number of spiral interleaves             */
    size_t Nc;             /* number of coils                          */
    size_t N;              /* number of reconstructed pixels           */

    /* cg parameters */
    float beta;
    size_t maxiter;

    /* pointers to interleaved complex data of a MATLAB array (mxArray) */
    mxComplexSingle *d       = NULL; /* Nk x Ni x Nc */
    mxComplexSingle *E       = NULL; /* Nk x N  x Nc */
    mxComplexSingle *csm     = NULL; /* N1 x N2 x Nc */
    mxComplexSingle *m_tilde = NULL; /* N x 1        */
    double          *delta   = NULL; /* maxiter x 1  */

    /* help text in the event of an input error */
    const char help_txt[] = 
"iterative_image_reconstruction_precompute -- higher order reconstruction using a conjugate gradient solver \n\
                                              with a precomputed E matrix                                   \n\
 Usage                                                                                                      \n\
   m_tilde = iterative_image_reconstruction_precompute(d, E, csm, beta, maxiter)                            \n\
 Inputs                                                                                                     \n\
   d          k-space data                      (float , complex) (Nk x Ni x Nc)                            \n\
   E          encoding matrix                   (float , complex) (Nk x N  x Ni)                            \n\
   csm        coil sensitivity maps             (float , complex) (N1 x N2 x Nc)                            \n\
   beta       Tikhonov regularization parameter (double, real   ) (scalar)                                  \n\
   maxiter    maximum number of CG iterations   (double, real   ) (scalar)                                  \n\
 Outputs                                                                                                    \n\
   m_tilde    reconstruction                                 (float , complex) (N x maxiter)                \n\
   delta      ||E^H * E * m - E^H * d||_l2 / ||E^H * d||_l2  (double, real   ) (maxiter x 1)                \n\
                                                                                                            \n\
 Description                                                                                                \n\
                                                                                                            \n\
 See Also                                                                                                   \n\
                                                                                                            \n\
Copyright (c) 2020. Nam Gyun Lee                                                                            \n\
Comments? e-mail namgyunl@usc.edu (or ggang56@gmail.com)                                                    \n";

    /*-------------------------------------------------------------------*/
    /* Check for proper number of arguments                              */
    /*-------------------------------------------------------------------*/
    if (nrhs != 5)
    {
        printf("%s", help_txt);
        mexErrMsgTxt("5 inputs are required.");
    }
    if (nlhs != 2)
    {
        printf("%s", help_txt);
        mexErrMsgTxt("2 outputs are required.");
    }

    /*-------------------------------------------------------------------*/
    /* Check the data type of inputs                                     */
    /*-------------------------------------------------------------------*/
    if (!mxIsSingle(prhs[0]) || !mxIsComplex(prhs[0]))
    {
        printf("%s", help_txt);
        mexErrMsgTxt("The 1st input, d, must be single and complex.");
    }
    if (!mxIsSingle(prhs[1]) || !mxIsComplex(prhs[1]))
    {
        printf("%s", help_txt);
        mexErrMsgTxt("The 2nd input, E, must be single and complex.");
    }
    if (!mxIsSingle(prhs[2]) || !mxIsComplex(prhs[2]))
    {
        printf("%s", help_txt);
        mexErrMsgTxt("The 3rd input, s, must be single and complex.");
    }
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]))
    {
        printf("%s", help_txt);
        mexErrMsgTxt("The 4th input, beta, must be double and real.");
    }
    if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]))
    {
        printf("%s", help_txt);
        mexErrMsgTxt("The 5th input, maxiter, must be double and real.");
    }

    /*-------------------------------------------------------------------*/
    /* Set size, dimension related parameters                            */
    /*-------------------------------------------------------------------*/
    d_size = mxGetDimensions(prhs[0]);
    Nk = d_size[0];

    if (mxGetNumberOfDimensions(prhs[0]) == 1)
    {
        Ni = 1;
        Nc = 1;
    } else if (mxGetNumberOfDimensions(prhs[0]) == 2)
    {
        Ni = d_size[1];
        Nc = 1;
    } else if (mxGetNumberOfDimensions(prhs[0]) == 3)
    {
        Ni = d_size[1];
        Nc = d_size[2];
    }
    N = mxGetN(prhs[1]) / Ni;

    /*-------------------------------------------------------------------*/
    /* Set a pointer to the interleaved complex data of a MATLAB array   */
    /*-------------------------------------------------------------------*/
    d   = mxGetComplexSingles(prhs[0]);
    E   = mxGetComplexSingles(prhs[1]);
    csm = mxGetComplexSingles(prhs[2]);

    /*-------------------------------------------------------------------*/
    /* Dereference a value directly from a pointer to mxArray            */
    /*-------------------------------------------------------------------*/
    beta    = (float)  *mxGetDoubles(prhs[3]);
    maxiter = (size_t) *mxGetDoubles(prhs[4]);

    #ifdef DEBUG
    printf("Nk = %lu, Ni = %lu, Nc = %lu, N = %lu\n", Nk, Ni, Nc, N);
    printf("beta = %g, maxiter = %lu\n", beta, maxiter);
    #endif

    /*-------------------------------------------------------------------*/
    /* Create output MATLAB arrays for the return argument               */
    /*-------------------------------------------------------------------*/
    plhs[0] = mxCreateNumericMatrix(N, maxiter, mxSINGLE_CLASS, mxCOMPLEX);
    plhs[1] = mxCreateNumericMatrix(maxiter, 1, mxDOUBLE_CLASS, mxREAL);
    m_tilde = mxGetComplexSingles(plhs[0]);
    delta   = mxGetDoubles(plhs[1]);

    /*-------------------------------------------------------------------*/
    /* Conjugate gradient algorithm                                      */
    /*-------------------------------------------------------------------*/
    iterative_image_reconstruction_precompute(d, E, csm, beta, maxiter, Nk, Nc, N, Ni, m_tilde, delta);
}

/*=======================================================================*/
/*  iterative_image_reconstruction_precompute                            */
/*=======================================================================*/
void iterative_image_reconstruction_precompute(mxComplexSingle *d,
                                               mxComplexSingle *E,
                                               mxComplexSingle *csm,
                                               float beta,
                                               size_t maxiter,
                                               size_t Nk,
                                               size_t Nc,
                                               size_t N,
                                               size_t Ni,
                                               mxComplexSingle *m_tilde,
                                               double *delta)
{
    mxComplexSingle *p       = NULL;
    mxComplexSingle *Ehv     = NULL;
    mxComplexSingle *Sp      = NULL;
    mxComplexSingle *ESp     = NULL;
    mxComplexSingle *ShEhESp = NULL;
    mxComplexSingle *EhESp   = NULL;
    mxComplexSingle *q       = NULL;
    mxComplexSingle *r       = NULL;
    mxComplexSingle *r_old   = NULL;

    size_t c;
    size_t i;
    size_t iter;
    size_t blocksize_N;  /* size of dynamically allocated N elements */
    size_t blocksize_Nk; /* size of dynamically allocated Nk elements */

    complex cbeta   = {beta,0.0};
    complex aHa     = {0.0,0.0};
    complex rHr     = {0.0,0.0};
    complex pHq     = {0.0,0.0};
    complex rHr_old = {0.0,0.0};
    complex ca      = {0.0,0.0};

    double numerator_real;
    double numerator_imag;
    double denominator;

    /* constant variables for BLAS call */
    complex cone  = {1.0,0.0};
    complex czero = {0.0,0.0};
    ptrdiff_t one = 1;

    clock_t start, stop;
    double elapsed_time;

    #ifdef DEBUG
    printf("DEBUG: beta = %g, maxiter = %lu\n", beta, maxiter);
    printf("DEBUG: cbeta.r = %g, cbeta.i = %g\n", cbeta.r, cbeta.i);
    printf("DEBUG: cone.r = %g, cone.i = %g\n", cone.r, cone.i);
    #endif

    /*-------------------------------------------------------------------*/
    /* Allocate aligned memory (32-byte aligned) on the CPU              */
    /*-------------------------------------------------------------------*/
    blocksize_N  = sizeof(mxComplexSingle) * N;
    blocksize_Nk = sizeof(mxComplexSingle) * Nk;
    p       = (mxComplexSingle *) _mm_malloc(blocksize_N , 32);
    Ehv     = (mxComplexSingle *) _mm_malloc(blocksize_N , 32);
    Sp      = (mxComplexSingle *) _mm_malloc(blocksize_N , 32);
    ESp     = (mxComplexSingle *) _mm_malloc(blocksize_Nk, 32);
    EhESp   = (mxComplexSingle *) _mm_malloc(blocksize_N , 32);
    ShEhESp = (mxComplexSingle *) _mm_malloc(blocksize_N , 32);
    q       = (mxComplexSingle *) _mm_malloc(blocksize_N , 32);
    r       = (mxComplexSingle *) _mm_malloc(blocksize_N , 32);
    r_old   = (mxComplexSingle *) _mm_malloc(blocksize_N , 32);

    /*-------------------------------------------------------------------*/
    /* Initialize dynamically allocated arrays to 0                      */
    /*-------------------------------------------------------------------*/
    memset(p      , '\0', blocksize_N);
    memset(Ehv    , '\0', blocksize_N);
    memset(Sp     , '\0', blocksize_N);
    memset(ESp    , '\0', blocksize_Nk);
    memset(EhESp  , '\0', blocksize_N);
    memset(ShEhESp, '\0', blocksize_N);
    memset(q      , '\0', blocksize_N);
    memset(r      , '\0', blocksize_N);
    memset(r_old  , '\0', blocksize_N);

    /*-------------------------------------------------------------------*/
    /* Calculate p = E^H(d) (N x 1)                                      */
    /* A_{i,c}^H * d_{i,c} = Sc^H * Ei^H * d_{i,c}                       */
    /*-------------------------------------------------------------------*/
    for (c = 0; c < Nc; c++)
    {
        memset(Ehv, '\0', blocksize_N);

        for (i = 0; i < Ni; i++)
        {
            start = clock();
            printf("Calculating p = E^H(d) (c=%2lu/%lu, i=%2lu/%2lu)... ", c+1, Nc, i+1, Ni);

            /*-----------------------------------------------------------*/
            /* cgemv: y <- alpha*A^H*x + beta*y                          */
            /* Ehv = Ehv + E(:,:,i)' * d(:,i,c);                         */
            /*-----------------------------------------------------------*/
            cgemv("C", &Nk, &N, (float *)&cone, (float *)(E + i * Nk * N), 
                  &Nk, (float *)(d + i * Nk + c * Nk * Ni), &one, (float *)&cone, (float *)Ehv, &one);

            stop = clock();
            elapsed_time = (double) (stop - start) / CLOCKS_PER_SEC; /* [sec] */
            printf("done! (%6.4f sec)\n", elapsed_time);
        }

        /*---------------------------------------------------------------*/
        /* z <- alpha * conj(x) * y + beta * z , i = 1,...,n             */
        /* p = p + reshape(conj(s(:,:,c)), [N 1]) .* Ehv;                */
        /*---------------------------------------------------------------*/
        chad("C", &N, (float *)&cone, (float *)(csm + c * N), (float *)Ehv, (float *)&cone, (float *)p);
    }

    /*-------------------------------------------------------------------*/
    /* r = p (N x 1)                                                     */
    /*-------------------------------------------------------------------*/
    ccopy(&N, (float *)p, &one, (float *)r, &one);

    /*-------------------------------------------------------------------*/
    /* Calculate a^H * a                                                 */
    /*-------------------------------------------------------------------*/
    #ifndef FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID
    aHa = cdotc(&N, (float *)p, &one, (float *)p, &one);
    #else
    cdotc((complex *)&aHa, &N, (float *)p, &one, (float *)p, &one);
    #endif

    /*-------------------------------------------------------------------*/
    /* main loop                                                         */
    /*-------------------------------------------------------------------*/
    for (iter = 0; iter < maxiter; iter++)
    {
        printf("MaxGIRF reconstruction using a conjugate gradient solver, iteration = %lu:\n", iter+1);

        /*---------------------------------------------------------------*/
        /* Calculate delta = r^H * r / a^H * a                           */
        /*---------------------------------------------------------------*/
        #ifndef FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID
        rHr = cdotc(&N, (float *)r, &one, (float *)r, &one);
        #else
        cdotc((complex *)&rHr, &N, (float *)r, &one, (float *)r, &one);
        #endif
        delta[iter] = (double) rHr.r / (double) aHa.r;

        /*---------------------------------------------------------------*/
        /* Set arrays to 0                                               */
        /*---------------------------------------------------------------*/
        memset(Sp     , '\0', blocksize_N);
        memset(ESp    , '\0', blocksize_Nk);
        memset(EhESp  , '\0', blocksize_N);
        memset(ShEhESp, '\0', blocksize_N);

        /*---------------------------------------------------------------*/
        /* Calculate ShEhESp (N x 1)                                     */
        /*                                                               */
        /* (A_{i,c}^H * A_{i,c}^H) * p = Sc^H * Ei^H * Ei * Sc * p       */
        /* = [           ][     ][           ][           ][ ]           */
        /*   [           ][     ][     Ei    ][           ][ ]           */
        /*   [    Sc^H   ][ Ei^H][           ][    Sc     ][p]           */ 
        /*   [           ][     ]             [           ][ ]           */
        /*   [           ][     ]             [           ][ ]           */
        /*---------------------------------------------------------------*/
        for (i = 0; i < Ni; i++)
        {
            for (c = 0; c < Nc; c++)
            {
                start = clock();
                printf("(it=%2lu/%2lu): Calculating ShEhESp (c=%2lu/%2lu, i=%2lu/%2lu)... ", iter+1, maxiter, c+1, Nc, i+1, Ni);\

                /*-------------------------------------------------------*/
                /* z[i] <- alpha*x[i]*y[i] + beta*z[i], i = 1,...,n      */
                /* Sp = reshape(s(:,:,c), [N 1]) .* p;                   */
                /*-------------------------------------------------------*/
                chad("N", &N, (float *)&cone, (float *)(csm + c * N), (float *)p, (float *)&czero, (float *)Sp);

                /*-------------------------------------------------------*/
                /* cgemv: y <- alpha*A*x + beta*y                        */
                /* ESp = E(:,:,i) * Sp;                                  */
                /* ESp: Nk x 1, E(:,:,i) in Nk x N, Sp in N x 1          */
                /*-------------------------------------------------------*/
                cgemv("N", &Nk, &N, (float *)&cone, (float *)(E + i * Nk * N),
                      &Nk, (float *)Sp, &one, (float *)&czero, (float *)ESp, &one);

                /*-------------------------------------------------------*/
                /* cgemv: y <- alpha*A^H*x + beta*y                      */
                /* EhESp = E(:,:,i)' * ESp;                              */
                /* EhESp in N x 1, E(:,:,i)' in N x Nk, ESp in Nk x 1    */
                /*-------------------------------------------------------*/
                cgemv("C", &Nk, &N, (float *)&cone, (float *)(E + i * Nk * N),
                      &Nk, (float *)ESp, &one, (float *)&czero, (float *)EhESp, &one);

                /*--------------------------------------------------------------*/
                /* z[i] <- alpha*conj(x[i])*y[i] + beta*z[i], i = 1,...,n       */
                /* ShEhESp = ShEhESp + reshape(conj(s(:,:,c)), [N 1]) .* EhESp; */
                /*--------------------------------------------------------------*/
                chad("C", &N, (float *)&cone, (float *)(csm + c * N), (float *)EhESp, (float *)&cone, (float *)ShEhESp);
                //chad_avx("C", N, (float *)(csm + c * N), (float *)EhESp, (float *)&cone, (float *)ShEhESp);

                stop = clock();
                elapsed_time = (double) (stop - start) / CLOCKS_PER_SEC; /* [sec] */
                printf("done! (%6.4f sec)\n", elapsed_time);
            }
        }

        /*---------------------------------------------------------------*/
        /* Calculate q (N x 1)                                           */
        /* q = ShEhESp + beta * p;                                       */
        /* ccopy: y <- x                                                 */
        /*        q <- ShEhESp                                           */
        /* caxpy: y <- a*x + y                                           */
        /*        q <- beta * p + q                                      */
        /*---------------------------------------------------------------*/
        ccopy(&N, (float *)ShEhESp, &one, (float *)q, &one);
        caxpy(&N, (float *)&cbeta, (float *)p, &one, (float *)q, &one);

        /*---------------------------------------------------------------*/
        /* Calculate m_tilde (N x 1)                                     */
        /* m_tilde = m_tilde + ((r' * r) / (p' * q)) * p;                */
        /*---------------------------------------------------------------*/
        #ifndef FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID
        pHq = cdotc(&N, (float *)p, &one, (float *)q, &one);
        #else
        cdotc((complex *)&pHq, &N, (float *)p, &one, (float *)q, &one);
        #endif

        /*---------------------------------------------------------------*/
        /* a = rHr.r, b = rHr.i, c = pHq.r, d = pHq.i                    */
        /* (a + jb) / (c + jd) = (a + jb)(c - jd) / (c^2 + d^2)          */
        /* = ((ac + bd) + j(bc - ad)) / (c^2 + d^2)                      */
        /*---------------------------------------------------------------*/
        numerator_real = ((double) rHr.r * (double) pHq.r + (double) rHr.i * (double) pHq.i);
        numerator_imag = ((double) rHr.i * (double) pHq.r - (double) rHr.r * (double) pHq.i);
        denominator    = ((double) pHq.r * (double) pHq.r + (double) pHq.i * (double) pHq.i);
        ca.r = (float) (numerator_real / denominator);
        ca.i = (float) (numerator_imag / denominator);

        caxpy(&N, (float *)&ca, (float *)p, &one, (float *)(m_tilde + (maxiter - 1) * N), &one);

        /*---------------------------------------------------------------*/
        /* Save the intermediate image                                   */
        /*---------------------------------------------------------------*/
        ccopy(&N, (float *)(m_tilde + (maxiter - 1) * N), &one, (float *)(m_tilde + iter * N), &one);

        /*---------------------------------------------------------------*/
        /* Calculate r_old (N x 1)                                       */
        /* r_old = r;                                                    */
        /*---------------------------------------------------------------*/
        ccopy(&N, (float *)r, &one, (float *)r_old, &one);

        /*---------------------------------------------------------------*/
        /* Calculate r (N x 1)                                           */
        /* r  = r - ((r' * r) / (p' * q)) * q;                           */
        /*---------------------------------------------------------------*/
        ca.r = -ca.r;
        ca.i = -ca.i;
        caxpy(&N, (float *)&ca, (float *)q, &one, (float *)r, &one);

        /*---------------------------------------------------------------*/
        /* Calculate p (N x 1)                                           */
        /* p = r + ((r' * r) / (r_old' * r_old)) * p;                    */
        /* cscal: x <- a*x                                               */
        /*        p <- ((r' * r) / (r_old' * r_old)) * p                 */
        /* caxpy: y <- a*x + y                                           */
        /*        p <- r + p                                             */
        /*---------------------------------------------------------------*/
        #ifndef FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID
        rHr     = cdotc(&N, (float *)r    , &one, (float *)r    , &one);
        rHr_old = cdotc(&N, (float *)r_old, &one, (float *)r_old, &one);
        #else
        cdotc((complex *)&rHr    , &N, (float *)r    , &one, (float *)r    , &one);
        cdotc((complex *)&rHr_old, &N, (float *)r_old, &one, (float *)r_old, &one);
        #endif
        ca.r = rHr.r / rHr_old.r;
        ca.i = 0;
        cscal(&N, (float *)&ca, (float *)p, &one);
        caxpy(&N, (float *)&cone, (float *)r, &one, (float *)p, &one);
    }

    /*-------------------------------------------------------------------*/
    /* Free CPU memory                                                   */
    /*-------------------------------------------------------------------*/
    _mm_free(p);
    _mm_free(Ehv);
    _mm_free(Sp);
    _mm_free(ESp);
    _mm_free(ShEhESp);
    _mm_free(EhESp);
    _mm_free(q);
    _mm_free(r);
    _mm_free(r_old);
}

/*-----------------------------------------------------------------------*/
/* Compute the Hadamard product of two vectors x and y, and store the    */
/* result in vector z.                                                   */
/* z[i] <- alpha * op(x[i]) * y[i] + beta * z[i] , i = 1,...,n           */
/*-----------------------------------------------------------------------*/
void chad(const char *trans, const size_t *n, const float *alpha, const float *x, const float *y, const float *beta, float *z)
{
    size_t i;
    size_t offset;
    float rx;
    float ix;
    float ry;
    float iy;
    float rz;
    float iz;
    float rxy;
    float ixy;

    size_t n_ = *n;

    float ra = *(alpha);
    float ia = *(alpha + 1);
    float rb = *(beta);
    float ib = *(beta + 1);

    /*-------------------------------------------------------------------*/
    /* (ra + j*ia) * (rx + j*ix) * (ry + j*iy)                           */
    /* = (ra + j*ia) * [(rx * ry - ix * iy) + j*(rx * iy + ix * ry)]     */
    /* = ra*(rx * ry - ix * iy) - ia*(rx * iy + ix * ry)                 */
    /*   + j*( ra*(rx * iy + ix * ry) + ia*(rx * ry - ix * iy) )         */
    /*-------------------------------------------------------------------*/
    if (strcmp(trans, "N") == 0)
    {
        for (i = 0; i < n_; i++)
        {
            offset = 2 * i;
            rx =  x[offset];
            ix =  x[offset + 1];
            ry =  y[offset];
            iy =  y[offset + 1];
            rz =  z[offset];
            iz =  z[offset + 1];

            rxy = rx * ry - ix * iy;
            ixy = rx * iy + ix * ry;

            z[offset]     = ra * rxy - ia * ixy + (rb * rz - ib * iz);
            z[offset + 1] = ra * ixy + ia * rxy + (rb * iz + ib * rz);
        }
    } else if (strcmp(trans, "C") == 0)
    {
        for (i = 0; i < n_; i++)
        {
            offset = 2 * i;
            rx =  x[offset];
            ix = -x[offset + 1];
            ry =  y[offset];
            iy =  y[offset + 1];
            rz =  z[offset];
            iz =  z[offset + 1];

            rxy = rx * ry - ix * iy;
            ixy = rx * iy + ix * ry;

            z[offset]     = ra * rxy - ia * ixy + (rb * rz - ib * iz);
            z[offset + 1] = ra * ixy + ia * rxy + (rb * iz + ib * rz);
        }
    }
}

/*=======================================================================*/
/*  H I S T O R Y                                                        */
/*=======================================================================*/
/* 2020-04-14 Namgyun Lee                                                */
/*      Introducing source-file into ITERATIVE_IMAGE_RECONSTRUCTION_BLAS */
/* 2020-04-14 Namgyun Lee                                                */
/*      Finished the initial version                                     */
