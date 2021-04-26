#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe.h"
#include "polynomials.h"
#include "lagrangehalfc_impl.h"
#include "lagrangehalfc_arithmetic.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"

using namespace std;

#define N 1024
#define PRINT_N 12

void print_T_poly(TorusPolynomial *const p)
{
    for (int i = 0; i < PRINT_N; i++)
    {
        printf("%d, ", p->coefsT[i]);
    }
    if (N > PRINT_N) printf("...");
}

void print_I_poly(IntPolynomial *const p)
{
    for (int i = 0; i < PRINT_N; i++)
    {
        printf("%d, ", p->coefs[i]);
    }
    if (N > PRINT_N) printf("...");
}

void print_L_poly(LagrangeHalfCPolynomial *const fp)
{
    LagrangeHalfCPolynomial_IMPL * __fp = (LagrangeHalfCPolynomial_IMPL *)(fp);
    for (int i = 0; i < PRINT_N/2; i++)
    {
        printf("%+.2f%+.2fi, ", __fp->coefsC[i].real(), __fp->coefsC[i].imag());
    }
    if (N > PRINT_N) printf("...");
}



// =============================================================================
//
//  MAIN
//

int32_t main(int32_t argc, char **argv)
{
    IntPolynomial   *p = new_IntPolynomial(N);
    TorusPolynomial *q = new_TorusPolynomial(N);
    TorusPolynomial *r = new_TorusPolynomial(N);

    LagrangeHalfCPolynomial * fp = new_LagrangeHalfCPolynomial(N);
    LagrangeHalfCPolynomial * fq = new_LagrangeHalfCPolynomial(N);
    LagrangeHalfCPolynomial * fr = new_LagrangeHalfCPolynomial(N);

    intPolynomialClear(p);
    p->coefs[0] = 1;
    p->coefs[3] = 2;
    torusPolynomialClear(q);
    q->coefsT[0] = -1;
    q->coefsT[5] = 3;

    //~ FFT_Processor_nayuki_FFNT *const ffnt = fp->proc;

    // FFNT input polynomials
    IntPolynomial_ifft(fp, p);   // here FFNT performs direct transform
    TorusPolynomial_ifft(fq, q);

    // dyadic multiplication
    LagrangeHalfCPolynomialMul(fr, fp, fq);

    // IFFNT output
    TorusPolynomial_fft(r, fr);   // here FFNT performs inverse transform

    // print outputs
    printf("p = (");
    print_I_poly(p);
    printf(")\n");

    printf("q = (");
    print_T_poly(q);
    printf(")\n");

    printf("---------------------------------------------------------------\n");

    printf("fp = (");
    print_L_poly(fp);
    printf(")\n");

    printf("fq = (");
    print_L_poly(fq);
    printf(")\n");

    printf("fr = (");
    print_L_poly(fr);
    printf(")\n");

    printf("---------------------------------------------------------------\n");

    printf("r = (");
    print_T_poly(r);
    printf(")\n");

    // cleanup
    destroy_LagrangeHalfCPolynomial(fr);
    destroy_LagrangeHalfCPolynomial(fq);
    destroy_LagrangeHalfCPolynomial(fp);
    destroy_TorusPolynomial(r);
    destroy_TorusPolynomial(q);
    destroy_IntPolynomial(p);

    return 0;
}
