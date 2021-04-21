#include <complex>
#include <polynomials.h>
#include "lagrangehalfc_impl.h"
#include "ffnt.h"
#include <cassert>
#include <cmath>

FFT_Processor_nayuki_FFNT::FFT_Processor_nayuki_FFNT(const int32_t N)
    : _2N(2*N), N(N), Ns2(N/2)
{
    // stores inputs and (intermediate) results
    real_inout = (double *)malloc(sizeof(double) * N);
    imag_inout = (double *)malloc(sizeof(double) * N);

    ffnt_tables_2N   = (void *)fft_init(_2N);
    ffnt_tables_N_2  = (void *)fft_init(Ns2);
    iffnt_tables_N_2 = (void *)ifft_init(Ns2);

    //TODO figure out what is the representation of ??? - 1 in FFNT
    omegaxminus1 = (cplx *)malloc(sizeof(cplx) * _2N);
    for (int32_t x = 0; x < _2N; x++)
    {
        omegaxminus1[x]=cplx(cos(x*M_PI/N)-1., sin(x*M_PI/N));
        // instead of cos(x*M_PI/N)-1. + sin(x*M_PI/N) * 1i
        // exp(i.x.pi/N)-1
    }
}

void FFT_Processor_nayuki_FFNT::execute_reverse_int(cplx * res,         // LagrangeHalfCPolynomial_IMPL.coefsC
                                                    const int32_t * a)  // IntPolynomial.coefs
{
    double * res_dbl = (double *)res;

    // cast integral input to double
    for (int32_t i = 0; i < N; i++)
        real_inout[i] = (double)(a[i]);

    // direct !!
    ffnt_transform((FftTables *)ffnt_tables_2N, (FftTables *)ffnt_tables_N_2, real_inout, imag_inout);

    //TODO this copying is completely useless
    for (int32_t i = 0; i < Ns2; i++)
    {
        res_dbl[2*i]   = real_inout[i];
        res_dbl[2*i+1] = imag_inout[i];
    }
}

void FFT_Processor_nayuki_FFNT::execute_reverse_torus32(cplx * res,         // LagrangeHalfCPolynomial_IMPL.coefsC
                                                        const Torus32 * a)  // TorusPolynomial.coefsT (int32_t)
{
    const double _2pm32 = (double)1.0 / (INT64_C(1) << 32);
    double * res_dbl = (double *)res;

    // scale Torus input
    for (int32_t i = 0; i < N; i++)
        real_inout[i] = _2pm32 * a[i];

    // direct !!
    ffnt_transform((FftTables *)ffnt_tables_2N, (FftTables *)ffnt_tables_N_2, real_inout, imag_inout);

    //TODO this copying is completely useless
    for (int32_t i = 0; i < Ns2; i++)
    {
        res_dbl[2*i] = real_inout[i];
        res_dbl[2*i+1] = imag_inout[i];
    }
}

void FFT_Processor_nayuki_FFNT::execute_direct_torus32(Torus32 * res,    // TorusPolynomial.coefsT
                                                       const cplx * a)   // LagrangeHalfCPolynomial_IMPL.coefsC
{
    const double _2p32 = (double)(INT64_C(1) << 32);
    double * a_dbl = (double *)a;

    // init FFNT arys
    for (int32_t i = 0; i < Ns2; i++)
    {
        real_inout[i] = a_dbl[2*i];
        imag_inout[i] = a_dbl[2*i+1];
    }

    // reverse !!
    iffnt_transform((FftTables *)ffnt_tables_2N, (FftTables *)iffnt_tables_N_2, real_inout, imag_inout);

    // scale result back to Torus
    for (int32_t i = 0; i < N; i++)
        res[i] = (Torus32)round(real_inout[i] * _2p32);
}

// OK
FFT_Processor_nayuki_FFNT::~FFT_Processor_nayuki_FFNT()
{
    tables_destroy(ffnt_tables_2N);
    tables_destroy(ffnt_tables_N_2);
    tables_destroy(iffnt_tables_N_2);

    free(real_inout);
    free(imag_inout);
    free(omegaxminus1);
}

// add other degrees
thread_local FFT_Processor_nayuki_FFNT fp1024_nayuki_FFNT(1024);

/**
 * Top-level FFT functions
 */
EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial * result,
                               const IntPolynomial * p)
{
    LagrangeHalfCPolynomial_IMPL * r = (LagrangeHalfCPolynomial_IMPL *) result;
    fp1024_nayuki_FFNT.execute_reverse_int(r->coefsC, p->coefs);
}

EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial * result,
                                 const TorusPolynomial * p)
{
    LagrangeHalfCPolynomial_IMPL * r = (LagrangeHalfCPolynomial_IMPL *) result;
    fp1024_nayuki_FFNT.execute_reverse_torus32(r->coefsC, p->coefsT);
}

EXPORT void TorusPolynomial_fft(TorusPolynomial * result,
                                const LagrangeHalfCPolynomial * p)
{
    LagrangeHalfCPolynomial_IMPL * r = (LagrangeHalfCPolynomial_IMPL *) p;
    fp1024_nayuki_FFNT.execute_direct_torus32(result->coefsT, r->coefsC);
}
