#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"

#define PI 4

using namespace std;

void die_soon(string message)
{
    cerr << "(!) " << message << "\n    Aborting ..." << endl;
    abort();
}


// -----------------------------------------------------------------------------
//  En/Decryption
//
void paral_sym_encr_priv(LweSample *ct,
                         const int32_t message,
                         const TFheGateBootstrappingSecretKeySet *sk)
{
    Torus32 _1s16 = modSwitchToTorus32(1, 1 << PI);

    // scale to [-8/16, 7/16]
    //                         ___/ 8 \___     _/ mask 1111 \_     ___/ 8 \___
    Torus32 mu = (((message + (1 << (PI-1))) & ((1 << PI) - 1)) - (1 << (PI-1))) * _1s16;
    double alpha = sk->params->in_out_params->alpha_min; //TODO: specify noise
    lweSymEncrypt(ct, mu, alpha, sk->lwe_key);
}

void paral_sym_encr(LweSample *ct,
                    const int32_t message,
                    const TFheGateBootstrappingSecretKeySet *sk)
{
    if ((message < -2) || (2 < message))
        die_soon("Out of the alphabet A = [-2 .. 2].");

    paral_sym_encr_priv(ct, message, sk);
}

int32_t paral_sym_decr(const LweSample *sample,
                       const TFheGateBootstrappingSecretKeySet *sk)
{
    Torus32 _1s16 = modSwitchToTorus32(1, 1 << PI);
    Torus32 mu = lwePhase(sample, sk->lwe_key);
    return (mu + (_1s16 >> 1)) >> (32 - PI);
}


// -----------------------------------------------------------------------------
//  LUT Bootstrapping: Threshold
//
void paral_bs_eq()
{
}

void paral_bs_gleq(LweSample *result,
                   const LweSample *sample,
                   const uint32_t thr,
                   const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweBootstrappingKeyFFT *bkFFT = bk->bkFFT;
    LweSample *tmp = new_LweSample(&bkFFT->accum_params->extracted_lweparams);
    const Torus32 MU = modSwitchToTorus32(1, 1 << PI);

    //~ tfhe_bootstrap_woKS_FFT(tmp, bkFFT, MU, sample);

            const TGswParams *bk_params = bkFFT->bk_params;
            const TLweParams *accum_params = bkFFT->accum_params;
            const LweParams *in_params = bkFFT->in_out_params;
            const int32_t N = accum_params->N;
            const int32_t Nx2 = 2 * N;
            const int32_t n = in_params->n;

            TorusPolynomial *testvect = new_TorusPolynomial(N);
            int32_t *bara = new int32_t[N];


            // Modulus switching
            int32_t barb = modSwitchFromTorus32(sample->b, Nx2);
            for (int32_t i = 0; i < n; i++) {
                bara[i] = modSwitchFromTorus32(sample->a[i], Nx2);
            }

            // the initial testvec = [mu,mu,mu,...,mu]
            for (int32_t i = 0; i < N/4; i++) testvect->coefsT[i] = 0*MU;
            for (int32_t i = N/4; i < 3*N/4; i++) testvect->coefsT[i] = 1*MU;
            for (int32_t i = 3*N/4; i < N; i++) testvect->coefsT[i] = 0*MU;

            // Bootstrapping rotation and extraction
            tfhe_blindRotateAndExtract_FFT(tmp, testvect, bkFFT->bkFFT, barb, bara, n, bk_params);
            //FUCKUP #02:   TFheGateBootstrappingCloudKeySet has member 'bkFFT',
            //              which is LweBootstrappingKeyFFT, which also has member 'bkFFT', which is TGswSampleFFT


            delete[] bara;
            delete_TorusPolynomial(testvect);

    // -------------------------------------------------------------------------

    // Key switching
    lweKeySwitch(result, bkFFT->ks, tmp);

    delete_LweSample(tmp);
}


// -----------------------------------------------------------------------------
//  Helper Functions
//
void paral_calc_qi(LweSample *qi,
                   const LweSample *w_i0,
                   const LweSample *w_i1,
                   const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    // aux variables
    LweSample *r1 = new_LweSample(io_lwe_params);
    //~ LweSample *r2 = new_LweSample(io_lwe_params);
    //~ LweSample *r3 = new_LweSample(io_lwe_params);
    LweSample *r23= new_LweSample(io_lwe_params);

    //TODO
    // r1   = w_i   <> +-3
    // r2   = w_i   == +-2
    // r3   = w_i-1 <> +-2
    // r23  = r2+r3 == +-2

    // q_i = r1 + r23
    lweNoiselessTrivial(qi, 0, io_lwe_params);
    lweAddTo(qi, r1, io_lwe_params);
    lweAddTo(qi, r23, io_lwe_params);
}


// -----------------------------------------------------------------------------
//  Parallel Addition
//
// !! there must be two more samples to the right in x and y !!
void paral_add(LweSample *z,
               const LweSample *x,
               const LweSample *y,
               const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    // alloc aux variable arrays for w_i, q_i
    LweSample *w = new_LweSample_array(3, io_lwe_params);
    LweSample *q = new_LweSample_array(2, io_lwe_params);

    // calc w_i = x_i + y_i   for i, i-1, i-2
    for (int i = 0; i < 3; i++)
    {
        lweNoiselessTrivial(w + i, 0, io_lwe_params);
        lweAddTo(w + i, x + i, io_lwe_params);
        lweAddTo(w + i, y + i, io_lwe_params);
    }

    // calc q_i and q_i-1
    paral_calc_qi(q + 0, w + 0, w + 1, bk);
    paral_calc_qi(q + 1, w + 1, w + 2, bk);

    // calculate result: z_i = w_i - 4q_i + q_i-1
    lweNoiselessTrivial(z, 0, io_lwe_params);
    lweAddTo(z, w, io_lwe_params);
    lweSubMulTo(z, 4, q, io_lwe_params);
    lweAddTo(z, q + 1, io_lwe_params);

    delete_LweSample_array(2, q);
    delete_LweSample_array(3, w);

    //TODO more like functional bootstrap: identity [-2 .. 2]
    //~ tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
}



// =============================================================================
//
//  MAIN
//

int32_t main(int32_t argc, char **argv)
{
#ifndef NDEBUG
    cout << "DEBUG MODE!" << endl;
#endif

    // roll dice
    srand(time(NULL));

    // generate TFHE params
    int32_t minimum_lambda = 100;
    //TODO generate appropriate params
    TFheGateBootstrappingParameterSet *tfhe_params = new_default_gate_bootstrapping_parameters(minimum_lambda);
    const LweParams *io_lwe_params = tfhe_params->in_out_params;
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *tfhe_keys = new_random_gate_bootstrapping_secret_keyset(tfhe_params);
    // alloc samples
    LweSample *x = new_LweSample(io_lwe_params);
    LweSample *y = new_LweSample(io_lwe_params);
    LweSample *r = new_LweSample(io_lwe_params);

    for (int32_t i = 0; i < 16; i++)
    {
        // encrypt
        paral_sym_encr_priv(x, i - 8, tfhe_keys);

        // bootstrap
        paral_bs_gleq(r, x, 3, &(tfhe_keys->cloud)); //FUCKUP #01: member 'cloud' is a struct, but not a pointer (unlike others)

        // decrypt
        int32_t x_plain = paral_sym_decr(x, tfhe_keys);
        int32_t r_plain = paral_sym_decr(r, tfhe_keys);

        printf("D[E(%+d)] = %+d .. bootstrapped to %+d\n", i-8, x_plain, r_plain);
    }

    //~ cout << "starting bootstrapping ..." << endl;

    //~ clock_t begin1 = clock();
    //~ bootsAND(r, x, y, &tfhe_keys->cloud);
    //~ clock_t end1 = clock();

    //~ cout << "finished bootstrapping, total time " << (end1 - begin1) << " [us]" << endl;

    // verify
    //~ int32_t x_plain = paral_sym_decr(x, tfhe_keys);
    //~ int32_t y_plain = paral_sym_decr(y, tfhe_keys);
    //~ int32_t r_plain = paral_sym_decr(r, tfhe_keys);

    //~ cout << "r = x + y ... " << r_plain << " = " << x_plain << " + " << y_plain << endl;

    // cleanup
    delete_LweSample(r);
    delete_LweSample(y);
    delete_LweSample(x);

    delete_gate_bootstrapping_secret_keyset(tfhe_keys);
    delete_gate_bootstrapping_parameters(tfhe_params);

    return 0;
}
