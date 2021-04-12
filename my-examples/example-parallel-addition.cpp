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
    cerr << message << endl;
    abort();
}

void paral_sym_encr_priv(LweSample *ct,
                         const int32_t message,
                         const TFheGateBootstrappingSecretKeySet *key)
{
    Torus32 _1s16 = modSwitchToTorus32(1, 1 << PI);

    // scale to [-8/16, 7/16]
    //                         ___/ 8 \___     _/ mask 1111 \_     ___/ 8 \___
    Torus32 mu = (((message + (1 << (PI-1))) & ((1 << PI) - 1)) - (1 << (PI-1))) * _1s16;
    double alpha = key->params->in_out_params->alpha_min; //TODO: specify noise
    lweSymEncrypt(ct, mu, alpha, key->lwe_key);
}

void paral_sym_encr(LweSample *ct,
                    const int32_t message,
                    const TFheGateBootstrappingSecretKeySet *key)
{
    if (message < -2) || (2 < message)
        die_soon("Out of the alphabet A = [-2 .. 2].\n");

    paral_sym_encr_priv(ct, message, key);
}

void paral_calc_qi(LweSample *qi,
                   const LweSample *w_i0,
                   const LweSample *w_i1,
                   const TFheGateBootstrappingSecretKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    //TODO calc aux variables
    LweSample *r1 = new_LweSample(io_lwe_params);   // w_i   <> +-3
    LweSample *r2 = new_LweSample(io_lwe_params);   // w_i   == +-2
    LweSample *r3 = new_LweSample(io_lwe_params);   // w_i-1 <> +-2
    LweSample *r23= new_LweSample(io_lwe_params);   // r2+r3 == +-2   (r23)

    // q_i = r1 + r23
    lweNoiselessTrivial(qi, 0, io_lwe_params);
    lweAddTo(qi, r1, io_lwe_params);
    lweAddTo(qi, r23, io_lwe_params);
}

// there must be two more samples to the right in x and y !!
void paral_add(LweSample *z,
               const LweSample *x,
               const LweSample *y,
               const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    // calc w_i = x_i + y_i   for i, i-1, i-2
    LweSample *w_i0 = new_LweSample(io_lwe_params);
    lweAddTo(w_i0, x, io_lwe_params);
    lweAddTo(w_i0, y, io_lwe_params);

    LweSample *w_i1 = new_LweSample(io_lwe_params);
    lweAddTo(w_i1, x + 1, io_lwe_params);
    lweAddTo(w_i1, y + 1, io_lwe_params);

    LweSample *w_i2 = new_LweSample(io_lwe_params);
    lweAddTo(w_i2, x + 2, io_lwe_params);
    lweAddTo(w_i2, y + 2, io_lwe_params);

    // calc q_i0 and q_i1
    paral_calc_qi(q_i0, w_i0, w_i1, io_lwe_params)
    paral_calc_qi(q_i1, w_i1, w_i2, io_lwe_params)

    // calculate result: z_i = w_i - 4q_i + q_i-1
    lweNoiselessTrivial(z, 0, io_lwe_params);
    lweAddTo(z, w_i0, io_lwe_params);
    lweSubMulTo(z, 4, q_i0, io_lwe_params)
    lweAddTo(z, q_i1, io_lwe_params);

    //TODO more like functional bootstrap: identity [-2 .. 2]
    //~ tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
}

int32_t paral_sym_decr(const LweSample *sample,
                       const TFheGateBootstrappingSecretKeySet *key)
{
    Torus32 _1s16 = modSwitchToTorus32(1, 1 << PI);
    Torus32 mu = lwePhase(sample, key->lwe_key);
    return (mu + (_1s16 >> 1)) >> (32 - PI);
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

    // encrypt
    for (int32_t i = 0; i < 16; i++)
    {
        paral_sym_encr(x, i - 8, tfhe_keys);
        int32_t x_plain = paral_sym_decr(x, tfhe_keys);
        cout << "D[E(" << i-8 << ")] = " << x_plain << endl;
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
