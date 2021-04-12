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

using namespace std;

void paral_sym_encr(LweSample *result,
                    const int32_t message,
                    const TFheGateBootstrappingSecretKeySet *key)
{
    Torus32 _1s16 = modSwitchToTorus32(1, 16);
    Torus32 mu = (((message + 8) & 0x0f) - 8) * _1s16;  // scales to [-8/16, 7/16], consider only [-2,2] interval?
    double alpha = key->params->in_out_params->alpha_min; //TODO: specify noise
    lweSymEncrypt(result, mu, alpha, key->lwe_key);
}

void paral_add(LweSample *result,
               const LweSample *ca,
               const LweSample *cb,
               const TFheGateBootstrappingCloudKeySet *bk)
{
    static const Torus32 MU = modSwitchToTorus32(1, 8);
    const LweParams *in_out_params = bk->params->in_out_params;

    LweSample *temp_result = new_LweSample(in_out_params);

    //compute: (0,-1/8) + ca + cb
    static const Torus32 AndConst = modSwitchToTorus32(-1, 8);
    lweNoiselessTrivial(temp_result, AndConst, in_out_params);
    lweAddTo(temp_result, ca, in_out_params);
    lweAddTo(temp_result, cb, in_out_params);

    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);

    delete_LweSample(temp_result);
}

int32_t paral_sym_decr(const LweSample *sample,
                       const TFheGateBootstrappingSecretKeySet *key)
{
    Torus32 _1s16 = modSwitchToTorus32(1, 16);
    Torus32 mu = lwePhase(sample, key->lwe_key);
    return (mu + (_1s16 >> 1)) / _1s16;
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
