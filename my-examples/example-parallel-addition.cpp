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
    TFheGateBootstrappingParameterSet *tfhe_params = new_default_gate_bootstrapping_parameters(minimum_lambda);
    const LweParams *io_lwe_params = tfhe_params->in_out_params;
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *tfhe_keys = new_random_gate_bootstrapping_secret_keyset(tfhe_params);
    // alloc samples
    LweSample *x = new_LweSample(io_lwe_params);
    LweSample *y = new_LweSample(io_lwe_params);
    LweSample *r = new_LweSample(io_lwe_params);

    // encrypt
    bootsSymEncrypt(x, rand() % 2, tfhe_keys);
    bootsSymEncrypt(y, rand() % 2, tfhe_keys);

    cout << "starting bootstrapping ..." << endl;

    clock_t begin1 = clock();
    bootsAND(r, x, y, &tfhe_keys->cloud);
    clock_t end1 = clock();

    cout << "finished bootstrapping, total time " << (end1 - begin1) << " [us]" << endl;

    // verify
    bool x_plain = bootsSymDecrypt(x, tfhe_keys);
    bool y_plain = bootsSymDecrypt(y, tfhe_keys);
    bool r_plain = bootsSymDecrypt(r, tfhe_keys);

    cout << "r = x AND y ... " << (int)r_plain << " = " << x_plain << " AND " << y_plain << endl;

    // cleanup
    delete_LweSample(r);
    delete_LweSample(y);
    delete_LweSample(x);

    delete_gate_bootstrapping_secret_keyset(tfhe_keys);
    delete_gate_bootstrapping_parameters(tfhe_params);

    return 0;
}
