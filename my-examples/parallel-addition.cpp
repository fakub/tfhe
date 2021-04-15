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

#include <parallel-addition-impl.h>

//~ #define BS_TEST
#define PA_TEST

using namespace std;


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
    //TODO roll TFHE dice
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *tfhe_keys = new_random_gate_bootstrapping_secret_keyset(tfhe_params);


#ifdef BS_TEST
    // -------------------------------------------------------------------------
    //
    //  Bootstrapping Test
    //

    // alloc samples
    LweSample *a  = new_LweSample(io_lwe_params);   // for BS testing
    LweSample *id = new_LweSample(io_lwe_params);
    LweSample *gl = new_LweSample(io_lwe_params);
    LweSample *eq = new_LweSample(io_lwe_params);

    // print table heading
    printf("--------------------------------------------------------------------------------\n");
    printf(" Encr -> Decr  | Id. | <=> 3 | == 2 | timing (Id.)\n");
    printf("--------------------------------------------------------------------------------\n");

    for (int32_t i = 0; i < (1 << PI); i++)
    {
        // encrypt
        paral_sym_encr_priv(a, i - (1 << (PI-1)), tfhe_keys);

        // bootstrap
        clock_t begin_id = clock();
        bs_id(id, a, &(tfhe_keys->cloud)); //FUCKUP #01: member 'cloud' is a struct, but not a pointer (unlike others)
        clock_t end_id = clock();

        // clock_t begin_gl = clock();
        bs_gleq(gl, a, 3, &(tfhe_keys->cloud));
        // clock_t end_gl = clock();

        // clock_t begin_eq = clock();
        bs_eq(eq, a, 2, &(tfhe_keys->cloud));
        // clock_t end_eq = clock();

        // decrypt
        int32_t a_plain     = paral_sym_decr(a,  tfhe_keys);
        int32_t id_plain    = paral_sym_decr(id, tfhe_keys);
        int32_t gl_plain    = paral_sym_decr(gl, tfhe_keys);
        int32_t eq_plain    = paral_sym_decr(eq, tfhe_keys);

        printf(" D[E(%+d)] = %+d |  %+d |   %+d  |  %+d  | %lu ms\n", i-8, a_plain,
                                    id_plain, gl_plain, eq_plain,
                                                            (end_id - begin_id) / 1000);
    }
    printf("\n");

    // cleanup
    delete_LweSample(eq);
    delete_LweSample(gl);
    delete_LweSample(id);
    delete_LweSample(a);
#endif


#ifdef PA_TEST
    // -------------------------------------------------------------------------
    //
    //  Parallel Addition Test
    //

    // alloc plaintexts & samples (n.b., opposite order, i.e., big endian (?))
    #define PA_WLEN 10
    int32_t x_plain[PA_WLEN] = {-2,+1,+2,+2,+1,-2,-1,+2,+0,+1,};   // LSB-first
    int32_t y_plain[PA_WLEN] = {+0,-2,+1,+2,+1,+0,-2,+2,+1,+2,};
    int32_t z_plain[PA_WLEN + 1];

    int64_t exp_sum = paral_eval(&x_plain[0], PA_WLEN) + paral_eval(&y_plain[0], PA_WLEN);

    LweSample *x = new_LweSample_array(PA_WLEN,     io_lwe_params);   // for parallel addition
    LweSample *y = new_LweSample_array(PA_WLEN,     io_lwe_params);
    LweSample *z = new_LweSample_array(PA_WLEN + 1, io_lwe_params);

    // encrypt
    for (int32_t i = 0; i < PA_WLEN; i++)
    {
        paral_sym_encr(x + i, x_plain[i], tfhe_keys);
        paral_sym_encr(y + i, y_plain[i], tfhe_keys);
    }

    // print inputs
    printf("\n------------------------------------------------------------\n");

    // x
    printf(" X  | +0 ");
    for (int32_t i = PA_WLEN - 1; i >= 0; i--)
        printf("| %+d ", x_plain[i]);
    printf("| %+9ld\n", paral_eval(&x_plain[0], PA_WLEN));

    // y
    printf(" Y  | +0 ");
    for (int32_t i = PA_WLEN - 1; i >= 0; i--)
        printf("| %+d ", y_plain[i]);
    printf("| %+9ld\n------------------------------------------------------------\n", paral_eval(&y_plain[0], PA_WLEN));

    // parallel addition
    parallel_add(z, x, y, PA_WLEN, &(tfhe_keys->cloud));

    // decrypt
    for (int32_t i = 0; i <= PA_WLEN; i++)
        z_plain[i] = paral_sym_decr(z + i, tfhe_keys);

    // print results
    // z
    printf(" Z  ");
    for (int32_t i = PA_WLEN; i >= 0; i--)
    {
        printf("| %+d ", z_plain[i]);
    }
    printf("| %+9ld   %s (exp. %+9ld)\n",
              paral_eval(&z_plain[0], PA_WLEN + 1),
                      exp_sum == paral_eval(&z_plain[0], PA_WLEN + 1) ? "\033[1;32mPASS\033[0m" : "\033[1;31mFAIL\033[0m",
                               exp_sum);
    printf("------------------------------------------------------------\n");

    // cleanup
    delete_LweSample_array(PA_WLEN,     x);
    delete_LweSample_array(PA_WLEN,     y);
    delete_LweSample_array(PA_WLEN + 1, z);
#endif


    // -------------------------------------------------------------------------
    //
    //  Cleanup
    //

    delete_gate_bootstrapping_secret_keyset(tfhe_keys);
    delete_gate_bootstrapping_parameters(tfhe_params);

    return 0;
}
