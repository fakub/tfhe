#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>

#include "tfhe.h"
#include "tfhe_garbage_collector.h"
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
    uint64_t seed = time(NULL);
    srand(seed);

    // w/o KS:
    // for  80-bit security @ pi = 4
    //      N =  512, n = 682, l = 20, mlogal = 23.745
    // for 128-bit security @ pi = 4
    //      N = 1024, n = 768, l = 20, mlogal = 24.331

    // with KS:
    //TODO generate appropriate TFHE params with KS

        static const int32_t N = 1024;
        static const int32_t k = 1;
        static const int32_t n = 630;
        static const int32_t bk_l = 3;
        static const int32_t bk_Bgbit = 7;
        static const int32_t ks_basebit = 2;
        static const int32_t ks_length = 8;
        static const double ks_stdev = pow(2.,-15); // standard deviation
        static const double bk_stdev = pow(2.,-25); // standard deviation
        static const double max_stdev = 0.012467;   // max standard deviation for a 1/4 msg space

        LweParams *params_in = new_LweParams(n, ks_stdev, max_stdev);
        TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
        TGswParams *params_bk = new_TGswParams(bk_l, bk_Bgbit, params_accum);

        TfheGarbageCollector::register_param(params_in);
        TfheGarbageCollector::register_param(params_accum);
        TfheGarbageCollector::register_param(params_bk);

        TFheGateBootstrappingParameterSet *tfhe_params = new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);


    const LweParams *io_lwe_params = tfhe_params->in_out_params;

    // seed TFHE
    tfhe_random_generator_setSeed((uint32_t*)(&seed), 2);
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
