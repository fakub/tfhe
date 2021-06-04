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

//~ #define SEQ_TEST
#define  PA_TEST
//~ #define  BS_TEST

using namespace std;


int32_t main(int32_t argc, char **argv)
{
#ifndef NDEBUG
    cout << "DEBUG MODE!" << endl;
#endif

    // roll dice
    uint64_t seed = time(NULL);
    srand(seed);
    // seed TFHE
    tfhe_random_generator_setSeed((uint32_t*)(&seed), 2);


#ifdef SEQ_TEST
    // -------------------------------------------------------------------------
    //
    //  Sequential Addition with Carry Test
    //
    if (SEQ_SCENARIO != SEQ_TFHE_PARAMS_INDEX) fprintf(stderr, "(w) TFHE parameters do not correspond with parallel addition scenario!\n");
    printf("\n    <<<<    Sequential Addition Test    >>>>\n\n");
    printf("Scenario: %c\n", 'A' + SEQ_SCENARIO - 1);
    printf("Params:   %c\n\n", 'A' + SEQ_TFHE_PARAMS_INDEX - 1);fflush(stdout);

    // setup TFHE params
    TFheGateBootstrappingParameterSet *seq_tfhe_params = NULL;
    setup_TFHE_params(SEQ_TFHE_PARAMS_INDEX, &seq_tfhe_params);
    const LweParams *seq_io_lwe_params = seq_tfhe_params->in_out_params;
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *seq_tfhe_keys = new_random_gate_bootstrapping_secret_keyset(seq_tfhe_params);

    // alloc plaintexts & samples (n.b., opposite order, i.e., big endian (?))
    //~ #define PA_WLEN 10
    //~ int32_t x_plain[PA_WLEN] = {-2,+1,+2,+2,+1,-2,-1,+2,+0,+1,};   // LSB-first
    //~ int32_t y_plain[PA_WLEN] = {+0,-2,+1,+2,+1,+0,-2,+2,+1,+2,};
    //~ int32_t z_plain[PA_WLEN + 1];

    //~ int64_t exp_sum = paral_eval(&x_plain[0], PA_WLEN) + paral_eval(&y_plain[0], PA_WLEN);

    //~ LweSample *x = new_LweSample_array(PA_WLEN,     seq_io_lwe_params);   // for parallel addition
    //~ LweSample *y = new_LweSample_array(PA_WLEN,     seq_io_lwe_params);
    //~ LweSample *z = new_LweSample_array(PA_WLEN + 1, seq_io_lwe_params);

    //~ // encrypt
    //~ for (int32_t i = 0; i < PA_WLEN; i++)
    //~ {
        //~ paral_sym_encr(x + i, x_plain[i], seq_tfhe_keys);
        //~ paral_sym_encr(y + i, y_plain[i], seq_tfhe_keys);
    //~ }

    //~ // print inputs
    //~ printf("------------------------------------------------------------\n");

    //~ // x
    //~ printf(" X  | +0 ");
    //~ for (int32_t i = PA_WLEN - 1; i >= 0; i--)
        //~ printf("| %+d ", x_plain[i]);
    //~ printf("| %+9ld\n", paral_eval(&x_plain[0], PA_WLEN));

    //~ // y
    //~ printf(" Y  | +0 ");
    //~ for (int32_t i = PA_WLEN - 1; i >= 0; i--)
        //~ printf("| %+d ", y_plain[i]);
    //~ printf("| %+9ld\n------------------------------------------------------------   ", paral_eval(&y_plain[0], PA_WLEN));

    //~ // parallel addition
    //~ parallel_add(z, x, y, PA_WLEN, &(tfhe_keys->cloud));

    //~ // decrypt
    //~ for (int32_t i = 0; i <= PA_WLEN; i++)
        //~ z_plain[i] = paral_sym_decr(z + i, tfhe_keys);

    //~ // print results
    //~ // z
    //~ printf(" Z  ");
    //~ for (int32_t i = PA_WLEN; i >= 0; i--)
    //~ {
        //~ printf("| %+d ", z_plain[i]);
    //~ }
    //~ printf("| %+9ld   %s (exp. %+9ld)\n",
              //~ paral_eval(&z_plain[0], PA_WLEN + 1),
                      //~ exp_sum == paral_eval(&z_plain[0], PA_WLEN + 1) ? "\033[1;32mPASS\033[0m" : "\033[1;31mFAIL\033[0m",
                               //~ exp_sum);
    //~ printf("------------------------------------------------------------\n");

    //~ // cleanup
    //~ delete_LweSample_array(PA_WLEN,     x);
    //~ delete_LweSample_array(PA_WLEN,     y);
    //~ delete_LweSample_array(PA_WLEN + 1, z);
#endif


#ifdef PA_TEST
    // -------------------------------------------------------------------------
    //
    //  Parallel Addition Test
    //
    if (PA_SCENARIO != PA_TFHE_PARAMS_INDEX) fprintf(stderr, "(w) TFHE parameters do not correspond with parallel addition scenario!\n");
    printf("\n    <<<<    Parallel Addition Test    >>>>\n\n");
    printf("Scenario: %c\n", 'A' + PA_SCENARIO - 1);
    printf("Params:   %c\n\n", 'A' + PA_TFHE_PARAMS_INDEX - 1);fflush(stdout);

    // setup TFHE params
    TFheGateBootstrappingParameterSet *pa_tfhe_params = NULL;
    setup_TFHE_params(PA_TFHE_PARAMS_INDEX, &pa_tfhe_params);
    const LweParams *pa_io_lwe_params = pa_tfhe_params->in_out_params;
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *pa_tfhe_keys = new_random_gate_bootstrapping_secret_keyset(pa_tfhe_params);

    // alloc plaintexts & samples (n.b., opposite order, i.e., big endian (?))
    #define PA_WLEN 10
    int32_t x_plain[PA_WLEN] = {-2,+1,+2,+2,+1,-2,-1,+2,+0,+1,};   // LSB-first
    int32_t y_plain[PA_WLEN] = {+0,-2,+1,+2,+1,+0,-2,+2,+1,+2,};
    int32_t z_plain[PA_WLEN + 1];

    int64_t exp_sum = paral_eval(&x_plain[0], PA_WLEN) + paral_eval(&y_plain[0], PA_WLEN);

    LweSample *x = new_LweSample_array(PA_WLEN,     pa_io_lwe_params);   // for parallel addition
    LweSample *y = new_LweSample_array(PA_WLEN,     pa_io_lwe_params);
    LweSample *z = new_LweSample_array(PA_WLEN + 1, pa_io_lwe_params);

    // encrypt
    for (int32_t i = 0; i < PA_WLEN; i++)
    {
        paral_sym_encr(x + i, x_plain[i], pa_tfhe_keys);
        paral_sym_encr(y + i, y_plain[i], pa_tfhe_keys);
    }

    // print inputs
    printf("------------------------------------------------------------\n");

    // x
    printf(" X  | +0 ");
    for (int32_t i = PA_WLEN - 1; i >= 0; i--)
        printf("| %+d ", x_plain[i]);
    printf("| %+9ld\n", paral_eval(&x_plain[0], PA_WLEN));

    // y
    printf(" Y  | +0 ");
    for (int32_t i = PA_WLEN - 1; i >= 0; i--)
        printf("| %+d ", y_plain[i]);
    printf("| %+9ld\n------------------------------------------------------------   ", paral_eval(&y_plain[0], PA_WLEN));

    // parallel addition
    parallel_add(z, x, y, PA_WLEN, &(pa_tfhe_keys->cloud));

    // decrypt
    for (int32_t i = 0; i <= PA_WLEN; i++)
        z_plain[i] = paral_sym_decr(z + i, pa_tfhe_keys);

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


#ifdef BS_TEST
    // -------------------------------------------------------------------------
    //
    //  Bootstrapping Test
    //
    printf("\n    <<<<    Bootstrapping Test    >>>>\n\n");fflush(stdout);

    // setup TFHE params
    TFheGateBootstrappingParameterSet *bs_tfhe_params = NULL;
    setup_TFHE_params(BS_TFHE_PARAMS_INDEX, &bs_tfhe_params);
    const LweParams *bs_io_lwe_params = bs_tfhe_params->in_out_params;
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *bs_tfhe_keys = new_random_gate_bootstrapping_secret_keyset(bs_tfhe_params);

    // alloc samples
    LweSample *a  = new_LweSample(bs_io_lwe_params);   // for BS testing
    LweSample *id = new_LweSample(bs_io_lwe_params);
    LweSample *gl = new_LweSample(bs_io_lwe_params);
    LweSample *eq = new_LweSample(bs_io_lwe_params);

    // print table heading
    printf("--------------------------------------------------------------------------------\n");
    printf(" Encr -> Decr  | Id. | <=> 3 | == 2 | timing (Id.)\n");
    printf("--------------------------------------------------------------------------------\n");

    for (int32_t i = 0; i < (1 << PI); i++)
    {
        // encrypt
        paral_sym_encr_priv(a, i - (1 << (PI-1)), bs_tfhe_keys);

        // bootstrap
        clock_t begin_id = clock();
        bs_id(id, a, &(bs_tfhe_keys->cloud)); //FUCKUP #01: member 'cloud' is a struct, but not a pointer (unlike others)
        clock_t end_id = clock();

        // clock_t begin_gl = clock();
        bs_gleq(gl, a, 3, &(bs_tfhe_keys->cloud));
        // clock_t end_gl = clock();

        // clock_t begin_eq = clock();
        bs_eq(eq, a, 2, &(bs_tfhe_keys->cloud));
        // clock_t end_eq = clock();

        // decrypt
        int32_t a_plain     = paral_sym_decr(a,  bs_tfhe_keys);
        int32_t id_plain    = paral_sym_decr(id, bs_tfhe_keys);
        int32_t gl_plain    = paral_sym_decr(gl, bs_tfhe_keys);
        int32_t eq_plain    = paral_sym_decr(eq, bs_tfhe_keys);

        printf(" D[E(%+d)] = %+d |  %+d |   %+d  |  %+d  | %lu ms\n", i - (1 << (PI-1)), a_plain,
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


    // -------------------------------------------------------------------------
    //
    //  Cleanup
    //

#ifdef SEQ_TEST
    delete_gate_bootstrapping_secret_keyset(seq_tfhe_keys);
    delete_gate_bootstrapping_parameters(   seq_tfhe_params);
#endif
#ifdef  PA_TEST
    delete_gate_bootstrapping_secret_keyset( pa_tfhe_keys);
    delete_gate_bootstrapping_parameters(    pa_tfhe_params);
#endif
#ifdef  BS_TEST
    delete_gate_bootstrapping_secret_keyset( bs_tfhe_keys);
    delete_gate_bootstrapping_parameters(    bs_tfhe_params);
#endif

    return 0;
}
