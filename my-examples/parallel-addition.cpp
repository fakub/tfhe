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

#define SEQ_TEST
#define  PA_TEST
#define  BS_TEST

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
    if (SEQ_SCENARIO != SEQ_TFHE_PARAMS_INDEX) fprintf(stderr, "(w) TFHE parameters do not correspond with sequential addition scenario!\n");
    printf("\n\n================================================================================\n");
    printf("\n    <<<<    Sequential Addition Test    >>>>\n\n");
    printf("Scenario: %c\n", 'A' + SEQ_SCENARIO - 1);
    printf("Params:   %c\n\n", 'A' + SEQ_TFHE_PARAMS_INDEX - 1);fflush(stdout);

    // setup TFHE params
    TFheGateBootstrappingParameterSet *seq_tfhe_params = NULL;
    setup_TFHE_params(SEQ_TFHE_PARAMS_INDEX, &seq_tfhe_params);
    const LweParams *seq_io_lwe_params = seq_tfhe_params->in_out_params;
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *seq_tfhe_keys = new_random_gate_bootstrapping_secret_keyset(seq_tfhe_params);

    // alloc plaintexts & samples (n.b., opposite order, i.e., LSB-first)
    #define SEQ_BITLEN 20
    int32_t x_seq_plain_bin[SEQ_BITLEN] = {0,1,0,0,0,1,0,1,1,1,0,1,0,1,1,0,0,0,1,0,};   // LSB-first
    int32_t y_seq_plain_bin[SEQ_BITLEN] = {0,0,0,1,0,0,0,1,1,0,0,0,1,1,1,0,1,0,0,1,};

    // base 4 scenario C
#if (SEQ_SCENARIO == C_CARRY_2_BIT)
    #define SEQ_WLEN ((SEQ_BITLEN) / 2)
    int32_t x_seq_plain[SEQ_WLEN];
    int32_t y_seq_plain[SEQ_WLEN];
    for (int32_t i = 0; i < SEQ_WLEN; i++)
    {
        x_seq_plain[i] = x_seq_plain_bin[2*i] + 2*x_seq_plain_bin[2*i+1];
        y_seq_plain[i] = y_seq_plain_bin[2*i] + 2*y_seq_plain_bin[2*i+1];
    }
    int64_t (*seq_eval)(const int32_t *const, const uint32_t) = &quad_eval;

    // base 2 scenarios A, B
#else
    #define SEQ_WLEN SEQ_BITLEN
    int32_t x_seq_plain[SEQ_WLEN];
    int32_t y_seq_plain[SEQ_WLEN];
    for (int32_t i = 0; i < SEQ_WLEN; i++)
    {
        x_seq_plain[i] = x_seq_plain_bin[i];
        y_seq_plain[i] = y_seq_plain_bin[i];
    }
    int64_t (*seq_eval)(const int32_t *const, const uint32_t) = &bin_eval;
#endif
    int32_t z_seq_plain[SEQ_WLEN + 1];

    int64_t exp_seq_sum = seq_eval(&x_seq_plain[0], SEQ_WLEN) + seq_eval(&y_seq_plain[0], SEQ_WLEN);

    LweSample *x_seq = new_LweSample_array(SEQ_WLEN,     seq_io_lwe_params);   // for sequential addition
    LweSample *y_seq = new_LweSample_array(SEQ_WLEN,     seq_io_lwe_params);
    LweSample *z_seq = new_LweSample_array(SEQ_WLEN + 1, seq_io_lwe_params);

    // encrypt
    for (int32_t i = 0; i < SEQ_WLEN; i++)
    {
#if (SEQ_SCENARIO == C_CARRY_2_BIT)
        seq_quad_sym_encr(x_seq + i, x_seq_plain[i], seq_tfhe_keys);
        seq_quad_sym_encr(y_seq + i, y_seq_plain[i], seq_tfhe_keys);
#else
        bin_sym_encr(x_seq + i, x_seq_plain[i], seq_tfhe_keys);
        bin_sym_encr(y_seq + i, y_seq_plain[i], seq_tfhe_keys);
#endif
    }

    // print inputs
    for (int32_t i = 0; i < SEQ_WLEN+2; i++) printf("-----");
    printf("\n");

    // x
    printf(" X  | +0 ");
    for (int32_t i = SEQ_WLEN - 1; i >= 0; i--)
        printf("| %+d ", x_seq_plain[i]);
    printf("| %+9ld\n", seq_eval(&x_seq_plain[0], SEQ_WLEN));

    // y
    printf(" Y  | +0 ");
    for (int32_t i = SEQ_WLEN - 1; i >= 0; i--)
        printf("| %+d ", y_seq_plain[i]);
    printf("| %+9ld\n", seq_eval(&y_seq_plain[0], SEQ_WLEN));
    // ------------------------------------------------------
    for (int32_t i = 0; i < SEQ_WLEN+2; i++) printf("-----");
    printf("   ");

    // sequential addition
    sequential_add(z_seq, x_seq, y_seq, SEQ_WLEN,
#ifdef DBG_OUT
                   seq_tfhe_keys,
#endif
                   &(seq_tfhe_keys->cloud));

    // decrypt
    for (int32_t i = 0; i <= SEQ_WLEN; i++)
#if (SEQ_SCENARIO == C_CARRY_2_BIT)
        z_seq_plain[i] = sym_decr(z_seq + i, seq_tfhe_keys);
#else
        z_seq_plain[i] = bin_sym_decr(z_seq + i, seq_tfhe_keys);
#endif

    // print results
    // z
    printf(" Z  ");
    for (int32_t i = SEQ_WLEN; i >= 0; i--)
    {
        printf("| %+d ", z_seq_plain[i]);
    }
    printf("| %+9ld   %s (exp. %+9ld)\n",
                            seq_eval(&z_seq_plain[0], SEQ_WLEN + 1),
                            exp_seq_sum == seq_eval(&z_seq_plain[0], SEQ_WLEN + 1) ? "\033[1;32mPASS\033[0m" : "\033[1;31mFAIL\033[0m",
                            exp_seq_sum);
    for (int32_t i = 0; i < SEQ_WLEN+2; i++) printf("-----");
    printf("\n");

    // cleanup
    delete_LweSample_array(SEQ_WLEN,     x_seq);
    delete_LweSample_array(SEQ_WLEN,     y_seq);
    delete_LweSample_array(SEQ_WLEN + 1, z_seq);
#endif


#ifdef PA_TEST
    // -------------------------------------------------------------------------
    //
    //  Parallel Addition Test
    //
    if (PA_SCENARIO != PA_TFHE_PARAMS_INDEX) fprintf(stderr, "(w) TFHE parameters do not correspond with parallel addition scenario!\n");
    printf("\n\n================================================================================\n");
    printf("\n    <<<<    Parallel Addition Test    >>>>\n\n");
    printf("Scenario: %c\n", 'A' + PA_SCENARIO - 1);
    printf("Params:   %c\n\n", 'A' + PA_TFHE_PARAMS_INDEX - 1);fflush(stdout);

    // setup TFHE params
    TFheGateBootstrappingParameterSet *pa_tfhe_params = NULL;
    setup_TFHE_params(PA_TFHE_PARAMS_INDEX, &pa_tfhe_params);
    const LweParams *pa_io_lwe_params = pa_tfhe_params->in_out_params;
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *pa_tfhe_keys = new_random_gate_bootstrapping_secret_keyset(pa_tfhe_params);

    // alloc plaintexts & samples (n.b., opposite order, i.e., LSB-first)
    #define PA_WLEN 10
    int32_t x_pa_plain[PA_WLEN] = {-2,+1,+2,+2,+1,-2,-1,+2,+0,+1,};   // LSB-first
    int32_t y_pa_plain[PA_WLEN] = {+0,-2,+1,+2,+1,+0,-2,+2,+1,+2,};
    int32_t z_pa_plain[PA_WLEN + 1];

    int64_t exp_pa_sum = quad_eval(&x_pa_plain[0], PA_WLEN) + quad_eval(&y_pa_plain[0], PA_WLEN);

    LweSample *x_pa = new_LweSample_array(PA_WLEN,     pa_io_lwe_params);   // for parallel addition
    LweSample *y_pa = new_LweSample_array(PA_WLEN,     pa_io_lwe_params);
    LweSample *z_pa = new_LweSample_array(PA_WLEN + 1, pa_io_lwe_params);

    // encrypt
    for (int32_t i = 0; i < PA_WLEN; i++)
    {
        paral_sym_encr(x_pa + i, x_pa_plain[i], pa_tfhe_keys);
        paral_sym_encr(y_pa + i, y_pa_plain[i], pa_tfhe_keys);
    }

    // print inputs
    for (int32_t i = 0; i < PA_WLEN+2; i++) printf("-----");
    printf("\n");

    // x
    printf(" X  | +0 ");
    for (int32_t i = PA_WLEN - 1; i >= 0; i--)
        printf("| %+d ", x_pa_plain[i]);
    printf("| %+9ld\n", quad_eval(&x_pa_plain[0], PA_WLEN));

    // y
    printf(" Y  | +0 ");
    for (int32_t i = PA_WLEN - 1; i >= 0; i--)
        printf("| %+d ", y_pa_plain[i]);
    printf("| %+9ld\n", quad_eval(&y_pa_plain[0], PA_WLEN));
    for (int32_t i = 0; i < PA_WLEN+2; i++) printf("-----");
    printf("   ");

    // parallel addition
    parallel_add(z_pa, x_pa, y_pa, PA_WLEN, &(pa_tfhe_keys->cloud));

    // decrypt
    for (int32_t i = 0; i <= PA_WLEN; i++)
        z_pa_plain[i] = sym_decr(z_pa + i, pa_tfhe_keys);

    // print results
    // z
    printf(" Z  ");
    for (int32_t i = PA_WLEN; i >= 0; i--)
    {
        printf("| %+d ", z_pa_plain[i]);
    }
    printf("| %+9ld   %s (exp. %+9ld)\n",
                            quad_eval(&z_pa_plain[0], PA_WLEN + 1),
                            exp_pa_sum == quad_eval(&z_pa_plain[0], PA_WLEN + 1) ? "\033[1;32mPASS\033[0m" : "\033[1;31mFAIL\033[0m",
                            exp_pa_sum);
    for (int32_t i = 0; i < PA_WLEN+2; i++) printf("-----");
    printf("\n");

    // cleanup
    delete_LweSample_array(PA_WLEN,     x_pa);
    delete_LweSample_array(PA_WLEN,     y_pa);
    delete_LweSample_array(PA_WLEN + 1, z_pa);
#endif


#ifdef BS_TEST
    // -------------------------------------------------------------------------
    //
    //  Bootstrapping Test
    //
    printf("\n\n================================================================================\n");
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
        sym_encr_priv(a, i - (1 << (PI-1)), bs_tfhe_keys);

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
        int32_t a_plain     = sym_decr(a,  bs_tfhe_keys);
        int32_t id_plain    = sym_decr(id, bs_tfhe_keys);
        int32_t gl_plain    = sym_decr(gl, bs_tfhe_keys);
        int32_t eq_plain    = sym_decr(eq, bs_tfhe_keys);

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
