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

//~ #define  SEQ_TEST
#define  PA_TEST_BIN
#define  PA_TEST_QUAD
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
        z_seq_plain[i] = sym_decr(z_seq + i, PI_S, seq_tfhe_keys);
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


#ifdef PA_TEST_BIN   //TODO
    // -------------------------------------------------------------------------
    //
    //  Parallel Addition Test
    //
    if (PA_SCENARIO_BIN != PAB_TFHE_PARAMS_INDEX) fprintf(stderr, "(w) TFHE parameters do not correspond with binary parallel addition scenario!\n");
    printf("\n\n================================================================================\n");
    printf("\n    <<<<    Parallel Addition Test -- Binary    >>>>\n\n");
    printf("Scenario: %c\n", 'A' + PA_SCENARIO_BIN - 1);
    printf("Params:   %c\n\n", 'A' + PAB_TFHE_PARAMS_INDEX - 1);fflush(stdout);

    // setup TFHE params
    TFheGateBootstrappingParameterSet *pa_bin_tfhe_params = NULL;
    setup_TFHE_params(PAB_TFHE_PARAMS_INDEX, &pa_bin_tfhe_params);
    const LweParams *pa_bin_io_lwe_params = pa_bin_tfhe_params->in_out_params;
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *pa_bin_tfhe_keys = new_random_gate_bootstrapping_secret_keyset(pa_bin_tfhe_params);

    // alloc plaintexts & samples (n.b., opposite order, i.e., LSB-first)
    #define PA_BIN_WLEN 20
    int32_t x_pa_bin_plain[PA_BIN_WLEN] = {0, 1,-1, 0, 1,-1, 0,-1,-1,-1, 0, 1, 0, 1, 1, 0,-1, 0, 1,-1, };   // LSB-first
    int32_t y_pa_bin_plain[PA_BIN_WLEN] = {0, 1, 0, 1, 1, 0, 1,-1,-1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, };
    int32_t z_pa_bin_plain[PA_BIN_WLEN + 1];

    int64_t exp_pa_bin_sum = bin_eval(&x_pa_bin_plain[0], PA_BIN_WLEN) + bin_eval(&y_pa_bin_plain[0], PA_BIN_WLEN);

    LweSample *x_pa_bin = new_LweSample_array(PA_BIN_WLEN,     pa_bin_io_lwe_params);   // for parallel addition
    LweSample *y_pa_bin = new_LweSample_array(PA_BIN_WLEN,     pa_bin_io_lwe_params);
    LweSample *z_pa_bin = new_LweSample_array(PA_BIN_WLEN + 1, pa_bin_io_lwe_params);

    // encrypt
    for (int32_t i = 0; i < PA_BIN_WLEN; i++)
    {
        paral_bin_sym_encr(x_pa_bin + i, x_pa_bin_plain[i], pa_bin_tfhe_keys);
        paral_bin_sym_encr(y_pa_bin + i, y_pa_bin_plain[i], pa_bin_tfhe_keys);
    }

    // print inputs
    for (int32_t i = 0; i < PA_BIN_WLEN+2; i++) printf("-----");
    printf("\n");

    // x
    printf(" X  | +0 ");
    for (int32_t i = PA_BIN_WLEN - 1; i >= 0; i--)
        printf("| %+d ", x_pa_bin_plain[i]);
    printf("| %+9ld\n", bin_eval(&x_pa_bin_plain[0], PA_BIN_WLEN));

    // y
    printf(" Y  | +0 ");
    for (int32_t i = PA_BIN_WLEN - 1; i >= 0; i--)
        printf("| %+d ", y_pa_bin_plain[i]);
    printf("| %+9ld\n", bin_eval(&y_pa_bin_plain[0], PA_BIN_WLEN));
    for (int32_t i = 0; i < PA_BIN_WLEN+2; i++) printf("-----");
    printf("   ");

    // parallel addition
    parallel_add_bin(z_pa_bin, x_pa_bin, y_pa_bin,
                     PA_BIN_WLEN,
#ifdef DBG_OUT
                     pa_bin_tfhe_keys,
#endif
                     &(pa_bin_tfhe_keys->cloud));

    // decrypt
    for (int32_t i = 0; i <= PA_BIN_WLEN; i++)
        z_pa_bin_plain[i] = sym_decr(z_pa_bin + i, PI_B, pa_bin_tfhe_keys);

    // print results
    // z
    printf(" Z  ");
    for (int32_t i = PA_BIN_WLEN; i >= 0; i--)
    {
        printf("| %+d ", z_pa_bin_plain[i]);
    }
    printf("| %+9ld   %s (exp. %+9ld)\n",
                            bin_eval(&z_pa_bin_plain[0], PA_BIN_WLEN + 1),
                            exp_pa_bin_sum == bin_eval(&z_pa_bin_plain[0], PA_BIN_WLEN + 1) ? "\033[1;32mPASS\033[0m" : "\033[1;31mFAIL\033[0m",
                            exp_pa_bin_sum);
    for (int32_t i = 0; i < PA_BIN_WLEN+2; i++) printf("-----");
    printf("\n");

    // cleanup
    delete_LweSample_array(PA_BIN_WLEN,     x_pa_bin);
    delete_LweSample_array(PA_BIN_WLEN,     y_pa_bin);
    delete_LweSample_array(PA_BIN_WLEN + 1, z_pa_bin);
#endif


#ifdef PA_TEST_QUAD
    // -------------------------------------------------------------------------
    //
    //  Parallel Addition Test
    //
    if (PA_SCENARIO_QUAD != PAQ_TFHE_PARAMS_INDEX) fprintf(stderr, "(w) TFHE parameters do not correspond with quad parallel addition scenario!\n");
    printf("\n\n================================================================================\n");
    printf("\n    <<<<    Parallel Addition Test -- Quad    >>>>\n\n");
    printf("Scenario: %c\n", 'A' + PA_SCENARIO_QUAD - 1);
    printf("Params:   %c\n\n", 'A' + PAQ_TFHE_PARAMS_INDEX - 1);fflush(stdout);

    // setup TFHE params
    TFheGateBootstrappingParameterSet *pa_quad_tfhe_params = NULL;
    setup_TFHE_params(PAQ_TFHE_PARAMS_INDEX, &pa_quad_tfhe_params);
    const LweParams *pa_quad_io_lwe_params = pa_quad_tfhe_params->in_out_params;
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *pa_quad_tfhe_keys = new_random_gate_bootstrapping_secret_keyset(pa_quad_tfhe_params);

    // alloc plaintexts & samples (n.b., opposite order, i.e., LSB-first)
    #define PA_QUAD_WLEN 10
    int32_t x_pa_plain[PA_QUAD_WLEN] = {-2,+1,+2,+2,+1,-2,-1,+2,+0,+1,};   // LSB-first
    int32_t y_pa_plain[PA_QUAD_WLEN] = {+0,-2,+1,+2,+1,+0,-2,+2,+1,+2,};
    int32_t z_pa_plain[PA_QUAD_WLEN + 1];

    int64_t exp_pa_quad_sum = quad_eval(&x_pa_plain[0], PA_QUAD_WLEN) + quad_eval(&y_pa_plain[0], PA_QUAD_WLEN);

    LweSample *x_pa_quad = new_LweSample_array(PA_QUAD_WLEN,     pa_quad_io_lwe_params);   // for parallel addition
    LweSample *y_pa_quad = new_LweSample_array(PA_QUAD_WLEN,     pa_quad_io_lwe_params);
    LweSample *z_pa_quad = new_LweSample_array(PA_QUAD_WLEN + 1, pa_quad_io_lwe_params);

    // encrypt
    for (int32_t i = 0; i < PA_QUAD_WLEN; i++)
    {
        paral_quad_sym_encr(x_pa_quad + i, x_pa_plain[i], pa_quad_tfhe_keys);
        paral_quad_sym_encr(y_pa_quad + i, y_pa_plain[i], pa_quad_tfhe_keys);
    }

    // print inputs
    for (int32_t i = 0; i < PA_QUAD_WLEN+2; i++) printf("-----");
    printf("\n");

    // x
    printf(" X  | +0 ");
    for (int32_t i = PA_QUAD_WLEN - 1; i >= 0; i--)
        printf("| %+d ", x_pa_plain[i]);
    printf("| %+9ld\n", quad_eval(&x_pa_plain[0], PA_QUAD_WLEN));

    // y
    printf(" Y  | +0 ");
    for (int32_t i = PA_QUAD_WLEN - 1; i >= 0; i--)
        printf("| %+d ", y_pa_plain[i]);
    printf("| %+9ld\n", quad_eval(&y_pa_plain[0], PA_QUAD_WLEN));
    for (int32_t i = 0; i < PA_QUAD_WLEN+2; i++) printf("-----");
    printf("   ");

    // parallel addition
    parallel_add_quad(z_pa_quad, x_pa_quad, y_pa_quad,
                      PA_QUAD_WLEN,
#ifdef DBG_OUT
                      pa_quad_tfhe_keys,
#endif
                      &(pa_quad_tfhe_keys->cloud));

    // decrypt
    for (int32_t i = 0; i <= PA_QUAD_WLEN; i++)
        z_pa_plain[i] = sym_decr(z_pa_quad + i, PI_Q, pa_quad_tfhe_keys);

    // print results
    // z
    printf(" Z  ");
    for (int32_t i = PA_QUAD_WLEN; i >= 0; i--)
    {
        printf("| %+d ", z_pa_plain[i]);
    }
    printf("| %+9ld   %s (exp. %+9ld)\n",
                            quad_eval(&z_pa_plain[0], PA_QUAD_WLEN + 1),
                            exp_pa_quad_sum == quad_eval(&z_pa_plain[0], PA_QUAD_WLEN + 1) ? "\033[1;32mPASS\033[0m" : "\033[1;31mFAIL\033[0m",
                            exp_pa_quad_sum);
    for (int32_t i = 0; i < PA_QUAD_WLEN+2; i++) printf("-----");
    printf("\n");

    // cleanup
    delete_LweSample_array(PA_QUAD_WLEN,     x_pa_quad);
    delete_LweSample_array(PA_QUAD_WLEN,     y_pa_quad);
    delete_LweSample_array(PA_QUAD_WLEN + 1, z_pa_quad);
#endif


#ifdef BS_TEST
    // -------------------------------------------------------------------------
    //
    //  Bootstrapping Test
    //
    printf("\n\n================================================================================\n");
    printf("\n    <<<<    Bootstrapping Test    >>>>\n\n");fflush(stdout);

    // setup pi
    const uint32_t pi = tfhe_params_store[BS_TFHE_PARAMS_INDEX].pi;
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

    //DBG
    printf("pi = %d\n", pi);

    // print table heading
    printf("--------------------------------------------------------------------------------\n");
    printf(" Encr -> Decr    | Id.  | <=> 3 | == 2 | timing (Id.)\n");
    printf("--------------------------------------------------------------------------------\n");

    for (int32_t i = 0; i < (1 << pi); i++)
    {
        // encrypt
        sym_encr_priv(a, i - (1 << (pi-1)), pi, bs_tfhe_keys);

        // bootstrap
        clock_t begin_id = clock();
        bs_id(id, a, pi, &(bs_tfhe_keys->cloud)); //FUCKUP #01: member 'cloud' is a struct, but not a pointer (unlike others)
        clock_t end_id = clock();

        // clock_t begin_gl = clock();
        bs_gleq(gl, a, pi, 3, &(bs_tfhe_keys->cloud));
        // clock_t end_gl = clock();

        // clock_t begin_eq = clock();
        bs_eq(eq, a, pi, 2, &(bs_tfhe_keys->cloud));
        // clock_t end_eq = clock();

        // decrypt
        int32_t a_plain     = sym_decr(a,  pi, bs_tfhe_keys);
        int32_t id_plain    = sym_decr(id, pi, bs_tfhe_keys);
        int32_t gl_plain    = sym_decr(gl, pi, bs_tfhe_keys);
        int32_t eq_plain    = sym_decr(eq, pi, bs_tfhe_keys);

        printf(" D[E(%+3d)] = %+3d |  %+3d |   %+d  |  %+d  | %lu ms\n", i - (1 << (pi-1)), a_plain,
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
#ifdef  PA_TEST_BIN
    delete_gate_bootstrapping_secret_keyset( pa_bin_tfhe_keys);
    delete_gate_bootstrapping_parameters(    pa_bin_tfhe_params);
#endif
#ifdef  PA_TEST_QUAD
    delete_gate_bootstrapping_secret_keyset( pa_quad_tfhe_keys);
    delete_gate_bootstrapping_parameters(    pa_quad_tfhe_params);
#endif
#ifdef  BS_TEST
    delete_gate_bootstrapping_secret_keyset( bs_tfhe_keys);
    delete_gate_bootstrapping_parameters(    bs_tfhe_params);
#endif

    return 0;
}
