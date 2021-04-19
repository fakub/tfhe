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

#define BS_TEST
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

    // orig. params
    //~ static const int32_t N = 1024;
    //~ static const int32_t k = 1;
    //~ static const int32_t n = 630;
    //~ static const int32_t bk_l = 3;
    //~ static const int32_t bk_Bgbit = 7;
    //~ static const int32_t ks_basebit = 2;
    //~ static const int32_t ks_length = 8;
    //~ static const double ks_stdev = pow(2.,-15); // standard deviation
    //~ static const double bk_stdev = pow(2.,-25); // standard deviation
    //~ static const double max_stdev = 0.012467;   // max standard deviation for a 1/4 msg space

    // w/o KS: DEPRECATED
    // for  80-bit security @ pi = 4
    //      N =  512, n = 682, l = 20, mlogal = 23.745
    // for 128-bit security @ pi = 4
    //      N = 1024, n = 768, l = 20, mlogal = 24.331

    // with KS:
    // -------------------------------------------------------------------------
    // pi = 2 ; delta2 = 1 ; n = 384 ; nn = 1024 ; gamma = 15
    //   # 128-bit sec., n_max 43k
    //   # t = 11; l = 1; -log(a_KS_n) = 13.31; -log(a_BK_N) = 31.17    (little more than 128-bit)

        //~ static const int32_t N = 1024;
        //~ static const int32_t k = 1;
        //~ static const int32_t n = 384;
        //~ static const int32_t bk_l = 1;
        //~ static const int32_t bk_Bgbit = 15;
        //~ static const int32_t ks_basebit = 1;
        //~ static const int32_t ks_length = 11;
        //~ static const double ks_stdev = pow(2.,-13); // standard deviation
        //~ static const double bk_stdev = pow(2.,-31); // standard deviation
        //~ static const double max_stdev = 0.04167;    // max standard deviation for a 1/4 msg space

    // pi = 2 ; delta2 = Math.log2(3) ; n = 384 ; nn = 1024 ; gamma = 16
    //   # ~128-bit sec.,  n_max 43k,
    //   # t = 11; l = 1; -log(a_KS_n) = 13.61; -log(a_BK_N) = 32.46    (close to 128-bit)

        //~ static const int32_t N = 1024;
        //~ static const int32_t k = 1;
        //~ static const int32_t n = 384;
        //~ static const int32_t bk_l = 1;
        //~ static const int32_t bk_Bgbit = 16;
        //~ static const int32_t ks_basebit = 1;
        //~ static const int32_t ks_length = 11;
        //~ static const double ks_stdev = pow(2.,-13); // standard deviation
        //~ static const double bk_stdev = pow(2.,-33); // standard deviation
        //~ static const double max_stdev = 0.04167;    // max standard deviation for a 1/4 msg space

                // it turns out that decreasing l is very important, however, this is only possible with larger N
                // pi = 2 ; delta2 = Math.log2(3) ; n = 384 ; nn = 512 ; gamma = 1
                //   # ~100-bit sec.,  n_max 11k,
                //   # t = 11; l = 15; -log(a_KS_n) = 13.11; -log(a_BK_N) = 18.92    (close to 100-bit)

                    //~ static const int32_t N = 512;
                    //~ static const int32_t k = 1;
                    //~ static const int32_t n = 384;
                    //~ static const int32_t bk_l = 15;
                    //~ static const int32_t bk_Bgbit = 1;
                    //~ static const int32_t ks_basebit = 1;
                    //~ static const int32_t ks_length = 11;
                    //~ static const double ks_stdev = pow(2.,-13); // standard deviation
                    //~ static const double bk_stdev = pow(2.,-19); // standard deviation
                    //~ static const double max_stdev = 0.04167;    // max standard deviation for a 1/4 msg space

    // -------------------------------------------------------------------------
    // pi = 4 ; delta2 = Math.log2(3)  ; n = 520 ; nn = 1024 ; gamma = 9
    //   # ~128-bit sec.,  n_max 2729.7
    //   # t = 13; l = 2; -log(a_KS_n) = 15.73; -log(a_BK_N) = 28.18    (more than 128-bit)

        //~ static const int32_t N = 1024;
        //~ static const int32_t k = 1;
        //~ static const int32_t n = 520;
        //~ static const int32_t bk_l = 2;
        //~ static const int32_t bk_Bgbit = 9;
        //~ static const int32_t ks_basebit = 1;
        //~ static const int32_t ks_length = 15;
        //~ static const double ks_stdev = pow(2.,-18); // standard deviation
        //~ static const double bk_stdev = pow(2.,-31); // standard deviation
        //~ static const double max_stdev = 0.01042;    // max standard deviation for a 1/4 msg space

    // pi = 4 ; delta2 = Math.log2(36)  ; n = 520 ; nn = 1024 ; gamma = 10
    //   # ~128-bit sec.,  n_max 2729.7
    //   # t = 15; l = 2; -log(a_KS_n) = 17.62; -log(a_BK_N) = 30.97    (little more than 128-bit)

                                                // orig. params
        //~ static const int32_t N = 1024;
        //~ static const int32_t k = 1;
        //~ static const int32_t n = 520;           // 630
        //~ static const int32_t bk_l = 2;         //   3
        //~ static const int32_t bk_Bgbit = 10;      //   7
        //~ static const int32_t ks_basebit = 1;
        //~ static const int32_t ks_length = 15;    //   8
        //~ static const double ks_stdev = pow(2.,-18); // standard deviation   // -15
        //~ static const double bk_stdev = pow(2.,-31); // standard deviation
        //~ static const double max_stdev = 0.01042;    // max standard deviation for a 1/4 msg space   // 0.012467

    // -------------------------------------------------------------------------
    // pi = 5 ; delta2 = Math.log2(36)  ; n = 550 ; nn = 1024 ; gamma = 11
    //   # ~128-bit sec.,  n_max 682
    //   # t = 16; l = 2; -log(a_KS_n) = 18.67; -log(a_BK_N) = 33.01    (very close to >128-bit)

        static const int32_t N = 1024;
        static const int32_t k = 1;
        static const int32_t n = 550;
        static const int32_t bk_l = 2;
        static const int32_t bk_Bgbit = 11;
        static const int32_t ks_basebit = 1;
        static const int32_t ks_length = 16;
        static const double ks_stdev = pow(2.,-19); // standard deviation
        static const double bk_stdev = pow(2.,-33); // standard deviation
        static const double max_stdev = 0.005208;   // max standard deviation for a 1/4 msg space

    // -------------------------------------------------------------------------
    // pi = 7 ; delta2 = Math.log2(74)  ; n = 682 ; nn = 4096 ; gamma = 24
    //   # ~128-bit sec.,  n_max 682
    //   # t = 20; l = 1; -log(a_KS_n) = 22.35; -log(a_BK_N) = 49.19    (more than 128-bit)
    //   !! might be problem with float precision !! only 53 bits

        //~ static const int32_t N = 4096;
        //~ static const int32_t k = 1;
        //~ static const int32_t n = 682;
        //~ static const int32_t bk_l = 1;
        //~ static const int32_t bk_Bgbit = 24;
        //~ static const int32_t ks_basebit = 1;
        //~ static const int32_t ks_length = 20;
        //~ static const double ks_stdev = pow(2.,-22); // standard deviation
        //~ static const double bk_stdev = pow(2.,-49); // standard deviation
        //~ static const double max_stdev = 0.001302;   // max standard deviation for a 1/4 msg space

    // -------------------------------------------------------------------------

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
    printf("\n    <<<<    Bootstrapping Test    >>>>\n\n");fflush(stdout);

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
    printf("\n    <<<<    Parallel Addition Test    >>>>\n\n");fflush(stdout);

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
