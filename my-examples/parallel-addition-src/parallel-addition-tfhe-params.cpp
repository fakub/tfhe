#include <parallel-addition-impl.h>


// =============================================================================
//
//  Extern variables
//

const tfhe_params_t tfhe_params_store[N_PARAM_SETS] =
{
    {   // TFHE_LIB
        .pi = 2,
        .N = 1024,
        .k = 1,
        .n = 630,
        .bk_l = 3,
        .bk_Bgbit = 7,
        .ks_basebit = 2,
        .ks_length = 8,
#ifdef NO_NOISE
        .ks_stdev = 0.0, .bk_stdev = 0.0,
#else
        .ks_stdev = pow(2.,-15),
        .bk_stdev = pow(2.,-25),
#endif
        .max_stdev = 0.012467,
    },
    {   // A_CARRY_5_BS_BIN                 55 ms
        .pi = 2,
        .N = 1024,
        .k = 1,
        .n = 400,
        .bk_l = 1,
        .bk_Bgbit = 15,
        .ks_basebit = 1,
        .ks_length = 11,
#ifdef NO_NOISE
        .ks_stdev = 0.0, .bk_stdev = 0.0,
#else
        .ks_stdev = pow(2.,-13.31),
        .bk_stdev = pow(2.,-31.20),
#endif
        .max_stdev = 0.04167,
    },
    {   // B_CARRY_2_BS_BIN                 60 ms
        .pi = 2,
        .N = 1024,
        .k = 1,
        .n = 420,
        .bk_l = 1,
        .bk_Bgbit = 16,
        .ks_basebit = 1,
        .ks_length = 11,
#ifdef NO_NOISE
        .ks_stdev = 0.0, .bk_stdev = 0.0,
#else
        .ks_stdev = pow(2.,-13.61),
        .bk_stdev = pow(2.,-32.53),
#endif
        .max_stdev = 0.04167,
    },
    {   // C_CARRY_2_BIT                    95 ms
        .pi = 3,
        .N = 1024,
        .k = 1,
        .n = 490,
        .bk_l = 2,
        .bk_Bgbit = 9,
        .ks_basebit = 1,
        .ks_length = 14,
#ifdef NO_NOISE
        .ks_stdev = 0.0, .bk_stdev = 0.0,
#else
        .ks_stdev = pow(2.,-16.11),
        .bk_stdev = pow(2.,-28.47),
#endif
        .max_stdev = 0.02084,
    },
    {   // D_PARALLEL_2                     98 ms
        .pi = 3,
        .N = 1024,
        .k = 1,
        .n = 480,
        .bk_l = 2,
        .bk_Bgbit = 9,
        .ks_basebit = 1,
        .ks_length = 13,
#ifdef NO_NOISE
        .ks_stdev = 0.0, .bk_stdev = 0.0,
#else
        .ks_stdev = pow(2.,-15.73),
        .bk_stdev = pow(2.,-28.12),
#endif
        .max_stdev = 0.02084,
    },
    {   // E_PARALLEL_2                     103 ms
        .pi = 4,
        .N = 1024,
        .k = 1,
        .n = 510,
        .bk_l = 2,
        .bk_Bgbit = 10,
        .ks_basebit = 1,
        .ks_length = 14,
#ifdef NO_NOISE
        .ks_stdev = 0.0, .bk_stdev = 0.0,
#else
        .ks_stdev = pow(2.,-16.78),
        .bk_stdev = pow(2.,-30.17),
#endif
        .max_stdev = 0.01042,
    },
    {   // F_PARALLEL_2                     111 ms
        .pi = 5,
        .N = 1024,
        .k = 1,
        .n = 560,
        .bk_l = 2,
        .bk_Bgbit = 10,
        .ks_basebit = 1,
        .ks_length = 16,
#ifdef NO_NOISE
        .ks_stdev = 0.0, .bk_stdev = 0.0,
#else
        .ks_stdev = pow(2.,-18.25),
        .bk_stdev = pow(2.,-31.60),
#endif
        .max_stdev = 0.005208,
    },
    {   // G_PARALLEL_4                     106 ms
        .pi = 4,
        .N = 1024,
        .k = 1,
        .n = 540,
        .bk_l = 2,
        .bk_Bgbit = 10,
        .ks_basebit = 1,
        .ks_length = 15,
#ifdef NO_NOISE
        .ks_stdev = 0.0, .bk_stdev = 0.0,
#else
        .ks_stdev = pow(2.,-17.62),
        .bk_stdev = pow(2.,-31.00),
#endif
        .max_stdev = 0.01042,
    },
    {   // H_PARALLEL_4                     113 ms
        .pi = 5,
        .N = 1024,
        .k = 1,
        .n = 570,
        .bk_l = 2,
        .bk_Bgbit = 11,
        .ks_basebit = 1,
        .ks_length = 16,
#ifdef NO_NOISE
        .ks_stdev = 0.0, .bk_stdev = 0.0,
#else
        .ks_stdev = pow(2.,-18.67),
        .bk_stdev = pow(2.,-33.04),
#endif
        .max_stdev = 0.005208,
    },
    {   // I_PARALLEL_4                     547 ms, but erroneous results !! (and the problem is not at N = 4096 FFT processor,
        //                                          probably 2^-49 is just too close to the precision of double)
        .pi = 7,
        .N = 4096,
        .k = 1,
        .n = 680,
        .bk_l = 1,
        .bk_Bgbit = 24,
        .ks_basebit = 1,
        .ks_length = 20,
#ifdef NO_NOISE
        .ks_stdev = 0.0, .bk_stdev = 0.0,
#else
        .ks_stdev = pow(2.,-22.35),
        .bk_stdev = pow(2.,-49.19),
#endif
        .max_stdev = 0.001302,
    },
};
