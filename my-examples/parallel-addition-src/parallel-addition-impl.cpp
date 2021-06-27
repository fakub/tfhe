#include <cstdio>

#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"
#include "tfhe_garbage_collector.h"

#include <bootstrapping-functions.h>
#include <parallel-addition-impl.h>


// =============================================================================
//
//  Static Function Prototypes
//

// -----------------------------------------------------------------------------
//  Bin
//
/**
 *  @brief          Description
 *
 */
static void paral_calc_qi_bin(LweSample *qi,
                              const LweSample *w_i0,
                              const LweSample *w_i1,
                              const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
static void paral_calc_zi_bin(LweSample *zi,
                              const LweSample *w_i,
                              const LweSample *q_i0,
                              const LweSample *q_i1,
                              const TFheGateBootstrappingCloudKeySet *bk);

// -----------------------------------------------------------------------------
//  Quad
//
/**
 *  @brief          Description
 *
 */
static void paral_calc_qi_quad(LweSample *qi,
                               const LweSample *w_i0,
                               const LweSample *w_i1,
                               const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
static void paral_calc_zi_quad(LweSample *zi,
                               const LweSample *w_i,
                               const LweSample *q_i0,
                               const LweSample *q_i1,
                               const TFheGateBootstrappingCloudKeySet *bk);


// =============================================================================
//
//  Extern variables
//

//  ----    TFHE params moved to separate file    ----

extern const uint32_t PI_SEQ = tfhe_params_store[SEQ_TFHE_PARAMS_INDEX].pi;
extern const uint32_t PI_PAB = tfhe_params_store[PAB_TFHE_PARAMS_INDEX].pi;
extern const uint32_t PI_PAQ = tfhe_params_store[PAQ_TFHE_PARAMS_INDEX].pi;
extern const uint32_t PI_SGN = tfhe_params_store[SGN_TFHE_PARAMS_INDEX].pi;


// =============================================================================
//
//  Function Implementations
//

// -----------------------------------------------------------------------------
//  TFHE Param Setup
//
void setup_TFHE_params(const int tfhe_params_index,
                       TFheGateBootstrappingParameterSet **const tfhe_params)
{
    // choose TFHE param set from store
    const tfhe_params_t *const tfhe_params_set = &tfhe_params_store[tfhe_params_index];

    // setup TFHE Lib params
    LweParams *params_in = new_LweParams(tfhe_params_set->n,
                                         tfhe_params_set->ks_stdev,
                                         tfhe_params_set->max_stdev);
    TLweParams *params_accum = new_TLweParams(tfhe_params_set->N,
                                              tfhe_params_set->k,
                                              tfhe_params_set->bk_stdev,
                                              tfhe_params_set->max_stdev);
    TGswParams *params_bk = new_TGswParams(tfhe_params_set->bk_l,
                                           tfhe_params_set->bk_Bgbit,
                                           params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    *tfhe_params = new TFheGateBootstrappingParameterSet(tfhe_params_set->ks_length,
                                                         tfhe_params_set->ks_basebit,
                                                         params_in, params_bk);
}

// -----------------------------------------------------------------------------
//  En/Decryption
//
void sym_encr_priv(LweSample *ct,
                   const int32_t message,
                   const uint32_t pi,
                   const TFheGateBootstrappingSecretKeySet *sk)
{
    Torus32 _1s16 = modSwitchToTorus32(1, 1 << pi);

    // scale to [-8/16, 7/16]
    //                         ___/ 8 \___     _/ mask 1111 \_     ___/ 8 \___
    Torus32 mu = (((message + (1 << (pi-1))) & ((1 << pi) - 1)) - (1 << (pi-1))) * _1s16;
    double alpha = sk->params->in_out_params->alpha_min;
    lweSymEncrypt(ct, mu, alpha, sk->lwe_key);
}

void bin_sym_encr(LweSample *ct,
                  const int32_t message,
                  const TFheGateBootstrappingSecretKeySet *sk)
{
    if ((message < 0) || (1 < message))
        die_soon("Out of the alphabet A = {0,1}.");

    bootsSymEncrypt(ct, message, sk);
}

void seq_quad_sym_encr(LweSample *ct,
                       const int32_t message,
                       const TFheGateBootstrappingSecretKeySet *sk)
{
    if ((message < 0) || (3 < message))
        die_soon("Out of the alphabet A = [0 .. 3].");

    sym_encr_priv(ct, message, PI_SEQ, sk);
}

void paral_quad_sym_encr(LweSample *ct,
                         const int32_t message,
                         const TFheGateBootstrappingSecretKeySet *sk)
{
    if ((message < -2) || (2 < message))
        die_soon("Out of the alphabet A = [-2 .. 2].");

    sym_encr_priv(ct, message, PI_PAQ, sk);
}

void paral_bin_sym_encr(LweSample *ct,
                        const int32_t message,
                        const TFheGateBootstrappingSecretKeySet *sk)
{
    if ((message < -1) || (1 < message))
        die_soon("Out of the alphabet A = [-1 .. 1].");

    sym_encr_priv(ct, message, PI_PAB, sk);
}

void sgn_sym_encr(LweSample *ct,
                  const int32_t message,
                  const TFheGateBootstrappingSecretKeySet *sk)
{
    if ((message <= -(1 << (PI_SGN-1))) || ((1 << (PI_SGN-1)) <= message))
        die_soon("Out of the interval = (-2^(pi-1) .. 2^(pi-1)).");

    sym_encr_priv(ct, message, PI_SGN, sk);
}

int32_t sym_decr(const LweSample *ct,
                 const uint32_t pi,
                 const TFheGateBootstrappingSecretKeySet *sk)
{
    Torus32 _1s16 = modSwitchToTorus32(1, 1 << pi);
    Torus32 mu = lwePhase(ct, sk->lwe_key);
    return (mu + (_1s16 >> 1)) >> (32 - pi);
}

int32_t bin_sym_decr(const LweSample *ct,
                     const TFheGateBootstrappingSecretKeySet *sk)
{
    return bootsSymDecrypt(ct, sk) ? 1 : 0;
}

void bin_noiseless(LweSample *ct,
                   const int32_t message,
                   const TFheGateBootstrappingCloudKeySet *bk)
{
    if ((message < 0) || (1 < message))
        die_soon("Out of the alphabet A = {0,1}.");

    const LweParams *io_lwe_params = bk->params->in_out_params;

    Torus32 _1s8 = modSwitchToTorus32(1, 8);
    Torus32 mu = (message == 1) ? _1s8 : -_1s8;
    lweNoiselessTrivial(ct, mu, io_lwe_params);
}

// -----------------------------------------------------------------------------
//  Sequential Addition
//
void sequential_add(LweSample *z,
                    const LweSample *x,
                    const LweSample *y,
                    const uint32_t wlen,
#ifdef DBG_OUT
                    const TFheGateBootstrappingSecretKeySet *sk,
#endif
                    const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

#ifdef DBG_OUT
    printf("\n");
#endif

    //  SCENARIO binary with 5 bootstraps (ala TFHE Library)
#if SEQ_SCENARIO == A_CARRY_5_BS_BIN
{
    // alloc aux arrays
    LweSample *carry    = new_LweSample_array(2, io_lwe_params);
    LweSample *temp     = new_LweSample_array(3, io_lwe_params);

    // init c = 0
    bin_noiseless(carry, 0, bk);

    // calc z's and carry
    for (uint32_t i = 0; i < wlen; i++)
    {
#ifndef DBG_OUT
        // progress bar ...
        printf("-");fflush(stdout);
#endif

        // z_i = x_i XOR y_i XOR c
        bootsXOR(temp, x + i, y + i, bk);   // w_i = x_i XOR y_i (temp_0)
        bootsXOR(z + i, temp, carry, bk);

        // c = (x_i AND y_i) XOR (c AND (x_i XOR y_i))
        bootsAND(temp + 1, x + i, y + i, bk);   // temp_1 = x_i AND y_i
        bootsAND(temp + 2, carry, temp,  bk);   // temp_2 = c AND w_i
        bootsXOR(carry + 1, temp + 1, temp + 2, bk);
        bootsCOPY(carry, carry + 1, bk);

#ifdef DBG_OUT
        printf("Position #%d\n", i);
        int32_t x_plain = bin_sym_decr(x + i, sk);
        int32_t y_plain = bin_sym_decr(y + i, sk);
        int32_t c_plain = bin_sym_decr(carry, sk);
        int32_t z_plain = bin_sym_decr(z + i, sk);
        printf("    x = %d, y = %d, z = %d, c = %d\n", x_plain, y_plain, z_plain, c_plain);fflush(stdout);
#endif
    }

    // i = wlen
    bootsCOPY(z + wlen, carry, bk);

    // cleanup
    delete_LweSample_array(3, temp);
    delete_LweSample_array(2, carry);
}

    //  SCENARIO binary with 2 bootstraps
#elif SEQ_SCENARIO == B_CARRY_2_BS_BIN
{
    // alloc carry
    LweSample *c = new_LweSample(io_lwe_params);

    // init c = 0
    bin_noiseless(c, 0, bk);

    // calc z's and carry
    for (uint32_t i = 0; i < wlen; i++)
    {
#ifndef DBG_OUT
        // progress bar ...
        printf("-");fflush(stdout);
#endif

        // z_i = XOR3(x_i, y_i, c)
        bootsXOR3(z + i, x + i, y + i, c, bk);
        // c = 2OF3(x_i, y_i, c)
        boots2OF3(c,     x + i, y + i, c, bk);

#ifdef DBG_OUT
        printf("Position #%d\n", i);
        int32_t x_plain = bin_sym_decr(x + i, sk);
        int32_t y_plain = bin_sym_decr(y + i, sk);
        int32_t c_plain = bin_sym_decr(c,     sk);
        int32_t z_plain = bin_sym_decr(z + i, sk);
        printf("    x = %d, y = %d, z = %d, c = %d\n", x_plain, y_plain, z_plain, c_plain);fflush(stdout);
#endif
    }

    // i = wlen
    bootsCOPY(z + wlen, c, bk);

    // cleanup
    delete_LweSample(c);
}

    //  SCENARIO quad with 5 bootstraps
#elif SEQ_SCENARIO == C_CARRY_2_BIT
{
    // alloc carry & w
    LweSample *c = new_LweSample(io_lwe_params);
    LweSample *w = new_LweSample(io_lwe_params);

    // init c = 0
    lweNoiselessTrivial(c, 0, io_lwe_params);

    // calc z's and carry
    for (uint32_t i = 0; i < wlen; i++)
    {
#ifndef DBG_OUT
        // progress bar ...
        printf("-");fflush(stdout);
#endif

        // w = x_i + y_i + c
        // init w = 0
        lweNoiselessTrivial(w, 0, io_lwe_params);
        // w += x_i
        lweAddTo(w, x + i, io_lwe_params);
        // w += y_i
        lweAddTo(w, y + i, io_lwe_params);
        // w += c
        lweAddTo(w,     c, io_lwe_params);

        // c = (0,0,0,0,1,1,1,1)[w]   (using trick with 1/2)
        bs_01(c, w, PI_SEQ, bk);

        // z_i = w_i - 4 c
        lweSubMulTo(w, 4, c, io_lwe_params);

        // bootstrap with identity (at positive interval)
        bs_pos_id(z + i, w, PI_SEQ, bk);

#ifdef DBG_OUT
        printf("Position #%d\n", i);
        int32_t x_plain = sym_decr(x + i, PI_SEQ, sk);
        int32_t y_plain = sym_decr(y + i, PI_SEQ, sk);
        int32_t c_plain = sym_decr(c,     PI_SEQ, sk);
        int32_t z_plain = sym_decr(z + i, PI_SEQ, sk);
        printf("    x = %d, y = %d, z = %d, c = %d\n", x_plain, y_plain, z_plain, c_plain);fflush(stdout);
#endif
    }

    // i = wlen
    lweCopy(z + wlen, c, io_lwe_params);

    // cleanup
    delete_LweSample(w);
    delete_LweSample(c);
}

#else
    #pragma message "Invalid sequential addition scenario " XSTR(SEQ_SCENARIO)
    #error "As stated above."
#endif

    // progress bar end
    printf("\n");
}

// -----------------------------------------------------------------------------
//  Parallel Addition (Bin)
//
void parallel_add_bin(LweSample *z,
                      const LweSample *x,
                      const LweSample *y,
                      const uint32_t wlen_bin,
#ifdef DBG_OUT
                      const TFheGateBootstrappingSecretKeySet *sk,
#endif
                      const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

#ifdef DBG_OUT
    printf("\n");
#endif

    // alloc aux arrays
    LweSample *w = new_LweSample_array(wlen_bin, io_lwe_params);
    LweSample *q = new_LweSample_array(wlen_bin, io_lwe_params);

    // calc w_i = x_i + y_i
    for (uint32_t i = 0; i < wlen_bin; i++)
    {
        lweCopy (w + i, x + i, io_lwe_params);
        lweAddTo(w + i, y + i, io_lwe_params);
    }

    // calc q_i's
    // i = 0
    paral_calc_qi_bin(q, w, NULL, bk);
    // i > 0
    for (uint32_t i = 1; i < wlen_bin; i++)
    {
        paral_calc_qi_bin(q + i, w + i, w + i - 1, bk);
    }

    // calc z_i's
    // i = 0
    paral_calc_zi_bin(z, w, q, NULL, bk);
    // i > 0
    for (uint32_t i = 1; i < wlen_bin; i++)
    {
        paral_calc_zi_bin(z + i, w + i, q + i, q + i - 1, bk);
    }
    // i = wlen_bin
    paral_calc_zi_bin(z + wlen_bin, NULL, NULL, q + wlen_bin - 1, bk);

#ifdef DBG_OUT
    for (uint32_t i = 0; i < wlen_bin; i++)
    {
        printf("Position #%d\n", i);
        int32_t x_plain = sym_decr(x + i, PI_PAB, sk);
        int32_t y_plain = sym_decr(y + i, PI_PAB, sk);
        int32_t z_plain = sym_decr(z + i, PI_PAB, sk);
        int32_t w_plain = sym_decr(w + i, PI_PAB, sk);
        int32_t q_plain = sym_decr(q + i, PI_PAB, sk);
        printf("    x = %d, y = %d, z = %d, w = %d, q = %d\n", x_plain, y_plain, z_plain, w_plain, q_plain);fflush(stdout);
    }
#endif

    // progress bar end
    printf("\n");

    // cleanup
    delete_LweSample_array(wlen_bin, q);
    delete_LweSample_array(wlen_bin, w);
}

// -----------------------------------------------------------------------------
//  Parallel Addition (Quad)
//
void parallel_add_quad(LweSample *z,
                       const LweSample *x,
                       const LweSample *y,
                       const uint32_t wlen_quad,
#ifdef DBG_OUT
                       const TFheGateBootstrappingSecretKeySet *sk,
#endif
                       const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

#ifdef DBG_OUT
    printf("\n");
#endif

    // alloc aux arrays
    LweSample *w = new_LweSample_array(wlen_quad, io_lwe_params);
    LweSample *q = new_LweSample_array(wlen_quad, io_lwe_params);

    // calc w_i = x_i + y_i
    for (uint32_t i = 0; i < wlen_quad; i++)
    {
        lweCopy (w + i, x + i, io_lwe_params);
        lweAddTo(w + i, y + i, io_lwe_params);
    }

    // calc q_i's
    // i = 0
    paral_calc_qi_quad(q, w, NULL, bk);
    // i > 0
    for (uint32_t i = 1; i < wlen_quad; i++)
    {
        paral_calc_qi_quad(q + i, w + i, w + i - 1, bk);
    }

    // calc z_i's
    // i = 0
    paral_calc_zi_quad(z, w, q, NULL, bk);
    // i > 0
    for (uint32_t i = 1; i < wlen_quad; i++)
    {
        paral_calc_zi_quad(z + i, w + i, q + i, q + i - 1, bk);
    }
    // i = wlen_quad
    paral_calc_zi_quad(z + wlen_quad, NULL, NULL, q + wlen_quad - 1, bk);

#ifdef DBG_OUT
    for (uint32_t i = 0; i < wlen_quad; i++)
    {
        printf("Position #%d\n", i);
        int32_t x_plain = sym_decr(x + i, PI_PAQ, sk);
        int32_t y_plain = sym_decr(y + i, PI_PAQ, sk);
        int32_t z_plain = sym_decr(z + i, PI_PAQ, sk);
        int32_t w_plain = sym_decr(w + i, PI_PAQ, sk);
        int32_t q_plain = sym_decr(q + i, PI_PAQ, sk);
        printf("    x = %d, y = %d, z = %d, w = %d, q = %d\n", x_plain, y_plain, z_plain, w_plain, q_plain);fflush(stdout);
    }
#endif

    // progress bar end
    printf("\n");

    // cleanup
    delete_LweSample_array(wlen_quad, q);
    delete_LweSample_array(wlen_quad, w);
}

// -----------------------------------------------------------------------------
//  Parallel Signum
//
void parallel_sgn(LweSample *sgn,
                  const LweSample *x,
                  const uint32_t wlen_sgn,
#ifdef DBG_OUT
                  const TFheGateBootstrappingSecretKeySet *sk,
#endif
                  const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

#ifdef DBG_OUT
    printf("\n");
#endif

    if (wlen_sgn == 1)
    {
        // bootstrap one more time (this is actually signum) & return
        bs_gleq(sgn, x, 1, 1, PI_SGN, bk);

#ifdef DBG_OUT
        printf("Last step\n");
        int32_t  x_plain = sym_decr(x,   PI_SGN, sk);
        int32_t  s_plain = sym_decr(sgn, PI_SGN, sk);
        printf("    sgn(%d) = %d\n", x_plain, s_plain);fflush(stdout);
#endif

        // progress bar end
        printf("\n");

        return;
    }

    const uint32_t gamma = PI_SGN - 1;
    // calc new length (ceiled wlen_sgn / gamma)
    const uint32_t new_wlen = (wlen_sgn-1) / gamma + 1;

#ifdef DBG_OUT
    printf("Size = %d\n", wlen_sgn);
#endif

    // alloc aux variables / arrays
    LweSample *b = new_LweSample_array(new_wlen, io_lwe_params);
    LweSample *s = new_LweSample(io_lwe_params);

    for (uint32_t j = 0; j < new_wlen; j++)
    {
#ifdef DBG_OUT
        //~ int32_t  x2i_plain = sym_decr(x+2*i, PI_SGN, sk);
        //~ sym_encr_priv(x+2*i, x2i_plain, PI_SGN, sk);
        //~ int32_t x2i1_plain = 0;
        //~ LweSample *z1    = new_LweSample(io_lwe_params);    // additive 1
        //~ LweSample *xi_an = new_LweSample(io_lwe_params);    // x_i for analysis
        //~ LweSample *r0_an = new_LweSample(io_lwe_params);    // r_0 for analysis
        //~ sym_encr_priv(z1, 1, PI_SGN, sk);
#else
        // progress bar ...
        printf("-");fflush(stdout);
#endif

        // init b + j = 0
        lweNoiselessTrivial(b + j, 0, io_lwe_params);

        for (uint32_t i = 0; i < gamma; i++)
        {
            //         s  = (x_(gamma*j+i) <=> +-1) * 2^i
            //TODO set to zero if out of range
            if (gamma*j+i >= wlen_sgn)
            {
                lweNoiselessTrivial(s, 0, io_lwe_params);
            }
            else
            {
                bs_gleq(   s,    x +gamma*j+i,       1,  (1<<i), PI_SGN, bk);
            }
            //       b_j += s
            lweAddTo(b+j,   s, io_lwe_params);
#ifdef DBG_OUT
            //~ int32_t s_plain = sym_decr(s,    PI_SGN, sk);
#endif
        }

#ifdef DBG_OUT
        printf("    Position #%d\n", j);
        int32_t  bj_plain = sym_decr(b + j, PI_SGN, sk);
        printf("        b_j = %+d\n", bj_plain);fflush(stdout);

        //~ for (int32_t k = -4; k <= 4; k++)
        //~ {
            //~ lweCopy(    xi_an, x + 2*i, io_lwe_params);
            //~ lweAddMulTo(xi_an, k, z1, io_lwe_params);
            //~ bs_gleq(r0_an, xi_an, 1, 2, PI_SGN, bk);
            //~ int32_t r0_an_plain = sym_decr(r0_an, PI_SGN, sk);
            //~ printf("            k = %+d: r_0 = %+d\n", k, r0_an_plain);
        //~ }
#endif
    }

    // call recurently
    parallel_sgn(sgn, b, new_wlen,
#ifdef DBG_OUT
                 sk,
#endif
                 bk);

    // cleanup
    delete_LweSample_array(new_wlen, b);
    delete_LweSample(s);
}

// -----------------------------------------------------------------------------
//  Misc
//
void die_soon(const char* message)
{
    fprintf(stderr, "(!) %s\n    Aborting ...\n", message);
    abort();
}

int64_t quad_eval(const int32_t *const x,
                  const uint32_t len)
{
    int64_t r = 0;

    for (uint32_t i = 0; i < len; i++)
    {
        r += x[i] * (1 << (2*i));
    }

    return r;
}

int64_t bin_eval(const int32_t *const x,
                 const uint32_t len)
{
    int64_t r = 0;

    for (uint32_t i = 0; i < len; i++)
    {
        r += x[i] * (1 << i);
    }

    return r;
}

int32_t sgn_eval(const int32_t *const x,
                 const uint32_t len)
{
    if (len == 0) return 0;

    int32_t i = len;

    do
    {
        i--;
    } while (i > 0 && x[i] == 0);

    return (x[i] > 0) ? 1 : (x[i] < 0 ? -1 : 0);
}


// =============================================================================
//
//  Static Function Implementations
//

// -----------------------------------------------------------------------------
//  Bin
//
static void paral_calc_qi_bin(LweSample *qi,
                              const LweSample *w_i0,
                              const LweSample *w_i1,
                              const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    // progress bar ...
    printf("-");fflush(stdout);

    //  SCENARIO IIa-D
#if PA_SCENARIO_BIN == D_PARALLEL_2
    {
        // aux variables
        LweSample *r1 = new_LweSample(io_lwe_params);
        LweSample *r2 = new_LweSample(io_lwe_params);
        LweSample *r3 = new_LweSample(io_lwe_params);
        LweSample *r23= new_LweSample(io_lwe_params);

        //      r1   =   w_i <=> +-2
        bs_gleq(r1,      w_i0,     2,   1, PI_PAB, bk);

        //      r2   =   w_i == +-1
        bs_eq  (r2,      w_i0,    1,    PI_PAB, bk);

        if (w_i1 == NULL)
        {
            //  zero
            lweNoiselessTrivial(r3, 0,  io_lwe_params);
        }
        else
        {
            //      r3   =   w_i-1 <=> +-1
            bs_gleq(r3,      w_i1,       1,   1, PI_PAB, bk);
        }

        //      r23: r2+r3 == +-2
        lweAddTo(    r2,r3,         io_lwe_params);
        bs_eq  (r23, r2,        2,  PI_PAB, bk);

        //      q_i   =   r1 + r23
        lweCopy (qi,      r1,       io_lwe_params);
        lweAddTo(qi,           r23, io_lwe_params);
    }

    //  SCENARIO IIa-E
#elif PA_SCENARIO_BIN == E_PARALLEL_2
    {
        // aux variables
        LweSample *r1 = new_LweSample(io_lwe_params);
        LweSample *r2 = new_LweSample(io_lwe_params);
        LweSample *w_i1__2_r2 = new_LweSample(io_lwe_params);
        LweSample *r23= new_LweSample(io_lwe_params);

        //      r1   =   w_i <=> +-2
        bs_gleq(r1,      w_i0,     2,   1, PI_PAB, bk);

        //      r2   =   w_i == +-1
        bs_eq  (r2,      w_i0,    1,    PI_PAB, bk);

        if (w_i1 == NULL)
        {
            lweNoiselessTrivial(w_i1__2_r2, 0,  io_lwe_params);             // w_i1__2_r2 = w_i-1
        }
        else
        {
            lweCopy(w_i1__2_r2, w_i1, io_lwe_params);                       // w_i1__2_r2 = w_i-1
        }
        //      r23   =   w_i-1     + 2 r2 <=> +-3
        lweAddMulTo(      w_i1__2_r2, 2,r2,         io_lwe_params);         // w_i1__2_r2 = w_i-1 + 2 r2
        bs_gleq(r23,      w_i1__2_r2,            3, 1, PI_PAB, bk);

        //      q_i   =   r1 + r23
        lweCopy (qi,      r1,       io_lwe_params);
        lweAddTo(qi,           r23, io_lwe_params);

        //TODO cleanup
    }

    //  SCENARIO IIa-F
#elif PA_SCENARIO_BIN == F_PARALLEL_2
    {
        // aux variables
        LweSample *r1 = new_LweSample(io_lwe_params);

        if (w_i1 == NULL)
        {
            lweNoiselessTrivial(r1, 0, io_lwe_params);              // r1 = w_i-1
        }
        else
        {
            lweCopy(r1, w_i1, io_lwe_params);                       // r1 = w_i-1
        }

        //      q_i   =   w_i-1 + 3 w_i <=> +-4
        lweAddMulTo(      r1,     3,w_i0,          io_lwe_params); // r1 = w_i-1 + 3 w_i
        bs_gleq(qi,       r1,                 4,   1, PI_PAB, bk);
    }

#else
    #pragma message "Invalid parallel addition scenario " XSTR(PA_SCENARIO_QUAD)
    #error "As stated above."
#endif
}

static void paral_calc_zi_bin(LweSample *zi,
                              const LweSample *w_i,
                              const LweSample *q_i0,
                              const LweSample *q_i1,
                              const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    // aux variables
    LweSample *tmpz = new_LweSample(io_lwe_params);

    // setup with w_i (or 0)
    if (w_i == NULL)
    {
        lweNoiselessTrivial(tmpz, 0, io_lwe_params);
    }
    else
    {
        lweCopy(tmpz, w_i, io_lwe_params);
    }

    // subtract 2 q_i
    if (q_i0 != NULL)
    {
        lweSubMulTo(tmpz, 2, q_i0, io_lwe_params);
    }

    // add q_i-1
    if (q_i1 != NULL)
    {
        lweAddTo(tmpz, q_i1, io_lwe_params);
    }

    // progress bar ...
    printf("-");fflush(stdout);

    // bootstrap to refresh noise
    bs_id(zi, tmpz, PI_PAB, bk);

    // cleanup
    delete_LweSample(tmpz);
}

// -----------------------------------------------------------------------------
//  Quad
//
static void paral_calc_qi_quad(LweSample *qi,
                               const LweSample *w_i0,
                               const LweSample *w_i1,
                               const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    // progress bar ...
    printf("-");fflush(stdout);

    //  SCENARIO IIb-G
#if PA_SCENARIO_QUAD == G_PARALLEL_4
    {
        // aux variables
        LweSample *r1 = new_LweSample(io_lwe_params);
        LweSample *r2 = new_LweSample(io_lwe_params);
        LweSample *r3 = new_LweSample(io_lwe_params);
        LweSample *r23= new_LweSample(io_lwe_params);

        //      r1   =   w_i <=> +-3
        bs_gleq(r1,      w_i0,     3,   1, PI_PAQ, bk);

        //      r2   =   w_i == +-2
        bs_eq  (r2,      w_i0,    2,    PI_PAQ, bk);

        if (w_i1 == NULL)
        {
            //  zero
            lweNoiselessTrivial(r3, 0,  io_lwe_params);
        }
        else
        {
            //      r3   =   w_i-1 <=> +-2
            bs_gleq(r3,      w_i1,       2,   1, PI_PAQ, bk);
        }

        //      r23: r2+r3 == +-2
        lweAddTo(    r2,r3,         io_lwe_params);
        bs_eq  (r23, r2,        2,  PI_PAQ, bk);

        //      q_i   =   r1 + r23
        lweCopy (qi,      r1,       io_lwe_params);
        lweAddTo(qi,           r23, io_lwe_params);
    }

    //  SCENARIO IIb-H
#elif PA_SCENARIO_QUAD == H_PARALLEL_4
    {
        // aux variables
        LweSample *r1 = new_LweSample(io_lwe_params);
        LweSample *r2 = new_LweSample(io_lwe_params);
        LweSample *w_i1__3_r2 = new_LweSample(io_lwe_params);
        LweSample *r23= new_LweSample(io_lwe_params);

        //      r1   =   w_i <=> +-3
        bs_gleq(r1,      w_i0,     3,   1, PI_PAQ, bk);

        //      r2   =   w_i == +-2
        bs_eq  (r2,      w_i0,    2,    PI_PAQ, bk);

        if (w_i1 == NULL)
        {
            lweNoiselessTrivial(w_i1__3_r2, 0,  io_lwe_params);             // w_i1__3_r2 = w_i-1
        }
        else
        {
            lweCopy(w_i1__3_r2, w_i1, io_lwe_params);                       // w_i1__3_r2 = w_i-1
        }
        //      r23   =   w_i-1     + 3 r2 <=> +-5
        lweAddMulTo(      w_i1__3_r2, 3,r2,         io_lwe_params);         // w_i1__3_r2 = w_i-1 + 3 r2
        bs_gleq(r23,      w_i1__3_r2,            5, 1, PI_PAQ, bk);

        //      q_i   =   r1 + r23
        lweCopy (qi,      r1,       io_lwe_params);
        lweAddTo(qi,           r23, io_lwe_params);
    }

    //  SCENARIO IIb-I
#elif PA_SCENARIO_QUAD == I_PARALLEL_4
    {
        // aux variables
        LweSample *r1 = new_LweSample(io_lwe_params);

        if (w_i1 == NULL)
        {
            lweNoiselessTrivial(r1, 0, io_lwe_params);              // r1 = w_i-1
        }
        else
        {
            lweCopy(r1, w_i1, io_lwe_params);                       // r1 = w_i-1
        }

        //      q_i   =   w_i-1 + 6 w_i <=> +-14
        lweAddMulTo(      r1,     6,w_i0,           io_lwe_params); // r1 = w_i-1 + 6 w_i
        bs_gleq(qi,       r1,                 14,   1, PI_PAQ, bk);
    }

#else
    #pragma message "Invalid parallel addition scenario " XSTR(PA_SCENARIO_QUAD)
    #error "As stated above."
#endif
}

static void paral_calc_zi_quad(LweSample *zi,
                               const LweSample *w_i,
                               const LweSample *q_i0,
                               const LweSample *q_i1,
                               const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    // aux variables
    LweSample *tmpz = new_LweSample(io_lwe_params);

    // setup with w_i (or 0)
    if (w_i == NULL)
    {
        lweNoiselessTrivial(tmpz, 0, io_lwe_params);
    }
    else
    {
        lweCopy(tmpz, w_i, io_lwe_params);
    }

    // subtract 4 q_i
    if (q_i0 != NULL)
    {
        lweSubMulTo(tmpz, 4, q_i0, io_lwe_params);
    }

    // add q_i-1
    if (q_i1 != NULL)
    {
        lweAddTo(tmpz, q_i1, io_lwe_params);
    }

#ifndef DBG_OUT
    // progress bar ...
    printf("-");fflush(stdout);
#endif

    // bootstrap to refresh noise
    bs_id(zi, tmpz, PI_PAQ, bk);

    // cleanup
    delete_LweSample(tmpz);
}
