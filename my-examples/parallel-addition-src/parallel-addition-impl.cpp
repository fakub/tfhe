#include <cstdio>

#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"

#include <parallel-addition-impl.h>


// =============================================================================
//
//  Static Function Prototypes
//

/**
 *  @brief          Description
 *
 */
static void bs_set_tv_identity(Torus32 *const tv,
                               const uint32_t N,
                               const Torus32 MU);

/**
 *  @brief          Description
 *
 */
static void bs_set_tv_gleq(Torus32 *const tv,
                           const uint32_t N,
                           const uint32_t thr,
                           const Torus32 MU);

/**
 *  @brief          Description
 *
 */
static void bs_set_tv_eq(Torus32 *const tv,
                         const uint32_t N,
                         const uint32_t thr,
                         const Torus32 MU);

/**
 *  @brief          Description
 *
 */
static void bs_with_tv(LweSample *result,
                       const LweSample *sample,
                       const TFheGateBootstrappingCloudKeySet *bk,
                       const TorusPolynomial *testvect);

/**
 *  @brief          Description
 *
 */
static void paral_calc_qi(LweSample *qi,
                          const LweSample *w_i0,
                          const LweSample *w_i1,
                          const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
static void paral_calc_zi(LweSample *zi,
                          const LweSample *w_i,
                          const LweSample *q_i0,
                          const LweSample *q_i1,
                          const TFheGateBootstrappingCloudKeySet *bk);


// =============================================================================
//
//  Extern variables
//

//  ----    TFHE params moved to separate file    ----


// =============================================================================
//
//  Function Implementations
//

// -----------------------------------------------------------------------------
//  LUT Bootstrapping: Identity, Threshold, Equality
//
void bs_id(LweSample *result,
           const LweSample *sample,
           const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << PI);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    bs_set_tv_identity(&testvect->coefsT[0], N, MU);

    // call BS priv
    bs_with_tv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}

void bs_gleq(LweSample *result,
             const LweSample *sample,
             const uint32_t thr,
             const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << PI);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    bs_set_tv_gleq(&testvect->coefsT[0], N, thr, MU);

    // call BS priv
    bs_with_tv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}

void bs_eq(LweSample *result,
           const LweSample *sample,
           const uint32_t thr,
           const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << PI);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    bs_set_tv_eq(&testvect->coefsT[0], N, thr, MU);

    // call BS priv
    bs_with_tv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}

// -----------------------------------------------------------------------------
//  En/Decryption
//
void paral_sym_encr_priv(LweSample *ct,
                         const int32_t message,
                         const TFheGateBootstrappingSecretKeySet *sk)
{
    Torus32 _1s16 = modSwitchToTorus32(1, 1 << PI);

    // scale to [-8/16, 7/16]
    //                         ___/ 8 \___     _/ mask 1111 \_     ___/ 8 \___
    Torus32 mu = (((message + (1 << (PI-1))) & ((1 << PI) - 1)) - (1 << (PI-1))) * _1s16;
    double alpha = sk->params->in_out_params->alpha_min; //TODO orig: specify noise
    lweSymEncrypt(ct, mu, alpha, sk->lwe_key);
}

void paral_sym_encr(LweSample *ct,
                    const int32_t message,
                    const TFheGateBootstrappingSecretKeySet *sk)
{
    if ((message < -2) || (2 < message))
        die_soon("Out of the alphabet A = [-2 .. 2].");

    paral_sym_encr_priv(ct, message, sk);
}

int32_t paral_sym_decr(const LweSample *sample,
                       const TFheGateBootstrappingSecretKeySet *sk)
{
    Torus32 _1s16 = modSwitchToTorus32(1, 1 << PI);
    Torus32 mu = lwePhase(sample, sk->lwe_key);
    return (mu + (_1s16 >> 1)) >> (32 - PI);
}

// -----------------------------------------------------------------------------
//  Parallel Addition
//
void parallel_add(LweSample *z,
                  const LweSample *x,
                  const LweSample *y,
                  const uint32_t wlen,
                  const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    // alloc aux arrays
    LweSample *w = new_LweSample_array(wlen, io_lwe_params);
    LweSample *q = new_LweSample_array(wlen, io_lwe_params);

    // calc w_i = x_i + y_i
    for (uint32_t i = 0; i < wlen; i++)
    {
        lweCopy (w + i, x + i, io_lwe_params);
        lweAddTo(w + i, y + i, io_lwe_params);
    }
            //~ // print w
            //~ printf(" W  | +0 ");
            //~ for (int32_t i = wlen - 1; i >= 0; i--)
            //~ {
                //~ plain = paral_sym_decr(w + i, sk);
                //~ printf("| %+d ", plain);
            //~ }
            //~ printf("|\n");

    // calc q_i's
    // i = 0
    paral_calc_qi(q, w, NULL, bk);
    // i > 0
    for (uint32_t i = 1; i < wlen; i++)
    {
        paral_calc_qi(q + i, w + i, w + i - 1, bk);
    }

    // calc z_i's
    // i = 0
    paral_calc_zi(z, w, q, NULL, bk);
    // i > 0
    for (uint32_t i = 1; i < wlen; i++)
    {
        paral_calc_zi(z + i, w + i, q + i, q + i - 1, bk);
    }
    // i = wlen
    paral_calc_zi(z + wlen, NULL, NULL, q + wlen - 1, bk);

    // progress bar end
    printf("\n");

    // cleanup
    delete_LweSample_array(wlen, q);
    delete_LweSample_array(wlen, w);
}

// -----------------------------------------------------------------------------
//  Misc
//
void die_soon(const char* message)
{
    fprintf(stderr, "(!) %s\n    Aborting ...\n", message);
    abort();
}


// =============================================================================
//
//  Static Function Implementations
//

static void bs_set_tv_identity(Torus32 *const tv,
                               const uint32_t N,
                               const Torus32 MU)
{
    uint32_t i, s;

    // identity for parallel addition needs to work on the interval [-2, 2] (can be extended to [-4,4])
    //          <-- to be set up manually ----->    <-- negacyclic extension ------>
    //   4    |                 ————             |                                   |
    //   3    |             ————    ————         |                                   |
    //   2    |         ————            ————     |                                   |
    //   1    |     ————                    ———— |    8   9  10  11  12  13  14  15  |
    //   0 ·· | ———— ··························  |  ———— ··························· | ····
    //  -1    |   0   1   2   3   4   5   6   7  |      ————                    ———— |
    //  -2    |                                  |          ————            ————     |
    //  -3    |                                  |              ————    ————         |
    //  -4    |                                  |                  ————             |
    //          <-- used ---------->                                <-- used ------>

    // make a 0-stair around zero (and 2^(PI-1))
    for (i = 0; i < (N >> PI); i++)     tv[i] = 0;
    for (i = N - (N >> PI); i < N; i++) tv[i] = 0;

    // make other stairs
    // from 1 to 2^(PI-2) - 1
    for (s = 1; s < (1 << (PI-2)); s++)
    {
        for (i = s * (N >> (PI-1)) - (N >> PI);  i < s * (N >> (PI-1)) + (N >> PI); i++)
            tv[i] = s * MU;
    }
    // from 2^(PI-2) to 2^(PI-1) - 1
    for (s = (1 << (PI-2)); s < (1 << (PI-1)); s++)
    {
        for (i = s * (N >> (PI-1)) - (N >> PI);  i < s * (N >> (PI-1)) + (N >> PI); i++)
            tv[i] = (-s + (1 << (PI-1))) * MU;
    }
}

static void bs_set_tv_gleq(Torus32 *const tv,
                           const uint32_t N,
                           const uint32_t thr,
                           const Torus32 MU)
{
    if ((thr > (1 << (PI - 2))) || thr == 0)
        die_soon("Threshold for bootstrapping too large or zero.");

    uint32_t i, s;

    // make a 0-stair around zero
    for (i = 0; i < (N >> PI); i++)     tv[i] = 0;
    for (i = N - (N >> PI); i < N; i++) tv[i] = 0;

    // make other stairs
    // 0's first part
    for (s = 1; s < thr; s++)
    {
        for (i = s * (N >> (PI-1)) - (N >> PI);  i < s * (N >> (PI-1)) + (N >> PI); i++)
            tv[i] = 0 * MU;
    }
    // 1's
    for (s = thr; s <= (1 << (PI-1)) - thr; s++)
    {
        for (i = s * (N >> (PI-1)) - (N >> PI);  i < s * (N >> (PI-1)) + (N >> PI); i++)
            tv[i] = 1 * MU;
    }
    // 0's second part
    for (s = (1 << (PI-1)) - thr + 1; s < (1 << (PI-1)); s++)
    {
        for (i = s * (N >> (PI-1)) - (N >> PI);  i < s * (N >> (PI-1)) + (N >> PI); i++)
            tv[i] = 0 * MU;
    }
}

static void bs_set_tv_eq(Torus32 *const tv,
                         const uint32_t N,
                         const uint32_t thr,
                         const Torus32 MU)
{
    if ((thr > (1 << (PI - 2))) || thr == 0)
        die_soon("Threshold for bootstrapping too large or zero.");

    uint32_t i, s;

    // make a 0-stair around zero
    for (i = 0; i < (N >> PI); i++)     tv[i] = 0;
    for (i = N - (N >> PI); i < N; i++) tv[i] = 0;

    // make other stairs
    // 0's elsewhere
    for (s = 1; s < (1 << (PI-1)); s++)
    {
        for (i = s * (N >> (PI-1)) - (N >> PI);  i < s * (N >> (PI-1)) + (N >> PI); i++)
            tv[i] = 0 * MU;
    }
    // 1's at specific positions
    s = thr;
    for (i = s * (N >> (PI-1)) - (N >> PI);  i < s * (N >> (PI-1)) + (N >> PI); i++)
        tv[i] = 1 * MU;
    s = (1 << (PI-1)) - thr;
    for (i = s * (N >> (PI-1)) - (N >> PI);  i < s * (N >> (PI-1)) + (N >> PI); i++)
        tv[i] = 1 * MU;
}

static void bs_with_tv(LweSample *result,
                       const LweSample *sample,
                       const TFheGateBootstrappingCloudKeySet *bk,
                       const TorusPolynomial *testvect)
{
    const LweBootstrappingKeyFFT *bkFFT = bk->bkFFT;
    LweSample *tmp = new_LweSample(&bkFFT->accum_params->extracted_lweparams);

    //  Bootstrapping BEGIN   --------------------------------------------------
    const TGswParams *bk_params = bkFFT->bk_params;
    const int32_t N = bkFFT->accum_params->N;
    const int32_t n = bkFFT->in_out_params->n;
    int32_t *bara = new int32_t[N];

    // Modulus switching
    int32_t barb = modSwitchFromTorus32(sample->b, 2*N);
    for (int32_t i = 0; i < n; i++)
        bara[i] = modSwitchFromTorus32(sample->a[i], 2*N);

    // Blind rotation and sample extraction
    tfhe_blindRotateAndExtract_FFT(tmp, testvect, bkFFT->bkFFT, barb, bara, n, bk_params);
    //FUCKUP #02:   TFheGateBootstrappingCloudKeySet has member 'bkFFT',
    //              which is LweBootstrappingKeyFFT, which also has member 'bkFFT', which is TGswSampleFFT

    delete[] bara;
    //  Bootstrapping END   ----------------------------------------------------

    // Key switching
    lweKeySwitch(result, bkFFT->ks, tmp);

    delete_LweSample(tmp);
}

static void paral_calc_qi(LweSample *qi,
                          const LweSample *w_i0,
                          const LweSample *w_i1,
                          const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    // aux variables
    LweSample *r1 = new_LweSample(io_lwe_params);
    LweSample *r2 = new_LweSample(io_lwe_params);
    LweSample *r3 = new_LweSample(io_lwe_params);
    LweSample *r23= new_LweSample(io_lwe_params);

    // progress bar ...
    printf("-");fflush(stdout);

    //            r1:  w_i <=> +-3
    bs_gleq(r1,  w_i0,     3,     bk);

    //            r2:  w_i == +-2
    bs_eq  (r2,  w_i0,    2,      bk);

    if (w_i1 == NULL)
    {
        lweNoiselessTrivial(r3, 0, io_lwe_params);
    }
    else
    {
        //            r3:  w_i-1 <=> +-2
        bs_gleq(r3,  w_i1,       2,   bk);
    }

    //            r23: r2+r3 == +-2
    lweAddTo(          r2,r3,           io_lwe_params);
    bs_eq  (r23, r2,        2,    bk);

    // q_i = r1 + r23
    lweCopy (qi, r1, io_lwe_params);
    lweAddTo(qi, r23, io_lwe_params);
}

static void paral_calc_zi(LweSample *zi,
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

    // progress bar ...
    printf("-");fflush(stdout);

    // bootstrap to refresh noise
    bs_id(zi, tmpz, bk);

    // cleanup
    delete_LweSample(tmpz);
}

int64_t paral_eval(const int32_t *const x,
                   const uint32_t len)
{
    int64_t r = 0;

    for (uint32_t i = 0; i < len; i++)
    {
        r += x[i] * (1 << (2*i));
    }

    return r;
}
