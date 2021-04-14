#include <cstdio>

#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"

#include <parallel-addition-impl.h>


// Static Function Prototypes   ------------------------------------------------

/**
 *  @brief          Description
 *
 */
static void paral_bs_set_tv_identity(Torus32 *const tv,
                                     const uint32_t N,
                                     const Torus32 MU);

/**
 *  @brief          Description
 *
 */
static void paral_bs_set_tv_gleq(Torus32 *const tv,
                                 const uint32_t N,
                                 const uint32_t thr,
                                 const Torus32 MU);

/**
 *  @brief          Description
 *
 */
static void paral_bs_set_tv_eq(Torus32 *const tv,
                               const uint32_t N,
                               const uint32_t thr,
                               const Torus32 MU);

/**
 *  @brief          Description
 *
 */
static void paral_bs_priv(LweSample *result,
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


// Function Implementations   --------------------------------------------------

void die_soon(const char* message)
{
    fprintf(stderr, "(!) %s\n    Aborting ...\n", message);
    abort();
}

void parallel_add(LweSample *z,
                  const LweSample *x,
                  const LweSample *y,
                  const uint32_t wlen,
                  const TFheGateBootstrappingCloudKeySet *bk)
{
    //TODO corner cases where there are supposed to be zeros
    for (uint32_t i = 2; i < wlen; i++)
    {
        paral_add(z + i, x + i, y + i, bk);
        printf("-");fflush(stdout);
    }
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
    double alpha = sk->params->in_out_params->alpha_min; //TODO: specify noise
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
//  LUT Bootstrapping: Threshold
//
static void paral_bs_set_tv_identity(Torus32 *const tv,
                                     const uint32_t N,
                                     const Torus32 MU)
{
    uint32_t i, s;

    // make a 0-stair around zero
    for (i = 0; i < (N >> PI); i++)     tv[i] = 0;
    for (i = N - (N >> PI); i < N; i++) tv[i] = 0;

    // make other stairs
    for (s = 1; s < (1 << (PI-1)); s++)   // due to negacyclicity, only half values are to be set
    {
        for (i = s * (N >> (PI-1)) - (N >> PI);  i < s * (N >> (PI-1)) + (N >> PI); i++)
            tv[i] = s * MU;
    }
}

static void paral_bs_set_tv_gleq(Torus32 *const tv,
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

static void paral_bs_set_tv_eq(Torus32 *const tv,
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

static void paral_bs_priv(LweSample *result,
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

void paral_bs_id(LweSample *result,
                 const LweSample *sample,
                 const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << PI);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    paral_bs_set_tv_identity(&testvect->coefsT[0], N, MU);

    // call BS priv
    paral_bs_priv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}

void paral_bs_gleq(LweSample *result,
                   const LweSample *sample,
                   const uint32_t thr,
                   const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << PI);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    paral_bs_set_tv_gleq(&testvect->coefsT[0], N, thr, MU);

    // call BS priv
    paral_bs_priv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}

void paral_bs_eq(LweSample *result,
                 const LweSample *sample,
                 const uint32_t thr,
                 const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << PI);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    paral_bs_set_tv_eq(&testvect->coefsT[0], N, thr, MU);

    // call BS priv
    paral_bs_priv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}


// -----------------------------------------------------------------------------
//  Helper Functions
//
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

    //            r1:  w_i <=> +-3
    paral_bs_gleq(r1,  w_i0,     3,     bk);

    //            r2:  w_i == +-2
    paral_bs_eq  (r2,  w_i0,    2,      bk);

    //            r3:  w_i-1 <=> +-2
    paral_bs_gleq(r3,  w_i1,       2,   bk);

    //            r23: r2+r3 == +-2
    lweAddTo(          r2,r3,           io_lwe_params);
    paral_bs_eq  (r23, r2,        2,    bk);

    // q_i = r1 + r23
    lweNoiselessTrivial(qi, 0, io_lwe_params);
    lweAddTo(qi, r1, io_lwe_params);
    lweAddTo(qi, r23, io_lwe_params);
}


// -----------------------------------------------------------------------------
//  Parallel Addition
//
// !! there must be two more samples to the right in x and y !!
void paral_add(LweSample *z,
               const LweSample *x,
               const LweSample *y,
               const TFheGateBootstrappingCloudKeySet *bk)
{
    const LweParams *io_lwe_params = bk->params->in_out_params;

    // alloc aux variable arrays for w_i, q_i
    LweSample *w = new_LweSample_array(3, io_lwe_params);
    LweSample *q = new_LweSample_array(2, io_lwe_params);
    LweSample *tmpz = new_LweSample(io_lwe_params);
                LweSample *tmpzz = new_LweSample(io_lwe_params);   // unless parameters allow addition of 12 values

    // calc w_i = x_i + y_i   for i, i-1, i-2
    for (int i = 0; i < 3; i++)
    {
        lweNoiselessTrivial(w + i, 0, io_lwe_params);
        lweAddTo(w + i, x - i, io_lwe_params);   // n.b., w_i-1 is at w + 1 position
        lweAddTo(w + i, y - i, io_lwe_params);
    }

    // calc q_i and q_i-1
    paral_calc_qi(q + 0, w + 0, w + 1, bk);
    paral_calc_qi(q + 1, w + 1, w + 2, bk);

    // calculate result: z_i = w_i - 4q_i + q_i-1
    lweNoiselessTrivial(tmpz, 0, io_lwe_params);
    lweAddTo(tmpz, w, io_lwe_params);
    //~ lweSubMulTo(tmpz, 4, q, io_lwe_params);   //TODO params to allow 12 additions
                // addition of 12 values is already too much for these parameters
                lweSubMulTo(tmpz, 2, q, io_lwe_params);
                paral_bs_id(tmpzz, tmpz, bk); lweCopy(tmpz, tmpzz, io_lwe_params);   // noise must be refreshed
                lweSubMulTo(tmpz, 2, q, io_lwe_params);

    lweAddTo(tmpz, q + 1, io_lwe_params);

    // bootstrap to refresh noise
    paral_bs_id(z, tmpz, bk);

    // cleanup
                delete_LweSample(tmpzz);   // unless parameters allow addition of 12 values
    delete_LweSample(tmpz);
    delete_LweSample_array(2, q);
    delete_LweSample_array(3, w);
}
