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

/**
 *  @brief          Identity around zero: [-2^pi / 4, 2^pi / 4]
 *
 */
static void bs_set_tv_identity(Torus32 *const tv,
                               const uint32_t N,
                               const uint32_t pi,
                               const Torus32 MU);

/**
 *  @brief          Identity at positive half: [0, 2^pi / 2)
 *
 */
static void bs_set_tv_mod(Torus32 *const tv,
                          const uint32_t N,
                          const uint32_t pi,
                          const uint32_t mod,
                          const Torus32 MU);

/**
 *  @brief          Description
 *
 */
static void bs_set_tv_gleq(Torus32 *const tv,
                           const uint32_t N,
                           const uint32_t pi,
                           const uint32_t thr,
                           const Torus32 MU);

/**
 *  @brief          Description
 *
 */
static void bs_set_tv_pos_gleq(Torus32 *const tv,
                               const uint32_t N,
                               const uint32_t pi,
                               const uint32_t thr,
                               const Torus32 MU);

/**
 *  @brief          Description
 *
 */
static void bs_set_tv_eq(Torus32 *const tv,
                         const uint32_t N,
                         const uint32_t pi,
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


// =============================================================================
//
//  Function Implementations
//

// -----------------------------------------------------------------------------
//  LUT Bootstrapping: Identity, Threshold, Equality
//
void bs_id(LweSample *result,
           const LweSample *sample,
           const uint32_t pi,
           const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << pi);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    bs_set_tv_identity(&testvect->coefsT[0], N, pi, MU);

    // call BS priv
    bs_with_tv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}

void bs_mod(LweSample *result,
            const LweSample *sample,
            const uint32_t mod,
            const uint32_t pi,
            const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << pi);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    bs_set_tv_mod(&testvect->coefsT[0], N, pi, mod, MU);

    // call BS priv
    bs_with_tv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}

void bs_gleq(LweSample *result,
             const LweSample *sample,
             const uint32_t thr,
             const uint32_t pi,
             const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << pi);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    bs_set_tv_gleq(&testvect->coefsT[0], N, pi, thr, MU);

    // call BS priv
    bs_with_tv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}

void bs_pos_gleq(LweSample *result,
                 const LweSample *sample,
                 const uint32_t thr,
                 const uint32_t pi,
                 const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << pi);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    bs_set_tv_pos_gleq(&testvect->coefsT[0], N, pi, thr, MU);

    // call BS priv
    bs_with_tv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}

void bs_eq(LweSample *result,
           const LweSample *sample,
           const uint32_t thr,
           const uint32_t pi,
           const TFheGateBootstrappingCloudKeySet *bk)
{
    // get N, MU
    const uint32_t N = (uint32_t)(bk->bkFFT->accum_params->N);
    const Torus32 MU = modSwitchToTorus32(1, 1 << pi);

    // init test vector with stair-case function
    TorusPolynomial *testvect = new_TorusPolynomial(N);
    bs_set_tv_eq(&testvect->coefsT[0], N, pi, thr, MU);

    // call BS priv
    bs_with_tv(result, sample, bk, testvect);

    // delete test vector
    delete_TorusPolynomial(testvect);
}


// =============================================================================
//
//  Static Function Implementations
//

static void bs_set_tv_identity(Torus32 *const tv,
                               const uint32_t N,
                               const uint32_t pi,
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

    // make a 0-stair around zero (and 2^(pi-1))
    for (i = 0; i < (N >> pi); i++)     tv[i] = 0;
    for (i = N - (N >> pi); i < N; i++) tv[i] = 0;

    // make other stairs
    // from 1 to 2^(pi-2) - 1
    for (s = 1u; s < (1u << (pi-2)); s++)
    {
        for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
            tv[i] = s * MU;
    }
    // from 2^(pi-2) to 2^(pi-1) - 1
    for (s = (1u << (pi-2)); s < (1u << (pi-1)); s++)
    {
        for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
            tv[i] = (-s + (1 << (pi-1))) * MU;
    }
}

static void bs_set_tv_mod(Torus32 *const tv,
                          const uint32_t N,
                          const uint32_t pi,
                          const uint32_t mod,
                          const Torus32 MU)
{
    if ((mod > (1u << (pi - 1))) || mod == 0)
    {
        char buff[100];
        sprintf(&buff[0], "Modulo for mod bootstrapping too large or zero: mod = %d > %d = 2^(pi-1).", mod, (1 << (pi - 1)));
        die_soon(&buff[0]);
    }

    uint32_t i, s;

    // make a 0-stair around zero (and 2^(pi-1))
    for (i = 0; i < (N >> pi); i++)     tv[i] = 0;
    for (i = N - (N >> pi); i < N; i++) tv[i] = 0;

    // make other stairs
    // from 1 to 2^(pi-1) - 1
    for (s = 1u; s < (1u << (pi-1)); s++)
    {
        for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
            tv[i] = (s % mod) * MU;
    }
}

static void bs_set_tv_gleq(Torus32 *const tv,
                           const uint32_t N,
                           const uint32_t pi,
                           const uint32_t thr,
                           const Torus32 MU)
{
    if ((thr > (1u << (pi - 2))) || thr == 0)
    {
        char buff[100];
        sprintf(&buff[0], "Threshold for around-zero bootstrapping too large or zero: thr = %d > %d = 2^(pi-2).", thr, (1 << (pi - 2)));
        die_soon(&buff[0]);
    }

    uint32_t i, s;

    // make a 0-stair around zero
    for (i = 0; i < (N >> pi); i++)     tv[i] = 0;
    for (i = N - (N >> pi); i < N; i++) tv[i] = 0;

    // make other stairs
    // 0's first part
    for (s = 1; s < thr; s++)
    {
        for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
            tv[i] = 0 * MU;
    }
    // 1's
    for (s = thr; s <= (1 << (pi-1)) - thr; s++)
    {
        for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
            tv[i] = 1 * MU;
    }
    // 0's second part
    for (s = (1u << (pi-1)) - thr + 1; s < (1u << (pi-1)); s++)
    {
        for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
            tv[i] = 0 * MU;
    }
}

static void bs_set_tv_pos_gleq(Torus32 *const tv,
                               const uint32_t N,
                               const uint32_t pi,
                               const uint32_t thr,
                               const Torus32 MU)
{
    if ((thr >= (1u << (pi - 1))) || thr == 0)
    {
        char buff[100];
        sprintf(&buff[0], "Threshold for positive bootstrapping too large or zero: thr = %d >= %d = 2^(pi-1).", thr, (1 << (pi - 1)));
        die_soon(&buff[0]);
    }

    uint32_t i, s;

    // make a 0-stair around zero
    for (i = 0; i < (N >> pi); i++)     tv[i] = 0;
    for (i = N - (N >> pi); i < N; i++) tv[i] = 0;

    // make other stairs
    // 0's first part
    for (s = 1; s < thr; s++)
    {
        for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
            tv[i] = 0 * MU;
    }
    // 1's
    for (s = thr; s < (1u << (pi-1)); s++)
    {
        for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
            tv[i] = 1 * MU;
    }
}

static void bs_set_tv_eq(Torus32 *const tv,
                         const uint32_t N,
                         const uint32_t pi,
                         const uint32_t thr,
                         const Torus32 MU)
{
    if ((thr > (1u << (pi - 2))) || thr == 0)
    {
        char buff[100];
        sprintf(&buff[0], "Threshold for bootstrapping too large or zero: thr = %d > %d = 2^(pi-2).", thr, (1 << (pi - 2)));
        die_soon(&buff[0]);
    }

    uint32_t i, s;

    // make a 0-stair around zero
    for (i = 0; i < (N >> pi); i++)     tv[i] = 0;
    for (i = N - (N >> pi); i < N; i++) tv[i] = 0;

    // make other stairs
    // 0's elsewhere
    for (s = 1; s < (1u << (pi-1)); s++)
    {
        for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
            tv[i] = 0 * MU;
    }
    // 1's at specific positions
    s = thr;
    for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
        tv[i] = 1 * MU;
    s = (1 << (pi-1)) - thr;
    for (i = s * (N >> (pi-1)) - (N >> pi);  i < s * (N >> (pi-1)) + (N >> pi); i++)
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
