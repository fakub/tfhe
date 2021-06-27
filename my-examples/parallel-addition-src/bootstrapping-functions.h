#ifndef BS_FUNCTIONS_H
#define BS_FUNCTIONS_H

#include <cstdint>

#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"

#include <parallel-addition-impl.h>

// -----------------------------------------------------------------------------
//  LUT Bootstrapping: Identity, Threshold, Equality
//

/**
 *  @brief          Identity around zero
 *
 */
void bs_id(LweSample *result,
           const LweSample *sample,
           const uint32_t pi,
           const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Identity around zero
 *
 */
void bs_pos_id(LweSample *result,
               const LweSample *sample,
               const uint32_t pi,
               const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Identity at positive half
 *
 */
void bs_mod(LweSample *result,
            const LweSample *sample,
            const uint32_t mod,
            const uint32_t pi,
            const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void bs_gleq(LweSample *result,
             const LweSample *sample,
             const uint32_t thr,
             const uint32_t mult,
             const uint32_t pi,
             const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void bs_pos_gleq(LweSample *result,
                 const LweSample *sample,
                 const uint32_t thr,
                 const uint32_t pi,
                 const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void bs_eq(LweSample *result,
           const LweSample *sample,
           const uint32_t thr,
           const uint32_t pi,
           const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void bs_01(LweSample *result,
           const LweSample *sample,
           const uint32_t pi,
           const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void bs_max_bin(LweSample *result,
                const LweSample *sample,
                const uint32_t pi,
                const TFheGateBootstrappingCloudKeySet *bk);

#endif // #ifndef BS_FUNCTIONS_H
