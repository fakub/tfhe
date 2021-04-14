#ifndef PARALLEL_ADDITION_IMPL_H
#define PARALLEL_ADDITION_IMPL_H

#include <cstdint>

#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"

#define PI 4

/**
 *  @brief          Lived fast -> die soon
 *
 *  @param[in]      Message to print
 *
 */
void die_soon(const char* message);

/**
 *  @brief          Main Parallel Addition Function
 *
 *  @param[out]     LWE Sample (length +wlen+ + 1)
 *  @param[in]      LWE Sample (length +wlen+)
 *  @param[in]      LWE Sample (length +wlen+)
 *  @param[in]      Length of LWE samples
 *  @param[in]      Bootstrapping Keys
 *
 */
void parallel_add(LweSample *z,
                  const LweSample *x,
                  const LweSample *y,
                  const uint32_t wlen,
                  const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void paral_sym_encr_priv(LweSample *ct,
                         const int32_t message,
                         const TFheGateBootstrappingSecretKeySet *sk);

/**
 *  @brief          Description
 *
 */
void paral_sym_encr(LweSample *ct,
                    const int32_t message,
                    const TFheGateBootstrappingSecretKeySet *sk);

/**
 *  @brief          Description
 *
 */
int32_t paral_sym_decr(const LweSample *sample,
                       const TFheGateBootstrappingSecretKeySet *sk);

/**
 *  @brief          Description
 *
 */
void paral_bs_id(LweSample *result,
                 const LweSample *sample,
                 const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void paral_bs_gleq(LweSample *result,
                   const LweSample *sample,
                   const uint32_t thr,
                   const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void paral_bs_eq(LweSample *result,
                 const LweSample *sample,
                 const uint32_t thr,
                 const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void paral_add(LweSample *z,
               const LweSample *x,
               const LweSample *y,
               const TFheGateBootstrappingCloudKeySet *bk);

#endif // #ifndef PARALLEL_ADDITION_IMPL_H
