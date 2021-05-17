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
#define N_WITH_CARRY_SCENARIOS 3
#define N_PARALLEL_SCENARIOS 3
#define N_PARAM_SETS (1 + (N_WITH_CARRY_SCENARIOS) + (N_PARALLEL_SCENARIOS))


// =============================================================================
//
//  Types
//

typedef enum addition_scenario_t
{
    TFHE_LIB                =  0,
    A_CARRY_2_GATE_TFHE     =  1,
    B_CARRY_3_GATE_2_BIT    =  2,
    C_CARRY_4_BIT           =  3,
    D_PARALLEL_SC_1         =  4,
    E_PARALLEL_SC_2         =  5,
    F_PARALLEL_SC_3         =  6,
} addition_scenario_t;

/*******************************************************************************
 *  structure to hold TFHE params
 * */
typedef struct tfhe_params_t
{
    int32_t N;
    int32_t k;
    int32_t n;
    int32_t bk_l;
    int32_t bk_Bgbit;
    int32_t ks_basebit;
    int32_t ks_length;
    double ks_stdev;    // standard deviation of KSK
    double bk_stdev;    // standard deviation of BK
    double max_stdev;   // max standard deviation for a 1/4 msg space
} tfhe_params_t;


// =============================================================================
//
//  Extern Variables
//

/*******************************************************************************
 *  store to hold TFHE params' structures
 *  n.b., currently only for 'with KS' scenario !!
 * */
extern const tfhe_params_t tfhe_params_store[N_PARAM_SETS];


// =============================================================================
//
//  Function Prototypes
//

// -----------------------------------------------------------------------------
//  LUT Bootstrapping: Identity, Threshold, Equality
//

/**
 *  @brief          Description
 *
 */
void bs_id(LweSample *result,
           const LweSample *sample,
           const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void bs_gleq(LweSample *result,
             const LweSample *sample,
             const uint32_t thr,
             const TFheGateBootstrappingCloudKeySet *bk);

/**
 *  @brief          Description
 *
 */
void bs_eq(LweSample *result,
           const LweSample *sample,
           const uint32_t thr,
           const TFheGateBootstrappingCloudKeySet *bk);

// -----------------------------------------------------------------------------
//  En/Decryption
//
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

// -----------------------------------------------------------------------------
//  Parallel Addition
//
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
void parallel_add(const addition_scenario_t scenario,
                  LweSample *z,
                  const LweSample *x,
                  const LweSample *y,
                  const uint32_t wlen,
                  const TFheGateBootstrappingCloudKeySet *bk);

// -----------------------------------------------------------------------------
//  Misc
//
/**
 *  @brief          Lived fast -> die soon
 *
 *  @param[in]      Message to print
 *
 */
void die_soon(const char* message);

/**
 *  @brief          Description
 *
 */
int64_t paral_eval(const int32_t *const x,
                   const uint32_t len);

#endif // #ifndef PARALLEL_ADDITION_IMPL_H
