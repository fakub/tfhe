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

#define TFHE_LIB              0
#define A_CARRY_2_GATE_TFHE   1
#define B_CARRY_3_GATE_2_BIT  2
#define C_CARRY_4_BIT         3
#define D_PARALLEL_SC_1       4
#define E_PARALLEL_SC_2       5
#define F_PARALLEL_SC_3       6

// choose parallel addition scenario (TFHE parameters are by default corresponding to this scenario)
#define  PA_SCENARIO    D_PARALLEL_SC_1

// choose sequential addition scenario (TFHE parameters are by default corresponding to this scenario)
#define SEQ_SCENARIO    A_CARRY_2_GATE_TFHE

// choose TFHE parameters for bootstrapping tests
#define BS_TFHE_PARAMS_INDEX    D_PARALLEL_SC_1

//  ----    do not edit    ----
#define  PA_TFHE_PARAMS_INDEX  PA_SCENARIO   // by default, use TFHE params derived for particular scenario
#define SEQ_TFHE_PARAMS_INDEX SEQ_SCENARIO   // by default, use TFHE params derived for particular scenario
#define N_WITH_CARRY_SCENARIOS 3
#define N_PARALLEL_SCENARIOS 3
#define N_PARAM_SETS (1 + (N_WITH_CARRY_SCENARIOS) + (N_PARALLEL_SCENARIOS))

#define XSTR(x) STR(x)
#define STR(x) #x


// =============================================================================
//
//  Types
//

/*******************************************************************************
 *  structure to hold TFHE params
 * */
typedef struct tfhe_params_t
{
    uint32_t pi;
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

extern const uint32_t PI;


// =============================================================================
//
//  Function Prototypes
//

// -----------------------------------------------------------------------------
//  TFHE Param Setup
//

/**
 *  @brief          Description
 *
 */
void setup_TFHE_params(const int tfhe_params_index,
                       TFheGateBootstrappingParameterSet **const tfhe_params);

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
void sym_encr_priv(LweSample *ct,
                   const int32_t message,
                   const TFheGateBootstrappingSecretKeySet *sk);

/**
 *  @brief          Description
 *
 */
void bin_sym_encr(LweSample *ct,
                  const int32_t message,
                  const TFheGateBootstrappingSecretKeySet *sk);

/**
 *  @brief          Description
 *
 */
void seq_quad_sym_encr(LweSample *ct,
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
int32_t sym_decr(const LweSample *sample,
                 const TFheGateBootstrappingSecretKeySet *sk);

/**
 *  @brief          Description
 *
 */
bool bin_sym_decr(LweSample *ct,
                  const TFheGateBootstrappingSecretKeySet *sk);

/**
 *  @brief          Description
 *
 */
void bin_noiseless(LweSample *ct,
                   const int32_t message,
                   const TFheGateBootstrappingCloudKeySet *bk);

// -----------------------------------------------------------------------------
//  Sequential Addition
//

/**
 *  @brief          Main Sequential Addition Function
 *
 *  @param[out]     LWE Sample (length +wlen+ + 1)
 *  @param[in]      LWE Sample (length +wlen+)
 *  @param[in]      LWE Sample (length +wlen+)
 *  @param[in]      Length of LWE samples
 *  @param[in]      Bootstrapping Keys
 *
 */

void sequential_add(LweSample *z,
                    const LweSample *x,
                    const LweSample *y,
                    const uint32_t wlen,
                    const TFheGateBootstrappingCloudKeySet *bk);

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
void parallel_add(LweSample *z,
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

/**
 *  @brief          Description
 *
 */
int64_t seq_eval(const int32_t *const x,
                 const uint32_t len);

#endif // #ifndef PARALLEL_ADDITION_IMPL_H
