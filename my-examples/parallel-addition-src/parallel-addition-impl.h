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

#include <bootstrapping-functions.h>

#define TFHE_LIB              0
#define A_CARRY_5_BS_BIN      1
#define B_CARRY_2_BS_BIN      2
#define C_CARRY_2_BIT         3
#define D_PARALLEL_2          4     // D_PARALLEL_SC_1
#define E_PARALLEL_2          5     // E_PARALLEL_SC_2
#define F_PARALLEL_2          6     // F_PARALLEL_SC_3
#define G_PARALLEL_4          7
#define H_PARALLEL_4          8
#define I_PARALLEL_4          9

// choose sequential addition scenario (TFHE parameters are by default corresponding to this scenario)
#define SEQ_SCENARIO            B_CARRY_2_BS_BIN

// choose parallel addition scenario (TFHE parameters are by default corresponding to this scenario)
#define  PA_SCENARIO_BIN        F_PARALLEL_2

// choose parallel addition scenario (TFHE parameters are by default corresponding to this scenario)
#define  PA_SCENARIO_QUAD       H_PARALLEL_4

// choose TFHE parameters for bootstrapping tests
#define BS_TFHE_PARAMS_INDEX    F_PARALLEL_2

#define DBG_OUT

//  ----    do not edit    ----
#define SEQ_TFHE_PARAMS_INDEX SEQ_SCENARIO          // by default, use TFHE params derived for particular scenario
#define PAB_TFHE_PARAMS_INDEX  PA_SCENARIO_BIN      // by default, use TFHE params derived for particular scenario
#define PAQ_TFHE_PARAMS_INDEX  PA_SCENARIO_QUAD     // by default, use TFHE params derived for particular scenario
#define N_WITH_CARRY_SCENARIOS 3
#define N_PARALLEL_SCENARIOS 6
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

extern const uint32_t PI_S;
extern const uint32_t PI_B;
extern const uint32_t PI_Q;


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
//  En/Decryption
//

/**
 *  @brief          Description
 *
 */
void sym_encr_priv(LweSample *ct,
                   const int32_t message,
                   const uint32_t pi,
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
void paral_quad_sym_encr(LweSample *ct,
                         const int32_t message,
                         const TFheGateBootstrappingSecretKeySet *sk);

/**
 *  @brief          Description
 *
 */
void paral_bin_sym_encr(LweSample *ct,
                        const int32_t message,
                        const TFheGateBootstrappingSecretKeySet *sk);

/**
 *  @brief          Description
 *
 */
int32_t sym_decr(const LweSample *ct,
                 const uint32_t pi,
                 const TFheGateBootstrappingSecretKeySet *sk);

/**
 *  @brief          Description
 *
 */
int32_t bin_sym_decr(const LweSample *ct,
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
#ifdef DBG_OUT
                    const TFheGateBootstrappingSecretKeySet *sk,
#endif
                    const TFheGateBootstrappingCloudKeySet *bk);

// -----------------------------------------------------------------------------
//  Parallel Addition (Bin)
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
void parallel_add_bin(LweSample *z,
                      const LweSample *x,
                      const LweSample *y,
                      const uint32_t wlen_bin,
#ifdef DBG_OUT
                      const TFheGateBootstrappingSecretKeySet *sk,
#endif
                      const TFheGateBootstrappingCloudKeySet *bk);

// -----------------------------------------------------------------------------
//  Parallel Addition (Quad)
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
void parallel_add_quad(LweSample *z,
                       const LweSample *x,
                       const LweSample *y,
                       const uint32_t wlen_quad,
#ifdef DBG_OUT
                       const TFheGateBootstrappingSecretKeySet *sk,
#endif
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
int64_t quad_eval(const int32_t *const x,
                  const uint32_t len);

/**
 *  @brief          Description
 *
 */
int64_t bin_eval(const int32_t *const x,
                 const uint32_t len);

#endif // #ifndef PARALLEL_ADDITION_IMPL_H
