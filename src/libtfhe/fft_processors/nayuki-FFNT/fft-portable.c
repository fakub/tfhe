/*
 * Fast Fourier transform (C)
 *
 * Copyright (c) 2016 Project Nayuki
 * https://www.nayuki.io/page/fast-fourier-transform-in-x86-assembly
 *
 * Extended by https://eprint.iacr.org/2021/480
 *
 * (MIT License)
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "fft.h"



// ===   Data types   ==========================================================

/*******************************************************************************
 * Type of tables.
 * */
typedef enum {
    FFT_e   = 0,
    IFFT_e  = 1,
} transform_t;

/*******************************************************************************
 * Structure that holds FFT tables.
 * */
typedef struct FftTables
{
    transform_t type;
    size_t n;
    size_t cos_denom;
    size_t tgt_size;
    size_t * bit_reversed;
    double * cos_table;
    double * sin_table;
    double * ct_path_tables;
    double * gs_path_tables;
} FftTables;



// ===   Private function prototypes   =========================================

/*******************************************************************************
 *  Implements FFT tables initialization.
 */
static FftTables * omegas_init(const size_t n,
                               const bool inverse);

/*******************************************************************************
 *  Implements FFT via Cooley-Tukey data path.
 */
static void dft_CT(const FftTables *const tables,
                   double *const real,
                   double *const imag,
                   const bool inverse,
                   const bool bitrev);

/*******************************************************************************
 *  Implements FFT via Gentleman-Sande data path.
 */
static void dft_GS(const FftTables *const tables,
                   double *const real,
                   double *const imag,
                   const bool inverse,
                   const bool bitrev);

/*******************************************************************************
 *  Inline function to permute +re+ and +im+ array, given a permutation.
 */
static inline void permute_ary(const size_t *const permutation,
                               double *const re,
                               double *const im,
                               const size_t n,
                               const bool inverse);
static double accurate_sine(uint64_t i, uint64_t n);
static int32_t floor_log2(size_t n);
static size_t reverse_bits(const size_t x, const uint32_t n);



// ===   Function implementations   ============================================

// ---   Init tables   ---------------------------------------------------------

FftTables * fft_init(const size_t n)
{
    return omegas_init(n, false);
}
FftTables * ifft_init(const size_t n)
{
    return omegas_init(n, true);
}

static FftTables * omegas_init(const size_t n,
                               const bool inverse)
{
    //TODO check these checks
    // check size
    if (n < 4 || (n & (n - 1)) != 0)
        return NULL;  // too small or not a power of 2
    if ((n / 2 + 2 * n - 8) > SIZE_MAX / sizeof(double) || n > SIZE_MAX / sizeof(size_t))
        return NULL;  // too large

    // alloc table structure
    FftTables * tables = (FftTables *)malloc(sizeof(FftTables));
    if (tables == NULL)
        return NULL;

    tables->n = n;
    tables->type = inverse ? IFFT_e : FFT_e;
    tables->cos_denom = n;
    tables->tgt_size = 2 * (tables->n - 4);

    // alloc arrays
    tables->bit_reversed = (size_t *)malloc(tables->n * sizeof(size_t));

    // cos/sin tables (needed for twisting in FFNT)
    tables->cos_table   = (double *)malloc(tables->cos_denom / 2 * sizeof(double));
    tables->sin_table   = (double *)malloc(tables->cos_denom / 2 * sizeof(double));
    tables->ct_path_tables = (double *)malloc(tables->tgt_size      * sizeof(double));
    tables->gs_path_tables = (double *)malloc(tables->tgt_size      * sizeof(double));

    // check allocation
    if (tables->bit_reversed    == NULL ||
        tables->cos_table       == NULL ||
        tables->sin_table       == NULL ||
        tables->ct_path_tables  == NULL ||
        tables->gs_path_tables  == NULL)
    {
        free(tables->bit_reversed);
        free(tables->cos_table);
        free(tables->sin_table);
        free(tables->ct_path_tables);
        free(tables->gs_path_tables);
        free(tables);
        return NULL;
    }

    // bit-reverse permutation
    size_t i;
    int32_t levels = floor_log2(tables->n);
    for (i = 0; i < tables->n; i++)
        tables->bit_reversed[i] = reverse_bits(i, levels);

    // cos/sin tables
    double angle;
    for (i = 0; i < tables->cos_denom / 2; i++)
    {
        angle = (inverse ? -1 : 1) * 2 * M_PI * i / (tables->cos_denom);
        tables->cos_table[i] = cos(angle);
        tables->sin_table[i] = sin(angle);
    }

    size_t size;
    size_t j, k;

    // trigonometric tables for each FFT internal level of the CT data path
    k = 0;
    for (size = 8; size <= tables->n; size <<= 1)
    {
        for (i = 0; i < size / 2; i += 4)
        {
            for (j = 0; j < 4; j++, k++)
            {
                tables->ct_path_tables[k] = accurate_sine(i + j + size / 4, size);  // Cosine
            }

            for (j = 0; j < 4; j++, k++)
            {
                tables->ct_path_tables[k] = (inverse ? -1 : 1) * accurate_sine(i + j, size);  // Sine
            }
        }
        // prevent overflow
        if (size == tables->n)
            break;
    }

    // trigonometric tables for each FFT internal level of the GS data path
    k = 0;
    for (size = tables->n; size >= 8; size >>= 1)
    {
        for (i = 0; i < size / 2; i += 4)
        {
            for (j = 0; j < 4; j++, k++)
            {
                tables->gs_path_tables[k] = accurate_sine(i + j + size / 4, size);  // Cosine, TODO
            }

            for (j = 0; j < 4; j++, k++)
            {
                tables->gs_path_tables[k] = (inverse ? -1 : 1) * accurate_sine(i + j, size);  // Sine, TODO
            }
        }
    }

    return tables;
}


// ---   FFT   -----------------------------------------------------------------

// in general, use the CT data path
void  fft_transform(const FftTables *const fft_tables,
                    double *const real,
                    double *const imag)
{
    fft_transform_CT(fft_tables, real, imag, true);
}
void ifft_transform(const FftTables *const ifft_tables,
                    double *const real,
                    double *const imag)
{
    ifft_transform_CT(ifft_tables, real, imag, true);
}

// call specific data path
void  fft_transform_CT(const FftTables *const fft_tables,
                       double *const real,
                       double *const imag,
                       const bool bitrev)
{
    dft_CT(fft_tables, real, imag, false, bitrev);
}
void ifft_transform_CT(const FftTables *const fft_tables,
                       double *const real,
                       double *const imag,
                       const bool bitrev)
{
    dft_CT(fft_tables, real, imag, true, bitrev);
}
void  fft_transform_GS(const FftTables *const fft_tables,
                       double *const real,
                       double *const imag,
                       const bool bitrev)
{
    dft_GS(fft_tables, real, imag, false, bitrev);
}
void ifft_transform_GS(const FftTables *const fft_tables,
                       double *const real,
                       double *const imag,
                       const bool bitrev)
{
    dft_GS(fft_tables, real, imag, true, bitrev);
}


// ---   FFNT   ----------------------------------------------------------------

void ffnt_transform(const FftTables *const ffnt_tables_2N,
                    const FftTables *const ffnt_tables_N_2,
                    double *const real,
                    double *const imag)
{
    const double *const fan_cos = ffnt_tables_2N->cos_table;
    const double *const fan_sin = ffnt_tables_2N->sin_table;
    size_t n_2 = ffnt_tables_N_2->n;

    // fan-style multiplication
    for (size_t i = 0; i < n_2; i++)
    {
        imag[i] = real[i] * fan_sin[i] + real[n_2+i] * fan_cos[i];
        real[i] = real[i] * fan_cos[i] - real[n_2+i] * fan_sin[i];
    }

    fft_transform_GS(ffnt_tables_N_2, real, imag, false);
}

void iffnt_transform(const FftTables *const ffnt_tables_2N,
                     const FftTables *const iffnt_tables_N_2,
                     double *const real,
                     double *const imag)
{
    const double *const fan_cos = ffnt_tables_2N->cos_table;
    const double *const fan_sin = ffnt_tables_2N->sin_table;
    size_t n_2 = iffnt_tables_N_2->n;

    ifft_transform_CT(iffnt_tables_N_2, real, imag, false);

    // fan-style multiplication (negative index)
    for (size_t i = 0; i < n_2; i++)
    {
        real[n_2+i] = -real[i] * fan_sin[i] + imag[i] * fan_cos[i];
        real[i]     =  real[i] * fan_cos[i] + imag[i] * fan_sin[i];
    }
}


// ---   Implementation of FFT via the Cooley-Tukey data path   ----------------

static void dft_CT(const FftTables *const tables,
                   double *const real,
                   double *const imag,
                   const bool inverse,
                   const bool bitrev)
{
    // check table type
    if (!((inverse && tables->type == IFFT_e) || (!inverse && tables->type == FFT_e))) return;

    const size_t n = tables->n;

    size_t i, j, k;
    double re, im, tpre, tpim;

    size_t halfsize, tablestep, size = 2;

    // Cooley-Tukey decimation-in-time
    const double * trigtables = tables->ct_path_tables;   // no *const since the pointer is iterated

    if (bitrev)
    {
#ifdef VERBOSE
    printf("\nBefore bit-rev:\n");
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

        const size_t *const bitreversed = tables->bit_reversed;
        permute_ary(bitreversed, real,
                    imag,
                    n, inverse);

#ifdef VERBOSE
    printf("\nAfter bit-rev:\n");
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE
    }
    else if (inverse)   // n.b., division is a part of bitreverse
    {
        for (i = 0; i < n; i++)
        {
            real[i] /= n;
            imag[i] /= n;
        }
    }

#ifdef VERBOSE
        printf("\nCT Round #%d: {\n", floor_log2(2));
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

    // size 2 merge (special)
    if (n >= 2)
    {
        for (i = 0; i < n; i += 2)
        {
            // addition and subtraction
            tpre = real[i];
            tpim = imag[i];

            real[i] += real[i + 1];
            imag[i] += imag[i + 1];
            real[i + 1] = tpre - real[i + 1];
            imag[i + 1] = tpim - imag[i + 1];
        }
    }

#ifdef VERBOSE
        printf("\nCT Round #%d: {\n", floor_log2(4));
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

    // size 4 merge (special)
    if (n >= 4)
    {
        if (!inverse)
        {
            for (i = 0; i < n; i += 4)
            {
                // even indices
                tpre = real[i + 2];   // y
                tpim = imag[i + 2];
                real[i + 2] = real[i] - tpre;   // y = x - y
                imag[i + 2] = imag[i] - tpim;
                real[i] += tpre;   // x = x + y
                imag[i] += tpim;

                // odd indices (times -i)
                tpre =  imag[i + 3];   // -iy
                tpim = -real[i + 3];
                real[i + 3] = real[i + 1] - tpre;   // y = x - iy
                imag[i + 3] = imag[i + 1] - tpim;
                real[i + 1] += tpre;   // x = x + iy
                imag[i + 1] += tpim;
            }
        }
        else
        {
            for (i = 0; i < n; i += 4)
            {
                // even indices
                tpre = real[i + 2];
                tpim = imag[i + 2];

                real[i + 2] = real[i] - tpre;
                imag[i + 2] = imag[i] - tpim;
                real[i] += tpre;
                imag[i] += tpim;

                // odd indices (times i)
                tpre = -imag[i + 3];
                tpim =  real[i + 3];

                real[i + 3] = real[i + 1] - tpre;
                imag[i + 3] = imag[i + 1] - tpim;
                real[i + 1] += tpre;
                imag[i + 1] += tpim;
            }
        }
    }

    for (size = 8; size <= n; size <<= 1)
    {
        halfsize = size >> 1;

#ifdef VERBOSE
        printf("\nCT Round #%d: {\n", floor_log2(size));
        for (int ii = 0; ii < n/2; ii++)
            printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

        for (i = 0; i < n; i += size)
        {
            for (j = 0, tablestep = 0;  j < halfsize;   j += 4, tablestep += 8)
            {
                for (k = 0; k < 4; k++)   // To simulate x86 AVX 4-vectors
                {
                    uint64_t vi = i + j + k;        // Vector index
                    uint64_t ti = tablestep + k;    // Table index

                    // in CT, calc
                    //  x = x + wy
                    //  y = x - wy

                    re = real[vi + halfsize];   // y
                    im = imag[vi + halfsize];
                    tpre = re * trigtables[ti] + im * trigtables[ti + 4];   // wy
                    tpim = im * trigtables[ti] - re * trigtables[ti + 4];

                    real[vi + halfsize] = real[vi] - tpre;   // y = x - wy
                    imag[vi + halfsize] = imag[vi] - tpim;
                    real[vi] += tpre;   // x = x + wy
                    imag[vi] += tpim;
                }
            }
        }
        if (size == n)   // Prevent overflow in 'size *= 2'
            break;

        trigtables += size;
    }

#ifdef VERBOSE
    printf("\nCT Result:\n");
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

}


// ---   Implementation of FFT via the Gentleman-Sande data path   -------------

static void dft_GS(const FftTables *const tables,
                   double *const real,
                   double *const imag,
                   const bool inverse,
                   const bool bitrev)
{
    // check table type
    if (!((inverse && tables->type == IFFT_e) || (!inverse && tables->type == FFT_e))) return;

    const size_t n = tables->n;

    size_t i, j, k;
    double re, im;

    size_t halfsize, tablestep, size;

    // Gentleman-Sande decimation-in-frequency
    const double * trigtables = tables->gs_path_tables;   // no *const since the pointer is iterated

    for (size = n; size >= 8; size >>= 1)
    {
        halfsize = size >> 1;

#ifdef VERBOSE
        printf("\nGS Round #%d: {\n", floor_log2(size));
        for (int ii = 0; ii < n/2; ii++)
            printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

        for (i = 0; i < n; i += size)
        {
            for (j = 0, tablestep = 0;  j < halfsize;   j += 4, tablestep += 8)
            {
                for (k = 0; k < 4; k++)   // To simulate x86 AVX 4-vectors
                {
                    uint64_t vi = i + j + k;        // Vector index
                    uint64_t ti = tablestep + k;    // Table index

                    // in GS, calc
                    //  x = x + y
                    //  y = w(x - y)
                    re = real[vi] - real[vi + halfsize];   // x - y
                    im = imag[vi] - imag[vi + halfsize];

                    real[vi] += real[vi + halfsize];   // x = x + y
                    imag[vi] += imag[vi + halfsize];
                    real[vi + halfsize] = re * trigtables[ti] + im * trigtables[ti + 4];   // y = w(x - y)
                    imag[vi + halfsize] = im * trigtables[ti] - re * trigtables[ti + 4];
                }
            }
        }

        trigtables += size;
    }

#ifdef VERBOSE
        printf("\nGS Round #%d: {\n", floor_log2(4));
        for (int ii = 0; ii < n/2; ii++)
            printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

    // size 4 merge (special)
    if (n >= 4)
    {
        if (!inverse)
        {
            for (i = 0; i < n; i += 4)
            {
                // even indices (same as CT)
                re = real[i + 2];   // y
                im = imag[i + 2];
                real[i + 2] = real[i] - re;   // y = x - y
                imag[i + 2] = imag[i] - im;
                real[i] += re;   // x = x + y
                imag[i] += im;

                // odd indices (times i)
                re =  (imag[i + 1] - imag[i + 3]);   // -i(x - y)
                im = -(real[i + 1] - real[i + 3]);
                real[i + 1] += real[i + 3];   // x = x + y
                imag[i + 1] += imag[i + 3];
                real[i + 3] = re;   // y = i(x - y)
                imag[i + 3] = im;
            }
        }
        else
        {
            for (i = 0; i < n; i += 4)
            {
                // even indices (same as CT)
                re = real[i + 2];   // y
                im = imag[i + 2];
                real[i + 2] = real[i] - re;   // y = x - y
                imag[i + 2] = imag[i] - im;
                real[i] += re;   // x = x + y
                imag[i] += im;

                // odd indices (times -i)
                re = -(imag[i + 1] - imag[i + 3]);   // i(x - y)
                im =  (real[i + 1] - real[i + 3]);
                real[i + 1] += real[i + 3];   // x = x + y
                imag[i + 1] += imag[i + 3];
                real[i + 3] = re;   // y = i(x - y)
                imag[i + 3] = im;
            }
        }
    }

#ifdef VERBOSE
        printf("\nGS Round #%d: {\n", floor_log2(2));
        for (int ii = 0; ii < n/2; ii++)
            printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

    // size 2 merge (special; same as CT)
    if (n >= 2)
    {
        for (i = 0; i < n; i += 2)
        {
            // addition and subtraction
            re = real[i];
            im = imag[i];

            real[i] += real[i + 1];   // x = x + y
            imag[i] += imag[i + 1];
            real[i + 1] = re - real[i + 1];   // y = x - y
            imag[i + 1] = im - imag[i + 1];
        }
    }

    if (bitrev)
    {
#ifdef VERBOSE
    printf("\nBefore bit-rev:\n");
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

        const size_t *const bitreversed = tables->bit_reversed;
        permute_ary(bitreversed, real,
                    imag,
                    n, inverse);
    }
    else if (inverse)   // n.b., division is a part of bitreverse
    {
        for (i = 0; i < n; i++)
        {
            real[i] /= n;
            imag[i] /= n;
        }
    }

#ifdef VERBOSE
    printf("\nGS Result:\n");
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

}


// ---   Helper functions   ----------------------------------------------------

static inline void permute_ary(const size_t *const permutation,
                               double *const re,
                               double *const im,
                               const size_t n,
                               const bool inverse)
{
    if (inverse)
    {
        for (size_t i = 0; i < n; i++)
        {
            size_t j = permutation[i];
            if (i <= j)
            {
                double tp0re = re[i];
                double tp1re = re[j];
                re[i] = tp1re / n;
                re[j] = tp0re / n;

                double tp0im = im[i];
                double tp1im = im[j];
                im[i] = tp1im / n;
                im[j] = tp0im / n;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < n; i++)
        {
            size_t j = permutation[i];
            if (i < j)
            {
                double tp0re = re[i];
                double tp1re = re[j];
                re[i] = tp1re;
                re[j] = tp0re;

                double tp0im = im[i];
                double tp1im = im[j];
                im[i] = tp1im;
                im[j] = tp0im;
            }
        }
    }
}

// n must be a multiple of 4
static double accurate_sine(uint64_t i, uint64_t n)
{
    if (n % 4 != 0)
    {
        return NAN;
    }
    else
    {
        int32_t neg = 0;  // Boolean
        // Reduce to full cycle
        i %= n;
        // Reduce to half cycle
        if (i >= n / 2)
        {
            neg = 1;
            i -= n / 2;
        }
        // Reduce to quarter cycle
        if (i >= n / 4)
            i = n / 2 - i;
        // Reduce to eighth cycle
        double val;
        if (i * 8 < n)
            val = sin(2 * M_PI * i / n);
        else
            val = cos(2 * M_PI * (n / 4 - i) / n);
        // Apply sign
        return neg ? -val : val;
    }
}

static int32_t floor_log2(size_t n)
{
    int32_t result = 0;
    for (; n > 1; n /= 2)
        result++;
    return result;
}

static size_t reverse_bits(size_t x, uint32_t n)
{
    size_t result = 0;
    uint32_t i;
    for (i = 0; i < n; i++, x >>= 1)
        result = (result << 1) | (x & 1);
    return result;
}


// ---   Destroy tables   ------------------------------------------------------

void tables_destroy(FftTables *const tables)
{
    if (tables == NULL)
        return;
    free(tables->bit_reversed);
    free(tables->cos_table);
    free(tables->sin_table);
    free(tables->ct_path_tables);
    free(tables->gs_path_tables);
    memset(tables, 0, sizeof(FftTables));
    free(tables);
}









/*---- Function implementations ----*/

// Returns a pointer to an opaque structure of FFT tables. n must be a power of 2.
void *fft_init(size_t n) {
	// Check size argument
	if (n <= 0 || (n & (n - 1)) != 0)
		return NULL;  // Error: Size is not a power of 2
	if (n / 2 > SIZE_MAX / sizeof(double) || n > SIZE_MAX / sizeof(size_t))
		return NULL;  // Error: Size is too large, which makes memory allocation impossible

	// Allocate structure
	struct FftTables *tables = malloc(sizeof(struct FftTables));
	if (tables == NULL)
		return NULL;
	tables->n = n;

	// Allocate arrays
	tables->bit_reversed = malloc(n * sizeof(size_t));
	tables->cos_table = malloc(n / 2 * sizeof(double));
	tables->sin_table = malloc(n / 2 * sizeof(double));
	if (tables->bit_reversed == NULL || tables->cos_table == NULL || tables->sin_table == NULL) {
		free(tables->bit_reversed);
		free(tables->cos_table);
		free(tables->sin_table);
		free(tables);
		return NULL;
	}

	// Precompute values and store to tables
	size_t i;
	int32_t levels = floor_log2(n);
	for (i = 0; i < n; i++)
		tables->bit_reversed[i] = reverse_bits(i, levels);
	for (i = 0; i < n / 2; i++) {
		double angle = 2 * M_PI * i / n;
		tables->cos_table[i] = cos(angle);
		tables->sin_table[i] = sin(angle);
	}
	return tables;
}


// Performs a forward FFT in place on the given arrays. The length is given by the tables struct.
void fft_transform(const void *tables, double *real, double *imag) {
	struct FftTables *tbl = (struct FftTables *)tables;
	size_t n = tbl->n;

	// Bit-reversed addressing permutation
	size_t *bitreversed = tbl->bit_reversed;
	size_t i;
	for (i = 0; i < n; i++) {
		size_t j = bitreversed[i];
		if (i < j) {
			double tp0re = real[i];
			double tp0im = imag[i];
			double tp1re = real[j];
			double tp1im = imag[j];
			real[i] = tp1re;
			imag[i] = tp1im;
			real[j] = tp0re;
			imag[j] = tp0im;
		}
	}

	// Cooley-Tukey decimation-in-time radix-2 FFT
	double *costable = tbl->cos_table;
	double *sintable = tbl->sin_table;
	size_t size;
	for (size = 2; size <= n; size *= 2) {
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (i = 0; i < n; i += size) {
			size_t j, k;
			for (j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				double tpre =  real[j+halfsize] * costable[k] + imag[j+halfsize] * sintable[k];
				double tpim = -real[j+halfsize] * sintable[k] + imag[j+halfsize] * costable[k];
				real[j + halfsize] = real[j] - tpre;
				imag[j + halfsize] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
}


// Deallocates the given structure of FFT tables.
void fft_destroy(void *tables) {
	if (tables == NULL)
		return;
	struct FftTables *tbl = (struct FftTables *)tables;
	free(tbl->bit_reversed);
	free(tbl->cos_table);
	free(tbl->sin_table);
	memset(tbl, 0, sizeof(struct FftTables));  // Prevent accidental memory reuse
	free(tbl);
}


// Returns the largest i such that 2^i <= n.
static int32_t floor_log2(size_t n) {
	int32_t result = 0;
	for (; n > 1; n /= 2)
		result++;
	return result;
}


// Returns the bit reversal of the n-bit unsigned integer x.
static size_t reverse_bits(size_t x, uint32_t n) {
	size_t result = 0;
	uint32_t i;
	for (i = 0; i < n; i++, x >>= 1)
		result = (result << 1) | (x & 1);
	return result;
}
