/*
 * Fast Fourier transform
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
#ifndef FFNT_H
#define FFNT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdbool.h>

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

/*******************************************************************************
 * Initialize FFT tables.
 * */
FftTables *  fft_init(const size_t n);
/*******************************************************************************
 * Initialize inverse FFT tables.
 * */
FftTables * ifft_init(const size_t n);

/*******************************************************************************
 * Run FFT (via Cooley-Tukey data path).
 * */
void  fft_transform(const FftTables *const fft_tables,
                    double *const real,
                    double *const imag);
/*******************************************************************************
 * Run inverse FFT (via Cooley-Tukey data path).
 * */
void ifft_transform(const FftTables *const ifft_tables,
                    double *const real,
                    double *const imag);

/*******************************************************************************
 * Run FFNT (does not apply bit-reverse).
 * */
void ffnt_transform(const FftTables *const ffnt_tables_2N,
                    const FftTables *const ffnt_tables_N_2,
                    double *const real,
                    double *const imag);
/*******************************************************************************
 * Run inverse FFNT (does not apply bit-reverse).
 * */
void iffnt_transform(const FftTables *const ffnt_tables_2N,
                     const FftTables *const iffnt_tables_N_2,
                     double *const real,
                     double *const imag);

/*******************************************************************************
 * Run FFT via Cooley-Tukey data path.
 *
 * Choose whether to use bit-reverse of the input using the +bitrev+ switch.
 *
 * */
void  fft_transform_CT(const FftTables *const fft_tables,
                       double *const real,
                       double *const imag,
                       const bool bitrev);
/*******************************************************************************
 * Run inverse FFT via Cooley-Tukey data path.
 *
 * Choose whether to use bit-reverse of the input using the +bitrev+ switch.
 *
 * */
void ifft_transform_CT(const FftTables *const fft_tables,
                       double *const real,
                       double *const imag,
                       const bool bitrev);
/*******************************************************************************
 * Run FFT via Gentleman-Sande data path.
 *
 * Choose whether to use bit-reverse of the output using the +bitrev+ switch.
 *
 * */
void  fft_transform_GS(const FftTables *const fft_tables,
                       double *const real,
                       double *const imag,
                       const bool bitrev);
/*******************************************************************************
 * Run inverse FFT via Gentleman-Sande data path.
 *
 * Choose whether to use bit-reverse of the output using the +bitrev+ switch.
 *
 * */
void ifft_transform_GS(const FftTables *const fft_tables,
                       double *const real,
                       double *const imag,
                       const bool bitrev);

/*******************************************************************************
 * Final cleanup.
 * */
void tables_destroy(void *const tables);

#ifdef __cplusplus
}
#endif

#endif   // FFNT_H
