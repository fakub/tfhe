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

#ifndef FFT_NAYUKI_FFNT_H
#define FFT_NAYUKI_FFNT_H


#ifdef __cplusplus
extern "C" {
#endif

void* fft_init(size_t n);

void* fft_init_reverse(size_t n);

void fft_transform(const void *tables, double *real, double *imag);

void fft_transform_reverse(const void *tables, double *real, double *imag);

void fft_destroy(void *tables);

#ifdef __cplusplus
}
#endif

#endif //FFT_NAYUKI_FFNT_H
