cmake_minimum_required(VERSION 3.0)

# This is the nayuki ffnt processor for the tfhe library
enable_language(ASM)

set(SRCS_PORTABLE
    fft-portable.c
    fft_processor_nayuki_FFNT.cpp
    lagrangehalfc_impl.cpp
    )

set(SRCS_AVX
    fft-x8664-avx-aux.c
    fft-x8664-avx.s
    fft-x8664-avx-reverse.s
    fft_processor_nayuki_FFNT.cpp
    lagrangehalfc_impl.cpp
    )

set(HEADERS
    ffnt.h
    lagrangehalfc_impl.h
    )

if (ENABLE_NAYUKI_FFNT_PORTABLE)
    add_library(tfhe-fft-nayuki-FFNT-portable OBJECT ${SRCS_PORTABLE} ${HEADERS})
    set_property(TARGET tfhe-fft-nayuki-FFNT-portable PROPERTY POSITION_INDEPENDENT_CODE ON)
endif (ENABLE_NAYUKI_FFNT_PORTABLE)

if (ENABLE_NAYUKI_FFNT_AVX)
    message( SEND_ERROR "Nayuki-FFNT processor not yet implemented." )
    #~ add_library(tfhe-fft-nayuki-FFNT-avx OBJECT ${SRCS_AVX} ${HEADERS})
    #~ set_property(TARGET tfhe-fft-nayuki-FFNT-avx PROPERTY POSITION_INDEPENDENT_CODE ON)
endif (ENABLE_NAYUKI_FFNT_AVX)
