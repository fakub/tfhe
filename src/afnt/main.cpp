#include <cstdio>
#include <cstdint>

#include "nayuki-fft.h"

int main(int argc, char ** argv)
{
#define P_SIZE 8
    int32_t p[P_SIZE] = {1,3,-2,0,   -4,-3,2,-1};
    double pior[P_SIZE], pioi[P_SIZE];
    //  FFT(p) = [(-4.0+0.0i), (8.536-0.95i), (-3.0-1.0i), (1.464-8.95i), (-2.0+0.0i), (1.464+8.95i), (-3.0+1.0i), (8.536+0.95i)]
    // AFNT(p) = [(3.015+6.006i), (2.588-6.996i), (5.069+1.004i), (-6.672-1.994i)]


    // ---   FFT(p)   ----------------------------------------------------------

    // pre-compute FFT tables (omegas)
    void * omegas = fft_init(P_SIZE);
    if (omegas == NULL) return -1;

    // copy/typecast input
    for (int i = 0; i < P_SIZE; i++)
    {
        pior[i] = (double)p[i];
        pioi[i] = (double)0.0;
    }

    // FFT
    fft_transform(omegas, &pior[0], &pioi[0]);

    // destroy tables
    fft_destroy(omegas);


    // ---   Print results   ---------------------------------------------------

    printf("    p  = [");
    for (int i = 0; i < P_SIZE-1; i++)
        printf("%d, ", p[i]);
    printf("%d]\n", p[P_SIZE-1]);

    printf("FFT(p) = [");
    for (int i = 0; i < P_SIZE-1; i++)
        printf("(%.3f+%.3fi), ", pior[i], pioi[i]);
    printf("(%.3f+%.3fi)]\n", pior[P_SIZE-1], pioi[P_SIZE-1]);

    return 0;
}
