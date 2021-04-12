#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"

using namespace std;

void dieDramatically(string message)
{
    cerr << message << endl;
    abort();
}

void full_adder_MUX(LweSample *sum,
                    const LweSample *x,
                    const LweSample *y,
                    const int32_t nb_bits,
                    const TFheGateBootstrappingSecretKeySet *tfhe_keys)
{
    const LweParams *io_lwe_params = tfhe_keys->params->in_out_params;
    LweSample *carry    = new_LweSample_array(2, io_lwe_params);
    LweSample *temp     = new_LweSample_array(2, io_lwe_params);

    // init c = 0
    bootsSymEncrypt(carry, 0, tfhe_keys);

    for (int32_t i = 0; i < nb_bits; ++i)
    {
        // z_i = x_i XOR y_i XOR c
        bootsXOR(temp, x + i, y + i, &tfhe_keys->cloud);   // w_i = x_i XOR y_i (temp_0)
        bootsXOR(sum + i, temp, carry, &tfhe_keys->cloud);

        // c = (x_i XOR y_i) ? c : (x_i AND y_i)
        bootsAND(temp + 1, x + i, y + i, &tfhe_keys->cloud); // temp_1 = x_i AND y_i
        bootsMUX(carry + 1, temp, carry, temp + 1, &tfhe_keys->cloud);

        // correctness check
        bool mess1 = bootsSymDecrypt(temp, tfhe_keys);
        bool mess2 = bootsSymDecrypt(carry, tfhe_keys);
        bool mess3 = bootsSymDecrypt(temp + 1, tfhe_keys);
        bool messmux = bootsSymDecrypt(carry + 1, tfhe_keys);

        if (messmux != (mess1 ? mess2 : mess3))
        {
            cout << "ERROR!!! " << i << " - ";
            cout << t32tod(lwePhase(temp, tfhe_keys->lwe_key)) << " - ";
            cout << t32tod(lwePhase(carry, tfhe_keys->lwe_key)) << " - ";
            cout << t32tod(lwePhase(temp + 1, tfhe_keys->lwe_key)) << " - ";
            cout << t32tod(lwePhase(carry + 1, tfhe_keys->lwe_key)) << endl;
        }

        bootsCOPY(carry, carry + 1, &tfhe_keys->cloud);
    }
    bootsCOPY(sum + nb_bits, carry, &tfhe_keys->cloud);

    delete_LweSample_array(2, temp);
    delete_LweSample_array(2, carry);
}

void full_adder(LweSample *sum,
                const LweSample *x,
                const LweSample *y,
                const int32_t nb_bits,
                const TFheGateBootstrappingSecretKeySet *tfhe_keys)
{
    const LweParams *io_lwe_params = tfhe_keys->params->in_out_params;
    LweSample *carry    = new_LweSample_array(2, io_lwe_params);
    LweSample *temp     = new_LweSample_array(3, io_lwe_params);

    // init c = 0
    bootsSymEncrypt(carry, 0, tfhe_keys);

    for (int32_t i = 0; i < nb_bits; ++i)
    {
        // z_i = x_i XOR y_i XOR c
        bootsXOR(temp, x + i, y + i, &tfhe_keys->cloud);   // w_i = x_i XOR y_i (temp_0)
        bootsXOR(sum + i, temp, carry, &tfhe_keys->cloud);

        // c = (x_i AND y_i) XOR (c AND (x_i XOR y_i))
        bootsAND(temp + 1, x + i, y + i, &tfhe_keys->cloud);   // temp_1 = x_i AND y_i
        bootsAND(temp + 2, carry, temp, &tfhe_keys->cloud);    // temp_2 = c AND w_i
        bootsXOR(carry + 1, temp + 1, temp + 2, &tfhe_keys->cloud);
        bootsCOPY(carry, carry + 1, &tfhe_keys->cloud);
    }
    bootsCOPY(sum + nb_bits, carry, &tfhe_keys->cloud);

    delete_LweSample_array(3, temp);
    delete_LweSample_array(2, carry);
}

void comparison_MUX(LweSample *comp, const LweSample *x, const LweSample *y, const int32_t nb_bits,
                    const TFheGateBootstrappingSecretKeySet *tfhe_keys) {
    const LweParams *io_lwe_params = tfhe_keys->params->in_out_params;
    LweSample *carry = new_LweSample_array(2, io_lwe_params);
    LweSample *temp = new_LweSample(io_lwe_params);

    // init c = 1
    bootsSymEncrypt(carry, 1, tfhe_keys);

    // calc c = (x <= y)
    for (int32_t i = 0; i < nb_bits; ++i)
    {
        // c = w_i ? y_i : c (propagates as the final result)
        bootsXOR(temp, x + i, y + i, &tfhe_keys->cloud);   // w_i = x_i XOR y_i (temp)
        bootsMUX(carry + 1, temp, y + i, carry, &tfhe_keys->cloud);
        bootsCOPY(carry, carry + 1, &tfhe_keys->cloud);
    }
    bootsCOPY(comp, carry, &tfhe_keys->cloud);

    delete_LweSample(temp);
    delete_LweSample_array(2, carry);
}

#ifndef NDEBUG
extern const TLweKey *debug_accum_key;
extern const LweKey *debug_extract_key;
extern const LweKey *debug_in_key;
#endif



// =============================================================================
//
//  MAIN
//

int32_t main(int32_t argc, char **argv)
{
#ifndef NDEBUG
    cout << "DEBUG MODE!" << endl;
#endif
    const int32_t nb_bits = 16;
    const int32_t nb_trials = 3;

    // generate TFHE params
    int32_t minimum_lambda = 100;
    TFheGateBootstrappingParameterSet *tfhe_params = new_default_gate_bootstrapping_parameters(minimum_lambda);
    const LweParams *io_lwe_params = tfhe_params->in_out_params;
    // generate TFHE secret keys
    TFheGateBootstrappingSecretKeySet *tfhe_keys = new_random_gate_bootstrapping_secret_keyset(tfhe_params);
    // alloc samples
    LweSample *x    = new_LweSample_array(nb_bits, io_lwe_params);
    LweSample *y    = new_LweSample_array(nb_bits, io_lwe_params);
    LweSample *sum  = new_LweSample_array(nb_bits + 1, io_lwe_params);
    LweSample *comp = new_LweSample(io_lwe_params);


    // -------------------------------------------------------------------------
    //  main loop -- trials
    //
    for (int32_t trial = 0; trial < nb_trials; ++trial)
    {
        // generate random-bit samples
        for (int32_t i = 0; i < nb_bits; ++i)
        {
            bootsSymEncrypt(x + i, rand() % 2, tfhe_keys);
            bootsSymEncrypt(y + i, rand() % 2, tfhe_keys);
        }


        // run MUX addition   --------------------------------------------------
        //
        cout << "starting bootstrapping " << nb_bits << "-bits addition circuit (FA in MUX version)...trial " << trial
             << endl;

        clock_t begin1 = clock();
        full_adder_MUX(sum, x, y, nb_bits, tfhe_keys);
        clock_t end1 = clock();

        cout << "finished bootstrappings " << nb_bits << "-bits addition circuit (FA in MUX version)" << endl;
        cout << "total time (microsecs)... " << (end1 - begin1) << endl;

        // verify
        {
            bool messCarry = 0;
            for (int32_t i = 0; i < nb_bits; ++i)
            {
                bool messX = bootsSymDecrypt(x + i, tfhe_keys);
                bool messY = bootsSymDecrypt(y + i, tfhe_keys);
                bool messSum = bootsSymDecrypt(sum + i, tfhe_keys);

                if (messSum != (messX ^ messY ^ messCarry))
                {
                    cout << "ERROR!!! " << trial << "," << i << " - ";
                    cout << t32tod(lwePhase(x + i, tfhe_keys->lwe_key)) << " - ";
                    cout << t32tod(lwePhase(y + i, tfhe_keys->lwe_key)) << " - ";
                    cout << t32tod(lwePhase(sum + i, tfhe_keys->lwe_key)) << endl;
                }

                messCarry = messCarry ? (messX || messY) : (messX && messY);
            }
            bool messSum = bootsSymDecrypt(sum + nb_bits, tfhe_keys);
            if (messSum != messCarry)
                cout << "ERROR!!! " << trial << "," << nb_bits << endl;
        }


        // run no-MUX addition   -----------------------------------------------
        //
        cout << "starting bootstrapping " << nb_bits << "-bits addition circuit (FA)...trial " << trial << endl;

        clock_t begin2 = clock();
        full_adder(sum, x, y, nb_bits, tfhe_keys);
        clock_t end2 = clock();

        cout << "finished bootstrappings " << nb_bits << "-bits addition circuit (FA)" << endl;
        cout << "total time (microsecs)... " << (end2 - begin2) << endl;

        // verify
        {
            bool messCarry = 0;
            for (int32_t i = 0; i < nb_bits; ++i)
            {
                bool messX = bootsSymDecrypt(x + i, tfhe_keys);
                bool messY = bootsSymDecrypt(y + i, tfhe_keys);
                bool messSum = bootsSymDecrypt(sum + i, tfhe_keys);

                if (messSum != (messX ^ messY ^ messCarry))
                    cout << "ERROR!!! " << trial << "," << i << endl;

                messCarry = messCarry ? (messX || messY) : (messX && messY);
            }
            bool messSum = bootsSymDecrypt(sum + nb_bits, tfhe_keys);
            if (messSum != messCarry)
                cout << "ERROR!!! " << trial << "," << nb_bits << endl;
        }


        // run comparison   -----------------------------------------------
        //
        cout << "starting bootstrapping " << nb_bits << "-bits comparison...trial " << trial << endl;

        clock_t begin3 = clock();
        comparison_MUX(comp, x, y, nb_bits, tfhe_keys);
        clock_t end3 = clock();

        cout << "finished bootstrappings " << nb_bits << "-bits comparison" << endl;
        cout << "total time (microsecs)... " << (end3 - begin3) << endl;

        // verify
        {
            bool messCarry = 1;
            for (int32_t i = 0; i < nb_bits; ++i)
            {
                bool messX = bootsSymDecrypt(x + i, tfhe_keys);
                bool messY = bootsSymDecrypt(y + i, tfhe_keys);

                messCarry = (messX ^ messY) ? messY : messCarry;
            }
            bool messComp = bootsSymDecrypt(comp, tfhe_keys);
            if (messComp != messCarry)
                cout << "ERROR!!! " << trial << "," << nb_bits << endl;
        }
    }

    // cleanup
    delete_LweSample(comp);
    delete_LweSample_array(nb_bits + 1, sum);
    delete_LweSample_array(nb_bits, y);
    delete_LweSample_array(nb_bits, x);

    delete_gate_bootstrapping_secret_keyset(tfhe_keys);
    delete_gate_bootstrapping_parameters(tfhe_params);

    return 0;
}
