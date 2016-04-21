#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "lwe.h"
#include "multiplication.h"
#include "polynomials.h"

#include "lwesamples.h"
#include "lweparams.h"

using namespace std;



// **********************************************************************************
// ********************************* MAIN *******************************************
// **********************************************************************************
double approxEquals(Torus32 a, Torus32 b) { return abs(a-b)<10; }

int main(int argc, char** argv) {
    
    LWEParams* params = new_LWEParams(512, 0.2, 0.5); //les deux alpha mis un peu au hasard
    int n = params->n;
    LWEKey* key = new_LWEKey(params);
    LWESample* cipher = new_LWESample(params);
    Torus32 mu = dtot32(0.5);
    //Attention, 1<<30 correspond au message 0.25!! Ila: t'as raison!
    double alpha = 0.0625;
    Torus32 phi;
    double message;  
    int Msize = 2;  

    lweKeyGen(key);
    lweSymEncrypt(cipher, mu, alpha, key);
    cout << "a = [";
    for (int i = 0; i < n-1; ++i) cout << t32tod(cipher->a[i]) << ", ";
    cout << t32tod(cipher->a[n-1]) << "]" << endl;
    cout << "b = " << t32tod(cipher->b) << endl;

    phi = lwePhase(cipher, key);
    cout << "phi = " << t32tod(phi) << endl;
    message = lweSymDecrypt(cipher, key, Msize);
    cout << "message = " << t32tod(message) << endl;

    //lwe crash test
    int failures = 0;
    int trials = 1000;
    for (int i=0; i<trials; i++) {
        Torus32 input = dtot32((i%3)/3.);
	lweKeyGen(key);
	lweSymEncrypt(cipher, input, 0.047, key); // Ila: attention au niveau de bruit!!! à voir (0.06 n'est pas le bon je crois, 0.047 marche parfaitement)
	phi = lwePhase(cipher, key);
	Torus32 decrypted = lweSymDecrypt(cipher, key, 3);
	if ( !approxEquals(input,decrypted) ) {
	    cerr << "WARNING: the msg " << t32tod(input) << " gave phase " << t32tod(phi) << " and was incorrectly decrypted to " << t32tod(decrypted) << endl;
	    failures++;
	}
    }
    cout << "There were " << failures << " failures out of " << trials << " trials" << endl;
    cout << "(it might be normal)" << endl;

    return 0;
}
