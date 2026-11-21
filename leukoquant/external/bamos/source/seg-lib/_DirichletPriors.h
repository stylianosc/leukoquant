#ifndef _SEG_DP_H
#define _SEG_DP_H

#include <math.h>
#define gamma 0.5772156649
//#define MaxSupport 15
#define DefaultProp 1E-6

using namespace std;

#include <iostream>
#include "_TreeEM_new.h"

inline float pow_int1(const float base,
                             int exp){
    if(exp==0){return 1;}
    float result = base;
    while (--exp){result *= base;}
    return result;
}

// Returns the value of the digamma function applied to x
float Digamma(float x);

// Returns the value taken by the derivative of the digamma function in x
float DigammaDer(float x);

// Returns the value of x such that psi(x) = y using the Newton method to invert the digamma function
float NewtonDigamma(float y);

// Returns the array of
float * AlphaCalculation(float * p);

float BandwidthCalculation(int * Counts,SEG_PARAMETERS * segment_param);

float GaussianKernel(float x1, float x2, float BW);

float * DirectKernelSmoothing(int * Counts,SEG_PARAMETERS * segment_param);

float * DerivativeLocalRegression(int i, int * Counts, int order, float BW, float * Coeffs);

float GradientDescentPoissonRegression(int i, int * Counts, int order,SEG_PARAMETERS * segment_param);

float * PoissonRegressionKernelSmoothing(int * Counts, int order, SEG_PARAMETERS* segment_param);

float * MultivariateJointHistogram(vector<int *> Counts,vector<int*> CountDistribution);

float * GaussianBlurring(float * CountHistogram, float gauss_std, vector<int> dim,bool NormalisedSum=1);



#endif
