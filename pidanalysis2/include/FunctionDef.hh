#ifndef FunctionDef_h
#define FunctionDef_h

#include <cmath>
#include <iostream>

#include "Rtypes.h"
#include "TF1.h"
#include "Vector3D.hh"

// calculate distance between ray and ditection point
double distance_func(double* xx, double* par);

// quantum efficiency for MPPC
double quantumEff(double* wlen, double* par);

// refractive index of Aerogel
double refractive_index(double* wlen, double* par);

// transmittance of Aerogel
double transmittance(double* wlen, double* par);
double transmittance2(double* wlen, double* par);

// calculate all detected photon 
double allDetectedPhoton(double mass, double momentum, double thickness, TF1* qeff, TF1* transmit,
                         TF1* r_index, double geo_eff, double wlenmin = 300.0, double wlenmax = 910.0);

// probability distribution of cherenkov angle
double angle_probDis(double theta, double theta_theory, double sigma);
//{
//    return std::exp(-0.5*std::pow((theta-theta_theory)/sigma,2.0))/(std::pow(2.0*M_PI,3.0/2)*sigma);
//}


#endif