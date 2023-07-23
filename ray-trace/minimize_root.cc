#include <iostream>
#include <vector>
#include <algorithm>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

double distance_func(const double* xx){
    const double x = xx[0];
    const double y = xx[1];
    const double a = xx[2];

    return (x-a)*(x-a) + y*y + 3.0;
}

void minimize_root(const char* minName="Minuit2", const char* algoName="", int randomSeed = -1){

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    min->SetMaxFunctionCalls(10000);
    min->SetMaxIterations(10000);
    min->SetTolerance(1e-5);
    min->SetPrintLevel(1);

    ROOT::Math::Functor func(&distance_func, 2);
    double step[2] = {0.01,0.01};

    double variable[2] = {-1.0, 1.2};

    min->SetFunction(func);

    min->SetVariable(0, "x", variable[0], step[0]);
    min->SetVariable(1, "y", variable[1], step[1]);

    min->SetFixedVariable(2, "a", 2.2);
    
    min->Minimize();

    const double* xs = min->X();

    std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
              << min->MinValue()  << std::endl;
    
    
}