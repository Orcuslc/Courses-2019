#include "1d_linear_hydrodynamics.h"
#include "cmath"
#include <valarray>
#include "string"
using namespace std;

double k0 = 2*M_PI;
double u0 = 0.1;
double sigma = 1.0;

valarray<double> tspan = {0., 1.};
valarray<double> xspan = {0., 1.};

int nx = 128;
double dx = (tspan[1]-tspan[0])/(double)nx;
double dt = 1./128.;

// output format: 1-3 row: rho, u, p at t0, etc.
string save_path = "data.txt";