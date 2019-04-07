#ifndef _1D_NONLINEAR_HYDRODYNAMICS
#define _1D_NONLINEAR_HYDRODYNAMICS
#endif

#include "cmath"
#include "string"
#include <fstream>
#include <valarray>

using namespace std;

// abbv for array
typedef valarray<double> arr;

// save path
extern string save_path;

// wave vector parameters
extern double k0, u1, sigma;

// equilibrium quantities
extern double rho0, u0, Gamma, p0;

// grid paramters
extern double dt;
extern int nx;
extern arr tspan, xspan;

// initial conditions
arr rho_initial(const arr& x);
arr u_initial(const arr& x);
arr p_initial(const arr& x);

// write data to file
void write(string path, arr& data);

// leapfrog scheme
void leapfrog(const arr& xspan, const arr& tspan, const int nx, const double dt);