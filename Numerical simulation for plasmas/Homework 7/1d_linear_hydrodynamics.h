#ifndef _1D_LINEAR_HYDRODYNAMICS
#define _1D_LINEAR_HYDRODYNAMICS
#endif

#include "cmath"
#include <valarray>
#include "string"
#include <fstream>
using namespace std;

// wave factor parameters
extern double k0, u0, sigma;

// initializing functions
inline valarray<double> rho_initial(valarray<double> x) {
	valarray<double> rho = sigma*u0*sin(k0*x);
	return rho;
}

inline valarray<double> u_initial(const valarray<double>& x) {
	valarray<double> u = u0*sin(k0*x);
	return u;
}

inline valarray<double> p_initial(const valarray<double>& x) {
	valarray<double> p = sigma*u0*sin(k0*x);
	return p;
}

// grid parameters
extern valarray<double> tspan, xspan;
extern double dt, dx;
extern int nx;

// A step of lax method with periodic boundary conditions for equation
// du/dt = -dp/dx
valarray<double> lax_step(const valarray<double>& u, const valarray<double>& p, const double dx, const double dt);

// save path
extern string save_path;
void write(string path, valarray<double>& data);

// lax method for the equations
// drho/dt = -du/dx
// du/dt = -dp/dx
// dp/dt = -du/dx
void lax(const valarray<double>& xspan, const valarray<double>& tspan, const int nx, const double dt); 