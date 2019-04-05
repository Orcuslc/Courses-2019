#include "1d_nonlinear_hydrodynamics.h"
#include "cmath"
#include <valarray>
#include "string"

using namespace std;

// wave vector parameters
double k0 = 2*M_PI;
double u1 = 0.08;
double sigma = 1.0;

// equilibrium quantities
double rho0 = 1.0;
double u0 = 0.0;
double Gamma = 5.0/3.0;
double p0 = 1.0/Gamma;

// grid parameters
arr tspan = {0., 2.};
arr xspan = {0., 1.};
int nx = 128;
double dt = 1./512.;
double dx = (xspan[1] - xspan[0])/(double)nx;

// initial functions
arr rho_initial(const arr& x) {
	arr rho = rho0 + sigma*u1*sin(k0*x);
	return rho;
}
arr u_initial(const arr& x) {
	arr u = u0 + u1*sin(k0*x);
	return u;
}
arr p_initial(const arr& x) {
	arr p = p0 + sigma*u1*sin(k0*x);
	return p;
}

string save_path = "data.txt";