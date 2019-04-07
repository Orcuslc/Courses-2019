#include "1d_nonlinear_hydrodynamics.h" 
using namespace std;

int main() {
	leapfrog(xspan, tspan, nx, dt);
	return 0;
}