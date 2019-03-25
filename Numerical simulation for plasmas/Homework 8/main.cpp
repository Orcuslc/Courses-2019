#include "1d_linear_hydrodynamics.h"
#include<valarray> 
#include<iostream>
using namespace std;

int main() {
	leapfrog(xspan, tspan, nx, dt);
	return 0; 
}