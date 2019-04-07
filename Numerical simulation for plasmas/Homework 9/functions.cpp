#include "1d_nonlinear_hydrodynamics.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <valarray>
#include <vector>
using namespace std;

arr central_difference(const arr& f, const double dx, const int boundary) {
    // compute derivative w.r.t. x, i.e., df/dx at each x point by central difference
    // boundary: 3 for periodic boundary conditions
    int len = f.size();
    arr dfdx(0., len);

    switch(boundary) {
        case 3: // periodic
            dfdx[slice(1, len-2, 1)] = (f[slice(2, len-1, 1)] - f[slice(0, len-3, 1)])/(2*dx);
            dfdx[len-1] = dfdx[0] = (f[1] - f[len-2])/(2*dx);
            break;
    }
    return dfdx;
}

arr lax_update(const arr& u, const arr& dp, const double dx, const double dt) {
    // lax method for du/dt = dp/dx;
    int len = u.size();
    arr u1(0., len);

    u1[slice(1, len-2, 1)] = 0.5*(u[slice(2, len-1, 1)]+u[slice(0, len-3, 1)]) + valarray<double>(dp[slice(1, len-2, 1)])*dt;
    u1[0] = 0.5*(u[1]+u[len-2]) + dp[0]*dt;
    u1[len-1] = 0.5*(u[1]+u[len-2]) + dp[len-1]*dt;
    return u1;
}

arr leapfrog_update(const arr& u, const arr& dp, const double dt) {
    // leapfrog update for du/dt = dp/dx;
    arr u1 = u+2*dt*dp;
    return u1;
}

void leapfrog(const arr& xspan, const arr& tspan, const int nx, const double dt) {
    // x nodes
    double dx = (xspan[1]-xspan[0])/(double)nx;
    arr x(0., nx+1);
    for(int i = 0; i < nx; i++) {
        x[i] = xspan[0]+(double)i*dx;
    }
    x[nx] = xspan[1];

    // initial
    arr rho = rho_initial(x), u = u_initial(x), p = p_initial(x);
    write(save_path, rho);
    write(save_path, u);
    write(save_path, p);

    // Compute the next step with Lax method
    // 1. compute central difference
    arr drhodx = central_difference(rho, dx, 3);
    arr dudx = central_difference(u, dx, 3);
    arr dpdx = central_difference(p, dx, 3);
    // 2. compute partial derivatives w.r.t. x
    arr drho = -(u*drhodx+rho*dudx);
    arr du = -(u*dudx+1./rho*dpdx);
    arr dp = -(u*dpdx+Gamma*p*dudx);
    // 3. update with lax step
    arr rho1 = lax_update(rho, drho, dx, dt);
    arr u1 = lax_update(u, du, dx, dt);
    arr p1 = lax_update(p, dp, dx, dt);
    write(save_path, rho1);
    write(save_path, u1);
    write(save_path, p1);

    // update with leapfrog method
    int nt = floor((tspan[1] - tspan[0])/dt);
    for(int n = 0; n < nt; n++) {
        // 1. compute central difference
        drhodx = central_difference(rho1, dx, 3);
        dudx = central_difference(u1, dx, 3);
        dpdx = central_difference(p1, dx, 3);

        // 2. compute partial derivatives w.r.t. x
        drho = -(u1*drhodx+rho1*dudx);
        du = -(u1*dudx+1./rho1*dpdx);
        dp = -(u1*dpdx+Gamma*p1*dudx);

        // 3. update with leapfrog step
        rho = leapfrog_update(rho, drho, dt);
        u = leapfrog_update(u, du, dt);
        p = leapfrog_update(p, dp, dt);

        // 4. write
        write(save_path, rho);
        write(save_path, u);
        write(save_path, p);

        // 5. swap
        rho.swap(rho1);
        u.swap(u1);
        p.swap(p1);
    } 
}

void write(string path, arr& data) {
    ofstream ofs (path, ios::app);  
    for(auto c: data) {
        ofs << c << ",";
    }
    ofs << "\n";
}