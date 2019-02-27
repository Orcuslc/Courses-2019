include("./approximation.jl");
include("./polynomial.jl");
include("./quadrature.jl");
using .Approximation, .Quadrature, .Polynomial


x = -1:.1:1. |> collect;
f = x -> Legendre(x, 2);

x, w = legendre_gauss_weights_nodes(3, -1., 1.);
print(f(x).*w)