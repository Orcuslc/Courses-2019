include("./approximation.jl");
include("./polynomial.jl");
include("./quadrature.jl");
using .Approximation, .Quadrature, .Polynomial


x = -1:1:1. |> collect;

f = x -> x;
pif = legendre_interpolation(f, x);

x1 = -1.:.2:1. |> collect;
print(pif(x1));