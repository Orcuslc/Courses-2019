include("./approximation.jl");
include("./polynomial.jl");
include("./quadrature.jl");
using .Approximation, .Quadrature, .Polynomial


x = -1:1:1. |> collect;
x2 = -1.:.5:1. |> collect;
f = x -> x;

pif = legendre_interpolation(x, f(x));
qif = legendre_least_square(x, f(x), 1);

x1 = -1.:.2:1. |> collect;
print(pif(x1));
print(qif(x1))