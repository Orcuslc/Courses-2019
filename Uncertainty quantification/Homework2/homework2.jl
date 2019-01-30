include("./polynomial.jl");
using .Polynomial;

# x = (-1.:.01:1.) |> collect;
# y1 = Legendre(x, 5);
# y2 = Legendre_Normalized(x, 5);

# using PyCall;
# @pyimport matplotlib.pyplot as plt
# for i in 1:6
# 	plt.plot(x, y2[:, i]);
# end
# plt.show();

f1 = x -> abs(sin(pi*x));
f2 = x -> abs(x);
f3 = x -> cos(pi*x);

e1 = zeros(21, 1);
e2 = zeros(21, 1);
e3 = zeros(21, 1);

for i in 1:21
	e1[i] = Legendre_Normalized_Projection_Error(f1, i-1);
	e2[i] = Legendre_Normalized_Projection_Error(f2, i-1);
	e3[i] = Legendre_Normalized_Projection_Error(f3, i-1);
end

x = 0:20;

using PyCall;
@pyimport matplotlib.pyplot as plt
plt.plot(x, e1, label = "|sin(pi*x)|")
plt.semilogy(x, e2, label = "|x|")
plt.legend()
plt.xticks(0:2:20)
plt.xlabel("N")
plt.ylabel("projection error")
plt.grid(true)
plt.show()

plt.figure()
plt.semilogy(x, e3, label = "cos(pi*x)")
plt.legend()
plt.xticks(0:2:20)
plt.xlabel("N")
plt.ylabel("projection error")
plt.grid(true)
plt.show()