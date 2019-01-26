import numpy as np
from matplotlib import pyplot as plt
from orthnet import *

x1 = np.linspace(-1, 1, 100).reshape((-1, 1))
x2 = np.linspace(0, 5, 1000).reshape((-1, 1))
x3 = np.linspace(-2, 2, 1000).reshape((-1, 1))
legendre = Legendre(x1, 5)
laguerre = Laguerre(x2, 5)
hermite = Hermite(x3, 5)
chebyshev = Chebyshev(x1, 5)

plt.figure(1)
for i in range(6):
	plt.plot(x1, legendre.list[i], label = "n="+str(i))
plt.legend()
plt.title("Legendre Polynomials")
plt.show()

plt.figure(2)
for i in range(6):
	plt.plot(x2, laguerre.list[i], label = "n="+str(i))
plt.legend()
plt.title("Laguerre Polynomials")
plt.show()

plt.figure(3)
for i in range(6):
	plt.plot(x3, hermite.list[i], label = "n="+str(i))
plt.legend()
plt.title("Hermite Polynomials")
plt.show()

plt.figure(4)
for i in range(6):
	plt.plot(x1, chebyshev.list[i], label = "n="+str(i))
plt.legend()
plt.title("Chebyshev Polynomials")
plt.show()