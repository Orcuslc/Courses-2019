import matplotlib.pyplot as plt
import numpy as np

# read data
with open("data.txt") as f:
	data = f.read().split('\n')[:-1]
to_num = lambda row: np.asarray(list(map(lambda x: float(x), row.split(',')[:-1])))
rho = list(map(to_num, data[0::3]))
u = list(map(to_num, data[1::3]))
p = list(map(to_num,data[2::3]))

# compute local sound speed
gamma = 5/3
c_s = [np.sqrt(gamma*p_/rho_) for (p_, rho_) in zip(p, rho)]

dt = 1./128
get_index = lambda x: int(x/dt)
x = np.linspace(0., 1., 129)

# # problem (f)
# for t in [0.25*i for i in range(5)]:
# 	plt.plot(x, u[get_index(t)], label = "t="+str(t))
# plt.legend()
# plt.show()

# Riemann invariants
J_plus = [ui + 2/(gamma-1)*c_si for (ui, c_si) in zip(u, c_s)]
J_minus = [ui - 2/(gamma-1)*c_si for (ui, c_si) in zip(u, c_s)]

# # problem (h)
# for t in [0.25*i for i in range(1, 3)]:
# 	plt.plot(x, J_plus[get_index(t)], label = "J+, t = "+str(t))
# 	plt.plot(x, J_minus[get_index(t)], label = "J-, t = "+str(t))
# plt.legend()
# plt.show()

# since J+ is not a constant, we need to propogate the x position of each variable at speed u+c_s. 

import math
advance_x = lambda x, t, u, c_s: (x+(u+c_s)*t) - np.floor((x+(u+c_s)*t))

u_analytical = [(J_plus_i+J_minus_i)/2 for (J_plus_i, J_minus_i) in zip(J_plus, J_minus)]

# # problem (i)
# for t in [0.25, 1.0]:
# 	plt.plot(x, u[get_index(t)], label = "u_numerical, t = "+str(t))

# 	x_adv = advance_x(x, t, u[get_index(t)], c_s[get_index(t)])
# 	# x_ana, u_ana = zip(*sorted(zip(x_adv, u_analytical[get_index(t)])))
# 	x_ana, u_ana = zip(*sorted(zip(x_adv, u_analytical[get_index(t)])))
# 	plt.plot(x_ana, u_ana, label = "u_analytical, t = "+str(t))
# plt.legend()
# plt.show()

# problem (j)
plt.plot(x, u[get_index(2.0)], label = "u_numerical, t = 2.0")
x_adv = advance_x(x, 2.0, u[get_index(2)], c_s[get_index(2)])
# x_ana, u_ana = zip(*sorted(zip(x_adv, u_analytical[get_index(t)])))
x_ana, u_ana = zip(*sorted(zip(x_adv, u_analytical[get_index(2)])))
plt.plot(x, u_ana, label = "u_analytical, t = 2.0")
plt.legend()
plt.show()