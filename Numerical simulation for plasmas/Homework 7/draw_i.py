import matplotlib.pyplot as plt
import numpy as np

with open("data_64.txt") as f:
	data = f.read().split('\n')[:-1]

to_num = lambda row: np.asarray(list(map(lambda x: float(x), row.split(',')[:-1])))
rho = list(map(to_num, data[0::3]))
u = list(map(to_num, data[1::3]))
p = list(map(to_num,data[2::3]))

x = np.linspace(0., 1., 129)
dt = 1./64

get_index = lambda x: int(x/dt)

for t in [0, 0.25, 0.5, 0.75, 0.8125]:
	plt.plot(x, u[get_index(t)], label = "t="+str(t))
plt.legend()
plt.show()