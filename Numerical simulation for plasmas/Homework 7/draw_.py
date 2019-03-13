import matplotlib.pyplot as plt
import numpy as np

with open("data_1024.txt") as f:
	data = f.read().split('\n')[:-1]

to_num = lambda row: np.asarray(list(map(lambda x: float(x), row.split(',')[:-1])))
rho = list(map(to_num, data[0::3]))
u = list(map(to_num, data[1::3]))
p = list(map(to_num,data[2::3]))

x = np.linspace(0., 1., 129)
dt = 1./1024

get_index = lambda x: int(x/dt)

for t in [0.01*i for i in range(10)]:
	plt.plot(x, u[get_index(t)], label = "t="+str(t))
plt.legend()
plt.show()