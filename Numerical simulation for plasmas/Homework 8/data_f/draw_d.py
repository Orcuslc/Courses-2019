import matplotlib.pyplot as plt
import numpy as np

for nt in [16, 24, 32, 64, 128, 256]:
	with open("data_{0}.txt".format(nt)) as f:
		data = f.read().split('\n')[:-1]

	to_num = lambda row: np.asarray(list(map(lambda x: float(x), row.split(',')[:-1])))
	rho = list(map(to_num, data[0::3]))
	u = list(map(to_num, data[1::3]))
	p = list(map(to_num,data[2::3]))

	x = np.linspace(0., 1., 17)

	plt.plot(x, u[-1], label = "t=100, n="+str(nt))
plt.legend()
plt.show()