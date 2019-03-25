from scipy.optimize import curve_fit as cfit
import numpy as np
import matplotlib.pyplot as plt

def func(x, delta):
	return 0.1*np.sin(2*np.pi*x+delta)

deltas = []
nts = [128, 192, 256, 512, 1024, 2048]

for nt in nts:
	with open("data_{0}.txt".format(nt)) as f:
		data = f.read().split('\n')[:-1]

	to_num = lambda row: np.asarray(list(map(lambda x: float(x), row.split(',')[:-1])))
	rho = list(map(to_num, data[0::3]))
	u = list(map(to_num, data[1::3]))
	p = list(map(to_num,data[2::3]))

	x = np.linspace(0., 1., 129)

	popt, _ = cfit(func, x, u[-1])
	deltas.append(popt[0])

plt.plot(128/np.asarray(nts), np.asarray(deltas)/np.pi, '.-')
plt.xlim(max(128/np.asarray(nts)), min(128/np.asarray(nts)))
plt.xticks(128/np.asarray(nts))
plt.xlabel("c Dt/Dx")
plt.ylabel("delta/pi")
plt.show()