from scipy.optimize import curve_fit as cfit
import numpy as np
import matplotlib.pyplot as plt

def func(x, uf, delta):
	return uf*np.sin(2*np.pi*x+delta)

deltas = []
ufs = []
nts = [16, 24, 32, 64, 128, 256]

for nt in nts:
	with open("data_{0}.txt".format(nt)) as f:
		data = f.read().split('\n')[:-1]

	to_num = lambda row: np.asarray(list(map(lambda x: float(x), row.split(',')[:-1])))
	rho = list(map(to_num, data[0::3]))
	u = list(map(to_num, data[1::3]))
	p = list(map(to_num,data[2::3]))

	x = np.linspace(0., 1., 17)

	popt, _ = cfit(func, x, u[-1])
	ufs.append(popt[0])
	deltas.append(popt[1])

plt.plot(16/np.asarray(nts), np.asarray(deltas)/np.pi, '.-', label='delta')
plt.xlim(max(16/np.asarray(nts)), min(16/np.asarray(nts)))
plt.xticks(16/np.asarray(nts))
plt.xlabel("c Dt/Dx")
plt.ylabel("delta/pi")
plt.show()

plt.plot(16/np.asarray(nts), np.asarray(ufs)/0.1, '.-')
plt.xlim(max(16/np.asarray(nts)), min(16/np.asarray(nts)))
plt.xticks(16/np.asarray(nts))
plt.xlabel("c Dt/Dx")
plt.ylabel("uf/u0")
plt.show()