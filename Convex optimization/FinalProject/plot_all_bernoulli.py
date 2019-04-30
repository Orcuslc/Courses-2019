import matplotlib.pyplot as plt
import h5py
import numpy as np
import sys

true_sol = h5py.File("bernoulli.mat", "r")
x_true = np.asarray(true_sol["x"]).reshape((-1, 1))
y_true = np.asarray(true_sol["y"]).reshape((1, -1))
z_true = np.asarray(true_sol["z"])[0]
w_true = np.asarray(true_sol["w"])

save_path = "bernoulli"

methods = ['SSG', 'APG', 'PG1', 'PG2']
xerr = {}
yerr = {}
zerr = {}
werr = {}

accuracies = {}
objectives = {}

for method in methods:
	filename = 'results/bernoulli/{0}.h5'.format(method)

	f = h5py.File(filename, 'r')
	results = f['results']

	accuracies[method] = np.asarray(results['accuracy']).reshape((-1, 1))
	objectives[method] = np.asarray(results['objective']).reshape((-1, 1))

	variables = f['variables']

	x = np.asarray([np.linalg.norm(xi.reshape((-1, 1)) - x_true)/np.linalg.norm(x_true) for xi in variables['x']])
	y = np.asarray([np.linalg.norm(yi - y_true)/np.linalg.norm(y_true) for yi in variables['y']])
	z = np.asarray([np.linalg.norm(zi - z_true)/np.linalg.norm(z_true) for zi in variables['z']])
	w = np.asarray([np.linalg.norm(wi[0] - w_true)/np.linalg.norm(w_true) for wi in variables['w']])

	xerr[method] = x
	yerr[method] = y
	zerr[method] = z
	werr[method] = w

K = range(0, 201, 2)

# plot objective
plt.figure(0)
lines = []
for method in methods:
	l, = plt.semilogy(K, objectives[method], linewidth = 2, label = method)
	lines.append(l)
plt.legend(lines, methods)
plt.grid(True)
plt.xlabel("Iterations")
plt.ylabel("Objective")
plt.savefig("plots/{0}/objective.jpg".format(save_path))

# plot accuracy
plt.figure(1)
lines = []
for method in methods:
	l, = plt.plot(K, accuracies[method], linewidth = 2, label = method)
	lines.append(l)
plt.legend(lines, methods)
plt.grid(True)
plt.xlabel("Iterations")
plt.ylabel("Accuracy (%)")
plt.savefig("plots/{0}/accuracy.jpg".format(save_path))

# plot subgradients
plt.figure(2)
lines = []
for method in methods:
	l, = plt.semilogy(K, xerr[method], linewidth = 2, label = method)
	lines.append(l)
plt.legend(lines, methods)
plt.grid(True)
plt.xlabel("Iterations")
plt.ylabel(r'$\frac{\Vert x-x_t\Vert_2}{\Vert x_t\Vert_2} $')
plt.savefig("plots/{0}/xerr.jpg".format(save_path))

plt.figure(3)
lines = []
for method in methods:
	l, = plt.semilogy(K, yerr[method], linewidth = 2, label = method)
	lines.append(l)
plt.legend(lines, methods)
plt.grid(True)
plt.xlabel("Iterations")
plt.ylabel(r"$\frac{\Vert y-y_t\Vert_2}{\Vert y_t\Vert_2} $")
plt.savefig("plots/{0}/yerr.jpg".format(save_path))

plt.figure(4)
lines = []
for method in methods:
	l, = plt.semilogy(K, zerr[method], linewidth = 2, label = method)
	lines.append(l)
plt.legend(lines, methods)
plt.grid(True)
plt.xlabel("Iterations")
plt.ylabel(r"$\frac{\Vert z-z_t\Vert_2}{\Vert z_t\Vert_2} $")
plt.savefig("plots/{0}/zerr.jpg".format(save_path))


plt.figure(5)
lines = []
for method in methods:
	l, = plt.semilogy(K, werr[method], linewidth = 2, label = method)
	lines.append(l)
plt.legend(lines, methods)
plt.grid(True)
plt.xlabel("Iterations")
plt.ylabel(r"$\frac{\Vert w-w_t\Vert_2}{\Vert w_t\Vert_2}$")
plt.savefig("plots/{0}/werr.jpg".format(save_path))


plt.show()

# if __name__ == '__main__':
# 	N = np.arange(101)

# 	# for objective and accuracy
# 	fig, ax1 = plt.subplots()
# 	color = 'tab:red'
# 	ax1.set_xlabel('Iterations (w.r.t. full batch)')
# 	ax1.set_ylabel('Objective', color=color)
# 	ax1.semilogy(objective, color=color)
# 	ax1.tick_params(axis='y', labelcolor=color)

# 	ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# 	color = 'tab:blue'
# 	ax2.set_ylabel('Accuracy', color=color)  # we already handled the x-label with ax1
# 	ax2.plot(accuracy, color=color)
# 	ax2.tick_params(axis='y', labelcolor=color)
# 	fig.tight_layout()  # otherwise the right y-label is slightly clipped
# 	plt.grid(True)
# 	plt.show()

# 	# for subgradients
# 	l1, = plt.semilogy(dx, label = "dx")
# 	l2, = plt.semilogy(dy, label = "dy")
# 	l3, = plt.semilogy(dz, label = "dz")
# 	l4, = plt.semilogy(dw, label = "dw")
# 	plt.xlabel('Iterations (w.r.t. full batch)')
# 	plt.ylabel('Norm of subgradients')
# 	plt.legend()
# 	plt.grid(True)
# 	plt.show()