import matplotlib.pyplot as plt
import h5py
import numpy as np
import sys

# save_path = "rc1.binary"
save_path = "real-sim"

objectives = {}
accuracies = {}
dxs = {}
dys = {}
dzs = {}
dws = {}

methods = ['SSG', 'APG', 'PG1', 'PG2']

for method in methods:
	filename = 'results/{1}/{0}.h5'.format(method, save_path)

	f = h5py.File(filename, 'r')
	results = f['results']
	variables = f['variables']

	objective = results['objective'][0]
	accuracy = results['accuracy'][0]

	dx = results['dx']
	dy = results['dy']
	dz = results['dz']
	dw = results['dw']

	# norm of each subgradient
	dx = [np.linalg.norm(dxi.reshape((-1, 1))) for dxi in dx]
	dy = [np.linalg.norm(dyi.reshape((-1, 1))) for dyi in dy]
	dz = [np.linalg.norm(dzi.reshape((-1, 1))) for dzi in dz]
	dw = [np.linalg.norm(dwi.reshape((-1, 1))) for dwi in dw]

	# store results in dictionary
	objectives[method] = objective
	accuracies[method] = accuracy
	dxs[method] = dx
	dys[method] = dy
	dzs[method] = dz
	dws[method] = dw 

K = range(0, 201, 2)

# plot objective
plt.figure(0)
lines = []
for method in methods:
	l, = plt.semilogy(K[:len(objectives[method])], objectives[method], linewidth = 2, label = method)
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
	l, = plt.plot(K[:len(objectives[method])], accuracies[method], linewidth = 2, label = method)
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
	l, = plt.semilogy(K[:len(objectives[method])][1:], dxs[method][1:], linewidth = 2, label = method)
	lines.append(l)
plt.legend(lines, methods)
plt.grid(True)
plt.xlabel("Iterations")
plt.ylabel(r'$\Vert \partial_x\Vert_2 $', fontsize = 14)
plt.savefig("plots/{0}/dx.jpg".format(save_path))

plt.figure(3)
lines = []
for method in methods:
	l, = plt.semilogy(K[:len(objectives[method])][1:], dys[method][1:], linewidth = 2, label = method)
	lines.append(l)
plt.legend(lines, methods)
plt.grid(True)
plt.xlabel("Iterations")
plt.ylabel(r"$\Vert \partial_y\Vert_2 $")
plt.savefig("plots/{0}/dy.jpg".format(save_path))

plt.figure(4)
lines = []
for method in methods:
	l, = plt.semilogy(K[:len(objectives[method])][1:], dzs[method][1:], linewidth = 2, label = method)
	lines.append(l)
plt.legend(lines, methods)
plt.grid(True)
plt.xlabel("Iterations")
plt.ylabel(r"$\Vert \partial_z\Vert_2 $")
plt.savefig("plots/{0}/dz.jpg".format(save_path))


plt.figure(5)
lines = []
for method in methods:
	l, = plt.semilogy(K[:len(objectives[method])][1:], dws[method][1:], linewidth = 2, label = method)
	lines.append(l)
plt.legend(lines, methods)
plt.grid(True)
plt.xlabel("Iterations")
plt.ylabel(r"$\Vert \partial_w\Vert_2 $")
plt.savefig("plots/{0}/dw.jpg".format(save_path))


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