import matplotlib.pyplot as plt
import h5py
import numpy as np
import sys

dataset = sys.argv[1]
method = sys.argv[2]

filename = 'results/{1}/{0}.h5'.format(method, dataset)

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

if __name__ == '__main__':
	N = np.arange(101)

	# for objective and accuracy
	fig, ax1 = plt.subplots()
	color = 'tab:red'
	ax1.set_xlabel('Iterations (w.r.t. full batch)')
	ax1.set_ylabel('Objective', color=color)
	ax1.semilogy(objective, color=color)
	ax1.tick_params(axis='y', labelcolor=color)

	ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
	color = 'tab:blue'
	ax2.set_ylabel('Accuracy', color=color)  # we already handled the x-label with ax1
	ax2.plot(accuracy, color=color)
	ax2.tick_params(axis='y', labelcolor=color)
	fig.tight_layout()  # otherwise the right y-label is slightly clipped
	plt.grid(True)
	plt.show()

	# for subgradients
	l1, = plt.semilogy(dx, label = "dx")
	l2, = plt.semilogy(dy, label = "dy")
	l3, = plt.semilogy(dz, label = "dz")
	l4, = plt.semilogy(dw, label = "dw")
	plt.xlabel('Iterations (w.r.t. full batch)')
	plt.ylabel('Norm of subgradients')
	plt.legend()
	plt.grid(True)
	plt.show()