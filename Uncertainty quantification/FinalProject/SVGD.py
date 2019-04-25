import numpy as np
from scipy.spatial.distance import pdist, squareform

class SVGD:
	def __init__(self):
		pass

	def gaussian_kernel(self, theta, h):
		pairwise_dist = squareform(pdist(theta, 'euclidean'))**2
		
		if h < 0: # use median trick: Page 183 of http://alex.smola.org/teaching/kernelcourse/day_2.pdf ? 
			h = np.median(pairwise_dist)
			h = np.sqrt(0.5*h/np.log(theta.shape[0]+1))

		# compute kernel: K(x_i, x_j) = exp(-1/(2*h^2) ||x_i - x_j||^2) 
		kxx = np.exp(-pairwise_dist/(2*h**2))

		# compute derivative of kernel w.r.t. x, i.e.
		# dkdx_i = sum_j K(x_i, x_j)(-1/h^2 (x_i - x_j))
		dkdx = 1/(h**2)*(-np.dot(kxx, theta) + theta*kxx.sum(axis = 1, keepdims = True))

		return kxx, dkdx



	def update(self, theta, dx_logp, n_iter = 1000, learning_rate = 1e-3, alpha = 0.9):
		# update by Adagrad with momentum
		fudge_factor = 1e-6
		historical_grad = 0.0

		for i in range(n_iter):
			# grad_x_logp(x)
			grad_logp = dx_logp(theta)

			# kernel matrix
			kxx, dkdx = self.gaussian_kernel(theta, h = -1)

			# gradient of theta
			phi = (np.dot(kxx, grad_logp) + dkdx)/theta.shape[0]

			# Adagrad
			if i == 0:
				historical_grad += phi**2
			else:
				historical_grad = alpha*historical_grad + (1-alpha)*phi**2
			theta += learning_rate*phi/(fudge_factor+np.sqrt(historical_grad))

		return theta