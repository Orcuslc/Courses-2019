import numpy as np
import torch
import torch.nn as nn

class SVGD:
	"""
		Stein Variational Gradient Descent training engine for Bayesian Neural Networks

		BNN: a list of Bayesian Neural Networks (nn.Module)
	"""
	def __init__(self, BNN):
		self.BNN = BNN

		# number of particles (models) to approx the distribution
		self.n_particles = len(BNN)


	def _RBF_kernel(self, theta):
		"""
		Compute the RBF Kernel K(x, x') and gradient w.r.t. x, 
		i.e., (8) in SVGD paper
		"""
		pdist = self._pairwise_dist(theta)

		# median trick, described in Section 5 of SVDG paper
		h_square = 0.5*pdist.median()/np.log(self.n_particles)

		# RBF kernel
		kernel = torch.exp(pdist/(2*h_square))

		# gradient of kernel w.r.t. x
		grad_kernel = 1/h_square*torch.mm(kernel.sum(1).diag()-kernel, theta)

		return kernel, grad_kernel


	def _pairwise_dist(self, theta):
		"""
		Compute pairwise distance of theta. Let x_i be the i^th sample of theta (i^th row), then D(x_i, x_j) = || x_i - x_j ||_2^2
		"""
		return torch.norm(theta[:, None] - theta, dim = 2, p = 2)**2

	def train(self, data, epochs):
		# set BNN to training mode
		self.BNN.train(mode = True)

		# main iteration
		for batch_index, (x, y) in enumerate(data):
			
