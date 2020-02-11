import torch
import torch.nn as nn
import numpy as np
from scipy.spatial.distance import pdist, squareform
import time
import copy

class BNN:
	"""
	Bayesian neural network:
		y = f(x, w) + n,
	`w` is the parameters (weights, biases) of the net, and `n` is the additional noise.

	Here we only consider the output-wise noise, i.e., 
		n = \sigma*\epsilon,
	where `\epsilon ~ N(0, 1)` is the gaussian noise. Define `\beta = 1/\sigma^2` be the noise precision, then
		n ~ N(0, \beta^-1).

	Prior:
		- p(w | \alpha) = N(w | 0, \alpha^-1)
		- p(\alpha) = Gamma(\alpha | a0, b0)
	=>	w ~ student T(w | 0, a0/b0, 2*a0)

		- p(n | \beta) = N(n | 0, \beta^-1)
		- p(\beta) = Gamma(\beta | a1, b1)
	
	Then
		- p(y | w, x, n) = \prod N(y_i | f(x_i, w), \beta^-1)

	Posterior:
		p(w, \alpha, \beta) = p(y | w, x, n) * p(w | \alpha) * p(\alpha) * p(\beta)
	
	Parameters:
		- model: the deterministic model (nn.Module)
		- N_particle: the number of particles (models) used to approx the posterior distribution
		- alpha_prior: [a0, b0]
		- beta_prior: [a1, b1]
	"""
	def __init__(self, model, N_particle, alpha_prior, beta_prior):
		self.N_particle = N_particle

		# initialize the particles
		models = []
		for i in range(N_particle):
			

	def _initialize(self, )