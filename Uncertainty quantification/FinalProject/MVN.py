import numpy as np
from SVGD import SVGD

class MVN:
	def __init__(self, mu, Sigma):
		self.mu = mu
		self.Sigma = Sigma

	def dx_logp(self, theta):
		return -1*np.dot(theta - np.tile(self.mu, (theta.shape[0], 1)), self.Sigma)

if __name__ == '__main__':
	Sigma = np.array([[0.2260, 0.1652], [0.1652, 0.6779]])
	mu = np.array([-0.6871, 0.8010])

	model = MVN(mu, Sigma)

	x0 = np.random.normal(0, 1, [10, 2]);
	theta = SVGD().update(x0, model.dx_logp, 1000, 0.01)

	print(mu)
	print(np.mean(theta, axis = 0))