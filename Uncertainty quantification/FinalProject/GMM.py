import numpy as np
from SVGD import SVGD

class GMM1d:
	def __init__(self, mu, sigma, p):
		self.mu = mu
		self.sigma = sigma
		self.p = p

	def grad_logp(self, theta):
		prob_for_each_model = 1./np.sqrt(2*np.pi*self.sigma**2)*np.exp(-(theta - self.mu)**2/(2*self.sigma**2))
		return np.sum(self.p*prob_for_each_model*(-(theta - self.mu)/(self.sigma**2)), axis = 1, keepdims = True)/np.sum(self.p*prob_for_each_model, axis = 1, keepdims = True)

if __name__ == '__main__':
	mu = np.array([-2, 2])
	sigma = np.array([1, 1])
	p = np.array([1/3, 2/3])

	GMM = GMM1d(mu, sigma, p)
	x0 = np.random.normal(-10, 1, [100, 1])

	theta = SVGD().update(x0, GMM.grad_logp, 1000, 0.1)


	print(theta)

	import matplotlib.pyplot as plt
	import seaborn as sb

	xx = np.linspace(-10, 10, 1000).reshape((-1, 1))
	f = lambda x: np.sum(p*1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2)), axis = 1)
	plt.plot(xx, f(xx), '-')
	sb.distplot(theta, hist = True, kde = True, bins = 20)
	plt.show()
