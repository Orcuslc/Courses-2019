import numpy as np


def sigma(x, w): # logistic sigmoid
	return 1/(1+np.exp(-x.dot(w)))

def expand_constant(X):
	return np.hstack([np.ones((X.shape[0], 1)), X])

def L(x, y, w): # loss
	si = sigma(x, w)
	return -np.mean(y*np.log(si) + (1-y)*np.log(1-si)), si

def dL(x, y, si, w): # gradient
	lam = 0
	return np.mean((si - y)*x, axis = 0).reshape((-1, 1)) - 2*lam*w

class Logistic_Regression:
	def __init__(self):
		self.w = None

	def fit(self, solver, X, y, w0, *, test_X = None, test_y = None, debug = False, **kwargs):
		if solver == "GD":
			self.w, hist = gradient_descent_for_logistic_regression(X, y, w0, test_X = test_X, test_y = test_y, debug = debug, **kwargs)
		elif solver == "SA":
			self.w, hist = simulated_annealing_for_logistic_regression(X, y, w0, test_X = test_X, test_y = test_y, debug = debug, **kwargs)
		return self.w, hist

	def predict(self, X):
		return sigma(expand_constant(X), self.w)		

def gradient_descent_for_logistic_regression(X, y, w0, *, num_iter, lr, test_X = None, test_y = None, debug = False, print_every = 1000):
	"""Gradient descent for logistic regression
	
	Arguments:
		X {np.ndarray} -- feature matrix, of shape [n_sample, n_feature]
		y {np.ndarray} -- label vector, of shape [n_sample, 1]
		w0 {np.ndarray} -- initial guess of weights, w = [w_0, w_1, ..., w_k]
		num_iter {int} -- number of iterations	
		lr {float} -- learning rate
		test_X, test_y {np.ndarray} -- test data (default: None)
	"""
	X = expand_constant(X) # stack X with a constant vector
	test_X = expand_constant(test_X)

	hist = {"train": [], "test": []}
	for i in range(num_iter):
		loss, si = L(X, y, w0) # loss
		hist["train"].append(loss)
		if test_X is not None and test_y is not None: # test loss
			test_loss, _ = L(test_X, test_y, w0)
			hist["test"].append(test_loss)

		w = w0 - lr*dL(X, y, si, w0) # GD
		if debug:
			if i % print_every == 0:
				print("Iteration: {0}, training loss: {1}".format(i, loss))
		w0 = w

	return w, hist

def neighbour(w, scale): # find one random neighbour
	# make sure delta in [-1, 1]
	delta = np.random.random(w.shape)*2-1
	return w + delta*scale

def acceptance(L0, L1, temperature):
	"""Yields the probability of acceptance
	
	Arguments:
		L0 {float} -- current loss
		L1 {float} -- new loss
		temperature {float} -- current temperature
	"""
	return 1 if L1 < L0 else np.exp(-(L1-L0)/temperature)

def linearly_decreasing_temperature(i, num_iter, T0):
	"""Yields current temperature (with linear decrease)
	
	Arguments:
		T0 {float} -- initial temperature
		i {int} -- current iteration
		num_iter {int} -- total iterations
	"""
	return T0*(1-i/num_iter)

def simulated_annealing_for_logistic_regression(X, y, w0, *, num_iter, temperature, neighbour, acceptance, test_X = None, test_y = None, debug = False, print_every = 1000):
	"""Simulated_annealing for logistic regression
	
	Arguments:
		X {np.ndarray} -- feature matrix, of shape [n_sample, n_feature]
		y {np.ndarray} -- label vector, of shape [n_sample, 1]
		w0 {np.ndarray} -- initial guess of weights, w = [w_0, w_1, ..., w_k].T
		num_iter {int} -- number of iterations	
		temperature {function} -- temperature (of i, num_iter)
		neighbour {function} -- random neighbour (of w)
		acceptance {function} -- acceptance of the new state (of L0, L1, T)
		test_X, test_y {np.ndarray} -- test data (default: None)
	"""
	X = expand_constant(X) # stack X with a constant vector
	test_X = expand_constant(test_X)

	hist = {"train": [], "test": []}
	loss, _ = L(X, y, w0)
	hist["train"].append(loss)
	if test_X is not None and test_y is not None:
		test_loss, _ = L(test_X, test_y, w0)
		hist["test"].append(loss)

	for i in range(num_iter):
		T = temperature(i, num_iter)
		w = neighbour(w0)
		new_loss, _ = L(X, y, w)
		acc = 0 # flag of acceptance
		if acceptance(loss, new_loss, T) > np.random.random():
			w0 = w
			loss = new_loss
			hist["train"].append(new_loss)
			
			if test_X is not None and test_y is not None: # test loss
				test_loss, _ = L(test_X, test_y, w0)
				hist["test"].append(test_loss)
			acc = 1

		if debug:
			if i % print_every == 0:
				print("Iteration: {0}, training loss: {1}, {2} new state".format(i, new_loss, "accept" if acc else "reject"))
	
	return w, hist
				
