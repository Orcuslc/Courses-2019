import numpy as np

def sqrt_l2_norm(x):
	return np.sqrt(np.sum(x**2, axis = 1, keepdims = True))

def l2_error(true, pred, relative = True):
	"""Compute error between the true vector and the predicted vector
	
	Arguments:
		true {np.ndarray} -- true values
		pred {np.ndarray} -- predicted values
	
	Keyword Arguments:
		relative {bool} -- compute absolute or relative errors (w.r.t. true values) (default: {True})
	
	Returns:
		error -- {np.ndarray}, the error for each prediction (absolute or relative)
	"""
	error = sqrt_l2_norm(true - pred)
	if relative:
		error /= sqrt_l2_norm(true)
	return error

def mse(true, pred):
	return np.mean(np.sum((true - pred)**2, axis = 1))