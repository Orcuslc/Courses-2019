import numpy as np
import matplotlib.pyplot as plt
from functools import wraps
import tensorflow as tf
import tensorflow.keras as keras


def cut(x, domain):
	return np.where((x >= domain[0]) & (x <= domain[1]))

def numerical_solver_PML(I, V, f, c, L1, L2, dx, T, dt, cmax):
	"""
	Solve u_tt = d(c^2(x)u_x)/dx + f(x, t) on (L1, L2)x(0, T].
		I: u(t = 0)
		V: du_dt(t = 0)
		f = f(x, t)
		c = c(x)
		cmax: max(c(x)) on [L1, L2]
	"""
	# compute CFL
	CFL = cmax*dt/dx
	print("Courant Number: {0}".format(CFL))
	if CFL > 1:
		raise ValueError("CFL unsatisfied")

	Nt = int(round(T/dt))
	t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
	Nx = int(round((L2-L1)/dx))
	x = np.linspace(L1, L2, Nx+1)       # Mesh points in space
	
	u = np.zeros((Nt+1, Nx+1))   # Solution array at new time level

	# Load initial condition into u_1
	u[0, :] = I(x)

	# Special formula for first time step
# 	u[1, 1:-1] = u[0, 1:-1] \
# 				+ dt*V(x[1:-1]) \
# 				+ 0.5*(dt/dx)**2 \
# 					*(c((x[2:]+x[1:-1])/2)**2*(u[0, 2:] - u[0, 1:-1]) \
# 						- c((x[1:-1]+x[:-2])/2)**2*(u[0, 1:-1] - u[0, :-2])) \
# 				+ 0.5*dt**2*f(x[1:-1], t[0])
	u[1, 1:-1] = u[0, 1:-1] \
				+ dt*V(x[1:-1]) \
				+ 0.5*(dt/dx)**2 \
					*((c(x[2:])**2+c(x[1:-1])**2)/2*(u[0, 2:] - u[0, 1:-1]) \
						- (c(x[1:-1])**2+c(x[:-2])**2)/2*(u[0, 1:-1] - u[0, :-2])) \
				+ 0.5*dt**2*f(x[1:-1], t[0])
	
	u[1, 0] = 0;  u[1, Nx] = 0
	
	for n in range(1, Nt):
		# Update all inner points at time t[n+1]
# 		u[n+1, 1:-1] = -u[n-1, 1:-1] + 2*u[n, 1:-1] \
# 						+ (dt/dx)**2 \
# 							*(c((x[2:]+x[1:-1])/2)**2*(u[n, 2:] - u[n, 1:-1]) \
# 								- c((x[1:-1]+x[:-2])/2)**2*(u[n, 1:-1] - u[n, :-2])) \
# 						+ dt**2*f(x[1:-1], t[n])
		u[n+1, 1:-1] = -u[n-1, 1:-1] + 2*u[n, 1:-1] \
						+ (dt/dx)**2 \
							*((c(x[2:])**2+c(x[1:-1])**2)/2*(u[n, 2:] - u[n, 1:-1]) \
								- (c(x[1:-1])**2+c(x[:-2])**2)/2*(u[n, 1:-1] - u[n, :-2])) \
						+ dt**2*f(x[1:-1], t[n])

		# Insert boundary conditions
		u[n+1, 0] = 0;  u[n+1, Nx] = 0

	return u, x, t
	

def numerical_solver(I, V, f, c, L1, L2, dt, C, T):
	"""Solve u_tt=c^2*u_xx + f on (L1,L2)x(0,T].
		I: u(t = 0)
		V: du_dt(t = 0)
	 	C: Courant number
	"""
	Nt = int(round(T/dt))
	t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
	dx = dt*c/float(C)
	Nx = int(round((L2-L1)/dx))
	x = np.linspace(L1, L2, Nx+1)       # Mesh points in space
	C2 = C**2                      # Help variable in the scheme

	u = np.zeros((Nt+1, Nx+1))   # Solution array at new time level

	# Load initial condition into u_1
	u[0, :] = I(x)

	# Special formula for first time step
	u[1, 1:-1] = u[0, 1:-1] + dt*V(x[1:-1]) + 0.5*C2*(u[0, :-2] - 2*u[0, 1:-1] + u[0, 2:]) + 0.5*dt**2*f(x[1:-1], t[0])
	u[1, 0] = 0;  u[1, Nx] = 0
	
	for n in range(1, Nt):
		# Update all inner points at time t[n+1]
		u[n+1, 1:-1] = -u[n-1, 1:-1] + 2*u[n, 1:-1] + C2*(u[n, :-2] - 2*u[n, 1:-1] + u[n, 2:]) + dt**2*f(x[1:-1], t[n])

		# Insert boundary conditions
		u[n+1, 0] = 0;  u[n+1, Nx] = 0

	return u, x, t

def dAlembert(I, c, x, t):
	"""
	Solve the Cauchy problem 
		u_tt = c^2 u_xx
	with initial condition
		u(x, t = 0) = I(x), u_t(x, t = 0) = 0.
	"""
	func_true = lambda x, t: 1/2*I(x.flatten() + c*t.flatten()) + 1/2*I(x.flatten() - c*t.flatten())
	u = np.zeros((len(t), len(x)))
	for i, ti in enumerate(t):
		u[i, :] = np.squeeze(func_true(x, np.ones_like(x)*ti))
	return u

def transform(x, L1, L2):
	return x*(L2-L1) + L1

def sqrt_l2_norm(x):
	return np.sqrt(np.sum(x**2, axis = 0, keepdims = True))

def l2_error(true, pred, relative = True):
	error = sqrt_l2_norm(true - pred)
	if relative:
		error /= sqrt_l2_norm(true)
	return error

def plot(model, x, t_indices, u_true, dt, col_index = 0, ncol = 4):
	if len(t_indices) < ncol:
		f, ax = plt.subplots(1, len(t_indices), figsize = (len(t_indices)*5, 5))
	else:
		f, ax = plt.subplots((len(t_indices)-1)//ncol + 1, ncol, figsize = (20, (20 // ncol)*5))
	f2, ax2 = plt.subplots(1, 2, figsize = (10, 5))
	
	t = []
	abs_err = []
	rel_err = []

	for (i, ti) in enumerate(t_indices):
		t.append(ti*dt)
		
		u_pred = model(np.hstack([x, np.ones_like(x)*t[-1]]))[:, col_index:col_index+1]
		u_true_ = u_true[ti, :].reshape((-1, 1))
		abs_err.append(l2_error(u_true_, u_pred, False).flatten()[0])
		rel_err.append(l2_error(u_true_, u_pred, True).flatten()[0])
		ax[i // ncol][i % ncol].plot(x, u_true_, label = "true")
		ax[i // ncol][i % ncol].plot(x, u_pred, label = "pred")
		ax[i // ncol][i % ncol].text(-1.8, -0.35, "abs_err: {0:.4f}\nrel_err: {1:.4f}".format(abs_err[-1], rel_err[-1]))
		ax[i // ncol][i % ncol].legend()
		ax[i // ncol][i % ncol].grid()
		ax[i // ncol][i % ncol].set_ylim([-0.5, 0.6])
		ax[i // ncol][i % ncol].set_title("t = {0:.2f}".format(t[-1]))
	
	ax2[0].plot(t, abs_err)
	ax2[0].set_title("abs. err.")
	ax2[0].set_xlabel("time")
	ax2[0].grid()
	ax2[1].plot(t, rel_err)
	ax2[1].set_title("rel. err.")
	ax2[1].set_xlabel("time")
	ax2[1].grid()
	plt.show()

def cast_to_tf_constant(dtype):
	def decorate(func):
		@wraps(func)
		def cast(*args):
			newargs = []
			for arg in args:
				newargs.append(tf.constant(arg, dtype = dtype))
			return func(*newargs)
		return cast
	return decorate