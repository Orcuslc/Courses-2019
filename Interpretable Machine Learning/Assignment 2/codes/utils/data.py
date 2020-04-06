import numpy as np

def tensor_grid(x):
	"""build tensor grid for multiple parameters
	
	Arguments:
		x {tuple or list of np.array} -- parameters
	
	Returns:
		grid {np.ndarray} -- tensor grids

	Example:
		>>> tensor_grid(([1, 2], [3, 4], [5, 6, 7]))
		>>> np.array([[1, 3, 5],
					[1, 3, 6],
					[1, 3, 7],
					[1, 4, 5],
					[1, 4, 6],
					[1, 4, 7],
					[2, 3, 5],
					[2, 3, 6],
					[2, 3, 7],
					[2, 4, 5],
					[2, 4, 6],
					[2, 4, 7]])
	"""
	return np.vstack(np.meshgrid(*x, indexing = 'ij')).reshape((len(x), -1)).T

def train_valid_index_split(all_index, train_size = None, valid_split = 0.3):
	"""Split train index and valid index, sampled from a large candidate set
	
	Arguments:
		all_index {np.array or int} -- Index of all training data; if an `int`, all_index = np.arange(all_index)
	
	Keyword Arguments:
		train_size {int} -- size of training set, which includes both training and validation sets (Default: {None} -- use the whole candidate set as training set)
		valid_split {float between 0 and 1} -- Fraction of validation data within training set (Default: {0.3})
	"""
	all_index = np.arange(all_index) if isinstance(all_index, int) else np.array(all_index)
	train_size = len(all_index) if train_size is None else train_size
	train_index_ = np.random.choice(all_index, train_size, replace = False)
	train_index, valid_index = np.split(train_index_, [int(train_size*(1-valid_split))])
	return train_index, valid_index

def train_valid_index_split_two_stage(all_index, train_size_1 = None, train_size_2 = None, valid_split = 0.3):
	"""Split train index and valid index for two-stage training; the training (validation) index for stage 2 are all contained in the training (validation) index for stage 1
	
	Arguments:
		all_index {np.array or int} -- Index of all training data; if an `int`, all_index = np.arange(all_index)
	
	Keyword Arguments:
		train_size_1 {int} -- size of training set in stage 1 (default: {None} -- use the whole set)
		train_size_2 {int} -- size of training set in stage 2 (default: {None} -- use the whole set)
		valid_split {float between 0 and 1} -- Fraction of validation data within training set (default: {0.3})
	"""
	all_index = np.arange(all_index) if isinstance(all_index, int) else np.array(all_index)

	train_size_2 = len(all_index) if train_size_2 is None else train_size_2
	train_index_2_ = np.random.choice(all_index, train_size_2, replace = False)
	train_index_2, valid_index_2 = np.split(train_index_2_, [int(train_size_2*(1-valid_split))])

	all_index = np.setdiff1d(all_index, train_index_2)
	train_index_1_ = np.random.choice(all_index, train_size_1-train_size_2, replace = False)
	train_index_1, valid_index_1 = np.split(train_index_1_, [int((train_size_1-train_size_2)*(1-valid_split))])
	train_index_1 = np.hstack([train_index_1, train_index_2])
	valid_index_1 = np.hstack([valid_index_1, valid_index_2])
	return train_index_1, valid_index_1, train_index_2, valid_index_2


def reshape_to_trajectory(data, trajectory_length):
	"""reshape dataset to trajectories
	
	Arguments:
		data {np.ndarray} -- snapshot matrix, of shape [n_sample, n_feature]
		trajectory_length {int} -- length of each trajectory, i.e., number of timesteps in each trajectory
	"""
	return data.reshape((-1, trajectory_length*data.shape[1]))