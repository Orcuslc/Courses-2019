import numpy as np
import scipy.io as sio
from six.moves import cPickle as pickle
from tensorflow.keras.models import load_model, model_from_json

def read_data(file):
	"""Read data
	
	Arguments:
		file {str} -- path to the file
	"""
	if file.endswith(".mat"):
		return read_mat(file)
	if file.endswith(".pkl"):
		return load_data(file)

def read_mat(file):
	"""Read a .MAT file and return corresponding numpy formats 
	(i.e., transpose [N_feature, N_sample] MATLAB format into [N_sample, N_feature] Python format)
	
	Arguments:
		file {string} -- path to the mat file

	Returns:
		res {dict} -- data in Numpy format
	"""
	data = sio.loadmat(file)
	res = {}
	for (k, v) in data.items():
		if k.startswith("__"): # meta informations
			continue
		res[k] = v.T
	return res

def save_data(file, data):
	"""save data to pickle file
	
	Arguments:
		file {str} -- path to the file
		data {*} -- whatever data can be sequentialized
	"""
	with open(file, "wb") as f:
		pickle.dump(data, f)

def load_data(file):
	"""load data from pickle file
	
	Arguments:
		file {str} -- path to the file
	
	Returns:
		data -- data stored in the file
	"""
	with open(file, "rb") as f:
		data = pickle.load(f)
	return data

def save_nn(path, nn):
	"""save nn to hdf5 file
	
	Arguments:
		path {str} -- path to the file (with filename): i.e., "./xxx"
		nn {keras.models.Model} -- trained neural network
	"""
	nn.save_weights("{0}_weights.h5".format(path))
	nn_json = nn.to_json()
	with open("{0}_model.json".format(path), "w") as f:
		f.write(nn_json)
	
def load_nn(path):
	"""load model from disk
	
	Arguments:
		path {str} -- path to the file (with filename): i.e., "./xxx"
	
	Returns:
		model {keras.models.Model} -- trained neural network
	"""
	with open("{0}_model.json".format(path), "r") as f:
		model = model_from_json(f.read())
	model.load_weights("{0}_weights.h5".format(path))
	return model