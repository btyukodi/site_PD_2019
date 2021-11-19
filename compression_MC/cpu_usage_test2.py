def test():

	import numpy as np 
	from scipy.linalg import expm

	axis = np.array([5.0, 8.0, -2.0])
	axis = axis/(np.sum(axis**2)**0.5)

	theta = 0.001


	for t in range(100000):
		cx = np.cross(np.eye(3), axis*theta)
		M0 = expm(cx)
		cx = np.dot(M0, cx) 

	print cx

