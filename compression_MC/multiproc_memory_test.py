from multiprocessing import Pool
import numpy as np

def worker(x):
	l = 10000000
	y = (np.arange(l)+x)/float(l)
	return sum(y*0.5)


pool = Pool(processes=10)


xval = np.arange(10)
for t in range(10000):
	results = []
	for x in xval:
		#print x
		results.append(pool.apply_async(worker, [x]))
		#results.append(worker(x))
	for ix, res in enumerate(results):
		xval[ix] = res.get()

pool.close()
pool.join()