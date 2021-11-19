from pymongo import MongoClient
def client():
	return MongoClient('xxx', 27017, username='xxx', password='xxx', 
			authSource='xxx', authMechanism='SCRAM-SHA-1')

data_root = '/workstation/assembly_MC/data/'#'/home/btyukodi/assembly_MC/data/'#'/scratch/btyukodi/data/'#'/home/btyukodi/assembly_MC/data/'
json_file= '/home/btyukodi/runs.staged-runs.json'

#submit_job and submit_process_job could also be defined here, per host
