from pymongo import MongoClient
import json
from bson.objectid import ObjectId 
import warnings
import config

class DummyCollection:
	"""Needs in case runs are loaded from file. It's a workaround until DB is working on hpcc"""
	def __init__(self, filename):
		data = []
		for line in open(filename, "r"):
			data.append(json.loads(line))
		for d in data:
			d['_id'] = ObjectId(d['_id']['$oid'])
			try:
				d['init_config']['id'] = ObjectId(d['init_config']['id']['$oid'])
			except:
				pass
		self.data = data
		
	def find(self, q):
		key = q.keys()[0]
		#return filter(lambda x: x[key]==q[key], self.data)
		return filter(lambda x: False if key not in x.keys() else x[key]==q[key], self.data)
	def find_one(self, q):
		return self.find(q)[0]
	def update_one(self, *args, **kwargs):
		pass


def get_db_collection():
	try:
		client = config.client()
		client.server_info()
	except:
		json_file = config.json_file#'/home/btyukodi/runs.staged-runs.json'
		warnings.warn("DB connection failed. Using Json dump "+json_file)
		return DummyCollection(json_file)
	db = client['runs']
	collection = db['staged-runs']
	return collection

#def get_data_root_config():
#	try:
#		client = mongo_config.client()
#		client.server_info()
#	except:
#		json_file = '/home/btyukodi/runs.data_roots.json'
#		warnings.warn("DB connection failed. Using Json dump "+json_file)
#		return DummyCollection(json_file)
#
#
#	db = client['runs']
#	collection = db['data_roots']
#	return collection


def init_globals(parameters):
	"""Should be called before every run. Sets and checks global variables."""
	global params, l_fusion, d_max, l_add
	params = parameters

	l_fusion=params['global_params']['l_fusion']
	d_max = params['global_params']['d_max']
	l_add = params['global_params']['l_add']                      


	#check excluder overlap symmetry
	Exclude_Lookup = {}
	for exc1 in params['excluders']:
		for exc2 in params['excluders'][exc1]['overlaps']:
			Exclude_Lookup[(exc1, exc2)] = params['excluders'][exc1]['overlaps'][exc2]                            


	for key in Exclude_Lookup.keys():
		if Exclude_Lookup[key[::-1]] <> Exclude_Lookup[key]:
			raise ValueError('Excluder interactions have to be symmetric')                        

	#check the same for edge angle, modulus, affinity symmetry
	for prop in ['equilibrium_angle', 'bending_modulus', 'binding_affinity']:
		Edge_Lookup = {}
		for e1 in params['edges']:
			for e2 in params['edges'][e1][prop]:
				Edge_Lookup[(e1, e2)] = params['edges'][e1][prop][e2]
		for key in Edge_Lookup.keys():
			if Edge_Lookup[key[::-1]] <> Edge_Lookup[key]:
				raise ValueError('Edge interactions have to be symmetric')  


def load_parameters_from_db(_id):
	collection = get_db_collection()
	return collection.find_one({"_id":_id})

def update_capsid_parameters(capsid, new_params):
	"""Given a capsid, its parameters can be updated according to new_params"""
	subunits = capsid.subunits
	edges = capsid.get_edges()

	for subunit in capsid.subunits:
		sub_type = subunit.subunit_type
		subunit.activity = new_params['subunits'][sub_type]['activity']
		subunit.number_of_rotational_configurations = new_params['subunits'][sub_type]['number_of_rotational_configurations']
		for excluder in subunit.excluders:
			exc_type = excluder.excluder_type
			excluder.overlaps = new_params['excluders'][exc_type]['overlaps']
			excluder.R = new_params['excluders'][exc_type]['radius']
	for edge in edges:
		edge_type = edge.edge_type
		for key in new_params['edges'][edge_type]:
			setattr(edge, key, new_params['edges'][edge_type][key])
