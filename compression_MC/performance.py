import umbrella_run as urun
import run as run
import gdbm
import shelve
import numpy as np
import data_process as dp
def make_clickable(val):
    return '<a href="http://localhost:5000/{}" target="_blank">{}</a>'.format(val,val)

A = 15.0
##tubule params
params = {u'edges': {u'type1': {u'bending_modulus': {u'type1': 300.0,
    u'type2': float('inf'),
    u'type3': float('inf')},
   u'binding_affinity': {u'type1': 0.0, u'type2': 0.0, u'type3': 0.0},
   u'elastic_modulus': 300.0,
   u'equilibrium_angle': {u'type1': 0.18, u'type2': 0.0, u'type3': 0.0},
   u'equilibrium_length': 1.0},
  u'type2': {u'bending_modulus': {u'type1': np.float64('inf'),
    u'type2': 300.0,
    u'type3': np.float64('inf')},
   u'binding_affinity': {u'type1': 0.0, u'type2': 0.0, u'type3': 0.0},
   u'elastic_modulus': 300.0,
   u'equilibrium_angle': {u'type1': 0.0, u'type2': -0.18, u'type3': 0.0},
   u'equilibrium_length': 1.0},
  u'type3': {u'bending_modulus': {u'type1': np.float64('inf'),
    u'type2': float('inf'),
    u'type3': 300.0},
   u'binding_affinity': {u'type1': 0.0, u'type2': 0.0, u'type3': 0.0},
   u'elastic_modulus': 300.0,
   u'equilibrium_angle': {u'type1': 0.0,
    u'type2': 0.0,
    u'type3': 0.57},#-0.29000000000000004},
   u'equilibrium_length': 1.0}},
 u'excluders': {u'com': {u'overlaps': {u'com': True}, u'radius': 0.2}},
 u'global_params': {u'd_max': 0.05,
  u'dtsave': 1000,
  u'ensemble': 1,
  u'l_add': 0.5,
  u'l_fusion': 0.5,
  u'tmax': 800},
 u'subunits': {u'sub1': {u'activity': -5.0,
   u'edges': {u'edge1': u'type1', u'edge2': u'type2', u'edge3': u'type3'},
   u'excluders': {u'exc1': u'com'},
   u'insertion_rate': 0.005,
   u'number_of_rotational_configurations': 3}}}
for key in params['edges'].keys():
        params['edges'][key]['binding_affinity'][key] = A
run_params = {'parameters':params,  'data_location':'data/APS/profiler9/', 'comments':'testrun'}


import gdbm_compat
mydb = gdbm_compat.open_compat('data/APS/profiler8/Capsid', 'r')
shelf = shelve.Shelf(mydb)
capsid = shelf['6000']

for sub in capsid.subunits:
    sub.normal = np.zeros(3)
