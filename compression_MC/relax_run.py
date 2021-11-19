import numpy as np
import pandas as pd
import data_process as dp
import parameters as pr
from pandas.io.json import json_normalize
import re
import sys
from bson.objectid import ObjectId
#sys.path.append('/share/apps/scisoft/ANACONDA/5_2/envs/py_2_7_15/lib/python2.7/lib-dynload')
import timeit
sys.setrecursionlimit(100000)


import monte_carlo as mc
import parameters as pr
from collections import OrderedDict
import Utils as ut
import shelve
import gdbm
import numpy as np   


def relax_run(run_params, data_root='', init_config=None, run_id=None):
    """Example run function. If _id provided, will update entries during run in db."""


 
    #subset of run_params concerning the system only, but not the metadata
    params = run_params['parameters']
    #initialize global parameters
    pr.init_globals(params)


    #number of steps to run    
    tmax = params['global_params']['tmax']
    #frequency of dumping data
    dtsave = params['global_params']['dtsave']
    #random number generator. Currently pointless to pass it as argument because of the use of sets
    RNG = np.random.RandomState()#mc.pr.seed)
    #number of subunits when to shoot it down
    N_subunits_approx= 200#4.0*np.pi/(3.0*3.0**0.5 * np.sin(0.5*mc.pr.equilibrium_angle)**2)
    #save the data here and create the directory tree
    path = data_root+run_params['data_location']#'data/test'#ut.create_data_folder(params)#root='/mnt/4TB_HDD/btyukodi/')
    ut.ensure_dir(path)
    #gdbm database/shelf for the capsid
    db = gdbm.open(path+'/Capsid','n')
    shelf = shelve.Shelf(db)

    #start with a single subunit capsid
    if not init_config:
        capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])#mc.create_n_gon(5)#mc.create_hex_sheet(1,1)
    else:
        capsid = init_config    
        #!!!need to update capsid parameters from params[]
        pr.update_capsid_parameters(capsid, params)


    if run_id:
        collection = pr.get_db_collection()


    #thermalize the capsid
    for t in range(10000):
        mc.attempt_move(capsid, RNG)

    capsid.kT=1.5
    #simulation steps
    for t in range(tmax):
            #if t%500==0:
                
                #params['edges']['type3']['equilibrium_angle']['type3']-=0.005
                #print "WTF", params['edges']['type3']['equilibrium_angle']['type3']
                #pr.update_capsid_parameters(capsid, params)        
                #print [(e.edge_type, e.equilibrium_angle) for e in capsid.subunits[0].edges]
        
        
            #print some garbage regularly and end the simulation if the capsid has grown too large
            if t%5000==0:
                n_surface_edges = len(filter(lambda x: x.is_surface_edge(), capsid.get_edges()))
                print t, len(capsid.subunits), RNG.rand(), n_surface_edges
           #dump the data into shelf
            if t%dtsave==0:
                shelf[str(t)]= capsid
                print "dump"
                print "kT", capsid.kT
            #vertex move attempt
            for i in range(len(capsid.hubs)):
                mc.attempt_move(capsid, RNG)

            if t>2000:    
                capsid.kT-=0.0001
                if capsid.kT<0.0:
                    break    
 

    shelf.close()
    return capsid


#run_params = pr.get_db_collection().find_one({"_id":ObjectId("5c2e5608946df694be826db5")})    
#run_params['data_location'] = 'cylinder_blow_test/'
#run_params['data_root'] = 'data/'
#run_params['parameters']['global_params']['dtsave'] = 50
#run_params['parameters']['global_params']['tmax']= 600

#init_config = dp.get_single_capsid_snapshot(run_id=ObjectId("5c2e5608946df694be826db5"), timestep=-1)

#relax_run(run_params, data_root='data/', init_config=init_config)