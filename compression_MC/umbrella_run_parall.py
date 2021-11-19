import monte_carlo as mc
import parameters as pr
from collections import OrderedDict
import Utils as ut
import shelve
import gdbm
import numpy as np    
import copy
import data_process as dp
from multiprocessing import Pool
from bson.objectid import ObjectId

import sys

#sys.path.append('/share/apps/scisoft/ANACONDA/5_2/envs/py_2_7_15/lib/python2.7/lib-dynload')
import timeit

def propagate(capsid, dt, umbrella_window, run_params, data_root):   


    params = run_params['parameters']
    dtsave = params['global_params']['dtsave']
    #random number generator. Currently pointless to pass it as argument because of the use of sets
    #mc.pr.seed)
    #number of subunits when to shoot it down
    N_subunits_approx= 300#4.0*np.pi/(3.0*3.0**0.5 * np.sin(0.5*mc.pr.equilibrium_angle)**2)
    #save the data here and create the directory tree
    path = data_root+run_params['data_location']


    RNG = np.random.RandomState()
    print "started..."+str(umbrella_window.left)

    db = gdbm.open(path+'/tmp/Capsid'+str(umbrella_window.left),'n')
    shelf = shelve.Shelf(db)

    dbu = gdbm.open(path+'/tmp/Umbrella_hist'+str(umbrella_window.left),'n')
    shelf_u = shelve.Shelf(dbu)     

    mc.update_neighbor_lists(capsid)
    #thermalize the capsid
    for t in range(10000):
        mc.attempt_move(capsid, RNG)

    connectivity_changed = True        
    

    connection_types = ['type1_fusion', 'type1_fission', 'type2_fusion', 'type2_fission']

    insertion_types = ['insertion', 'removal', 'wedge_insertion', 'wedge_removal']

    #simulation steps
    for t in xrange(1,dt):
            #print some garbage regularly and end the simulation if the capsid has grown too large
            if t%5000==0:
                n_surface_edges = len(filter(lambda x: x.is_surface_edge(), capsid.get_edges()))
                print t, len(capsid.subunits), RNG.rand(), n_surface_edges
                if n_surface_edges==0:
                    break
                n_subunits = len(capsid.subunits)
                #overgrowth
                if n_subunits>210:#2*N_subunits_approx:
                    break
                #if run_id:
                #    collection.update_one({"_id":run_id},{"$set":{'status':{'nsteps': t, 'n_subunits':n_subunits}}}, upsert = True)  


                #check full overlap from time to time
                if mc.check_overlap(capsid.subunits, capsid):
                    raise ValueError("Overlap has happened! Check your neighbor list settings!")    

            #dump the data into shelf
            if t%dtsave==0:
                shelf[str(t)]= capsid

                P_N = np.histogram(umbrella_window.values, bins=np.arange(0,200))
                shelf_u[str(umbrella_window.left)] = P_N                      
            #vertex move attempt
            for i in range(len(capsid.hubs)):
                mc.attempt_move(capsid, RNG)

            if t%50==0:    
                mc.update_neighbor_lists(capsid)


            attempt_type = np.random.choice(connection_types)                
            #fusion/fission attempts
            if (attempt_type=='type1_fusion'):
                if connectivity_changed:
                    type1_fusion_pairs = mc.get_type1_fusion_pairs(capsid)               
                if (mc.attempt_type1_fusion(capsid, RNG, type1_fusion_pairs)):
                    connectivity_changed = True
                    print "type1 fusion"


            if (attempt_type=='type2_fusion'): 
                if connectivity_changed:
                    type2_fusion_pairs = mc.get_type2_fusion_pairs(capsid)                 
                if (mc.attempt_type2_fusion(capsid, RNG, type2_fusion_pairs)):
                    connectivity_changed = True
                    print "type2 fusion"   


            if (attempt_type=='type1_fission'):     
                if connectivity_changed:
                    type1_fission_pairs = mc.get_type1_fission_pairs(capsid)                               
                if (mc.attempt_type1_fission(capsid, RNG, type1_fission_pairs)):
                    connectivity_changed = True
                    print "*** type1 fission, removable: "#, mc.get_removable_subunits(capsid)    
                    print "*** timestep ", t


            if (attempt_type=='type2_fission'):
                if connectivity_changed:
                    type2_fission_triplets = mc.get_type2_fission_triplets(capsid)              
                if (mc.attempt_type2_fission(capsid, RNG, type2_fission_triplets)):
                    connectivity_changed = True
                    print "type2 fission" 

            attempt_type = np.random.choice(insertion_types)                     

            Nsub = len(capsid.subunits)
            if (attempt_type=='insertion'):        
                
                if (Nsub<umbrella_window.right-1):
                    #insertion attempts; iterate over the subunit types and try inserting them according to their rates
                    #structure_changed is just a small optimization: don't look for surface edges again if nothing was
                    #inserted before
                    #structure_changed = True
                    for subunit_type in params['subunits']:
                        k_insertion = params['subunits'][subunit_type]['insertion_rate']
                        if connectivity_changed:
                            surface_edges = capsid.get_surface_edges()
                            #structure_changed=False
                        p_propose = k_insertion*len(surface_edges)
                        rinsert = RNG.rand()
                        if rinsert<p_propose:
                            if (mc.attempt_insertion(capsid, surface_edges, RNG, subunit_type)):
                                connectivity_changed = True
                                print "insertion ",subunit_type   
            Nsub = len(capsid.subunits)
            if (attempt_type=='removal'):
                if (Nsub>umbrella_window.left):
                    #removal attempts; get_removable_subunits() returns the subunits that are attached along a single edge (2 hubs) to the capsid
                    #since now we have multiple subunit types, this has to be filtered to the current type
                    #structure_changed = True
                    for subunit_type in params['subunits']:
                        k_insertion = params['subunits'][subunit_type]['insertion_rate'] 
                        if connectivity_changed:                   
                            removable_subunits = mc.get_removable_subunits(capsid)
                            #structure_changed=False
                        this_type_removable_subunits = filter(lambda x: x.subunit_type==subunit_type, removable_subunits)
                        p_propose = k_insertion*len(removable_subunits) #???? len(this_type_removable_subunits) ???
                        if RNG.rand()<p_propose:
                            if (mc.attempt_removal(capsid, this_type_removable_subunits, RNG)):
                                connectivity_changed = True
                                print "removal ", "subunit_type"
            Nsub = len(capsid.subunits)
            if (attempt_type=='wedge_insertion'):
                if (Nsub<umbrella_window.right-1):
                    #wedge insertion attempt
                    #structure_changed = True
                    for subunit_type in params['subunits']:
                        k_insertion = params['subunits'][subunit_type]['insertion_rate']
                        if connectivity_changed:
                            wedges = mc.get_open_wedge_pairs(capsid)
                            #structure_changed=False
                        p_propose = k_insertion*len(wedges)
                        rinsert = RNG.rand()
                        if rinsert<p_propose:
                            #print "WEDGE insertion attempting..."
                            if (mc.attempt_wedge_insertion(capsid, wedges, RNG, subunit_type)):
                                connectivity_changed = True
                                print "@@@ WEDGE insertion ",subunit_type 
            Nsub = len(capsid.subunits)
            if (attempt_type=='wedge_removal'):                                               
                if (Nsub>umbrella_window.left):
                    #wedge removal attempt
                    #structure_changed = True
                    for subunit_type in params['subunits']:
                        k_insertion = params['subunits'][subunit_type]['insertion_rate'] 
                        if connectivity_changed: 
                            #print "WEDGE REM ATT 1"                  
                            removable_wedge_subunits = mc.get_removable_wedge_subunits(capsid)
                            #structure_changed=False
                        this_type_removable_subunits = filter(lambda x: x.subunit_type==subunit_type, removable_wedge_subunits)
                        p_propose = k_insertion*len(removable_wedge_subunits) #???? len(this_type_removable_subunits) ???
                        #print "********** ", len(removable_subunits) 
                        if RNG.rand()<p_propose:
                            #print "WEDGE REM ATT"
                            if (mc.attempt_wedge_removal(capsid, this_type_removable_subunits, RNG)):
                                connectivity_changed = True
                                print "@@@ WEDGE removal ", "subunit_type", len(removable_wedge_subunits)
                                #print len(removable_subunits)

            if connectivity_changed:
                edge_fusion_hubs = mc.get_edge_fusion_hubs(capsid)
            #do these every step                                
            if (mc.attempt_edge_fusion(capsid, RNG, edge_fusion_hubs)):
                    print "EDGE FUSION"
                    connectivity_changed = True

            if connectivity_changed:
                edge_fission_hubs = mc.get_edge_fission_hubs(capsid)                     
            if (mc.attempt_edge_fission(capsid, RNG, edge_fission_hubs)):
                    connectivity_changed = True                
                    print "EDGE FISSION"



            #verify
            Nsub = len(capsid.subunits)
            if (Nsub>=umbrella_window.right) or (Nsub<umbrella_window.left):
                print Nsub, umbrella_window.right, umbrella_window.left, t
                raise ValueError("Jackass, you're out of your umbrella window")

            #update init_conf for the next window
            #middle = (umbrella_window.right+umbrella_window.left)/2                
            #if (Nsub>=middle):
            #    umbrella_window.init_conf = copy.deepcopy(capsid) 
            #    #print "xxxx ", middle, len(umbrella_window.init_conf.subunits)  
            umbrella_window.add_value(Nsub)


            #if t%umbrella_nsteps==0:
            #    #only shift window if initial configuration is available for the next window
            #    if umbrella_window.init_conf is not None:
            #        print "cccc ", Nsub, umbrella_window.right, umbrella_window.left, t, umbrella_window.init_conf
            #        print "dddd ", len(umbrella_window.init_conf.subunits), len(capsid.subunits)
            #        print "shifted..."
            #        P_N = np.histogram(umbrella_window.values, bins=np.arange(0,200))
            #        print sum(P_N[0]>0)
            #        shelf_u[str(umbrella_window.left)] = P_N    
            #        umbrella_window.shift((umbrella_window.right-umbrella_window.left)/2)
            #        capsid = copy.deepcopy(umbrella_window.init_conf)
            #        connectivity_changed = True
            #        #may need equilibration here
            #        umbrella_window.init_conf = None
            #        umbrella_window.values = []
            #    else:
            #        print "not yet... ",t
            #        #continue


            if connectivity_changed:
                surface_edges = capsid.get_surface_edges()
                removable_subunits = mc.get_removable_subunits(capsid)
                wedges = mc.get_open_wedge_pairs(capsid)
                removable_wedge_subunits = mc.get_removable_wedge_subunits(capsid)                
                type1_fusion_pairs = mc.get_type1_fusion_pairs(capsid)   
                type2_fusion_pairs = mc.get_type2_fusion_pairs(capsid)  
                type1_fission_pairs = mc.get_type1_fission_pairs(capsid)     
                type2_fission_triplets = mc.get_type2_fission_triplets(capsid)
                edge_fusion_hubs = mc.get_edge_fusion_hubs(capsid)
                edge_fission_hubs = mc.get_edge_fission_hubs(capsid)
                connectivity_changed = False
    shelf_u.close()
    shelf.close()                    
    return capsid




def run(run_params, data_root='', init_config=None, run_id=None):
    """Example run function. If _id provided, will update entries during run in db."""

    sys.setrecursionlimit(100000)


    #subset of run_params concerning the system only, but not the metadata
    params = run_params['parameters']
    #seed_capsids_id = ObjectId('5cfa910d946df6282a4b4f26')

    #print "FK", dp.get_data_location(seed_capsids_id)[1]


    #initialize global parameters
    pr.init_globals(params)


    #seed structures for umbrella sampling; should be a run_id where a seed_capsids[n_sub] shelf is located
    seed_capsids_id = params['seed_capsids']
    print seed_capsids_id 
    seed_path = dp.get_data_location(seed_capsids_id)[1]
    dbs = gdbm.open(seed_path+'/seed_capsids','r')
    seed_capsid_shelf = shelve.Shelf(dbs)
    seed_capsids = {}
    for ns in seed_capsid_shelf.keys():
        seed_capsids[int(ns)] = seed_capsid_shelf[ns]
    seed_capsid_shelf.close()





    #number of steps to run    
    tmax = params['global_params']['tmax']
    #frequency of dumping data
    dtsave = params['global_params']['dtsave']
    #random number generator. Currently pointless to pass it as argument because of the use of sets
    #mc.pr.seed)
    #number of subunits when to shoot it down
    N_subunits_approx= 300#4.0*np.pi/(3.0*3.0**0.5 * np.sin(0.5*mc.pr.equilibrium_angle)**2)
    #save the data here and create the directory tree
    path = data_root+run_params['data_location']#'data/test'#ut.create_data_folder(params)#root='/mnt/4TB_HDD/btyukodi/')
    ut.ensure_dir(path)
    ut.ensure_dir(path+'/tmp/')
 


    pool = Pool(processes=28)

    umbrella_nsteps = 12000 #number of steps in an umbrella window

    results = []
    windows_lefts = np.arange(2,32)
    for ns in windows_lefts:
        capsid = seed_capsids[ns]    

        umbrella_window = mc.UmbrellaWindow(ns,ns+4,1)
        #!!!need to update capsid parameters from params[]
        pr.update_capsid_parameters(capsid, params)
        #print ns
        
        results.append(pool.apply_async(propagate, [capsid,umbrella_nsteps, umbrella_window, run_params, data_root]))
        #pool.apply_async(f, (ns,))
        #propagate(capsid,umbrella_nsteps, ns, run_params, data_root)
    pool.close()
    pool.join()

    #force the error if something went wrong
    for res in results:
        c = res.get()

    #join results here, delete tmp
    #for parallel tempering (eventually) leave everything in tmp and with another script, create new database run entries for those and move to folders
    if run_id:
            collection = pr.get_db_collection()

    db = gdbm.open(path+'/Capsid','n')
    shelf = shelve.Shelf(db)

    dbu = gdbm.open(path+'/Umbrella_hist','n')
    shelf_u = shelve.Shelf(dbu)

    for wl in windows_lefts:
        db0 = gdbm.open(path+'/tmp/Capsid'+str(wl),'r')
        shelf0 = shelve.Shelf(db0)

        if len(shelf.keys())==0:
            index_delta=1
        else:
            index_delta = np.max(np.array(shelf.keys()).astype(int))
        for t in np.array(shelf0.keys()).astype(int):
            capsid = shelf0[str(t)]
            shelf[str(index_delta+t)] = capsid   
        shelf0.close()             
    
        dbu0 = gdbm.open(path+'/tmp/Umbrella_hist'+str(wl),'r')
        shelf_u0 = shelve.Shelf(dbu0)

        umbrella_hist = shelf_u0[str(wl)]
        shelf_u[str(wl)] = umbrella_hist
        shelf_u0.close()
    shelf.close()
    shelf_u.close()        

    #remove tmp folders here
            
##from mpl_toolkits.mplot3d import Axes3D
##import matplotlib.pyplot as plt
##def plot_capsid(capsid):
##    cx, cy, cz = 0.0, 0.0, 0.0
##    for hub in capsid.hubs:
##        if len(hub.vertices)>0:
##            cx+=hub.vertices[0].x
##            cy+=hub.vertices[0].y
##            cz+=hub.vertices[0].z
##    cx/=len(capsid.hubs)
##    cy/=len(capsid.hubs)
##    cz/=len(capsid.hubs)
##    lim =5
##    fig = plt.figure()
##    ax = fig.gca(projection='3d')
##    #sub.plot(ax)
##    #sub2.plot(ax)
##    capsid.plot(ax)
##  
##    ax.set_xlim([cx-lim,cx+lim])
##    ax.set_ylim([cy-lim,cy+lim])
##    ax.set_zlim([cz-lim,cz+lim])
##    ax.set_aspect('equal')
##    plt.show()   
##
##plot_capsid(capsid)
