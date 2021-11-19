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


def propagate(capsid, kT, t_init, t_final, umbrella_window, run_params, RNG, path):


    #path = run_params['data_location']

    db = gdbm.open(path+'/tmp/Capsid_kT'+("%.9f" % kT),'w')
    shelf = shelve.Shelf(db)
    

    dbu = gdbm.open(path+'/tmp/Umbrella_hist_kT'+("%.9f" % kT),'w')
    shelf_u = shelve.Shelf(dbu)
      

    params = run_params['parameters']

    dtsave = params['global_params']['dtsave']

    connectivity_changed = True        
    connection_types = ['type1_fusion', 'type1_fission', 'type2_fusion', 'type2_fission']
    insertion_types = ['insertion', 'removal', 'wedge_insertion', 'wedge_removal']


    mc.update_neighbor_lists(capsid)

    #simulation steps --> keep only this in propagate, from t1 to t2. Window shifting and swapping happens outside, in series
    for t in xrange(t_init,t_final):
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
                        p_propose = k_insertion*len(this_type_removable_subunits) #???? len(this_type_removable_subunits) ???
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
                        p_propose = k_insertion*len(this_type_removable_subunits)
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
            middle = (umbrella_window.right+umbrella_window.left)/2                
            if (Nsub>=middle) & (umbrella_window.init_conf is None):
                umbrella_window.init_conf = copy.deepcopy(capsid) 
                #for i in range(1000):
                #    print "ADDED"
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
    return capsid, umbrella_window