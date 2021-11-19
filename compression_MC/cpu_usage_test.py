def run(run_params, data_root='', init_config=None, run_id=None):
    """Example run function. If _id provided, will update entries during run in db."""
    print "baa"
    import sys

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
    import copy
    import geometry as g
    import helpers as h
    import energy as e
 
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
    N_subunits_approx= 300#4.0*np.pi/(3.0*3.0**0.5 * np.sin(0.5*mc.pr.equilibrium_angle)**2)
    #save the data here and create the directory tree
    path = data_root+run_params['data_location']#'data/test'#ut.create_data_folder(params)#root='/mnt/4TB_HDD/btyukodi/')
    ut.ensure_dir(path)
    #gdbm database/shelf for the capsid
    db = gdbm.open(path+'/Capsid','n')
    shelf = shelve.Shelf(db)
    print "PATH: ",path

    dbu = gdbm.open(path+'/Umbrella_hist','n')
    shelf_u = shelve.Shelf(dbu)

    #start with a single subunit capsid
    if not init_config:
        capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])#mc.create_n_gon(5)#mc.create_hex_sheet(1,1)
    else:
        capsid = init_config    
        #!!!need to update capsid parameters from params[]
        pr.update_capsid_parameters(capsid, params)


    capsid.kT = params['global_params']['kT']
    if run_id:
        collection = pr.get_db_collection()

    mc.update_neighbor_lists(capsid)
    #thermalize the capsid
    for t in range(10000):
        mc.attempt_move(capsid, RNG)

    connectivity_changed = True        


    umbrella_window = mc.UmbrellaWindow(0,111,1)
    umbrella_nsteps = 7500000*30 #number of steps in an umbrella window

    connection_types = ['type1_fusion', 'type1_fission', 'type2_fusion', 'type2_fission']

    insertion_types = ['insertion', 'removal', 'wedge_insertion', 'wedge_removal']


    for t in range(1000):

            surface_edges = capsid.get_surface_edges()
            mc.attempt_insertion(capsid, surface_edges, RNG, 'sub1')

    print "n=", len(capsid.subunits)           

    #simulation steps
    for t in range(1,tmax):
            #print some garbage regularly and end the simulation if the capsid has grown too large
            if t%500000==0:
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
                #print 'nsub = ', len(capsid.subunits)
            #vertex move attempt
            #for i in range(len(capsid.hubs)):
            #    mc.attempt_move(capsid, RNG)

            #if t%50==0:    
            mc.update_neighbor_lists(capsid)
             


            attempt_type = np.random.choice(insertion_types)                     

            Nsub = len(capsid.subunits)

            surface_edges = capsid.get_surface_edges()


            surface_edges = capsid.get_surface_edges()
            mc.attempt_insertion(capsid, surface_edges, RNG, 'sub1')

            for kk in range(1000):
                #h.update_excluders_position(capsid.subunits)
            #mc.attempt_insertion(capsid, surface_edges, RNG, 'sub1')
                #Nn = g.get_triangle_normal(capsid.subunits[0])
                #O = e.check_overlap(capsid.subunits, capsid)
                #E = e.elastic_energy(capsid.subunits)
                #Nnn = g.cross3D([3,6.0,8.0], [5,1.0,8.0])
                z = np.dot([3,6.0,8.0], [5,1.0,8.0])
                #g.rotate_vertex(capsid.subunits[0].vertices[0],,theta)
##            if (attempt_type=='insertion'):        
##                
##                if (Nsub<umbrella_window.right-1):
##                    #insertion attempts; iterate over the subunit types and try inserting them according to their rates
##                    #structure_changed is just a small optimization: don't look for surface edges again if nothing was
##                    #inserted before
##                    #structure_changed = True
##                    for subunit_type in params['subunits']:
##                        k_insertion = params['subunits'][subunit_type]['insertion_rate']
##                        if connectivity_changed:
##                            surface_edges = capsid.get_surface_edges()
##                            #structure_changed=False
##                        p_propose = k_insertion*len(surface_edges)
##                        rinsert = RNG.rand()
##                        if rinsert<p_propose:
##                            if (mc.attempt_insertion(capsid, surface_edges, RNG, subunit_type)):
##                                connectivity_changed = True
                        ##        print "insertion ",subunit_type   
##            Nsub = len(capsid.subunits)
##            if (attempt_type=='removal'):
##                if (Nsub>umbrella_window.left):
##                    #removal attempts; get_removable_subunits() returns the subunits that are attached along a single edge (2 hubs) to the capsid
##                    #since now we have multiple subunit types, this has to be filtered to the current type
##                    #structure_changed = True
##                    for subunit_type in params['subunits']:
##                        k_insertion = params['subunits'][subunit_type]['insertion_rate'] 
##                        if connectivity_changed:                   
##                            removable_subunits = mc.get_removable_subunits(capsid)
##                            #structure_changed=False
##                        this_type_removable_subunits = filter(lambda x: x.subunit_type==subunit_type, removable_subunits)
##                        p_propose = k_insertion*len(removable_subunits) #???? len(this_type_removable_subunits) ???
##                        if RNG.rand()<p_propose:
##                            if (mc.attempt_removal(capsid, this_type_removable_subunits, RNG)):
##                                connectivity_changed = True
##                                print "removal ", "subunit_type"



            if connectivity_changed:
                surface_edges = capsid.get_surface_edges()
                #removable_subunits = mc.get_removable_subunits(capsid)

                connectivity_changed = False


    shelf_u.close()
    shelf.close()
    return capsid



import parameters as pr
from bson.objectid import ObjectId
params = pr.load_parameters_from_db(ObjectId('5cf81833946df6e286164bf5'))
run(params, data_root="test/cpu/")