
# cython: profile=False

def generate_seeds(run_params, data_root='', init_config=None, run_id=None):
    """Example run function. If _id provided, will update entries during run in db."""

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
 
    #subset of run_params concerning the system only, but not the metadata
    params = run_params['parameters']

    #boost activities here? or binding strengths?
    params['subunits']['sub1']['activity'] = 10.0


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
    db = gdbm.open(path+'/seed_capsids','n')
    shelf = shelve.Shelf(db)

    #start with a single subunit capsid
    if not init_config:
        capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])#mc.create_n_gon(5)#mc.create_hex_sheet(1,1)
    else:
        capsid = init_config    
        #!!!need to update capsid parameters from params[]
        pr.update_capsid_parameters(capsid, params)


    mc.update_neighbor_lists(capsid)
    #thermalize the capsid
    for t in range(10000):
        mc.attempt_move(capsid, RNG)

    connectivity_changed = True

    n_subunits_prev = 0
    #simulation steps
    for t in xrange(tmax):
       
            #print some garbage regularly and end the simulation if the capsid has grown too large
            if t%5000==0:
                n_surface_edges = len(filter(lambda x: x.is_surface_edge(), capsid.get_edges()))
                print t, len(capsid.subunits), RNG.rand(), n_surface_edges
                if n_surface_edges==0:
                    break
                n_subunits = len(capsid.subunits)
        		#overgrowth
                if n_subunits>100:
                    break
                if run_id:
                    collection.update_one({"_id":run_id},{"$set":{'status':{'nsteps': t, 'n_subunits':n_subunits}}}, upsert = True)    

                #check full overlap from time to time
                if mc.check_overlap(capsid.subunits, capsid):
                    raise ValueError("Overlap has happened! Check your neighbor list settings!")    
            #dump the data into shelf
            n_subunits = len(capsid.subunits)
            if n_subunits > n_subunits_prev:
                shelf[str(n_subunits)]= capsid
                n_subunits_prev = n_subunits
            #vertex move attempt
            for i in range(len(capsid.hubs)):
                mc.attempt_move(capsid, RNG)
            #print "----- MOVED  "
            if t%50==0:    
                mc.update_neighbor_lists(capsid)

            #fusion/fission attempts
            if connectivity_changed:
                type1_fusion_pairs = mc.get_type1_fusion_pairs(capsid)
            if (mc.attempt_type1_fusion(capsid, RNG, type1_fusion_pairs)):
                connectivity_changed = True
                print "type1 fusion"

            if connectivity_changed:
                type2_fusion_pairs = mc.get_type2_fusion_pairs(capsid)                
            if (mc.attempt_type2_fusion(capsid, RNG, type2_fusion_pairs)):
                connectivity_changed = True                
                print "type2 fusion"   

            if connectivity_changed:
                type1_fission_pairs = mc.get_type1_fission_pairs(capsid)    
            if (mc.attempt_type1_fission(capsid, RNG, type1_fission_pairs)):
                connectivity_changed = True                
                print "*** type1 fission"
                print "*** timestep ", t

            if connectivity_changed:
                type2_fission_triplets = mc.get_type2_fission_triplets(capsid)
            if (mc.attempt_type2_fission(capsid, RNG, type2_fission_triplets)):
                connectivity_changed = True                
                print "type2 fission" 

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
                        #structure_changed = True
                        connectivity_changed = True                        
                        print "insertion ",subunit_type       

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
                        #structure_changed = True
                        connectivity_changed = True
                        print "removal ", "subunit_type"


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
                        #structure_changed = True
                        connectivity_changed = True
                        print "@@@ WEDGE insertion ",subunit_type                

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
                        #structure_changed = True
                        connectivity_changed = True                        
                        print "@@@ WEDGE removal ", "subunit_type", len(removable_wedge_subunits)
                        #print len(removable_subunits)



            if connectivity_changed:
                edge_fusion_hubs = mc.get_edge_fusion_hubs(capsid)
            if (mc.attempt_edge_fusion(capsid, RNG, edge_fusion_hubs)):
                print "EDGE FUSION"
                connectivity_changed = True

            if connectivity_changed:
                edge_fission_hubs = mc.get_edge_fission_hubs(capsid)        
            if (mc.attempt_edge_fission(capsid, RNG, edge_fission_hubs)):
                print "EDGE FISSION"
                connectivity_changed = True
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
            if n_subunits>100:
                break 
    shelf.close()
    return capsid



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
