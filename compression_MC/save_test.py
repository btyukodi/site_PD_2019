    for wedge_pair in new_wedge_pairs:
            next_capsid = fuse_copy(next_capsid, wedge_pair)
            relax(next_capsid, RNG)
            configurations.append(next_capsid)
            generate_next_state(next_capsid)
            #relax(next_capsid, RNG)
        #get_open_wedge_pairs(next_capsid)
        #for type1_fusion_pairs in ...
            #next_capsid = fuse_copy()
            #generate_next_state
generate_next_state(capsid)
886/272: len(configurations)
886/273:
pent_conf = [conf for conf in configurations if len(conf.subunits)==5]
pent_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==5 ]
886/274: dp.plot_capsid(pent_conf2[0])
886/275: dp.plot_capsid(pent_conf[0])
886/276:
elastic_modulus=50.0
kT = 1.0
bending_modulus = 300.0
theta = 0.1

#!!! NO EXCLUDERS FOR CONNECTIVITY

params= {'global_params':{'ensemble':2, 'l_fusion':1.0, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':30000000, 'dtsave':2000, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                             'excluders':{'exc1':'com'},
                             'activity':-3.0,
                             'insertion_rate':0.1,
                             'number_of_rotational_configurations':1}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':theta},
                            'bending_modulus':{'type1':bending_modulus},
                            'binding_affinity':{'type1':5.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },                
                            
                 },
         'excluders':{'com':{'overlaps':{'com':False}, 'radius':0.0}
                     }
                  }
pr.init_globals(params)
886/277:
##GENERATE AND EQUILIBRATE WITHOUT EXCLUDERS. ADD EXCLUDERS AFTERWARDS, REMOVE CONFIGURATIONS WHERE THEY OVERLAP
#assumes each connectivity has a minimum and not a glassy landscape
#compute 3 moments of inertia to identify degenerate structures
#includes l_fusion for fusion, otherwise runs in geometry/relaxation problems
#won't generate too high energy confs due to l_fusion
capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])
configurations = []
def generate_next_state(capsid):
    if len(capsid.subunits)>4:
        return
    print capsid, len(capsid.subunits)
    surface_edges = capsid.get_surface_edges()
    for surface_edge in surface_edges:
        next_capsid = insert_copy(capsid, surface_edge, 'sub1')
        #should equilibrate capsid here to drive it towards the global minimum of the potential energy landscape!!
        #more important for frustrated structures to equilibrate here, rather than at the end
        #step1 add temperature for a few timesteps
        #step2 athermal equilibrate. This should avoid numerical problems

        relax(next_capsid, RNG)

        configurations.append(next_capsid)
        generate_next_state(next_capsid)
        
        #relax(next_capsid, RNG)
        
        new_wedge_pairs = get_new_wedge_pairs(next_capsid, next_capsid.subunits[-1])
        print new_wedge_pairs, len(next_capsid.subunits)
        for wedge_pair in new_wedge_pairs:
            next_capsid = fuse_copy(next_capsid, wedge_pair)
            relax(next_capsid, RNG)
            configurations.append(next_capsid)
            generate_next_state(next_capsid)
            #relax(next_capsid, RNG)
        #get_open_wedge_pairs(next_capsid)
        #for type1_fusion_pairs in ...
            #next_capsid = fuse_copy()
            #generate_next_state
generate_next_state(capsid)
886/278: len(configurations)
886/279:
pent_conf = [conf for conf in configurations if len(conf.subunits)==5]
pent_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==5 ]
886/280: dp.plot_capsid(pent_conf[0])
886/281: dp.plot_capsid(pent_conf2[0])
886/282:
elastic_modulus=50.0
kT = 1.0
bending_modulus = 300.0
theta = 0.2

#!!! NO EXCLUDERS FOR CONNECTIVITY

params= {'global_params':{'ensemble':2, 'l_fusion':1.0, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':30000000, 'dtsave':2000, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                             'excluders':{'exc1':'com'},
                             'activity':-3.0,
                             'insertion_rate':0.1,
                             'number_of_rotational_configurations':1}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':theta},
                            'bending_modulus':{'type1':bending_modulus},
                            'binding_affinity':{'type1':5.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },                
                            
                 },
         'excluders':{'com':{'overlaps':{'com':False}, 'radius':0.0}
                     }
                  }
pr.init_globals(params)
886/283:
##GENERATE AND EQUILIBRATE WITHOUT EXCLUDERS. ADD EXCLUDERS AFTERWARDS, REMOVE CONFIGURATIONS WHERE THEY OVERLAP
#assumes each connectivity has a minimum and not a glassy landscape
#compute 3 moments of inertia to identify degenerate structures
#includes l_fusion for fusion, otherwise runs in geometry/relaxation problems
#won't generate too high energy confs due to l_fusion
capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])
configurations = []
def generate_next_state(capsid):
    if len(capsid.subunits)>4:
        return
    print capsid, len(capsid.subunits)
    surface_edges = capsid.get_surface_edges()
    for surface_edge in surface_edges:
        next_capsid = insert_copy(capsid, surface_edge, 'sub1')
        #should equilibrate capsid here to drive it towards the global minimum of the potential energy landscape!!
        #more important for frustrated structures to equilibrate here, rather than at the end
        #step1 add temperature for a few timesteps
        #step2 athermal equilibrate. This should avoid numerical problems

        relax(next_capsid, RNG)

        configurations.append(next_capsid)
        generate_next_state(next_capsid)
        
        #relax(next_capsid, RNG)
        
        new_wedge_pairs = get_new_wedge_pairs(next_capsid, next_capsid.subunits[-1])
        print new_wedge_pairs, len(next_capsid.subunits)
        for wedge_pair in new_wedge_pairs:
            next_capsid = fuse_copy(next_capsid, wedge_pair)
            relax(next_capsid, RNG)
            configurations.append(next_capsid)
            generate_next_state(next_capsid)
            #relax(next_capsid, RNG)
        #get_open_wedge_pairs(next_capsid)
        #for type1_fusion_pairs in ...
            #next_capsid = fuse_copy()
            #generate_next_state
generate_next_state(capsid)
886/284: len(configurations)
886/285:
pent_conf = [conf for conf in configurations if len(conf.subunits)==5]
pent_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==5 ]
886/286: pent_conf2
886/287:
elastic_modulus=50.0
kT = 1.0
bending_modulus = 300.0
theta = 0.3

#!!! NO EXCLUDERS FOR CONNECTIVITY

params= {'global_params':{'ensemble':2, 'l_fusion':1.0, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':30000000, 'dtsave':2000, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                             'excluders':{'exc1':'com'},
                             'activity':-3.0,
                             'insertion_rate':0.1,
                             'number_of_rotational_configurations':1}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':theta},
                            'bending_modulus':{'type1':bending_modulus},
                            'binding_affinity':{'type1':5.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },                
                            
                 },
         'excluders':{'com':{'overlaps':{'com':False}, 'radius':0.0}
                     }
                  }
pr.init_globals(params)
886/288: len(configurations)
886/289:
elastic_modulus=50.0
kT = 1.0
bending_modulus = 300.0
theta = 0.3

#!!! NO EXCLUDERS FOR CONNECTIVITY

params= {'global_params':{'ensemble':2, 'l_fusion':1.0, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':30000000, 'dtsave':2000, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                             'excluders':{'exc1':'com'},
                             'activity':-3.0,
                             'insertion_rate':0.1,
                             'number_of_rotational_configurations':1}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':theta},
                            'bending_modulus':{'type1':bending_modulus},
                            'binding_affinity':{'type1':5.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },                
                            
                 },
         'excluders':{'com':{'overlaps':{'com':False}, 'radius':0.0}
                     }
                  }
pr.init_globals(params)
886/290:
##GENERATE AND EQUILIBRATE WITHOUT EXCLUDERS. ADD EXCLUDERS AFTERWARDS, REMOVE CONFIGURATIONS WHERE THEY OVERLAP
#assumes each connectivity has a minimum and not a glassy landscape
#compute 3 moments of inertia to identify degenerate structures
#includes l_fusion for fusion, otherwise runs in geometry/relaxation problems
#won't generate too high energy confs due to l_fusion
capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])
configurations = []
def generate_next_state(capsid):
    if len(capsid.subunits)>4:
        return
    print capsid, len(capsid.subunits)
    surface_edges = capsid.get_surface_edges()
    for surface_edge in surface_edges:
        next_capsid = insert_copy(capsid, surface_edge, 'sub1')
        #should equilibrate capsid here to drive it towards the global minimum of the potential energy landscape!!
        #more important for frustrated structures to equilibrate here, rather than at the end
        #step1 add temperature for a few timesteps
        #step2 athermal equilibrate. This should avoid numerical problems

        relax(next_capsid, RNG)

        configurations.append(next_capsid)
        generate_next_state(next_capsid)
        
        #relax(next_capsid, RNG)
        
        new_wedge_pairs = get_new_wedge_pairs(next_capsid, next_capsid.subunits[-1])
        print new_wedge_pairs, len(next_capsid.subunits)
        for wedge_pair in new_wedge_pairs:
            next_capsid = fuse_copy(next_capsid, wedge_pair)
            relax(next_capsid, RNG)
            configurations.append(next_capsid)
            generate_next_state(next_capsid)
            #relax(next_capsid, RNG)
        #get_open_wedge_pairs(next_capsid)
        #for type1_fusion_pairs in ...
            #next_capsid = fuse_copy()
            #generate_next_state
generate_next_state(capsid)
886/291: len(configurations)
886/292:
pent_conf = [conf for conf in configurations if len(conf.subunits)==5]
pent_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==5 ]
886/293: pent_conf2
886/294:
RNG = np.random.RandomState()

#insert new subunit, returns copy
#overlap should not count in establishing connectivities
def insert_copy(capsid, surface_edge, subunit_type, nrotations=0):

    edge = surface_edge   
    v_from = edge.vertex_from
    v_to = edge.vertex_to

    v1 = Vertex(v_to.x,v_to.y,v_to.z, None, None)
    v2 = Vertex(v_from.x,v_from.y,v_from.z, None, None) 
    #compute position of new vertex
    r0 = 0.5*np.array([v1.x+v2.x, v1.y+v2.y, v1.z+v2.z])
    n = get_triangle_normal(edge.parent_subunit)
    a = np.array([v1.x, v1.y, v1.z])-r0
    nxa = np.array(cross3D(n,a)).astype(float)
    v3_pos = -nxa/norm3D(nxa)*edge.equilibrium_length*np.sqrt(3.0)/2.0
    #print "v3 ", v3_pos

    v3 = Vertex(v3_pos[0],v3_pos[1],v3_pos[2], None, None)

    sub = create_subunit_from_vertices(v1, v2, v3, subunit_type)
     
    if nrotations==0:
        pass
    if nrotations>0:
        #rotate once, then rename
        xtmp, ytmp, ztmp = v3.x, v3.y, v3.z
        v3.x, v3.y, v3.z = v2.x, v2.y, v2.z
        v2.x, v2.y, v2.z = v1.x, v1.y, v1.z
        v1.x, v1.y, v1.z = xtmp, ytmp, ztmp
        vtmp = v2
        v2 = v3
        v3 = v1
        v1 = vtmp
    if nrotations>1:
        #rotate again, then rename
        xtmp, ytmp, ztmp = v3.x, v3.y, v3.z
        v3.x, v3.y, v3.z = v2.x, v2.y, v2.z
        v2.x, v2.y, v2.z = v1.x, v1.y, v1.z
        v1.x, v1.y, v1.z = xtmp, ytmp, ztmp
        vtmp = v2
        v2 = v3
        v3 = v1
        v1 = vtmp

    edge1 = edge
    edge2 = get_edges_between_hubs(v1.parent_hub, v2.parent_hub)[0] 
    #!!!!!!!!
    #!!!!!!!!
    # check for incompatible edges here

    theta = edge1.equilibrium_angle[edge2.edge_type]#pr.equilibrium_angle
    axis = a/norm3D(a)
    rotate_vertex(v3,axis,theta)
    
    v3.x+=r0[0]
    v3.y+=r0[1]
    v3.z+=r0[2]   


    #update_excluders_position([sub])
    #if check_overlap([sub], capsid):
    #    return None
    capsid.add_new_subunit(sub)
    capsid.connect_subunit(v_to.parent_hub, v_from.parent_hub, v1.parent_hub, v2.parent_hub)
    new_capsid = copy.deepcopy(capsid)
    capsid.disconnect_subunit(sub)
    capsid.remove_subunit(sub)
    
    return new_capsid



def fuse_copy(capsid, wedge_pair):
    #get_open_wedge_pairs(capsid)
    
    hub1, hub2 = wedge_pair
    hub3 = get_common_neighbor_hubs(hub1, hub2)[0]
    edge1 = get_edges_between_hubs(hub1, hub3)[0]
    edge2 = get_edges_between_hubs(hub2, hub3)[0]

    affected_subunits = list(set([vertex.parent_subunit for vertex in hub1.vertices]+
                                 [vertex.parent_subunit for vertex in hub2.vertices]))

    energy_before = elastic_energy(affected_subunits)
   
    
    dx = 0.5*(hub1.vertices[0].x - hub2.vertices[0].x)
    dy = 0.5*(hub1.vertices[0].y - hub2.vertices[0].y)
    dz = 0.5*(hub1.vertices[0].z - hub2.vertices[0].z)
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!! This should be a rigid body rotation for merge
    # otherwise theta=0 structures are too compressed and never relax
    # Rotate at least the newly inserted hub
    #!!!!!!!!!!!!!!!!!!!
    
    hub1.move(-dx, -dy, -dz)
    hub2.move(dx, dy, dz)

    #update_excluders_position(affected_subunits)
    #if check_overlap(affected_subunits, capsid):
    #    hub1.move(dx, dy, dz)
    #    hub2.move(-dx, -dy, -dz)
    #    update_excluders_position(affected_subunits)
    #    return False
    #print "merge start", number_of_edges_between_hubs(hub1, hub2)
    hub1_initial_vertices = hub1.vertices #potential source of bug: deep vs shallow copy
    hub1.merge_into(hub2)
    
    energy_after = elastic_energy(affected_subunits)

    dE = energy_after - energy_before
    

    new_capsid = copy.deepcopy(capsid)
        
    hub1 = hub2.split_from(hub1_initial_vertices, hub1)
    hub1.move(dx, dy, dz)
    hub2.move(-dx, -dy, -dz)
    #update_excluders_position(affected_subunits)
        #limit to avoid dumb fusions
    if dE>10.0:
        return None
    return new_capsid

def get_new_wedge_pairs(capsid, new_sub):
    surface_hubs = capsid.get_surface_hubs()
    new_sub_hubs = [v.parent_hub for v in new_sub.vertices]
    
    hub2 = filter(lambda h: len(h.vertices)==1, new_sub_hubs)
    if len(hub2)==0:
        return []
    else:
        hub2 = hub2[0]
    
    wedge_pairs = []
    for hub1 in surface_hubs:
            if hub1==hub2:
                continue
            #should check if hub1 and hub2 are not connected neighbors!!
            if number_of_edges_between_hubs(hub1, hub2)>0:
                continue
            #type 1 fusion is wedge closure; need a single common neighbor 
            common_neighbors = get_common_neighbor_hubs(hub1, hub2)
            if len(common_neighbors)<>1:
                continue

            common_neighbor = common_neighbors[0]
            #both edges have to be surface edges
            if number_of_edges_between_hubs(hub1, common_neighbor)>1:
                continue
            if number_of_edges_between_hubs(hub2, common_neighbor)>1:
                continue                
            
            #if distance(hub1.vertices[0], hub2.vertices[0])<pr.l_fusion:
            #if within_cube(hub1.vertices[0], hub2.vertices[0], pr.l_fusion):
            wedge_pairs.append([hub1, hub2])
    fusion_pairs = wedge_pairs#filter(lambda h: within_cube(h[0].vertices[0], h[1].vertices[0], pr.l_fusion) ,wedge_pairs)
    return fusion_pairs

def relax(capsid, RNG):
    for sub in capsid.subunits:
        sub.neighbor_list = []
    for i in range(5*len(capsid.subunits)):
        mc.attempt_move(capsid, RNG)
886/295:
elastic_modulus=50.0
kT = 1.0
bending_modulus = 300.0
theta = 0.3

#!!! NO EXCLUDERS FOR CONNECTIVITY

params= {'global_params':{'ensemble':2, 'l_fusion':1.0, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':30000000, 'dtsave':2000, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                             'excluders':{'exc1':'com'},
                             'activity':-3.0,
                             'insertion_rate':0.1,
                             'number_of_rotational_configurations':1}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':theta},
                            'bending_modulus':{'type1':bending_modulus},
                            'binding_affinity':{'type1':5.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },                
                            
                 },
         'excluders':{'com':{'overlaps':{'com':False}, 'radius':0.0}
                     }
                  }
pr.init_globals(params)
886/296:
##GENERATE AND EQUILIBRATE WITHOUT EXCLUDERS. ADD EXCLUDERS AFTERWARDS, REMOVE CONFIGURATIONS WHERE THEY OVERLAP
#assumes each connectivity has a minimum and not a glassy landscape
#compute 3 moments of inertia to identify degenerate structures
#includes l_fusion for fusion, otherwise runs in geometry/relaxation problems
#won't generate too high energy confs due to l_fusion
capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])
configurations = []
def generate_next_state(capsid):
    if len(capsid.subunits)>4:
        return
    print capsid, len(capsid.subunits)
    surface_edges = capsid.get_surface_edges()
    for surface_edge in surface_edges:
        next_capsid = insert_copy(capsid, surface_edge, 'sub1')
        #should equilibrate capsid here to drive it towards the global minimum of the potential energy landscape!!
        #more important for frustrated structures to equilibrate here, rather than at the end
        #step1 add temperature for a few timesteps
        #step2 athermal equilibrate. This should avoid numerical problems

        relax(next_capsid, RNG)

        configurations.append(next_capsid)
        generate_next_state(next_capsid)
        
        #relax(next_capsid, RNG)
        
        new_wedge_pairs = get_new_wedge_pairs(next_capsid, next_capsid.subunits[-1])
        print new_wedge_pairs, len(next_capsid.subunits)
        for wedge_pair in new_wedge_pairs:
            next_capsid = fuse_copy(next_capsid, wedge_pair)
            if next_capsid is None:
                continue
            relax(next_capsid, RNG)
            configurations.append(next_capsid)
            generate_next_state(next_capsid)
            #relax(next_capsid, RNG)
        #get_open_wedge_pairs(next_capsid)
        #for type1_fusion_pairs in ...
            #next_capsid = fuse_copy()
            #generate_next_state
generate_next_state(capsid)
886/297: len(configurations)
886/298:
pent_conf = [conf for conf in configurations if len(conf.subunits)==5]
pent_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==5 ]
886/299: pent_conf2
886/300: np.exp(-10)
886/301:
RNG = np.random.RandomState()

#insert new subunit, returns copy
#overlap should not count in establishing connectivities
def insert_copy(capsid, surface_edge, subunit_type, nrotations=0):

    edge = surface_edge   
    v_from = edge.vertex_from
    v_to = edge.vertex_to

    v1 = Vertex(v_to.x,v_to.y,v_to.z, None, None)
    v2 = Vertex(v_from.x,v_from.y,v_from.z, None, None) 
    #compute position of new vertex
    r0 = 0.5*np.array([v1.x+v2.x, v1.y+v2.y, v1.z+v2.z])
    n = get_triangle_normal(edge.parent_subunit)
    a = np.array([v1.x, v1.y, v1.z])-r0
    nxa = np.array(cross3D(n,a)).astype(float)
    v3_pos = -nxa/norm3D(nxa)*edge.equilibrium_length*np.sqrt(3.0)/2.0
    #print "v3 ", v3_pos

    v3 = Vertex(v3_pos[0],v3_pos[1],v3_pos[2], None, None)

    sub = create_subunit_from_vertices(v1, v2, v3, subunit_type)
     
    if nrotations==0:
        pass
    if nrotations>0:
        #rotate once, then rename
        xtmp, ytmp, ztmp = v3.x, v3.y, v3.z
        v3.x, v3.y, v3.z = v2.x, v2.y, v2.z
        v2.x, v2.y, v2.z = v1.x, v1.y, v1.z
        v1.x, v1.y, v1.z = xtmp, ytmp, ztmp
        vtmp = v2
        v2 = v3
        v3 = v1
        v1 = vtmp
    if nrotations>1:
        #rotate again, then rename
        xtmp, ytmp, ztmp = v3.x, v3.y, v3.z
        v3.x, v3.y, v3.z = v2.x, v2.y, v2.z
        v2.x, v2.y, v2.z = v1.x, v1.y, v1.z
        v1.x, v1.y, v1.z = xtmp, ytmp, ztmp
        vtmp = v2
        v2 = v3
        v3 = v1
        v1 = vtmp

    edge1 = edge
    edge2 = get_edges_between_hubs(v1.parent_hub, v2.parent_hub)[0] 
    #!!!!!!!!
    #!!!!!!!!
    # check for incompatible edges here

    theta = edge1.equilibrium_angle[edge2.edge_type]#pr.equilibrium_angle
    axis = a/norm3D(a)
    rotate_vertex(v3,axis,theta)
    
    v3.x+=r0[0]
    v3.y+=r0[1]
    v3.z+=r0[2]   


    #update_excluders_position([sub])
    #if check_overlap([sub], capsid):
    #    return None
    capsid.add_new_subunit(sub)
    capsid.connect_subunit(v_to.parent_hub, v_from.parent_hub, v1.parent_hub, v2.parent_hub)
    new_capsid = copy.deepcopy(capsid)
    capsid.disconnect_subunit(sub)
    capsid.remove_subunit(sub)
    
    return new_capsid



def fuse_copy(capsid, wedge_pair):
    #get_open_wedge_pairs(capsid)
    
    hub1, hub2 = wedge_pair
    hub3 = get_common_neighbor_hubs(hub1, hub2)[0]
    edge1 = get_edges_between_hubs(hub1, hub3)[0]
    edge2 = get_edges_between_hubs(hub2, hub3)[0]

    affected_subunits = list(set([vertex.parent_subunit for vertex in hub1.vertices]+
                                 [vertex.parent_subunit for vertex in hub2.vertices]))

    energy_before = elastic_energy(affected_subunits)
   
    
    dx = 0.5*(hub1.vertices[0].x - hub2.vertices[0].x)
    dy = 0.5*(hub1.vertices[0].y - hub2.vertices[0].y)
    dz = 0.5*(hub1.vertices[0].z - hub2.vertices[0].z)
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!! This should be a rigid body rotation for merge
    # otherwise theta=0 structures are too compressed and never relax
    # Rotate at least the newly inserted hub
    #!!!!!!!!!!!!!!!!!!!
    
    hub1.move(-dx, -dy, -dz)
    hub2.move(dx, dy, dz)

    #update_excluders_position(affected_subunits)
    #if check_overlap(affected_subunits, capsid):
    #    hub1.move(dx, dy, dz)
    #    hub2.move(-dx, -dy, -dz)
    #    update_excluders_position(affected_subunits)
    #    return False
    #print "merge start", number_of_edges_between_hubs(hub1, hub2)
    hub1_initial_vertices = hub1.vertices #potential source of bug: deep vs shallow copy
    hub1.merge_into(hub2)
    
    energy_after = elastic_energy(affected_subunits)

    dE = energy_after - energy_before
    

    new_capsid = copy.deepcopy(capsid)
        
    hub1 = hub2.split_from(hub1_initial_vertices, hub1)
    hub1.move(dx, dy, dz)
    hub2.move(-dx, -dy, -dz)
    #update_excluders_position(affected_subunits)
        #limit to avoid dumb fusions
    if dE>20.0:
        return None
    return new_capsid

def get_new_wedge_pairs(capsid, new_sub):
    surface_hubs = capsid.get_surface_hubs()
    new_sub_hubs = [v.parent_hub for v in new_sub.vertices]
    
    hub2 = filter(lambda h: len(h.vertices)==1, new_sub_hubs)
    if len(hub2)==0:
        return []
    else:
        hub2 = hub2[0]
    
    wedge_pairs = []
    for hub1 in surface_hubs:
            if hub1==hub2:
                continue
            #should check if hub1 and hub2 are not connected neighbors!!
            if number_of_edges_between_hubs(hub1, hub2)>0:
                continue
            #type 1 fusion is wedge closure; need a single common neighbor 
            common_neighbors = get_common_neighbor_hubs(hub1, hub2)
            if len(common_neighbors)<>1:
                continue

            common_neighbor = common_neighbors[0]
            #both edges have to be surface edges
            if number_of_edges_between_hubs(hub1, common_neighbor)>1:
                continue
            if number_of_edges_between_hubs(hub2, common_neighbor)>1:
                continue                
            
            #if distance(hub1.vertices[0], hub2.vertices[0])<pr.l_fusion:
            #if within_cube(hub1.vertices[0], hub2.vertices[0], pr.l_fusion):
            wedge_pairs.append([hub1, hub2])
    fusion_pairs = wedge_pairs#filter(lambda h: within_cube(h[0].vertices[0], h[1].vertices[0], pr.l_fusion) ,wedge_pairs)
    return fusion_pairs

def relax(capsid, RNG):
    for sub in capsid.subunits:
        sub.neighbor_list = []
    for i in range(5*len(capsid.subunits)):
        mc.attempt_move(capsid, RNG)
886/302:
elastic_modulus=50.0
kT = 1.0
bending_modulus = 300.0
theta = 0.3

#!!! NO EXCLUDERS FOR CONNECTIVITY

params= {'global_params':{'ensemble':2, 'l_fusion':1.0, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':30000000, 'dtsave':2000, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                             'excluders':{'exc1':'com'},
                             'activity':-3.0,
                             'insertion_rate':0.1,
                             'number_of_rotational_configurations':1}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':theta},
                            'bending_modulus':{'type1':bending_modulus},
                            'binding_affinity':{'type1':5.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },                
                            
                 },
         'excluders':{'com':{'overlaps':{'com':False}, 'radius':0.0}
                     }
                  }
pr.init_globals(params)
886/303:
##GENERATE AND EQUILIBRATE WITHOUT EXCLUDERS. ADD EXCLUDERS AFTERWARDS, REMOVE CONFIGURATIONS WHERE THEY OVERLAP
#assumes each connectivity has a minimum and not a glassy landscape
#compute 3 moments of inertia to identify degenerate structures
#includes l_fusion for fusion, otherwise runs in geometry/relaxation problems
#won't generate too high energy confs due to l_fusion
capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])
configurations = []
def generate_next_state(capsid):
    if len(capsid.subunits)>4:
        return
    print capsid, len(capsid.subunits)
    surface_edges = capsid.get_surface_edges()
    for surface_edge in surface_edges:
        next_capsid = insert_copy(capsid, surface_edge, 'sub1')
        #should equilibrate capsid here to drive it towards the global minimum of the potential energy landscape!!
        #more important for frustrated structures to equilibrate here, rather than at the end
        #step1 add temperature for a few timesteps
        #step2 athermal equilibrate. This should avoid numerical problems

        relax(next_capsid, RNG)

        configurations.append(next_capsid)
        generate_next_state(next_capsid)
        
        #relax(next_capsid, RNG)
        
        new_wedge_pairs = get_new_wedge_pairs(next_capsid, next_capsid.subunits[-1])
        print new_wedge_pairs, len(next_capsid.subunits)
        for wedge_pair in new_wedge_pairs:
            next_capsid = fuse_copy(next_capsid, wedge_pair)
            if next_capsid is None:
                continue
            relax(next_capsid, RNG)
            configurations.append(next_capsid)
            generate_next_state(next_capsid)
            #relax(next_capsid, RNG)
        #get_open_wedge_pairs(next_capsid)
        #for type1_fusion_pairs in ...
            #next_capsid = fuse_copy()
            #generate_next_state
generate_next_state(capsid)
886/304: len(configurations)
886/305:
pent_conf = [conf for conf in configurations if len(conf.subunits)==5]
pent_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==5 ]
886/306: pent_conf2
886/307:
elastic_modulus=50.0
kT = 1.0
bending_modulus = 50.0
theta = 0.3

#!!! NO EXCLUDERS FOR CONNECTIVITY

params= {'global_params':{'ensemble':2, 'l_fusion':1.0, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':30000000, 'dtsave':2000, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                             'excluders':{'exc1':'com'},
                             'activity':-3.0,
                             'insertion_rate':0.1,
                             'number_of_rotational_configurations':1}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':theta},
                            'bending_modulus':{'type1':bending_modulus},
                            'binding_affinity':{'type1':5.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },                
                            
                 },
         'excluders':{'com':{'overlaps':{'com':False}, 'radius':0.0}
                     }
                  }
pr.init_globals(params)
886/308:
##GENERATE AND EQUILIBRATE WITHOUT EXCLUDERS. ADD EXCLUDERS AFTERWARDS, REMOVE CONFIGURATIONS WHERE THEY OVERLAP
#assumes each connectivity has a minimum and not a glassy landscape
#compute 3 moments of inertia to identify degenerate structures
#includes l_fusion for fusion, otherwise runs in geometry/relaxation problems
#won't generate too high energy confs due to l_fusion
capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])
configurations = []
def generate_next_state(capsid):
    if len(capsid.subunits)>4:
        return
    print capsid, len(capsid.subunits)
    surface_edges = capsid.get_surface_edges()
    for surface_edge in surface_edges:
        next_capsid = insert_copy(capsid, surface_edge, 'sub1')
        #should equilibrate capsid here to drive it towards the global minimum of the potential energy landscape!!
        #more important for frustrated structures to equilibrate here, rather than at the end
        #step1 add temperature for a few timesteps
        #step2 athermal equilibrate. This should avoid numerical problems

        relax(next_capsid, RNG)

        configurations.append(next_capsid)
        generate_next_state(next_capsid)
        
        #relax(next_capsid, RNG)
        
        new_wedge_pairs = get_new_wedge_pairs(next_capsid, next_capsid.subunits[-1])
        print new_wedge_pairs, len(next_capsid.subunits)
        for wedge_pair in new_wedge_pairs:
            next_capsid = fuse_copy(next_capsid, wedge_pair)
            if next_capsid is None:
                continue
            relax(next_capsid, RNG)
            configurations.append(next_capsid)
            generate_next_state(next_capsid)
            #relax(next_capsid, RNG)
        #get_open_wedge_pairs(next_capsid)
        #for type1_fusion_pairs in ...
            #next_capsid = fuse_copy()
            #generate_next_state
generate_next_state(capsid)
886/309: len(configurations)
886/310:
pent_conf = [conf for conf in configurations if len(conf.subunits)==5]
pent_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==5 ]
886/311: pent_conf2
886/312:
pent_conf = [conf for conf in configurations if len(conf.subunits)==5]
pent_conf2 = [conf for conf in pent_conf if  len(conf.get_surface_edges())==5 ]
886/313: pent_conf2
886/314: dp.plot_capsid(pent_conf2[0])
886/315: dp.plot_capsid(pent_conf2[10])
886/316: dp.plot_capsid(pent_conf2[30])
886/317:
RNG = np.random.RandomState()

#insert new subunit, returns copy
#overlap should not count in establishing connectivities
def insert_copy(capsid, surface_edge, subunit_type, nrotations=0):

    edge = surface_edge   
    v_from = edge.vertex_from
    v_to = edge.vertex_to

    v1 = Vertex(v_to.x,v_to.y,v_to.z, None, None)
    v2 = Vertex(v_from.x,v_from.y,v_from.z, None, None) 
    #compute position of new vertex
    r0 = 0.5*np.array([v1.x+v2.x, v1.y+v2.y, v1.z+v2.z])
    n = get_triangle_normal(edge.parent_subunit)
    a = np.array([v1.x, v1.y, v1.z])-r0
    nxa = np.array(cross3D(n,a)).astype(float)
    v3_pos = -nxa/norm3D(nxa)*edge.equilibrium_length*np.sqrt(3.0)/2.0
    #print "v3 ", v3_pos

    v3 = Vertex(v3_pos[0],v3_pos[1],v3_pos[2], None, None)

    sub = create_subunit_from_vertices(v1, v2, v3, subunit_type)
     
    if nrotations==0:
        pass
    if nrotations>0:
        #rotate once, then rename
        xtmp, ytmp, ztmp = v3.x, v3.y, v3.z
        v3.x, v3.y, v3.z = v2.x, v2.y, v2.z
        v2.x, v2.y, v2.z = v1.x, v1.y, v1.z
        v1.x, v1.y, v1.z = xtmp, ytmp, ztmp
        vtmp = v2
        v2 = v3
        v3 = v1
        v1 = vtmp
    if nrotations>1:
        #rotate again, then rename
        xtmp, ytmp, ztmp = v3.x, v3.y, v3.z
        v3.x, v3.y, v3.z = v2.x, v2.y, v2.z
        v2.x, v2.y, v2.z = v1.x, v1.y, v1.z
        v1.x, v1.y, v1.z = xtmp, ytmp, ztmp
        vtmp = v2
        v2 = v3
        v3 = v1
        v1 = vtmp

    edge1 = edge
    edge2 = get_edges_between_hubs(v1.parent_hub, v2.parent_hub)[0] 
    #!!!!!!!!
    #!!!!!!!!
    # check for incompatible edges here

    theta = edge1.equilibrium_angle[edge2.edge_type]#pr.equilibrium_angle
    axis = a/norm3D(a)
    rotate_vertex(v3,axis,theta)
    
    v3.x+=r0[0]
    v3.y+=r0[1]
    v3.z+=r0[2]   


    #update_excluders_position([sub])
    #if check_overlap([sub], capsid):
    #    return None
    capsid.add_new_subunit(sub)
    capsid.connect_subunit(v_to.parent_hub, v_from.parent_hub, v1.parent_hub, v2.parent_hub)
    new_capsid = copy.deepcopy(capsid)
    capsid.disconnect_subunit(sub)
    capsid.remove_subunit(sub)
    
    return new_capsid



def fuse_copy(capsid, wedge_pair):
    #get_open_wedge_pairs(capsid)
    
    hub1, hub2 = wedge_pair
    hub3 = get_common_neighbor_hubs(hub1, hub2)[0]
    edge1 = get_edges_between_hubs(hub1, hub3)[0]
    edge2 = get_edges_between_hubs(hub2, hub3)[0]

    affected_subunits = list(set([vertex.parent_subunit for vertex in hub1.vertices]+
                                 [vertex.parent_subunit for vertex in hub2.vertices]))

    energy_before = elastic_energy(affected_subunits)
   
    
    dx = 0.5*(hub1.vertices[0].x - hub2.vertices[0].x)
    dy = 0.5*(hub1.vertices[0].y - hub2.vertices[0].y)
    dz = 0.5*(hub1.vertices[0].z - hub2.vertices[0].z)
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!! This should be a rigid body rotation for merge
    # otherwise theta=0 structures are too compressed and never relax
    # Rotate at least the newly inserted hub
    #!!!!!!!!!!!!!!!!!!!
    
    hub1.move(-dx, -dy, -dz)
    hub2.move(dx, dy, dz)

    #update_excluders_position(affected_subunits)
    #if check_overlap(affected_subunits, capsid):
    #    hub1.move(dx, dy, dz)
    #    hub2.move(-dx, -dy, -dz)
    #    update_excluders_position(affected_subunits)
    #    return False
    #print "merge start", number_of_edges_between_hubs(hub1, hub2)
    hub1_initial_vertices = hub1.vertices #potential source of bug: deep vs shallow copy
    hub1.merge_into(hub2)
    
    energy_after = elastic_energy(affected_subunits)

    dE = energy_after - energy_before
    

    new_capsid = copy.deepcopy(capsid)
        
    hub1 = hub2.split_from(hub1_initial_vertices, hub1)
    hub1.move(dx, dy, dz)
    hub2.move(-dx, -dy, -dz)
    #update_excluders_position(affected_subunits)
        #limit to avoid dumb fusions
    #if dE>20.0:
    #    return None
    return new_capsid

def get_new_wedge_pairs(capsid, new_sub):
    surface_hubs = capsid.get_surface_hubs()
    new_sub_hubs = [v.parent_hub for v in new_sub.vertices]
    
    hub2 = filter(lambda h: len(h.vertices)==1, new_sub_hubs)
    if len(hub2)==0:
        return []
    else:
        hub2 = hub2[0]
    
    wedge_pairs = []
    for hub1 in surface_hubs:
            if hub1==hub2:
                continue
            #should check if hub1 and hub2 are not connected neighbors!!
            if number_of_edges_between_hubs(hub1, hub2)>0:
                continue
            #type 1 fusion is wedge closure; need a single common neighbor 
            common_neighbors = get_common_neighbor_hubs(hub1, hub2)
            if len(common_neighbors)<>1:
                continue

            common_neighbor = common_neighbors[0]
            #both edges have to be surface edges
            if number_of_edges_between_hubs(hub1, common_neighbor)>1:
                continue
            if number_of_edges_between_hubs(hub2, common_neighbor)>1:
                continue                
            
            #if distance(hub1.vertices[0], hub2.vertices[0])<pr.l_fusion:
            #if within_cube(hub1.vertices[0], hub2.vertices[0], pr.l_fusion):
            wedge_pairs.append([hub1, hub2])
    fusion_pairs = filter(lambda h: within_cube(h[0].vertices[0], h[1].vertices[0], pr.l_fusion) ,wedge_pairs)
    return fusion_pairs

def relax(capsid, RNG):
    for sub in capsid.subunits:
        sub.neighbor_list = []
    for i in range(5*len(capsid.subunits)):
        mc.attempt_move(capsid, RNG)
886/318:
elastic_modulus=50.0
kT = 1.0
bending_modulus = 50.0
theta = 0.3

#!!! NO EXCLUDERS FOR CONNECTIVITY

params= {'global_params':{'ensemble':2, 'l_fusion':1.0, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':30000000, 'dtsave':2000, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                             'excluders':{'exc1':'com'},
                             'activity':-3.0,
                             'insertion_rate':0.1,
                             'number_of_rotational_configurations':1}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':theta},
                            'bending_modulus':{'type1':bending_modulus},
                            'binding_affinity':{'type1':5.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },                
                            
                 },
         'excluders':{'com':{'overlaps':{'com':False}, 'radius':0.0}
                     }
                  }
pr.init_globals(params)
886/319:
##GENERATE AND EQUILIBRATE WITHOUT EXCLUDERS. ADD EXCLUDERS AFTERWARDS, REMOVE CONFIGURATIONS WHERE THEY OVERLAP
#assumes each connectivity has a minimum and not a glassy landscape
#compute 3 moments of inertia to identify degenerate structures
#includes l_fusion for fusion, otherwise runs in geometry/relaxation problems
#won't generate too high energy confs due to l_fusion
capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])
configurations = []
def generate_next_state(capsid):
    if len(capsid.subunits)>4:
        return
    print capsid, len(capsid.subunits)
    surface_edges = capsid.get_surface_edges()
    for surface_edge in surface_edges:
        next_capsid = insert_copy(capsid, surface_edge, 'sub1')
        #should equilibrate capsid here to drive it towards the global minimum of the potential energy landscape!!
        #more important for frustrated structures to equilibrate here, rather than at the end
        #step1 add temperature for a few timesteps
        #step2 athermal equilibrate. This should avoid numerical problems

        relax(next_capsid, RNG)

        configurations.append(next_capsid)
        generate_next_state(next_capsid)
        
        #relax(next_capsid, RNG)
        
        new_wedge_pairs = get_new_wedge_pairs(next_capsid, next_capsid.subunits[-1])
        print new_wedge_pairs, len(next_capsid.subunits)
        for wedge_pair in new_wedge_pairs:
            next_capsid = fuse_copy(next_capsid, wedge_pair)
            if next_capsid is None:
                continue
            relax(next_capsid, RNG)
            configurations.append(next_capsid)
            generate_next_state(next_capsid)
            #relax(next_capsid, RNG)
        #get_open_wedge_pairs(next_capsid)
        #for type1_fusion_pairs in ...
            #next_capsid = fuse_copy()
            #generate_next_state
generate_next_state(capsid)
886/320: len(configurations)
886/321:
pent_conf = [conf for conf in configurations if len(conf.subunits)==5]
pent_conf2 = [conf for conf in pent_conf if  len(conf.get_surface_edges())==5 ]
886/322: pent_conf2
886/323: dp.plot_capsid(pent_conf2[30])
886/324: dp.plot_capsid(pent_conf2[10])
886/325: hex_conf = [conf for conf in configurations if len(conf.subunits)==6]
886/326: hex_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==8 ]
886/327: hex_conf2
886/328:
RNG = np.random.RandomState()

#insert new subunit, returns copy
#overlap should not count in establishing connectivities
def insert_copy(capsid, surface_edge, subunit_type, nrotations=0):

    edge = surface_edge   
    v_from = edge.vertex_from
    v_to = edge.vertex_to

    v1 = Vertex(v_to.x,v_to.y,v_to.z, None, None)
    v2 = Vertex(v_from.x,v_from.y,v_from.z, None, None) 
    #compute position of new vertex
    r0 = 0.5*np.array([v1.x+v2.x, v1.y+v2.y, v1.z+v2.z])
    n = get_triangle_normal(edge.parent_subunit)
    a = np.array([v1.x, v1.y, v1.z])-r0
    nxa = np.array(cross3D(n,a)).astype(float)
    v3_pos = -nxa/norm3D(nxa)*edge.equilibrium_length*np.sqrt(3.0)/2.0
    #print "v3 ", v3_pos

    v3 = Vertex(v3_pos[0],v3_pos[1],v3_pos[2], None, None)

    sub = create_subunit_from_vertices(v1, v2, v3, subunit_type)
     
    if nrotations==0:
        pass
    if nrotations>0:
        #rotate once, then rename
        xtmp, ytmp, ztmp = v3.x, v3.y, v3.z
        v3.x, v3.y, v3.z = v2.x, v2.y, v2.z
        v2.x, v2.y, v2.z = v1.x, v1.y, v1.z
        v1.x, v1.y, v1.z = xtmp, ytmp, ztmp
        vtmp = v2
        v2 = v3
        v3 = v1
        v1 = vtmp
    if nrotations>1:
        #rotate again, then rename
        xtmp, ytmp, ztmp = v3.x, v3.y, v3.z
        v3.x, v3.y, v3.z = v2.x, v2.y, v2.z
        v2.x, v2.y, v2.z = v1.x, v1.y, v1.z
        v1.x, v1.y, v1.z = xtmp, ytmp, ztmp
        vtmp = v2
        v2 = v3
        v3 = v1
        v1 = vtmp

    edge1 = edge
    edge2 = get_edges_between_hubs(v1.parent_hub, v2.parent_hub)[0] 
    #!!!!!!!!
    #!!!!!!!!
    # check for incompatible edges here

    theta = edge1.equilibrium_angle[edge2.edge_type]#pr.equilibrium_angle
    axis = a/norm3D(a)
    rotate_vertex(v3,axis,theta)
    
    v3.x+=r0[0]
    v3.y+=r0[1]
    v3.z+=r0[2]   


    #update_excluders_position([sub])
    #if check_overlap([sub], capsid):
    #    return None
    capsid.add_new_subunit(sub)
    capsid.connect_subunit(v_to.parent_hub, v_from.parent_hub, v1.parent_hub, v2.parent_hub)
    new_capsid = copy.deepcopy(capsid)
    capsid.disconnect_subunit(sub)
    capsid.remove_subunit(sub)
    
    return new_capsid



def fuse_copy(capsid, wedge_pair):
    #get_open_wedge_pairs(capsid)
    
    hub1, hub2 = wedge_pair
    hub3 = get_common_neighbor_hubs(hub1, hub2)[0]
    edge1 = get_edges_between_hubs(hub1, hub3)[0]
    edge2 = get_edges_between_hubs(hub2, hub3)[0]

    affected_subunits = list(set([vertex.parent_subunit for vertex in hub1.vertices]+
                                 [vertex.parent_subunit for vertex in hub2.vertices]))

    #@energy_before = elastic_energy(affected_subunits)
   
    
    dx = 0.5*(hub1.vertices[0].x - hub2.vertices[0].x)
    dy = 0.5*(hub1.vertices[0].y - hub2.vertices[0].y)
    dz = 0.5*(hub1.vertices[0].z - hub2.vertices[0].z)
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!! This should be a rigid body rotation for merge
    # otherwise theta=0 structures are too compressed and never relax
    # Rotate at least the newly inserted hub
    #!!!!!!!!!!!!!!!!!!!
    
    hub1.move(-dx, -dy, -dz)
    hub2.move(dx, dy, dz)

    #update_excluders_position(affected_subunits)
    #if check_overlap(affected_subunits, capsid):
    #    hub1.move(dx, dy, dz)
    #    hub2.move(-dx, -dy, -dz)
    #    update_excluders_position(affected_subunits)
    #    return False
    #print "merge start", number_of_edges_between_hubs(hub1, hub2)
    hub1_initial_vertices = hub1.vertices #potential source of bug: deep vs shallow copy
    hub1.merge_into(hub2)
    
    #@energy_after = elastic_energy(affected_subunits)

    #@dE = energy_after - energy_before
    

    new_capsid = copy.deepcopy(capsid)
        
    hub1 = hub2.split_from(hub1_initial_vertices, hub1)
    hub1.move(dx, dy, dz)
    hub2.move(-dx, -dy, -dz)
    #update_excluders_position(affected_subunits)
        #limit to avoid dumb fusions
    #if dE>20.0:
    #    return None
    return new_capsid

def get_new_wedge_pairs(capsid, new_sub):
    surface_hubs = capsid.get_surface_hubs()
    new_sub_hubs = [v.parent_hub for v in new_sub.vertices]
    
    hub2 = filter(lambda h: len(h.vertices)==1, new_sub_hubs)
    if len(hub2)==0:
        return []
    else:
        hub2 = hub2[0]
    
    wedge_pairs = []
    for hub1 in surface_hubs:
            if hub1==hub2:
                continue
            #should check if hub1 and hub2 are not connected neighbors!!
            if number_of_edges_between_hubs(hub1, hub2)>0:
                continue
            #type 1 fusion is wedge closure; need a single common neighbor 
            common_neighbors = get_common_neighbor_hubs(hub1, hub2)
            if len(common_neighbors)<>1:
                continue

            common_neighbor = common_neighbors[0]
            #both edges have to be surface edges
            if number_of_edges_between_hubs(hub1, common_neighbor)>1:
                continue
            if number_of_edges_between_hubs(hub2, common_neighbor)>1:
                continue                
            
            #if distance(hub1.vertices[0], hub2.vertices[0])<pr.l_fusion:
            #if within_cube(hub1.vertices[0], hub2.vertices[0], pr.l_fusion):
            wedge_pairs.append([hub1, hub2])
    fusion_pairs = filter(lambda h: within_cube(h[0].vertices[0], h[1].vertices[0], pr.l_fusion) ,wedge_pairs)
    return fusion_pairs

def relax(capsid, RNG):
    for sub in capsid.subunits:
        sub.neighbor_list = []
    for i in range(5*len(capsid.subunits)):
        mc.attempt_move(capsid, RNG)
886/329:
elastic_modulus=50.0
kT = 1.0
bending_modulus = 50.0
theta = 0.3

#!!! NO EXCLUDERS FOR CONNECTIVITY

params= {'global_params':{'ensemble':2, 'l_fusion':1.0, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':30000000, 'dtsave':2000, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                             'excluders':{'exc1':'com'},
                             'activity':-3.0,
                             'insertion_rate':0.1,
                             'number_of_rotational_configurations':1}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':theta},
                            'bending_modulus':{'type1':bending_modulus},
                            'binding_affinity':{'type1':5.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },                
                            
                 },
         'excluders':{'com':{'overlaps':{'com':False}, 'radius':0.0}
                     }
                  }
pr.init_globals(params)
886/330:
##GENERATE AND EQUILIBRATE WITHOUT EXCLUDERS. ADD EXCLUDERS AFTERWARDS, REMOVE CONFIGURATIONS WHERE THEY OVERLAP
#assumes each connectivity has a minimum and not a glassy landscape
#compute 3 moments of inertia to identify degenerate structures
#includes l_fusion for fusion, otherwise runs in geometry/relaxation problems
#won't generate too high energy confs due to l_fusion
capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])
configurations = []
def generate_next_state(capsid):
    if len(capsid.subunits)>5:
        return
    print capsid, len(capsid.subunits)
    surface_edges = capsid.get_surface_edges()
    for surface_edge in surface_edges:
        next_capsid = insert_copy(capsid, surface_edge, 'sub1')
        #should equilibrate capsid here to drive it towards the global minimum of the potential energy landscape!!
        #more important for frustrated structures to equilibrate here, rather than at the end
        #step1 add temperature for a few timesteps
        #step2 athermal equilibrate. This should avoid numerical problems

        relax(next_capsid, RNG)

        configurations.append(next_capsid)
        generate_next_state(next_capsid)
        
        #relax(next_capsid, RNG)
        
        new_wedge_pairs = get_new_wedge_pairs(next_capsid, next_capsid.subunits[-1])
        print new_wedge_pairs, len(next_capsid.subunits)
        for wedge_pair in new_wedge_pairs:
            next_capsid = fuse_copy(next_capsid, wedge_pair)
            if next_capsid is None:
                continue
            relax(next_capsid, RNG)
            configurations.append(next_capsid)
            generate_next_state(next_capsid)
            #relax(next_capsid, RNG)
        #get_open_wedge_pairs(next_capsid)
        #for type1_fusion_pairs in ...
            #next_capsid = fuse_copy()
            #generate_next_state
generate_next_state(capsid)
886/331: len(configurations)
886/332: hex_conf = [conf for conf in configurations if len(conf.subunits)==6]
886/333: hex_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==8 ]
886/334: hex_conf2
886/335: dp.plot_capsid(pent_conf2[10])
886/336: dp.plot_capsid(hex_conf2[10])
886/337: hex_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==6 ]
886/338: hex_conf2
886/339: dp.plot_capsid(hex_conf2[10])
886/340: hex_conf = [conf for conf in configurations if len(conf.subunits)==6]
886/341: hex_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==6 ]
886/342: dp.plot_capsid(hex_conf2[10])
886/343: dp.plot_capsid(hex_conf2[15])
886/344: dp.plot_capsid(hex_conf2[25])
886/345: dp.plot_capsid(hex_conf2[55])
886/346: dp.plot_capsid(hex_conf2[45])
886/347: dp.plot_capsid(hex_conf2[0])
886/348: dp.plot_capsid(hex_conf2[-1])
886/349: dp.plot_capsid(hex_conf2[100])
886/350:
elastic_modulus=50.0
kT = 1.0
bending_modulus = 50.0
theta = 0.7

#!!! NO EXCLUDERS FOR CONNECTIVITY

params= {'global_params':{'ensemble':2, 'l_fusion':1.0, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':30000000, 'dtsave':2000, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                             'excluders':{'exc1':'com'},
                             'activity':-3.0,
                             'insertion_rate':0.1,
                             'number_of_rotational_configurations':1}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':theta},
                            'bending_modulus':{'type1':bending_modulus},
                            'binding_affinity':{'type1':5.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },                
                            
                 },
         'excluders':{'com':{'overlaps':{'com':False}, 'radius':0.0}
                     }
                  }
pr.init_globals(params)
886/351:
##GENERATE AND EQUILIBRATE WITHOUT EXCLUDERS. ADD EXCLUDERS AFTERWARDS, REMOVE CONFIGURATIONS WHERE THEY OVERLAP
#assumes each connectivity has a minimum and not a glassy landscape
#compute 3 moments of inertia to identify degenerate structures
#includes l_fusion for fusion, otherwise runs in geometry/relaxation problems
#won't generate too high energy confs due to l_fusion
capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])
configurations = []
def generate_next_state(capsid):
    if len(capsid.subunits)>5:
        return
    print capsid, len(capsid.subunits)
    surface_edges = capsid.get_surface_edges()
    for surface_edge in surface_edges:
        next_capsid = insert_copy(capsid, surface_edge, 'sub1')
        #should equilibrate capsid here to drive it towards the global minimum of the potential energy landscape!!
        #more important for frustrated structures to equilibrate here, rather than at the end
        #step1 add temperature for a few timesteps
        #step2 athermal equilibrate. This should avoid numerical problems

        relax(next_capsid, RNG)

        configurations.append(next_capsid)
        generate_next_state(next_capsid)
        
        #relax(next_capsid, RNG)
        
        new_wedge_pairs = get_new_wedge_pairs(next_capsid, next_capsid.subunits[-1])
        print new_wedge_pairs, len(next_capsid.subunits)
        for wedge_pair in new_wedge_pairs:
            next_capsid = fuse_copy(next_capsid, wedge_pair)
            if next_capsid is None:
                continue
            relax(next_capsid, RNG)
            configurations.append(next_capsid)
            generate_next_state(next_capsid)
            #relax(next_capsid, RNG)
        #get_open_wedge_pairs(next_capsid)
        #for type1_fusion_pairs in ...
            #next_capsid = fuse_copy()
            #generate_next_state
generate_next_state(capsid)
886/352: len(configurations)
886/353: hex_conf = [conf for conf in configurations if len(conf.subunits)==6]
886/354: hex_conf2 = [conf for conf in hex_conf if  len(conf.get_surface_edges())==6 ]
886/355:
trimers = [conf for conf in configurations if len(conf.subunits)==3]
pyr = [conf for conf in trimers if  len(conf.get_surface_edges())==3 ]
886/356: len(pyr)
886/357:
pent_conf = [conf for conf in configurations if len(conf.subunits)==5]
pent_conf2 = [conf for conf in pent_conf if  len(conf.get_surface_edges())==5 ]
886/358: dp.plot_capsid(hex_conf2[100])
886/359: dp.plot_capsid(pent_conf2[100])
886/360:
Graphs = []
for conf in configurations:
    n_hubs = len(conf.hubs)
    A = np.zeros(shape=(n_hubs, n_hubs))
    for ix,hub1 in enumerate(conf.hubs):
        for jx,hub2 in enumerate(conf.hubs):
            A[ix, jx] = len(get_edges_between_hubs(hub1, hub2))
    Graphs.append(A)
886/361: Graphs[-1]
886/362:
import monte_carlo as mc
import parameters as pr
from helpers import *
import copy
import pylab as pl
import Utils as ut
import data_process as dp
from energy import *
import networkx as nx
886/363:
Graphs = []
for conf in configurations:
    n_hubs = len(conf.hubs)
    A = np.zeros(shape=(n_hubs, n_hubs))
    for ix,hub1 in enumerate(conf.hubs):
        for jx,hub2 in enumerate(conf.hubs):
            A[ix, jx] = len(get_edges_between_hubs(hub1, hub2))
    Graphs.append(nx.from_numpy_matrix(A))
886/364:
for g1 in Graphs:
    for g2 in Graphs:
        nx.is_isomorphic(g1, g2)
