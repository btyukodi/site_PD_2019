import monte_carlo as mc
import parameters as pr
from helpers import *
import copy
import pylab as pl
import Utils as ut
import data_process as dp
from energy import *
import networkx as nx
import pandas as pd
import energy as en


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
    if (edge1.bending_modulus[edge2.edge_type]==float('inf')):
        return None

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
    #@if dE>20.0:
    #@    return None
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

def thermalize(capsid, RNG):
    for sub in capsid.subunits:
        sub.neighbor_list = []
    for i in range(5*len(capsid.subunits)):
        mc.attempt_move(capsid, RNG)


        
def moment_of_inertia(capsid):
    posx = np.array([h.vertices[0].x for h in capsid.hubs])
    posy = np.array([h.vertices[0].y for h in capsid.hubs])
    posz = np.array([h.vertices[0].z for h in capsid.hubs])
    comx, comy, comz = np.mean(posx), np.mean(posy), np.mean(posz)
    posx-=comx
    posy-=comy
    posz-=comz
    Ix, Iy, Iz = np.sum(posx**2), np.sum(posy**2), np.sum(posz**2)
    return Ix+Iy+Iz



##GENERATE AND EQUILIBRATE WITHOUT EXCLUDERS. ADD EXCLUDERS AFTERWARDS, REMOVE CONFIGURATIONS WHERE THEY OVERLAP
#assumes each connectivity has a minimum and not a glassy landscape
#compute 3 moments of inertia to identify degenerate structures
#includes l_fusion for fusion, otherwise runs in geometry/relaxation problems
#won't generate too high energy confs due to l_fusion
def generate_most_configurations(params, n_sub_max):
	RNG = np.random.RandomState()
	pr.init_globals(params)
	capsid = mc.create_single_subunit_capsid(subunit_type=params['subunits'].keys()[0])
	configurations = []
	def generate_next_state(capsid):
	    if len(capsid.subunits)>n_sub_max:
	        return
	    #print capsid, len(capsid.subunits)
	    surface_edges = capsid.get_surface_edges()
	    for surface_edge in surface_edges:
	        next_capsid = insert_copy(capsid, surface_edge, 'sub1')
	        if next_capsid is None:
	                continue        
	        #should equilibrate capsid here to drive it towards the global minimum of the potential energy landscape!!
	        #more important for frustrated structures to equilibrate here, rather than at the end
	        #step1 add temperature for a few timesteps
	        #step2 athermal equilibrate. This should avoid numerical problems

	        thermalize(next_capsid, RNG)

	        configurations.append(next_capsid)
	        generate_next_state(next_capsid)
	        
	        #thermalize(next_capsid, RNG)
	        
	        new_wedge_pairs = get_new_wedge_pairs(next_capsid, next_capsid.subunits[-1])
	        #print new_wedge_pairs, len(next_capsid.subunits)
	        for wedge_pair in new_wedge_pairs:
	            next_capsid = fuse_copy(next_capsid, wedge_pair)
	            if next_capsid is None:
	                continue
	            thermalize(next_capsid, RNG)
	            configurations.append(next_capsid)
	            generate_next_state(next_capsid)
	            #thermalize(next_capsid, RNG)
	        #get_open_wedge_pairs(next_capsid)
	        #for type1_fusion_pairs in ...
	            #next_capsid = fuse_copy()
	            #generate_next_state
	generate_next_state(capsid)	            
	return configurations

def get_distinct_structures(configurations):
	#construct connection Graphs from adjacency matrices
	Graphs = []
	for conf in configurations:
	    n_hubs = len(conf.hubs)
	    A = np.zeros(shape=(n_hubs, n_hubs))
	    for ix,hub1 in enumerate(conf.hubs):
	        for jx,hub2 in enumerate(conf.hubs):
	            A[ix, jx] = len(get_edges_between_hubs(hub1, hub2))
	    Graphs.append(nx.from_numpy_matrix(A))

	#test for isomorphism between the graphs
	isomorphic = np.zeros(shape=(len(Graphs), len(Graphs)))
	for ix in range(len(Graphs)):
	    g1 = Graphs[ix]
	    hub_degree_ix = sorted([len(h.vertices) for h in configurations[ix].hubs])    
	    for jx in range(ix+1, len(Graphs)):
	        g2 = Graphs[jx]     
	        hub_degree_jx = sorted([len(h.vertices) for h in configurations[jx].hubs])
	        
	        filter_condition = (len(configurations[ix].hubs)==len(configurations[jx].hubs)) &\
	                            (hub_degree_ix==hub_degree_jx)
	                            #(len(configurations[ix].get_vertices())==len(configurations[jx].get_vertices())) &\
	                           #(len(configurations[ix].get_surface_hubs())==len(configurations[jx].get_surface_hubs()))
	        #filter_condition = nx.could_be_isomorphic(g1, g2)

	        if (filter_condition):
	            isomorphic[ix, jx]=int(nx.is_isomorphic(g1, g2))
	            isomorphic[jx, ix]=isomorphic[ix, jx]
	            #??????????? OK ????
	            if (isomorphic[ix, jx]==1):
	                #print ix, jx
	                break    

	#construct a graph using the isomorphic matrix as adjacency matrix
	isomorphism_graph = nx.from_numpy_matrix(isomorphic)
	#connected subgraphs here are all isomorphic
	distinct_structures = list(nx.connected_component_subgraphs(isomorphism_graph))	      

	return [list(structure.nodes) for structure in distinct_structures]          


def cool_slowly(capsid, n_cool_steps=10000):
	RNG = np.random.RandomState()
	kT = capsid.kT
	dT = kT/float(n_cool_steps)
	#cool off slowly
	for t in range(n_cool_steps):
	    pr.d_max = 0.01
	    for i in range(5*len(capsid.hubs)):
	        mc.attempt_move(capsid, RNG)
	    capsid.kT-=dT	

	#refine to get to the local minimum
	for dm in range(1,10):
	    pr.d_max/=(dm*10)
	    for t in range(10000):
	        for i in range(len(capsid.hubs)):
	            mc.attempt_move_athermal(capsid, RNG)


def construct_hessian(capsid, dr=1e-5):
	n_hubs = len(capsid.hubs)
	H = np.zeros(shape=(3*n_hubs, 3*n_hubs))
	H_diag = np.zeros(shape=(3*n_hubs, 3*n_hubs))
	H_test2 = np.zeros(shape=(3*n_hubs, 3*n_hubs))
	h_s = []
	#will need to loop over neighbor pairs
	hub_i = capsid.hubs[0]
	hub_j = hub_i#hub_i.get_neighbor_hubs()[2]
	#dr = 1e-5
	n_dr = 10

	#can be made faster by looping through the neighboring hubs only which "interact"
	for ix, hub_i in enumerate(capsid.hubs):
	    for jx, hub_j in enumerate(capsid.hubs):
	        if hub_i <> hub_j:
	            #print ix, jx, len(get_common_neighbor_hubs(hub_i, hub_j))
	            #H elements zero for these
	            if len(get_common_neighbor_hubs(hub_i, hub_j))==0:
	                continue
	            if len(get_common_neighbor_hubs(hub_i, hub_j))==1:
	                if hub_i not in hub_j.get_neighbor_hubs():
	                    continue                
	        
	        h = np.zeros(shape=(3,3))
	        h_test = np.zeros(shape=(3,3))
	        h_std = np.zeros(shape=(3,3))
	        #will need to loop over cartesian coordinates
	        for alpha in [0,1,2]:
	            for beta in [0,1,2]:
	                #diagonal elements of H
	                if (hub_i==hub_j) & (alpha==beta):
	                    continue
	                E = []
	                xi_val = []
	                xj_val = []
	                for x_i in np.arange(-n_dr, n_dr+1)*dr:
	                    xi, yi, zi = np.roll([x_i, 0.0, 0.0], alpha)
	                    hub_i.move(xi, yi, zi)
	                    for x_j in np.arange(-n_dr, n_dr+1)*dr:
	                        xj, yj, zj = np.roll([x_j, 0.0, 0.0], beta)
	                        hub_j.move(xj, yj, zj)
	                        E.append(en.elastic_energy(capsid.subunits))
	                        xi_val.append(x_i)
	                        xj_val.append(x_j)
	                        hub_j.move(-xj, -yj, -zj)    
	                    hub_i.move(-xi, -yi, -zi)    
	                xi_val = np.array(xi_val)
	                xj_val = np.array(xj_val)

	                Er = np.array(E).reshape(2*n_dr+1, 2*n_dr+1)
	                Dx, Dy = np.gradient(Er, dr)
	                D2x, D2xy = np.gradient(Dx, dr)
	                D2yx, D2y = np.gradient(Dy, dr)
	                #if hub_i<>hub_j:
	                    #if alpha<>beta:
	                h[alpha, beta]=np.mean(D2xy[2:-2,:2:-2])
	                
	                #h_std[alpha, beta]=np.std(D2xy) / np.mean(D2xy)
	                #DOUBLE TRIPLE CHECK THIS
	                if hub_i==hub_j:
	                    H_diag[3*ix+alpha, 3*ix+alpha] = np.mean(D2x[2:-2,:])
	                    #print 0.5*np.mean(D2x[2:-2,:]), 3*ix+alpha, 3*ix+alpha
	                    #H_diag[3*ix+alpha, 3*ix+alpha] = np.mean(D2y[:,2:-2])
	        H[3*ix:3*(ix+1), 3*jx:3*(jx+1)]=h
	        h_s.append(h)
	            #else:
	            #    h[alpha, beta]=np.mean(D2x[2:-2,:])
	            #    h_test[alpha,beta] = np.mean(D2y[:,2:-2])
	            #    h_std[alpha, beta]=np.std(D2x[2:-2,:]) / np.mean(D2x[2:-2,:])                  
	        #else:
	        #    h[alpha, beta]=np.mean(D2x) ##may
	H[np.diag_indices_from(H)] = np.diag(H_diag)
	return H

#def get_eigen