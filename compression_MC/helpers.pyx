# cython: profile=False

#import numpy as np
from Capsid import *
from numpy import cross, eye, dot
import parameters as pr
import numpy as np
#from math import exp
#from geometry import *
cdef extern from "math.h":
    double sqrt(double x)
    double fabs(double x) 


def com_position(subunit):
    """Returns the center of mass coordinates of a given subunit computed as the average position of vertices"""
    comx = np.mean([vertex.x for vertex in subunit.vertices])
    comy = np.mean([vertex.y for vertex in subunit.vertices])
    comz = np.mean([vertex.z for vertex in subunit.vertices])
    return (comx, comy, comz)

#faster than the generic com_position
def com_position_triangle(subunit):
    """Same as com_position(), but faster and only works for triangle subunits"""
    v0, v1, v2 = subunit.vertices
    cdef double comx = 1.0/3.0*(v0.x+v1.x+v2.x)
    cdef double comy = 1.0/3.0*(v0.y+v1.y+v2.y)
    cdef double comz = 1.0/3.0*(v0.z+v1.z+v2.z)
    return (comx, comy, comz)

#once it works, change to a general subunit [v1, v2, ... vn]
def create_triangle_subunit_from_vertices(v1, v2, v3, subunit_type):
    """Returns a triangle subunit initialized by vertices given as parameters. 
    subunit_type should be listed in the global parameter dictionary which serves as a prototype for the subunit."""
    edges = pr.params['subunits'][subunit_type]['edges']
    number_of_edges = len(edges)
    e12 = Edge(v1,v2, 1.0)
    e12.edge_type = edges['edge1']
    v1.edge_out = e12
    v2.edge_in = e12
    e23 = Edge(v2,v3, 1.0)
    e23.edge_type = edges['edge2']
    v2.edge_out = e23
    v3.edge_in = e23
    e31 = Edge(v3,v1, 1.0)
    e31.edge_type = edges['edge3']
    v3.edge_out = e31
    v1.edge_in = e31
    sub = Subunit()
    sub.vertices = [v1, v2, v3]
    sub.edges = [e12, e23, e31]

    #!!!other type of excluders based on the subunit_type should be added here
    ##exc = Excluder(R=pr.excluder_radius, get_position=com_position_triangle, excluder_type='com')
    ##sub.excluders = [exc]
    sub.excluders = []
    for excluder in pr.params['subunits'][subunit_type]['excluders']:
        excluder_type = pr.params['subunits'][subunit_type]['excluders'][excluder]
                                                                          #!!! will need a fix to get position function from parameters or lookup      
        exc = Excluder(R=pr.params['excluders'][excluder_type]['radius'], get_position=com_position_triangle, excluder_type=excluder_type)
        exc.overlaps = pr.params['excluders'][excluder_type]['overlaps']
        exc.parent_subunit = sub
        exc.update_position()
        sub.excluders.append(exc)
    for v in [v1, v2, v3]:
        v.parent_subunit = sub
        h = Hub()
        h.add_vertex(v)
    for e in [e12, e23, e31]:
        e.parent_subunit = sub
    for edge in sub.edges:
        for key in pr.params['edges'][edge.edge_type]:
            setattr(edge, key, pr.params['edges'][edge.edge_type][key])

    sub.subunit_type = subunit_type
    sub.mu = pr.params['subunits'][subunit_type]['mu']
    #should be inferred from pr.params['subunits'][subunit_type]['edges'] ?
    sub.number_of_rotational_configurations = pr.params['subunits'][subunit_type]['number_of_rotational_configurations']#3

    return sub    


#wrapper
def create_subunit_from_vertices(v1, v2, v3, subunit_type):
    """Wrapper function to include other subunits than triangles at some point"""
    return create_triangle_subunit_from_vertices(v1, v2, v3, subunit_type)


def create_n_gon(n, R=1.0, subunit_type='sub1'):
    """Returns a capsid formed of n subunits sharing a common vertex and the others lying on a circle of radius R
    n=6 gives a hexagon, n=5 a pentagon, etc."""
    x = []
    y = []
    da = 2.0*np.pi/n
    for i in xrange(n):
        x.append(R*np.cos(i*da))
        y.append(R*np.sin(i*da))
    
    vc = Vertex(0, 0, 0, None, None)
    v0 = Vertex(x[0], y[0], 0, None, None)
    v1 = Vertex(x[1], y[1], 0, None, None)
    v0_first = v0
    vc_first = vc
    sub1 = create_triangle_subunit_from_vertices(vc, v1, v0, subunit_type)
    capsid = Capsid()
    capsid.add_new_subunit(sub1)
    
    for i in xrange(1, n):
        v1_prev = v1
        vc_prev = vc
        vc = Vertex(0, 0, 0, None, None)
        v0 = Vertex(x[i], y[i], 0, None, None)
        v1 = Vertex(x[(i+1)%n], y[(i+1)%n], 0, None, None)
        sub = create_triangle_subunit_from_vertices(vc, v1, v0, subunit_type)
        capsid.add_new_subunit(sub)
        capsid.connect_subunit(vc_prev.parent_hub, v1_prev.parent_hub, vc.parent_hub, v0.parent_hub)
    
    #capsid.connect_subunit(vc.parent_hub, v1.parent_hub, vc_first.parent_hub, v0_first.parent_hub)    
    v1.parent_hub.merge_into(v0_first.parent_hub)
    return capsid


def create_wedge_template():
    """Returns a capsid formed of 5 subunits with a wedge"""
    n=6
    R=1.0
    x = []
    y = []
    da = 2.0*np.pi/n
    for i in xrange(n):
        x.append(R*np.cos(i*da))
        y.append(R*np.sin(i*da))
    
    vc = Vertex(0, 0, 0, None, None)
    v0 = Vertex(x[0], y[0], 0, None, None)
    v1 = Vertex(x[1], y[1], 0, None, None)
    v0_first = v0
    vc_first = vc
    sub1 = create_triangle_subunit_from_vertices(vc, v1, v0, subunit_type="sub0")
    capsid = Capsid()
    capsid.add_new_subunit(sub1)
    
    for i in xrange(1, n-2):
        v1_prev = v1
        vc_prev = vc
        vc = Vertex(0, 0, 0, None, None)
        v0 = Vertex(x[i], y[i], 0, None, None)
        v1 = Vertex(x[(i+1)%n], y[(i+1)%n], 0, None, None)
        sub = create_triangle_subunit_from_vertices(vc, v1, v0, subunit_type="sub2")
        capsid.add_new_subunit(sub)
        capsid.connect_subunit(vc_prev.parent_hub, v1_prev.parent_hub, vc.parent_hub, v0.parent_hub)
    
    for i in xrange(n-2, n-1):
        v1_prev = v1
        vc_prev = vc
        vc = Vertex(0, 0, 0, None, None)
        v0 = Vertex(x[i], y[i], 0, None, None)
        v1 = Vertex(x[(i+1)%n], y[(i+1)%n], 0, None, None)
        sub = create_triangle_subunit_from_vertices(vc, v1, v0, subunit_type="sub0")
        capsid.add_new_subunit(sub)
        capsid.connect_subunit(vc_prev.parent_hub, v1_prev.parent_hub, vc.parent_hub, v0.parent_hub)


    #capsid.connect_subunit(vc.parent_hub, v1.parent_hub, vc_first.parent_hub, v0_first.parent_hub)    
    #v1.parent_hub.merge_into(v0_first.parent_hub)
    return capsid


def create_wedge_template2():
    """Returns a capsid formed of 6 subunits not closed at the end"""
    n=6
    R=1.0
    x = []
    y = []
    da = 2.0*np.pi/n
    for i in xrange(n):
        x.append(R*np.cos(i*da))
        y.append(R*np.sin(i*da))
    
    vc = Vertex(0, 0, 0, None, None)
    v0 = Vertex(x[0], y[0], 0, None, None)
    v1 = Vertex(x[1], y[1], 0, None, None)
    v0_first = v0
    vc_first = vc
    sub1 = create_triangle_subunit_from_vertices(vc, v1, v0, subunit_type="sub0")
    capsid = Capsid()
    capsid.add_new_subunit(sub1)
    
    for i in xrange(1, n-2):
        v1_prev = v1
        vc_prev = vc
        vc = Vertex(0, 0, 0, None, None)
        v0 = Vertex(x[i], y[i], 0, None, None)
        v1 = Vertex(x[(i+1)%n], y[(i+1)%n], 0, None, None)
        sub = create_triangle_subunit_from_vertices(vc, v1, v0, subunit_type="sub2")
        capsid.add_new_subunit(sub)
        capsid.connect_subunit(vc_prev.parent_hub, v1_prev.parent_hub, vc.parent_hub, v0.parent_hub)
    
    for i in xrange(n-2, n-1):
        v1_prev = v1
        vc_prev = vc
        vc = Vertex(0, 0, 0, None, None)
        v0 = Vertex(x[i], y[i], 0, None, None)
        v1 = Vertex(x[(i+1)%n], y[(i+1)%n], 0, None, None)
        sub = create_triangle_subunit_from_vertices(vc, v1, v0, subunit_type="sub0")
        capsid.add_new_subunit(sub)
        capsid.connect_subunit(vc_prev.parent_hub, v1_prev.parent_hub, vc.parent_hub, v0.parent_hub)

    for i in xrange(n-1, n):
        v1_prev = v1
        vc_prev = vc
        vc = Vertex(0, 0, 0, None, None)
        v0 = Vertex(x[i], y[i], 0, None, None)
        v1 = Vertex(x[(i+1)%n], y[(i+1)%n], 0, None, None)
        sub = create_triangle_subunit_from_vertices(vc, v1, v0, subunit_type="sub1")
        capsid.add_new_subunit(sub)
        capsid.connect_subunit(vc_prev.parent_hub, v1_prev.parent_hub, vc.parent_hub, v0.parent_hub)


    #capsid.connect_subunit(vc.parent_hub, v1.parent_hub, vc_first.parent_hub, v0_first.parent_hub)    
    #v1.parent_hub.merge_into(v0_first.parent_hub)
    return capsid


def create_single_subunit_capsid(subunit_type='sub1'):
    """Returns a capsid formed of a single subunit."""
    v0 = Vertex(0, 0, 0, None, None)
    v1 = Vertex(0, 1, 0, None, None)   
    v2 = Vertex(3.0**0.5/2.0, 0.5, 0, None, None)    
    sub = create_subunit_from_vertices(v0,v1,v2, subunit_type)
    capsid = Capsid()
    capsid.add_new_subunit(sub)
    return capsid    

def create_hex_sheet(size1, size2):
    """Create a flat sheet of triangle subunits arranged in a triangular lattice."""
    hubs = []
    hexagons = []
    for k in xrange(size1):
        for l in xrange(size2):
            hexagon = create_n_gon(6)
            for hub in hexagon.hubs:
                for vertex in hub.vertices:
                    vertex.x+=0.5*(3*(k+l))
                    vertex.y+=0.5*(3.0**0.5*(k-l))
            hubs = hubs + hexagon.hubs
            hexagons.append(hexagon)
    capsid = Capsid()
    for hexagon in hexagons:
        for subunit in hexagon.subunits:
                capsid.add_new_subunit(subunit)
    #hubs were added multiple times
    capsid.hubs = list(set(capsid.hubs))
    #has to be a while loop; while there are close hubs, keep merging
    is_there_to_merge = True
    while is_there_to_merge:
        is_there_to_merge = False
        for hub1 in hubs:
            for hub2 in hubs:
                if hub1==hub2:
                    continue
                if distance(hub1.vertices[0], hub2.vertices[0])<1e-4:
                    v2 = hub2.vertices[0]
                    hub2.merge_into(hub1)
                    hubs.remove(hub2)
                    hubs.append(v2.parent_hub)
                    #print (hub2 in capsid.hubs)
                    is_there_to_merge = True
                    break
            if (is_there_to_merge):
                break
    return capsid


def number_of_edges_between_hubs(hub1, hub2):
    """Given two hubs, returns the number of edges between them. This should be 0, 1 or 2"""
    n_edges_between_hubs = 0
    for vertex in hub1.vertices:
        n_edges_between_hubs +=int(vertex.edge_out.vertex_to.parent_hub ==hub2)
    for vertex in hub2.vertices:
        n_edges_between_hubs +=int(vertex.edge_out.vertex_to.parent_hub == hub1)
    return n_edges_between_hubs

def get_edges_between_hubs(hub1, hub2):
    """Given two hubs, returns the edges between them. There should be 0,1 or 2 edges between any two hubs"""
    edges_between_hubs = []
    for vertex in hub1.vertices:
        if (vertex.edge_out.vertex_to.parent_hub ==hub2):
            edges_between_hubs.append(vertex.edge_out)
    for vertex in hub2.vertices:
        if (vertex.edge_out.vertex_to.parent_hub == hub1):
            edges_between_hubs.append(vertex.edge_out)            
    return edges_between_hubs

def get_common_neighbor_hubs(hub1, hub2):
    """Given two hubs, returns a list of hubs that are neighbors (i.e. connected to) of both hubs."""
    hub1_neighbors = hub1.get_neighbor_hubs()
    hub2_neighbors = hub2.get_neighbor_hubs()
    return list(set(hub1_neighbors) & set(hub2_neighbors))

def get_type1_fusion_pairs(capsid):
    """Type 1 fusion is wedge fusion. Function returns the possible hub pairs that can be fused.
    These are the hubs that share a single common neighbor hub and both of them are on the surface."""
    #wedge closure fusion
    #a hub is on surface if it has two surface edges
    ##surface_hubs = capsid.get_surface_hubs()
    ##fusion_pairs = []
    ##for hub1 in surface_hubs:
    ##    for hub2 in surface_hubs:
    ##        if hub1==hub2:
    ##            continue
    ##        #should check if hub1 and hub2 are not connected neighbors!!
    ##        if number_of_edges_between_hubs(hub1, hub2)>0:
    ##            continue
    ##        #type 1 fusion is wedge closure; need a single common neighbor 
    ##        common_neighbors = get_common_neighbor_hubs(hub1, hub2)
    ##        if len(common_neighbors)<>1:
    ##            continue
##
    ##        common_neighbor = common_neighbors[0]
    ##        #both edges have to be surface edges
    ##        if number_of_edges_between_hubs(hub1, common_neighbor)>1:
    ##            continue
    ##        if number_of_edges_between_hubs(hub2, common_neighbor)>1:
    ##            continue                
    ##        
    ##        #if distance(hub1.vertices[0], hub2.vertices[0])<pr.l_fusion:
    ##        if within_cube(hub1.vertices[0], hub2.vertices[0], pr.l_fusion):
    ##            fusion_pairs.append([hub1, hub2])

    #might be a little slower this way because has to build the full wedge list (as opposed to the actual pairs only)
    wedge_pairs = get_open_wedge_pairs(capsid)

    fusion_pairs = filter(lambda h: within_cube(h[0].vertices[0], h[1].vertices[0], pr.l_fusion) ,wedge_pairs)

    #otherwise subunit will fuse into plane
    fusion_pairs = filter(lambda h: len(h[2].vertices)>2 ,fusion_pairs)

    #DB! ---
    #fusion_pairs = filter(lambda h: len(h[2].vertices)>5 ,fusion_pairs)
    #--- DB
    return fusion_pairs


def get_open_wedge_pairs(capsid):    
    """Basically the same as get_type1_fusion_pairs, except that it doesn't test for the distance between the fusion candidate hubs.
    Function not entirely correct: in the case of 2 triangle subunits, there are two closeable wedges and returns none because can't choose between the two."""
    surface_hubs = capsid.get_surface_hubs()
    wedge_pairs = []
    #for hub1 in surface_hubs:
    #    for hub2 in surface_hubs:
    for i1 in range(len(surface_hubs)):
        for i2 in range(i1, len(surface_hubs)):
            hub1 = surface_hubs[i1]
            hub2 = surface_hubs[i2]
            if hub1==hub2:
                continue
            #should check if hub1 and hub2 are not connected neighbors!!
            if number_of_edges_between_hubs(hub1, hub2)>0:
                continue
            #type 1 fusion is wedge closure; need a single common neighbor 
            common_neighbors = get_common_neighbor_hubs(hub1, hub2)
            #print "CN: ", common_neighbors
            if len(common_neighbors)==0:
                continue
            if len(common_neighbors)<>1:
                continue
            common_neighbor = None
            for cn in common_neighbors:
                if ((number_of_edges_between_hubs(hub1, cn)==1) & (number_of_edges_between_hubs(hub2, cn)==1)):
                    common_neighbor = cn	 

            #@common_neighbor = common_neighbors[0]
            #both edges have to be surface edges
            #@if number_of_edges_between_hubs(hub1, common_neighbor)>1:
            #@    continue
            #@if number_of_edges_between_hubs(hub2, common_neighbor)>1:
            #@    continue                
            
            #if distance(hub1.vertices[0], hub2.vertices[0])<pr.l_fusion:
            #if within_cube(hub1.vertices[0], hub2.vertices[0], pr.l_fusion):
            if common_neighbor is not None:
                wedge_pairs.append([hub1, hub2, common_neighbor])
    
    return wedge_pairs




def get_type2_fusion_pairs(capsid):
    """Type 2 fusion is crack fusion. Function returns the possible hub pairs that can be fused.
    These are surface hubs sharing two common neighbors."""
    #crack closing fusion
    #a hub is on surface if it has two surface edges
    surface_hubs = capsid.get_surface_hubs()
    fusion_pairs = []
    for hub1 in surface_hubs:
        for hub2 in surface_hubs:
            if hub1==hub2:
                continue
            #should check if hub1 and hub2 are not connected neighbors!!
            if number_of_edges_between_hubs(hub1, hub2)>0:
                continue
            #type 2 fusion is crack closure; need 2 common neighbors
            common_neighbors = get_common_neighbor_hubs(hub1, hub2)
            if len(common_neighbors)<>2:
                continue
            #also all involved edges have to be surface edges       
            nedges1 = sum(map(lambda x: number_of_edges_between_hubs(hub1, x), common_neighbors))
            nedges2 = sum(map(lambda x: number_of_edges_between_hubs(hub2, x), common_neighbors))
            if nedges1<>2:
                continue
            if nedges2<>2:
                continue

            #if distance(hub1.vertices[0], hub2.vertices[0])<pr.l_fusion:
            if within_cube(hub1.vertices[0], hub2.vertices[0], pr.l_fusion):
                fusion_pairs.append([hub1, hub2])
    
    return fusion_pairs

def get_type1_fission_pairs(capsid):
    """Type 1 fission is wedge opening. Function returns pairs of hubs [(hub1, hub2),...] defining edges that can undergo type 1 fission. 
    hub1 is assumed to be the surface hub which can break, opening a wedge."""
    #wedge opening
    fission_pairs = []
    #fission moves opening up a wedge
    surface_hubs = capsid.get_surface_hubs()
    for hub in surface_hubs:
        if len(hub.vertices)==1:
            continue
        neighbor_hubs = hub.get_neighbor_hubs()
        for neighbor in neighbor_hubs:
            if neighbor in surface_hubs:
                continue
           #DB!!! -----
           # if 'sub1' not in [v.parent_subunit.subunit_type for v in hub.vertices]:
           #     continue
           #DB ------- 
            fission_pairs.append((hub, neighbor))
    return fission_pairs
    
        
def get_type2_fission_triplets(capsid):
    """Type 2 fission is crack opening. Function returns triplets of hubs [(hub1, hub2, hub3), ...] defining 2 edges that can undergo
    type 2 fission. hub2 is assumed to split, opening a (hub2->hub1) and a (hub2->hub3) crack."""

    #crack opening
    fission_triplets = []
    surface_hubs = capsid.get_surface_hubs()
    for hub in capsid.hubs:
        #crack cannot open up from surface hubs
        if hub in surface_hubs:
            continue
        neighbors = hub.get_neighbor_hubs()
        #non-surface neighbors
        deep_neighbors = filter(lambda x: x not in surface_hubs, neighbors)
        for neigh1 in deep_neighbors:
            for neigh2 in deep_neighbors:
                if neigh1 <> neigh2:
                    fission_triplets.append((neigh1, hub, neigh2))
    #each move is added twice; should remove symmetric occurences
    return fission_triplets

def get_edge_fusion_hubs(capsid):
    """Edge fusion is two unconnected edges merging. Function returns list of four hubs [[hub1A, hub1B, hub2A, hub2B],...]
    defining edges A and B. hub1A will merge into hub1B, and hub2A will merge into hub2B"""
    edge_fusion_hubs = []
    surface_edges = filter(lambda x: x.is_surface_edge(), capsid.get_edges())
    for edgeA in surface_edges:
        hub1A = edgeA.vertex_from.parent_hub
        hub2A = edgeA.vertex_to.parent_hub
        for edgeB in surface_edges:
            hub2B = edgeB.vertex_from.parent_hub
            hub1B = edgeB.vertex_to.parent_hub 
            if (hub1A==hub1B):
                continue
            if (hub2A==hub2B):
                continue           


            common_hubs1 = get_common_neighbor_hubs(hub1A, hub1B)
            common_hubs2 = get_common_neighbor_hubs(hub2A, hub2B)
            if len(common_hubs1)+len(common_hubs2)>0:
                continue #then it's a wedge or crack fusion
            else:
                #check distances
                #d1 = distance(hub1A.vertices[0], hub1B.vertices[0])
                #d2 = distance(hub2A.vertices[0], hub2B.vertices[0])
                #if ((d1<pr.l_fusion) & (d2<pr.l_fusion)):
                d1 = within_cube(hub1A.vertices[0], hub1B.vertices[0], pr.l_fusion)
                d2 = within_cube(hub2A.vertices[0], hub2B.vertices[0], pr.l_fusion)
                if (d1 & d2):
                    edge_fusion_hubs.append([hub1A, hub1B, hub2A, hub2B])

    return edge_fusion_hubs


def get_edge_fission_hubs(capsid):
    """Reverse of edge fusion. Returns hub pairs defining edges that can be splitted by this move."""
    #edge fission can happen on edges which are not on surface, 
    #but both their hubs are
    edge_fission_hubs = []
    non_surface_edges = filter(lambda x: not x.is_surface_edge(), capsid.get_edges())
    for edge in non_surface_edges:
        hub1 = edge.vertex_from.parent_hub
        hub2 = edge.vertex_to.parent_hub
        if (hub1.is_surface_hub() & hub2.is_surface_hub()):
            if (hub2, hub1) not in edge_fission_hubs:
                edge_fission_hubs.append((hub1, hub2))
    disconnecting_moves = []
    #now have to make sure that edge fission would not disconnect the graph
    for hub1, hub2 in edge_fission_hubs:
        next_hub = hub1
        for i in xrange(1000): #should be while loop
            edges_out = [vertex.edge_out for vertex in next_hub.vertices]
            surface_edges_out = filter(lambda x: x.is_surface_edge(),edges_out)
            if len(surface_edges_out)<>1:
                raise ValueError('Edge fission screwed up #1!')
            if i>900:
                raise ValueError('Edge fission screwed up #2!')                
            next_hub = surface_edges_out[0].vertex_to.parent_hub
            #could make a loop to arrive back to hub1
            if next_hub==hub1:
                break
            #could make a loop to arrive to hub2; disconnected domain
            if next_hub==hub2:
                disconnecting_moves.append((hub1, hub2))
                break
    return [(hub1, hub2) for (hub1, hub2) in edge_fission_hubs
            if (hub1, hub2) not in disconnecting_moves]


def get_edge_fission_hubs_breakup_allowed(capsid):
    """Reverse of edge fusion. Returns hub pairs defining edges that can be splitted by this move."""
    #edge fission can happen on edges which are not on surface, 
    #but both their hubs are
    edge_fission_hubs = []
    non_surface_edges = filter(lambda x: not x.is_surface_edge(), capsid.get_edges())
    for edge in non_surface_edges:
        hub1 = edge.vertex_from.parent_hub
        hub2 = edge.vertex_to.parent_hub
        if (hub1.is_surface_hub() & hub2.is_surface_hub()):
            if (hub2, hub1) not in edge_fission_hubs:
                edge_fission_hubs.append((hub1, hub2))
    #moves that would break off single subunits
    forbidden_moves = []
    #now have to make sure that edge fission would not disconnect single subunits
    
    for hub1, hub2 in edge_fission_hubs:
        visited_subunits = []
        next_hub = hub1
        for i in xrange(1000): #should be while loop
            edges_out = [vertex.edge_out for vertex in next_hub.vertices]
            surface_edges_out = filter(lambda x: x.is_surface_edge(),edges_out)
            visited_subunits.append(surface_edges_out[0].parent_subunit)
            if len(surface_edges_out)<>1:
                raise ValueError('Edge fission screwed up #1!')
            if i>900:
                raise ValueError('Edge fission screwed up #2!')                
            next_hub = surface_edges_out[0].vertex_to.parent_hub
            #could make a loop to arrive back to hub1
            if next_hub==hub1:
                break
            #could make a loop to arrive to hub2; disconnected domain
            if next_hub==hub2:
                #disconnecting_moves.append((hub1, hub2))
                #a piece will break off
                visited_subunits = list(set(visited_subunits))
                if len(visited_subunits)==1:
                    forbidden_moves.append((hub1, hub2))
                if (len(capsid.subunits)-len(visited_subunits))==1:   
                    forbidden_moves.append((hub1, hub2))     

                #don't allow for breakoff until after nucleation. this has to go in run.py
                if len(capsid.subunits)<30:
                    forbidden_moves.append((hub1, hub2))                                     
                break
    return [(hub1, hub2) for (hub1, hub2) in edge_fission_hubs
            if (hub1, hub2) not in forbidden_moves]


def get_removable_subunits(capsid):
    """Returns subunits which are only connected along one edge and thus can be removed. Does not select between subunit types."""
    #find subunits that can be removed. these are the ones that are attached along one edge only
    removable_subunits = []
    for subunit in capsid.subunits:
        edges = subunit.edges
        vertices = subunit.vertices
        number_of_surface_edges = len(filter(lambda x: x.is_surface_edge(), edges))
        if (number_of_surface_edges == (len(vertices)-1)):
            removable_subunits.append(subunit)
    if len(capsid.subunits)<2: #don't remove last subunit
        return []
    return removable_subunits

def get_removable_wedge_subunits(capsid):
    removable_subunits = []
    for subunit in capsid.subunits:
        edges = subunit.edges
        vertices = subunit.vertices
        surface_edges = filter(lambda x: x.is_surface_edge(), edges)
        number_of_surface_edges = len(surface_edges)
        if (number_of_surface_edges == 1):
            #have to test to make sure structure doesn't fall apart by removing the wedge, eg. /\/\/\/\ zigzag stripe
            #looks like this is only possible if wedge hub is a non-surface hub
            # a triangle subunit with a single free edge thus can have at most 1 inner hub

            if sum([(not v.parent_hub.is_surface_hub()) for v in vertices])==1:

                removable_subunits.append(subunit)
    if len(capsid.subunits)<2: #don't remove last subunit
        return []
    return removable_subunits



def update_excluders_position(affected_subunits):
    """Excluder positions are updated after every move for the affected subunits. This saves time from recalculating all of them unless subunit vertices were moved."""
    excluders = map(lambda x: x.excluders, affected_subunits)
    excluders_flattened = [y for x in excluders for y in x]
    map(lambda x: x.update_position(), excluders_flattened)
    return 



def update_neighbor_lists(capsid):
    for sub in capsid.subunits:
        sub.neighbor_list = []
    for sub1 in capsid.subunits:
        for sub2 in capsid.subunits:
            if sub1==sub2:
                continue
            xA, yA, zA = com_position_triangle(sub1)
            xB, yB, zB = com_position_triangle(sub2)
            L = 1.2*sub1.edges[0].equilibrium_length              
            if ((fabs(xA-xB)<L) & (fabs(yA-yB)<L) & (fabs(zA-zB)<L)):
                sub1.neighbor_list.append(sub2)


class UmbrellaWindow(object):
    #no idea what delta was for
    def __init__(self, left=0, right=4, delta=1):
        self.left=left
        self.right = right
        self.delta = delta
        self.values = []
        self.init_conf = None
    def add_value(self, v):
        self.values.append(v)
    def shift(self,dx):
        self.left+=dx
        self.right+=dx
    #forward compatibility    
    def bias_potential(self, N):
        return 0.0


class HarmonicUmbrellaWindow(object):
    def __init__(self, spring_const, N0):
        self.spring_const=spring_const
        self.N0 = N0
        self.values = []

        #for backward compatibility with UmbrellaWindow
        self.left = 0
        self.right = 10000
        self.init_conf=None

    def bias_potential(self, N):
        return 0.5*self.spring_const*(N-self.N0)**2
            
    def add_value(self, v):
        self.values.append(v)


class Wall(object):
    def __init__(self, hardness, z_position, amplitude=1.0):
        self.hardness = hardness
        self.z_position = z_position
        self.amplitude = amplitude
    def potential_energy(self, z):
        """Second wall assumed at z=0"""
        #return np.exp(-self.hardness*(z-self.z_position))/abs((self.z_position-z)) + np.exp(-self.hardness*(z))/abs(z)
        return self.amplitude*(np.exp(self.hardness*(z-self.z_position)) + np.exp(-self.hardness*z))
    def move(self, dz):
        self.z_position+=dz


def create_tubule(Nperim, L):
    import copy
    l0=1.0
    theta = 2.0*np.pi / Nperim    
    R = l0/(2.0*np.sin(theta/2))
    dz = l0*(3.0**0.5)/2
    capsid = Capsid()

    def create_first_crown():

        #first crown    
        xn = []
        yn = []

        xm = []
        ym = []

        
        for n in range(Nperim):
            theta_n = n*theta
            xn.append(R*np.cos(theta_n))
            yn.append(R*np.sin(theta_n))

            theta_m = (n-1)*theta + theta/2
            xm.append(R*np.cos(theta_m))
            ym.append(R*np.sin(theta_m))    

        for n in range(1, Nperim+1):

            nx = n % Nperim
            v1 = Vertex(xn[nx-1], yn[nx-1], 0.0 , None, None)
            v2 = Vertex(xn[nx], yn[nx], 0.0 , None, None)    
            v3 = Vertex(xm[nx], ym[nx], dz , None, None)   
            sub = create_subunit_from_vertices(v1, v2, v3, 'sub1')
            capsid.add_new_subunit(sub)


    def create_second_crown():        
        #second (flipped) crown
        xn = []
        yn = []

        xm = []
        ym = []


        for n in range(Nperim):
            theta_n = n*theta - theta/2
            xn.append(R*np.cos(theta_n))
            yn.append(R*np.sin(theta_n))

            theta_m = (n-1)*theta
            xm.append(R*np.cos(theta_m))
            ym.append(R*np.sin(theta_m))    

        for n in range(1, Nperim+1):

            nx = n % Nperim
            v1 = Vertex(xn[nx-1], yn[nx-1], dz , None, None)
            v2 = Vertex(xn[nx], yn[nx], dz , None, None)    
            v3 = Vertex(xm[nx], ym[nx], 0.0 , None, None)   
            #check orientation in second crown
            sub = create_subunit_from_vertices(v2, v1, v3, 'sub1')
            capsid.add_new_subunit(sub)

    def create_third_crown():    

        #third(shifted) crown    
        xn = []
        yn = []

        xm = []
        ym = []


        for n in range(Nperim):
            theta_n = n*theta + theta/2
            xn.append(R*np.cos(theta_n))
            yn.append(R*np.sin(theta_n))

            theta_m = (n)*theta
            xm.append(R*np.cos(theta_m))
            ym.append(R*np.sin(theta_m))    

        for n in range(1, Nperim+1):

            nx = n % Nperim
            v1 = Vertex(xn[nx-1], yn[nx-1], dz , None, None)
            v2 = Vertex(xn[nx], yn[nx], dz , None, None)    
            v3 = Vertex(xm[nx], ym[nx], 2*dz , None, None)   
            sub = create_subunit_from_vertices(v1, v2, v3, 'sub1')
            capsid.add_new_subunit(sub)    

    def create_fourth_crown():    

        #fourth(shifted and flipped) crown    
        xn = []
        yn = []

        xm = []
        ym = []


        for n in range(Nperim):
            theta_n = n*theta #- theta/2
            xn.append(R*np.cos(theta_n))
            yn.append(R*np.sin(theta_n))

            theta_m = (n-1)*theta + theta/2
            xm.append(R*np.cos(theta_m))
            ym.append(R*np.sin(theta_m))    

        for n in range(1, Nperim+1):

            nx = n % Nperim
            v1 = Vertex(xn[nx-1], yn[nx-1], 2*dz , None, None)
            v2 = Vertex(xn[nx], yn[nx], 2*dz , None, None)    
            v3 = Vertex(xm[nx], ym[nx], dz , None, None)   
            sub = create_subunit_from_vertices(v2, v1, v3, 'sub1')
            capsid.add_new_subunit(sub)     
   
    def merge_hubs():
        hubs = capsid.hubs
        #has to be a while loop; while there are close hubs, keep merging
        is_there_to_merge = True
        while is_there_to_merge:
            hubs = capsid.hubs
            is_there_to_merge = False
            for hub1 in hubs:
                for hub2 in hubs:
                    if hub1==hub2:
                        continue
                    if distance(hub1.vertices[0], hub2.vertices[0])<1e-4:
                        v2 = hub2.vertices[0]
                        hub2.merge_into(hub1)
                        is_there_to_merge = True
                        break
                if (is_there_to_merge):
                    break

    
    create_first_crown()
    create_second_crown()
    if (L<>1):
        create_third_crown()
        create_fourth_crown()
    
    new_subs = []
    if L%2<>0:
        rn = range(L/2)
    else:
        rn = range(L/2-1)
    for l in rn:
        subs = capsid.subunits

        for sub in subs:
            new_sub = copy.deepcopy(sub)
            for v in new_sub.vertices:
                v.z+=2*dz*(l+1)
            new_subs.append(new_sub)    

    if L%2<>0:
        new_subs = new_subs[:-2*Nperim]
    for new_sub in new_subs:
        capsid.add_new_subunit(new_sub)
    merge_hubs()    
    return capsid    


def create_trumpet(Nperim, L, parameters):
    import monte_carlo as mc
    RNG = np.random.RandomState()
    capsid = create_tubule(Nperim, L)
    pr.update_capsid_parameters(capsid, parameters)    
    E_elast = []
    #nsteps = 1000#0
    #ndrop = 100#0
    mc.update_neighbor_lists(capsid)
    nsteps = parameters['global_params']['nsteps_quench']
    ndrop = nsteps/10
    pr.d_max = parameters['global_params']['d_max']
    for s in range(2*nsteps):  
        #pr.d_max=0.05#(2*1.0/50.0)**0.5
        #@if s%ndrop==0:
        #@    pr.d_max/=1.2
        accepted = []
        for i in range(len(capsid.hubs)):
            accepted.append(mc.attempt_move_athermal(capsid, RNG))
        #adaptive step; if too many rejections, decrease d_max   
        if float(sum(accepted)) / float(len(accepted)) < 0.5:
            pr.d_max/=1.1
        else:
            pr.d_max*=1.1    

        E_elast.append(mc.elastic_energy(capsid.subunits))
    return capsid, E_elast  


def create_trumpet_thermalize_quench(Nperim, L, parameters):
    import monte_carlo as mc    
    nsteps = parameters['global_params']['nsteps_quench']
    RNG = np.random.RandomState()
    capsid = create_tubule(Nperim, L)
    pr.update_capsid_parameters(capsid, parameters)    
    E_elast = []
    #nsteps = 1000#0
    #ndrop = 100#0
    mc.update_neighbor_lists(capsid)

    #thermalize
    capsid.kT = 1.0
    for i in range(10*len(capsid.hubs)):
        mc.attempt_move(capsid, RNG) 

    #quench           
    dTdrop = (capsid.kT-0.05)/nsteps
    for s in range(nsteps):
        for i in range(len(capsid.hubs)):
            mc.attempt_move(capsid, RNG) 
        capsid.kT-=dTdrop

    ndrop = nsteps/10
    pr.d_max = parameters['global_params']['d_max']
    for s in range(2*nsteps):  
        #pr.d_max=0.05#(2*1.0/50.0)**0.5
        if s%ndrop==0:
            pr.d_max/=1.2
        for i in range(len(capsid.hubs)):
            mc.attempt_move_athermal(capsid, RNG)
        E_elast.append(mc.elastic_energy(capsid.subunits))
    return capsid, E_elast      


def crack_and_quench(capsid, parameters):
    import monte_carlo as mc     
    nsteps = parameters['global_params']['nsteps_quench']
    RNG = np.random.RandomState()
    crack_edge(capsid)
    E_elast = []
    mc.update_neighbor_lists(capsid)
    capsid.kT = 1.0
    for i in range(50*len(capsid.hubs)):
        mc.attempt_move(capsid, RNG) 

    #quench           
    dTdrop = (capsid.kT-0.05)/nsteps
    for s in range(nsteps):
        for i in range(len(capsid.hubs)):
            mc.attempt_move(capsid, RNG) 
        capsid.kT-=dTdrop

    ndrop = nsteps/10
    pr.d_max = parameters['global_params']['d_max']
    for s in range(2*nsteps):  
        #pr.d_max=0.05#(2*1.0/50.0)**0.5
        if s%ndrop==0:
            pr.d_max/=1.2
        for i in range(len(capsid.hubs)):
            mc.attempt_move_athermal(capsid, RNG)
        E_elast.append(mc.elastic_energy(capsid.subunits))
    return capsid, E_elast        



def create_tubule_hanging_wedges(Nperim, L):
    """Hanging wedges are added to both rims to keep double edges at the rim for uniform stretching moduli"""
    import copy
    l0=1.0
    theta = 2.0*np.pi / Nperim    
    R = l0/(2.0*np.sin(theta/2))
    dz = l0*(3.0**0.5)/2
    capsid = Capsid()

    def create_first_crown():

        #first crown    
        xn = []
        yn = []

        xm = []
        ym = []

        
        for n in range(Nperim):
            theta_n = n*theta
            xn.append(R*np.cos(theta_n))
            yn.append(R*np.sin(theta_n))

            theta_m = (n-1)*theta + theta/2
            xm.append(R*np.cos(theta_m))
            ym.append(R*np.sin(theta_m))    

        for n in range(1, Nperim+1):

            nx = n % Nperim
            v1 = Vertex(xn[nx-1], yn[nx-1], 0.0 , None, None)
            v2 = Vertex(xn[nx], yn[nx], 0.0 , None, None)    
            v3 = Vertex(xm[nx], ym[nx], dz , None, None)   
            sub = create_subunit_from_vertices(v1, v2, v3, 'sub1')
            capsid.add_new_subunit(sub)


    def create_second_crown():        
        #second (flipped) crown
        xn = []
        yn = []

        xm = []
        ym = []


        for n in range(Nperim):
            theta_n = n*theta - theta/2
            xn.append(R*np.cos(theta_n))
            yn.append(R*np.sin(theta_n))

            theta_m = (n-1)*theta
            xm.append(R*np.cos(theta_m))
            ym.append(R*np.sin(theta_m))    

        for n in range(1, Nperim+1):

            nx = n % Nperim
            v1 = Vertex(xn[nx-1], yn[nx-1], dz , None, None)
            v2 = Vertex(xn[nx], yn[nx], dz , None, None)    
            v3 = Vertex(xm[nx], ym[nx], 0.0 , None, None)   
            #check orientation in second crown
            sub = create_subunit_from_vertices(v2, v1, v3, 'sub1')
            capsid.add_new_subunit(sub)

    def create_third_crown():    

        #third(shifted) crown    
        xn = []
        yn = []

        xm = []
        ym = []


        for n in range(Nperim):
            theta_n = n*theta + theta/2
            xn.append(R*np.cos(theta_n))
            yn.append(R*np.sin(theta_n))

            theta_m = (n)*theta
            xm.append(R*np.cos(theta_m))
            ym.append(R*np.sin(theta_m))    

        for n in range(1, Nperim+1):

            nx = n % Nperim
            v1 = Vertex(xn[nx-1], yn[nx-1], dz , None, None)
            v2 = Vertex(xn[nx], yn[nx], dz , None, None)    
            v3 = Vertex(xm[nx], ym[nx], 2*dz , None, None)   
            sub = create_subunit_from_vertices(v1, v2, v3, 'sub1')
            capsid.add_new_subunit(sub)    

    def create_fourth_crown():    

        #fourth(shifted and flipped) crown    
        xn = []
        yn = []

        xm = []
        ym = []


        for n in range(Nperim):
            theta_n = n*theta #- theta/2
            xn.append(R*np.cos(theta_n))
            yn.append(R*np.sin(theta_n))

            theta_m = (n-1)*theta + theta/2
            xm.append(R*np.cos(theta_m))
            ym.append(R*np.sin(theta_m))    

        for n in range(1, Nperim+1):

            nx = n % Nperim
            v1 = Vertex(xn[nx-1], yn[nx-1], 2*dz , None, None)
            v2 = Vertex(xn[nx], yn[nx], 2*dz , None, None)    
            v3 = Vertex(xm[nx], ym[nx], dz , None, None)   
            sub = create_subunit_from_vertices(v2, v1, v3, 'sub1')
            capsid.add_new_subunit(sub)     
   
    def merge_hubs():
        hubs = capsid.hubs
        #has to be a while loop; while there are close hubs, keep merging
        is_there_to_merge = True
        while is_there_to_merge:
            hubs = capsid.hubs
            is_there_to_merge = False
            for hub1 in hubs:
                for hub2 in hubs:
                    if hub1==hub2:
                        continue
                    if distance(hub1.vertices[0], hub2.vertices[0])<1e-4:
                        v2 = hub2.vertices[0]
                        hub2.merge_into(hub1)
                        is_there_to_merge = True
                        break
                if (is_there_to_merge):
                    break

    
    create_first_crown()
    create_second_crown()
    if (L<>1):
        create_third_crown()
        create_fourth_crown()
    
    new_subs = []
    if L%2<>0:
        rn = range(L/2)
    else:
        rn = range(L/2-1)
    for l in rn:
        subs = capsid.subunits

        for sub in subs:
            new_sub = copy.deepcopy(sub)
            for v in new_sub.vertices:
                v.z+=2*dz*(l+1)
            new_subs.append(new_sub)    

    if L%2<>0:
        new_subs = new_subs[:-2*Nperim]
    for new_sub in new_subs:
        capsid.add_new_subunit(new_sub)
        
    Ns = len(capsid.subunits)
    if L%2<>0:
        #print "bb"
        create_third_crown()
        for sub in capsid.subunits[Ns:]:
            for v in sub.vertices:
                l = rn[-1]
                v.z+=2*dz*(l+1)
        Ns = len(capsid.subunits)        
        create_fourth_crown()
        for sub in capsid.subunits[Ns:]:
            for v in sub.vertices:
                l = rn[-1]
                v.z+=-2*dz#0*2*dz*(l+1)         
                
    Ns = len(capsid.subunits)                
    if L%2==0:
        print "aa"
        create_first_crown()
        for sub in capsid.subunits[Ns:]:
            for v in sub.vertices:
                if L==2:
                    l=-1
                else:
                    l = rn[-1]
                v.z+=2*dz*(l+2)
        Ns = len(capsid.subunits)        
        create_fourth_crown()
        for sub in capsid.subunits[Ns:]:
            for v in sub.vertices:
                if L==2:
                    l=0
                else:
                    l = rn[-1]                
                v.z+=-2*dz#0*2*dz*(l+1)                   
            

    merge_hubs()    
    return capsid   


def create_trumpet_hanging_wedges(Nperim, L, parameters):
    import monte_carlo as mc
    RNG = np.random.RandomState()
    capsid = create_tubule_hanging_wedges(Nperim, L)
    pr.update_capsid_parameters(capsid, parameters)    
    E_elast = []
    #nsteps = 1000#0
    #ndrop = 100#0
    mc.update_neighbor_lists(capsid)
    nsteps = parameters['global_params']['nsteps_quench']
    ndrop = nsteps/10
    pr.d_max = 0.1
    for s in range(2*nsteps):  
        #pr.d_max=0.05#(2*1.0/50.0)**0.5
        if s%ndrop==0:
            pr.d_max/=1.2
        for i in range(len(capsid.hubs)):
            mc.attempt_move_athermal(capsid, RNG)
        E_elast.append(mc.elastic_energy(capsid.subunits))
    return capsid, E_elast  


def crack_edge(capsid):
    hub = capsid.get_surface_hubs()[0]
    vertices = hub.vertices
    for vertex in vertices:
        if sum([e.is_surface_edge() for e in vertex.parent_subunit.edges]) == 1:
            hub.split_from([vertex])        
            break    

#converts a capsid object to faces.dat, edges.dat, vertices.dat for the C++ code
def convert_capsid_to_cpp_dat(capsid, outfolder):
    vertices = [h.vertices[0] for h in capsid.hubs]
    faces = capsid.subunits
    of = open(outfolder+"/vertices.dat", "w")
    for ix, v in enumerate(vertices):
        print ix+1, v.x, v.y, v.z
        of.write(str(ix+1)+" "+str(v.x)+" "+ str(v.y)+" "+str(v.z)+"\n")
    of.close()
    
    of = open(outfolder+"/faces.dat", "w")
    for ix, f in enumerate(faces):
        e1, e2, e3 = f.edges
        v1 = e1.vertex_from
        v2 = e1.vertex_to
        v3 = e1.vertex_to.edge_out.vertex_to

        v1i = vertices.index(v1.parent_hub.vertices[0])
        v2i = vertices.index(v2.parent_hub.vertices[0])
        v3i = vertices.index(v3.parent_hub.vertices[0])    
        print ix+1, v1i+1, v2i+1, v3i+1, f.subunit_type.replace('sub','')
        of.write(str(ix+1)+" "+str(v1i+1)+" "+str(v2i+1)+" "+str(v3i+1)+" "+f.subunit_type.replace('sub','')+"\n")
    of.close()    
    
    of = open(outfolder+"/edges.dat", "w")
    for ix, f in enumerate(faces):
        for ixe, e in enumerate(f.edges):
            v1 = e.vertex_from
            v2 = e.vertex_to

            v1i = vertices.index(v1.parent_hub.vertices[0])
            v2i = vertices.index(v2.parent_hub.vertices[0])

            print ix*3+ixe+1, v1i+1, v2i+1, e1.edge_type.replace('type', '')
            of.write(str(ix*3+ixe+1)+" "+ str(v1i+1)+" "+ str(v2i+1)+" "+ e.edge_type.replace('type', '')+"\n")
    of.close()
