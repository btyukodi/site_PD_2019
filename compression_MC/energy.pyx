# cython: profile=False

#from Capsid import *
from geometry import *
from helpers import *
import numpy as np
import math

cimport numpy as np
import cython
from cpython cimport array

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

cdef extern from "math.h":
    double sqrt(double x)
    double fabs(double x) 
    double cos(double x)
    double acos(double x)

def edge_stretch_energy(edge):
    """Computes the stretching energy on an edge"""
    v1 = edge.vertex_from
    v2 = edge.vertex_to
    l = distance(v1, v2)
    return 0.5*edge.elastic_modulus*(l - edge.equilibrium_length)**2

def subunit_energy(subunit):
    """Computes the stretching energy of a subunit"""
    #return sum(map(lambda e: edge_stretch_energy(e), subunit.edges))
    return sum([edge_stretch_energy(e) for e in subunit.edges])

def subunit_energy_old(subunit):
    """Computes the stretching energy of a subunit; could be rewritten to use edge_stretch_energy"""
    E = 0.0
    for v in subunit.vertices:
        edge = v.edge_out
        vnext = edge.vertex_to
        l = distance(v, vnext)
        E += 0.5*edge.elastic_modulus*(l - edge.equilibrium_length)**2
    return E    

def bending_energy(edge):
    """Computes the bending energy along an edge. Note that this energy is shared between the two edges connecting the common hubs."""
    hub1 = edge.vertex_from.parent_hub
    hub2 = edge.vertex_to.parent_hub
    edges = get_edges_between_hubs(hub1, hub2)
    if len(edges)==1:
        return 0 #surface hub, bears no energy
    sub1 = edges[0].parent_subunit
    sub2 = edges[1].parent_subunit
    cdef double bending_modulus = edges[0].bending_modulus[ edges[1].edge_type]#0.5*(edges[0].bending_modulus + edges[1].bending_modulus) #
    cdef double equilibrium_angle = edges[0].equilibrium_angle[ edges[1].edge_type]#0.5*(edges[0].equilibrium_angle + edges[1].equilibrium_angle) #
    cdef double norm1[3] 
    norm1 = get_triangle_normal(sub1)
    cdef double norm2[3] 
    norm2 = get_triangle_normal(sub2)
    #update_subunit_normal(sub1)
    #update_subunit_normal(sub2)
    #norm1 = sub1.normal
    #norm2 = sub2.normal

    e0 = edges[0]
    v1 = e0.vertex_from
    v2 = e0.vertex_to

    cdef double edge_vector[3]
    edge_vector[:] = [v2.x-v1.x, v2.y-v1.y, v2.z-v1.z]
    cdef double n1xn2[3] 
    n1xn2 = [norm1[1]*norm2[2] - norm1[2]*norm2[1], norm1[2]*norm2[0] - norm1[0]*norm2[2], norm1[0]*norm2[1] - norm1[1]*norm2[0]]
    cdef bint convex = (n1xn2[0]*edge_vector[0] + n1xn2[1]*edge_vector[1] + n1xn2[2]*edge_vector[2])<0
    #convex = np.dot(cross3D(norm1, norm2), edge_vector)<0
    #theta = math.arcos(np.dot(norm1, norm2))
    #print norm1, norm2, (norm1[0]*norm2[0] + norm1[1]*norm2[1]+norm1[2]*norm2[2])
    dotprod = norm1[0]*norm2[0] + norm1[1]*norm2[1]+norm1[2]*norm2[2]
    if dotprod>1.0:
        print acos(dotprod)
        raise ValueError("dotprod > 1")
    if dotprod<-1.0:
        print acos(dotprod)
        raise ValueError("dotprod < -1")
    cdef double theta = acos(dotprod)
    if convex:
        theta = -theta
    return bending_modulus*(1.0-cos(theta - equilibrium_angle))


def elastic_energy(affected_subunits):
    """Computes the stretching+bending energies of a list of subunits. Function used to calculate pre- and after- move energies"""
    E = 0.0
    already_counted = []
    for subunit in affected_subunits:
        E+=subunit_energy(subunit)
        for edge in subunit.edges:
            hub1 = edge.vertex_from.parent_hub
            hub2 = edge.vertex_to.parent_hub
            if set([hub1, hub2]) in already_counted: #avoid double counting of edge
                continue
            already_counted.append(set([hub1, hub2]))
            E+=bending_energy(edge)
    return E

#includes elastic, binding and -mu*N
def full_energy(capsid):
    E = elastic_energy(capsid.subunits)
    K = 0.0
    mu_N = 0.0
    already_counted = []
    for hub1 in capsid.hubs:
        for hub2 in capsid.hubs:
            if set([hub1, hub2]) in already_counted: #avoid double counting of edge
                continue            
            edges = get_edges_between_hubs(hub1, hub2)
            if len(edges)==2:
                edge1, edge2 = edges
                K+=edge1.binding_affinity[edge2.edge_type]
    for subunit in capsid.subunits:
        mu_N+=subunit.activity
    return E+K-mu_N  

#DEBUG
def elastic_energy_DB(affected_subunits):
    #print affected_subunits[0].vertices[0]
    capsid = affected_subunits[0].vertices[0].parent_hub.parent_capsid  
    N = len(capsid.subunits)
    return 0.0#-0.05*N  

def edge_energy(edge):
    """0.5*Bending+stretching energy of an edge. Bending energy is shared between adjacent bonds, thus the 0.5 factor. Used in data processing."""
    return 0.5*bending_energy(edge) + edge_stretch_energy(edge)    

#checks overlap of a list of subunits with a capsid
def check_overlap(subunits, capsid):
    """Checks overlap of a list of subunits with a capsid. Call helpers.update_excluders_position() 
    if excluder positions are not up-to-date"""
    #check for overlaps; this sucks because it scales as number_of_subunits^2
    #update_excluders_position(subunits)
    for subunit1 in capsid.subunits:
        for subunit2 in subunits:
            if subunit1.do_overlap(subunit2):
                return True
    return False

def check_neighbor_overlap(subunits):
    return (True in [sub1.do_overlap(sub2) for sub1 in subunits for sub2 in sub1.neighbor_list])
    #for sub1 in subunits:
    #    if True in [sub1.do_overlap(sub2) for sub2 in sub1.neighbor_list]:
    #        return True
        #for sub2 in sub1.neighbor_list:
        #    if sub1.do_overlap(sub2):
        #        return True
    #return False

#eventually one could implement hard constraints on maximum 
#edge length and bending angle to make sure that vertices never leave
#the l_fusion^3 and l_add^3 volumes
def hard_constraints(affected_subunits):
    pass
