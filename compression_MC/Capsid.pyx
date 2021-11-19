# cython: profile=False

#import numpy as np
import pylab as pl
#import dummy0 as pl
#from itertools import product
import itertools
from numpy import cross, eye, dot, einsum, zeros
#from scipy.linalg import expm, norm
from geometry import *
import parameters as pr


class Edge(object):
    def __init__(self, vertex_from, vertex_to, equilibrium_length ,
                 edge_type='type1', edge_id=None, parent_subunit=None,
                elastic_modulus=1.0, bending_modulus=1.0, equilibrium_angle=0.0,
                binding_affinity=1.0):
        self.vertex_from = vertex_from
        self.vertex_to = vertex_to
        self.equilibrium_length = equilibrium_length
        self.edge_type = edge_type
        self.edge_id = edge_id
        self.parent_subunit = parent_subunit
        self.elastic_modulus = elastic_modulus
        self.bending_modulus = bending_modulus #this could come from a lookup between diff edge types
        self.equilibrium_angle = equilibrium_angle #this could come from a lookup between diff edge types
        self.binding_affinity = binding_affinity #this could come from a lookup between diff edge types
    def is_surface_edge(self):
        #need more careful check here based on neighbor parent hubs
        hub_from = self.vertex_from.parent_hub
        hub_to = self.vertex_to.parent_hub
        n_edges_between_hubs = 0
        #!for vertex in hub_from.vertices:
        #!    n_edges_between_hubs +=int(vertex.edge_out.vertex_to.parent_hub ==hub_to)
        n_edges_between_hubs +=sum([vertex.edge_out.vertex_to.parent_hub ==hub_to
        							for vertex in hub_from.vertices])

        #!for vertex in hub_to.vertices:
        #!    n_edges_between_hubs +=int(vertex.edge_out.vertex_to.parent_hub == hub_from)
        n_edges_between_hubs +=sum([vertex.edge_out.vertex_to.parent_hub ==hub_from 
        							for vertex in hub_to.vertices])

        #edges_between_hubs = filter(lambda vertex: vertex.edge_out.vertex_to.parent_hub ==hub_to, hub_from.vertices) +\
        #                    filter(lambda vertex: vertex.edge_out.vertex_to.parent_hub ==hub_from, hub_to.vertices)


        #if len(edges_between_hubs)==1:
        #    return True
        if n_edges_between_hubs==1:
            return True
        return False

        
    def plot(self,ax):
        #ax.quiver([self.vertex_from.x],
        #        [self.vertex_from.y],
        #        [self.vertex_from.z],
        #        [self.vertex_to.x-self.vertex_from.x],
        #        [self.vertex_to.y-self.vertex_from.y],
        #        [self.vertex_to.z-self.vertex_from.z], pivot='tail', length=distance(self.vertex_from, self.vertex_to))
        ax.plot([self.vertex_from.x, self.vertex_to.x],
               [self.vertex_from.y, self.vertex_to.y],
               [self.vertex_from.z, self.vertex_to.z],'-' )
        
class Hub(object):
    def __init__(self):
        self.vertices = []
        parent_capsid=None
    def add_vertex(self,vertex):
        self.vertices.append(vertex)
        vertex.parent_hub = self
    def remove_vertex(self,vertex):
        self.vertices.remove(vertex)
        vertex.parent_hub = None
    def merge_into(self, hub):
        for vertex in self.vertices:
            hub.add_vertex(vertex)
        self.vertices=[]
        self.parent_capsid.hubs.remove(self)
        self.parent_capsid=None
        #return current instance hub which is deleted
        return self
        #delete current instance if needed
    # if new_hub given, that will be used as the second hub    
    def split_from(self, vertices_to_split, new_hub=None):
        if len(self.vertices)==1:
            return
        if not new_hub:    
            hub = Hub()
        else:
            hub = new_hub            
        for vertex in vertices_to_split:
            hub.add_vertex(vertex)
        self.vertices = [vertex for vertex in self.vertices if vertex not in vertices_to_split]
        self.parent_capsid.hubs.append(hub)
        hub.parent_capsid = self.parent_capsid
        return hub
    def is_surface_hub(self):
        edges = [vertex.edge_out for vertex in self.vertices] +\
                [vertex.edge_in for vertex in self.vertices]
        surface_edges_in_hub = filter(lambda x: x.is_surface_edge(), edges)
        if len(surface_edges_in_hub)>0:
            return True
        return False
               
    def get_neighbor_hubs(self):
        #neighbor_hubs = []
        #for vertex in self.vertices:
        #    neighbor_hubs.append(vertex.edge_out.vertex_to.parent_hub)
        #    neighbor_hubs.append(vertex.edge_in.vertex_from.parent_hub)
        neighbor_hubs =  [vertex.edge_out.vertex_to.parent_hub for vertex in self.vertices ] +  [vertex.edge_in.vertex_from.parent_hub for vertex in self.vertices ]
        return list(set(neighbor_hubs))
    
    def move(self, dx, dy, dz):
        for vertex in self.vertices:
            vertex.x+=dx
            vertex.y+=dy
            vertex.z+=dz
        
    def plot(self, ax):
        x, y, z = self.vertices[0].x, self.vertices[0].y, self.vertices[0].z
        ax.plot([x], [y], [z], 'o', markersize=5, color='g')
        

class Vertex(object):
    #@cpdef public double x, y, z
    #@cdef public object edge_in, edge_out, parent_subunit, parent_hub
    def __init__(self, x, y, z, edge_in, edge_out, vertex_type=0, vertex_id=None, parent_subunit=None):
    #@def __init__(self, double x, double y, double z, edge_in, edge_out, parent_subunit=None):        
        self.x = x
        self.y = y
        self.z = z
        #$self.vertex_type = vertex_type
        #$self.vertex_id = vertex_id
        self.parent_subunit = parent_subunit
        #doubly linked list via edge_in, edge_out
        self.edge_in = edge_in
        self.edge_out = edge_out
        #hub the vertex belongs to
        self.parent_hub = None
    #@cpdef void move(self, double dx, double dy, double dz):
    def move(self, dx, dy, dz):
        self.x+=dx
        self.y+=dy
        self.z+=dz

    def plot(self, ax):
        ax.plot([self.x], [self.y], [self.z], 'o', color='b')

class Excluder(object):
    def __init__(self, R, get_position, parent_subunit=None, excluder_type=None, overlaps=None):
        self.R=R
        self.parent_subunit = parent_subunit
        self.get_pos = get_position #pointer to function which computes its position based on vertex positions
        self.excluder_type=excluder_type
        self.overlaps = overlaps
    def get_position(self):
        return self.get_pos(self.parent_subunit)
    def update_position(self):
        self.x, self.y, self.z = self.get_position()
    def plot(self,ax):
        x, y, z = self.x, self.y, self.z#self.get_position()
        ax.plot([x], [y], [z], 'o', markersize=2, color='r') #change this to plotting sphere eventually
        
class Subunit(object):
    def __init__(self, activity=1.0):
        #vertex list in the subunit; could use single vertex only, but easier to loop through
        self.vertices = []
        self.edges = []
        self.excluders = []
        self.activity = activity #exp(\beta \mu) / \lambda^9
        #number of distinguishable rotational configurations
        self.number_of_rotational_configurations = 1
        self.normal = zeros(3)
    def translate(self, dx, dy, dz):
        map(lambda v: v.move(dx, dy, dz), self.vertices)
    def do_overlap(self, subunit):
        #nothing overlaps with itself
        if subunit==self:
            return False
##        for excluder1 in self.excluders:
##            for excluder2 in subunit.excluders:
##                #if not (pr.Exclude_Lookup[excluder1.excluder_type, excluder2.excluder_type]):
##                if not excluder1.overlaps[excluder2.excluder_type]:
##                    continue
##                #x1, y1, z1 = excluder1.get_position()
##                #x2, y2, z2 = excluder2.get_position()
##                x1, y1, z1 = excluder1.x, excluder1.y, excluder1.z
##                x2, y2, z2 = excluder2.x, excluder2.y, excluder2.z
##                #dr = [x1-x2, y1-y2, z1-z2]
##                #dr2 = einsum('i,i->', dr, dr)
##                R12 = excluder1.R + excluder2.R
##                #if dr2 < R12**2:
##                if ((x1-x2)**2 + (y1-y2)**2+(z1-z2)**2)<R12*R12:
##                    return True

        return True in [(excluder1.x-excluder2.x)**2 + (excluder1.y-excluder2.y)**2+(excluder1.z-excluder2.z)**2 < (excluder1.R + excluder2.R)**2 
                        for excluder1 in self.excluders for excluder2 in subunit.excluders]            

##        return False
    def plot(self, ax):
        for edge in self.edges:
            edge.plot(ax)
        for vertex in self.vertices:
            vertex.plot(ax)
        for excluder in self.excluders:
            excluder.plot(ax)
        
    #rotate around point and axis, deform along m0, m1, m2
    #def rotate, deform, translate
    
    
class Capsid(object):
    def __init__(self):
        self.subunits = []
        self.hubs = []
        self.kT = 1.0
    def get_vertices(self):
        vertices_nested = [subunit.vertices for subunit in self.subunits]
        return list(itertools.chain(*vertices_nested))
    def get_edges(self):
        edges_nested = [subunit.edges for subunit in self.subunits]
        return list(itertools.chain(*edges_nested))    
    def get_surface_edges(self):
        return filter(lambda x: x.is_surface_edge(), self.get_edges())
    def get_surface_vertices(self):
        surface_edges = self.get_surface_edges()
        surface_vertices = list(set([edge.vertex_from for edge in surface_edges]+
                                [edge.vertex_to for edge in surface_edges]))
        return surface_vertices
    def get_surface_hubs(self):
        surface_vertices = self.get_surface_vertices()
        surface_hubs = list(set([vertex.parent_hub for vertex in surface_vertices]))
        return surface_hubs
    
    def add_new_subunit(self, subunit):
        self.subunits.append(subunit)
        for vertex in subunit.vertices:
            self.hubs.append(vertex.parent_hub)
            vertex.parent_hub.parent_capsid = self
        
    def remove_subunit(self, subunit):
        self.subunits.remove(subunit)
        for vertex in subunit.vertices:
            self.hubs.remove(vertex.parent_hub)
            vertex.parent_hub.parent_capsid = None       
                    
  
    def connect_subunit(self, capsid_h1, capsid_h2, subunit_h1, subunit_h2):
        if subunit_h1.vertices[0].parent_subunit not in self.subunits:
            raise ValueError('Add subunit to capsid before connecting')
        if (distance(capsid_h1.vertices[0], subunit_h1.vertices[0])+
            distance(capsid_h2.vertices[0], subunit_h2.vertices[0])>1e-3):
                raise ValueError('Vertices too far to connect')
        if capsid_h1.vertices[0].parent_subunit.do_overlap(subunit_h1.vertices[0].parent_subunit):
            raise ValueError('Cannot connect: excluders overlap')
        subunit_h1.merge_into(capsid_h1)
        subunit_h2.merge_into(capsid_h2)

    def disconnect_subunit(self, subunit):
        for vertex in subunit.vertices:
            vertex.parent_hub.split_from([vertex])

                    
    def plot(self, ax):
        for subunit in self.subunits:
            subunit.plot(ax)
        for hub in self.hubs:
            hub.plot(ax)

        
