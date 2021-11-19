# cython: profile=False

import numpy as np
#from Capsid import *
from numpy import eye, dot, cross
from scipy.linalg import expm#, norm
#from math import fabs

cimport numpy as np
import cython
from cpython cimport array

#import math

cdef extern from "math.h":
    double sqrt(double x)
    double fabs(double x)    

def cross3D(U, V):
    return [U[1]*V[2] - U[2]*V[1], U[2]*V[0] - U[0]*V[2], U[0]*V[1] - U[1]*V[0]]
    #return np.array([U[1]*V[2] - U[2]*V[1], U[2]*V[0] - U[0]*V[2], U[0]*V[1] - U[1]*V[0]]).astype(float)

cdef array.array cross3D_c(array.array U, array.array  V):
    return array.array('d', [U[1]*V[2] - U[2]*V[1], U[2]*V[0] - U[0]*V[2], U[0]*V[1] - U[1]*V[0]])

def norm3D(U):
    return float(U[0]**2 + U[1]**2 + U[2]**2)**0.5

def distance(vertex1, vertex2):
        cdef double xA, yA, zA 
        xA, yA, zA = vertex1.x, vertex1.y, vertex1.z
        cdef double xB, yB, zB 
        xB, yB, zB = vertex2.x, vertex2.y, vertex2.z
        return sqrt((xA-xB)*(xA-xB) + (yA-yB)*(yA-yB) + (zA-zB)*(zA-zB))

def within_cube(vertex1, vertex2, double L):
        cdef xA, yA, zA 
        xA, yA, zA = vertex1.x, vertex1.y, vertex1.z
        cdef xB, yB, zB 
        xB, yB, zB= vertex2.x, vertex2.y, vertex2.z 
        return (fabs(xA-xB)<L) & (fabs(yA-yB)<L) & (fabs(zA-zB)<L) 
    
def M(axis, theta):
    return expm(cross(eye(3), axis/norm3D(axis)*theta))

#rotation as per https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
def rotate_vertex(vertex, axis, theta):
    M0 = M(axis, theta)
    vertex.x, vertex.y, vertex.z = dot(M0, [vertex.x, vertex.y, vertex.z]) 
def rotate_hub(hub, axis, theta):
    M0 = M(axis, theta)
    for vertex in hub.vertices:
        vertex.x, vertex.y, vertex.z = dot(M0, [vertex.x, vertex.y, vertex.z]) 
   
       
cpdef get_triangle_normal(subunit):
    p1, p2, p3 = subunit.vertices
    if (p1.edge_out.vertex_to <> p2):
        p1, p2, p3 = subunit.vertices[::-1]

    #U = [p2.x - p1.x, p2.y - p1.y, p2.z - p1.z]
    #V = [p3.x - p1.x, p3.y - p1.y, p3.z - p1.z]

    #normal = cross3D(U,V)

    cdef double normal[3]
    normal[:] = [(p2.y - p1.y)*(p3.z - p1.z) - (p2.z - p1.z)*(p3.y - p1.y), (p2.z - p1.z)*(p3.x - p1.x) - (p2.x - p1.x)*(p3.z - p1.z), (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x)]

    #U[1]*V[2] - U[2]*V[1] = (p2.y - p1.y)*(p3.z - p1.z) - (p2.z - p1.z)*(p3.y - p1.y)
    #U[2]*V[0] - U[0]*V[2] = (p2.z - p1.z)*(p3.x - p1.x) - (p2.x - p1.x)*(p3.z - p1.z)
    #U[0]*V[1] - U[1]*V[0] = (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x)
    cdef double N = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
    normal[0]/=N
    normal[1]/=N
    normal[2]/=N

    return normal
    #return normal/norm3D(normal)# faster than scipy's norm()

#cpdef normal_from_ordered_coordinates(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3):
#    cdef double normal[3]
#    normal[:] = [(y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1), (z2 - z1)*(x3 - x1) - (x2 - x1)*(z3 - z1), (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1)]
#
#        #U[1]*V[2] - U[2]*V[1] = (p2.y - p1.y)*(p3.z - p1.z) - (p2.z - p1.z)*(p3.y - p1.y)
#        #U[2]*V[0] - U[0]*V[2] = (p2.z - p1.z)*(p3.x - p1.x) - (p2.x - p1.x)*(p3.z - p1.z)
#        #U[0]*V[1] - U[1]*V[0] = (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x)
#    cdef double N = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
#    normal[0]/=N
#    normal[1]/=N
#    normal[2]/=N
#
#    return normal

def get_triangle_normal_old(subunit):
    p1, p2, p3 = subunit.vertices
    if (p1.edge_out.vertex_to <> p2):
        p1, p2, p3 = subunit.vertices[::-1]

    #U = [p2.x - p1.x, p2.y - p1.y, p2.z - p1.z]
    #V = [p3.x - p1.x, p3.y - p1.y, p3.z - p1.z]

    #normal = cross3D(U,V)

    normal = [(p2.y - p1.y)*(p3.z - p1.z) - (p2.z - p1.z)*(p3.y - p1.y),
        (p2.z - p1.z)*(p3.x - p1.x) - (p2.x - p1.x)*(p3.z - p1.z),
        (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x)]

    N = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
    normal[0]/=N
    normal[1]/=N
    normal[2]/=N

    return normal
    #return normal/norm3D(normal)# faster than scipy's norm()    


def update_subunit_normal(subunit):
    p1, p2, p3 = subunit.vertices
    U = [p2.x - p1.x, p2.y - p1.y, p2.z - p1.z] 
    V = [p3.x - p1.x, p3.y - p1.y, p3.z - p1.z]
    subunit.normal[[0,1,2]] = [U[1]*V[2] - U[2]*V[1], U[2]*V[0] - U[0]*V[2], U[0]*V[1] - U[1]*V[0]]
    subunit.normal/=norm3D(subunit.normal)