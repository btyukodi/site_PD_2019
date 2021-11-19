# cython: profile=False
#from Capsid import *
from helpers import *
from energy import *
#from geometry import *
import parameters as pr
import math

#DB
import shelve
import gdbm


def attempt_move(capsid, RNG, wall=None):
    #choose with equal probability from hubs to move
    candidates = capsid.hubs
    which = RNG.randint(len(candidates))
    candidate = candidates[which]

    #DB!!! -------- don't move vertices of sub2 type subunits
    #if 'sub2' in [v.parent_subunit.subunit_type for v in candidate.vertices]:
    #    return False
    #END DB -------

    vertex = candidate.vertices[0]
    affected_subunits = [v.parent_subunit for v in candidate.vertices]
    
    energy_before = elastic_energy(affected_subunits)
    if wall is not None:
        energy_before+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))


    #for now move in a box volume dV = d_max^3
    dx, dy, dz = (0.5-RNG.rand(3))*2.0*pr.d_max

    candidate.move(dx, dy, dz)
    #check for overlaps; this sucks because it scales as number_of_subunits^2

    update_excluders_position(affected_subunits)
    #if check_overlap(affected_subunits, capsid):
    if check_neighbor_overlap(affected_subunits):
        candidate.move(-dx, -dy, -dz)
        update_excluders_position(affected_subunits)
        #pr.rej+=1
        return False

    #!!! should check if vertex left v_add volume
            
    energy_after = elastic_energy(affected_subunits)
    if wall is not None:
        energy_after+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))


    dE = energy_after-energy_before

    #DB ----
    #if dE < 1000.0:
    #    dE=0.0
    #--DB      

    if (math.exp(-dE/capsid.kT) > RNG.rand()):
        return True
    else:
        candidate.move(-dx, -dy, -dz)
        update_excluders_position(affected_subunits)
        return False


#rename this to attempt_free_insertion()  and write an attempt_insertion() wrapper in which you decide what kind of insertion to make
#do the same for removal and rewrite helpers.get_removable subunits to return all subunits at the edge    
def attempt_insertion(capsid, surface_edges, RNG, subunit_type='sub1', harmonic_umbrella_window=None):

    if len(surface_edges)==0:
        return False
    #pick an edge
    which = RNG.randint(len(surface_edges))
    edge = surface_edges[which]    
    v_from = edge.vertex_from
    v_to = edge.vertex_to
    affected_subunits = [edge.parent_subunit]
    energy_before = elastic_energy(affected_subunits)

    if harmonic_umbrella_window is not None:
        energy_before+=harmonic_umbrella_window.bias_potential(len(capsid.subunits))

    v1 = Vertex(v_to.x,v_to.y,v_to.z, None, None)
    v2 = Vertex(v_from.x,v_from.y,v_from.z, None, None) 
    #compute position of new vertex
    r0 = 0.5*np.array([v1.x+v2.x, v1.y+v2.y, v1.z+v2.z])
    n = get_triangle_normal(edge.parent_subunit)
    a = np.array([v1.x, v1.y, v1.z])-r0
    nxa = np.array(cross3D(n,a)).astype(float)
    v3_pos = -nxa/norm3D(nxa)*edge.equilibrium_length*np.sqrt(3.0)/2.0
    #print "xxx", v1.x, v1.y, v1.z
    
    dx, dy, dz = (0.5-RNG.rand(3))*pr.l_add
    #print "###", dx, dy, dz
    v3 = Vertex(v3_pos[0],v3_pos[1],v3_pos[2], None, None)
    #!!! need to rotate if equilibrium edge angle is not 0
    #theta = edge.equilibrium_angle #pr.equilibrium_angle for now #might need to symmetrize:0.8=5* + sub.edges[0].equilibrium_angle)

    #sub = create_curved_triangle_subunit_from_vertices(v1, v2, v3)
    sub = create_subunit_from_vertices(v1, v2, v3, subunit_type)
     
    #shift edges/rotate subunit; will need to be rewritten properly for other shapes

    nrotations = RNG.randint(sub.number_of_rotational_configurations)
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
 

    theta = edge1.equilibrium_angle[edge2.edge_type]#pr.equilibrium_angle
    axis = a/norm3D(a)
    rotate_vertex(v3,axis,theta)
    v3.x+=r0[0]+dx
    v3.y+=r0[1]+dy
    v3.z+=r0[2]+dz

    update_excluders_position([sub])
    if check_overlap([sub], capsid):
        return False
    capsid.add_new_subunit(sub)
    capsid.connect_subunit(v_to.parent_hub, v_from.parent_hub, v1.parent_hub, v2.parent_hub)
    affected_subunits.append(sub)
    #could be made faster by setting dE = subunit_energy(sub) + bending_energy(edge)
    energy_after = elastic_energy(affected_subunits)
    if harmonic_umbrella_window is not None:
        energy_after+=harmonic_umbrella_window.bias_potential(len(capsid.subunits))



    dE = energy_after - energy_before
    #DB ----
    #if dE < 1000.0:
    #    dE=0.0
    #--DB    

    #@zs = sub.activity#pr.activity
    #@K = edge1.binding_affinity[edge2.edge_type]#pr.affinity
    #@p = pr.l_add**3 *sub.number_of_rotational_configurations*math.exp((-dE+zs+K)/capsid.kT)
    mu = sub.mu
    e_b = edge1.e_b[edge2.edge_type]
    #this should probably be va**2
    p =  (pr.va * pr.l_add**3 / pr.v0**3) *sub.number_of_rotational_configurations*math.exp(-(dE-mu+e_b)/capsid.kT)
    #p =  (pr.va * pr.l_add**3 / pr.v0**3) *math.exp(-(dE-mu+e_b)/capsid.kT)

    #p=math.exp((zs)/capsid.kT)
    #print "p=",p
    if (RNG.rand() < p):
        #print p, which, v3.x, v3.y, v3.z
        update_neighbor_lists(capsid)
        return True
    capsid.disconnect_subunit(sub)
    capsid.remove_subunit(sub)
    return False
    

def attempt_wedge_insertion(capsid, wedges, RNG, subunit_type='sub1', harmonic_umbrella_window=None):
    #wedges are type1 fusion pairs
    #wedges = get_open_wedge_pairs(capsid)
    which = RNG.randint(len(wedges))
    if len(wedges)==0:
        return False    

    #@hub1, hub2 = wedges[which]
    #@hub3 = get_common_neighbor_hubs(hub1, hub2)
    hub1, hub2, hub3 = wedges[which]
    #@if len(hub3)<>1:
    #@    raise ValueError("Wedge insertion screwed up")
    #@hub3 = hub3[0]
    
    edge1 = get_edges_between_hubs(hub1, hub3)[0]
    edge2 = get_edges_between_hubs(hub2, hub3)[0]

    affected_subunits = [edge1.parent_subunit, edge2.parent_subunit]
    energy_before = 0.0#elastic_energy(affected_subunits)
    if harmonic_umbrella_window is not None:
        energy_before+=harmonic_umbrella_window.bias_potential(len(capsid.subunits))    


    if edge1.vertex_to in hub3.vertices:

        vA = edge1.vertex_to
        vB = edge1.vertex_from
        vC = edge2.vertex_to


    if edge1.vertex_from in hub3.vertices:

        vA = edge1.vertex_from
        vB = edge1.vertex_to
        vC = edge2.vertex_from

    v1 = Vertex(vA.x,vA.y,vA.z, None, None)
    v2 = Vertex(vB.x,vB.y,vB.z, None, None) 
    v3 = Vertex(vC.x,vC.y,vC.z, None, None)       


    if edge1.vertex_to in hub3.vertices:
        sub = create_subunit_from_vertices(v1, v2, v3, subunit_type)

    if edge1.vertex_from in hub3.vertices:
        sub = create_subunit_from_vertices(v1, v3, v2, subunit_type)        
    #shift edges/rotate subunit; will need to be rewritten properly for other shapes
    nrotations = RNG.randint(sub.number_of_rotational_configurations)
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

    update_excluders_position([sub])
    if check_overlap([sub], capsid):
        return False
    capsid.add_new_subunit(sub)
    #

    #print "CONN1..."
    capsid.connect_subunit(vA.parent_hub, vB.parent_hub, v1.parent_hub, v2.parent_hub)
    v3.parent_hub.merge_into(vC.parent_hub)
    #print "CONN2..."
    #capsid.connect_subunit(v1.parent_hub, v3.parent_hub, vA.parent_hub, vC.parent_hub)  
    #print "CONN3"  
    affected_subunits.append(sub)
    energy_after = subunit_energy(sub) + bending_energy(edge1) + bending_energy(edge2)#elastic_energy(affected_subunits)
    if harmonic_umbrella_window is not None:
        energy_after+=harmonic_umbrella_window.bias_potential(len(capsid.subunits))    
    #capsid.connect_subunit(v_to.parent_hub, v_from.parent_hub, v1.parent_hub, v2.parent_hub)
    dE = energy_after - energy_before
    #DB ----
    #if dE < 1000.0:
    #    dE=0.0
    #--DB

    #!zs = sub.activity #set sub.activity somewhere!!!
    #!K = 0.5*(edge.binding_affinity + sub.edges[0].binding_affinity)
    #!p=2.0*sub.activity*K*K*pr.l_add**3 *np.exp(-dE)
    #@zs = sub.activity#pr.activity

    edges1 = get_edges_between_hubs(vA.parent_hub, vB.parent_hub)
    edges2 = get_edges_between_hubs(vA.parent_hub, vC.parent_hub)    

    #@K1 = edges1[0].binding_affinity[edges1[1].edge_type]
    #@K2 = edges2[0].binding_affinity[edges2[1].edge_type]
    
    #@K = K1+K2#pr.affinity
    #K = edge.binding_affinity[new_edge_type]#pr.affinity
    #p=zs*K*pr.l_add**3 *np.exp(-dE)
    #now zs is ln(zs) and K is ln(K)
    #@p = sub.number_of_rotational_configurations*math.exp((-dE+zs+K)/capsid.kT)


    mu = sub.mu    
    e_b1 = edges1[0].e_b[edges1[1].edge_type]
    e_b2 = edges2[0].e_b[edges2[1].edge_type]


    e_b = e_b1 + e_b2
    p =  (pr.va**2 / pr.v0**3) *sub.number_of_rotational_configurations*math.exp(-(dE-mu+e_b)/capsid.kT)
    #p =  (pr.va**2 / pr.v0**3) *math.exp(-(dE-mu+e_b)/capsid.kT)

    #print "p=",p
    if (RNG.rand() < p):
        #print p, which, v3.x, v3.y, v3.z
        update_neighbor_lists(capsid)
        return True
    capsid.disconnect_subunit(sub)
    capsid.remove_subunit(sub)
    return False


def attempt_wedge_removal(capsid, removable_wedge_subunits, RNG, harmonic_umbrella_window=None):
    if len(removable_wedge_subunits)==0:
        return False
    which = RNG.randint(len(removable_wedge_subunits))
    subunit = removable_wedge_subunits[which]
    subunit_surface_edge = filter(lambda x: x.is_surface_edge(), subunit.edges)
    #non-surface edges of subunit
    edge_subunit = [edge for edge in subunit.edges 
                    if edge not in subunit_surface_edge]
    hub1 = edge_subunit[0].vertex_from.parent_hub
    hub2 = edge_subunit[0].vertex_to.parent_hub

    hub3 = edge_subunit[1].vertex_from.parent_hub
    hub4 = edge_subunit[1].vertex_to.parent_hub


    edges1 = get_edges_between_hubs(hub1, hub2)
    edges2 = get_edges_between_hubs(hub3, hub4)

    neighbor1 =[s.parent_subunit for s in edges1]
    neighbor1.remove(subunit)
    neighbor1 = neighbor1[0]
    neighbor2 =[s.parent_subunit for s in edges2]
    neighbor2.remove(subunit)
    neighbor2 = neighbor2[0]

    affected_subunits = [subunit, neighbor1, neighbor2]    
    #@energy_before = elastic_energy(affected_subunits)
    #@energy_after = elastic_energy([neighbor1, neighbor2])#subunit_energy(neighbor1) + subunit_energy(neighbor2)    

    #@dE = energy_after - energy_before

    dE = -(subunit_energy(subunit) + bending_energy(edge_subunit[0]) + bending_energy(edge_subunit[1]))
    if harmonic_umbrella_window is not None:
        dE+=harmonic_umbrella_window.bias_potential(len(capsid.subunits)-1) - harmonic_umbrella_window.bias_potential(len(capsid.subunits))

    #DB ----
    #if dE < 1000.0:
    #    dE=0.0
    #--DB    

    #@zs = subunit.activity#pr.activity
    


    #@K1 = edges1[0].binding_affinity[edges1[1].edge_type]
    #@K2 = edges2[0].binding_affinity[edges2[1].edge_type]

    #@K = K1 + K2
    #p=1.0/(zs*K*pr.l_add**3) *np.exp(-dE)
    #@p=1.0/(subunit.number_of_rotational_configurations) *math.exp((-dE-zs-K)/capsid.kT)

    mu = subunit.mu    
    e_b1 = edges1[0].e_b[edges1[1].edge_type]
    e_b2 = edges2[0].e_b[edges2[1].edge_type]


    e_b = e_b1 + e_b2
    p =  1.0/((pr.va**2 / pr.v0**3) *subunit.number_of_rotational_configurations)*math.exp(-(dE+mu-e_b)/capsid.kT)
    #p =  1.0/((pr.va**2 / pr.v0**3) )*math.exp(-(dE+mu-e_b)/capsid.kT)

    #print "p=",p
    if (p>RNG.rand()):
        capsid.disconnect_subunit(subunit)
        capsid.remove_subunit(subunit)
        update_neighbor_lists(capsid)
        return True
    return False    


def attempt_removal(capsid, removable_subunits, RNG, harmonic_umbrella_window=None):
    #removable_subunits = get_removable_subunits(capsid)
    if len(removable_subunits)==0:
        #print "nothing to remove"
        return False
    which = RNG.randint(len(removable_subunits))
    subunit = removable_subunits[which]
    subunit_surface_edges = filter(lambda x: x.is_surface_edge(), subunit.edges)

    #non-surface edge of subunit
    edge_subunit = [edge for edge in subunit.edges 
                    if edge not in subunit_surface_edges][0]
    hub1 = edge_subunit.vertex_from.parent_hub
    hub2 = edge_subunit.vertex_to.parent_hub
    edge_neighbor =get_edges_between_hubs(hub1, hub2)

    edge_neighbor.remove(edge_subunit)
    edge_neighbor = edge_neighbor[0]

    neighbor_subunit = edge_neighbor.parent_subunit
    affected_subunits = [subunit, neighbor_subunit]
    #@neighbor_stretch_energy = subunit_energy(neighbor_subunit)
    #@energy_before = subunit_energy(subunit) + neighbor_stretch_energy + bending_energy(edge_subunit)#elastic_energy(affected_subunits)
    #print "BEND TEST: ", bending_energy(edge_subunit), bending_energy(edge_neighbor)
    #bending_energy(edge_subunit) == bending_energy(edge_neighbor)
    #@energy_after = neighbor_stretch_energy#elastic_energy([neighbor_subunit])#!!!!!!DEBUGsubunit_energy(neighbor_subunit)

    #@dE = energy_after - energy_before


    dE = -(subunit_energy(subunit) + bending_energy(edge_subunit))
    if harmonic_umbrella_window is not None:
        dE+=harmonic_umbrella_window.bias_potential(len(capsid.subunits)-1) - harmonic_umbrella_window.bias_potential(len(capsid.subunits))    
    #DB ----
    #if dE < 1000.0:
    #    dE=0.0
    #--DB

    #print "---- DE="+str(dE)
    #!zs = subunit.activity
    #!K = 0.5*(edge_subunit.binding_affinity + edge_neighbor.binding_affinity)
    #!p=1.0/(2.0*zs*K*K*pr.l_add**3) *np.exp(-dE)
    #@zs = subunit.activity#pr.activity
    #K = edge_subunit.binding_affinity[edge_neighbor.edge_type]#pr.affinity
    #@K = edge_subunit.binding_affinity[edge_neighbor.edge_type]#pr.affinity
    #p=1.0/(zs*K*pr.l_add**3) *np.exp(-dE)
    #@p=1.0/(pr.l_add**3*subunit.number_of_rotational_configurations) *math.exp((-dE-zs-K)/capsid.kT)

    mu = subunit.mu
    e_b = edge_subunit.e_b[edge_neighbor.edge_type]
    #this should probably be va**2
    p =  1.0/((pr.va * pr.l_add**3 / pr.v0**3) *subunit.number_of_rotational_configurations)*math.exp(-(dE+mu-e_b)/capsid.kT) 
    #p =  1.0/((pr.va * pr.l_add**3 / pr.v0**3) )*math.exp(-(dE+mu-e_b)/capsid.kT)    
    #p=math.exp((-zs)/capsid.kT)
    #print "p=",p
    if (p>RNG.rand()):
        capsid.disconnect_subunit(subunit)
        capsid.remove_subunit(subunit)
        update_neighbor_lists(capsid)
        return True
    return False

def attempt_type1_fusion(capsid, RNG, fusion_pairs, wall=None):
    #fusion_pairs = get_type1_fusion_pairs(capsid)
    if len(fusion_pairs)==0:
        return False
    which = RNG.randint(len(fusion_pairs))


    hub1, hub2, hub3 = fusion_pairs[which]
    #@hub3 = get_common_neighbor_hubs(hub1, hub2)[0]
    edge1 = get_edges_between_hubs(hub1, hub3)[0]
    edge2 = get_edges_between_hubs(hub2, hub3)[0]

    affected_subunits = list(set([vertex.parent_subunit for vertex in hub1.vertices]+
                                 [vertex.parent_subunit for vertex in hub2.vertices]))
    #E1 = elastic_energy(affected_subunits)
    #affected_subunits = 
    energy_before = elastic_energy(affected_subunits)
    if wall is not None:
        energy_before+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))

    
    dx = 0.5*(hub1.vertices[0].x - hub2.vertices[0].x)
    dy = 0.5*(hub1.vertices[0].y - hub2.vertices[0].y)
    dz = 0.5*(hub1.vertices[0].z - hub2.vertices[0].z)
    hub1.move(-dx, -dy, -dz)
    hub2.move(dx, dy, dz)

    update_excluders_position(affected_subunits)
    if check_overlap(affected_subunits, capsid):
        hub1.move(dx, dy, dz)
        hub2.move(-dx, -dy, -dz)
        update_excluders_position(affected_subunits)
        return False
    #print "merge start", number_of_edges_between_hubs(hub1, hub2)
    hub1_initial_vertices = hub1.vertices #potential source of bug: deep vs shallow copy
    hub1.merge_into(hub2)
    energy_after = elastic_energy(affected_subunits)
    if wall is not None:
        energy_after+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))        
    #E2 = elastic_energy(affected_subunits)

    dE = energy_after - energy_before

    #print "QQQQ ", E2-E1, dE, dE-(E2-E1)
    #if abs(dE-(E2-E1))>1e-8:
    #    print "RRRRR ", E2-E1, dE, dE-(E2-E1)

    #--DB 
    #if dE < 1000.0:
    #    dE=0.0
    #--DB
  


    #will need change depending on edge/fusion type
    #@K = edge1.binding_affinity[edge2.edge_type]#pr.affinity
    #no: K = pr.va*edge1.e_b[edge2.edge_type]
    #verify acceptance probabilities!

    #p=(K/pr.l_fusion**3 *np.exp(-dE))
    e_b = edge1.e_b[edge2.edge_type]
    #print "e_b fuse: ",e_b
    #print "dE fuse: ", dE
    #@p=(1.0/pr.l_fusion**3 *math.exp((-dE+K)/capsid.kT))
    p=1.0/((pr.l_fusion)**3) *pr.va * math.exp(-(dE+e_b)/capsid.kT)
    if (p>RNG.rand()):
        #print p
        return True
    hub1 = hub2.split_from(hub1_initial_vertices, hub1)
    hub1.move(dx, dy, dz)
    hub2.move(-dx, -dy, -dz)
    update_excluders_position(affected_subunits)
    return False

def attempt_type2_fusion(capsid, RNG, fusion_pairs, wall=None):
    #fusion_pairs = get_type2_fusion_pairs(capsid)
    if len(fusion_pairs)==0:
        return False
    which = RNG.randint(len(fusion_pairs))

    #will need to check if fusion is type1 or type2 for the energy increase
    hub1, hub2 = fusion_pairs[which]
    hub3, hub4 = get_common_neighbor_hubs(hub1, hub2)

    affected_subunits = list(set([vertex.parent_subunit for vertex in hub1.vertices]+
                                 [vertex.parent_subunit for vertex in hub2.vertices]))

    energy_before = elastic_energy(affected_subunits)
    if wall is not None:
        energy_before+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))    
    
    dx = 0.5*(hub1.vertices[0].x - hub2.vertices[0].x)
    dy = 0.5*(hub1.vertices[0].y - hub2.vertices[0].y)
    dz = 0.5*(hub1.vertices[0].z - hub2.vertices[0].z)
    hub1.move(-dx, -dy, -dz)
    hub2.move(dx, dy, dz)

    update_excluders_position(affected_subunits)
    if check_overlap(affected_subunits, capsid):
        hub1.move(dx, dy, dz)
        hub2.move(-dx, -dy, -dz)
        update_excluders_position(affected_subunits)
        return False
    #print "merge start", number_of_edges_between_hubs(hub1, hub2)
    hub1_initial_vertices = hub1.vertices #potential source of bug: deep vs shallow copy
    hub1.merge_into(hub2)
    energy_after = elastic_energy(affected_subunits)
    if wall is not None:
        energy_after+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))        

    edges23 = get_edges_between_hubs(hub2, hub3)
    edges24 = get_edges_between_hubs(hub2, hub4) 
    #@K1 = edges23[0].binding_affinity[edges23[1].edge_type]   
    #@K2 = edges24[0].binding_affinity[edges24[1].edge_type]   

    e_b1 = edges23[0].e_b[edges23[1].edge_type]   
    e_b2 = edges24[0].e_b[edges24[1].edge_type] 
    e_b = e_b1 + e_b2

    dE = energy_after - energy_before
    #will need change depending on edge/fusion type
    #K = pr.affinity
    #verify acceptance probabilities!

    #p=(K*K/pr.l_fusion**3 *np.exp(-dE))   
    #p = (1.0/pr.l_fusion**3 *np.exp(-dE+2*K)) 
    p = 1.0/pr.l_fusion**3 *pr.va**2 * math.exp(-(dE+e_b)/capsid.kT)
    if (p>RNG.rand()):
        return True
    hub1 = hub2.split_from(hub1_initial_vertices, hub1)
    hub1.move(dx, dy, dz)
    hub2.move(-dx, -dy, -dz)
    update_excluders_position(affected_subunits)
    return False
    
   
#wedge fission
def attempt_type1_fission(capsid, RNG, fission_pairs, wall=None):
    #fission_pairs = get_type1_fission_pairs(capsid)
    if len(fission_pairs)==0:
        return False
    which = RNG.randint(len(fission_pairs))
    #hub1 should be the surface hub to break
    hub1, hub2 = fission_pairs[which]
    edges = get_edges_between_hubs(hub1, hub2)
    #affected_subunits = list(set([vertex.parent_subunit for vertex in hub1.vertices]))
    affected_subunits = list(set([vertex.parent_subunit for vertex in hub1.vertices]+
                                 [vertex.parent_subunit for vertex in hub2.vertices]))
    #E1 = elastic_energy(affected_subunits)
    #affected_subunits = capsid.subunits
    energy_before = elastic_energy(affected_subunits)
    if wall is not None:
        energy_before+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))


    vertices_left = []

    next_hub = hub2
    for i in xrange(1000): #should be while loop; safety to avoid infinite loops
        if i>900:
            raise ValueError('Type 1 fission screwed up')
        edges12 = get_edges_between_hubs(hub1, next_hub)
        #when hit the nearest surface edge
        if len(edges12)==1:
            break
        subunit_left = edges12[0].parent_subunit
        subunit_right = edges12[1].parent_subunit        
        #find subunit_left vertex that's in hub1
        vertices_left.append(list(set(subunit_left.vertices) & set(hub1.vertices))[0])        
        #there should be two edges from/to this vertex
        edge1 = vertices_left[-1].edge_in
        edge2 = vertices_left[-1].edge_out
        possible_hubs = [edge1.vertex_from.parent_hub, edge2.vertex_to.parent_hub]
        possible_hubs.remove(next_hub)
        next_hub = possible_hubs[0]

    new_hub = hub1.split_from(vertices_left)
    dx, dy, dz = (0.5-RNG.rand(3))*pr.l_fusion#*0.5
    #dx, dy = (0.5-RNG.rand(2))*pr.l_fusion
    #dz = 
    hub1.move(dx, dy, dz)
    new_hub.move(-dx, -dy, -dz)
    #check overlap and energetics
    update_excluders_position(affected_subunits)
    if check_overlap(affected_subunits, capsid):
        hub1.move(-dx, -dy, -dz)
        new_hub.move(dx, dy, dz)
        #hub1.merge_into(new_hub)
        new_hub.merge_into(hub1)
        update_excluders_position(affected_subunits)
        return False

    #E2 = elastic_energy(affected_subunits)
    energy_after = elastic_energy(affected_subunits)
    if wall is not None:
        energy_after+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))    
    dE = energy_after - energy_before

    #print "xQQQQ ", E2-E1, dE, dE-(E2-E1)
    #if abs(dE-(E2-E1))>1e-8:
    #    print "xRRRRR ", E2-E1, dE, dE-(E2-E1)    
    #print "dE fiss:", dE, energy_after, energy_before
    #DB ----
    #if dE < 1000.0:
    #    dE=0.0
    #if dE>100:
        #dE=-10.0
    #    db = gdbm.open('/home/btyukodi/assembly_MC/DEBUG/Capsid','n')
    #    shelf = shelve.Shelf(db)
    #    shelf['0'] = capsid
    #    shelf.close()
    #-----

    #@K = edges[0].binding_affinity[edges[1].edge_type]#pr.affinity
    e_b = edges[0].e_b[edges[1].edge_type]
    #print "e_b fiss: ",e_b
    #p = 1.0/(K/pr.l_fusion**3) * np.exp(-dE)
    #@p = 1.0/(1.0/pr.l_fusion**3) * math.exp((-dE-K)/capsid.kT)
    p = 1.0/(1.0/((pr.l_fusion)**3) * pr.va) * math.exp(-(dE-e_b)/capsid.kT)
    if (RNG.rand()<p):
        print "type1 fission p=",p
        print "E_after, E_before", energy_after, energy_before
        return True #accept
    else:
        hub1.move(-dx, -dy, -dz)
        new_hub.move(dx, dy, dz)
        #hub1.merge_into(new_hub)
        new_hub.merge_into(hub1)
        update_excluders_position(affected_subunits)
        return False
 

#crack fission    
def attempt_type2_fission(capsid, RNG, fission_triplets, wall=None):
    #fission_triplets = get_type2_fission_triplets(capsid)
    if len(fission_triplets)==0:
        return False
    which = RNG.randint(len(fission_triplets))
    hub1, hub2, hub3 = fission_triplets[which]
    edges12 = get_edges_between_hubs(hub1, hub2)
    edges23 = get_edges_between_hubs(hub2, hub3)
    affected_subunits = list(set([vertex.parent_subunit for vertex in hub1.vertices]+
                                 [vertex.parent_subunit for vertex in hub2.vertices]+
                                 [vertex.parent_subunit for vertex in hub3.vertices]))
    energy_before = elastic_energy(affected_subunits)
    if wall is not None:
        energy_before+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))    
    vertices_left = []

    next_hub = hub1
    for i in xrange(1000): #should be while loop; safety to avoid infinite loops
        if i>900:
            raise ValueError('Type 2 fission screwed up')
        edges12 = get_edges_between_hubs(hub2, next_hub)

        subunit_left = edges12[0].parent_subunit
        subunit_right = edges12[1].parent_subunit        
        #find subunit_left vertex that's in hub1
        vertices_left.append(list(set(subunit_left.vertices) & set(hub2.vertices))[0])        
        #there should be two edges from/to this vertex
        edge1 = vertices_left[-1].edge_in
        edge2 = vertices_left[-1].edge_out
        possible_hubs = [edge1.vertex_from.parent_hub, edge2.vertex_to.parent_hub]
        possible_hubs.remove(next_hub)
        #when hit hub3
        next_hub = possible_hubs[0]
        if next_hub==hub3:
            break


    new_hub = hub2.split_from(vertices_left)
    dx, dy, dz = (0.5-RNG.rand(3))*pr.l_fusion #*0.5
    hub2.move(dx, dy, dz)
    new_hub.move(-dx, -dy, -dz)
    #check overlap and energetics
    update_excluders_position(affected_subunits)
    if check_overlap(affected_subunits, capsid):
        hub2.move(-dx, -dy, -dz)
        new_hub.move(dx, dy, dz)
        #hub2.merge_into(new_hub)
        new_hub.merge_into(hub2)
        update_excluders_position(affected_subunits)
        return False

    energy_after = elastic_energy(affected_subunits)
    if wall is not None:
        energy_after+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))      
    dE = energy_after - energy_before
    #K = pr.affinity
    #@K1 = edges12[0].binding_affinity[edges12[1].edge_type]
    #@K2 = edges23[0].binding_affinity[edges23[1].edge_type] 
    e_b1 = edges12[0].e_b[edges12[1].edge_type]
    e_b2 = edges23[0].e_b[edges23[1].edge_type] 
    e_b = e_b1 + e_b2

    #type 2 probabilities need to be computed
    #p = 1.0/(K*K/pr.l_fusion**3) * np.exp(-dE)
    #@p = 1.0/(1.0/pr.l_fusion**3) * math.exp((-dE-K1-K2)/capsid.kT)
    p = 1.0/(1.0/pr.l_fusion**3 * pr.va**2) * math.exp(-(dE-e_b)/capsid.kT)
    #print '^^^', p, dE, e_b

    if (RNG.rand()<p):
        return True #accept
    else:
        hub2.move(-dx, -dy, -dz)
        new_hub.move(dx, dy, dz)
        #hub2.merge_into(new_hub)
        new_hub.merge_into(hub2)
        update_excluders_position(affected_subunits)
        return False
   

def attempt_edge_fusion(capsid, RNG, edge_fusion_hubs, wall=None):
    #edge_fusion_hubs = get_edge_fusion_hubs(capsid)
    if len(edge_fusion_hubs)==0:
        return False
    which = RNG.randint(len(edge_fusion_hubs))
    hub1A, hub1B, hub2A, hub2B = edge_fusion_hubs[which]
    edgeA = get_edges_between_hubs(hub1A, hub2A)[0] #should be of length 1
    edgeB = get_edges_between_hubs(hub1B, hub2B)[0] #should be of length 1
    affected_subunits = list(set([vertex.parent_subunit for vertex in hub1A.vertices]+
                                 [vertex.parent_subunit for vertex in hub2A.vertices]+
                                 [vertex.parent_subunit for vertex in hub1B.vertices]+
                                 [vertex.parent_subunit for vertex in hub2B.vertices]))
    energy_before = elastic_energy(affected_subunits)
    if wall is not None:
        energy_before+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))        
    #fuse 1A and 2A
    dx1 = 0.5*(hub1A.vertices[0].x - hub1B.vertices[0].x)
    dy1 = 0.5*(hub1A.vertices[0].y - hub1B.vertices[0].y)
    dz1 = 0.5*(hub1A.vertices[0].z - hub1B.vertices[0].z)
    hub1A.move(-dx1, -dy1, -dz1)
    hub1B.move(dx1, dy1, dz1) 

    #fuse 1B and 2B
    dx2 = 0.5*(hub2A.vertices[0].x - hub2B.vertices[0].x)
    dy2 = 0.5*(hub2A.vertices[0].y - hub2B.vertices[0].y)
    dz2 = 0.5*(hub2A.vertices[0].z - hub2B.vertices[0].z)
    hub2A.move(-dx2, -dy2, -dz2)
    hub2B.move(dx2, dy2, dz2) 

    update_excluders_position(affected_subunits)
    if check_overlap(affected_subunits, capsid):
        hub1A.move(dx1, dy1, dz1)
        hub1B.move(-dx1, -dy1, -dz1)
        hub2A.move(dx2, dy2, dz2)
        hub2B.move(-dx2, -dy2, -dz2)
        update_excluders_position(affected_subunits)
        return False  

    hub1A_initial_vertices = hub1A.vertices #potential source of bug: deep vs shallow copy
    hub1A.merge_into(hub1B) 

    hub2A_initial_vertices = hub2A.vertices #potential source of bug: deep vs shallow copy
    hub2A.merge_into(hub2B) 

    energy_after = elastic_energy(affected_subunits)
    if wall is not None:
        energy_after+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))        
    dE = energy_after-energy_before
    #K = pr.affinity
    #@K = edgeA.binding_affinity[edgeB.edge_type]
    e_b = edgeA.e_b[edgeB.edge_type]
    #p=(K/pr.l_add**6 *np.exp(-dE))
    #@p = (1.0/pr.l_add**6 *math.exp((-dE+K)/capsid.kT))
    p = 1.0/pr.l_fusion**6 * pr.va * math.exp(-(dE+e_b)/capsid.kT)  
    if (RNG.rand()<p):
        return True

    hub1A = hub1B.split_from(hub1A_initial_vertices, hub1A)
    hub2A = hub2B.split_from(hub2A_initial_vertices, hub2A)

    hub1A.move(dx1, dy1, dz1)
    hub1B.move(-dx1, -dy1, -dz1)
    hub2A.move(dx2, dy2, dz2)
    hub2B.move(-dx2, -dy2, -dz2)  
    update_excluders_position(affected_subunits)  
    #print "EF: ", dE, p
    return False

    #hub1 = hub2.split_from(hub1_initial_vertices)
    #hub1.move(dx, dy, dz)
    #hub2.move(-dx, -dy, -dz)    
#----------------------  

def attempt_edge_fission(capsid, RNG, edge_fission_hubs, wall=None):
    #edge_fission_hubs = get_edge_fission_hubs(capsid)
    if len(edge_fission_hubs)==0:
        return False
    which = RNG.randint(len(edge_fission_hubs))
    hub1, hub2 = edge_fission_hubs[which]
    edges = get_edges_between_hubs(hub1, hub2)
    affected_subunits = list(set([vertex.parent_subunit for vertex in hub1.vertices]+
                            [vertex.parent_subunit for vertex in hub2.vertices]))

    energy_before = elastic_energy(affected_subunits)
    if wall is not None:
        energy_before+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))        

    vertices_left_hub1 = []
    next_hub = hub2
    for i in xrange(1000): #should be while loop; safety to avoid infinite loops
        if i>900:
            raise ValueError('Edge fission screwed up #1')
        edges12 = get_edges_between_hubs(hub1, next_hub)
        #when hit the nearest surface edge
        if len(edges12)==1:
            break
        subunit_left = edges12[0].parent_subunit
        subunit_right = edges12[1].parent_subunit        
        #find subunit_left vertex that's in hub1
        vertices_left_hub1.append(list(set(subunit_left.vertices) & set(hub1.vertices))[0])        
        #there should be two edges from/to this vertex
        edge1 = vertices_left_hub1[-1].edge_in
        edge2 = vertices_left_hub1[-1].edge_out
        possible_hubs = [edge1.vertex_from.parent_hub, edge2.vertex_to.parent_hub]
        possible_hubs.remove(next_hub)
        next_hub = possible_hubs[0]

    vertices_left_hub2 = []
    next_hub = hub1
    for i in xrange(1000): #should be while loop; safety to avoid infinite loops
        if i>900:
            raise ValueError('Edge fission screwed up #2')
        edges12 = get_edges_between_hubs(hub2, next_hub)
        #when hit the nearest surface edge
        if len(edges12)==1:
            break
        subunit_left = edges12[0].parent_subunit
        subunit_right = edges12[1].parent_subunit        
        #find subunit_left vertex that's in hub1
        vertices_left_hub2.append(list(set(subunit_left.vertices) & set(hub2.vertices))[0])        
        #there should be two edges from/to this vertex
        edge1 = vertices_left_hub2[-1].edge_in
        edge2 = vertices_left_hub2[-1].edge_out
        possible_hubs = [edge1.vertex_from.parent_hub, edge2.vertex_to.parent_hub]
        possible_hubs.remove(next_hub)
        next_hub = possible_hubs[0]              

    new_hub1 = hub1.split_from(vertices_left_hub1)
    dx1, dy1, dz1 = (0.5-RNG.rand(3))*pr.l_fusion #*0.5
    hub1.move(dx1, dy1, dz1)
    new_hub1.move(-dx1, -dy1, -dz1)

    new_hub2 = hub2.split_from(vertices_left_hub2)
    dx2, dy2, dz2 = (0.5-RNG.rand(3))*pr.l_fusion #*0.5
    hub2.move(dx2, dy2, dz2)
    new_hub2.move(-dx2, -dy2, -dz2)

    #check overlap and energetics
    update_excluders_position(affected_subunits)
    if check_overlap(affected_subunits, capsid):
        hub1.move(-dx1, -dy1, -dz1)
        new_hub1.move(dx1, dy1, dz1)
        #hub1.merge_into(new_hub1)
        new_hub1.merge_into(hub1)

        hub2.move(-dx2, -dy2, -dz2)
        new_hub2.move(dx2, dy2, dz2)
        #hub2.merge_into(new_hub2)
        new_hub2.merge_into(hub2)
        update_excluders_position(affected_subunits)
        return False

    energy_after = elastic_energy(affected_subunits)
    if wall is not None:
        energy_after+=sum(wall.potential_energy(np.array([v.z for sub in affected_subunits for v in sub.vertices])))        
    dE = energy_after - energy_before
    #K = pr.affinity
    #@K = edges[0].binding_affinity[edges[1].edge_type]
    e_b = edges[0].e_b[edges[1].edge_type]
    #p = 1.0/(K/pr.l_fusion**6) * np.exp(-dE)
    #@p = 1.0/(1.0/pr.l_fusion**6) * math.exp((-dE-K)/capsid.kT)
    p = 1.0/(1.0/pr.l_fusion**6 * pr.va) * math.exp(-(dE-e_b)/capsid.kT)
    if (RNG.rand()<p):
        print "Edge fission p=",p
        print "E_after, E_before", energy_after, energy_before
        return True #accept
   
    hub1.move(-dx1, -dy1, -dz1)
    new_hub1.move(dx1, dy1, dz1)
    #hub1.merge_into(new_hub1)
    new_hub1.merge_into(hub1)

    hub2.move(-dx2, -dy2, -dz2)
    new_hub2.move(dx2, dy2, dz2)
#    hub2.merge_into(new_hub2)
    new_hub2.merge_into(hub2)
    update_excluders_position(affected_subunits)
    return False



def attempt_move_athermal(capsid, RNG):
    #choose with equal probability from hubs to move
    candidates = capsid.hubs
    which = RNG.randint(len(candidates))
    candidate = candidates[which]

    vertex = candidate.vertices[0]
    #!!!!!!!1
    #!!!!!!!!
    #fix it, otherwise will be slow!!!!
    affected_subunits = [v.parent_subunit for v in candidate.vertices]
    #affected_subunits = capsid.subunits
    
    energy_before = elastic_energy(affected_subunits)

    #for now move in a box volume dV = d_max^3
    dx, dy, dz = (0.5-RNG.rand(3))*2.0*pr.d_max

    candidate.move(dx, dy, dz)
    #check for overlaps; this sucks because it scales as number_of_subunits^2

    update_excluders_position(affected_subunits)
    #if check_overlap(affected_subunits, capsid):
    if check_neighbor_overlap(affected_subunits):
        candidate.move(-dx, -dy, -dz)
        update_excluders_position(affected_subunits)
        #pr.rej+=1
        return False

    #!!! should check if vertex left v_add volume
            
    energy_after = elastic_energy(affected_subunits)
    dE = energy_after-energy_before
    if dE<0:
        #pr.acc+=1
        return True
    else:
        candidate.move(-dx, -dy, -dz)
        update_excluders_position(affected_subunits)
        #pr.rej+=1
        return False
