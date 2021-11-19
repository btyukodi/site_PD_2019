import jobs
from run import run
import time
import numpy as np
group_id = ('%.10f' % time.time()).replace('.','')
print "group_id ", group_id

#define a template indicating parameters that will be changed
params_template = {'global_params':{'ensemble':1, 'l_fusion':0.5, 'd_max':0.05, 'l_add':0.5, 'tmax':1000000, 'dtsave':2000},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type2', 'edge3':'type3'},
                             'excluders':{'exc1':'com'},
                             'activity':-7,
                             'insertion_rate':0.0005,
                             'number_of_rotational_configurations':3}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':'#THETA1', 'type2':0.0, 'type3':0.0},
                            'bending_modulus':{'type1':300.0, 'type2':float('inf'), 'type3':float('inf')},
                            'binding_affinity':{'type1':'#AFFINITY1', 'type2':0.0, 'type3':0.0},
                            'elastic_modulus':300.0, 
                            'equilibrium_length':1.0 },
                  'type2': {'equilibrium_angle':{'type1':0.0, 'type2':'#THETA2', 'type3':0.0},
                            'bending_modulus':{'type1':float('inf'), 'type2':300.0, 'type3':float('inf')},
                            'binding_affinity':{'type1':0.0, 'type2':'#AFFINITY2', 'type3':0.0},
                            'elastic_modulus':300.0, 
                            'equilibrium_length':1.0 },
                  'type3': {'equilibrium_angle':{'type1':0.0, 'type2':0.0, 'type3':'#THETA3'},
                            'bending_modulus':{'type1':float('inf'), 'type2':float('inf'), 'type3':300.0},
                            'binding_affinity':{'type1':0.0, 'type2':0.0, 'type3':'#AFFINITY3'},
                            'elastic_modulus':300.0, 
                            'equilibrium_length':1.0 }                        
                            
                 },
         'excluders':{'com':{'overlaps':{'com':True}, 'radius':0.2}
                     }
                  }



AL = np.linspace(0,15,16)
AV = np.linspace(0,15,16)
A = []

for al in AL:
    for av in AV:
        if (al<8.0) and (av<8.0):
            continue
        else:
            A.append((al,av))

#------- rings --------------------------
theta = 0.28
params = []

for al, av in A:
    params.append((theta, theta, -theta/2, av, av, al))

#give the values for the parameters to be changed
params_values = {'names':['#THETA1','#THETA2' ,'#THETA3', '#AFFINITY1', '#AFFINITY2', '#AFFINITY3'],
                 'values':params}
                             
#stage the set
set_id = jobs.stage_set(params_template, params_values, run, comments='Stevens cylinder ring: (theta, theta, -theta/2). theta=0.28 (N=13 in this geometry) and different binding affinities', group_id=group_id)


#------- zigzag --------------------------
theta = 0.483
params = []

for al, av in A:
    params.append((0.0, 0.0, theta, al, al, av))

#give the values for the parameters to be changed
params_values = {'names':['#THETA1','#THETA2' ,'#THETA3', '#AFFINITY1', '#AFFINITY2', '#AFFINITY3'],
                 'values':params}
                             
#stage the set
set_id = jobs.stage_set(params_template, params_values, run, comments='Stevens cylinder zigzag: (0, 0, theta). theta=0.483 (N=13 in this geometry) and different binding affinities', group_id=group_id)

print "group_id: ", group_id
#jobs.submit_group(group_id)
