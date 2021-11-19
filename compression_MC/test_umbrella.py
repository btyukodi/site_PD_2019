import jobs
from umbrella_run import run
import time
import numpy as np
group_id = ('%.10f' % time.time()).replace('.','')
print "group_id ", group_id

#define a template indicating parameters that will be changed
params_template = {'global_params':{'ensemble':2, 'l_fusion':0.5, 'd_max':0.05, 'l_add':0.5, 'tmax':10000000, 'dtsave':2000},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type2', 'edge3':'type3'},
                             'excluders':{'exc1':'com'},
                             'activity':'#ACTIVITY',
                             'insertion_rate':0.005,
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



A = np.linspace(4.5,6.5,11)
zs = [-3.0, -2.0, -1.0]
#------- rings --------------------------
theta = 0.367
params = []

for a in A:
    for z in zs:
      for e in range(2):
          params.append((theta, theta, -2*theta, a, a, a, z))
params_values = {'names':['#THETA1','#THETA2' ,'#THETA3', '#AFFINITY1', '#AFFINITY2', '#AFFINITY3', '#ACTIVITY'],
                 'values':params}
set_id = jobs.stage_set(params_template, params_values, run, comments='Umbrella test: Chips: (theta, theta, -2theta)', group_id=group_id)
##for a in A:
##    for z in zs:
##    	for e in range(2):
##        	params.append((theta, theta, -theta, a, a, a, z))
##
###give the values for the parameters to be changed
##params_values = {'names':['#THETA1','#THETA2' ,'#THETA3', '#AFFINITY1', '#AFFINITY2', '#AFFINITY3', '#ACTIVITY'],
##                 'values':params}
##                             
###stage the set
##set_id = jobs.stage_set(params_template, params_values, run, comments='Umbrella test: Trumpet ring: (theta, theta, -theta/2+delta). theta=0.533 (N=7 in this geometry)', group_id=group_id)
##
##
##params = []
##
##for a in A:
##    for z in zs:
##    	for e in range(2):    	
##        	params.append((theta, theta, -theta/2, a, a, a, z))
##
###give the values for the parameters to be changed
##params_values = {'names':['#THETA1','#THETA2' ,'#THETA3', '#AFFINITY1', '#AFFINITY2', '#AFFINITY3', '#ACTIVITY'],
##                 'values':params}
##                             
###stage the set
##set_id = jobs.stage_set(params_template, params_values, run, comments='Umbrella test: cylinder ring: (theta, theta, -theta/2+delta). theta=0.533 (N=7 in this geometry)', group_id=group_id)



print "group_id: ", group_id
#jobs.submit_group(group_id)
