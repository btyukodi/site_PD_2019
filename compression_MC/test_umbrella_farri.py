import jobs
from umbrella_run import run
import time
import numpy as np
group_id = ('%.10f' % time.time()).replace('.','')
print "group_id ", group_id

elastic_modulus=10.0
kT = 0.5
bending_modulus = 300.0
#define a template indicating parameters that will be changed
params_template = {'global_params':{'ensemble':2, 'l_fusion':0.5, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                                    'l_add':3*(2*kT/min(elastic_modulus,bending_modulus))**0.5, 'tmax':10000000, 'dtsave':100, 'kT':kT,},
         'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type2', 'edge3':'type3'},
                             'excluders':{'exc1':'com'},
                             'activity':'#ACTIVITY',
                             'insertion_rate':0.0005,
                             'number_of_rotational_configurations':3}
                    },
         'edges':{'type1': {'equilibrium_angle':{'type1':'#THETA1', 'type2':0.0, 'type3':0.0},
                            'bending_modulus':{'type1':bending_modulus, 'type2':float('inf'), 'type3':float('inf')},
                            'binding_affinity':{'type1':'#AFFINITY1', 'type2':0.0, 'type3':0.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },
                  'type2': {'equilibrium_angle':{'type1':0.0, 'type2':'#THETA2', 'type3':0.0},
                            'bending_modulus':{'type1':float('inf'), 'type2':bending_modulus, 'type3':float('inf')},
                            'binding_affinity':{'type1':0.0, 'type2':'#AFFINITY2', 'type3':0.0},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 },
                  'type3': {'equilibrium_angle':{'type1':0.0, 'type2':0.0, 'type3':'#THETA3'},
                            'bending_modulus':{'type1':float('inf'), 'type2':float('inf'), 'type3':bending_modulus},
                            'binding_affinity':{'type1':0.0, 'type2':0.0, 'type3':'#AFFINITY3'},
                            'elastic_modulus':elastic_modulus, 
                            'equilibrium_length':1.0 }                        
                            
                 },
         'excluders':{'com':{'overlaps':{'com':True}, 'radius':0.2}
                     }
                  }



#A = np.linspace(4.5,6.5,11)
#A = np.linspace(0.5,4.5,21)
#A = np.linspace(3.0,7.0, 11)
A = np.linspace(2.0,15.0, 61)
zs =[-1.0]# [-4.0, -3.0, -2.5, -1.0]#, -2.0]
#------- rings --------------------------
theta = 0.367
params = []

for a in A:
    for z in zs:
      for e in range(2):
          #params.append((theta, theta, -theta, 4*a, 4*a, a, z))
          params.append((theta, theta, -theta, 4*a, 4*a, a, z))          
params_values = {'names':['#THETA1','#THETA2' ,'#THETA3', '#AFFINITY1', '#AFFINITY2', '#AFFINITY3', '#ACTIVITY'],
                 'values':params}
set_id = jobs.stage_set(params_template, params_values, run, comments='Umbrella test: Farri trumpet fixed: (4*aff, 4*aff, aff), (theta, theta, -theta)', group_id=group_id, host='login-00')


##for a in [10.0]:
##    for z in [-1.0]:
##      for e in [1]:
##          params.append((theta, theta, -theta, 4*a, 4*a, a, z))
##params_values = {'names':['#THETA1','#THETA2' ,'#THETA3', '#AFFINITY1', '#AFFINITY2', '#AFFINITY3', '#ACTIVITY'],
##                 'values':params}
##set_id = jobs.stage_set(params_template, params_values, run, comments='dummy', group_id=group_id, host='WS')




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
