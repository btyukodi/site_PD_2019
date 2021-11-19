import jobs
from run import run
import time
group_id = ('%.10f' % time.time()).replace('.','')
print "group_id ", group_id

#define a template indicating parameters that will be changed
params_template = {'global_params':{'ensemble':1, 'l_fusion':0.5, 'd_max':0.05, 'l_add':0.5, 'tmax':1000000, 'dtsave':500},
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

t13 = 0.28 
params = []

a1 = 24.0
a2 = 0.0
params.append((t13, t13, -t13/2.0, a1, a1, a2))



a1 = 0.0
a2 = 24.0
params.append((t13, t13, -t13/2.0, a1, a1, a2))




#affinity_values = [11.0, 13.0]#[11.0, 11.5, 12.0, 12.5, 13.0, 13.5]#[11.0, 12.0, 13.0, 14.0, 15.0]
#for a1 in affinity_values:
#    for a2 in affinity_values:
#        for i in range(15):
#            params.append((t13, t13, -t13/2.0, a1, a1, a2))


#give the values for the parameters to be changed
params_values = {'names':['#THETA1','#THETA2' ,'#THETA3', '#AFFINITY1', '#AFFINITY2', '#AFFINITY3'],
                 'values':params}
                             
#stage the set
set_id = jobs.stage_set(params_template, params_values, run, comments='TEST', group_id=group_id)



#submit the set for run
#jobs.submit_set(set_id)
#jobs.submit_group(group_id)
