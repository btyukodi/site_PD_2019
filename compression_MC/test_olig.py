import jobs
from oligomerization_run_DEBUG import run
import time
import numpy as np
group_id = ('%.10f' % time.time()).replace('.','')
print "group_id ", group_id

theta = 0.0#1.2

##for kT in [0.5, 1.0, 1.5]:
##    for mu in [-3.0, -2.0, -1.0]:
##        for e_b in [4.0, 5.0, 6.0]:
##            for elastic_modulus in [50.0, 100.0]:#, 200.0]:
##                for bending_modulus in [50.0, 100.0]:#, 200.0]:
##                    for ens in range(5):

for kT in [1.0, 1.5]:
    for mu in [-3.0, -2.0, -1.0]:
        for e_b in [3.0, 4.0, 5.0]:#[4.0, 5.0, 6.0]:
            for elastic_modulus in [50.0, 100.0]:#, 200.0]:
                for bending_modulus in [50.0, 100.0]:#, 200.0]:
                    for ens in range(5):                        

                        params = {'global_params':{'ensemble':ens, 'l_fusion':0.5, 'd_max':(2*kT/(elastic_modulus+bending_modulus))**0.5, 'l_add':3*(2*kT/(elastic_modulus+bending_modulus))**0.5, 
                        'tmax':750000/4, 'dtsave':5000, 'kT':kT},
                                 'subunits':{'sub1':{'edges':{'edge1':'type1', 'edge2':'type1', 'edge3':'type1'},
                                                     'excluders':{'exc1':'com'},
                                                     'activity':mu,
                                                     'insertion_rate':0.1,
                                                     'number_of_rotational_configurations':3}
                                            },
                                 'edges':{'type1': {'equilibrium_angle':{'type1':theta, 'type2':0.0, 'type3':0.0},
                                                    'bending_modulus':{'type1':bending_modulus, 'type2':float('inf'), 'type3':float('inf')},
                                                    'binding_affinity':{'type1':e_b, 'type2':0.0, 'type3':0.0},
                                                    'elastic_modulus':elastic_modulus, 
                                                    'equilibrium_length':1.0 },
                                          'type2': {'equilibrium_angle':{'type1':0.0, 'type2':theta, 'type3':0.0},
                                                    'bending_modulus':{'type1':float('inf'), 'type2':bending_modulus, 'type3':float('inf')},
                                                    'binding_affinity':{'type1':0.0, 'type2':e_b, 'type3':0.0},
                                                    'elastic_modulus':elastic_modulus, 
                                                    'equilibrium_length':1.0 },
                                          'type3': {'equilibrium_angle':{'type1':0.0, 'type2':0.0, 'type3':theta},
                                                    'bending_modulus':{'type1':float('inf'), 'type2':float('inf'), 'type3':bending_modulus},
                                                    'binding_affinity':{'type1':0.0, 'type2':0.0, 'type3':e_b},
                                                    'elastic_modulus':elastic_modulus, 
                                                    'equilibrium_length':1.0 }                        
                                                    
                                         },
                                 'excluders':{'com':{'overlaps':{'com':True}, 'radius':0.2}
                                             }
                                          }
                        jobs.stage_runs_in_db(params, run, "oligo test - isotropic F=0", group_id, host='btyukodi-MS-7B09')
#print 3*(2*kT/(elastic_modulus+bending_modulus))**0.5
#jobs.stage_runs_in_db(params, run, "oligo test", group_id)
jobs.submit_group(group_id)
print group_id