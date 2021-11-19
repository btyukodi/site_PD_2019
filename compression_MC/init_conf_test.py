import jobs
from run import run
import time
import numpy as np
from bson.objectid import ObjectId
import parameters as pr
group_id = ('%.10f' % time.time()).replace('.','')
print "group_id ", group_id


#ring
init_conf_id = ObjectId("5c2e9d27946df6b86dedf275")
#init_conf_id = ObjectId('5c421c15946df6e447399c99')#ObjectId('5c3f5899946df64cf998a702')
params = pr.get_db_collection().find_one({"_id":init_conf_id})
params['parameters']['global_params']['dtsave'] = 200


for timestep in [46000]:
    for A in [13.0]:
        for zs in [-7.0]:
            for key in params['parameters']['edges'].keys():
                params['parameters']['edges'][key]['binding_affinity'][key] = A

            #params['parameters']['subunits']['sub1']['insertion_rate'] = 0.005
            params['parameters']['subunits']['sub1']['activity'] = zs
            for ens in range(1):
                init_config = {"id": init_conf_id, "timestep":timestep}
                jobs.stage_runs_in_db(params=params['parameters'], run_function=run, 
                    comments='Trumpet ring seeded WEDGE TEST: (theta, theta, -theta/2+delta). theta=0.28 (N=13 in this geometry)',
                    group_id=group_id, init_config=init_config)


#zigzag
#init_conf_id = ObjectId("5c2e9d27946df6b86dedf2b9")
##init_conf_id = ObjectId('5c4227e4946df6e447399c9d')
##params = pr.get_db_collection().find_one({"_id":init_conf_id})
##params['parameters']['global_params']['dtsave'] = 200
##
##
##for timestep in [40000, 62000, 84000]:
##    for A in np.linspace(1,10,10):
##
##        for key in params['parameters']['edges'].keys():
##            params['parameters']['edges'][key]['binding_affinity'][key] = A
##        for ens in range(1):
##            init_config = {"id": init_conf_id, "timestep":timestep}
##            jobs.stage_runs_in_db(params=params['parameters'], run_function=run, 
##                comments='Trumpet zigzag seeded DISC: (delta, delta, theta). theta=0.483 (N=13 in this geometry)',
##                group_id=group_id, init_config=init_config)        