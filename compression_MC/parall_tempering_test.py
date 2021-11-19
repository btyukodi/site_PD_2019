import data_process as dp
import parameters as pr
from bson.objectid import ObjectId
from umbrella_run_parall_tempering import run
import time
import jobs
from bson.objectid import ObjectId

group_id = ('%.10f' % time.time()).replace('.','')
print "group_id ", group_id

#params = pr.load_parameters_from_db(ObjectId('5cfa910d946df6282a4b4f26'))
params = pr.load_parameters_from_db(ObjectId('5d41ea01946df6283be90e62'))


print params.keys()
params['parameters']['global_params']['kT_values'] = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2]


_id = jobs.stage_runs_in_db(params['parameters'], run, "parall tempering test", group_id, host='btyukodi-MS-7B09')
jobs.submit_job(ObjectId(_id))



#run(params, data_root='data/parall_tempering/')
