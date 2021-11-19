import data_process as dp
import parameters as pr
from bson.objectid import ObjectId
import umbrella_run_parall

seed_capsids_id = ObjectId('5cfa910d946df6282a4b4f26')
seed_path = dp.get_data_location(seed_capsids_id)[1]
print seed_path

params = pr.load_parameters_from_db(seed_capsids_id)


params['parameters']['seed_capsids']=ObjectId('5cfa910d946df6282a4b4f26')
#params['parameters']['global_params']['dtsave']=50


umbrella_run_parall.run(params, data_root='data/')
