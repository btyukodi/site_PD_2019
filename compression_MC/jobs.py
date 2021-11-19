import platform
import parameters as pr
from multiprocessing import Pool, current_process
import os
import timeit
from bson.objectid import ObjectId
import data_process as dp
import Utils as ut
import config

def stage_runs_in_db(params, run_function, comments='', group_id=None, set_id=None, init_config=None, host=None):
    """Adds an entry to the database, i.e. stages a run that can be launched later on"""
    collection =pr. get_db_collection()
    run_entry = {'parameters':params, 'status':'staged'}
    _id = collection.insert_one(run_entry).inserted_id
    data_location = '/'.join(list(map(''.join, zip(*[iter(str(_id))]*4))))
    collection.update_one({"_id":_id},{"$set":{'data_location':data_location, 'init_config':init_config}})
    #from run_MS import run
    import inspect
    source = inspect.getsource(run_function)
    collection.update_one({"_id":_id},{"$set":{'run_function':run_function.__name__,'source':source, 'comments':comments, 'group_id':group_id, 'set_id':set_id, 'host':host}}) 

    return _id

#run one simulation
def run_one(run_id):
    """Runs a single simulation staged in the database"""
    collection = pr.get_db_collection()
    staged_host = collection.find_one({"_id":run_id})['host']

    #data_root_conf = pr.get_data_root_config()
    data_root = config.data_root#data_root_conf.find_one({"host_from":staged_host, "host_to":staged_host})['data_root']    


    
    run_source = collection.find_one({"_id":run_id})['source']
    init_config_entry = collection.find_one({"_id":run_id})['init_config']
    if init_config_entry:
        init_config_id = init_config_entry['id']
        init_config_timestep = init_config_entry['timestep']
        #should check here if old parameter's key structure matches the new one
        init_config = dp.get_single_capsid_snapshot(run_id=init_config_id, timestep=init_config_timestep)
    else:
        init_config=None

    #define the run() function from the database
    exec(run_source)
    start = timeit.default_timer()
    #@collection.update_one({"_id":run_id},{"$set":{'status':{'started':start}, 'host':platform.node(), 'data_root':data_root}}, upsert = True)
    run_params = pr.load_parameters_from_db(run_id)
    run_function = run_params['run_function']
    run = locals()[run_function]

    run(run_params=run_params, data_root=data_root, init_config=init_config, run_id=run_id)
    stop = timeit.default_timer()
    #@collection.update_one({"_id":run_id},{"$set":{'status':{'started':start, 'ended':stop, 'DeltaT':(stop-start)/60}}}, upsert = True)  
    print 'Time: ', stop - start, (stop-start)/60  

#submit a single job
def submit_job(run_id):
    """Submits a single simulation from the database as a job. Rules for various hosts have to be added here."""


    host = platform.node()
    if host=='btyukodi-MS-7B09':
        run_one(run_id)
    if host=='hpcc.brandeis.edu':
        #data_root = '/work/btyukodi/data/'
        data_root = config.data_root
        f = open('sbatch_job_template.sh','r')
        sbatch_script = f.read()
        f.close()
        #sbatch_script = create_SBATCH_script(run_id)
        sbatch_script = sbatch_script.replace('DATAROOT', data_root)
        sbatch_script = sbatch_script.replace('RUNID', str(run_id))        
        sbatch_script = sbatch_script.replace('JOBNAME', str(run_id))
        sbatch_script = sbatch_script.replace('OUTPUTFILE', data_root+'jobs/'+str(run_id)+'.out')
        sbatch_script = sbatch_script.replace('ERRORFILE', data_root+'jobs/'+str(run_id)+'.err')
        #maybe save the jobfile itself as well?
        #os.system('sbatch < '+sbatch_script)
        ut.ensure_dir(data_root+'jobs/')
        f = open(data_root+'jobs/'+str(run_id)+'.sh','w')
        f.write(sbatch_script)
       	f.close()
       	os.system('sbatch '+data_root+'jobs/'+str(run_id)+'.sh')

    if host in ['login-00','login-01']:
    	#data_root = '/scratch/btyukodi/data/'
        #data_root_conf = pr.get_data_root_config()
        data_root = config.data_root#data_root_conf.find_one({"host_from":host, "host_to":host})['data_root']  


    	f = open('discovery_job_template.sh','r')
        sbatch_script = f.read()
        f.close()
        #sbatch_script = create_SBATCH_script(run_id)
        #sbatch_script = sbatch_script.replace('DATAROOT', data_root)
        sbatch_script = sbatch_script.replace('RUNID', str(run_id))        
        sbatch_script = sbatch_script.replace('JOBNAME', str(run_id))
        sbatch_script = sbatch_script.replace('OUTPUTFILE', data_root+'jobs/'+str(run_id)+'.out')
        sbatch_script = sbatch_script.replace('ERRORFILE', data_root+'jobs/'+str(run_id)+'.err')
        #maybe save the jobfile itself as well?
        #os.system('sbatch < '+sbatch_script)
        ut.ensure_dir(data_root+'jobs/')
        f = open(data_root+'jobs/'+str(run_id)+'.sh','w')
        f.write(sbatch_script)
       	f.close()
       	os.system('sbatch '+data_root+'jobs/'+str(run_id)+'.sh')    	


def submit_process_job(run_id,  frequency=1, process_functions=[]):
    """Submits data processing jobs. process_functions is a list of functions to be called at each snapshot or #freq snapshot. Should be merged with submit_jobs """       	
    host = platform.node()
    if host=='btyukodi-MS-7B09':
    	dp.process(run_id=run_id)
    if host=='hpcc.brandeis.edu':
        data_root = config.data_root#'/work/btyukodi/data/'
        f = open('sbatch_proc_template.sh','r')
        sbatch_script = f.read()
        f.close()
        #sbatch_script = create_SBATCH_script(run_id)
        sbatch_script = sbatch_script.replace('DATAROOT', data_root)
        sbatch_script = sbatch_script.replace('RUNID', str(run_id))        
        sbatch_script = sbatch_script.replace('JOBNAME', str(run_id))
        sbatch_script = sbatch_script.replace('OUTPUTFILE', data_root+'jobs/proc_'+str(run_id)+'.out')
        sbatch_script = sbatch_script.replace('ERRORFILE', data_root+'jobs/proc_'+str(run_id)+'.err')
        #maybe save the jobfile itself as well?
        #os.system('sbatch < '+sbatch_script)
        ut.ensure_dir(data_root+'jobs/')
        f = open(data_root+'jobs/proc_'+str(run_id)+'.sh','w')
        f.write(sbatch_script)
       	f.close()
       	os.system('sbatch '+data_root+'jobs/proc_'+str(run_id)+'.sh')

    if host in ['login-00','login-01']:
    	#data_root = '/scratch/btyukodi/data/'
        #data_root_conf = pr.get_data_root_config()
        #---- host_to should be the one from staged-runs['host']
        data_root = config.data_root#data_root_conf.find_one({"host_from":host, "host_to": host})['data_root']  

    	f = open('discovery_proc_template.sh','r')
        sbatch_script = f.read()
        f.close()
        #sbatch_script = create_SBATCH_script(run_id)
        #sbatch_script = sbatch_script.replace('DATAROOT', data_root)
        sbatch_script = sbatch_script.replace('RUNID', str(run_id))        
        sbatch_script = sbatch_script.replace('JOBNAME', str(run_id))
        sbatch_script = sbatch_script.replace('OUTPUTFILE', data_root+'jobs/proc_'+str(run_id)+'.out')
        sbatch_script = sbatch_script.replace('ERRORFILE', data_root+'jobs/proc_'+str(run_id)+'.err')
        #maybe save the jobfile itself as well?
        #os.system('sbatch < '+sbatch_script)
        ut.ensure_dir(data_root+'jobs/')
        f = open(data_root+'jobs/proc_'+str(run_id)+'.sh','w')
        f.write(sbatch_script)
       	f.close()
       	os.system('sbatch '+data_root+'jobs/proc_'+str(run_id)+'.sh')  

#submit a group
def submit_group(group_id):
    """Submits a group of simulations from the database as jobs. Rules for various hosts have to be added here."""
    host = platform.node()

    if host=='btyukodi-MS-7B09':
        pool = Pool(processes=28)
        cursor = pr.get_db_collection().find({'group_id':group_id})
        for document in cursor:
            pool.apply_async(run_one, [document['_id']])
        pool.close()
        pool.join()
        #apply async run_one(run_id, 'data/')
    if host=='hpcc.brandeis.edu':  
        cursor = pr.get_db_collection().find({'group_id':group_id})
        for document in cursor:
            submit_job(document['_id'])
    if host in ['login-00','login-01']:
        cursor = pr.get_db_collection().find({'group_id':group_id})
        for document in cursor:
            submit_job(document['_id'])            

#submit a set
def submit_set(set_id):
    """Submits a set of simulations from the database as jobs. Rules for various hosts have to be added here."""
    host = platform.node()

    if host=='btyukodi-MS-7B09':
        pool = Pool()
        cursor = pr.get_db_collection().find({'set_id':set_id})
        for document in cursor:
            pool.apply_async(run_one, [document['_id']])
        pool.close()
        pool.join()
        #apply async run_one(run_id, 'data/')
    if host=='hpcc.brandeis.edu':  
        cursor = pr.get_db_collection().find({'set_id':set_id})
        for document in cursor:
            submit_job(document['_id'])

#stage a set of runs
#a set is a bunch of runs with most parameters kept fixed and a few varying across the runs
def stage_set(params_template, params_values, run_function, comments='', group_id=None, host=None):
    """Stage a set of runs.
       A set is a bunch of runs with most parameters kept fixed and a few varying across the runs.
       params_template is like a params dictionary, except that some values are replaced by '#PARAM1_VALUE' strings.
       params_values then sets values for these parameters like params_values={'names':['#PARAM1_VALUE', 'PARAM2_VALUE'], values:[(2, 5), (2,7), ...]}"""
    run_collection = pr.get_db_collection()
    db = run_collection.database
    set_collection = db['staged-sets']

    set_id = set_collection.insert_one({'params_template':params_template, 'params_values':params_values, 'comments':comments}).inserted_id
    param_names = params_values['names']
    rids = []
    for values in params_values['values']:
        params = params_template
        for ix in range(len(values)):
            params = replace_in_template(params, param_names[ix], values[ix])
        rid = stage_runs_in_db(params, run_function, comments, group_id, set_id, host=host)
        rids.append(rid)
    set_collection.update_one({'_id':set_id}, {'$set':{'run_ids':rids}})
    return set_id    

def remove_run(run_id):
	pass

def remove_set(set_id):
	pass

def remove_group(group_id):
	pass


def replace_in_template(template, to_replace, replacement):
    """Given a dictionary template, function replaces all values to_replace to the value replacement"""
    def replace_in_template_rec(dict_data, to_replace, replacement):
        for key, value in dict_data.items():
            if isinstance(value, dict):
                replace_in_template_rec(value, to_replace, replacement)
            else:
                if (dict_data[key]==to_replace):
                    dict_data[key] = replacement
                    #print key,replacement
    import copy
    dict_data = copy.deepcopy(template)
    replace_in_template_rec(dict_data, to_replace, replacement)
    return dict_data      
