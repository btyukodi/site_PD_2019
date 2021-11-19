import os, errno
from collections import OrderedDict

def ensure_dir(dirname):
    """
    Ensure that a named directory exists; if it does not, attempt to create it.
    """
    try:
        os.makedirs(dirname)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise


def create_data_folder(params, root=''):
	folder = root+'data/'+'/'.join([k+str(v) for k, v in params.iteritems()])
	ensure_dir(folder)
	return folder
