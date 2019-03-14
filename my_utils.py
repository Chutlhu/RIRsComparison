import numpy as np
import pickle as pkl
import hdf5storage
import json
import progressbar as pb

def save_to_pickle(obj, filename ):
    with open(filename, 'wb') as f:
        pkl.dump(obj, f, pkl.HIGHEST_PROTOCOL)

def write_to_json(json_obj, filename):
    with open(filename, 'w') as f:
        json.dump(json_obj, f)

def load_from_pickle(filename):
    with open(filename, 'rb') as f:
        return pkl.load(f)

def load_matfile(filename, var = None):
    print('Loading: ' + filename + ' ...')
    mat = hdf5storage.loadmat(filename)
    if not var is None:
        mat = mat[var]
    print('done')
    return mat


#define progress timer class
class progress_timer:

    def __init__(self, n_iter, description="Something"):
        self.n_iter         = n_iter
        self.iter           = 0
        self.description    = description + ': '
        self.timer          = None
        self.initialize()

    def initialize(self):
        #initialize timer
        widgets = [
            self.description, pb.Percentage(),
            ' ', pb.Bar(marker='#', left='[', right=']'),
            ' ', pb.ETA(),
        ]
        self.timer = pb.ProgressBar(widgets=widgets, maxval=self.n_iter).start()

    def update(self, q=1):
        #update timer
        self.timer.update(self.iter)
        self.iter += q

    def finish(self):
        #end timer
        self.timer.finish()
