import pickle

from cobra.core.gene import GPR

_orig_init = GPR.__init__

def patched_init(self, *args, **kwargs):
    # Allow extra positional arguments
    if len(args) > 1:
        # Python 3.13 restored a bogus tuple of args
        # we only want the first one (expression)
        args = (args[0],)
    _orig_init(self, *args, **kwargs)

GPR.__init__ = patched_init

def load_pickle_me_model(path):
    with open(path, "rb") as infile:
        return pickle.load(infile)

def save_pickle_me_model(me, path):
    with open(path, "wb") as outfile:
        pickle.dump(me, outfile)
