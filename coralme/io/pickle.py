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

# Custom unpickler to remap old references
class FixDefaultParamsUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "coralme.core.parameters" and name == "DefaultParameters":
            from coralme.core.extended_classes import DefaultParameters
            return DefaultParameters
        # retrocompatibility with pandas <= 1.5.1
        elif module == "pandas.core.indexes.numeric":
            from pandas import Index
            return Index
        return super().find_class(module, name)

def load_pickle_me_model(path, repair = False):
    with open(path, "rb") as infile:
        me = FixDefaultParamsUnpickler(infile).load()
    if repair:
        me.repair()
    return me

def save_pickle_me_model(me, path):
    with open(path, "wb") as outfile:
        pickle.dump(me, outfile)
