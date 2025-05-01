import pytest
from importlib.resources import files
from coralme.io.pickle import load_pickle_me_model,save_pickle_me_model
# dir = str(files("coralme"))
def test_save_pickle_me_model(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    save_pickle_me_model(builder1.me_model,"./tests/data/base_model/{}.pkl".format(builder1.me_model.id))

def test_load_pickle_me_model(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    model = load_pickle_me_model("./tests/data/base_model/{}.pkl".format(builder1.me_model.id))
    model.optimize()
    assert model.solution is not None