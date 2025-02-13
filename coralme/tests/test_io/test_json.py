import pytest
from importlib.resources import files
from coralme.io.json import load_json_me_model,save_json_me_model

dir = str(files("coralme"))
def test_save_json_me_model(shared_builder):
    save_json_me_model(shared_builder.me_model,"{}/tests/data/{}.json".format(dir,shared_builder.me_model.id))

def test_load_json_me_model(shared_builder):
    model = load_json_me_model("{}/tests/data/{}.json".format(dir,shared_builder.me_model.id))
    model.optimize()
    assert model.solution is not None

def test_old_json_me_model():
    model = load_json_me_model("{}/tests/data/VERSION-1.0.0-EXAMPLE.json".format(dir))
    model.optimize()
    assert model.solution is not None
