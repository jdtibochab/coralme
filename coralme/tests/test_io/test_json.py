import pytest
from importlib.resources import files
from coralme.io.json import load_json_me_model,save_json_me_model

dir = str(files("coralme"))
def test_save_json_me_model(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    save_json_me_model(builder1.me_model,"{}/tests/data/base_model/{}.json".format(dir,builder1.me_model.id))

def test_load_json_me_model(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    model = load_json_me_model("{}/tests/data/base_model/{}.json".format(dir,builder1.me_model.id))
    model.optimize()
    assert model.solution is not None

def test_old_json_me_model():
    model = load_json_me_model("{}/tests/data/VERSION-1.0.0-EXAMPLE.json".format(dir))
    model.optimize()
    assert model.solution is not None
