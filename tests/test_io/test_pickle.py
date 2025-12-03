import pytest
from coralme.io.pickle import load_pickle_me_model, save_pickle_me_model

def test_save_and_load_pickle_me_model(shared_builder_troubleshooted, tmp_path):
    builder1 = pytest.shared_builder_troubleshooted
    file_path = tmp_path / "{}.pkl".format(builder1.me_model.id)

    # Save
    save_pickle_me_model(builder1.me_model, file_path)
    assert file_path.exists()

    # Load
    model = load_pickle_me_model(file_path)
    model.optimize()
    assert model.solution is not None
