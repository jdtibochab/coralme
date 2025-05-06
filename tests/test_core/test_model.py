import pytest

@pytest.mark.parametrize("attr,expected",
                         [
                            ('optimize',True),
                            ('optimize_windows',True),
                            ('feasibility',True),
                            ('feas_gurobi',True),
                            ('feas_cplex',True),
                            ('construct_lp_problem',True),
                            ('fva',True)
                             ])
def test_model_optimize(attr, expected, shared_builder):
    builder = pytest.shared_builder
    model = builder.me_model
    assert hasattr(model, attr) == expected, "Model does not have {} method".format(attr)
    
@pytest.mark.parametrize("attr,expected",
                         [
                            ('optimize',True),
                            ('optimize_windows',True),
                            ('feasibility',True),
                            ('feas_gurobi',True),
                            ('feas_cplex',True),
                            ('construct_lp_problem',True),
                            ('fva',True)
                             ])
def test_model_optimize_troubleshooted(attr, expected, shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    model = builder.me_model
    assert hasattr(model, attr) == expected, "Model does not have {} method".format(attr)