import pytest
import sys
from coralme.core.model import MEModel
from coralme.core.optimization import optimize
from coralme.core.optimization import feas_gurobi
from coralme.core.optimization import feasibility
from coralme.core.optimization import fva
from coralme.core.optimization import _get_evaluated_nlp
from coralme.core.optimization import compute_solution_error
from coralme.core.optimization import construct_lp_problem

def test_optimize(shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    model = builder.me_model
    optimize(model,solver = "qminos")
    assert model.solution is not None
    try:
        import gurobipy
        optimize(model,solver = "gurobi")
        assert model.solution is not None
    except:
        pass
    if sys.version_info == (3,10):
        model.optimize(solver = 'cplex') # needs CPLEX runtime installed manually (only python 3.10?)
        assert model.solution is not None

def test_feasibility(shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    model = builder.me_model
    feasibility(model)
    assert model.solution is not None
    try:
        import gurobipy
        feas_gurobi(model)
        assert model.solution is not None
    except:
        pass
    if sys.version_info == (3,10):
        model.feas_cplex() # needs CPLEX runtime installed manually (only python 3.10?)
        assert model.solution is not None

def test_fva(shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    model = builder.me_model
    reaction_list = model.reactions.query("FWD")[:10]
    fva(model,
        reaction_list, 0.99, mu_fixed = None, objective = 'biomass_dilution',
        max_mu = 2.8100561374051836, min_mu = 0., maxIter = 100, lambdify = True,
        tolerance = 1e-6, precision = 'quad', verbose = True)

def test_optimize_m_model(shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    model = builder.me_model.gem
    optimize(model,solver = "qminos")
    assert model.solution is not None
    try:
        import gurobipy
        optimize(model,solver = "gurobi")
        assert model.solution is not None
    except:
        pass
    if sys.version_info == (3,10):
        model.optimize(solver = 'cplex') # needs CPLEX runtime installed manually (only python 3.10?)
        assert model.solution is not None

def test_optimize_m_model_from_cobra(shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    model = MEModel.from_cobra(builder.me_model.gem)
    optimize(model,solver = "qminos")
    assert model.solution is not None
    try:
        import gurobipy
        optimize(model,solver = "gurobi")
        assert model.solution is not None
    except:
        pass
    if sys.version_info == (3,10):
        model.optimize(solver = 'cplex') # needs CPLEX runtime installed manually (only python 3.10?)
        assert model.solution is not None

def test_get_evaluated_nlp(shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    _get_evaluated_nlp(builder.me_model, keys = { builder.me_model.mu.magnitude : builder.me_model.solution.objective_value })

def test_compute_solution_error(shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    builder.me_model.optimize(solver = 'qminos')
    compute_solution_error(builder.me_model, keys = { builder.me_model.mu.magnitude : builder.me_model.solution.objective_value })

def test_construct_lp_problem_me_model(shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    construct_lp_problem(builder.me_model, lambdify = False, per_position = False, as_dict = False, statistics = False)
    construct_lp_problem(builder.me_model, lambdify = True, per_position = False, as_dict = False, statistics = False)
    construct_lp_problem(builder.me_model, lambdify = True, per_position = True, as_dict = False, statistics = False)
    construct_lp_problem(builder.me_model, lambdify = True, per_position = True, as_dict = True, statistics = False)
    construct_lp_problem(builder.me_model, lambdify = True, per_position = True, as_dict = True, statistics = True)

def test_construct_lp_problem_m_model(shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    construct_lp_problem(builder.me_model.gem, lambdify = False, per_position = False, as_dict = False, statistics = False)
    construct_lp_problem(builder.me_model.gem, lambdify = True, per_position = False, as_dict = False, statistics = False)
    construct_lp_problem(builder.me_model.gem, lambdify = True, per_position = True, as_dict = False, statistics = False)
    construct_lp_problem(builder.me_model.gem, lambdify = True, per_position = True, as_dict = True, statistics = False)
    construct_lp_problem(builder.me_model.gem, lambdify = True, per_position = True, as_dict = True, statistics = True)

def test_construct_lp_problem_m_model_from_cobra(shared_builder_troubleshooted):
    builder = pytest.shared_builder_troubleshooted
    model = MEModel.from_cobra(builder.me_model.gem)
    construct_lp_problem(model, lambdify = False, per_position = False, as_dict = False, statistics = False)
    construct_lp_problem(model, lambdify = True, per_position = False, as_dict = False, statistics = False)
    construct_lp_problem(model, lambdify = True, per_position = True, as_dict = False, statistics = False)
    construct_lp_problem(model, lambdify = True, per_position = True, as_dict = True, statistics = False)
    construct_lp_problem(model, lambdify = True, per_position = True, as_dict = True, statistics = True)
