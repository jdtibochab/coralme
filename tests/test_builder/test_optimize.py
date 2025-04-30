import sys
import coralme
import pytest

def test_optimize_me_model(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    builder1.me_model.optimize(solver = 'qminos')
    assert builder1.me_model.solution is not None
    builder1.me_model.optimize(solver = 'gurobi')
    assert builder1.me_model.solution is not None
    if sys.version_info == (3,10):
        builder1.me_model.optimize(solver = 'cplex') # needs CPLEX runtime installed manually (only python 3.10?)
        assert builder1.me_model.solution is not None

def test_optimize_m_model(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    coralme.core.optimization.optimize(builder1.me_model.gem, solver = 'qminos')
    assert builder1.me_model.solution is not None
    coralme.core.optimization.optimize(builder1.me_model.gem, solver = 'gurobi')
    assert builder1.me_model.solution is not None
    if sys.version_info == (3,10):
        coralme.core.optimization.optimize(builder1.me_model.gem, solver = 'cplex') # needs CPLEX runtime installed manually (only python 3.10?)
        assert builder1.me_model.solution is not None

def test_optimize_m_model_from_cobra(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    builder1.me_model.gem = coralme.core.model.MEModel.from_cobra(builder1.me_model.gem)
    coralme.core.optimization.optimize(builder1.me_model.gem, solver = 'qminos')
    assert builder1.me_model.solution is not None
    coralme.core.optimization.optimize(builder1.me_model.gem, solver = 'gurobi')
    assert builder1.me_model.solution is not None
    if sys.version_info == (3,10):
        bcoralme.core.optimization.optimize(builder1.me_model.gem, solver = 'cplex') # needs CPLEX runtime installed manually (only python 3.10?)
        assert builder1.me_model.solution is not None

def test_get_evaluated_nlp(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    coralme.core.optimization._get_evaluated_nlp(builder1.me_model, keys = { builder1.me_model.mu.magnitude : builder1.me_model.solution.objective_value })

def test_compute_solution_error(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    builder1.me_model.optimize(solver = 'qminos')
    coralme.core.optimization.compute_solution_error(builder1.me_model, keys = { builder1.me_model.mu.magnitude : builder1.me_model.solution.objective_value })

def test_construct_lp_problem_me_model(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    coralme.core.optimization.construct_lp_problem(builder1.me_model, lambdify = False, per_position = False, as_dict = False, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model, lambdify = True, per_position = False, as_dict = False, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model, lambdify = True, per_position = True, as_dict = False, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model, lambdify = True, per_position = True, as_dict = True, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model, lambdify = True, per_position = True, as_dict = True, statistics = True)

def test_construct_lp_problem_m_model(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    coralme.core.optimization.construct_lp_problem(builder1.me_model.gem, lambdify = False, per_position = False, as_dict = False, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model.gem, lambdify = True, per_position = False, as_dict = False, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model.gem, lambdify = True, per_position = True, as_dict = False, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model.gem, lambdify = True, per_position = True, as_dict = True, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model.gem, lambdify = True, per_position = True, as_dict = True, statistics = True)

def test_construct_lp_problem_m_model_from_cobra(shared_builder_troubleshooted):
    builder1 = pytest.shared_builder_troubleshooted
    builder1.me_model.gem = coralme.core.model.MEModel.from_cobra(builder1.me_model.gem)
    coralme.core.optimization.construct_lp_problem(builder1.me_model.gem, lambdify = False, per_position = False, as_dict = False, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model.gem, lambdify = True, per_position = False, as_dict = False, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model.gem, lambdify = True, per_position = True, as_dict = False, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model.gem, lambdify = True, per_position = True, as_dict = True, statistics = False)
    coralme.core.optimization.construct_lp_problem(builder1.me_model.gem, lambdify = True, per_position = True, as_dict = True, statistics = True)
