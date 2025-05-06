"""
| solver↴ model type→ | coralme ME-model | coralme M-model | COBRApy M-model |
|---------------------|------------------|-----------------|-----------------|
| qminos              | Correct          | Correct         | Correct         |
| gurobi              | Correct          | Correct         | Correct         |
| cplex               | Correct          | Correct         | Correct         |
|---------------------|------------------|-----------------|-----------------|
"""

import scipy
import numpy
import pint
import sympy
import sys
import pandas

import cobra
import coralme

# simulation helpers and other functions
def _check_options(model = None, keys = dict(), tolerance = 1e-6, precision = 'quad'):
	# check options
	tolerance = tolerance if tolerance >= 1e-15 else 1e-6
	precision = precision if precision in [ 'quad', 'double', 'dq', 'dqq' ] else 'quad'

	if isinstance(model, coralme.core.model.MEModel) and not model.notes.get('from cobra', False):
		if len(keys.items()) == 0.:
			keys = { model.mu.magnitude : 0.01 }

		for key in list(keys.keys()):
			if isinstance(key, pint.Quantity):
				keys[key.magnitude] = keys.pop(key)
			elif isinstance(key, sympy.Symbol):
				pass
			else:
				keys[sympy.Symbol(key, positive = True)] = keys.pop(key)
	elif isinstance(model, coralme.core.model.MEModel) and model.notes['from cobra'] or isinstance(model, cobra.core.model.Model):
		keys = {} # M-models do not require setting up a key dictionary (mu symbols and their values)

	return keys, tolerance, precision

def _get_evaluated_nlp(model = None, keys = dict(), **kwargs):
	keys, _, _ = _check_options(model, keys = keys) # tolerance and precision not used here

	# populate with stoichiometry with replacement of mu's (Sf contains Se)
	# for single evaluations of the LP problem, direct replacement is faster than lambdify
	Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = kwargs['lp'] if 'lp' in kwargs else construct_lp_problem(model, lambdify = False)

	if lambdas is None:
		Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, Se, lb, ub, keys, atoms)
	else:
		Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, lambdas, lb, ub, keys, atoms)

	return Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm

def compute_solution_error(model, keys = dict()):
	errors = {}

	if not hasattr(model, 'solution'):
		model.optimize()

	lp = construct_lp_problem(model, as_dict = False, lambdify = False)
	Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = _get_evaluated_nlp(keys = keys, **{ 'lp' : lp })

	rows, cols, values = zip(*[ (i, j, v) for (i, j), v in Sf.items() ])
	Sf = scipy.sparse.coo_matrix((values, (rows, cols)), shape = (max(rows) + 1, max(cols) + 1))
	x = numpy.array(list(model.solution.fluxes.values()))

	err = abs(Sf * x)
	errors['max_error'] = err.max()
	errors['sum_error'] = err.sum()
	ub_err = min(ub - x)
	errors['upper_bound_error'] = abs(ub_err) if ub_err < 0 else 0
	lb_err = min(x - lb)
	errors['lower_bound_error'] = abs(lb_err) if lb_err < 0 else 0
	return errors

def construct_lp_problem(model, lambdify = False, per_position = False, as_dict = False, statistics = False) -> tuple:
	"""
	lambdify
		Returns lambda functions for each symbolic stoichiometric coefficient

	per_position
		Returns a list of lambda functions instead of a single 'vectorized' lambda function.
		The lambdify and evaluation is slower, but it allows direct manipulation of the LP.

	as_dict
		Returns a dictionary with keys matching the input of coralme.solver.solver.ME_NLP

	Output:
		A tuple of 11 elements or a dictionary with 11 keys:
			Dictionary with numeric stoichiometric coefficients: { (met, rxn) : float }
			Dictionary with symbolic stoichiometric coefficients: { (met, rxn) : symbol }
			List of lower bounds (numeric and symbolic)
			List of upper bounds (numeric and symbolic)
			List of metabolic bounds (see metabolites._bound property)
			List of objectives (see reaction.objective_coefficient property)
			List of constraint senses (always 'E')
			Set of atoms (i.e., free symbols in symbolic stoichiometric coefficients)
			Dictionary of lambda functions for each symbolic stoichiometry coefficient
			List of reaction IDs (useful when only the LP exists)
			List of metabolites IDs (useful when only the LP exists)
	"""

	# populate empty dictionaries with stoichiometry
	Sf = dict() # floats
	Se = dict() # expressions
	Lr = [ x.id for x in model.reactions ] # reaction identifiers
	Lm = [ x.id for x in model.metabolites ] # metabolite identifiers

	# check how many variables are in the ME-model
	atoms = [] # SymPyDeprecationWarning

	for idx, rxn in enumerate(model.reactions):
		# metabolites derives from symbolic_stoichiometry, replacing everything except model.mu
		for met, value in rxn.metabolites.items():
			met_index = model.metabolites.index(met)
			if hasattr(value, 'subs'):
				# atoms.add(list(value.free_symbols)[0])
				# atoms.update(list(value.free_symbols))
				# TODO: if two or more ME-models are merged, detect if 'mu' is unique or not
				free_symbols = list(value.free_symbols)[0] # only mu
				if free_symbols not in atoms:
					atoms.append(free_symbols)
				Se[met_index, idx] = value
			else:
				Sf[met_index, idx] = value

	if isinstance(model, coralme.core.model.MEModel) and not model.notes.get('from cobra', False):
		if hasattr(model.mu, 'magnitude'):
			lb, ub = zip(*[ (rxn.lower_bound.magnitude, rxn.upper_bound.magnitude) if rxn.functional() else (0., 0.) for rxn in model.reactions ])
		else:
			lb, ub = zip(*[ (rxn.lower_bound, rxn.upper_bound) if rxn.functional() else (0., 0.) for rxn in model.reactions ])
		# evaluate bounds (e.g., DNA_replication)
		lb, ub = zip(*[ (lb.subs(model.default_parameters) if hasattr(lb, 'subs') else lb, ub.subs(model.default_parameters) if hasattr(ub, 'subs') else ub) for lb, ub in zip(lb, ub) ])
	elif isinstance(model, coralme.core.model.MEModel) and model.notes['from cobra']:
		lb, ub = zip(*[ (rxn.lower_bound.magnitude, rxn.upper_bound.magnitude) if rxn.functional() else (0., 0.) for rxn in model.reactions ])
	else:
		lb, ub = zip(*[rxn.bounds for rxn in model.reactions ]) # COBRApy M-models

	b = [ m._bound for m in model.metabolites ] # accumulation
	c = [ r.objective_coefficient for r in model.reactions ]
	# constraint sense eventually will be in the metabolite object
	cs = [ 'E' for m in model.metabolites ]

	if lambdify and isinstance(model, coralme.core.model.MEModel) and not model.notes.get('from cobra', False):
		# 2-3x faster than lambdas = { k:v for k,v in zip(Se.keys(), fn(list(Se.values()))) }
		kwargs = {"docstring_limit":None} if sys.version_info >= (3,8) else {}
		if per_position:
			fn = numpy.vectorize(lambda x: sympy.lambdify(atoms, x, **kwargs))
			lb = [ x for x in fn(lb) ]
			ub = [ x for x in fn(ub) ]
			lambdas = { k:v for k,v in zip(Se.keys(), fn(list(Se.values()))) }
		else:
			lb = sympy.lambdify(atoms, lb, **kwargs) # 5x faster than [ x for x in fn(lb) ]
			ub = sympy.lambdify(atoms, ub, **kwargs) # 5x faster than [ x for x in fn(lb) ]
			lambdas = (list(Se.keys()), sympy.lambdify(atoms, list(Se.values()),**kwargs))
	else:
		lambdas = None

	if statistics:
		print('Sf has {:d} non-zero coefficients ({:.2%})'.format(len(Sf), len(Sf) / (len(Lm)*len(Lr)) ))
		print('Se has {:d} non-zero coefficients ({:.2%})'.format(len(Se), len(Se) / (len(Lm)*len(Lr)) ))

	#TODO: can't pickle attribute lookup _lambdifygenerated on __main__ failed
	#model.lp_full_symbolic = Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm

	lb = list(lb) if isinstance(lb, tuple) else lb
	ub = list(ub) if isinstance(ub, tuple) else ub
	if as_dict:
		return {
			'Sf' : Sf,
			'Se' : Se,
			'xl' : lb,
			'xu' : ub,
			'b' : b,
			'c' : c,
			'cs' : cs,
			'mu' : atoms,
			'lambdas' : lambdas,
			'Lr' : Lr, # list of reaction IDs
			'Lm' : Lm, # list of metabolite IDs
			}
	else:
		return Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm

def rank(model, mu = 0.001):
	Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = construct_lp_problem(model)
	Sp = scipy.sparse.dok_matrix((len(b), len(c)))

	for idx, idj in Sf.keys():
		Sp[idx, idj] = Sf[idx, idj]

	for idx, idj in Se.keys():
		Sp[idx, idj] = float(Se[idx, idj].subs({ model.mu.magnitude : mu }))

	return numpy.linalg.matrix_rank(Sp.todense())

def _solver_solution_to_cobrapy_solution(model, muopt, xopt, yopt, zopt, stat, solver = 'qminos'):
	if hasattr(model, 'reactions'):
		Lr = [ x.id for x in model.reactions ]
		Lm = [ x.id for x in model.metabolites ]
	else:
		Lr, Lm = model

	if solver in ['qminos', 'gurobi']:
		#f = sum([ rxn.objective_coefficient * xopt[idx] for idx, rxn in enumerate(model.reactions) ])
		#x_primal = xopt[ 0:len(model.reactions) ]   # The remainder are the slacks
		x_dict = { rxn : float(xopt[idx]) for idx, rxn in enumerate(Lr) }
		y_dict = { met : float(yopt[idx]) for idx, met in enumerate(Lm) }
		z_dict = { rxn : float(zopt[idx]) for idx, rxn in enumerate(Lr) }
	elif solver == 'cplex':
		#x_primal =
		x_dict = { rxn: float(xopt[rxn].solution_value) for idx, rxn in enumerate(Lr) }
		y_dict = { met: float(yopt[met].dual_value) for idx, met in enumerate(Lm) }
		z_dict = { rxn: float(zopt[rxn].reduced_cost) for idx, rxn in enumerate(Lr) }
	else:
		raise ValueError('solver output not compatible.')

	#model.me.solution = Solution(f, x_primal, x_dict, y, y_dict, 'qminos', time_elapsed, status)
	return cobra.core.Solution(
		objective_value = muopt,
		status = stat,
		fluxes = x_dict, # x_primal is a numpy.array with only fluxes info
		reduced_costs = z_dict,
		shadow_prices = y_dict,
		)

def _set_gurobi_params(gpModel, precision = 'quad', method = 0, ncpus = 1):
	# gpModel.Params.Threads = ncpus
	gpModel.Params.OutputFlag = 0
	gpModel.Params.Presolve = 0
	if precision == 'quad':
		gpModel.Params.Quad = 1
	gpModel.Params.NumericFocus = 3
	gpModel.Params.FeasibilityTol = 1e-9
	gpModel.Params.IntFeasTol = 1e-9
	gpModel.Params.OptimalityTol = 1e-9
	gpModel.Params.Method = method
	gpModel.Params.BarQCPConvTol = 1e-9
	gpModel.Params.BarConvTol = 1e-10
	gpModel.Params.BarHomogeneous = -1
	gpModel.Params.BarCorrectors = 1
	gpModel.Params.Crossover = 4

def _make_gpModel(Sf, lb, ub, c, Lr, Lm, precision = 'quad', method = 2, ncpus = 1):
	import gurobipy as gp
	from gurobipy import GRB
	gpModel = gp.Model()

	# Set params
	_set_gurobi_params(gpModel, precision = precision, method = method, ncpus = ncpus)

	# Define decision variables
	# x = {}
	# for idx, rxn in enumerate(model.reactions):
	# 	x[idx] = gpModel.addVar(lb = lb[idx], ub = ub[idx], name = rxn.id, vtype = GRB.CONTINUOUS)
	x = gpModel.addVars(range(0, len(Lr)), lb = lb, ub = ub, vtype = GRB.CONTINUOUS) # 2x faster than addVar

	# Set objective function
	# lst = [ x[idx] for idx, rxn in enumerate(model.reactions) if rxn.objective_coefficient != 0 ]
	lst = [ x[idx] for idx, obj in enumerate(c) if obj != 0 ] # 4x faster
	gpModel.setObjective(gp.quicksum(lst), gp.GRB.MAXIMIZE)

	# Add constraints for system of linear equations
	for jdx, met in enumerate(Lm):
		lhs = gp.LinExpr()
		for idx, rxn in enumerate(Lr):
			if (jdx, idx) in Sf: # Sf is a dictionary
				lhs += Sf[(jdx, idx)] * x[idx]
		gpModel.addConstr(lhs == 0)

	return gpModel

def _make_mpModel(Sf, lb, ub, c, Lr, Lm):
	# create a cplex model
	from docplex.mp.model import Model
	mpModel = Model(float_precision = 17, cts_by_name = True)

	# Define decision variables
	x = {}
	for idx, rxn in enumerate(Lr):
		x[idx] = mpModel.continuous_var(lb = lb[idx], ub = ub[idx], name = rxn)

	# Set objective function
	lst = [ x[idx] for idx, obj in enumerate(c) if obj != 0 ]
	mpModel.maximize(mpModel.sum(lst))

	# Add constraints for system of linear equations
	for jdx, met in enumerate(Lm):
		lhs = mpModel.linear_expr()
		for idx, rxn in enumerate(Lr):
			if (jdx, idx) in Sf: # Sf is a dictionary
				lhs += Sf[(jdx, idx)] * x[idx]
		mpModel.add_constraint(lhs == 0, ctname = met)

	return mpModel

# Based on Maxwell Neal's work
def _guess_basis(model, keys = dict(), tolerance = 1e-6, precision = 'quad', method = 2, ncpus = 1, **kwargs):
	keys, tolerance, precision = _check_options(model, keys, tolerance, precision)

	lp = construct_lp_problem(model, as_dict = False, lambdify = False)
	Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = _get_evaluated_nlp(keys = keys, **{ 'lp' : lp })

	# make model
	gpModel = _make_gpModel(Sf, lb, ub, c, Lr, Lm, precision = precision, method = method, ncpus = ncpus)

	# optimize
	gpModel.optimize()

	# get basis
	import gurobipy as gp
	if gpModel.status == gp.GRB.Status.OPTIMAL:
		gbasis = numpy.array(gpModel.vbasis)
		basis_guess = numpy.zeros(len(model.reactions) + len(model.metabolites) + 1)
		basis_guess[0:len(model.reactions)][gbasis == +0] = 3
		basis_guess[0:len(model.reactions)][gbasis == -1] = 0
		basis_guess[0:len(model.reactions)][gbasis == -2] = 1
		basis_guess[0:len(model.reactions)][gbasis == -3] = 2

		basis_guess[-1] = 3 # TODO: add comment
		basis_guess = numpy.int32(basis_guess)

		return basis_guess
	else:
		raise ValueError('Optimization failed. Please choose another value for the growth rate.')

# simulation methods: fva, optimize, feasibility, optimize_windows and feas_windows := { feas_gurobi, feas_cplex }
def fva(model,
	reaction_list, fraction_of_optimum, mu_fixed = None, objective = 'biomass_dilution',
	max_mu = 2.8100561374051836, min_mu = 0., maxIter = 100, lambdify = True,
	tolerance = 1e-6, precision = 'quad', verbose = True):

	"""
	Determine the minimum and maximum flux value for each reaction constrained
	to a fraction of the current growth rate (default = 1.0)

	Parameters
	----------
	reaction_list : list of cobra.Reaction or str, optional
		List of reactions IDs and/or reaction objects
	fraction_of_optimum : float, optional
		Must be <= 1.0. Requires that the objective value is at least the
		fraction times maximum objective value. A value of 0.85 for instance
		means that the objective has to be at least at 85% percent of its
		maximum (default 1.0).
	mu_fixed : float, optional
		Set it to avoid the optimization of a ME-model. The growth rate must
		be feasible. If not, the ME-model will be optimized with the following
		options:

		max_mu : float, optional
			Maximum growth rate for initializing the growth rate binary search (GRBS).
		min_mu : float, optional
			Minimum growth rate for initializing GRBS.
		maxIter : int
			Maximum number of iterations for GRBS.
		lambdify : bool
			If True, returns a dictionary of lambda functions for each symbolic
			stoichiometric coefficient
		tolerance : float
			Tolerance for the convergence of GRBS.
		precision : str, {"quad", "double", "dq", "dqq"}
			Precision (quad or double precision) for the GRBS

	verbose : bool
		If True, allow printing.
	"""

	# max_mu is constrained by the fastest-growing bacterium (14.8 doubling time)
	# https://www.nature.com/articles/s41564-019-0423-8

	# check options
	_, tolerance, precision = _check_options(model = model, tolerance = tolerance, precision = precision)
	fraction_of_optimum = fraction_of_optimum if fraction_of_optimum <= 1.0 and fraction_of_optimum >= 0.0 else 1.0
	if isinstance(reaction_list, str):
		reaction_list = [reaction_list]

	# populate with stoichiometry, no replacement of mu's
	if hasattr(model, 'construct_lp_problem'):
		# check if the ME-model has a solution
		if mu_fixed is not None and not hasattr(model, 'solution'):
			model.optimize(max_mu = max_mu, min_mu = min_mu, maxIter = maxIter, lambdify = lambdify,
				tolerance = tolerance, precision = precision, verbose = verbose)

		# set mu_fixed for replacement in a ME-model.
		mu_fixed = model.solution.fluxes.get(objective, mu_fixed) * fraction_of_optimum

		# get mathematical representation
		Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = construct_lp_problem(model, lambdify = lambdify)
	else:
		# not a ME-model, and objective bounds usually are (0, 1000)
		if model.reactions.has_id(objective):
			model.reactions.get_by_id(objective).lower_bound = mu_fixed * fraction_of_optimum
			model.reactions.get_by_id(objective).upper_bound = mu_fixed
		else:
			raise ValueError('Objective reaction \'{:s}\' not in the M-model.'.format(objective))

		# get mathematical representation
		Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = construct_lp_problem(model)

	if verbose:
		print('Running FVA for {:d} reactions. Maximum growth rate fixed to {:g}'.format(len(reaction_list), mu_fixed))

	me_nlp = coralme.solver.solver.ME_NLP(Sf, Se, b, c, lb, ub, cs, atoms, lambdas)

	# We need only reaction objects
	rxns_fva = []
	for rxn in reaction_list:
		if isinstance(rxn, str) and model.reactions.has_id(rxn):
			rxns_fva.append(model.reactions.get_by_id(rxn))
		else:
			rxns_fva.append(rxn)

	obj_inds0 = [ model.reactions.index(rxn) for rxn in rxns_fva for j in range(0, 2) ]
	obj_coeffs = [ ci for rxn in rxns_fva for ci in (1.0, -1.0) ]

	# varyME is a specialized method for multiple min/maximization problems
	obj_inds0, nVary, obj_vals = me_nlp.varyme(mu_fixed, obj_inds0, obj_coeffs, basis = None, verbosity = verbose)

	# Return result consistent with cobrapy FVA
	fva_result = {
		(model.reactions[obj_inds0[2*i]].id): {
			'maximum':obj_vals[2*i],
			'minimum':obj_vals[2*i+1]
			} for i in range(0, nVary//2) }

	return pandas.DataFrame(fva_result).T

def optimize(model,
	max_mu = 2.8100561374051836, min_mu = 0., maxIter = 100, lambdify = True, basis = None,
	tolerance = 1e-6, precision = 'quad', verbose = True, get_reduced_costs = False, solver = "qminos"):

	"""Solves the NLP problem to obtain reaction fluxes for a ME-model.

	Parameters
	----------
	max_mu : float
		Maximum growth rate for initializing the growth rate binary search (GRBS).
	min_mu : float
		Minimum growth rate for initializing GRBS.
	maxIter : int
		Maximum number of iterations for GRBS.
	lambdify : bool
		If True, returns a dictionary of lambda functions for each symbolic
		stoichiometric coefficient.
	basis : list
		qminos basis to hot start optimization. Use `_guess_basis` to get a basis
	tolerance : float
		Tolerance for the convergence of GRBS.
	precision : str, {"quad", "double", "dq", "dqq"}
		Precision (quad or double precision) for the GRBS
	verbose : bool
		If True, allow printing.
	get_reduced_costs : bool
		If True, re-optimize but changing the objective function to 'biomass_dilution'
		and its bounds. New reduced costs and shadow prices will be returned.
	solver : str, { "qminos", "gurobi", "cplex" }
		Solver to use for optimization
	"""

	if isinstance(model, coralme.core.model.MEModel) and not model.notes.get('from cobra', False):
		if solver != "qminos":
			return optimize_windows(model, max_mu = max_mu, min_mu = min_mu, maxIter = maxIter, lambdify = lambdify, tolerance = tolerance, precision = precision, verbose = verbose, solver = solver)

	# check options
	# max_mu is constrained by the fastest-growing bacterium (14.8 min, doubling time)
	# https://www.nature.com/articles/s41564-019-0423-8
	min_mu = min_mu if min_mu >= 0. else 0.
	max_mu = max_mu if max_mu <= 2.8100561374051836 else 2.8100561374051836

	keys, tolerance, precision = _check_options(model, keys = dict(), tolerance = tolerance, precision = precision)

	assert get_reduced_costs == False or get_reduced_costs == lambdify == True, "get_reduced_costs requires lambdify=True"
	per_position = bool(get_reduced_costs)

	if hasattr(model, 'troubleshooting') and not model.troubleshooting or not hasattr(model, 'troubleshooting'):
		print('The MINOS and quad MINOS solvers are a courtesy of Prof Michael A. Saunders. Please cite Ma, D., Yang, L., Fleming, R. et al. Reliable and efficient solution of genome-scale models of Metabolism and macromolecular Expression. Sci Rep 7, 40863 (2017). https://doi.org/10.1038/srep40863\n')

	# populate with stoichiometry, no replacement of mu's
	if hasattr(model, 'construct_lp_problem') and not model.notes.get('from cobra', False):
		Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = construct_lp_problem(model, lambdify = lambdify, per_position = per_position)
	else:
		Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = construct_lp_problem(model, per_position = per_position)
		me_nlp = coralme.solver.solver.ME_NLP(Sf, Se, b, c, lb, ub, cs, atoms, lambdas)
		xopt, yopt, zopt, stat, basis = me_nlp.solvelp(.1, None, precision, probname = 'lp')

		if stat == 'optimal':
			muopt = float(sum([ x*c for x,c in zip(xopt, c) if c != 0 ]))
			model.solution = _solver_solution_to_cobrapy_solution(model, muopt, xopt, yopt, zopt, stat)
			return True
		else:
			if hasattr(model, 'solution'):
				del model.solution
			return False

	if len(atoms) > 1:
		print('Use `me_model.map_feasibility()` to obtain the boundary of feasible solutions.')
		print('Optimization will proceed replacing all growth keys with the same value.')

	me_nlp = coralme.solver.solver.ME_NLP(Sf, Se, b, c, lb, ub, cs, atoms, lambdas)
	muopt, xopt, yopt, zopt, basis, stat = me_nlp.bisectmu(
			mumax = max_mu,
			mumin = min_mu,
			maxIter = maxIter,
			basis = basis,
			tolerance = tolerance,
			precision = precision,
			verbose = verbose)

	if stat == 'optimal':
		# Adapted from Maxwell Neal, 2024
		if get_reduced_costs:
			rxn_idx =  {rxn.id : idx for idx, rxn in enumerate(model.reactions)}
			# Open biomass dilution bounds
			me_nlp.xl[rxn_idx["biomass_dilution"]] = lambda mu : 0.
			me_nlp.xu[rxn_idx["biomass_dilution"]] = lambda mu : 1000.
			# Set new objective coefficient
			me_nlp.c = [1.0 if r=="biomass_dilution" else 0.0 for r in rxn_idx]
			# Solve at muopt
			_xopt, yopt, zopt, _stat, _basis = me_nlp.solvelp(muf = muopt, basis = basis, precision = precision)

		model.solution = _solver_solution_to_cobrapy_solution(model, muopt, xopt, yopt, zopt, stat)
		model.basis = basis
		return True
	else:
		if hasattr(model, 'solution'):
			del model.solution
		if hasattr(model, 'basis'):
			model.basis = None
		return False

# WARNING: Experimental. We could not compile qminos under WinOS, and qminos has a licence restriction for its source code
def optimize_windows(model,
	max_mu = 2.8100561374051836, min_mu = 0., maxIter = 100, lambdify = True,
	tolerance = 1e-6, precision = 'quad', verbose = True, solver = 'gurobi'):

	"""Solves the NLP problem to obtain reaction fluxes for a ME-model. This
	method is used when setting a solver other than qMINOS. It allows to
	use coralME in other OS than Linux.

	Parameters
	----------
	max_mu : float
		Maximum growth rate for initializing the growth rate binary search (GRBS).
	min_mu : float
		Minimum growth rate for initializing GRBS.
	maxIter : int
		Maximum number of iterations for GRBS.
	lambdify : bool
		If True, returns a dictionary of lambda functions for each symbolic
		stoichiometric coefficient
	tolerance : float
		Tolerance for the convergence of GRBS.
	precision : str, {"quad", "double", "dq", "dqq"}
		Precision (quad or double precision) for the GRBS
	verbose : bool
		If True, allow printing.
	"""

	# check options
	keys, tolerance, precision = _check_options(model, keys = dict(), tolerance = tolerance, precision = precision)
	solver = solver if solver in [ 'gurobi', 'cplex' ] else 'gurobi'

	if solver == 'gurobi':
		model.check_feasibility = feas_gurobi
	elif solver == 'cplex':
		model.check_feasibility = feas_cplex
	else:
		print('The \'solver\' must be \'gurobi\' or \'cplex\'.')

	# populate with stoichiometry with replacement of mu's (Sf contains Se)
	# for multiple evaluations of the LP problem, replacement in lambdify'ed Se is faster overall
	Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = construct_lp_problem(model, lambdify = lambdify, per_position = True, as_dict = False)

	# test max_mu
	model.check_feasibility(model, keys = { model.mu.magnitude:max_mu }, precision = 'quad', **{ 'lp' : [Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm] })
	if hasattr(model, 'solution') and model.solution.status == 'optimal':
		return True
	else:
		for idx in range(1, maxIter + 1):
			# Just a sequence of feasibility checks
			muf = (min_mu + max_mu) / 2.
			model.check_feasibility(model, keys = { model.mu.magnitude:muf }, precision = 'quad', **{ 'lp' : [Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm] })

			if hasattr(model, 'solution') and model.solution.status == 'optimal':
				stat_new = 'optimal'
				min_mu = muf
			else:
				stat_new = 1
				max_mu = muf

			if verbose:
				print('{:s}\t{:.16f}\t{:s}'.format(str(idx).rjust(9), muf, 'Not feasible' if stat_new == 1 else stat_new.capitalize()))

			if abs(max_mu - min_mu) <= tolerance and stat_new == 'optimal':
				return True

			if max_mu <= tolerance:
				return False

def feas_windows(model, solver = 'gurobi'):
	if solver == 'gurobi':
		return feas_gurobi
	elif solver == 'cplex':
		return feas_cplex
	else:
		print('The \'solver\' must be \'gurobi\' or \'cplex\'.')
		return None

# WARNING: Experimental. We could not compile qminos under WinOS, and qminos has a licence restriction for its source code
def feas_cplex(model, keys = dict(), **kwargs):
	keys, tolerance, precision = _check_options(model, keys = keys, tolerance = 1e-6, precision = 'double')
	Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = _get_evaluated_nlp(model, keys = keys, **kwargs)

	# make cplex model and optimize
	mpModel = _make_mpModel(Sf, lb, ub, c, Lr, Lm)
	mpModel.solve()

	# output solution
	if mpModel.solve_details.status == 'optimal':
		# WARNING: the objective value is not the objective function flux, but rather the biomass_dilution flux
		muopt = mpModel._vars_by_name['biomass_dilution'].solution_value
		model.solution = _solver_solution_to_cobrapy_solution(model, muopt, mpModel._vars_by_name, mpModel._cts_by_name, mpModel._vars_by_name, stat = 'optimal', solver = 'cplex')
		return True
	else:
		if hasattr(model, 'solution'):
			del model.solution
		return False

# WARNING: Experimental. We could not compile qminos under WinOS, and qminos has a licence restriction for its source code
def feas_gurobi(model, keys = dict(), precision = 'quad', **kwargs):
	keys, tolerance, precision = _check_options(model, keys = keys, tolerance = 1e-6, precision = precision)
	Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = _get_evaluated_nlp(model, keys = keys, **kwargs)

	# make gurobi model and optimize
	gpModel = _make_gpModel(Sf, lb, ub, c, Lr, Lm, precision = precision, method = 2, ncpus = 1)
	gpModel.optimize()

	# output solution
	import gurobipy as gp
	if gpModel.status == gp.GRB.OPTIMAL:
		# WARNING: the objective value is not the objective function flux, but rather the biomass_dilution flux
		muopt = gpModel.x[0]
		model.solution = _solver_solution_to_cobrapy_solution(model, muopt, gpModel.x, gpModel.pi, gpModel.RC, stat = 'optimal', solver = 'gurobi')
		return True
	else:
		if hasattr(model, 'solution'):
			del model.solution
		return False

def feasibility(model, keys = dict(), tolerance = 1e-6, precision = 'quad', basis = None, **kwargs):
	keys, tolerance, precision = _check_options(model = model, keys = keys, tolerance = tolerance, precision = precision)
	Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = _get_evaluated_nlp(model, keys = keys, **kwargs)

	#me_nlp = ME_NLP(me)
	me_nlp = coralme.solver.solver.ME_NLP(Sf, dict(), b, c, lb, ub, cs, set(keys.keys()), None)
	muopt, xopt, yopt, zopt, basis, stat = me_nlp.bisectmu(
			mumax = 1., # mu was already replaced and maxIter is one, so a value here doesn't matter
			mumin = 0.,
			maxIter = 1,
			basis = basis,
			tolerance = tolerance,
			precision = precision,
			verbose = False)

	if stat == 'optimal':
		if len(model.reactions) > 1 and len(model.metabolites) > 1:
			model.solution = _solver_solution_to_cobrapy_solution(model, list(keys.values())[0], xopt, yopt, zopt, stat)
		else:
			x_primal = xopt[ 0:len(Lr) ]   # The remainder are the slacks
			x_dict = { rxn : xopt[idx] for idx, rxn in enumerate(Lr) }
			y_dict = { met : yopt[idx] for idx, met in enumerate(Lm) }
			z_dict = { rxn : zopt[idx] for idx, rxn in enumerate(Lr) }
			model.solution = cobra.core.Solution(
				objective_value = muopt,
				status = stat,
				fluxes = x_dict, # x_primal is a numpy.array with only fluxes info
				reduced_costs = z_dict,
				shadow_prices = y_dict,
				)

		model.basis = basis
		return True
	else:
		if hasattr(model, 'solution'):
			del model.solution
		if hasattr(model, 'basis'):
			model.basis = None
		return False

def map_feasibility(model, keys = { sympy.Symbol('mu', positive = True) : 1. }, tolerance = 1e-6, precision = 'quad'):
	return NotImplemented
