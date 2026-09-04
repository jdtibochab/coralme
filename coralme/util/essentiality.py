import copy
import cobra
import coralme
import tqdm

# Written originally by Rodrigo Santibanez for coralME models and COBRApy models
def perform_gene_knockouts(model, genes, mets_to_test = []):
	if isinstance(genes, (str, coralme.core.component.TranscribedGene)):
		genes = set([genes])

	if isinstance(model, coralme.core.model.MEModel) and model.notes.get('from cobra', False) is False:
		test = model.copy()
		for gene in genes:
			gene = gene.id if isinstance(gene, coralme.core.component.TranscribedGene) else gene
			gene = 'RNA_{:s}'.format(gene) if not gene.startswith('RNA_') else gene # only valid for ME-models
			if model.metabolites.has_id(gene):
				for TU in test.transcription_data:
					data = test.transcription_data.get_by_id(TU.id)
					test.transcription_data.get_by_id(TU.id).RNA_products = data.RNA_products.difference([gene])
				for rxn in test.reactions.query('transcription_'):
					rxn.update()
			else:
				raise AttributeError('Gene ID \'{:s}\' is not in the model.'.format(gene))
	else:
		test = copy.deepcopy(model)
		for gene in genes:
			gene = gene.id if isinstance(gene, cobra.core.gene.Gene) else gene
			# similar to cobra.manipulation.delete.remove_genes(model, genes, remove_reactions = True)
			if test.genes.has_id(gene):
				test.genes.get_by_id(gene).knock_out()
			else:
				raise AttributeError('Gene ID \'{:s}\' is not in the model.'.format(gene))

	for met in mets_to_test:
		if test.metabolites.has_id(met) and not test.reactions.has_id('SK_{:s}'.format(met)):
			test.add_boundary(test.metabolites.get_by_id(met), type = 'sink', lb = 0., ub = 0.)
		else:
			raise AttributeError('Metabolite ID \'{:s}\' is not in the model or a sink reaction already exists.'.format(gene))

	return test

def check_knockout_using_qminos(model, genes, optTol = 1e-15, feasTol = 1e-15):
	# test = copy.deepcopy(model)
	test = perform_gene_knockouts(model, genes)

	nlp = coralme.core.optimization.construct_lp_problem(test, lambdify = True, as_dict = True, per_position = False)
	solver = coralme.solver.solver.ME_NLP(**nlp)
	solver.opt_realdict['lp']['Optimality tol'] = optTol
	solver.opt_realdict['lp']['Feasibility tol'] = feasTol

	if isinstance(test, coralme.core.model.MEModel):
		muopt, xopt, yopt, zopt, basis, stat = solver.bisectmu()
	else:
		xopt, yopt, zopt, stat, basis = solver.solvelp(1., None, 'double')
		muopt = float(sum([ x*c for x,c in zip(xopt, nlp['c']) if c != 0 ]))

	solution = coralme.core.optimization._solver_solution_to_cobrapy_solution(test, muopt, xopt, yopt, zopt, stat)
	return solution

def create_ko_model_in_lp_format(model, genes, growth_rate, mets_to_test, *args):
	if len(args) == 4:
		model, genes, growth_rate, mets_to_test = args

	test = perform_gene_knockouts(model, genes, mets_to_test)
	nlp = test.construct_lp_problem(lambdify = False, as_dict = True, per_position = True)
	nlp['Sf'], nlp['Se'], nlp['xl'], nlp['xu'] = coralme.builder.helper_functions.evaluate_lp_problem(nlp['Sf'], nlp['Se'], nlp['xl'], nlp['xu'], { test.mu.magnitude : growth_rate }, nlp['mu'])

	indexes = { met:(test.reactions._dict['SK_{:s}'.format(met)], test.metabolites._dict[met]) for met in mets_to_test }

	return nlp, indexes

def check_many_mets_at_a_time(args):
	for mid, (rxn, met) in args[1].items():
		# print(args[0]['xl'][rxn], args[0]['xu'][rxn])
		args[0]['xl'][rxn] = -1000.
		args[0]['xu'][rxn] = +1000.
		# print(args[0]['xl'][rxn], args[0]['xu'][rxn])

	nlp = args[0]
	xopt, yopt, zopt, stat, basis = coralme.solver.solver.ME_NLP(**nlp).solvelp(muf = None, basis = None, precision = 'quad')
	muopt = [ x*c for x,c in zip(xopt, nlp['c']) if c != 0 ][0]
	sol = coralme.core.optimization._solver_solution_to_cobrapy_solution((nlp['Lr'], nlp['Lm']), muopt, xopt, yopt, zopt, stat)
	return sol

def check_all_mets_at_a_time(nlp, indexes):
	return check_many_mets_at_a_time((nlp, indexes))

def get_reduced_costs_from_nlp(nlp, objective_value = 0.1):
	xopt, yopt, zopt, stat, basis = coralme.solver.solver.ME_NLP(**nlp).solvelp(muf = objective_value, basis = None, precision = 'quad')
	muopt = [ x*c for x,c in zip(xopt, nlp['c']) if c != 0 ][0]
	sol = coralme.core.optimization._solver_solution_to_cobrapy_solution((nlp['Lr'], nlp['Lm']), muopt, xopt, yopt, zopt, stat)
	return sol.reduced_costs

def get_reduced_costs_from_model(model, objective_value = 0.1, target_reaction = 'biomass_dilution'):
	if not model.reactions.has_id(target_reaction):
		raise AttributeError('Model has no reaction \'{:s}\''.format(target_reaction))

	nlp = model.construct_lp_problem(lambdify = False, as_dict = True, per_position = True)

	# change objective function and its bounds
	rxn_id = { x:idx for idx,x in enumerate(nlp['Lr']) }

	# remove objective function
	nlp['c'] = [0.]*len(nlp['c'])

	# change target reaction
	nlp['xl'][rxn_id[target_reaction]] = 0.
	nlp['xu'][rxn_id[target_reaction]] = 1000.
	nlp['c'][rxn_id[target_reaction]] = 1.

	return get_reduced_costs_from_nlp(nlp, objective_value)

def revert_gene_knockouts(model, genes):
	raise NotImplementedError

def single_gene_deletion(model, gene, threshold = 0.01, solver = 'qminos'):
	if solver not in [ 'gurobi', 'qminos' ]:
		raise Exception('The solver argument should be \'qminos\' or \'gurobi\'.')
	if isinstance(gene, (list, set)):
		raise Exception('The method is limited to one gene only. Use model.perform_gene_knockouts(), followed by model.optimize() or model.feasibility().')

	test = perform_gene_knockouts(model, gene)

	if isinstance(model, coralme.core.model.MEModel):
		if test.feasibility({ test.mu.magnitude : threshold }):
			return gene, False # gene is not essential
		else:
			return gene, True # gene is essential
	else:
		if solver == 'qminos':
			# feasibility developed to work with ME-models only
			# output is True or False; if True, test.solution is created
			coralme.core.optimization.optimize(test)
		elif solver in ['gurobi']:
			test.solver = solver
			test.solution = test.optimize()
		else:
			test.solver = 'gurobi'

	# if sol.status != 'optimal' or sol.objective_value < threshold:
	if hasattr(test, 'solution'):
		if test.solution.status == 'infeasible':
			return gene, True # gene is essential
		elif test.solution.status == 'optimal' and test.solution.objective_value < threshold:
			return gene, True # gene is essential
		else:
			return gene, False # gene is not essential (over the threshold)
	else:
		return gene, True # gene is essential

def test_cofactor_essentiality(model, cofactors, threshold = 0.01):
	if isinstance(cofactors, str):
		if not model.metabolites.has_id(cofactors):
			raise AttributeError('Cofactor ID \'{:s}\' is not in the model.'.format(cofactors))
		else:
			cofactors = model.metabolites.get_by_id(cofactors)
	if isinstance(cofactors, coralme.core.component.Metabolite):
		cofactors = set([ cofactors ])

	test = copy.deepcopy(model)
	# test._mu = sympy.Symbol('mu', positive = True)
	# test._mu_old = test.mu

	# get proteins associated to cofactors
	proteins = set()
	for cofactor in cofactors:
		for rxn in [ x for x in cofactor.reactions if x.id.startswith('formation_') ]:
			for reactant in rxn.reactants:
				if isinstance(reactant, (coralme.core.component.ProcessedProtein, coralme.core.component.TranslatedGene)):
					proteins.add(reactant.id)

	if len(proteins) == 0:
		return None,None  

	# delete multiple formation reactions at once
	rxns = set()
	for protein in proteins:
		tmp = [ x.id for x in test.metabolites.get_by_id(protein).reactions 
		 	if x.id.startswith('formation') and 
			cofactor.id.removesuffix('_' + cofactor.compartment) in x.id 
			]
		for rxn in tmp:
			test.reactions.get_by_id(rxn).bounds = (0., 0.)
			rxns.add(rxn)
			
	if test.feasibility({ model.mu.magnitude : threshold }):
		return rxns,'Can grow without cofactor'
	else:
		return rxns,'Can\'t grow without cofactor'

def single_cofactor_essentiality_analysis(model, threshold = 0.01):
	if isinstance(model, coralme.core.model.MEModel) and model.notes.get('from cobra', False) is True:
		raise Exception('The model must be a coralME ME-model.')

	results = {}
	cofactors = model.get_cofactors
	for cofactor in tqdm.tqdm(cofactors):
		rxns, status = test_cofactor_essentiality(model, [cofactor], threshold = threshold)
		results[cofactor] = {'rxns': rxns, 'status': status}
	return results

# Originally developed by Diego Tec-Campos, UCSD, 2026
# Modified from COBRApy's single_gene_deletion functions for compatibility with coralME M-models
import pandas as pd
from typing import Iterable, List, Optional, Sequence, Tuple, Union
CofactorLike = Union[str, object]
def _reaction_cofactor_ids(reaction) -> set[str]:
    """Return cofactor identifiers occurring in a reaction cofactor rule."""
    rule = getattr(reaction, "cofactors", None)
    if rule is None:
        return set()
    return {str(cofactor) for cofactor in rule.genes}

def _model_cofactor_ids(model) -> List[str]:
    """Return sorted cofactor identifiers associated with model reactions."""
    cofactors = set()

    # Prefer a native CORALME cofactor collection if present.
    model_cofactors = getattr(model, "get_cofactors", None)
    if model_cofactors is not None:
        try:
            for cofactor in model_cofactors:
                cofactor_id = getattr(cofactor, "id", cofactor)
                cofactors.add(str(cofactor_id))
        except TypeError:
            pass

    # reaction.cofactors is the canonical fallback and is available for
    # models loaded through FromExcel -> MEModel.from_cobra.
    for reaction in model.reactions:
        cofactors.update(_reaction_cofactor_ids(reaction))

    return sorted(cofactors)

def _entity_ids(entities: Iterable[CofactorLike]) -> List[str]:
    """Return identifiers from strings or objects exposing an ``id`` field."""
    result = []

    for entity in entities:
        identifier = getattr(entity, "id", entity)
        result.append(str(identifier))

    return result

def _normalize_cofactor_list(
    model,
    cofactor_list: Optional[Iterable[CofactorLike]],
) -> List[str]:
    """Normalize and validate a requested cofactor list."""
    available = _model_cofactor_ids(model)
    available_set = set(available)

    if cofactor_list is None:
        return available

    if isinstance(cofactor_list, str) or hasattr(cofactor_list, "id"):
        cofactor_list = [cofactor_list]

    requested = _entity_ids(cofactor_list)

    # Remove duplicates while retaining user-supplied order.
    requested = list(dict.fromkeys(requested))

    missing = [cofactor for cofactor in requested if cofactor not in available_set]
    if missing:
        raise KeyError(
            "Cofactor(s) not associated with model reactions: "
            + ", ".join(missing)
        )

    return requested

def _get_coralme_growth(
    model: coralme.core.model.MEModel,
    **kwargs,
) -> Tuple[float, str]:
    """Return objective value and status for a CORALME model."""
    options = dict(kwargs)
    options.setdefault("verbose", False)

    # success = model.optimize(**options)
    success = model.feasibility(**options) # faster, we need a yes, no answer at a single growth rate

    if (
        success
        and getattr(model, "solution", None) is not None
        and model.solution.objective_value is not None
    ):
        return (
            float(model.solution.objective_value),
            str(model.solution.status),
        )

    return float("nan"), "not_optimal"

def _cofactor_deletion(
    model,
    cofactor_ids: Sequence[str],
    method: str = "fba",
    solution=None,
    **kwargs,
) -> Tuple[List[str], float, str]:
    """Perform one cofactor deletion simulation."""
    test_model, _ = perform_cofactor_knockouts(model, cofactor_ids)

    if isinstance(test_model, coralme.core.model.MEModel):
        if method != "fba":
            raise NotImplementedError(
                "MOMA and ROOM are not currently implemented for "
                "coralme.core.model.MEModel cofactor deletions."
            )

        growth, status = _get_coralme_growth(test_model, **kwargs)
        return list(cofactor_ids), growth, status

    if isinstance(test_model, cobra.Model):
        if "moma" in method:
            add_moma(
                test_model,
                solution=solution,
                linear="linear" in method,
            )
        elif "room" in method:
            add_room(
                test_model,
                solution=solution,
                linear="linear" in method,
                **kwargs,
            )

        growth, status = _get_cobra_growth(test_model)
        return list(cofactor_ids), growth, status

    raise TypeError(
        "model must be a cobra.Model or coralme.core.model.MEModel."
    )

def single_cofactor_deletion(
    model,
    cofactor_list: Optional[List[CofactorLike]] = None,
    method: str = "fba",
    solution=None,
    processes: Optional[int] = None,
    **kwargs,
) -> pd.DataFrame:
    """Knock out each cofactor from ``cofactor_list``.

    Parameters
    ----------
    model
        Cofactor-aware ``cobra.Model`` or ``coralme.core.model.MEModel``.
        Reaction cofactor requirements must be stored in
        ``reaction.cofactors`` as ``cobra.core.GPR`` objects.
    cofactor_list : list of str or cofactor objects, optional
        Cofactors to delete individually. If not passed, all cofactors
        associated with reaction cofactor rules are used.
    method : {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict growth (default ``"fba"``). For CORALME
        ``MEModel`` objects, the current implementation supports FBA only.
        COBRApy models additionally support MOMA and ROOM.
    solution : cobra.Solution, optional
        Previous solution used as reference for (linear) MOMA or ROOM.
        Ignored for FBA.
    processes : int, optional
        Included for API compatibility with COBRApy deletion functions.
        Parallel execution is not yet enabled for CORALME cofactor deletion.
        ``None`` and ``1`` execute sequentially.
    **kwargs
        For CORALME models, keyword arguments are forwarded to
        ``MEModel.optimize``. For COBRApy ROOM simulations, keyword arguments
        are forwarded to ``add_room``.

    Returns
    -------
    pandas.DataFrame
        One row per single-cofactor deletion with columns:

        ``ids``
            Set containing the deleted cofactor identifier.
        ``growth``
            Predicted objective value after deletion.
        ``status``
            Optimization status.

    Notes
    -----
    Only cofactors that occur in at least one reaction cofactor rule are
    included by default. Cofactors with zero reaction associations are not
    considered part of the functional deletion set.
    """
    method = str(method).lower().strip()

    valid_methods = {
        "fba",
        "moma",
        "linear moma",
        "room",
        "linear room",
    }

    if method not in valid_methods:
        raise ValueError(
            f"Unknown method '{method}'. "
            f"Expected one of {sorted(valid_methods)}."
        )

    if isinstance(model, coralme.core.model.MEModel) and method != "fba":
        raise NotImplementedError(
            "CORALME MEModel cofactor deletion currently supports FBA only."
        )

    if processes not in (None, 1):
        raise NotImplementedError(
            "Parallel cofactor deletion is not yet implemented. "
            "Use processes=None or processes=1."
        )

    cofactor_ids = _normalize_cofactor_list(model, cofactor_list)

    results = []

    for cofactor_id in cofactor_ids:
        ids, growth, status = _cofactor_deletion(
            model,
            [cofactor_id],
            method=method,
            solution=solution,
            **kwargs,
        )

        results.append(
            (
                set(ids),
                growth,
                status,
            )
        )

    return pd.DataFrame(
        results,
        columns=["ids", "growth", "status"],
    )

def _candidate_reactions(model, cofactor_ids: set[str]):
    """Return reactions containing at least one selected cofactor."""
    candidates = set()

    # Use a native cofactor.reactions interface when available. This keeps
    # the implementation compatible with CORALME's emerging cofactor API.
    model_cofactors = getattr(model, "get_cofactors", None)

    if model_cofactors is not None and hasattr(model_cofactors, "get_by_id"):
        for cofactor_id in cofactor_ids:
            try:
                cofactor = model_cofactors.get_by_id(cofactor_id)
            except (KeyError, ValueError):
                continue

            reactions = getattr(cofactor, "reactions", None)
            if reactions is not None:
                candidates.update(reactions)

    # GPR-based fallback used by the current public CORALME model loaded via
    # FromExcel -> MEModel.from_cobra.
    if not candidates:
        for reaction in model.reactions:
            if _reaction_cofactor_ids(reaction).intersection(cofactor_ids):
                candidates.add(reaction)

    return candidates

def constrained_reactions(
    model,
    cofactor_ids: Union[str, Iterable[str]],
) -> List[str]:
    """Return reactions disabled by deletion of one or more cofactors.

    Parameters
    ----------
    model
        Cofactor-aware COBRA or CORALME model.
    cofactor_ids
        Cofactor identifier or iterable of identifiers considered unavailable.

    Returns
    -------
    list of str
        Sorted reaction identifiers whose complete Boolean cofactor rule
        evaluates to ``False``.
    """
    if isinstance(cofactor_ids, str):
        knockout_ids = {cofactor_ids}
    else:
        knockout_ids = {str(cofactor) for cofactor in cofactor_ids}

    constrained = []

    for reaction in _candidate_reactions(model, knockout_ids):
        rule = getattr(reaction, "cofactors", None)

        if rule is None or not rule.to_string().strip():
            continue

        if not rule.eval(knockout_ids):
            constrained.append(reaction.id)

    return sorted(constrained)

def perform_cofactor_knockouts(
    model,
    cofactors: Union[CofactorLike, Iterable[CofactorLike]],
):
    """Return a copy of ``model`` after deleting selected cofactors.

    Reactions are constrained to zero only when their complete Boolean
    cofactor rule becomes false.

    Parameters
    ----------
    model
        Cofactor-aware COBRA or CORALME model.
    cofactors
        Cofactor identifier/object or iterable of identifiers/objects.

    Returns
    -------
    knockout_model
        Copy of the input model with affected reactions constrained to zero.
    constrained_rxns : list of str
        Reaction identifiers constrained by the deletion.
    """
    if isinstance(cofactors, str) or hasattr(cofactors, "id"):
        cofactors = [cofactors]

    cofactor_ids = _entity_ids(cofactors)

    knockout_model = model.copy()
    constrained_rxns = constrained_reactions(knockout_model, cofactor_ids)

    for reaction_id in constrained_rxns:
        knockout_model.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)

    return knockout_model, constrained_rxns