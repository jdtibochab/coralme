import copy
import coralme

# Written originally by Rodrigo Santibanez for ME model
def perform_gene_knockouts(model, genes):
	if isinstance(genes, str):
		genes = set([genes])

	if isinstance(model, coralme.core.model.MEModel):
		test = model.copy()
		for gene in genes:
			gene = 'RNA_{:s}'.format(gene) if not gene.startswith('RNA_') else gene # only valid for ME-models
			if model.metabolites.has_id(gene):
				for TU in test.transcription_data:
					data = test.transcription_data.get_by_id(TU.id)
					test.transcription_data.get_by_id(TU.id).RNA_products = data.RNA_products.difference([gene])
				for rxn in test.reactions.query('transcription_'):
					rxn.update()
			else:
				raise AttributeError('Gene ID \'{:s}\' not in the model.'.format(gene))
	else:
		test = copy.deepcopy(model)
		for gene in genes:
			gene = gene.replace('RNA_', '') if gene.startswith('RNA_') else gene # only valid for M-models
			# similar to cobra.manipulation.delete.remove_genes(model, genes, remove_reactions = True)
			test.genes.get_by_id(gene).knock_out()
	return test

def check_knockout_using_qminos(model, genes, optTol = 1e-15, feasTol = 1e-15):
	# test = copy.deepcopy(model)
	test = perform_gene_knockouts(model, genes)

	nlp = coralme.core.model.MEModel.construct_lp_problem(test, lambdify = True, as_dict = True)
	solver = coralme.solver.solver.ME_NLP(**nlp)
	solver.opt_realdict['lp']['Optimality tol'] = optTol
	solver.opt_realdict['lp']['Feasibility tol'] = feasTol

	if isinstance(test, coralme.core.model.MEModel):
		muopt, xopt, yopt, zopt, basis, stat = solver.bisectmu()
	else:
		xopt, yopt, zopt, stat, basis = solver.solvelp(1., None, 'double')
		muopt = float(sum([ x*c for x,c in zip(xopt, nlp['c']) if c != 0 ]))

	solution = coralme.core.model.MEModel._solver_solution_to_cobrapy_solution(test, muopt, xopt, yopt, zopt, stat)

	#df = test.solution.to_frame()
	# df[(df['fluxes'] < 0) & ~df.index.str.contains('^EX_')]

	return solution

def get_reduced_costs(model, objective_value = 0.1, target_reaction = 'biomass_dilution'):
	if not model.reactions.has_id(target_reaction):
		raise AttributeError('Model has no reaction \'{:s}\''.format(target_reaction))

	nlp = model.construct_lp_problem(as_dict = True, per_position = True)

	# change objective function and its bounds
	rxn_id = { x:idx for idx,x in enumerate(nlp['Lr']) }

	# remove objective function
	nlp['c'] = [0.]*len(nlp['c'])

	# change target reaction
	nlp['xl'][rxn_id[target_reaction]] = 0.
	nlp['xu'][rxn_id[target_reaction]] = 1000.
	nlp['c'][rxn_id[target_reaction]] = 1.

	solver = coralme.solver.solver.ME_NLP(**nlp)
	xopt, yopt, zopt, stat, basis = solver.solvelp(muf = objective_value, basis = None, precision = 'quad')
	muopt = [ x*c for x,c in zip(xopt, nlp['c']) if c != 0 ][0]
	sol = coralme.core.model.MEModel._solver_solution_to_cobrapy_solution(model, muopt, xopt, yopt, zopt, stat)

	return sol.reduced_costs

def revert_gene_knockouts(model, genes):
	raise NotImplementedError

def single_gene_deletion(model, gene, threshold = 0.01, solver = 'qminos'):
	if solver not in [ 'gurobi', 'qminos' ]:
		raise Exception('The solver argument should be \'qminos\' or any valid for COBRApy models, such as \'gurobi\'.')
	if isinstance(gene, (list, set)):
		raise Exception('The method is limited to one gene only. Use model.perform_gene_knockouts(), followed by model.optimize() or model.feasibility().')

	test = perform_gene_knockouts(model, gene)

	if isinstance(model, coralme.core.model.MEModel):
		if test.feasibility({ test.mu : threshold }):
			return gene, False # gene is not essential
		else:
			return gene, True # gene is essential
	else:
		if solver == 'qminos':
			# feasibility developed to work with ME-models only
			# output is True or False; if True, test.solution is created
			coralme.core.model.MEModel.optimize(test)
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

if __name__ == '__main__':
	single_gene_deletion(model, gene, threshold = 0.01, solver = 'qminos')
