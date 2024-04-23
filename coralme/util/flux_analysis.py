#!/usr/bin/python3
import cobra
import coralme
import pandas
import sympy
import tqdm


def exchange_single_model(me, flux_dict = 0, solution=0):
	"""
	Returns a summary of exchange reactions and fluxes
	"""
	complete_dict = {'id':[],'name':[],'reaction':[],'lower_bound':[],'upper_bound':[],'flux':[]}

	if solution:
		flux_dict = solution.fluxes
	elif not flux_dict:
		flux_dict = me.solution.fluxes

	for rxn in me.reactions:
		try:
			if rxn.reactants and rxn.products:
				continue
		except:
			continue
		flux = flux_dict[rxn.id]

		if not flux:
			continue
		rxn_name = rxn.name
		reaction = rxn.reaction
		lb = rxn.lower_bound
		ub = rxn.upper_bound

		complete_dict['id'].append(rxn.id)
		complete_dict['name'].append(rxn_name)
		complete_dict['reaction'].append(reaction)
		complete_dict['lower_bound'].append(lb)
		complete_dict['upper_bound'].append(ub)
		complete_dict['flux'].append(flux)


	df = pandas.DataFrame(complete_dict).set_index('id')
	return df

def get_met_coeff(stoich,growth_rate,growth_key='mu'):
	"""
	Returns a float stoichiometric coefficient of a metabolite
	in a reaction. If the coefficient is a sympy expression,
	it substitutes a growth rate to get a float.
	"""
	if hasattr(stoich, 'subs'):
		try:
			return float(stoich.subs(growth_key,growth_rate))
		except:
			return None
	return stoich

def summarize_reactions(model,met_id,only_types=(),ignore_types = ()):
	"""
	Returns a summary of reactions and their fluxes in a model
	"""
	reactions = get_reactions_of_met(model,met_id,only_types=only_types,
								 ignore_types=ignore_types,verbose=False)
	d = {}
	for r in reactions:
		if r.bounds == (0,0):
			continue
		d[r.id] = {
			'name':r.name,
			'gene_reaction_rule':r.gene_reaction_rule,
			'reaction':r.reaction,
			'notes':r.notes if r.notes else ''
		}
	df = pandas.DataFrame.from_dict(d).T
	return df[['name','gene_reaction_rule','reaction','notes']] if not df.empty else 'No reaction found'

def flux_based_reactions(model,
						 met_id,
						 only_types=(),
						 ignore_types = (),
						 threshold = 0.,
						 flux_dict=0,
						 solution = None,
						 keffs=False,
						 verbose=False):
	"""
	Returns a summary of the mass balance of a metabolite in a
	flux distribution.
	"""

	if not flux_dict:
		#flux_dict = model.solution.x_dict
		if not hasattr(model,'solution') or not model.solution:
			if solution is not None:
				flux_dict = solution.fluxes
			else:
				print('No solution in model object')
				flux_dict = {r.id:0. for r in model.reactions}
		else:
			flux_dict = model.solution.fluxes
	mu = model.mu if hasattr(model,'mu') else ''
	reactions = get_reactions_of_met(model,met_id,only_types=only_types,
									 ignore_types=ignore_types,verbose=False,growth_key=mu)
	if len(reactions) == 0:
		print('No reactions found for {}'.format(met_id))
		return

	met = model.metabolites.get_by_id(met_id)
	result_dict = {}
	g = flux_dict.get('biomass_dilution',None)
	for rxn in (tqdm.tqdm(reactions) if verbose else reactions):
		f = flux_dict[rxn.id]
		result_dict[rxn.id] = {}
		if f:
			coeff = get_met_coeff(rxn.metabolites[met],
								  g,
								growth_key=model.mu if hasattr(model,"mu") else None)

		else:
			coeff = 0
		if coeff is None:
			print('Could not convert expression to float in {}'.format(rxn.id))
			continue
		try:
			result_dict[rxn.id]['lb'] = rxn.lower_bound if isinstance(rxn.lower_bound, sympy.Symbol) else float(rxn.lower_bound)
			result_dict[rxn.id]['ub'] = rxn.upper_bound if isinstance(rxn.upper_bound, sympy.Symbol) else float(rxn.upper_bound)
		except:
			print('Could not convert bounds to float in {}'.format(rxn.id))

		result_dict[rxn.id]['rxn_flux'] = f
		result_dict[rxn.id]['met_flux'] = f*coeff
		result_dict[rxn.id]['reaction'] = rxn.reaction
		if keffs:
			result_dict[rxn.id]['keff'] = rxn.keff if hasattr(rxn,'keff') else ''
	df = pandas.DataFrame.from_dict(result_dict).T

	df['rxn_flux'] = df['rxn_flux'].astype(float)
	df['met_flux'] = df['met_flux'].astype(float)

	df = df.loc[df['met_flux'].abs().sort_values(ascending=False).index]
	return df#[df['ub'] != 0]

def get_reactions_of_met(me,met,s = 0, ignore_types = (),only_types = (), verbose = False,growth_key='mu'):
	"""
	Returns the reactions of a metabolite. If directionality is not set (s=0),
	the behavior is analogous to met.reactions. However, setting s=1 or s=-1,
	returns the reactions that produce or consume it, respectively.
	"""

	met_stoich = 0
	if only_types:
		only_reaction_types = tuple([getattr(coralme.core.reaction,i) for i in only_types])
	elif ignore_types:
		ignore_reaction_types = tuple([getattr(coralme.core.reaction,i) for i in ignore_types])
	reactions = []

	if not hasattr(me.metabolites,met):
		return reactions
	for rxn in me.metabolites.get_by_id(met).reactions:
		if only_types and not isinstance(rxn, only_reaction_types):
			continue
		elif ignore_types and isinstance(rxn, ignore_reaction_types):
			continue
		try:
			met_obj = me.metabolites.get_by_id(met)
			pos = 1 if get_met_coeff(rxn.metabolites[met_obj],0.1,growth_key=growth_key) > 0 else -1
			rev = 1 if rxn.lower_bound < 0 else 0
			fwd = 1 if rxn.upper_bound > 0 else 0
		except:
			if verbose:
				print(rxn.id, ' could not parse')
			else:
				pass
		try:
			if not s:
				reactions.append(rxn)
				if verbose:
					print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)

			elif s == pos*fwd or s == -pos*rev:
				reactions.append(rxn)
				if verbose:
					print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)

		except:
			if verbose:
				print(rxn.id, 'no reaction')
			else:
				pass
	return reactions
