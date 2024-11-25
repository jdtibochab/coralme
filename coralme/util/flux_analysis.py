#!/usr/bin/python3
import cobra
import coralme
import pandas
import sympy
import tqdm
import re


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
						 include_zeroes=True,
						 solution = None,
						 keffs=False,
						 verbose=False):
	"""
	Returns a summary of the mass balance of a metabolite in a
	flux distribution.
	"""
	if flux_dict:
		pass
	elif solution:
		flux_dict = solution.fluxes
	elif hasattr(model,'solution') and model.solution:
		flux_dict = model.solution.fluxes
	else:
		print('No solution in model object')
		flux_dict = {r.id:0. for r in model.reactions}
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
	if include_zeroes:
		return df#[df['ub'] != 0]
	else:
		return df[df['rxn_flux'] != 0.]

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

from coralme.builder.helper_functions import substitute_value,get_next_from_type
def get_immediate_partitioning(p):
	"""
	This function calculates the partitioning of a metabolite
	to its different immediate products across reactions.
	"""
	tmp = flux_based_reactions(p._model,p.id)["met_flux"]
	tmp = tmp[tmp<0]
	dct = tmp.div(tmp.sum()).to_dict()
	return {p._model.get(k):v for k,v in dct.items()}

def get_partitioning(m, seen = set(),final_fraction=1.0):
	"""
	This is a modified function from find_complexes, which keeps
	track of the partitioning of proteins according to flux
	distributions contained in a model object.
	"""
	if not m:
		return set()
	if m in seen:
		return set()
	if final_fraction == 0:
		return set()

	seen.add(m)

	# Reaction objects
	if isinstance(m,coralme.core.reaction.PostTranslationReaction):
		return get_partitioning(get_next_from_type(m.metabolites,coralme.core.component.ProcessedProtein), seen=seen,final_fraction=final_fraction)
	if isinstance(m,coralme.core.reaction.ComplexFormation):
		return get_partitioning(get_next_from_type(m.metabolites,coralme.core.component.Complex), seen=seen,final_fraction=final_fraction)
	if isinstance(m,coralme.core.reaction.GenericFormationReaction):
		return get_partitioning(get_next_from_type(m.metabolites,coralme.core.component.GenericComponent), seen=seen,final_fraction=final_fraction)
	if isinstance(m,coralme.core.reaction.tRNAChargingReaction):
		return get_partitioning(get_next_from_type(m.metabolites,coralme.core.component.GenerictRNA), seen=seen,final_fraction=final_fraction)
	if isinstance(m,coralme.core.reaction.MetabolicReaction):
		return get_partitioning(get_next_from_type(m.metabolites,coralme.core.component.Complex), seen=seen,final_fraction=final_fraction) | \
				get_partitioning(get_next_from_type(m.metabolites,coralme.core.component.GenericComponent), seen=seen,final_fraction=final_fraction)

	if isinstance(m,coralme.core.reaction.SummaryVariable):
		return set()

	partitioning = get_immediate_partitioning(m)

	# Metabolite objects
	if isinstance(m,coralme.core.component.TranslatedGene):
		cplxs = set()
		for r,fraction in partitioning.items():
			if substitute_value(m,r.metabolites[m] > 0):
				continue
			cplxs = cplxs | get_partitioning(r, seen=seen,final_fraction=final_fraction*fraction)
		return cplxs
	if isinstance(m,coralme.core.component.TranscribedGene):
		translated_protein = m.id.replace('RNA_','protein_')
		if translated_protein in m._model.metabolites:
			return get_partitioning(m._model.metabolites.get_by_id(translated_protein), seen=seen,final_fraction=final_fraction)
		cplxs = set()
		for r,fraction in partitioning.items():
			if substitute_value(m,r.metabolites[m] > 0):
				continue
			cplxs = cplxs | get_partitioning(r, seen=seen,final_fraction=final_fraction*fraction)
		return cplxs
	if isinstance(m,coralme.core.component.ProcessedProtein):
		cplxs = set()
		for r,fraction in partitioning.items():
			if substitute_value(m,r.metabolites[m] > 0):
				continue
			cplxs = cplxs | get_partitioning(r, seen=seen,final_fraction=final_fraction*fraction)
		return cplxs

	if isinstance(m,coralme.core.component.Complex) or isinstance(m,coralme.core.component.GenericComponent) or isinstance(m,coralme.core.component.GenerictRNA):
		other_formations = [(r,fraction) for r,fraction in partitioning.items() if (isinstance(r,coralme.core.reaction.ComplexFormation) or isinstance(r,coralme.core.reaction.GenericFormationReaction)) and substitute_value(m,r.metabolites[m]) < 0]
		cplxs = set([(m,final_fraction)])
		if other_formations:
			cplxs = set()
			for r,fraction in other_formations:
				cplxs = cplxs | get_partitioning(r, seen=seen,final_fraction=final_fraction*fraction)
		# print(3,cplxs)
		return cplxs

	return set()

def get_reduced_costs(nlp,muopt,rxn_idx,basis=None,precision=1e-6):
	# Adapted from Maxwell Neal, 2024
	# Open biomass dilution bounds
	nlp.xl[rxn_idx["biomass_dilution"]] = lambda mu : 0
	nlp.xu[rxn_idx["biomass_dilution"]] = lambda mu : 1000
	# Set new objective coefficient
	nlp.c = [1.0 if r=="biomass_dilution" else 0.0 for r in rxn_idx]
	# Solve at muopt
	_xopt, yopt, zopt, _stat, _basis = nlp.solvelp(muf = muopt, basis = basis, precision = precision)
	return _xopt, yopt, zopt, _stat, _basis

elements = ["C","H","O","N","P","S","Mn","Cr","Ni","Cu","Zn","Sb","Ca","H","Co","K","As","Cd","Mg","Mo","Fe","X","Hg","Pb","Ag","Se","Cl","Na","W","R"]
def get_biomass_formula(me,fluxes):
	if isinstance(fluxes,dict):
		fluxes = pandas.Series(fluxes)
	pattern = "^EX_|^DM_|^SK_|^TS"
	fluxes = fluxes[fluxes.index.str.contains(pattern)]
	fluxes = fluxes[fluxes!=0]
	fluxes = fluxes.drop([r for r in fluxes.index if me.get("lipid_biomass") in me.get(r).products])

	Formulas = {}
	for e in elements:
		Formulas[e] = {m:me.get(re.split(pattern,m)[1]).elements.get(e,0) for m in fluxes.index}
	Formulas = pandas.DataFrame.from_dict(Formulas)
	AbsoluteFormula = -fluxes.multiply(Formulas.T,axis=0).sum(axis=1)
	RelativeFormula = AbsoluteFormula/AbsoluteFormula["C"]
	return RelativeFormula[RelativeFormula>1e-3]
