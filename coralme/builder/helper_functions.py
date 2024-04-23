#!/usr/bin/python3
import re
import copy
import pandas
import tqdm
import cobra
import coralme


# from cobrame without changes
def get_base_complex_data(model, complex_id):
	"""
	If a complex is modified in a metabolic reaction it will not
	have a formation reaction associated with it. This function returns
	the complex data of the "base" complex, which will have the subunit
	stoichiometry of that complex
	"""

	# First try unmodified complex id
	try_1 = complex_id.split('_mod_')[0]
	if try_1 in model.process_data:
		return model.process_data.get_by_id(try_1)

	count = 0
	try_2 = complex_id.split('_mod_')[0] + '_'
	for i in model.process_data.query(try_2):
		if isinstance(i, coralme.core.processdata.ComplexData):
			count += 1
			data = i

	if count == 0:
		raise UserWarning('No base complex found for \'{:s}\'.'.format(complex_id))

	elif count > 1:
		raise UserWarning('More than one possible base complex found for \'{:s}\'.'.format(complex_id))

	return data

# from cobra with changes
def parse_composition(tmp_formula) -> dict:
	element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")
	"""Break the chemical formula down by element."""
	#tmp_formula = self.formula
	# commonly occuring characters in incorrectly constructed formulas
	if "*" in tmp_formula:
		warn(f"invalid character '*' found in formula '{self.formula}'")
		tmp_formula = self.formula.replace("*", "")
	if "(" in tmp_formula or ")" in tmp_formula:
		warn(f"parenthesis found in formula '{self.formula}'")
		return
	composition = {}
	parsed = element_re.findall(tmp_formula)
	for element, count in parsed:
		if count == "":
			count = 1
		else:
			try:
				count = float(count)
				int_count = int(count)
				if count == int_count:
					count = int_count
				else:
					warn(f"{count} is not an integer (in formula {self.formula})")
			except ValueError:
				warn(f"failed to parse {count} (in formula {self.formula})")
				#self.elements = {}
				return
		if element in composition:
			composition[element] += count
		else:
			composition[element] = count
	return composition

## Network traversing
def substitute_value(m,coeff):
	"""Substitute coefficient with a value"""
	if not hasattr(coeff,'subs'):
		return coeff
	return coeff.subs([(m._model.mu, 1e-6)])

def get_next_from_type(l,t):
	"""Get the first object that matches a type"""
	try:
		if isinstance(l,dict):
			return next(i for i,v in l.items() if isinstance(i,t) and substitute_value(i,v) > 0)
		return next(i for i in l if isinstance(i,t))
	except:
		return set()

def find_complexes(m, seen = set()):
	if not m:
		return set()
#	 print(m.id)
	if m in seen:
		return set()
#	 print(type(m))

	seen.add(m)

	# Metabolite objects
	if isinstance(m,coralme.core.component.TranslatedGene):
		cplxs = set()
		for r in m.reactions:
			if substitute_value(m,r.metabolites[m] > 0):
				continue
			cplxs = cplxs | find_complexes(r, seen=seen)
		return cplxs
	if isinstance(m,coralme.core.component.TranscribedGene):
		translated_protein = m.id.replace('RNA_','protein_')
		if translated_protein in m._model.metabolites:
			return find_complexes(m._model.metabolites.get_by_id(translated_protein), seen=seen)
		cplxs = set()
		for r in m.reactions:
			if substitute_value(m,r.metabolites[m] > 0):
				continue
			cplxs = cplxs | find_complexes(r, seen=seen)
		return cplxs
	if isinstance(m,coralme.core.component.ProcessedProtein):
		cplxs = set()
		for r in m.reactions:
			if substitute_value(m,r.metabolites[m] > 0):
				continue
			cplxs = cplxs | find_complexes(r, seen=seen)
		return cplxs

	if isinstance(m,coralme.core.component.Complex) or isinstance(m,coralme.core.component.GenericComponent) or isinstance(m,coralme.core.component.GenerictRNA):
		other_formations = [r for r in m.reactions if (isinstance(r,coralme.core.reaction.ComplexFormation) or isinstance(r,coralme.core.reaction.GenericFormationReaction)) and substitute_value(m,r.metabolites[m]) < 0]
#		 print(other_formations)
		cplxs = set([m])
		if other_formations:
			for r in other_formations:
				cplxs = cplxs | find_complexes(r, seen=seen)
		return cplxs
#		 return set([m])

	# Reaction objects
	if isinstance(m,coralme.core.reaction.PostTranslationReaction):
		return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.ProcessedProtein), seen=seen)
	if isinstance(m,coralme.core.reaction.ComplexFormation):
		return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.Complex), seen=seen)
	if isinstance(m,coralme.core.reaction.GenericFormationReaction):
		return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.GenericComponent), seen=seen)
	if isinstance(m,coralme.core.reaction.tRNAChargingReaction):
		return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.GenerictRNA), seen=seen)
	if isinstance(m,coralme.core.reaction.MetabolicReaction):
		return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.Complex), seen=seen) | \
				find_complexes(get_next_from_type(m.metabolites,coralme.core.component.GenericComponent), seen=seen)
	return set()

def get_functions(cplx):
	functions = set()
	for r in cplx.reactions:
		if isinstance(r,coralme.core.reaction.MetabolicReaction) and hasattr(r,'subsystem'):
			if r.subsystem:
				functions.add('Metabolic:' + r.subsystem)
				continue
			functions.add('Metabolic: No subsystem')
			continue
		if isinstance(r,coralme.core.reaction.TranslationReaction):
			functions.add('Translation')
		elif isinstance(r,coralme.core.reaction.TranscriptionReaction):
			functions.add('Transcription')
		elif isinstance(r,coralme.core.reaction.tRNAChargingReaction):
			functions.add('tRNA-Charging')
		elif isinstance(r,coralme.core.reaction.PostTranslationReaction):
			functions.add('Post-translation')
		elif isinstance(r,coralme.core.reaction.SummaryVariable):
			functions.add('Biomass')
	return functions

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

# Other
# Originally developed by JDTB, UCSD (2022)
def close_sink_and_solve(rxn_id):
	global _model
	global _muf
	model = copy.deepcopy(_model)
	model.reactions.get_by_id(rxn_id).bounds = (0, 0)
	model.optimize(max_mu = _muf, maxIter = 1)
	if not model.solution:
		return (rxn_id, False)
	else:
		return (rxn_id, True)

def change_reaction_id(model,old_id,new_id):
	old_rxn = model.reactions.get_by_id(old_id)
	rxn = cobra.Reaction(new_id)
	model.add_reactions([rxn])
	for k,v in old_rxn.metabolites.items():
		rxn.add_metabolites({k.id:v})
	rxn.bounds = old_rxn.bounds
	rxn.name = old_rxn.name
	rxn.subsystem = old_rxn.subsystem
	rxn.notes = old_rxn.notes
	rxn.gene_reaction_rule = old_rxn.gene_reaction_rule
	model.remove_reactions([old_rxn])

def get_metabolites_from_pattern(model,pattern):
	met_list = []
	for met in model.metabolites:
		if pattern in met.id:
			met_list.append(met.id)
	return met_list

def evaluate_lp_problem(Sf, Se, lb, ub, keys, atoms):
	lb = [ x(*[ keys[x] for x in list(atoms) ]) if hasattr(x, '__call__') else float(x.xreplace(keys)) if hasattr(x, 'subs') else x for x in lb ]
	ub = [ x(*[ keys[x] for x in list(atoms) ]) if hasattr(x, '__call__') else float(x.xreplace(keys)) if hasattr(x, 'subs') else x for x in ub ]
	Se = { k:x(*[ keys[x] for x in list(atoms) ]) if hasattr(x, '__call__') else float(x.xreplace(keys)) if hasattr(x, 'subs') else x for k,x in Se.items() }

	Sf.update(Se)

	return Sf, Se, lb, ub

def is_producible(me,met,growth_key_and_value):
	if met not in me.metabolites:
		return False
	r = add_exchange_reactions(me, [met], prefix = 'DM_')[0]
	r.bounds = (1e-16,1e-16)
	if me.check_feasibility(keys = growth_key_and_value):
		r.remove_from_model()
		return True
	else:
		r.remove_from_model()
		return False


def get_transport_reactions(model,met_id,comps=['e','c']):
	from_met = re.sub('_[a-z]$','_'+comps[0],met_id)
	to_met = re.sub('_[a-z]$','_'+comps[1],met_id)

	if isinstance(model,coralme.core.model.MEModel):
		reaction_type = ['MetabolicReaction']
	else:
		reaction_type = 0
	prod_rxns = [rxn.id for rxn in get_reactions_of_met(model,to_met,s=1,verbose=0,only_types=reaction_type)]
	cons_rxns = [rxn.id for rxn in get_reactions_of_met(model,from_met,s=-1,verbose=0,only_types=reaction_type)]
	transport_rxn_ids = list(set(prod_rxns)&set(cons_rxns))

	return [model.reactions.get_by_id(rxn_id) for rxn_id in transport_rxn_ids]

def get_all_transport_of_model(model):
	transport_reactions = []
	for r in tqdm.tqdm(model.reactions):
		comps = r.get_compartments()
		if len(comps) > 1:
			transport_reactions.append(r.id)
	return list(set(transport_reactions))


def format_kcats_from_DLKcat(df):
# df = pandas.read_csv("./bacillus/building_data/bacillus_rxn_kcats.tsv",sep='\t',index_col=0).set_index("reaction")

	df2 = pandas.DataFrame(columns=["direction","complex","mods","keff"])
	for r,keff in df['Kcat value (1/s)'].items():
		d = {}
		rid,cplx = re.split("_FWD_|_REV_",r)
		base_cplx = cplx.split("_mod_")[0]
		mods = cplx.split(base_cplx)[1]
		mods = " AND ".join(mods.split("_mod_")[1:])
		direc = "FWD" if "FWD" in r else "REV"
		d[rid] = {
			"direction":direc,
			"complex":base_cplx,
			"mods":mods,
			"keff":keff
		}
		df2 = pandas.concat([df2,pandas.DataFrame.from_dict(d).T])
	df2.index.name = 'reaction'
	return df2
#	 df2.to_csv("./bacillus/building_data/keffs.txt",sep='\t')

# Flux analysis migration
def exchange_single_model(me, flux_dict = 0, solution=0):
	return coralme.util.flux_analysis.exchange_single_model(me,
														 flux_dict = flux_dict,
														 solution=solution)

def get_met_coeff(stoich,growth_rate,growth_key='mu'):
	return coralme.util.flux_analysis.get_met_coeff(stoich,
												 growth_rate,
												 growth_key=growth_key)

def summarize_reactions(model,met_id,only_types=(),ignore_types = ()):
	return coralme.util.flux_analysis.summarize_reactions(model,
													   met_id,
													   only_types=only_types,
													   ignore_types = ignore_types)

def flux_based_reactions(model,
						 met_id,
						 only_types=(),
						 ignore_types = (),
						 threshold = 0.,
						 flux_dict=0,
						 solution = None,
						 keffs=False,
						 verbose=False):
	return coralme.util.flux_analysis.flux_based_reactions(model,
						 met_id,
						 only_types=only_types,
						 ignore_types = ignore_types,
						 threshold = threshold,
						 flux_dict=flux_dict,
						 solution = solution,
						 keffs=keffs,
						 verbose=verbose)

def get_reactions_of_met(me,met,s = 0, ignore_types = (),only_types = (), verbose = False,growth_key='mu'):
	return coralme.util.flux_analysis.get_reactions_of_met(me,
														met,
														s = s,
														ignore_types = ignore_types,
														only_types = only_types,
														verbose = verbose,
														growth_key = growth_key)
