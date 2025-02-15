import re

import logging
log = logging.getLogger(__name__)

import coralme
from collections import Counter

def get_remaining_complex_elements(model, met, modification_formulas):
	# get base complex and modifications, without compartment ID
	regex = r'_mod_([A-Za-z0-9_]*\(\d*\))'
	regex = re.compile(regex)

	tmp_met = coralme.core.component.Metabolite('tmp_met')
	base_complex = met.id.split('_mod_')[0]
	components = re.findall(regex, met.id)
	elements = Counter()

	# If a the completely unmodified complex is present in the model
	# has a formula, initialize the elements dictionary with that
	mets = model.metabolites
	if base_complex in mets and mets.get_by_id(base_complex).formula:
		elements.update(mets.get_by_id(base_complex).elements)

	# get component and stoichiometry per modification
	regex = r'([A-Za-z0-9_]*)\((\d*)\)'
	regex = re.compile(regex)
	for component in components:
		new_elements = elements.copy()
		new_complex = '_mod_'.join([base_complex, component])

		if new_complex in mets and mets.get_by_id(new_complex).formula:
			# default to new_complex elements if both new and old exist
			if base_complex in mets and mets.get_by_id(base_complex).formula:
				new_elements = Counter()
			formula = mets.get_by_id(new_complex).formula
			tmp_met.formula = formula
			new_elements.update(tmp_met.elements)

		## Net effect of an SH modification is adding a Sulfur to elements
		## SH should be in the me_mets file
		#elif 'SH(1)' in component:
			#new_elements['S'] += 1

		elif 'glycyl(1)' in component:
			new_elements['H'] -= 1

		# modifies O- to SH
		elif component == 'cosh(1)':
			new_elements['O'] -= 1
			new_elements['S'] += 1
			new_elements['H'] += 1

		elif component in modification_formulas: # Previously, a cofactor stoichiometry of 1 could be omitted
			formula = modification_formulas[component]
			tmp_met.formula = formula
			new_elements.update(tmp_met.elements)

		elif 'Oxidized(1)' in component and 'FLAVODOXIN' not in base_complex:
			new_elements.update({'H': -2})

		elif '(' in component: # FIXED: Check new complex ID convention
			component, value = re.findall(regex, component)[0]
			if component in modification_formulas:
				tmp_met.formula = modification_formulas[component]
			elif component + '_c' in mets:
				tmp_met.formula = mets.get_by_id(component + '_c').formula
			else:
				logging.warning('No formula found for modification \'{:s}\' either in me.metabolites nor metabolites.txt input.'.format(component))
				continue

			for e, v in tmp_met.elements.items():
				new_elements[e] += v * float(value)

		if elements == new_elements and 'FLAVODOXIN' not in base_complex:
			logging.warning('The stoichiometry of \'{:s}\' is identical to \'{:s}\'. Check if it is the correct behavior.'.format(met.id, base_complex))

		base_complex = '_mod_'.join([base_complex, component]) # DO NOT DELETE
		elements = new_elements.copy()

	return elements

def add_remaining_complex_formulas(model, modification_formulas):
	"""
	Add formula to complexes that are not formed from a complex formation
	reaction (i.e., complexes involved in metabolic reactions)
	"""
	element_dict = {}

	# Reset all formulas first
	complex_list = []
	for met in model.metabolites:
		# If met is not a complex or not formed by a complex formation reaction, do not reset
		if not isinstance(met, coralme.core.component.Complex) or met.id in model.process_data:
			continue

		for r in met.reactions:
			if hasattr(r, 'update'):
				r.update()

		met.formula = None
		met.elements = {}
		complex_list.append(met)

	# Get formulas only for complexes without complex formation reaction
	for met in complex_list:
		element_dict[met] = get_remaining_complex_elements(model, met, modification_formulas)

	# Adding elements for complexes dynamically can change function output
	# Update all of them after
	for met, elements in element_dict.items():
		coralme.util.massbalance.elements_to_formula(met, elements)
