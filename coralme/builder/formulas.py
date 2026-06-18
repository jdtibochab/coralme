import re
import tqdm
import logging
log = logging.getLogger(__name__)

import coralme
from coralme.core.extended_classes import log_format, bar_format
import collections
import copy

# WARNING: Originally on main.py
def add_remaining_complex_formulas(me, df_mets):
	"""
	Add formula to complexes that are not formed from a complex formation
	reaction (i.e., complexes involved in metabolic reactions)
	"""
	modification_formulas = df_mets[df_mets['type'].str.match('COFACTOR|MOD|MODIFICATION')]
	modification_formulas = dict(zip(modification_formulas['me_id'], modification_formulas['formula']))
	me.global_info['modification_formulas'] = modification_formulas

	modification_charges = {}
	if 'charge' in df_mets.columns:
		modification_charges = df_mets[df_mets['type'].str.match('COFACTOR|MOD|MODIFICATION')]
		modification_charges = dict(zip(modification_charges['me_id'], modification_charges['charge']))
		modification_charges = { k:float(v) for k,v in modification_charges.items() if v != '' }
		me.global_info['modification_charges'] = modification_charges

	for met in [ x for x in me.metabolites if '_mod_' in x.id and isinstance(x, coralme.core.component.Complex)]:
		met.formula = None
		met.elements = {}

		base_complex = met.id.split('_mod_')[0]
		if not me.metabolites.has_id(base_complex) and not hasattr(met, 'base_complex_elements'):
			logging.warning('WARNING: Base complex \'{:s}\' does not exist. Reload an unpruned coralME model.'.format(base_complex))
			continue # Without continue, add elemental contributions of the modifications
		if not hasattr(met, 'base_complex_elements'):
			met.base_complex_elements = collections.Counter(me.metabolites.get_by_id(base_complex).elements)
		base_complex_elements = collections.Counter(met.base_complex_elements)
		
		for mod in met.id.split('_mod_')[1:]:
			#for num in range(int(mod.rstrip(')').split('(')[1])):
			mod_elements = None
			mod_name = mod.split('(')[0]

			if mod_name in modification_formulas:
				mod_elements = coralme.builder.helper_functions.parse_composition(modification_formulas[mod_name])
				# 2fe2s_c and 4fe4s_c appear as free metabolites in reactions and need to have formula for correct mass balance determination
				if me.metabolites.has_id(mod_name + '_c') and me.metabolites.get_by_id(mod_name + '_c').formula is None:
					me.metabolites.get_by_id(mod_name + '_c').formula = modification_formulas[mod_name]
					logging.warning('WARNING: New formula for \'{:s}\' was updated using me_mets.txt file.'.format(mod_name + '_c'))
					me.metabolites.get_by_id(mod_name + '_c').charge = modification_charges.get(mod_name, 0)
					logging.warning('WARNING: New charge for \'{:s}\' was updated using me_mets.txt file.'.format(mod_name + '_c'))
				logging.warning('INFO: Elemental contribution for \'{:s}\' calculated from me_mets.txt file.'.format(mod_name))
			
			elif me.metabolites.has_id(mod_name + '_c') and me.metabolites.get_by_id(mod_name + '_c').formula is not None:
				mod_elements = me.metabolites.get_by_id(mod_name + '_c').elements
				logging.warning('INFO: Elemental contribution for \'{:s}\' calculated from metabolite formula.'.format(mod_name))
				
			# WARNING: electron carriers can, assuming they are neutral, transfer also protons
			# WARNING: Ferredoxins only transfer electrons; thioredoxins and others transfer protons and electrons.
			elif mod.startswith('Oxidized'):
				if base_complex in me.global_info['electron_transfers'].get('ferredoxins', []):
					mod_elements = {'H': 0}
					logging.warning('INFO: Elemental contribution in ferredoxin-type Complex \'{:s}\' calculated manually.'.format(met.id))
				elif base_complex in me.global_info['electron_transfers'].get('cytochromes', []):
					mod_elements = {'H': 0}
					logging.warning('INFO: Elemental contribution in cytochrome-type Complex \'{:s}\' calculated manually.'.format(met.id))
				elif 'Oxidized(1)' == mod and base_complex not in me.global_info['electron_transfers'].get('flavodoxins', ['FLAVODOXIN']):
					mod_elements = {'H': -2}
					logging.warning('INFO: Elemental contribution in thioredoxin-type Complex \'{:s}\' calculated manually.'.format(met.id))
				elif 'Oxidized(2)' == mod and base_complex not in me.global_info['electron_transfers'].get('flavodoxins', ['FLAVODOXIN']):
					mod_elements = {'H': -4}
					logging.warning('INFO: Elemental contribution in thioredoxin-type Complex \'{:s}\' calculated manually.'.format(met.id))
				# TODO: is the fmn cofactor in flavodoxin neutral?
				# WARNING: flavodoxin homologs might have a different base_complex ID compared to the ecolime model
				elif 'Oxidized(1)' == mod and base_complex in me.global_info['electron_transfers'].get('flavodoxins', ['FLAVODOXIN']):
					mod_elements = {'H': 0}
					logging.warning('INFO: Elemental contribution in flavodoxin-type Complex \'{:s}\' calculated manually.'.format(met.id))
				else:
					logging.warning('WARNING: Elemental contribution in \'{:s}\' could not be determined.'.format(met.id))
					logging.warning('INFO: Please check configuration file and add the base complex into the \'electron_transfers\' key.')

			# WARNING: Negative elemental contributions cannot be set in the metabolites.txt input file
			elif 'glycyl(1)' == mod:
				mod_elements = {'H': -1}
				logging.warning('INFO: Elemental contribution for \'{:s}\' calculated manually.'.format(mod_name))
			elif 'cosh(1)' == mod:
				mod_elements = {'H': +1, 'O': -1, 'S': +1}
				logging.warning('INFO: Elemental contribution for \'{:s}\' calculated manually.'.format(mod_name))
			else:
				if not me.metabolites.has_id(mod_name + '_c'):
					logging.warning('WARNING: Metabolite does not exist in M-model. Elemental contribution for \'{:s}\' could not be determined.'.format(mod_name))				
				else:
					logging.warning('WARNING: Elemental contribution for \'{:s}\' could not be determined. Please check me_mets.txt file.'.format(mod_name))

			if mod_elements:
				mod_elements = collections.Counter(mod_elements)
				mod_elements = { k:v * int(re.findall(r'\((\d+)\)', mod)[0]) for k,v in mod_elements.items() }
				base_complex_elements.update(mod_elements)
			else:
				logging.warning('WARNING: Attempt to calculate a corrected formula for \'{:s}\' failed. Please check if it is the correct behaviour, or if the modification \'{:s}_c\' exists as a metabolite in the ME-model or a formula is included in the me_mets.txt file.'.format(met.id, mod_name))
		
		complex_elements = { k:base_complex_elements[k] for k in sorted(base_complex_elements) if base_complex_elements[k] != 0 }
		met.formula = ''.join([ '{:s}{:d}'.format(k, v) for k,v in complex_elements.items() ])
		met.elements = coralme.builder.helper_functions.parse_composition(met.formula)
		logging.warning('INFO: Setting new formula for \'{:s}\' to \'{:s}\' successfully.'.format(met.id, met.formula))

	# Update a second time to incorporate all of the metabolite formulas correctly
	for data in tqdm.tqdm(me.subreaction_data.query(r'(?!^\w\w\w_addition_at_\w\w\w$)'), 'Recalculation of the elemental contribution in SubReactions...', bar_format = bar_format):
		data._element_contribution = data.calculate_element_contribution()

	# Update reactions affected by formula update
	for r in tqdm.tqdm(me.reactions.query('^formation_'), 'Updating all FormationReactions...', bar_format = bar_format):
		r.update()

	for r in tqdm.tqdm(me.reactions.query('_mod_lipoyl'), 'Updating FormationReactions involving a lipoyl prosthetic group...', bar_format = bar_format):
		r.update()

	for r in tqdm.tqdm(me.reactions.query('_mod_glycyl'), 'Updating FormationReactions involving a glycyl radical...', bar_format = bar_format):
		r.update()

	return None
