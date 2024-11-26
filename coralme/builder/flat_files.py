import tqdm
import json
import cobra
import pandas

import logging
log = logging.getLogger(__name__)

import coralme
from collections import defaultdict, Counter

def fix_id(id_str):
	return id_str.replace('_DASH_', '__')

def read(filename) -> pandas.DataFrame:
	return pandas.read_csv(filename, sep = '\t', index_col = None, comment = '#', skip_blank_lines = True, dtype = str, keep_default_na = False)

def get_complex_subunit_stoichiometry(complex_stoichiometry, rna_components = set()) -> dict:
	"""Returns dictionary of prot/prot, prot/rna complexes: {stoichiometry: {locus_tag: count}}
	"""
	complex_stoichiometry_dict = {}
	complex_stoichiometry.columns = ['Name', 'Stoichiometry']

	for key, row in complex_stoichiometry.iterrows():
		if key in complex_stoichiometry_dict.keys():
			logging.warning('Complex \'{:s}\' is present twice or more times in the DataFrame and the repeated item was ignored.'.format(key))
		else:
			complex_stoichiometry_dict[key] = {}

		for bnums in row['Stoichiometry'].split(' AND '):
			bnum, num = bnums.rstrip(')').split('(')
			stoichiometry = float(num) if not num == '' else 1.
			# add prefix to complex component
			#prefix = 'protein_' if bnum not in rna_components else 'RNA_'
			prefix = 'RNA_' if bnum in rna_components else 'protein_' if 'generic' not in bnum else ''
			complex_stoichiometry_dict[key][prefix + bnum] = stoichiometry

	return complex_stoichiometry_dict

def get_reaction_matrix_dict(matrix_df, compartments = {}, complex_set = set()) -> dict:
	"""Return dictionary representation of the metabolic reaction matrix.
	Updates metabolite id with compartment if not contained in complex_list
	"""
	#matrix_df.columns = ['Reaction', 'Metabolites', 'Compartment', 'Stoichiometry']
	matrix_df.columns = ['Reaction', 'Metabolites', 'Stoichiometry']
	#matrix_df.replace({'No_Compartment': 'Cytoplasm'}, inplace = True)

	# old code; replaced to get compartments from the me-model object
	#compartments = {'Cytosol': 'c', 'Periplasm': 'p', 'Extra-organism': 'e'}
	metabolic_reaction_dict = defaultdict(dict)
	for idx, row in matrix_df.iterrows():
		reaction = fix_id(row['Reaction'])
		metabolite = fix_id(row['Metabolites'])
		stoichiometry = row['Stoichiometry']
		#if compartments.get(row['Compartment']) is None:
			#compartment_id = ''
		#else:
			#compartment_id = '_{:s}'.format(compartments.get(row['Compartment']))

		# use compartment to append appropriate suffix
		#if metabolite.split('_mod_')[0] not in complex_set:
			#metabolite += compartment_id
		metabolic_reaction_dict[reaction][metabolite] = float(stoichiometry)

	return metabolic_reaction_dict

def get_complex_modifications(reaction_matrix, protein_complexes, complex_mods, compartments = set()) -> dict:
	"""
	"""
	#complex_dct = get_complex_subunit_stoichiometry(protein_complexes)
	#complex_set = set(complex_dct.keys())

	#ignored_complexes = set()
	#for met_stoich in coralme.builder.flat_files.get_reaction_matrix_dict(reaction_matrix).values():
		#for met, value in met_stoich.items():
			#if len(met.split('_mod_')[1:]) >= 1:
				#ignored_complexes.add(met)

	## correct ignored complexes
	#ignored_complexes = [ x for x in ignored_complexes
		#if x.split('_mod_')[1:] != ['pydx5p(1)'] and x.split('_mod_')[1:] != ['pydx5p(2)'] and x.split('_mod_')[1:] != ['4fe4s(2)']]

	new_mod_dict = {}
	for key, value in complex_mods.T.to_dict().items():
		# QueG_mono_mod_4fe4s(2), CPLX0-782_mod_4fe4s(2), CPLX0-246_CPLX0-1342_mod_pydx5p(1) and IscS_mod_pydx5p(2) must be exceptions
# 		if key in ignored_complexes:
# 			continue
		key = key.replace('_DASH_', '__')
		new_mod_dict[key] = {}
		new_mod_dict[key]['core_enzyme'] = value['Core_enzyme']
		new_mod_dict[key]['modifications'] = {}
		for mods in value['Modifications'].split(' AND '):
			mod, num_mods = mods.rstrip(')').split('(')
			if num_mods == '':
				num_mods = 1.
			else:
				num_mods = float(num_mods)

			mod = mod.replace('_DASH_', '__')
			new_mod_dict[key]['modifications'][mod + '_c'] = -num_mods

	return new_mod_dict

def get_reaction_to_complex(m_model, enz2rxn, modifications = True):
	"""anything not in this dict is assumed to be an orphan"""
	enz2rxn.columns = ['Complexes']
	enz2rxn = enz2rxn.apply(lambda x: x.replace('DASH', ''))

	rxn_to_complex_dict = defaultdict(set)
	for reaction, complexes in enz2rxn.itertuples():
		for cplx in complexes.split(' OR '):
			if modifications:
				rxn_to_complex_dict[reaction].add(cplx)
			else:
				rxn_to_complex_dict[reaction].add(cplx.split('_mod_')[0])

	for reaction in m_model.reactions:
		if 's0001' in reaction.gene_reaction_rule:
			rxn_to_complex_dict[reaction.id].add(None)

	return rxn_to_complex_dict

def remove_compartment(id_str):
	# original ecolime remove_compartment function
	#return id_str.replace('_c', '').replace('_p', '').replace('_e', '')
	return '_'.join(id_str.split('_')[:-1]) # compartment ID follows the last underscore

def process_reaction_matrix_dict(reaction_matrix, cplx_data, me_compartments = set()):
	reaction_matrix['Reaction'] = reaction_matrix['Reaction'].str.split(',')
	reaction_matrix = reaction_matrix.explode('Reaction')
	complex_dct = get_complex_subunit_stoichiometry(cplx_data)
	reaction_matrix_dict = get_reaction_matrix_dict(reaction_matrix, compartments = me_compartments, complex_set = set(complex_dct.keys()))

	return reaction_matrix_dict

def process_m_model(
	m_model, # mets_data, rxns_data, # all data
	mets_data, rxns_data, reaction_matrix, cplx_data, # only new data
	me_compartments = set(), defer_to_rxn_matrix = list(), repair = True):

	# copy the M-model object
	m_model = m_model.copy()

	# old code, known bug
	#for rxn in m_model.reactions:
		#if rxn.id.startswith('EX_') or rxn.id.startswith('DM_') or rxn.id.startswith('SK_') or rxn.id.startswith('BIOMASS_'):
			#continue

		## old code
		##if rxn.id not in reaction_matrix_dict.keys() or rxn.id in defer_to_rxn_matrix:
		#if rxn.id in defer_to_rxn_matrix:
			#rxn.remove_from_model(remove_orphans = True)
			#logging.warning('The MetabolicReaction \'{:s}\' (using \'defer_to_rxn_matrix\') was removed from the M-model metabolic network.'.format(rxn.id))

	#m_model.remove_reactions([ m_model.reactions.get_by_id(rxn) for rxn in defer_to_rxn_matrix ])
	rxns_to_remove = [ m_model.reactions.get_by_id(rxn) for rxn in defer_to_rxn_matrix if m_model.reactions.has_id(rxn) ]
	m_model.remove_reactions(rxns_to_remove)
# 	mets_to_remove = [ m for m in m_model.metabolites if len(m.reactions) == 0 ]
# 	m_model.remove_metabolites(mets_to_remove)

	# correct stoichiometry if rxn.check_mass_balance() == {'charge': +1.0, 'H': +1.0} or integer multiples of it
	count = 0
	for rxn in m_model.reactions:
		if rxn.id.startswith(('EX_', 'DM_', 'SK_')):
			continue

		check = rxn.check_mass_balance()
		if check == {}:
			continue
		elif set(check.keys()) == {'H', 'charge'} and check['charge'] == check['H']:
			compt = rxn.get_compartments()
			if len(compt) == 1:
				metabolites = Counter(m_model.reactions.get_by_id(rxn.id).metabolites)
				metabolites.update(Counter({ m_model.metabolites.get_by_id('h_{:s}'.format(compt[0])) : -1 * check['H'] }))

				m_model.reactions.get_by_id(rxn.id)._metabolites = { k:v for k,v in metabolites.items() if v != 0. }
				logging.warning('Stoichiometry for \'{:s}\' was corrected to mass balance protons.'.format(rxn.id))
			else:
				logging.warning('Stoichiometry for \'{:s}\' was not corrected due to more than one compartment detected in the reaction. Please check and correct if needed.'.format(rxn.id))
			count += 1
		else:
			logging.warning('Reaction \'{:s}\' is not mass/charge balanced. Please check and correct if needed.'.format(rxn.id))
			count += 1
	if count == 0:
		logging.warning('No mass/charge imbalance was detected in M-model (excluding prefixed EX, DM, and SK reactions).')

	# remove unused genes
	cobra.manipulation.delete.remove_genes(m_model, [ x for x in m_model.genes if len(x.reactions) == 0 ], remove_reactions = False)
	cobra.manipulation.delete.prune_unused_reactions(m_model) # reactions without metabolites
	# cobra.manipulation.delete.prune_unused_metabolites(m_model) # orphan metabolites

	# met_data DataFrame
	mets_data = mets_data[mets_data['type'].isin(['ADD', 'REPLACE'])]
	mets_data.columns = ['me_id', 'name', 'formula', 'compartment', 'type']
	# WARNING: Do we need to rename the index?
	#mets_data.rename(index = lambda x: x.replace('_DASH_', '__'), inplace = True)

	# process protein_complex DataFrame, and...
	complex_dct = get_complex_subunit_stoichiometry(cplx_data)
	# ...process reaction_matrix DataFrame
	reaction_matrix_dict = get_reaction_matrix_dict(reaction_matrix, compartments = me_compartments, complex_set = set(complex_dct.keys()))

	for rxn_id in reaction_matrix_dict:
		if rxn_id in m_model.reactions:
			logging.warning('The MetabolicReaction \'{:s}\' (using \'reaction_matrix.txt\') is already present in the M-model and it will be replaced.'.format(rxn_id))
			m_model.remove_reactions([rxn_id], remove_orphans = False)

		# Metabolites need to be added into the M-model first
		rxn_stoichiometry = reaction_matrix_dict[rxn_id]
		for met in rxn_stoichiometry:
			try:
				met_obj = m_model.metabolites.get_by_id(met)
			except KeyError:
				if met.startswith('RNA_'):
					met_obj = coralme.core.component.TranscribedGene(str(met), 'tRNA', '')
				elif met.startswith('generic_tRNA'):
					met_obj = coralme.core.component.GenerictRNA(str(met))
				else:
					met_obj = coralme.core.component.Metabolite(str(met))
				m_model.add_metabolites([met_obj])

			met_id = remove_compartment(met_obj.id)
			if met_id in mets_data.index and not met_obj.formula:
				met_obj.name = mets_data.loc[met_id, 'name']
				met_obj.formula = mets_data.loc[met_id, 'formula']

		# Add new reactions into the M-model
		rxn = coralme.core.reaction.MEReaction(rxn_id)
		m_model.add_reactions([rxn])
		rxn.add_metabolites(rxn_stoichiometry)
		if rxn_id in rxns_data.index:
			reversible = rxns_data.loc[rxn_id, 'is_reversible']
			subsystem = 'Not Determined' if rxns_data.loc[rxn_id, 'subsystems'] == 'False' else rxns_data.loc[rxn_id, 'subsystems']
		else:
			reversible = 'False'
			logging.warning('Unable to determine MetabolicReaction \'{:s}\' reversibility. Default value is \'False\'.'.format(rxn_id))
			subsystem = 'Not Determined'
			logging.warning('Unable to determine MetabolicReaction \'{:s}\' subsystem. Default value is \'Not Determined\'.'.format(rxn_id))

		rxn.lower_bound = -1000.0 if reversible.lower() == 'true' else 0.0
		logging.warning('The MetabolicReaction \'{:s}\' was created into the M-model (using \'reaction_matrix.txt\').'.format(rxn_id))
		logging.warning('MetabolicReaction \'{:s}\' reversibility is set to \'{:s}\' (using \'reaction_matrix.txt\').'.format(rxn_id, reversible))
		logging.warning('MetabolicReaction \'{:s}\' subsystem is set to \'{:s}\' (using \'reaction_matrix.txt\').'.format(rxn_id, subsystem))
		rxn.subsystem = subsystem

	# m_to_me_map DataFrame
	#m_to_me_map.columns = ['me_name']
	m_to_me_map = mets_data[mets_data['me_id'].notnull()]
	#m_to_me_map.columns = ['me_id', 'name', 'formula', 'compartment', 'data_source']
	m_to_me_map.rename(lambda x: x.replace('_DASH_', '__'), inplace = True)

	m_model.add_metabolites([ coralme.core.component.Complex(id = x) for x in set(m_to_me_map['me_id']) ])

	for rxn in m_model.reactions:
		#met_id = remove_compartment(met.id)
		# old code. mets_data contains all the metabolites (m+me)
		#if met_id not in mets_data.index and met_id in m_to_me_map.index:
		#if met_id in mets_data.index:
			# old code. m_to_me_map contains the new metabolite ID that maps the M-model metabolite
			#met_id = m_to_me_map.loc[met.id, 'me_name']
		#if mets_data.loc[met_id, 'me_id'] != '' and mets_data.loc[met_id, 'me_id'] != 'N/A':
		for met, coeff in rxn.metabolites.items():
			if met.id in m_to_me_map.index and m_to_me_map.loc[met.id, 'me_id'] != '':
				#met.id = mets_data.loc[met_id, 'me_id']
				#met_obj = coralme.core.component.Complex(id = m_to_me_map.loc[met.id, 'me_id'])
				#met_obj.name = m_to_me_map.loc[met.id, 'name']
				#met_obj.formula = m_to_me_map.loc[met.id, 'formula']

				new_id = m_model.metabolites.get_by_id(m_to_me_map.loc[met.id, 'me_id'])
				old_stoich = rxn.metabolites[m_model.metabolites.get_by_id(met.id)]
				rxn.add_metabolites({ new_id : old_stoich })
				logging.warning('Metabolite \'{:s}\' was replaced with \'{:s}\' in MetabolicReaction \'{:s}\'.'.format(met.id, m_to_me_map.loc[met.id, 'me_id'], rxn.id))

	m_model.remove_metabolites([ m_model.metabolites.get_by_id(x) for x in m_to_me_map[m_to_me_map['type'].str.fullmatch('REPLACE|REMOVE')].index if m_model.metabolites.has_id(x) ])

	# Add new metabolites (ME-metabolites) with properties into the "M-model"
	for m_met_id in m_to_me_map.index:
		me_met_id = m_to_me_map.loc[m_met_id, 'me_id']
		if m_model.metabolites.has_id(me_met_id):
			met_obj = m_model.metabolites.get_by_id(me_met_id)
		else:
			met_obj = coralme.core.component.Metabolite(str(me_met_id))

		met_obj.name = m_to_me_map.loc[m_met_id, 'name']
		met_obj.formula = m_to_me_map.loc[m_met_id, 'formula'] if m_to_me_map.loc[m_met_id, 'formula'] != '' else None
		met_obj.compartment = m_met_id.split('_')[-1]
		# WARNING: This "fails" silently if the met_obj ID already exists
		m_model.add_metabolites([met_obj])

	if repair:
		m_model.repair()

	return m_model

def get_m_model(
	modelID = 'M-model-from-ME-model',
	metabolites = 'metabolites.txt',
	reactions = 'reactions.txt',
	reaction_matrix = 'reaction_matrix.txt',
	prot_complexes = 'protein_complexes.txt',
	exchanges = 'exchange_bounds.txt',
	sources_and_sinks = 'reaction_matrix_sources_and_sinks.txt',
	compartments = {}
	) -> cobra.core.model.Model:

	m = cobra.core.model.Model(modelID)

	# match compartments in ME-model
	m.compartments = compartments
	compartment_lookup = {v:k for k,v in m.compartments.items()}

	met_info = pandas.read_csv(metabolites, sep = '\t', index_col = None, comment = '#', skip_blank_lines = True)
	met_info.columns = ['id', 'name', 'formula', 'compartment', 'data_source']
	met_info = met_info.set_index('id')

	complex_set = set(get_complex_subunit_stoichiometry(filename = prot_complexes).keys())

	for met_id in met_info.index:
		fixed_id = fix_id(met_id)
		for compartment in met_info.compartment[met_id].split('AND'):
			compartment = compartment.strip()
			if compartment == 'No_Compartment':
				print('Assigned {:s} to c'.format(met_id))
				compartment = m.compartments['c']
			new_met = cobra.Metabolite(
				fixed_id + '_' + compartment_lookup[compartment])
			new_met.name = met_info.name[met_id]
			new_met.formula = met_info.formula[met_id]
			m.add_metabolites(new_met)

	rxn_info = get_reaction_info_frame(reactions)
	rxn_dict = get_reaction_matrix_dict(filename = reaction_matrix, compartments = compartments, complex_set = complex_set)
	for rxn_id in rxn_info.index:
		reaction = coralme.core.reaction.MEReaction(rxn_id)
		reaction.name = rxn_info.description[rxn_id]
		for met_id, amount in rxn_dict[rxn_id].items():
			try:
				metabolite = m.metabolites.get_by_id(met_id)
			except KeyError:
				metabolite = cobra.Metabolite(met_id)
			reaction.add_metabolites({metabolite: amount})
		reaction.lower_bound = -1000. if rxn_info.is_reversible[rxn_id].lower() == 'true' else 0.
		reaction.upper_bound = 1000.
		if rxn_info.is_spontaneous[rxn_id].lower() == 'true':
			reaction.gene_reaction_rule = 's0001'
		m.add_reaction(reaction)

	sources_sinks = pandas.read_csv(sources_and_sinks, sep = '\t', index_col = None, comment = '#', skip_blank_lines = True)
	sources_sinks.columns = ['rxn_id', 'met_id', 'compartment', 'stoic']
	sources_sinks = sources_sinks.set_index('met_id')
	sources_sinks.index = [fix_id(i) for i in sources_sinks.index]

	source_amounts = pandas.read_csv(exchanges, sep = '\t', index_col = None, comment = '#', skip_blank_lines = True)
	source_amounts.columns = ['met_id', 'amount']
	source_amounts = source_amounts.set_index('met_id')
	source_amounts.index = [fix_id(i) for i in source_amounts.index]

	for met in sources_sinks.index:
		met_id = met + '_' + compartment_lookup[sources_sinks.compartment[met]]
		# EX_ or DM_ + met_id
		reaction_id = sources_sinks.rxn_id[met][:3] + met_id
		reaction = coralme.core.reaction.MEReaction(reaction_id)
		m.add_reaction(reaction)
		reaction.add_metabolites({m.metabolites.get_by_id(met_id): -1})
		# set bounds on exchanges
		if reaction.id.startswith('EX_') and met in source_amounts.index:
			reaction.lower_bound = -source_amounts.amount[met]

	return m

def get_trna_modification_targets(trna_mod) -> dict:
	trna_mod.columns = ['bnum', 'modification', 'positions', 'enzymes']

	trna_mod_dict = defaultdict(dict)
	for idx, mod in trna_mod.iterrows():
		mod_loc = '{:s}_at_{:s}'.format(mod['modification'], mod['positions'])
		trna_mod_dict[mod['bnum']][mod_loc] = 1.

	return trna_mod_dict

def get_reaction_keffs(me, keffs = 'keffs.json', verbose = True):
	def log(*args, **kwargs):
		if verbose:
			print(*args, **kwargs)

	with open(keffs, 'r') as infile:
		keffs = json.load(infile)

	new_keffs = {}
	for met_rxn in me.reactions:
		# skip spontaneous reactions
		if getattr(met_rxn, 'complex_data', None) is None:
			continue

		if isinstance(met_rxn, coralme.core.reaction.MetabolicReaction) and r.complex_data.id != 'CPLX_dummy':
			key = met_rxn.id.replace('-', '_DASH_').replace('__', '_DASH_').replace(':', '_COLON_')
			# key = met_rxn.id
			key = 'keff_' + key.replace('_FWD_', '_').replace('_REV_', '_')

			matches = [i for i in keffs if key in i]
			# get the direction
			if met_rxn.reverse:
				matches = [i for i in matches if i.endswith('_reverse_priming_keff')]
			else:
				matches = [i for i in matches if i.endswith('_forward_priming_keff')]

			if len(matches) == 1:
				new_keffs[met_rxn.id] = keffs[matches[0]]
			elif len(matches) > 0:
				if len(matches) == len([i for i in matches if key + '_mod_']):
					new_keffs[met_rxn.id] = keffs[matches[0]]
				else:
					log(key, len(matches))
			else:  # len(matches) == 0
				log('no keff found for ' + key)

	return new_keffs

def get_m_to_me_metabolite_mapping(filename = 'm_to_me_mets.txt', m_mets = 'm_name', me_mets = 'me_name') -> dict:
	"""returns a mapping from m metabolites to me metabolites"""
	df = pandas.read_csv(filename, sep = '\t', index_col = None, comment = '#', skip_blank_lines = True)
	return df.set_index(m_mets)[me_mets].dropna().to_dict()

def get_folding_data(filename, index_name = 'gene name', column_prefix = 'Dill_Keq_', column_suffix = 'C') -> pandas.DataFrame:
	"""returns formatted data length-based approximation of protein folding and aggregation data."""
	df = pandas.read_csv(filename, sep = '\t', index_col = None, comment = '#', skip_blank_lines = True)
	df = df.rename(columns = lambda x: x.replace(column_prefix, '').replace(column_suffix, ''))
	return df.set_index(index_name)
