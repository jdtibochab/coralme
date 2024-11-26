import Bio
import logging
log = logging.getLogger(__name__)

import coralme

#from functools import reduce
#def reducer(accumulator, element):
	#for key, value in element.items():
		#accumulator[key] = accumulator.get(key, 0) + value
	#return accumulator

from collections import Counter

def add_subreactions_to_model(me_model, subreactions):
	if not me_model.process_data.has_id('DEF'):
		stoichiometry = { 'h2o_c': -1.0, 'for_c': +1.0 }
		coralme.util.building.add_subreaction_data(
			me_model, modification_id = 'DEF', modification_stoichiometry = stoichiometry, modification_enzyme = 'CPLX_dummy')

	if not me_model.process_data.has_id('MAP'):
		stoichiometry = { 'h2o_c': -1.0, 'h_c': +1.0, 'met__L_c': +1.0 }
		coralme.util.building.add_subreaction_data(
			me_model, modification_id = 'MAP', modification_stoichiometry = stoichiometry, modification_enzyme = 'CPLX_dummy')

	if not me_model.process_data.has_id('FMETTRS'):
		stoichiometry = {
			'generic_tRNA_START_met__L_c': -1.0,
			me_model.metabolites.query('^10fthf')[0].id: -1.0,
			me_model.metabolites.query('^thf')[0].id: +1.0
			}
		coralme.util.building.add_subreaction_data(
			me_model, modification_id = 'FMETTRS', modification_stoichiometry = stoichiometry, modification_enzyme = 'CPLX_dummy')

	if not me_model.process_data.has_id('atp_hydrolysis_groel'):
		stoichiometry = { 'atp_c': -7.0, 'h2o_c': -7.0, 'adp_c': +7.0, 'h_c': +7.0, 'pi_c': +7.0 }
		coralme.util.building.add_subreaction_data(
			me_model, modification_id = 'atp_hydrolysis_groel', modification_stoichiometry = stoichiometry, modification_enzyme = None)

	# add subreactions associated with translation initiation, elongation, termination and postprocessing
	for subreaction in subreactions:
		for rxn, info in subreaction.items():
			data = coralme.core.processdata.SubreactionData(rxn, me_model)
			data.enzyme = info['enzymes']
			if me_model.global_info['translation_subreactions'][rxn] == '':
				data.stoichiometry = {}
			else:
				data.stoichiometry = me_model.process_data.get_by_id(me_model.global_info['translation_subreactions'][rxn]).stoichiometry
				data._element_contribution = data.calculate_element_contribution() #info.get('element_contribution', {})

	#old code
	#for rxn, info in termination_subreactions.items():
		#data = coralme.core.processdata.SubreactionData(rxn, me_model)
		#data.enzyme = info['enzymes']
		#data.stoichiometry = info['stoich']
		#data._element_contribution = info.get('element_contribution', {})

def add_charged_trna_subreactions(me_model, organelle = 'c', transl_table = set([11]), translation_stop_dict = {}, selenocysteine_enzymes = []):
	"""
	Create subreaction for each codon. this will be used to model
	the addition of charged tRNAs to the elongating peptide
	"""
	#for codon in coralme.util.dogma.codon_table:
	for codon in me_model.global_info['stop_codons']:
		#if coralme.util.dogma.codon_table[codon] == '*':
		stop_codon = codon.replace('T', 'U')
		stop_enzyme = translation_stop_dict.get(stop_codon, 'CPLX_dummy')
		me_model.add_metabolites([coralme.core.component.Complex(stop_enzyme)])
		subreaction_data = coralme.core.processdata.SubreactionData(stop_codon + '_' + stop_enzyme + '_mediated_termination_' + organelle, me_model)
		subreaction_data.enzyme = [stop_enzyme]
		subreaction_data.stoichiometry = {}
		#else:

	codon_table = Bio.Data.CodonTable.generic_by_id[list(transl_table)[0]]

	# check if the defined amino_acid_loader (default generic_Tuf) exists in the model
	if not me_model.metabolites.has_id(me_model.global_info['amino_acid_loader']):
		me_model.global_info['amino_acid_loader'] = 'CPLX_dummy'

	#for codon, aa in { k:v for k,v in me_model.global_info['codon_table'].forward_table.items() if 'U' not in k }.items():
	dct = { k:v for k,v in codon_table.forward_table.items() if 'U' not in k }
	if me_model.global_info.get('genetic_recoding', {}):
		# codon_table is codon to one-letter code amino acid
		for codon, tRNA_to_aa in me_model.global_info['genetic_recoding'].items():
			for aa, tRNA in tRNA_to_aa.items():
				dct.update({ codon : aa })
				logging.warning('Genetic recoding dictionary added to internal codon table.')

	for codon, aa in dct.items():
		if not me_model.process_data.has_id('atp_hydrolysis_trna_loading'):
			stoichiometry = { 'atp_c': -1.0, 'h2o_c': -1.0, 'amp_c': +1.0, 'h_c': +1.0, 'ppi_c': +1.0 }
			coralme.util.building.add_subreaction_data(
				me_model, modification_id = 'atp_hydrolysis_trna_loading', modification_stoichiometry = stoichiometry, modification_enzyme = None)

		if not me_model.process_data.has_id('gtp_hydrolysis'):
			stoichiometry = { 'gtp_c': -1.0, 'h2o_c': -1.0, 'gdp_c': +1.0, 'h_c': +1.0, 'pi_c': +1.0 }
			coralme.util.building.add_subreaction_data(
				me_model, modification_id = 'gtp_hydrolysis', modification_stoichiometry = stoichiometry, modification_enzyme = None)

		#full_aa = coralme.util.dogma.amino_acids[coralme.util.dogma.codon_table[codon]]
		if aa in coralme.util.dogma.amino_acids:
			full_aa = coralme.util.dogma.amino_acids[aa] # now without the compartment code (e.g., full_aa = gly; not gly_c)
		else:
			full_aa = aa.replace('_c', '') # do not add compartment (e.g., full_aa = ala__L; not ala__L_c)

		if not me_model.metabolites.has_id(full_aa + '_' + organelle):
			logging.warning('The amino acid \'{:s}\' in not in the ME-model.'.format(full_aa + '_' + organelle))
			continue

		#full_aa = full_aa.split('_')[0]
		subreaction_id = full_aa.split('_')[0] + '_addition_at_' + codon.replace('T', 'U')
		#subreaction_id = full_aa + '_' + organelle + '_addition_at_' + codon.replace('T', 'U')

		subreaction_data = coralme.core.processdata.SubreactionData(subreaction_id, me_model)
		trna = 'generic_tRNA_' + codon.replace('T', 'U') + '_' + full_aa + '_' + organelle
		# Default AA loader enzyme
		subreaction_data.enzyme = [me_model.global_info['amino_acid_loader']]
		# Accounts for GTP hydrolyzed by EF-TU and the ATP hydrolysis to AMP required to add the amino acid to the tRNA
		#subreaction_data.stoichiometry = reduce(
			#reducer, [{trna : -1}, me_model.global_info['atp_trna_loading'], me_model.global_info['gtp_hydrolysis']], {})
		subreaction_data.stoichiometry = Counter({trna : -1})
		subreaction_data.stoichiometry.update(me_model.process_data.get_by_id('atp_hydrolysis_trna_loading').stoichiometry)
		subreaction_data.stoichiometry.update(me_model.process_data.get_by_id('gtp_hydrolysis').stoichiometry)
		subreaction_data._element_contribution = me_model.metabolites.get_by_id(full_aa + '_' + organelle).elements # subreaction_data.calculate_element_contribution()
		logging.warning('A SubreactionData with ID \'{:s}\' was added to the ME-model.'.format(subreaction_id))

	# Add subreactions for selenocysteine
	if me_model.metabolites.has_id('generic_tRNA_UGA_ser__L_c'):
		subreaction_data = coralme.core.processdata.SubreactionData('sec_addition_at_UGA', me_model)
		subreaction_data.enzyme = selenocysteine_enzymes # Do not add brackets [] here, it is already a list
		subreaction_data.stoichiometry = { 'generic_tRNA_UGA_ser__L_c': -1.0, 'selnp_c': -1.0, 'h_c': +1.0, 'pi_c': +1.0 }
		subreaction_data._element_contribution = subreaction_data.calculate_element_contribution()
		logging.warning('A SubreactionData with ID \'{:s}\' was added to the ME-model.'.format('sec_addition_at_UGA'))
