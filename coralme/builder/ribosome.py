import logging
log = logging.getLogger(__name__)

import coralme

#def add_ribosome(me_model, ribosome_stoich, ribosome_subreactions, rrna_mods, modification_info, verbose = True):
def add_ribosome(me_model, ribosome_stoich, ribosome_subreactions, rrna_mods, verbose = True):
	lst = ['generic_Era', 'generic_RbfA', 'generic_RimM']
	for ComplexData in lst:
		if len(me_model.process_data.query(ComplexData)) == 0:
			logging.warning('ComplexData \'{:s}\' not in the ME-model. The ME-model would be unfeasible.'.format(ComplexData))
			logging.warning('Check if the behavior is correct or if \'{:s}\' is present in the \'Generic Complex ID\' column of the organism-specific matrix.'.format(ComplexData.replace('generic_', '')))

	ribosome_complex = coralme.core.processdata.ComplexData(me_model.global_info['ribosome_id'], me_model)

	for idx, mod_data in rrna_mods.iterrows():
		for position in mod_data.positions.split(','):
			rrna_mod = coralme.core.processdata.SubreactionData('{:s}_at_{:s}'.format(mod_data.modification, position), me_model)
			rrna_mod.enzyme = mod_data.enzymes.split(' AND ') if mod_data.enzymes != 'No_Machine' else ['CPLX_dummy']
			#rrna_mod.stoichiometry = modification_info[mod_data.modification]['metabolites']
			rrna_mod.stoichiometry = me_model.process_data.get_by_id(mod_data.modification).stoichiometry
			rrna_mod.keff = 65. # iOL uses 65. for all RNA mods

			# Add element contribution from modification to rRNA
			#rrna_mod._element_contribution = modification_info[mod_data.modification]['elements']
			rrna_mod._element_contribution = me_model.process_data.get_by_id(mod_data.modification).calculate_element_contribution()

			#if 'carriers' in modification_info[mod_data.modification].keys():
				#for carrier, stoich in modification_info[mod_data.modification]['carriers'].items():
					#if stoich < 0:
						#rrna_mod.enzyme += [carrier]
					#rrna_mod.stoichiometry[carrier] = stoich
			ribosome_complex.subreactions[rrna_mod.id] = 1.

	for subreaction_id in ribosome_subreactions:
		if not me_model.process_data.has_id('gtp_hydrolysis_era'):
			stoichiometry = { 'gtp_c': -2.0, 'h2o_c': -2.0, 'gdp_c': 2.0, 'h_c': 2.0, 'pi_c': 2.0 }
			coralme.util.building.add_subreaction_data(
				me_model, modification_id = 'gtp_hydrolysis_era', modification_stoichiometry = stoichiometry, modification_enzyme = None)

		# add subreaction to model
		if me_model.process_data.has_id(subreaction_id):
			subreaction = me_model.process_data.get_by_id(subreaction_id)
		else:
			subreaction = coralme.core.processdata.SubreactionData(subreaction_id, me_model)
			#subreaction.stoichiometry = ribosome_subreactions[subreaction_id]['stoich']
			#reaction_id = me_model.global_info['translation_subreactions'][subreaction_id]
			reaction_id = me_model.global_info['translation_subreactions'].get(subreaction_id,None)
			if reaction_id == '' or reaction_id is None:
				# If reaction_id is not in global_info it must have been defined in subreaction_matrix
				subreaction.stoichiometry = {}
			else:
				subreaction.stoichiometry = me_model.process_data.get_by_id(reaction_id).stoichiometry

		subreaction.enzyme = ribosome_subreactions[subreaction_id]['enzymes']
		# account for subreactions in complex data. num_mods is always 1
		ribosome_complex.subreactions[subreaction.id] = 1.

	# Ribosomes in iOL1650 contain 171 mg2 ions
	ribosome_complex.subreactions['mod_mg2_c'] = me_model.global_info['mg2_per_ribosome']
	ribosome_components = ribosome_complex.stoichiometry
	for process in ribosome_stoich:
		for protein, amount in ribosome_stoich[process]['stoich'].items():
			ribosome_components[protein] += amount

	ribosome_complex.create_complex_formation(verbose = verbose)
