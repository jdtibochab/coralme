import coralme

def add_subreactions_to_model(me_model, subreactions):
	if not me_model.process_data.has_id('atp_hydrolysis_rho'):
		stoichiometry = {'atp_c': -3.0, 'h2o_c': -3.0, 'adp_c': +3.0, 'h_c': +3.0, 'pi_c': +3.0}
		coralme.util.building.add_subreaction_data(
			me_model, modification_id = 'atp_hydrolysis_rho', modification_stoichiometry = stoichiometry, modification_enzyme = None)

	# add subreactions associated with transcription termination
	for subreaction in subreactions:
		for rxn, info in subreaction.items():
			data = coralme.core.processdata.SubreactionData(rxn, me_model)
			data.enzyme = info['enzymes']
			if me_model.global_info['transcription_subreactions'][rxn] == '':
				data.stoichiometry = {}
			else:
				data.stoichiometry = me_model.process_data.get_by_id(me_model.global_info['transcription_subreactions'][rxn]).stoichiometry
			data._element_contribution = data.calculate_element_contribution() #info.get('element_contribution', {})

def add_rna_polymerase_complexes(me_model, rna_polymerase_id_by_sigma_factor, verbose = True):
	for cplx, components in rna_polymerase_id_by_sigma_factor.items():
		if me_model.metabolites.has_id(components['sigma_factor']) and me_model.metabolites.has_id(components['polymerase']):
			rnap_complex = coralme.core.processdata.ComplexData(cplx, me_model)
			sigma_factor = components['sigma_factor']
			polymerase = components['polymerase']

			rnap_components = rnap_complex.stoichiometry
			rnap_components[sigma_factor] = 1.
			rnap_components[polymerase] = 1.

			rnap_complex.create_complex_formation(verbose = verbose)

def add_rna_excision_machinery(me_model, excision_type, stoichiometry):
	#for excision_type in excision_machinery:
	complex_data = coralme.core.processdata.ComplexData(excision_type + '_excision_machinery', me_model)

	#for machine in excision_machinery[excision_type]:
		#complex_data.stoichiometry[machine] = 1
	complex_data.stoichiometry = stoichiometry

	complex_data.create_complex_formation()
	modification = coralme.core.processdata.SubreactionData(excision_type + '_excision', me_model)
	modification.enzyme = [complex_data.id]

def add_rna_splicing(me_model):
	# Loop through transcription reactions and add appropriate splicing
	# machinery based on RNA types and number of splices required
	for data in me_model.transcription_data:
		#n_excised = sum(data.excised_bases.values())
		#n_cuts = len(data.RNA_products) * 2
		n_overlapping = data.n_overlapping
		n_excised = data.n_excised
		n_cuts = data.n_cuts

		#if n_excised == 0 or (n_excised + n_overlapping) == 0 or n_cuts == 0:
			#continue

		rna_types = list(data.RNA_types)
		n_trna = rna_types.count('tRNA')

		if 'rRNA' in set(rna_types):
			data.subreactions['rRNA_containing_excision'] = n_cuts
		elif n_trna == 1:
			data.subreactions['monocistronic_excision'] = n_cuts
		elif n_trna > 1:
			data.subreactions['polycistronic_wout_rRNA_excision'] = n_cuts
		else: # only applies to rnpB (RNase P catalytic RNA component)
			data.subreactions['monocistronic_excision'] = n_cuts

		# The non functional RNA segments need degraded back to nucleotides
		# TODO check if RNA_degradation requirement is per nucleotide
		data.subreactions['RNA_degradation_machine'] = n_cuts
		data.subreactions['RNA_degradation_atp_requirement'] = n_excised + n_overlapping
