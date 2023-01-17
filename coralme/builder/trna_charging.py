import coralme

#def add_trna_modification_procedures(me_model, trna_mods, modification_info):
def add_trna_modification_procedures(me_model, trna_mods):
	# remove duplications
	trna_mods = trna_mods.drop_duplicates(['modification', 'positions'], keep = 'first')
	for idx, mod_data in trna_mods.iterrows():
		for position in mod_data['positions'].split(','):
			#if mod_data.type == 'met_tRNA':
			if mod_data['bnum'] in me_model.global_info['START_tRNA']:
				#name = '{:s}_at_{:s}_in_{:s}'.format(mod_data.modification, position, mod_data.type)
				name = '{:s}_at_{:s}_in_{:s}'.format(mod_data['modification'], position, mod_data['bnum'])
			else:
				name = '{:s}_at_{:s}'.format(mod_data['modification'], position)

			trna_mod = coralme.core.processdata.SubreactionData(name, me_model)
			trna_mod.enzyme = mod_data['enzymes'].split(' AND ') if mod_data['enzymes'] != 'No_Machine' else None
			#trna_mod.stoichiometry = modification_info[mod_data.modification]['metabolites']
			try:
				trna_mod.stoichiometry = me_model.process_data.get_by_id(mod_data.modification).stoichiometry
			except:
				trna_mod.stoichiometry = {}
			trna_mod.keff = 65.  # iOL uses 65 for all tRNA mods

			# TODO: Modify to identify 'carriers' in SubReaction data
			#if 'carriers' in modification_info[mod_data.modification].keys():
				#for carrier, stoich in modification_info[mod_data.modification]['carriers'].items():
					#if stoich < 0:
						#trna_mod.enzyme += [carrier]
					#trna_mod.stoichiometry[carrier] = stoich

		# Add element contribution from modification to tRNA
		#trna_mod._element_contribution = modification_info[mod_data.modification]['elements']
		try:
			trna_mod._element_contribution = me_model.process_data.get_by_id(mod_data['modification']).calculate_element_contribution()
		except:
			trna_mod._element_contribution = {}