import logging
log = logging.getLogger(__name__)
import coralme

def add_iron_sulfur_modifications(me_model):
	#generic_fes_transfer_complexes = me_model.global_info['complex_cofactors']['generic_fes_transfer_complexes']
	fes_transfers = me_model.global_info['complex_cofactors']['fes_transfers']

	#identify if the user added MetabolicReactions to create modified FeS transfers
	#if not, FormationReaction's are modified during update(),
	#but missing 2fe2s and 4fe4s will be added as troubleshooted reactions
	for fes in ['2fe2s', '4fe4s']:
		# WARNING: If the organism doesn't have iron-sulfur loaders, CPLX_dummy will be added as an iron-sulfur loader
		# WARNING: The troubleshooter will add sink reactions for 4fe4s/2fe2s, making necessary formation reactions
		# See main.py, lines 1931-1936
		components = [ '{:s}_mod_{:s}(1)'.format(x, fes) for x in set(fes_transfers) if x != '' and x != 'CPLX_dummy' ]
		for component in components:
			query = me_model.metabolites.get_by_id(component).reactions
			query = [ x for x in query if x.metabolites[me_model.metabolites.get_by_id(component)] > 0 ]
			rtypes = [ type(x) for x in query ]
			# two type of reactions if user added the correct formation reaction in reaction_matrix.txt
			if coralme.core.reaction.MetabolicReaction in rtypes and coralme.core.reaction.ComplexFormation in rtypes:
				# we need to remove subreactions from process_data
				to_modify = [ x for x in query if isinstance(x, coralme.core.reaction.ComplexFormation)][0]
				me_model.process_data.get_by_id(to_modify.id.replace('formation_', '')).subreactions.pop('mod_{:s}_c'.format(fes))
				# and remove the formation reaction of the modified complex only
				me_model.process_data.remove(to_modify.id.replace('formation_', ''))
				me_model.remove_reactions([to_modify])

	for fes in ['2fe2s', '4fe4s']:
		name = 'generic_{:s}_transfer_complex'.format(fes)
		components = [ '{:s}_mod_{:s}(1)'.format(x, fes) for x in set(fes_transfers) if x != '' ]
		for component in components:
			if not me_model.metabolites.has_id(component):
				me_model.add_metabolites([coralme.core.component.Complex(component)])

		generic_fes_transfer = coralme.core.processdata.GenericData(name, me_model, components)
		generic_fes_transfer.create_reactions()

		# add 2fe2s or 4fe4s to the model if it doesn't exist
		me_model.add_metabolites([coralme.core.component.Metabolite(fes + '_c')])

		# create unloading reactions
		for name in set(fes_transfers):
			if name != '':
				rxn = coralme.core.reaction.MEReaction('_'.join([name, fes, 'unloading']))
				me_model.add_reactions([rxn])
				rxn.add_metabolites({ name + '_mod_' + fes + '(1)': -1, fes + '_c': 1, name: 1 })

	# add fes transfer enzymes to proper modification data
	if me_model.process_data.has_id('mod_2fe2s_c'):
		mod_2fe2s = me_model.process_data.mod_2fe2s_c
		mod_2fe2s.enzyme = ['generic_2fe2s_transfer_complex']
		mod_2fe2s.stoichiometry = { '2fe2s_c': -1 }
		mod_2fe2s._element_contribution = { 'Fe': 2, 'S': 2 }

	if me_model.process_data.has_id('mod_4fe4s_c'):
		mod_4fe4s = me_model.process_data.mod_4fe4s_c
		mod_4fe4s.enzyme = ['generic_4fe4s_transfer_complex']
		mod_4fe4s.stoichiometry = { '4fe4s_c': -1 }
		mod_4fe4s._element_contribution = { 'Fe': 4, 'S': 4 }

	if me_model.process_data.has_id('mod_3fe4s_c'):
		mod_3fe4s = me_model.process_data.mod_3fe4s_c
		mod_3fe4s.enzyme = ['generic_4fe4s_transfer_complex']
		mod_3fe4s.stoichiometry = { '4fe4s_c': -1, 'fe2_c': 1 }
		mod_3fe4s._element_contribution = { 'Fe': 3, 'S': 4 }

	fes_chaperones = me_model.global_info['complex_cofactors'].get('fes_chaperones', {})
	for chaperone in set(fes_chaperones.values()):
		new_mod = coralme.core.processdata.SubreactionData('mod_2fe2s_c_' + chaperone, me_model)
		new_mod.enzyme = [ chaperone, 'generic_2fe2s_transfer_complex' ]
		new_mod.stoichiometry = { '2fe2s_c': -1 }

	if me_model.process_data.has_id('mod_2fe2s_c'):
		for cplx_data in me_model.process_data.get_by_id('mod_2fe2s_c').get_complex_data():
			cplx_id = cplx_data.id.split('_mod')[0]
			if cplx_id in fes_chaperones:
				cplx_data.subreactions[ 'mod_2fe2s_c_' + fes_chaperones[cplx_id] ] = cplx_data.subreactions.pop('mod_2fe2s_c')

def _replace_modification(dct, me_model):
	modification = list(dct.keys())[0] # current modification ID in me.process_data
	new_subrxn_id = list(dct.values())[0][0] # subreaction ID in process_data to replace modification ID

	if me_model.process_data.has_id(new_subrxn_id):
		new_mod_data = me_model.process_data.get_by_id(new_subrxn_id)
	else: # create SubreactionData with new ID (e.g., mod_btn_c now biotin_ligase)
		new_mod_data = coralme.core.processdata.SubreactionData(new_subrxn_id, me_model)
		try:
			new_mod_data._element_contribution = new_mod_data.calculate_element_contribution()
		except:
			logging.warning('All metabolites in SubreactionData \'{:s}\' must have a formula to determine their elemental contribution.'.format(new_mod_data))

	# This function does not replace information in process_data; it rather changes information in subreactions
	for data in me_model.process_data.get_by_id(modification).get_complex_data():
		cplx_data = me_model.process_data.get_by_id(data.id)
		# copy the original data
		cplx_data.complex_id = data.complex_id
		cplx_data.stoichiometry = data.stoichiometry
		cplx_data.subreactions = data.subreactions.copy()
		# e.g.: remove mod_btn_c and replace it with the new subreaction id
		cplx_data.subreactions[new_mod_data.id] = cplx_data.subreactions.pop(modification)
		#cplx_data.create_complex_formation() # the cplx_data already exists

def add_btn_modifications(me_model):
	# dct = { "mod_btn_c" : [ "biotin_ligase" ] }
	dct = me_model.global_info['complex_cofactors'].get('biotin_subreactions', {})
	if bool(dct):
		_replace_modification(dct, me_model)

def add_2tpr3dpcoa_modifications(me_model):
	# dct = { "mod_2tpr3dpcoa_c" : [ "citx_transfer_to_citd" ] }
	dct = me_model.global_info['complex_cofactors'].get('citx_subreactions', {})
	if bool(dct):
		_replace_modification(dct, me_model)

def add_glycyl_modifications(me_model):
	# dct = { "mod_glycyl_c" : [ "gre_activation" ] }
	dct = me_model.global_info['complex_cofactors'].get('glycyl_subreactions', {})

	if dct:
		# WARNING: Do we need to correct all the modification.enzyme properties (biotin, 2tpr3dpcoa, pan4p)?
		data = me_model.subreaction_data.get_by_id(list(dct.values())[0][0])
		for met, stoich in data.stoichiometry.items():
			if isinstance(me_model.metabolites.get_by_id(met), coralme.core.component.Complex) and stoich < 0:
				data.enzyme.add(met)

	if bool(dct):
		_replace_modification(dct, me_model)

def add_pan4p_modifications(me_model):
	# dct = { "mod_pan4p_c" : [ "acpP_activation" ] }
	dct = me_model.global_info['complex_cofactors'].get('acps_subreactions', {})
	if bool(dct):
		_replace_modification(dct, me_model)

def add_FeFe_and_NiFe_modifications(me_model):
	# dct = { "mod_FeFe_cofactor_c" : [ "base_complex" ], "mod_NiFe_cofactor_c" : [ "base_complex" ] }
	dct = me_model.global_info['complex_cofactors'].get('FeFe/NiFe', {})

	for mod, base_complexes in dct.items():
		for base_complex in base_complexes:
			if base_complex != '' and me_model.process_data.has_id(mod):
				complex_data = list(me_model.process_data.get_by_id(mod).get_complex_data())
				if len(complex_data) > 0:
					for data in complex_data:
						cplx_data = me_model.process_data.get_by_id(data.id)
						cplx_data.complex_id = data.complex_id
						cplx_data.stoichiometry = { base_complex + '_' + mod[:-2] + '(1)' : 1 }
						# The modifications occurs now as follow:
						# base_complex + FeFe/NiFe => base_complex_mod_FeFe/NiFe
						# base_complex_mod_FeFe/NiFe + other cofactors => final modified complex
						cplx_data.subreactions[mod] = 0
				else:
					logging.warning('The ID \'{:s}\' in the configuration file has no base complexes assigned to it.'.format(mod))
			else:
				logging.warning('The ID \'{:s}\' in the configuration file does not exist in the ME-model.'.format(mod))

#def add_lipoate_modifications(me_model, lipoate_modifications):
def add_lipoyl_modifications(me_model):
	# Two different reactions can add a lipoate modification.
	# We create a separate SubreactionData for each one
	#for key, mod in lipoate_modifications):
	lipoate_modifications = me_model.global_info['complex_cofactors']['lipoate_subreactions']
	for mod in list(lipoate_modifications.values())[0]:
		if mod in me_model.process_data:
			mod_data = me_model.process_data.get_by_id(mod)
		else:
			# create SubreactionData
			mod_data = coralme.core.processdata.SubreactionData(mod, me_model)

		#info = me_model.process_data.get_by_id(mod)
		#mod_data.enzyme = mod_data.enzyme
		#mod_data.stoichiometry = mod_data.stoichiometry
		# element count for lipoate modifications
		try:
			mod_data._element_contribution = mod_data.calculate_element_contribution()
		except:
			logging.warning('All metabolites in SubreactionData \'{:s}\' must have a formula to determine their elemental contribution.'.format(mod))

	lipoate = me_model.process_data.get_by_id(list(lipoate_modifications.keys())[0])
	#alt_lipo = me_model.process_data.get_by_id('mod_lipo_c_alt')
	for data in lipoate.get_complex_data():
		#alt_cplx_data = coralme.core.processdata.ComplexData(data.id + '_alt', me_model)
		#alt_cplx_data.complex_id = data.complex_id
		#alt_cplx_data.stoichiometry = data.stoichiometry
		#alt_cplx_data.subreactions = data.subreactions.copy()
		#alt_cplx_data.subreactions[alt_lipo.id] = alt_cplx_data.subreactions.pop(lipoate.id)
		#alt_cplx_data.create_complex_formation()

		# We need to replace the `mod_lipoyl_c` (the key in lipoate_modifications) by the values in lipoate_modifications
		for mod in list(lipoate_modifications.values())[0]:
			new_cplx_data = coralme.core.processdata.ComplexData(data.id + '_' + mod, me_model)
			# copy the original data
			new_cplx_data.complex_id = data.complex_id
			new_cplx_data.stoichiometry = data.stoichiometry
			new_cplx_data.subreactions = data.subreactions.copy()
			# remove mod_lipoyl_c and replace it with the new subreaction id
			new_cplx_data.subreactions[mod] = new_cplx_data.subreactions.pop(lipoate.id)
			new_cplx_data.create_complex_formation()

	# Delete old ComplexFormation data
	# TODO: remove entry from process_data
	lst = list(me_model.process_data.get_by_id('mod_lipoyl_c').get_complex_data())
	lst = [ x.formation for x in lst ]
	me_model.remove_reactions(lst)
	#me_model.process_data.remove('mod_lipoyl_c')

def add_bmocogdp_chaperones(me_model):
	bmocogdp_chaperones = me_model.global_info['complex_cofactors']['bmocogdp_chaperones']
	for chaperone in set(bmocogdp_chaperones.values()):
		new_mod = coralme.core.processdata.SubreactionData('mod_bmocogdp_c_' + chaperone, me_model)
		new_mod.enzyme = [chaperone]
		new_mod.stoichiometry = {'bmocogdp_c': -1}

	for cplx_data in me_model.process_data.get_by_id('mod_bmocogdp_c').get_complex_data():
		cplx_id = cplx_data.id.split('_mod')[0]
		if cplx_id in bmocogdp_chaperones:
			cplx_data.subreactions['mod_bmocogdp_c_' + bmocogdp_chaperones[cplx_id]] = cplx_data.subreactions.pop('mod_bmocogdp_c')

# OLD CODE not used anymore
def add_modification_procedures(me_model):
	# add SubreactionData for iron sulfur clusters
	add_iron_sulfur_modifications(me_model)

	# lipoate modifications can be accomplished using two different mechanisms
	add_lipoate_modifications(me_model)

	# bmocogdp modifications have multiple selective chaperones that transfer
	# the metabolite to the target complexes
	add_bmocogdp_modifications(me_model)

